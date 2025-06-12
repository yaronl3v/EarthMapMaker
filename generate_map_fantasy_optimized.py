#!/usr/bin/env python3
"""
Real-continent Voronoi world-map generator
– removes small islands
– ~5 000 even-ish polygons
– curvy shared borders
– merges slivers
"""

import random, math, hashlib, json, pathlib
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon, mapping
from shapely.ops import unary_union
from shapely.errors import TopologicalError
from scipy.spatial import Voronoi
from shapely.strtree import STRtree

# Utility: build spatial index -------------------------------------------------

def build_spatial_index(polys: list[Polygon | MultiPolygon]):
    """Return (STRtree, wkb→index mapping) for fast spatial queries that survives geometry copying."""
    tree = STRtree(polys)
    wkb_map = {g.wkb: i for i, g in enumerate(polys)}
    return tree, wkb_map

# ---------- tweakables ----------
SEED              = 42
ISLAND_MIN_KM2    = 8_000     # keep UK=244 k, drop Japan=378 k? set lower if you want them
N_SEEDS           = 60_000
LLOYD_ITER        = 0            # 0 → keep more irregular seed layout
MAX_SEG_LEN       = 0.20         # ° before we insert extra vertices
MAX_OFFSET        = 0.08         # ° max border wiggle
USE_CURVES        = False        # toggle edge warping on/off
SLIVER_FACTOR     = 0.3          # merge cells smaller than avg*factor before starting merging
SLIVER_FACTOR2    = 0.5          # merge cells smaller than avg*factor in the final step
OUT               = "world_cells.geojson"
GROUP_MIN         = 3            # min polygons per funny group
GROUP_MAX         = 6            # max polygons per funny group
MERGE_PROBS       = {
    "random": 0.50,
    "vertical": 0.17,
    "horizontal": 0.17,
    "diagonal": 0.16,
}
# ---------------------------------

random.seed(SEED); np.random.seed(SEED)

# 1. load land, drop Antarctica, drop tiny islands
world = gpd.read_file("raw map/ne_110m_admin_0_countries.shp")
world = world[world["CONTINENT"] != "Antarctica"]
world["area_km2"] = world.to_crs(8857).area / 1e6          # equal-area
world = world[world["area_km2"] >= ISLAND_MIN_KM2]
land = unary_union(world.geometry)

# 2. random seed points --------------------------------------------------
def sample(poly, n, batch=20_000):
    minx, miny, maxx, maxy = poly.bounds
    pts = []
    while len(pts) < n:
        xs = np.random.uniform(minx, maxx, batch)
        ys = np.random.uniform(miny, maxy, batch)
        for x, y in zip(xs, ys):
            if poly.contains(Point(x, y)):
                pts.append((x, y))
                if len(pts) == n: break
    return np.array(pts)

seeds = sample(land, N_SEEDS)

# 3. Lloyd relaxation ----------------------------------------------------
def relax(pts, mask, iters):
    for _ in range(iters):
        vor = Voronoi(pts)
        new = []
        for r in vor.point_region:
            reg = vor.regions[r]
            if not reg or -1 in reg: continue
            poly = Polygon(vor.vertices[reg]).intersection(mask)
            if poly.is_valid and not poly.is_empty:
                c = poly.centroid; new.append((c.x, c.y))
        pts = np.array(new)
    return pts

seeds = relax(seeds, land, LLOYD_ITER)

# 4. Voronoi + clip ------------------------------------------------------
vor = Voronoi(seeds)
cells = []
for r in vor.point_region:
    reg = vor.regions[r]
    if not reg or -1 in reg: continue

    # Clip voronoi cell to land mass. .buffer(0) helps fix invalid geometries.
    cell = Polygon(vor.vertices[reg]).intersection(land).buffer(0)

    if cell.is_empty: continue

    # Explode MultiPolygons into individual polygons. This ensures that when a
    # voronoi cell is split by a coastline, each resulting piece is treated
    # as a separate entity, allowing us to find and merge tiny slivers.
    if cell.geom_type == 'Polygon':
        cells.append(cell)
    elif cell.geom_type == 'MultiPolygon':
        cells.extend([p for p in cell.geoms if p.area > 0])

# 5. edge wiggle ---------------------------------------------------------
def key(a,b): return tuple(sorted((a,b)))
edge_cache={}
def warp(a,b):
    dx,dy=b[0]-a[0],b[1]-a[1]; L=math.hypot(dx,dy)
    n=max(1,int(L/MAX_SEG_LEN)); h=hashlib.sha256(str(key(a,b)).encode()).digest()
    rng=random.Random(int.from_bytes(h[:4],"little"))
    phase=rng.uniform(0,2*math.pi); amp=rng.uniform(0.4,1)*MAX_OFFSET
    ux,uy=(-dy/L,dx/L) if L else (0,0); pts=[a]
    for i in range(1,n):
        t=i/n; pts.append((a[0]+dx*t+ux*amp*math.sin(phase+t*math.pi),
                           a[1]+dy*t+uy*amp*math.sin(phase+t*math.pi)))
    pts.append(b); return pts

def curvy(poly):
    if poly.geom_type=="Polygon":
        ring=list(poly.exterior.coords)
        new=[]
        for p,q in zip(ring,ring[1:]):
            seg=edge_cache.setdefault(key(p,q),warp(p,q))
            if seg[0]!=p: seg=list(reversed(seg))
            new.extend(seg[:-1])
        new.append(ring[-1])
        return Polygon(new)
    else:
        return MultiPolygon([curvy(p) for p in poly.geoms])

if USE_CURVES:
    cells = [curvy(c) for c in cells]

# After warping edges, clip again to the landmass to ensure no segment spills across the coastline
reclipped_cells = []
for cell in cells:
    # Ensure geometry validity before attempting costly operations
    if not cell.is_valid:
        cell = cell.buffer(0)

    try:
        clipped = cell.intersection(land).buffer(0)  # buffer(0) fixes potential topology errors
    except TopologicalError:
        # Skip problematic geometries that still fail after cleaning
        continue

    if clipped.is_empty:
        continue
    if clipped.geom_type == "Polygon":
        reclipped_cells.append(clipped)
    elif clipped.geom_type == "MultiPolygon":
        reclipped_cells.extend([p for p in clipped.geoms if p.area > 0])

cells = reclipped_cells  # replace with coast‐compliant versions

# 6. merge slivers -------------------------------------------------------
MERGE_PASSES = 5
print(f"Starting merge process for up to {MERGE_PASSES} passes.")
initial_polygon_count = len(cells)

def find_touching(i_cell: int, skip: set[int], tree: STRtree, id_map: dict[int, int]) -> list[int]:
    """Return indices of polygons that touch the target cell and are not in *skip*."""
    target_geom = cells[i_cell]
    candidates = tree.query(target_geom)
    idxs: list[int] = []
    for c in candidates:
        # Shapely ≥2 returns integer indexes; older versions return geometry objects
        if isinstance(c, (int, np.integer)):
            idxs.append(int(c))
        else:
            idx = id_map.get(c.wkb)
            if idx is not None:
                idxs.append(idx)

    return [j for j in idxs if j is not None and j != i_cell and j not in skip and target_geom.touches(cells[j])]

def safe_union(source_idx: int, target_idx: int) -> bool:
    """Try to union *source* into *target*.

    Returns True on success, False if the union raised a *TopologicalError*."""
    try:
        source = cells[source_idx].buffer(0)
        target = cells[target_idx].buffer(0)
        cells[target_idx] = source.union(target).buffer(0)
        return True
    except TopologicalError:
        return False

# Compute a fixed average cell area once, and keep the threshold constant across all merge passes
avg_area_fixed = sum(c.area for c in cells) / len(cells)
threshold = avg_area_fixed * SLIVER_FACTOR
print(f"Fixed average cell area: {avg_area_fixed:.4f} → merge threshold (fixed): {threshold:.4f}")

for i_pass in range(MERGE_PASSES):
    print(f"\n--- Merge pass {i_pass + 1}/{MERGE_PASSES} ---")
    small_cells = [idx for idx, c in enumerate(cells) if c.area < threshold]

    if not small_cells:
        print("No small cells left – merge complete.")
        break

    print(f"Small cells to process this pass: {len(small_cells)}")

    # Build spatial index for this pass
    tree, id_map = build_spatial_index(cells)

    removed: set[int] = set()
    merged_any = False

    for s_idx in small_cells:
        if s_idx in removed:
            continue  # already merged away this pass

        neighbours = find_touching(s_idx, removed, tree, id_map)
        if not neighbours:
            # fall back to closest polygon by centroid distance if none touch
            distances = [(j, cells[s_idx].distance(cells[j])) for j in range(len(cells)) if j != s_idx and j not in removed]
            neighbours = [min(distances, key=lambda x: x[1])[0]] if distances else []

        for n_idx in neighbours:
            if n_idx in removed:
                continue
            if safe_union(s_idx, n_idx):
                removed.add(s_idx)
                merged_any = True
                break  # success – go to next small cell

    print(f"Merged {len(removed)} small polygons in this pass.")

    cells = [c for k, c in enumerate(cells) if k not in removed]

    if not merged_any:
        # Could not merge any remaining small cells – likely due to topology issues
        print("No successful merges this pass. Stopping early to avoid infinite loop.")
        break

# 7. funny group merge ---------------------------------------------------
print("\nStarting funny grouping phase …")

rng = random.Random(SEED)

centroids = [c.centroid for c in cells]


def pick_merge_type() -> str:
    roll = rng.random()
    accum = 0.0
    for k, p in MERGE_PROBS.items():
        accum += p
        if roll <= accum:
            return k
    return "random"  # fallback


def direction_ok(base_idx: int, cand_idx: int, kind: str) -> bool:
    if kind == "random":
        return True
    dx = centroids[cand_idx].x - centroids[base_idx].x
    dy = centroids[cand_idx].y - centroids[base_idx].y
    if kind == "vertical":
        return abs(dy) > abs(dx)
    if kind == "horizontal":
        return abs(dx) > abs(dy)
    # diagonal
    return abs(abs(dx) - abs(dy)) < max(abs(dx), abs(dy)) * 0.3  # roughly diagonal


assigned: set[int] = set()
funny_cells: list[Polygon | MultiPolygon] = []

# Build spatial index once for neighbour queries during funny grouping
funny_tree, funny_id_map = build_spatial_index(cells)


def neighbours_of(i: int) -> list[int]:
    tgt = cells[i]
    cand_geoms = funny_tree.query(tgt.buffer(1e-5))  # slight expand to catch edge neighbours
    idxs: list[int] = []
    for c in cand_geoms:
        if isinstance(c, (int, np.integer)):
            idxs.append(int(c))
        else:
            idx = funny_id_map.get(c.wkb)
            if idx is not None:
                idxs.append(idx)
    return [j for j in idxs if j is not None and j != i and j not in assigned and (
        tgt.distance(cells[j]) < 1e-6)
    ]

for idx in range(len(cells)):
    if idx in assigned:
        continue
    target_size = rng.randint(GROUP_MIN, GROUP_MAX)
    merge_kind = pick_merge_type()
    group = {idx}
    assigned.add(idx)

    while len(group) < target_size:
        # collect neighbour candidates respecting direction
        cand = []
        for g_idx in list(group):
            for n_idx in neighbours_of(g_idx):
                if n_idx in assigned:
                    continue
                if direction_ok(g_idx, n_idx, merge_kind):
                    cand.append(n_idx)
        if not cand:  # fallback to any neighbour
            for g_idx in list(group):
                cand.extend(neighbours_of(g_idx))
            cand = [c for c in cand if c not in assigned]
        if not cand:
            break  # no more neighbours available
        pick = rng.choice(cand)
        group.add(pick)
        assigned.add(pick)

    # union geometries safely
    try:
        merged_poly = unary_union([cells[i].buffer(0) for i in group]).buffer(0)
        funny_cells.append(merged_poly)   
    except TopologicalError:
        # on failure, fall back to individual cells
        funny_cells.extend([cells[i] for i in group])

print(f"Formed {len(funny_cells)} funny groups from {len(cells)} base cells.")

final = funny_cells

# 8. final sliver cleanup -----------------------------------------------
print("\nFinal sliver cleanup …")

cells_cleanup = final  # work on a mutable reference
avg_area_after = sum(p.area for p in cells_cleanup) / len(cells_cleanup)
cleanup_threshold = avg_area_after * SLIVER_FACTOR2
print(f"Average area after funny grouping: {avg_area_after:.4f} → cleanup threshold: {cleanup_threshold:.4f}")

# Build spatial index once for this cleanup phase (much faster than full scans)
cleanup_tree, cleanup_id_map = build_spatial_index(cells_cleanup)

removed: set[int] = set()

def touch_indices(i: int) -> list[int]:
    """Return indices of polygons touching cell *i* and not yet *removed*."""
    tgt = cells_cleanup[i]
    cand_geoms = cleanup_tree.query(tgt)
    idxs: list[int] = []
    for g in cand_geoms:
        if isinstance(g, (int, np.integer)):
            j = int(g)
        else:
            j = cleanup_id_map.get(g.wkb)
        if j is None or j == i or j in removed:
            continue
        if tgt.touches(cells_cleanup[j]):
            idxs.append(j)
    return idxs

for idx, poly in enumerate(cells_cleanup):
    if idx in removed:
        continue
    if poly.area >= cleanup_threshold:
        continue
    neigh = touch_indices(idx)
    if not neigh:
        continue  # isolated tiny piece – leave it
    # merge with first neighbour (could be improved with area sorting)
    target_idx = neigh[0]
    try:
        merged = cells_cleanup[target_idx].union(poly).buffer(0)
        cells_cleanup[target_idx] = merged
        removed.add(idx)
    except TopologicalError:
        continue  # skip on failure

final = [p for k, p in enumerate(cells_cleanup) if k not in removed]
print(f"Removed {len(removed)} leftover small polygons.")

# 10. coastline enforcement --------------------------------------------
print("\nClipping final polygons to coastline boundary …")
clean_final: list[Polygon | MultiPolygon] = []
for poly in final:
    try:
        clipped = poly.intersection(land).buffer(0)
    except TopologicalError:
        continue
    if clipped.is_empty:
        continue
    if clipped.geom_type == "Polygon":
        clean_final.append(clipped)
    elif clipped.geom_type == "MultiPolygon":
        clean_final.extend([p for p in clipped.geoms if p.area > 0])

final = clean_final
print(f"After coastline clip: {len(final)} polygons remain.")

# 11. prune isolated tiny polygons ------------------------------------
print("\nPruning isolated tiny polygons …")

if final:  # guard against empty list
    avg_area_final = sum(p.area for p in final) / len(final)
    island_threshold = avg_area_final * SLIVER_FACTOR2  # reuse secondary sliver factor
    print(f"Average area after coastline clip: {avg_area_final:.4f} → island drop threshold: {island_threshold:.4f}")

    # Build spatial index for neighbour queries
    iso_tree, iso_id_map = build_spatial_index(final)

    kept: list[Polygon | MultiPolygon] = []
    dropped = 0

    for idx, poly in enumerate(final):
        # Check if *poly* touches any other polygon
        neighbours_geoms = iso_tree.query(poly)
        touches_any = False
        for g in neighbours_geoms:
            if isinstance(g, (int, np.integer)):
                j = int(g)
            else:
                j = iso_id_map.get(g.wkb)
            if j is None or j == idx:
                continue
            if poly.touches(final[j]):
                touches_any = True
                break

        # Drop if isolated AND smaller than threshold
        if (not touches_any) and (poly.area < island_threshold):
            dropped += 1
            continue

        kept.append(poly)

    final = kept
    print(f"Dropped {dropped} isolated small polygons.")
else:
    print("No polygons left to prune.")

# 9. save ---------------------------------------------------------------
json_path=pathlib.Path(OUT)
json_path.write_text(json.dumps({
    "type":"FeatureCollection",
    "features":[{"type":"Feature","geometry":mapping(g),"properties":{}} for g in final]
}))
print(f"Saved {len(final)} polygons after merge → {OUT}")

