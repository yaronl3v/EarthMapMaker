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

# ---------- tweakables ----------
SEED              = 42
ISLAND_MIN_KM2    = 50_000     # keep UK=244 k, drop Japan=378 k? set lower if you want them
N_SEEDS           = 5_000
LLOYD_ITER        = 2
MAX_SEG_LEN       = 0.20       # ° before we insert extra vertices
MAX_OFFSET        = 0.08       # ° max border wiggle
SLIVER_FACTOR     = 0.5        # merge cells smaller than avg*factor
OUT               = "world_cells.geojson"
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

cells=[curvy(c) for c in cells]

# 6. merge slivers -------------------------------------------------------
MERGE_PASSES = 5
print(f"Starting merge process for up to {MERGE_PASSES} passes.")
initial_polygon_count = len(cells)

def find_touching(i_cell: int, skip: set[int]) -> list[int]:
    """Return indices of polygons that touch the target cell and are not in *skip*."""
    target = cells[i_cell]
    return [j for j, other in enumerate(cells) if j != i_cell and j not in skip and target.touches(other)]

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

for i_pass in range(MERGE_PASSES):
    print(f"\n--- Merge pass {i_pass + 1}/{MERGE_PASSES} ---")
    avg_area = sum(c.area for c in cells) / len(cells)
    threshold = avg_area * SLIVER_FACTOR
    small_cells = [idx for idx, c in enumerate(cells) if c.area < threshold]

    if not small_cells:
        print("No small cells left – merge complete.")
        break

    print(f"Average cell area: {avg_area:.4f} → merge threshold: {threshold:.4f}")
    print(f"Small cells to process this pass: {len(small_cells)}")

    removed: set[int] = set()
    merged_any = False

    for s_idx in small_cells:
        if s_idx in removed:
            continue  # already merged away this pass

        neighbours = find_touching(s_idx, removed)
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

final = cells
print(f"\nFinished merging. Started with {initial_polygon_count}, ended with {len(final)} polygons.")

# 7. save ---------------------------------------------------------------
json_path=pathlib.Path(OUT)
json_path.write_text(json.dumps({
    "type":"FeatureCollection",
    "features":[{"type":"Feature","geometry":mapping(g),"properties":{}} for g in final]
}))
print(f"Saved {len(final)} polygons after merge → {OUT}")

"""
Adjust:

What you want	Parameter(s)
Keep UK & Japan	lower ISLAND_MIN_KM2
More / fewer cells	change N_SEEDS
Straighter / wigglier borders	tweak MAX_SEG_LEN, MAX_OFFSET
Merge more small cells	raise SLIVER_FACTOR
Flatter cell sizes	raise LLOYD_ITER

That's everything—use, modify, repeat.
""" 