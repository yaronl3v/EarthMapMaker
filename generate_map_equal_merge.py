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
from shapely.errors import TopologicalError, GEOSException
from scipy.spatial import Voronoi

# ---------- tweakables ----------
SEED              = 42
ISLAND_MIN_KM2    = 8_000      # drop only very small islands (keeps Israel, Denmark, etc.)
N_SEEDS           = 20_000
LLOYD_ITER        = 2          # fewer relaxations → less uniform seed spacing
MAX_SEG_LEN       = 0.30       # ° before we insert extra vertices – increase for smoother
MAX_OFFSET        = 0.08       # ° max border wiggle – gentler amplitude
SLIVER_FACTOR     = 0.6        # merge cells smaller than avg*factor
OUT               = "world_cells.geojson"
JITTER_STD        = 0.       # stronger jitter for more varied cell shapes/sizes
MIN_WIGGLES       = 1          # minimum full sine waves per border segment
MAX_WIGGLES       = 2          # maximum full sine waves per border segment
WARP_STYLE        = "angular"     # "sine" for smooth curves, "angular" for faceted edges
NOISE_FREQ        = 0.0       # slightly higher frequency for more local variation
NOISE_MIN_PROB    = 0.0       # baseline chance – sharper contrast between sparse and dense areas
MERGE_NEIGH_MIN   = 3          # each polygon will merge this many neighbours (min)
MERGE_NEIGH_MAX   = 6          # ... and at most this many
MERGE_SKIP_LARGE_AREA = 6.5    # deg² – polygons larger than this won't be merged into (but can absorb others)
SLIVER_MIN_AREA = 0.4           # deg² – absolute cutoff for sliver merging
SMALL_FIRST_AREA   = 1.0        # deg² – always process cells <= this area first
SMOOTH_ITER        = 0          # 0 = off; >0 = apply Chaikin smoothing iterations to all polygons
# ---------------------------------

random.seed(SEED); np.random.seed(SEED)

# 1. load land, drop Antarctica, drop tiny islands
world = gpd.read_file("raw map/ne_110m_admin_0_countries.shp")
world = world[world["CONTINENT"] != "Antarctica"]
world["area_km2"] = world.to_crs(8857).area / 1e6          # equal-area
world = world[world["area_km2"] >= ISLAND_MIN_KM2]
land = unary_union(world.geometry)

# 2. random seed points --------------------------------------------------
def _pseudo_noise(x: float, y: float, freq: float = NOISE_FREQ) -> float:
    """Cheap deterministic pseudo-Perlin-like noise in 1…1 for (lon, lat)."""
    return (
        math.sin(x * freq) +
        math.sin(y * freq * 1.3) +
        math.sin((x + y) * freq * 0.7)
    ) / 3.0


def sample(poly, n, batch: int = 20_000) -> np.ndarray:
    """Rejection-sample *n* points inside *poly* with probability skewed by noise.

    Areas where \_pseudo_noise is near +1 get almost every candidate accepted,
    while areas near −1 still keep at least *NOISE_MIN_PROB* of candidates.
    This varies seed density → variable-sized Voronoi cells."""
    minx, miny, maxx, maxy = poly.bounds
    pts: list[tuple[float, float]] = []
    while len(pts) < n:
        xs = np.random.uniform(minx, maxx, batch)
        ys = np.random.uniform(miny, maxy, batch)
        for x, y in zip(xs, ys):
            if not poly.contains(Point(x, y)):
                continue

            noise_val = _pseudo_noise(x, y)
            keep_prob = NOISE_MIN_PROB + (noise_val + 1) * 0.5 * (1 - NOISE_MIN_PROB)
            if random.random() > keep_prob:
                continue

            pts.append((x, y))
            if len(pts) == n:
                break
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

# NEW – add slight random jitter to seeds to break perfect regularity ----------------

def jitter_points(pts: np.ndarray, mask: Polygon | MultiPolygon, std: float) -> np.ndarray:
    """Return *pts* with Gaussian noise added (σ = *std*°).

    Each moved point is kept only if it remains inside *mask*; otherwise the
    original position is retained. This maintains roughly equal seed counts
    while introducing size variation between the resulting Voronoi cells."""
    jittered: list[tuple[float, float]] = []
    for x, y in pts:
        for _ in range(10):  # try a few times to keep point on land
            dx, dy = np.random.normal(scale=std, size=2)
            nx, ny = x + dx, y + dy
            if mask.contains(Point(nx, ny)):
                jittered.append((nx, ny))
                break
        else:
            # could not find on–land jitter – keep original
            jittered.append((x, y))
    return np.array(jittered)

seeds = jitter_points(seeds, land, JITTER_STD)

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
    h=hashlib.sha256(str(key(a,b)).encode()).digest()
    rng=random.Random(int.from_bytes(h[:4],"little"))

    if WARP_STYLE == "sine":
        cycles=rng.randint(MIN_WIGGLES,MAX_WIGGLES)
        n=max(2,int(L/MAX_SEG_LEN*cycles))
        phase=rng.uniform(0,2*math.pi)
        amp=rng.uniform(0.6,1.0)*MAX_OFFSET
        ux,uy=(-dy/L,dx/L) if L else (0,0); pts=[a]
        for i in range(1,n):
            t=i/n
            sine=math.sin(phase+t*math.pi*cycles)
            base_x=a[0]+dx*t
            base_y=a[1]+dy*t
            offset=amp*sine
            px=base_x+ux*offset
            py=base_y+uy*offset
            if not land.contains(Point(px,py)):
                # shrink offset until point lies on land (or give up)
                tmp=offset
                while abs(tmp) > 1e-6:
                    tmp*=0.5
                    px=base_x+ux*tmp; py=base_y+uy*tmp
                    if land.contains(Point(px,py)):
                        break
                else:
                    px,py=base_x,base_y  # fallback straight line point
            pts.append((px,py))
        pts.append(b)
        return pts

    # --- angular style -------------------------------------------------
    n=max(3,int(L/MAX_SEG_LEN))  # number of vertices along edge
    ux,uy=(-dy/L,dx/L) if L else (0,0)
    pts=[a]
    for i in range(1,n):
        t=i/n
        taper=math.sin(math.pi*t)          # 0 at ends, 1 at middle
        amp_off=rng.uniform(0.4,1.0)*MAX_OFFSET*taper
        sign=rng.choice((-1,1))
        base_x=a[0]+dx*t
        base_y=a[1]+dy*t
        px=base_x+ux*amp_off*sign
        py=base_y+uy*amp_off*sign
        if not land.contains(Point(px,py)):
            tmp=amp_off
            while abs(tmp) > 1e-6:
                tmp*=0.5
                px=base_x+ux*tmp*sign; py=base_y+uy*tmp*sign
                if land.contains(Point(px,py)):
                    break
            else:
                px,py=base_x,base_y
        pts.append((px,py))
    pts.append(b)
    return pts

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

# compute average area once (for information only) and use absolute threshold
global_avg_area = sum(c.area for c in cells) / len(cells)
threshold = SLIVER_MIN_AREA
print(f"Global average cell area: {global_avg_area:.4f} → sliver threshold (fixed): {threshold:.4f}")

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
    # decide which cells are considered *small* using fixed threshold
    small_cells = [idx for idx, c in enumerate(cells) if c.area < threshold]

    if not small_cells:
        print("No cells below sliver threshold – merge complete.")
        break

    print(f"Small cells this pass: {len(small_cells)} (threshold {threshold:.4f})")

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

# 6A. deterministic neighbour merges for shape variety ------------------
print("Deterministic neighbour merging …")
removed_rand: set[int] = set()
processed: set[int] = set()   # pickers we have already handled

# build processing order: smallest first, then others shuffled
small_indices = sorted([i for i, c in enumerate(cells) if c.area <= SMALL_FIRST_AREA], key=lambda i: cells[i].area)
other_indices = [i for i in range(len(cells)) if i not in small_indices]
random.shuffle(other_indices)
order = small_indices + other_indices

for idx in order:
    if idx in removed_rand or idx in processed:
        continue  # already merged or already acted as picker

    # skip if this polygon is larger than skip limit -> acts as anchor absorbing others but not absorbed
    if cells[idx].area > MERGE_SKIP_LARGE_AREA:
        continue

    # gather neighbours eligible to merge
    neighbours = [j for j in range(len(cells))
                  if j != idx and j not in removed_rand and j not in processed
                  and cells[idx].touches(cells[j]) and cells[j].area <= MERGE_SKIP_LARGE_AREA]

    if not neighbours:
        continue

    target_merge_count = random.randint(MERGE_NEIGH_MIN, MERGE_NEIGH_MAX)
    random.shuffle(neighbours)

    merged_count = 0
    current_area = cells[idx].area

    for j in neighbours:
        if merged_count >= target_merge_count:
            break
        if j in removed_rand:
            continue

        # skip neighbour if combined area would exceed limit
        if current_area + cells[j].area > MERGE_SKIP_LARGE_AREA:
            continue

        try:
            if safe_union(j, idx):
                removed_rand.add(j)
                merged_count += 1
                current_area = cells[idx].area  # update area after union
        except (TopologicalError, GEOSException):
            continue

    if merged_count:
        print(f"Merged {merged_count} neighbours into cell {idx} (area now {current_area:.2f} deg²)")

    # mark this cell as processed so it won't be absorbed later
    processed.add(idx)

if removed_rand:
    cells = [c for k, c in enumerate(cells) if k not in removed_rand]
    print(f"Deterministic merge: merged {len(removed_rand)} cells → {len(cells)} remaining.")

final = cells
print(f"\nFinished merging. Started with {initial_polygon_count}, ended with {len(final)} polygons.")

# 7. save ---------------------------------------------------------------
json_path=pathlib.Path(OUT)
json_path.write_text(json.dumps({
    "type":"FeatureCollection",
    "features":[{"type":"Feature","geometry":mapping(g),"properties":{}} for g in final]
}))
print(f"Saved {len(final)} polygons after merge → {OUT}")

# optional corner smoothing ------------------------------------------------

def chaikin(coords: list[tuple[float,float]]) -> list[tuple[float,float]]:
    """One iteration of Chaikin's corner-cutting algorithm on a closed ring."""
    new=[]
    for i in range(len(coords)-1):
        p=coords[i]; q=coords[i+1]
        new.append((0.75*p[0]+0.25*q[0], 0.75*p[1]+0.25*q[1]))
        new.append((0.25*p[0]+0.75*q[0], 0.25*p[1]+0.75*q[1]))
    new.append(new[0])
    return new

def smooth(poly, iters:int=1):
    if iters<=0:
        return poly
    if poly.geom_type=="Polygon":
        ring=list(poly.exterior.coords)
        for _ in range(iters):
            ring=chaikin(ring)
        return Polygon(ring)
    else:
        return MultiPolygon([smooth(p,iters) for p in poly.geoms])

if SMOOTH_ITER>0:
    print(f"Smoothing polygons with Chaikin (iterations={SMOOTH_ITER}) …")
    cells=[smooth(c, SMOOTH_ITER) for c in cells]

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