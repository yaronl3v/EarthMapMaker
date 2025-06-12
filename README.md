# EarthMapMaker

A Python toolkit for creating detailed, coastline-aware fantasy world maps. The primary script, **`generate_map_fantasy_optimized.py`**, carves the Earth's landmass into thousands of visually pleasing Voronoi regions and outputs GeoJSON.

## What it does

The main generator (`generate_map_fantasy_optimized.py`) samples up to 60 000 seed points over land, applies multi-pass merging, edge warping, and rigorous sliver cleanup to build a polygon mosaic that:

- Respects real coastlines – no cell spills into the sea
- Retains geographic quirks by dropping tiny islands below a configurable size
- Produces natural-looking, slightly irregular borders with optional curvature
- Removes micro-slivers and isolated specks for a tidy final layer
- Writes standard GeoJSON ready for any GIS or web-mapping stack

## Features (fantasy generator)

- **High seed density**: 60 k seeds → fine-grained regional detail
- **Curvy or straight edges**: Toggle border warping via `USE_CURVES`
- **Smart sliver merging**: Two-phase merge with size thresholds and spatial indexes
- **Funny group merging**: Randomised neighbour unions create organic region shapes
- **Efficient spatial indexing**: STRtree accelerates neighbour look-ups
- **Isolated-island pruning**: Final pass drops tiny standalone polygons

### Alternative generators

The project also ships two simpler generators for specific needs:

| Script | Purpose |
| ------ | ------- |
| `generate_map_eqaul_shapes_round.py` | Creates a perfectly uniform equal-area grid – ideal for statistical aggregation |
| `generate_map_equal_merge.py` | Builds variable-size cells with deterministic neighbour merges and fewer curves |

## Setup

1.  **Create a virtual environment:**
    ```bash
    python -m venv venv
    ```

2.  **Activate the virtual environment:**
    -   On Windows:
        ```bash
        .\\venv\\Scripts\\activate
        ```
    -   On macOS and Linux:
        ```bash
        source venv/bin/activate
        ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Running the project

After installing dependencies, pick a generator script:

- **Fantasy optimized map** (recommended)
  ```bash
  python generate_map_fantasy_optimized.py
  ```
- **Uniform equal-area grid**
  ```bash
  python generate_map_eqaul_shapes_round.py
  ```
- **Smart merged polygons**
  ```bash
  python generate_map_equal_merge.py
  ```

Each script writes a `world_cells.geojson` file in the project root.

## Output

The script produces a GeoJSON file with:
- Equal-area polygons covering the globe
- Consistent geometric shapes
- Standard GeoJSON format for easy integration with mapping tools
- Ready to use for data visualization and geographic analysis 