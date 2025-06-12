# EarthMapMaker

A Python tool that generates Earth maps with equal-sized polygons, creating clean geometric shapes without weird distortions and returning GeoJSON output.

## What it does

This project creates a world map divided into equal-sized polygons that cover the entire Earth's surface. Unlike traditional map projections that can create distorted shapes, this tool generates uniform geometric cells that maintain consistent size and shape across the globe. The output is a clean GeoJSON file that can be used for data visualization, geographic analysis, or any application requiring standardized geographic regions.

## Features

- **Equal-sized polygons**: All map cells have the same area, ensuring fair geographic representation
- **No weird shapes**: Clean, geometric polygons without projection distortions
- **GeoJSON output**: Standard format compatible with most mapping libraries and GIS tools
- **Full Earth coverage**: Complete coverage of the planet's surface
- **Easy to use**: Simple Python script with minimal dependencies

### Advanced generator (`generate_map_equal_merge.py`)

In addition to the equal-area generator above, the project includes an alternative script that focuses on visual variety and performance:

- **Smart polygon merging**: Deterministic neighbour merges create natural-looking regions without manual tuning
- **Efficient spatial indexing**: Uses STRtree indexes for lightning-fast neighbour and touch queries even with thousands of cells
- **Accurate coastline clipping**: All borders are clipped to the real coastline, guaranteeing that no geometry spills into the sea
- **Zero sliver leftovers**: Aggressive sliver removal passes ensure tiny polygons are merged away, leaving a clean dataset

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

After setting up the environment and installing the dependencies, choose one of the generator scripts:

- **Equal-sized polygons**
  ```bash
  python generate_map_eqaul_shapes_round.py
  ```
- **Smart merged polygons**
  ```bash
  python generate_map_equal_merge.py
  ```

Each script writes a `world_cells.geojson` file in the project root.

- The **equal-sized** mode produces a perfectly uniform grid ideal for statistical analysis.
- The **smart merged** mode produces variable-size but aesthetically pleasing polygons that hug coastlines closely and contain no micro-slivers.

## Output

The script produces a GeoJSON file with:
- Equal-area polygons covering the globe
- Consistent geometric shapes
- Standard GeoJSON format for easy integration with mapping tools
- Ready to use for data visualization and geographic analysis 