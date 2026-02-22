# create_WLC_indices
Python tools and tutorial for building raster-based suitability indices using Weighted Linear Combination (WLC) / Weighted Overlay Analysis (WOA).

This repository provides a reproducible workflow for combining multiple geospatial raster layers into composite indices commonly used in environmental, hydrological, and land-suitability studies.

## Overview 

Many spatial decision-making problems require integrating multiple environmental variables — such as soils, topography, vegetation, permeability, or accessibility — into a single index representing suitability, vulnerability, or potential.

This project implements that workflow using Python and GeoTIFF rasters.

The core script:

`SuitabilityMappingTools.py`

provides functions to:

    - Inspect raster datasets

    - Ensure spatial compatibility between layers

    - Assign weights to variables

    - Generate weighted composite indices

    - Visualize results

    - Produce basic diagnostic outputs

The tutorial notebook walks through the full process from raw input rasters to final suitability maps.

> This repository contains the index creation workflow only.

It does NOT include the full sensitivity analysis framework used in the revised version of the associated research paper resubmitted to Journal of Hydrology Regional Studies in 2026, which includes more advanced methods such as:

    - Probabilistic weight sampling

    - Monte Carlo approaches

    - Formal uncertainty quantification

    - Robustness analysis across parameter spaces

    - Scenario exploration beyond simple permutations or removal tests

The tools included here are intended for:

✔ Reproducible index generation
✔ Method transparency
✔ Educational purposes
✔ Rapid prototyping
✔ Demonstrating WLC methodology


What is Weighted Linear Combination (WLC)?

WLC combines standardized raster layers into a composite score:

$$
\text{Index} = \sum_{i=1}^{n} w_i \, x_i
$$

$$
\text{where:}
\quad
S = \text{composite suitability score (index)}, \\
w_i = \text{weight assigned to criterion } i, \\
x_i = \text{standardized value of criterion } i, \\
n = \text{number of criteria}, \\
\sum_{i=1}^{n} w_i = 1
$$


## Tutorial: What It Walks You Through

The tutorial notebook demonstrates a complete end-to-end workflow for creating raster-based suitability indices using a Weighted Linear Combination (WLC). It guides the user from raw input layers through preprocessing, weighting, index generation, and interpretation of results.

---

### 1) Preparing Input Raster Layers

The workflow begins by loading and inspecting GeoTIFF datasets representing the criteria to be combined (e.g., soils, slope, permeability, vegetation).

This step includes:

- Loading raster layers from a directory
- Inspecting metadata (projection, extent, resolution)
- Checking value ranges
- Identifying missing or NoData values
- Visualizing each layer

Functions demonstrated:

- `plot_layers()`
- `get_layer_info()`
- `get_dir_info()`

---

### 2) Ensuring Spatial Compatibility

All rasters used in a WLC must share identical spatial properties. The tutorial demonstrates how to ensure compatibility by resampling layers to match a reference raster.

Key requirements:

- Same coordinate reference system (CRS)
- Same cell size
- Same grid alignment
- Same spatial extent

Function demonstrated:

- `resample_layers()`

---

### 3) Selecting Variables

Users choose which raster layers will be included in the index calculation.

The selection tool:

- Lists available `.tif` files
- Allows interactive selection by index
- Retrieves metadata for selected layers
- Assigns equal weights initially

Function demonstrated:

- `select_tif_files()`

---

### 4) Defining Weights

Weights represent the relative importance of each variable in the composite index.

The tutorial shows how to:

- Use equal weighting
- Specify custom weights
- Verify that weights sum to one
- Inspect assigned weights

Function demonstrated:

- `define_weights()`

---

### 5) Creating the Weighted Index

This is the core step of the workflow.

The script:

- Reads each input raster
- Applies weights to each layer
- Handles NoData values appropriately
- Computes the weighted sum
- Writes the output as a new GeoTIFF
- Saves metadata documenting inputs and weights
- Optionally generates histograms and summary statistics

Function demonstrated:

- `create_weighted_index()`

Outputs include:

- Composite suitability raster
- Metadata file describing inputs and weights
- Histogram plots of index values
- Text summaries of distribution statistics

---

### 6) Visualizing and Interpreting Results

The tutorial demonstrates how to explore the resulting index:

- Viewing output rasters
- Understanding spatial patterns
- Interpreting index values
- Examining the distribution of suitability scores

Histogram outputs include:

- Pixel counts
- Percent area per value class
- Mean values
- Area estimates

---

## Optional: Basic Sensitivity Tests Included

Although this repository does not contain the full sensitivity analysis framework from the revised paper, it includes simple exploratory tests to assess how results change under alternative assumptions.

---

### Weight Permutation Test

This test swaps weights among variables while keeping the same values, allowing users to examine how sensitive results are to weight assignment.

Function:

- `sensitivity_analysis_weights()`

Outputs:

- Multiple index rasters for each permutation
- Summary statistics for each scenario
- CSV tables of results

---

### Variable Removal Test (One-at-a-Time)

Each variable is removed individually and the index recalculated using the remaining variables with equal weights. This helps identify influential variables.

Function:

- `sensitivity_analysis_removal()`

Outputs:

- Scenario rasters for each removal
- Control scenario using all variables
- Comparative statistics

---

### Summary Products for Multiple Scenarios

When multiple index maps are generated, the script can produce summary rasters that characterize variability across scenarios.

Generated products include:

- Mean index raster
- Standard deviation raster
- Range raster
- Coefficient of variation raster
- Markdown summaries of statistics
- CSV tables of results

Functions:

- `generate_summary_rasters()`
- `generate_summary_markdown()`
- `generate_markdown_from_csv()`

---

## Typical Use Cases

This workflow is suitable for studies involving:

- Groundwater recharge potential mapping
- Habitat suitability analysis
- Environmental vulnerability assessments
- Land capability evaluation
- Infrastructure siting
- Resource prioritization
- Spatial decision support systems
- Multi-criteria decision analysis (MCDA)

---

## Requirements

Python packages required:
    numpy
    pandas
    matplotlib
    rasterio
    glob
    itertools


Install via pip:

```
pip install rasterio numpy pandas matplotlib
```



---

## Input Data Requirements

Input rasters should:

- Represent meaningful criteria for the study objective
- Be standardized or scaled to comparable ranges
- Use consistent units where applicable
- Be spatially aligned and compatible
- Have clearly defined NoData values

---

## Limitations

This implementation assumes:

- Linear trade-offs between criteria
- Static weights across space
- Deterministic inputs
- No explicit uncertainty propagation

It does not implement:

- Nonlinear aggregation methods
- Fuzzy logic operators
- Bayesian approaches
- Data-driven weight calibration
- Machine learning models

---

## Intended Audience

This repository is designed for:

- Researchers learning WLC and MCDA methods
- Students in GIS, hydrology, or environmental science
- Practitioners seeking transparent workflows
- Analysts prototyping suitability indices
- Users transitioning from GUI-based GIS tools to Python

---

## Reproducibility

All outputs are generated programmatically and include metadata documenting:

- Input layers used
- Weights applied
- Processing date

This supports transparent and reproducible workflows.

---

## Future Extensions (Not Included Here)

The full research framework extends beyond this repository and may include:

- Formal uncertainty quantification
- Advanced sensitivity analysis methods
- Scenario ensembles
- Hierarchical index construction
- Automated processing pipelines
- Integration with decision models

---

## Citation

If using this workflow in academic work, cite this repository

