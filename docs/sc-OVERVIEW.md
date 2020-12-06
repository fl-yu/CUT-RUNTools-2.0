## The overview of single-cell module of CUT&RUNTools 2.0

<div align=center> <img src="../images/scCRtools.png" width="680" height="318"> </div> 

<p align="center">The schematic view of the single-cell data processing module of CUT&RUNTools</p>

The single-cell module of CUT&RUNTools generally includes four main processing steps for single-cell data analysis: 
-  (i) raw data processing and quality assessment
-  (ii) feature-by-cell matrix construction
-  (iii) dimensionality reduction and clustering analysis
-  (iv) cell-type-specific pseudo-bulk data analysis

The raw data of individual cells are trimmed and then aligned to the reference genome. Peak calling is performed for the data of cell aggregation. A comprehensive QC report and several diagnostic plots are generated to help to evaluate experiment quality. Three options of input features (peak-by-cell, bin-by-cell and customFeature-by-cell) are available for count matrix construction. The cell annotation from the clustering analysis will be generated. The cells from the same cluster are merged into a pseudo-bulk profile and the corresponding genome track files for both individual cells and the pooled signal are automatically generated. The bulk-level analysis including direct binding and TF footprinting analysis provided by the original version of CUT&RUNTools, and several newly developed downstream functional analysis such as Gene Ontology (GO) enrichment and genomic elements annotation analysis, (co-)regulatory TF discovery are also available to facilitate the interpretation of regulatory mechanism for different cell populations.

The main steps in the raw data processing and feature-by-cell matrix construction parts are performed in parallel to make full use of the available computational resource. For each run, only a configuration file specifying the options and parameters is required to be provided for the analysis. Users are able to easily run the entire workflow or an individual single step.

See the following links for more details:

- [Installation Instructions](./sc-INSTALL.md)
- [Quick Start](./sc-QUICK.md)
- [Usage Page](./sc-USAGE.md)
- [Directory Structure](./sc-DIRECTORY.md)
