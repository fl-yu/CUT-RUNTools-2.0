# CUT&RUNTools 2.0

#### *A tool for bulk and single-cell CUT&RUN/CUT&Tag data processing and analysis*

## Overview

CUT&RUNTools 2.0 is a major update of CUT&RUNTools [(link)](https://bitbucket.org/qzhudfci/cutruntools/), including a set of new features specially designed for CUT&RUN and CUT&Tag experiments. Both of the bulk and single-cell data can be processed, analyzed and interpreted using CUT&RUNTools 2.0.

## Main Features

*Single-cell data analysis*

- raw data processing in parallel
- quality assessment and visualization
- feature-by-cell matrix construction  
- dimensionality reduction and clustering analysis  
- cell-type-specific pseudo-bulk data analysis  
 

*Bulk data analysis update* 

- supporting spike-in sequence alignment and data normalization
- options of experiment for CUT&RUN or CUT&Tag data analysis  
- flexiable option for fragments selections (>120bp) 
- supporting different peak calling strategies 
- new functions for peaks annotation 
- compatiable with more computational platforms  


## Installation and Requirement
Please refer to [Installation Instructions](docs/INSTALL.md).

## Bulk data analysis

To learn new features for CUT&RUNTools 2.0 on the bulk data analysis, please see [New features](docs/bulk-news.md).


Once the package is installed, see [Data Example](docs/bulk-QUICK.md) for examples.

See [Usage Page](docs/bulk-USAGE.md) about how to set up the pipeline and detailed usage. 

See [Directory Structure](docs/bulk-DIRECTORY.md) for basic structure of input and output folders.  


--------


## Single-cell data analysis

To get started, please see [Overview](docs/sc-OVERVIEW.md).  

See [Data Example](docs/sc-QUICK.md) for a most excellent starting point to get familiar with the tool.  

See [Usage Page](docs/sc-USAGE.md) about how to set up the pipeline and detailed usage.  

See [Directory Structure](docs/sc-DIRECTORY.md) for basic structure of input and output folders.

<div align=center> <img src="images/scCRtools.png" width="680" height="318"> </div> 

<p align="center">The schematic view of the single-cell data processing module of CUT&RUNTools</p>  


If you run into issues and would like to report them, you can use the "Issues" tab on the left hand side.  
Alternatively, you can contact authors: fulong.yu{at}childrens.harvard.edu, or gcyuan{at}jimmy.harvard.edu.  
## Read our paper:  
*CUT&RUNTools 2.0: A pipeline for single-cell and bulk-level CUT&RUN and CUT&Tag data analysis*; bioRxiv 615179; doi: https://doi.org/10.1101/615179  
[*CUT&RUNTools: a flexible pipeline for CUT&RUN processing and footprint analysis*](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4); Genome Biol 20, 192; 2019

