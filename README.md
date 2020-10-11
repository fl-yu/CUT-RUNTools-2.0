# CUT&RUNTools 2.0

#### *A tool for bulk and single-cell CUT&RUN/CUT&Tag data processing and analysis*

## Overview

CUT&RUNTools 2.0 is a major update of CUT&RUNTools [(link)](https://bitbucket.org/qzhudfci/cutruntools/), including a set of new features specially designed for CUT&RUN and CUT&Tag experiments. In this release, both of the bulk and single-cell data can be processed, analyzed and interpreted using CUT&RUNTools 2.0.

## New Features

*Bulk data analysis* 

- new settings for CUT&Tag data analysis  
- spike-in alignment  
- different normalization method  
- support for analysis involving larger fragments (>120bp)  
- new functions of peak calling and annotation   
 
*Single-cell data analysis*

- raw data processing and quality assessment  
- feature-by-cell matrix construction  
- dimensionality reduction and clustering analysis  
- cell-type-specific pseudo-bulk data analysis  


## Bulk data analysis

To get started, please see [Installation Instructions](docs/bulk-INSTALL.md). 

Once the package is installed, see [Quick Start](docs/bulk-QUICK.md) for examples.

See [Usage Page](docs/bulk-USAGE.md) about how to set up the pipeline and detailed usage. 

See [Directory Structure](docs/bulk-DIRECTORY.md) for basic structure of input and output folders.  

Comparing to the last version for bulk data analysis, the requirement of SLURM system was removed to compatiable with more computational platforms. 

## Single-cell data analysis

To get started, please see [Overview](docs/sc-OVERVIEW.md) and [Installation Instructions](docs/sc-INSTALL.md).  

See [Quick Start](docs/sc-QUICK.md) for a most excellent starting point to get familiar with the tool.  

See [Usage Page](docs/sc-USAGE.md) about how to set up the pipeline and detailed usage.  

See [Directory Structure](docs/sc-DIRECTORY.md) for basic structure of input and output folders.

### The full tutorial of CUT&RUNTools 2.0 can be found [here](docs/2.0-TUTORIAL.md).

--------

If you run into issues and would like to report them, you can use the "Issues" tab on the left hand side.  
Alternatively, you can contact authors: fulong.yu{at}childrens.harvard.edu, or gcyuan{at}jimmy.harvard.edu.  
Read our [paper](). Sign up for CUT&RUNTools Google group mailing list to ask questions, receive updates.




