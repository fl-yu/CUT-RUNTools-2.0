## Pre-requisite software 
A few software need to be installed before running the single-cell module on the basis of [pre-requisites of bulk data analysis](). We encourage users to use conda to install and manage software.

Install GNU parallel (a software program developed for running several jobs in parallel to make full use of the computational resource (more information can be found [here](https://www.gnu.org/software/parallel))  

```
conda install -c conda-forge parallel
```

Install deeptools; get the path using the command: *which deeptools*

```
conda install -c bioconda deeptools
```

Install tabix (it is the a generic tool that indexes position sorted files in TAB-delimited formats); get the path using the command: *which tabix*

```
conda install -c bioconda tabix
```
    
Three python modules (umap, leidenalg, and igraph) are required to be installed on your system  
To install umap 
    
```
$pythonbin/pip install umap-learn leidenalg igraph
```
    
To test if the package is installed successfully, for example:
    
```
$pythonbin/pip list | grep umap
```

    
Additionally, R packages including reticulate("1.15"), leiden("0.3.3"), data.table("1.11.6"), Matrix("1.2.18"), irlba("2.3.3"), Rtsne("0.15"), RANN("2.6.1"), igraph("1.2.4.2"), uwot("0.1.8"), rGREAT("1.14.0"), ggplot2("3.3.0") need to be installed in the library of R you specified in the JSON configuration file. 



