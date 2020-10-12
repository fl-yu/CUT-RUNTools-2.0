## Pre-requisite software 
A few software need to be installed before running the single-cell module on the basis of [pre-requisites of bulk data analysis](./bulk-INSTALL.md). We encourage users to use conda to install and manage software.

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

## Configuration file

Similar to the bulk data analysis, CUT&RUNTools needs a JSON configuration file (named sc-config.json) to be written.

<!-- end list -->

``` json
"software_config" {
    	"Rscriptbin": "/n/app/R/3.3.3/bin",
    	"pythonbin": "/n/app/python/2.7.12/bin",
    	"perlbin": "/n/app/perl/5.24.0/bin",
    	"javabin": "/n/app/java/jdk-1.8u112/bin",
    	"trimmomaticbin": "/n/app/trimmomatic/0.36/bin",
    	"trimmomaticjarfile": "trimmomatic-0.36.jar",
   	"bowtie2bin": "/n/app/bowtie2/2.2.9/bin",
    	"samtoolsbin": "/n/app/samtools/1.3.1/bin",
    	"adapterpath": "/home/qz64/cutrun_pipeline/adapters", 
    	"picardbin": "/n/app/picard/2.8.0/bin",
    	"picardjarfile": "picard-2.8.0.jar",
    	"macs2bin": "/n/app/macs2/2.1.1.20160309/bin",
    	"macs2pythonlib": "/n/app/macs2/2.1.1.20160309/lib/python2.7/site-packages",
    	"kseqbin": "/home/qz64/cutrun_pipeline", 
    	"memebin": "/n/app/meme/4.12.0/bin", 
    	"bedopsbin": "/n/app/bedops/2.4.30", 
    	"bedtoolsbin": "/n/app/bedtools/2.26.0/bin",
    	"makecutmatrixbin": "/home/qz64/.local/bin",
    	"bt2idx": "/n/groups/shared_databases/bowtie2_indexes",
    	"genome_sequence": "/home/qz64/chrom.hg19/hg19.fa",
    	"extratoolsbin": "/home/qz64/cutrun_pipeline", 
    	"path_parallel": "/path/to/parallel", 
	"path_deeptools": "/path/to/deeptools",
	"path_tabix": "/path/to/tabix", 
    }
"sc_parameters": {
	"single_cell": "TRUE", 
	"fastq_directory": "/path/to/fastq", 
	"workdir": "/path/to/workdir", 
	"genome": "hg38", 
	"chrome_sizes_file": "/path/to/hg38.chrom.sizes",
	"cores": "8", 
	"percentage_rip": "30", 
	"num_reads_threshold": "10000", 
	"peak_caller": "macs2", 	
	"peak_type": "narrow", 
	"matrix_type": "peak_by_cell", 
	"bin_size": "5000", 
	"feature_file": "/path/to/feature_file", 
	"experiment_type": "CUT&Tag", 
	"cluster_resolution": "0.8", 
	"cluster_pc": "30", 
	"experiment_name": "scCUT&Tag", 
    }
"run_pipeline": {
	"entire_pipeline": "TRUE", 
	"individual_step": "NULL", 
	"step2_bamfile_dir": "$workdir/sc_aligned.aug10/dup.marked.clean", 
	"step2_output_dir": "$workdir/sc_countMatrix", 
	"step2_qc_pass_file": "$workdir/sc_qc/report/statistics_QCpassed.txt", 
	"step3_count_matrix": "$workdir/sc_countMatrix/feature-by-cell_matrix.txt", 
	"step3_output_dir": "$workdir/sc_cluster", 
	"step4_bamfile_dir": "$workdir/sc_aligned.aug10/dup.marked.clean", 
	"step4_cell_anno_file": "$workdir/sc_cluster/leiden_cluster_annotation.txt", 
	"step4_output_dir": "$workdir/sc_pseudoBulk", 
    }
      
```

The **software_config** section (the first 24 lines) concerns the software installation.  
Similar to the configure.json of the bulk data processing, all the requirements of software can be defined here. Three more paths (*path_parallel*, *path_deeptools* and *path_tabix*) should be specified compared to the bulk configure.json file. The rest is related to an actual analysis (explained in [Usage Page](./sc-USAGE.md)). 

## Validate prerequisites

To check if the paths are correct and if the softwares in these paths indeed work:

``` shell
./validate.py sc-config.json --ignore-input-output --software

```
This script checks that your configuration file is correct and all paths are correct. You will get an empty line if the validate.py script runs without errors.
***

See the following links for more details:

- [Overview](./sc-OVERVIEW.md)
- [Quick Start](./sc-QUICK.md)
- [Usage Page](./sc-USAGE.md)
- [Directory Structure](./sc-DIRECTORY.md)
- [The full tutorial of CUT&RUNTools 2.0](./2.0-TUTORIAL.md)

