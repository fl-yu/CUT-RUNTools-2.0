## Pre-requisite software 
A few software (list in the following) need to be installed before running the single-cell module on the basis of [pre-requisites of bulk data analysis](./bulk-INSTALL.md) (including Trimmomatic, Bowtie2, Samtools, Picard, MACS2, MEME, Bedops, Bedtools, deeptools and Download genome assemblies). We encourage users to use conda to install and manage software.

Install GNU parallel (a software program developed for running several jobs in parallel to make full use of the computational resource (more information can be found [here](https://www.gnu.org/software/parallel))  

```
conda install -c conda-forge parallel
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

```
{
	"software_config": {
		"Rscriptbin": "/path/to/R/bin",
		"pythonbin": "/path/to/python/bin",
		"perlbin": "/path/to/perl/bin",
		"javabin": "/path/to/java/bin",
		"trimmomaticbin": "/path/to/trimmomatic/bin",
		"trimmomaticjarfile": "trimmomatic-0.36.jar",
		"bowtie2bin": "/path/to/bowtie2/bin",
		"samtoolsbin": "/path/to/samtools/bin",
		"adapterpath": "/path/to/cutrun_pipeline/adapters", 
		"picardbin": "/path/to/picard/bin",
		"picardjarfile": "picard-2.8.0.jar",
		"macs2bin": "/path/to/macs2/bin",
		"macs2pythonlib": "/path/to/macs2/2.1.1.20160309/lib/python2.7/site-packages",
		"kseqbin": "/path/to/cutrun_pipeline", 
		"memebin": "/path/to/meme/bin", 
		"bedopsbin": "/path/to/bedops/bin", 
		"bedtoolsbin": "/path/to/bedtools/bin",
		"makecutmatrixbin": "/home/user/.local/bin",
		"bt2idx": "/path/to/bowtie2_indexes",
		"genome_sequence": "/path/to/chrom.hg19/hg19.fa",
		"extratoolsbin": "/home/user/cutrun_pipeline", 
		"path_parallel": "/path/to/parallel/bin", 
		"path_deeptools": "/path/to/deeptools/bin",
		"path_tabix": "/path/to/tabix/bin", 
    },
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
    },
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


