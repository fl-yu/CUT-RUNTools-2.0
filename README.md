## To do
1.bulk function code and tutorial  
2.install tutorial (if there is any script for software install) update and sc usage  
3.bulk code test  
4.tutorial update for integration  
5.News: what is the update of this version? New update of bulk normalization and spike-in processing



## The single-cell CUT&RUNTools module
**A new module of CUTRUNTools is provided for single-cell epigenome profiling data (scCUT&RUN and scCUT&Tag etc.) processing and analysis.** 

<div align=center> <img src="./images/scCRtools.png" width="680" height="318"> </div> 

<p align="center">The schematic view of the single-cell data processing module of CUT&RUNTools</p>


The single-cell module of CUT&RUNTools generally includes four main processing steps for single-cell data analysis: (i) raw data processing and quality assessment, (ii) feature-by-cell matrix construction, (iii) dimensionality reduction and clustering analysis, (iv) cell-type-specific pseudo-bulk data analysis. 

The raw data of individual cells are trimmed and then aligned to the reference genome. Peak calling is performed for the data of cell aggregation. A comprehensive QC report and several diagnostic plots are generated to help to evaluate experiment quality. Three options of input features (peak-by-cell, bin-by-cell and customFeature-by-cell) are available for count matrix construction. The cell annotation from the clustering analysis will be generated. The cells from the same cluster are merged into a pseudo-bulk profile and the corresponding genome track files for both individual cells and the pooled signal are automatically generated. The bulk-level analysis including direct binding and TF footprinting analysis provided by the original version of CUT&RUNTools, and several newly developed downstream functional analysis such as Gene Ontology (GO) enrichment and genomic elements annotation analysis, (co-)regulatory TF discovery are also available to facilitate the interpretation of regulatory mechanism for different cell populations.  
The main steps in the raw data processing and feature-by-cell matrix construction parts are performed in parallel to make full use of the available computational resource. For each run, only a configuration file specifying the options and parameters is required to be provided for the analysis. Users are able to easily run the entire workflow or an individual single step.



## Pre-requisite software in addition to bulk-data analysis
A few software need to be installed before running the single-cell module on the basis of pre-requisites of bulk data analysis. We encourage users to use conda to install and manage software.

Install GNU parallel (a software program developed for running several jobs in parallel to make full use of the computational resource, more information: [https://www.gnu.org/software/parallel] (https://www.gnu.org/software/parallel))  

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



## Usage

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

- The **software_config** section defines the software paths for various sequencing tools.  
Similar to the configure.json of the bulk data processing, all the requirements of software can be defined here. Three more paths (*path_parallel*, *path_deeptools* and *path_tabix*) should be specified compared to the bulk configure.json file.

- The **sc_parameters** section contains the general parameters for the data processing.

``` 
PARAMETERS
----------

	[single_cell]: To perform the single-cell data processing 
			options: TRUE (default), FALSE
						
	[fastq_directory]: A folder containing all the fastq files of barcode cells, pattern with ***_1.fastq.gz
		and ***_2.fastq.gz 
			options: NO default 
			 			
	[workdir]: A folder containg all the output of data processing and analysis
			options: NO default 
	
	[genome]: Reference genome version 
			options: hg38 (default), hg19, mm10 or mm9
			 			
	[chrome_sizes_file]: A file to define the chromosome lengths for a given genome
			options: NO default 
			 			
	[cores]: Number of processors to use
			options: 8 (default) 
			 			
	[percentage_rip]: A criterion to filter barcode cells based on the signal-to-noise (percentage of 
	      reads in peaks). The barcode cells will be filtered out if they having reads percentage in the
	      aggregation peaks less than the threshold
			options: 30 (default) or a "numeric" between 0 and 100. Set to 0 if you do not 
				 want to filter the cells based on this measurements
			 			
	[num_reads_threshold]: A criterion to filter barcode cells. The barcode cells will be filtered out 
	      if they having the number of properly paired reads less than the threshold. Set to 0 if you do
	      not want to filter the cells based on this measurements  
			options: 10000 (default) 
			 			
	[peak_caller]: To specify statistical significant peaks called from macs2 or SEACR method used for 
	      QC (i.e. cell filtering), countmatrix construction and differential peaks analysis 
			options: macs2 (default) or SEACR
			 			
	[peak_type]: If peak_caller is specified as macs2, peak_type will be activated. It should be chosen 
	      closely associated with the experimental types used. The narrow peak mode is suitable for a  
	      typical TF, while the broad peak mode is good to use for the chromatin factors generally having 
	      a large binding domain or epigenetic modification such as H3K27me3 or H3K9me3 having dispersed 
	      localization 
			options: narrow (default) or broad 
	
	[matrix_type]: This parameter control the feature used for count matrix generation (see Tutorial for 
	      further information) 
			options: bin_by_cell (default), peak_by_cell or customFeature_by_cell 
			 			
	[bin_size]: If matrix_type is specified as bin_by_cell, the feature file of genome-wide bins with 
	      resolution of bin_size is automatically generated and the feature-by-cell matrix will be 
	      constructed
			options: 5000 (default) 
			 			
	[feature_file]: If matrix_type is specified as customFeature_by_cell, user can use a customized  
	      feature file to construct the feature-by-cell matrix. This is useful when the peaks derived  
	      from bulk data or other informative data are desired as features to analyze the single-cell data.  
			options: NO default 
			 			
	[experiment_type]: Technology used for the single-cell experiment. It is important for filtering   
	      the peak, accurate cut point matrix generation and, motif and footprinting analysis 
			options: CUT&Tag (default) or CUT&RUN 
			 			
	
	[cluster_resolution]: The resolution parameter in the Leiden clustering algorithm to determine the 
	      granularity level at which communities are detected. The smaller of the value , the smaller  
	      number of clusters you will obtain in the clustering analysis 
			options: 0.8 (default) 
			 			
	[cluster_pc]: The number of PCs (principal components) used for construction of shared neighbor  
	      network (SNN) graph 
			options: 30 (default) 
			 			
	[experiment_name]: Name of the your analysis which is useful to name the name of output files
			options: scCUT&Tag (default)
			 			

```

- The **run_pipeline** section includes the parameters controlling how to process the data (run the entire pipeline or an individual analysis step).  

``` 
PARAMETERS
----------

	[entire_pipeline]: Whether to process the data with the entrie pipeline. If this paramter is set   
	      as " TRUE", the other parameters in this section will not be used 
			options: TRUE (default), FALSE
						
	[individual_step]: If entire_pipeline is specified as "FALSE", this parameter will be activated   
	      to run an individual step of the data processing pipeline. 
			options: It can be specified as one of fastq2peak (step 1 raw data  
				 processing and quality assessment), cells2count_matrix (step 2 
				 feature-by-cell matrix construction), count_matrix2clustering (step 3
				 dimensionality reduction and clustering analysis) and cells2psuedoBulk 
				 (step 4 cell-type-specific pseudo-bulk data analysis). NO default 
			 			 
	[step2_bamfile_dir]: A folder containing individual bam files of single cells which are used for   
	      feature-by-cell matrix construction. If individual_step is specified as cells2count_matrix,
	      this parameter will be activated. 
			options: NO default 	
			 		
	[step2_output_dir]: A folder for output. If individual_step is specified as cells2count_matrix,   
	      this parameter will be activated.
			options: NO default

	[step2_qc_pass_file]: A file containing the cell names which are used for cell filtration (see   
	       Tutorial for further information). It is an optional parameter and the user can set it as  
	       "FALSE" if all the cells are used for feature-by-cell matrix construction. If individual_step  
	       is specified as cells2count_matrix, this parameter will be activated.
			options: FALSE

	[step3_count_matrix]: The input file of feature-by-cell matrix. If individual_step is specified   
	      as count_matrix2clustering, this parameter will be activated.
			options: NO default 

	[step3_output_dir]: A folder for output. If individual_step is specified as count_matrix2clustering,  
	       this parameter will be activated. 
			options: NO default

	[step4_bamfile_dir]: A folder containing individual bam files of single cells that are used for   
	       cell aggregation. If individual_step is specified as cells2psuedoBulk, this parameter will  
	       be activated.  
			options: NO default

	[step4_cell_anno_file]: A file containing the cell names and their corresponding labels which are  
	       used for cell aggregation (see Tutorial for further information) 
			options: NO default

	[step4_output_dir]: A folder for output. If individual_step is specified as cells2psuedoBulk, this  
	       parameter will be activated.  
			options: NO default 
			 			
	

```
**NOTES**: If an individual analysis step is performed, the parameters in the sc\_parameters section should be tuned accordingly as needed. For step1, no parameters in the run\_pipeline section need to be changed except for *entire\_pipeline* and *individual_step*. For step2, *the matrix\_type* only accepts one of bin\_by\_cell and customFeature\_by\_cell.



## Quickstart
### Easy 1-step execution

We used a published single-cell CUT&Tag dataset which detected the H3K27me3 signal of individual H1 and K562 cells ([Kaya-Okur, et al. Nat Commun 10, 1930 (2019)](https://doi.org/10.1038/s41467-019-09982-5)), as an example to illustrate the utility and functionality of our tool. To walk through the pipeline quickly, we make an small example dataset only containing data from chr11 of 200 cells (can be found in the folder exampleData), by randomly selecting 100 cells for the two cell types, respectively. It will take roughly less than half an hour to run the entire pipeline from raw sequencing data using a laptop with 8 cores.

Uncompress the exampleData.tar.gz

```
tar -zxvf exampleData.tar.gz
```

To configure the exampleData_configure.sh file: setting the *workdir*; specify the *fastqdir* parameter as the folder of exampleData, then you can simply run the following command to perform the entire pipeline.

```
chmod +x ./run_scModule.sh   
./run_scModule.sh /path/to/scCutRunTools/src /path/to/exampleData_configure.sh
```

What happens next is that CUT&RUNTools will sequentially run the analysis pipeline (for details see the Usage and Tutorial below). 


## Tutorial

### **Step 1. raw data processing and quality assessment**
***Read Triming***

All the paired fastq files in the specified *fastq_directory* are trimmed via the two-step read trimming process (details see the bulk data trimming) to identify and trim any adapter sequence. All the fastq files are trimed in parallel on the available cores.

***Read Mapping***

The trimmed paired-end reads are aligned to a reference genome (hg38, hg19, mm10 or mm9) using bowtie2 software. To download proper indexes of bowtie2, just use the either the downloads on the [Bowtie2 homepage] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or the [Illumina iGenomes] (https://support.illumina.com/sequencing/sequencing_software/igenome.html).  All the fastq files are aligned in parallel on the available cores. The resulting data can be found in the folder *sc_trimmed* and *sc_trimmed3*.

***Read Filtering***

For each single cell experiment, unmapped, unpaired reads are removed.
The CUT&RUN and CUT&Tag technologies naturally have relatively low background signal comparing to traditional ChIP-seq, high sequencing depths for individual cells could result in severely redundant sequencing of PCR duplicates. For each cell, the fragments with identical start and end positions are counted once. 
Only Read pairs with mapping quality (MAPQ) > 30 and properly paired are retained for each barcode cell. (proper paired means reads are in Read1 forward, Read2 reverse orientation or Read1 reverse, Read2 forward orientation). All the bam files are processed and filtered in parallel on the available cores. The resulting data can be found in the folder *sc_aligned.aug10*.

***Data Aggregation and Peak Calling***

For single-cell epigenome data analysis, cell aggregation and the following pseudo-bulk peak calling are common methods to assess the data quality of the experiment and gain insight into finding potential biological functional regions. CUT&RUNTools aggregate all the barcode cells by merging all the single-cell bam files. Then we performed the peak calling analysis with the resulting pseudo-bulk bam. For the user convenience to get a suitable peak set for the functional interpretation of their data, the CUT&RUNTools will perform the peak calling using macs2 software with narrow peak mode (default) or broad peak mode for the pseudo-bulk bam files for either all the cell or cells from a specific cluster. Owing to the unique character of low background signal for CUT&RUN or CUT&Tag data comparing to the traditional ChIP-seq data, we also provide the option for the peak calling with SEACR method with stringent mode, which is specifically designed for the peak calling of such data. Users can change the setting in the configuration JSON file.
The resulting data could be found in the folder *sc_pseudoBulk*.

***QC Report, Diagnositic Figures and Cell Filtering***

To facilitate the user to quickly assess the data quality and general information for the experiment, CUT&RUNTools will summarize the mapping statistics and produce a quality control (QC) report for this experiment.  
The QC report (*qc_report.txt*) contains two parts including overall statistics and average statistics for all the barcode cells. The following QC report shows the general information of our example dataset.
 
| Overall Experimental Summary    |  |
| :----------------------- | :--- |
| Total Barcode Cells | 1,373 |
| Total Reads | 192,427,780 |
| Reads per Barcode Cells | 140,151 |
| Insert Size on Average   <img width=160/> | 230.3     <img width=80/>         |
| **Statistics of Barcode Cells (Median)** |  |
| Alignment Ratio | 99.4% |
| Reads Properly Paired (%) | 98,282 (92.1%) |
| Reads High Quality (%) | 91509 (87.4%) |
| Reads Duplicated (%) | 301 (0.3%) |
| Reads Nuclear (%) | 104,834 (99.5%) |
| Reads in Peaks (%) | 48,475 (45.3%) |


To intitutive show the overall QC metrics across all barcode cells, CUT&RUNTools also generate several diagnositic figures. 
<br/>


<div align=left> <img src="./images/QC-statistics-1.png" width="600" height="300"> </div> 

The insert/fragment size distribution is useful to inspect the nucleosomal periodicity for different TF factors and histone modification. Distribution of insert size of the total unique fragments for all the barcode cells will be showed.  

<div align=left> <img src="./images/QC-statistics-3-01.png" width="450" height="350"> </div> 

Strict quality control of single-cell data is essential to remove the contribution of low-quality cells. To ensure high-quality chromatin profiles used for the following analysis, we use two important metrics (*percentage\_rip* and *num\_reads_threshold*) to filter the cells. Cells with very few levels of reads found in the peaks may need to be excluded due to extreme low sequencing depth, technical artifacts and doublets in the single-cell sequencing. Generally, the peak from pseudo-bulk analysis of cell aggregation should be representative across the majority of cell types in the experiment, which can be a direct measurement for the signal-to-noise ratio of single-cell experiments. We recommend cells should be removed when the number less than a certain level (suggest 30% as default). Users can easily change the parameters in the configuration JSON file. If the users do not want to filter any cells, they just set both the parameters of *percentage\_rip* and *num\_reads_threshold* as 0. A scatter plot of the fraction of unique fragments in peaks versus the total number of unique fragments for each barcode cell. The default thresholds for the two parameters are indicated as dash lines in the figure.  

<div align=left> <img src="./images/QC-statistics-2-01.png" width="450" height="350"> </div> 

Also, if you want to see the detailed QC metrics information of all the barcodes, you can find it in the file *statistics\_QCpassed.txt* and *statistics\_QCfailed.txt*, which were generated based on the above two metrics. The file *statistics\_QCpassed.txt* is important to guide which cells could be used in the following analysis (based on the column Barcode\_name).
Both the files contain similar information for individual cells successfully ann failed to pass the QC filter criteria. Each row represents a barcode cell and each column represents the statics of QC metrics. The column name includes 

- Barcode\_name
- Total\_reads
- Overall\_alignment\_ratio
- Reads\_properly\_paired
- Reads\_properly\_paired\_percentage
- Reads\_duplicated
- Reads\_duplicated\_percentage
- Insert\_size\_average
- Insert\_size\_sd
- Reads\_nuclear
- Reads\_nuclear\_percentage
- Reads\_in\_peak
- Reads\_in\_peak\_percentage
- Reads\_MAPQ30
- Reads\_MAPQ30\_percentage
 
The resulting data can be found in the folder *sc_qc*.

### **Step 2. feature-by-cell matrix construction**

Recover informative features from sparse single-cell data is important to assess variability between cells in single-cell epigenome data analyses. To enable comprehensively capture biological relevant features to understand cell heterogeneity, we provide three optional methods (*peak\_by\_cell*, *bin\_by\_cell*, and *customFeature\_by\_cell*) for informative features generation to build the feature-by-cell matrix.

***bin\_by\_cell***

By default, genome-wide bins are taken as the input features which consider all the read coverage for each cell. One major merit of using genome-wide bins as input features is that it considers all the read coverage from each cell, not just only a fraction of reads occurs within peaks identified from bulk data or the aggregation data. Without any prior-knowledge for datasets analyzed, it is an unbiased method to build the count matrix.  
The genome-wide bins are automatically generated, and will be used as the features to build the count matrix. The *bin_size* parameter (5kb as default) is activated and could be changed by the user to generate bins with a different resolution. 

***peak\_by\_cell***

The peak file generated from the aggregation of all the cells is used when *peak\_by\_cell* is specified. This is straight forward manner to integrate prior knowledge to extract informative genomic regions with high signal-to-noise ratios. The peak file used as input of count matrix construction is closely associated with parameters *peak\_caller* and *peak\_type*. Note that peaks called from the pseudo-bulk data can serve as the useful features only when a sufficient number of single-cell profiles in the experiment are analyzed so that the pooled data can generate the robust aggregate signal. It is advisable to use to fast analyzing the data to gain the initial insights into major cell types but not focusing on the rare cell types. 

***customFeature\_by\_cell***

The *feature_file* should be specified by the user to construct the count matrix when the custom feature is used. It can typically be pre-defined genomic regions of peaks or biological relevant genomic regions derived from bulk data. The user-interested regions can sever as an important extended for different scenarios that cannot be covered by the above-mentioned strategies. For example, the user may interest using the peaks from the well-established bulk data from ENCODE, the functional genomic regions after feature selection, or the k-mer content and motif sequences of an interested TF or other co-factors. Also, the peak sets of pseudo-bulk data from particular cell clusters can be used as custom features for data re-analysis. 
Generally, a typical bed format data can be used for input of custom features. Here is an example:

```
chr1 10973  11152  
chr1 20032  20494  
chr1 121028 123022  
...
```

***Feature Filtration***

A number of blacklist regions (including ENCODE blacklisted regions, mitochondrial and other uninterst chromosomes) are collected for feature filtration. The [ENCODE blacklisted genomic regions](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/) are showed to prone to have artifact signal in next-generation sequencing experiments independent of cell line or experiment. 
The features overlapped with blacklist regions will be removed from downstream analysis. The removal of these regions eliminated significant background noise and potential artifacts, which will result in a more biologically meaningful interpretation of the single-cell data.  
The CUT&RUN experiment using a Protein A-Protein G-Micrococcal Nuclease (pAG-MNase) to cleave DNA surrounding binding/modification sites which has a sequence-dependent cut bias of ATAT or TATA. Many CUT&RUN data have reported the contamination of undetermined origin by sequences of (TA)n. If *experiment_type* is set to *CUT&RUN*, the feature will be further filtered based on (TA)n repeat sequences. 

***Count Matrix Construction***

Once the feature filtering is completed, all the bam files of individual cell are counted against the feature file to generate the feature-by-cell matrix in parallel on the available cores using software bedtools. The resulting files can be found in the *sc_countMatrix* folder.

### **Step 3. dimensionality reduction and clustering analysis**
***Binarization and Filtration***

Owing to the sparsity of count matrix of single-cell epigenomic data, CUT&RUNTools converts the feature-by-cell matrix into the SparseMatrix format and processing it by using the *‘Matrix’* package, which allows more efficiency of memory usage and computation. Then, the sparse matrix will be binarized by specifying a feature with value “1” when one or more reads fall in, and the value “0” otherwise. Two types of features are identified including rare ones, which defined as features detected only by very few cells, i.e. the maximum of 2 cells and 1/1000 of total cell numbers, and ubiquitous features, which are defined as features detected by most of the cells, i.e. 80% of total cells. They are removed in the binarized matrix as they perhaps due to alignment error or be random noise in the experiment and contribute very little to the dimensionality reduction analysis. This operation is also very useful to effectively reduce the dimensions of the feature-by-cell matrix and thus improves the efficiency of our pipeline.

***Normalization of Binarized Matrix***

The term frequency-inverse document frequency (TF-IDF) method is used for data normalization of binarized matrix with sequencing depth per cell and the total number of unique cell coverage per feature. 

***Dimensionality Reduction***

The latent semantic indexing (LSI) method is used for dimension reduction analysis, which conducted the principal component analysis (PCA) method against the TF-IDF normalized data to project the cells into a space of fewer dimensions. This is commonly used for single-cell epigenome data procession. CUT&RUNTools visualizes the results employing three commonly used methods for dimensionality reduction including Principal Component Analysis (PCA), Uniform Manifold Approximation and Projection (UMAP) and t-distributed stochastic neighbor embedding (t-SNE). The corresponding plots of the first two dimensions will be generated and saved, respectively. 

***Clustering Analysis***

The graph-based Leiden algorithm implemented with python package *leiden*, which is used as the unsupervised clustering method. Using the top 30 principal components, we built a Shared neighbor network (SNN) graph by considering each cell as a node and further finding its k-nearest neighbors according to the Euclidian distance. Leiden algorithm will be performed with resolution 0.8 (important for the determination of cluster number) as default. Users can use different thresholds of the PCs and resolution for the clustering analysis with different parameters of *cluster_pc* and *cluster_resolution*, respectively. 

***Cell Annotation and Data Visualization***

The file named *leiden_cluster_annotation.txt* could be produced stores cell annotations gained from the clustering analysis. It is important to guide the following cell aggregation analysis. Here is an example:

```
cell_name	leiden_cluster
SRX5204736.bam	1
SRX5204737.bam	1
SRX5204738.bam	2
SRX5204739.bam	2
...
```

Also, CUT&RUNTools will output three figures (both of the pdf and png versions) to show the results of dimensionality reduction and clustering analysis using PCA, tSNE and UMAP methods, respectively. A plot of two-dimensional projection of the scCUT&Tag data using t-SNE method:  
<br/>
<br/>

<div align=left> <img src="./images/tSNE-plot.png" width="450" height="350"> </div> 

The TF-IDF normalized count matrix for the experiment is saved and could be used for the further analysis. The resulting data can be found in the *sc_cluster* folder. 

### **Step 4. cell-type-specific pseudo-bulk data analysis**

***Cell Aggregation and Peak Calling***

The qualified profiles of individual cells from the same cell cluster (according to the cell annotation from clustering analysis) will be merged to generate the ensemble pseudo-bulk bam. The peak calling is performed with aggregated bam file and all the statistically significant peaks are generated for visualization and downstream analysis.

***The Single-cell/Pseudo-bulk Genomic Tracks***

To simultaneously visualize the genomic binding profiles of hundreds of thousands of single-cells, we use a new method to integrate information of interested cells into a single track file. The [qBED format](https://www.biorxiv.org/content/10.1101/2020.04.27.060061v2.full) is a text-based standard originally for encoding and analyzing calling card data and the qBED files can be visualized in the WashU Epigenome Browser. We found it can also efficiently illustrate the quantitative information of coverage information of hundreds of thousands of single cells simultaneously.  
For each single cell, we generated the bedgraph file to represent the coverage profiles and convert them to qBED format. The qBED files of cells grouped/annotated with the same cell cluster are sorted and combined into a new pooled qBED file. Then the pooled qBED is compressed using software tabix.  

Users can visualize the cell-type-specific signals of single cells by uploading both the qBED.gz and qBED.tabi files to the [WashU Epigenome Browser](https://epigenomegateway.wustl.edu/) via Local Tracks. The following information can be copied and pasted when uploading the files to configure the tracks. User can adjust the height/color/opacity according to cell number for visualization. 

```
{  
    "color":"#D12134", 
    "height":40,  
    "opacity":[80],
    "show":"ALL",
    "markerSize":1,
    }
 
```

The genome tracks of pseudo-bulk samples for each cell cluster are generated from the corresponding bam files. The aggregation data is normalized by the counts per million (CPM) method and the bigwig format files are generated using software deeptools. The resulting files of pooled signal (files suffix with .bw and .bedgraph) and the corresponding DNA binding sites (files suffix with .bed, .narrowPeak, or .broadPeak) can be uploaded together with single-cell qBED tracks to the genome browser for visualization and comparison.
The file suffix with .bed for each cell storing the genome coverage information is also generated for customizing visualization. The resulting file could be found in the *sc_pseudoBulk* folder.  
A screenshot of genome browser view is displayed.  
<br/>

<div align=left> <img src="./images/gb-chr17-48513628-48855563.png" width="1600" height="300"> </div> 


<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>

### **Downstream analysis of cell-cluster-specific peaks uisng Bulk-level analysis tools**
> ***Differential binding peaks between cell clusters***

Users can identify differential features among different cell groups based on the annotation from meta_table provided by the user or unsupervised analysis. The signals across all the cells for top features can be visualized. This function is implemented via FindAllmarker functions in Seurat. Wilcoxon test and top 10 for each cell clusters are retained.
TODO 
For each pseudo-bulk cluster, it finds the positive peaks of a single cluster (specified in cluster.pos), comparing to the other pseudogroups. This function is implemented via software diffbind. If by default cluster.neg=NULL, findDAR will look for a group of background cells closest to the positive cells. identify master regulators that are enriched in the differentially accessible regions (DARs)


> ***Motif discovery of (co-)regulatory TF***

We performed motif discovery using homer to infer the potential master regulator that controls for gene expression in each of the cell types.
If experiment_type is set to CUT&RUN, the feature will be further filtered based on (TA)n repeat sequences. 

We also filter the peak with blacklist regions before TF motif discovery and footprinting in the pseudo-bulk analysis. We also showed that simple repeat regions removal is essential to prevent from potential artifacts in the TF motif analysis for multiple single-cell cut&run datasets. 

> ***TF footprinting analysis***



> ***Functional annotation***


Motif scanning is carried out by using the FIMO package of the MEME suite (58) using the standard threshold (P < 1E-4). The background distribution is estimated by randomly selecting regions with matching C+G content. Known TF motifs are obtained from the JASPAR (35), TRANSFAC (34), and FACTORBOOK (36) databases. Statistical significance of motif enrichment is determined using Fisher’s exact test




## Directory structure of input and output folders 

```
	+ fastq_directory
	+ workdir
	    + sc_trimmed
	    + sc_trimmed3
	    + sc_log
	    	+ trim
	    	+ bowtie2
	    	+ bamSort
	    	+ bamMarkdup
	    + sc_aligned.aug10
	    	+ sorted
	    	+ dup.marked
	    	+ dup.marked.clean
	    + sc_countMatrix
	    + sc_cluster
	    + sc_qc
	    	+ report
	    	+ figure
	    + sc_pseudoBulk
	        + groups_aggregation
	            + single_cells
	            + pseudo_bulk_data
	                +  macs2
	                +  SEACR
	        + group_1
	            + single_cells
	            + pseudo_bulk_data
	        + group_2
	        	+ single_cells
	        	+ pseudo_bulk_data
  
```



Install haystack\_bio (it can highlight enriched TF motifs in variable and cell-type specific regions, more detail information:     https://github.com/pinellolab/haystack\_bio); get the path using the command: *which haystack_bio*

```
conda install haystack_bio
```
    
