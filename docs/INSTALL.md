# Software installation

CUT&RUNTools 2.0 requires **Python** (>=3.6), **R** (>=3.3), **Java** (>=1.8) and **Perl**. Installation also requires **GCC** to compile some C-based source code. Additionally, the following required tools should be already installed before running the setup.   

## Prerequisites
**In order to install the dependencies in Part 1 & 3, you can run the following commands to create conda environments with CUT&RUNTools's dependencies:**  
```
conda env create --file environment1.yml
conda env create --file environment2.yml
```  
**Otherwise, you could intall the dependencies separately as follow.**

**Part 1.**  
CUT&RUNTools 2.0 may work with different version of each tool.

* Bowtie2 [link](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Samtools [link](http://samtools.sourceforge.net/)
* MACS2 [link](https://github.com/taoliu/MACS)
* MEME [link](http://meme-suite.org/tools/meme)
* Bedops [link](https://bedops.readthedocs.io/en/latest/)
* Bedtools [link](https://bedtools.readthedocs.io/en/latest/)
* Deeptools [link](https://deeptools.readthedocs.io/en/develop/)
* GNU parallel [link](https://www.gnu.org/software/parallel/)
* tabix [link](https://davetang.org/muse/2013/02/22/using-tabix/)


We recommend the user to install and manage the software using [conda system](https://docs.conda.io/en/latest/) which could incorporate the dependencies for each software. 

1. Create and activate the conda environment to make a tidy environment to mange the software

```
conda create -n cutruntools2.1 python=3.6
source activate cutruntools2.1
```
2. Install the required software in the conda environment, you can skip any software which have been installed

```
conda install -y -c bioconda bowtie2
conda install -y -c bioconda samtools
conda install -y -c bioinfo macs2
conda install -y -c bioconda bedops
conda install -y -c bioconda bedtools
conda install -y -c bioconda deeptools
conda install -y -c conda-forge parallel
conda install -y -c bioconda tabix
conda install -y -c conda-forge python-igraph
```  
As the installation of MEME will frequently conflict with other software, we will create a new conda environment named `meme` to install this independently.
```
conda create -n meme python=3.6
source activate meme
conda install -y -c bioconda/label/cf201901 meme
```

On a Linux system, paths of software (e.g. /homes6/fulong/miniconda3/envs/cutruntools2) from the anaconda environment might be obtained in the following: 

```
conda env list
# conda environments:
#
base                        /homes6/fulong/miniconda3
cutruntools2             *  /homes6/fulong/miniconda3/envs/cutruntools2.1
meme                        /homes6/fulong/miniconda3/envs/meme
```

Software version is able to be obtained in the following:

```
conda env list
# conda software and version in the environment:
#
bedops                    2.4.39               hc9558a2_0    bioconda
bedtools                  2.29.2               hc088bd4_0    bioconda
bowtie2                   2.4.2            py36h5202f60_1    bioconda
deeptools                 3.5.0                      py_0    bioconda
deeptoolsintervals        0.1.9            py36h4c5857e_2    bioconda
macs2                     2.2.7.1          py36h4c5857e_1    bioconda
parallel                  20201122             ha770c72_0    conda-forge
samtools                  1.7                           1    bioconda
tabix                     0.2.6                ha92aebf_0    bioconda
```

**Part 2.**   
R packages including reticulate ("1.15"), leiden ("0.3.3"), data.table ("1.11.6"), Matrix ("1.2.18"), irlba ("2.3.3"), Rtsne ("0.15"), RANN ("2.6.1"), igraph ("1.2.4.2"), uwot ("0.1.8"), rGREAT ("1.14.0"), ggplot2 ("3.3.0"), CENTIPEDE ("1.2"), viridis, need to be installed in the library of R you specified in the JSON configuration file. The user can install these R packages automatically with the script ***r-pkgs-install.r*** in the `install` folder provided by our package.

```
# in R
source("r-pkgs-install.r")
```

**Part 3.**  
Three python modules (umap, leidenalg, and igraph) are required to be installed on your system. As igraph package was deprecated, we use conda to install igraph above.  
To install them  

```
pip install umap-learn leidenalg
```  

To test if the package is installed successfully, for example:

```
pip list | grep umap
```

**Part 4.**  

We also provided special notes for **Atactk** and **kseq**, **the installation files were already included in the packages**.

Two patches of `make_cut_matrix.patch` and `metrics.py.patch` for Atactk [(link)](https://github.com/ParkerLab/atactk) were provided to accurately estimate the cut frequency at single-base resolution. Install the patched version of the package by:

```
source atactk.install.sh
```
This will use pip to install the patched Atactk (*make_cut_matrix*) to the user's home directory (~/.local/bin).

Another software `kseq` for a special trimmer we wrote, which can further trim the reads by 6 nt to get rid of the problem of possible adapter run-through. To install:

```
source make_kseq_test.sh
```

**Part 5.**  
**[optional]** We provided script or method to download the reference genome and bowtie2 indexes. The genome sequence of a specific organism build (such as hg19, hg38) is required for genome alignment and motif discovery. We provide a script `assemblies.install` to download this automatically from UCSC. Two parameters were needed to be specified by the user, genome assembly (hg38, hg19, mm10 or mm9) and the path of software bedtools. 

```
chmod +x assemblies.install.sh
./assemblies.install.sh hg19 /path/to/bedtools
```
To download proper indexes of bowtie2, just use the either the downloads on the [Bowtie2 homepage] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or the [Illumina iGenomes] (https://support.illumina.com/sequencing/sequencing_software/igenome.html).  

**Part 6.**  
**[optional]** In this part, some necessary tools and data were already contained in CUT&RUNTools, so no installation is required.  
* SEACR (1.3) [link](https://github.com/FredHutch/SEACR)  
* picard (0.1.8) [link](http://broadinstitute.github.io/picard/command-line-overview.html)  
* trimmomatic (0.36) [link](https://github.com/timflutre/trimmomatic)  

Files frequently used which were already contained in CUT&RUNTools 2.0:
* genome size files (`assemblies` folder contains genome size files)
* blacklist regions (`blacklist` folder contains blacklist sequences)
* adaptor files (`adapters` folder contains Illumina Truseq3-PE adapter sequences)
* example fastq data (`exampleData` folder contains example data for the data analysis)


## Configuration file

The configuration file tells CutRunTools where to locate the prerequisite tools. This is a [JSON](http://www.json.org/) file. The sample JSON files `bulk-config.json` and `sc-config.json` for are provided below for bulk- and single-cell data processing and analysis, respectively. Filling in the information should be pretty easy: in most cases we need to provide the path to the `bin` directory of each tool.

**The sample JSON file of `bulk-config.json`**

```json
{
    "software_config": {
        "Rscriptbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "pythonbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "perlbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin",
        "javabin": "/homes6/fulong/miniconda3/envs/dfci1/bin",
        "bowtie2bin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "samtoolsbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "macs2bin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "memebin": "/homes6/fulong/miniconda3/envs/py3/bin", 
        "bedopsbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "bedtoolsbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "path_deeptools": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin",
        "path_parallel": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "path_tabix": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin",
        "bt2idx": "/gcdata/gcproj/fulong/Data/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index", 
        "genome_sequence": "/gcdata/gcproj/fulong/Data/Genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa", 
        "spike_in_bt2idx": "/gcdata/gcproj/fulong/Data/Genomes/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/Bowtie2Index", 
        "spike_in_sequence": "/gcdata/gcproj/fulong/Data/Genomes/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/WholeGenomeFasta/genome.fa", 
        "extratoolsbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "extrasettings": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install",
        "kseqbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "adapterpath": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/adapters", 
        "trimmomaticbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "picardbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "picardjarfile": "picard-2.8.0.jar", 
        "trimmomaticjarfile": "trimmomatic-0.36.jar", 
        "makecutmatrixbin": "/homes6/fulong/.local/bin"
    }, 
    "input_output": {
        "fastq_directory": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/exampleData", 
        "workdir": "/gcdata/gcproj/fulong/Software/cutrun-test/bulk-example-test", 
        "fastq_sequence_length": "42", 
        "organism_build": "hg19",
        "spike_in": "FALSE",
        "spike_in_norm": "FALSE",
        "spikein_scale": "10000",
        "frag_120": "TRUE",
        "peak_caller": "macs2",
        "dup_peak_calling": "FALSE",
        "cores": "8",
        "experiment_type": "CUT&RUN"
    }, 
    "motif_finding": {
        "num_bp_from_summit": "100", 
        "num_peaks": "1000", 
        "total_peaks": "2000", 
        "motif_scanning_pval": "0.0005", 
        "num_motifs": "10"
    }
}


```

The `software_config` section (the first 26 lines) concerns the software installation and required data, all the requirements of software can be defined here. The rest is related to an actual analysis (explained in [USAGE page]).


**The sample JSON file of `sc-config.json`** 

```json
{
    "software_config": {
        "Rscriptbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "pythonbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "perlbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin",
        "javabin": "/homes6/fulong/miniconda3/envs/dfci1/bin",
        "bowtie2bin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "samtoolsbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "macs2bin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "memebin": "/homes6/fulong/miniconda3/envs/py3/bin", 
        "bedopsbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "bedtoolsbin": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "path_deeptools": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin",
        "path_parallel": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin", 
        "path_tabix": "/homes6/fulong/miniconda3/envs/cutruntools2.1/bin",
        "bt2idx": "/gcdata/gcproj/fulong/Data/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index", 
        "genome_sequence": "/gcdata/gcproj/fulong/Data/Genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa", 
        "spike_in_bt2idx": "/gcdata/gcproj/fulong/Data/Genomes/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/Bowtie2Index", 
        "spike_in_sequence": "/gcdata/gcproj/fulong/Data/Genomes/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/WholeGenomeFasta/genome.fa", 
        "extratoolsbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "extrasettings": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install",
        "kseqbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "adapterpath": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/adapters", 
        "trimmomaticbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "picardbin": "/gcdata/gcproj/fulong/Software/cutrun-test/CUT-RUNTools-2.0/install", 
        "picardjarfile": "picard-2.8.0.jar", 
        "trimmomaticjarfile": "trimmomatic-0.36.jar", 
        "makecutmatrixbin": "/homes6/fulong/.local/bin"
    },
    "input_output": {
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

The `software_config` section (the first 26 lines) concerns the software installation and required data (same with the bulk version).
Similar to the configure.json of the bulk data processing, all the requirements of software can be defined here. Three more paths (path_parallel, path_deeptools and path_tabix) should be specified compared to the bulk configure.json file. The rest is related to an actual analysis (explained in [USAGE page]).



## Validate prerequisites

Briefly, a user writes a `bulk-config.json` or `sc-config.json` configuration file for a new analysis. Then, CutRunTools 2.0 can be run directly from the directory containing the CutRunTools scripts. To check if the paths are correct and if the softwares in these paths indeed work:

```
./validate.py bulk-config.json --ignore-input-output --software
```

OR

```
./validate.py sc-config.json --ignore-input-output --software
```

This script checks that your configuration file is correct and all paths are correct. You will get an empty line if the validate.py script runs without errors.








