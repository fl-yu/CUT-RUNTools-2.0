# Software installation of bulk data module

CUT&RUNTools 2.0 requires Python (3.6), R (3.6), Java and Perl. Installation also requires GCC to compile some C-based source code. 

## Prerequisites

Additionally, the following tools should be already installed. Check the corresponding website to see how to install them.  CUT&RUNTools 2.0 may work with a lower version of each tool. In the bracket is the version we have.

* Trimmomatic (0.36) [link](http://www.usadellab.org/cms/?page=trimmomatic)
* Bowtie2 (2.2.9) [link](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Samtools (1.3.1) [link](http://samtools.sourceforge.net/)
* MACS2 (2.1.1) [link](https://github.com/taoliu/MACS)
* MEME (4.12.0) [link](http://meme-suite.org/tools/meme)
* Bedops (2.4.30) [link](https://bedops.readthedocs.io/en/latest/)
* Bedtools (2.26.0) [link](https://bedtools.readthedocs.io/en/latest/)
* Deeptools (3.5.0) [link](https://deeptools.readthedocs.io/en/develop/)

Alternatively, we recommend the user to install the software using [conda system](https://docs.conda.io/en/latest/). On a Linux system, paths of software from the anaconda environment might be obtained in the following:

```
conda info -e
# conda environments:
#
base                     *  /home/fyu/anaconda2
bedGraphToBigWig            /home/fyu/anaconda2/envs/bedGraphToBigWig
bedops                      /home/fyu/anaconda2/envs/bedops
bedtools                    /home/fyu/anaconda2/envs/bedtools
bowtie2                     /home/fyu/anaconda2/envs/bowtie2
deeptools                   /home/fyu/anaconda2/envs/deeptools
macs2                       /home/fyu/anaconda2/envs/macs2
samtools                    /home/fyu/anaconda2/envs/samtools
```

We also provided special notes for **Atactk** and **UCSC-tools**, see below.
* Atactk [link](https://github.com/ParkerLab/atactk)- we provide special install instructions since we need to patch a source file.
* UCSC-tools [link](http://hgdownload.soe.ucsc.edu/admin/exe/)- we provide special install instructions.

Other tools already contained in CUT&RUNTools:

* Samblaster [link](https://github.com/GregoryFaust/samblaster)
* SEACR [link](https://github.com/FredHutch/SEACR)
* picard (0.1.8) [link](http://broadinstitute.github.io/picard/command-line-overview.html)
* trimmomatic (0.36) [link](https://github.com/timflutre/trimmomatic)

Files already contained in CUT&RUNTools:
* genome size files
* blacklist regions
* adaptor files
* example fastq data

### Atactk installation

This is a [python package](https://github.com/ParkerLab/atactk) that determines the enzyme cut frequency matrix. Originally for Tn5 transposase in ATAC-seq, the logic of the tool also applies to other digestions like CUT&RUN. 
However since original implementation is designed for ATAC-seq it needs to be patched in one small place in the code to make it suitable for CUT&RUN analysis. Otherwise, we estimate that the cut frequency is sometime off by 1 bp in its calculation.

We provide a patch for this problem. Install the patched version of the package by reading [`atactk.install.sh`](atactk.install.sh). The patches [`make_cut_matrix.patch`](make_cut_matrix.patch) and [`metrics.py.patch`](metrics.py.patch). Then install by:

```
source atactk.install.sh
```
This will use pip to install atactk to the user's home directory (~/.local/bin).

*  `adapterpath` contains Illumina Truseq3-PE adapter sequences (we provide them). 

### UCSC-tools installation

We need two specific tools bedGraph2BigWig and fetchChromSizes, and provide a script [`ucsc-tools.install`](ucsc-tools.install) to automatically download and install them from the UCSC genome browser:
```
source ucsc-tools.install
```
This will download the two executables in the current directory.

### kseq installation

We first install a special trimmer we wrote `kseq`. This tool further trims the reads by 6 nt to get rid of the problem of possible adapter run-through. To install:
```
source make_kseq_test.sh
```

## Download genome assemblies

The genome sequence of a specific organism build (such as hg19, hg38) is required for motif finding. We provide a script [`assemblies.install`](assemblies.install) to download this automatically from NCBI. We specifically require repeat-**masked** version of genome sequence file. See [`assemblies.install`](assemblies.install) for more details. 
```
source assemblies.install
```

## Install CENTIPEDE R packages for footprinting analysis

In R, run:
```
install.packages("CENTIPEDE", repos="http://R-Forge.R-project.org")
```


### A note about other Linux systems

If you are using Ubuntu Linux 18.04 or above, obtaining and installing the pre-requisites should be in most cases pretty simple. 
You can either use the `apt-get install` command to install things like `bedops`, `bedtools`, `bowtie2`. For other packages that are not in the Ubuntu package manager, 
they can be installed by going to the official channel, downloading the source tar file, extracting it, and doing `./configure --prefix=/usr/local` followed by `make && make install`.


## Configuration file

The configuration file tells CutRunTools where to locate the prerequisite tools. This is a [JSON](http://www.json.org/) file. A sample JSON file is provided below ([`bulk-config.json`](bulk-config.json)). Filling in the information should be pretty easy: in most cases we need to provide the path to the `bin` directory of each tool.

```json
{
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
	"macs2pythonlib": "/path/to/macs2/lib/python2.7/site-packages",
	"kseqbin": "/path/to/cutrun_pipeline", 
	"memebin": "/path/to/meme/bin", 
	"bedopsbin": "/path/to/bedops/bin", 
	"bedtoolsbin": "/path/to/bedtools/bin",
	"makecutmatrixbin": "/path/to/.local/bin",
	"bt2idx": "/path/to/bowtie2_indexes",
	"genome_sequence": "/path/to/chrom.hg19/hg19.fa",
	"extratoolsbin": "/path/to/cutrun_pipeline", 
	"extrasettings": "/path/to/cutrun_pipeline", 
	"input/output": {
		"fastq_directory": "/path/to/user/fastq_directory",
		"workdir": "/path/to/user/workdir",
		"fastq_sequence_length": 42,
		"organism_build": "hg19"
	},
	"motif_finding": {
		"num_bp_from_summit": 150,
		"num_peaks": 5000,
		"total_peaks": 15000,
		"motif_scanning_pval": 0.0005,
		"num_motifs": 20
	}
}
```
Pay attention to the **the first 20 lines of this file** which concern the software installation. The rest is related to an actual analysis (explained in [USAGE.md](USAGE.md)). 


### Validate prerequisites

To check if the paths are correct and if the softwares in these paths indeed work:
```
./validate.py bulk-config.json --ignore-input-output --software
```
This script checks that your configuration file is correct and all paths are correct. You will get an empty line if the validate.py script runs without errors.

## Install



CutRunTools does not need to be installed to any special location such as /usr/bin, or /usr/local/bin. 

It can be run directly from the directory containing the CutRunTools scripts.

The main program consists of `create_scripts.py`, `validate.py`, and a set of scripts in `aligned.aug10`, and in `macs2.narrow.aug18` which perform the important motif finding and footprinting analyses.

We suggest placing the CUT&RUNTools scripts to a more permanent location (like `~/cutruntools-scripts` or `/usr/local/cutruntools-scripts`). For example:

```
#suppose this is where cutruntools scripts are downloaded to:
pwd
/path/to/Downloads/qzhudfci-cutruntools-1806e65e5b28
#copy it in a more permanent location like home directory
mkdir ~/cutruntools-scripts
cp -r /path/to/Downloads/qzhudfci-cutruntools-1806e65e5b28/* ~/cutruntools-scripts/.
```

See [USAGE.md](USAGE.md) for details. Briefly, a user writes a `bulk-config.json` configuration file for a new analysis. CutRunTools uses this to generate a set of slurm-based job-submission scripts, customized to the user's samples and environment. These job-submission scripts can be directly used by the user to start analyzing his/her Cut&Run samples.

```
#Example
#Use the bulk-config.json from ~/cutruntools-scripts as a template for a new analysis
cp ~/cutruntools-scripts/config.json config.jan4.json
#modify config.jan4.json to your samples and needs
~/cutruntools-scripts/create_scripts.py config.jan4.json
#the SLURM-based job submission scripts would now have been created. See USAGE.md how to continue.
```




### Common errors flagged by validate.py

#### Validation of macs2, make_cut_matrix failed
This could be due to that the versions of Python are inconsistent. The version of Python indicated in the field `pythonbin` should be consistent with the version of Python used to install MACS2, and version used to install make_cut_matrix (`makecutmatrixbin`). 
To confirm this, on our system, we enter
```bash
head -n 1 /path/to/macs2/2.1.1.20160309/bin/macs2
#!/path/to/python/2.7.12/bin/python
head -n 1 /path/to/.local/bin/make_cut_matrix
#!/path/to/python/2.7.12/bin/python
```
Use `/path/to/python/2.7.12/bin` for the `pythonbin` field in `config.json`. 
Similarly for Perl (`perlbin`) and meme-chip (`memebin`). Version of Perl used to install meme should be used for `perlbin`:
```bash
head -n 1 /path/to/meme/bin/meme-chip
#!/path/to/perl/5.24.0/bin/perl
```
Use `/path/to/perl/5.24.0/bin` as value for `perlbin` in `config.json`.

#### Validation of meme, meme-chip failed, XML/Simple.pm not found
Validate script may flag an error with MEME installation: XML/Simple.pm is not found.
You may encounter this error if the XML PERL module is installed to a custom directory (i.e. home rather than in `/usr/` or `/usr/local`). 
If this is the case, the solution is to set up the Perl library environment variables.
On our system, the custom Perl modules are installed to `/path/to/perl5-O2/lib`, so we enter:
```bash
PERL5LIB="/path/to/perl5-O2/lib/perl5":$PERL5LIB
PERL_LOCAL_LIB_ROOT="/path/to/perl5-O2":$PERL_LOCAL_LIB_ROOT
export PERL5LIB
export PERL_LOCAL_LIB_ROOT
```
Then try validate script again.



