# CUT&RUNTools Installation

CUT&RUNTools requires Python 2.7, R 3.3, Java, Perl 5 and [Slurm](https://slurm.schedmd.com/) job submission environment. It is designed to run on a cluster set up.

Installation also requires GCC to compile some C-based source code. 

## Prerequisites

The following tools should be already installed. Check the corresponding website to see how to install them if not. Special notes for Atactk and UCSC-tools, see below.

In the bracket is the version we have. CUT&RUNTools may work with a lower version of each tool. 

* Trimmomatic (0.36) [link](http://www.usadellab.org/cms/?page=trimmomatic)
* Bowtie2 (2.2.9) [link](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Samtools (1.3.1) [link](http://samtools.sourceforge.net/)
* Picard (2.8.0) [link](https://broadinstitute.github.io/picard/)
* MACS2 (2.1.1) [link](https://github.com/taoliu/MACS)
* MEME (4.12.0) [link](http://meme-suite.org/tools/meme)
* Bedops (2.4.30) [link](https://bedops.readthedocs.io/en/latest/)
* Bedtools (2.26.0) [link](https://bedtools.readthedocs.io/en/latest/)
* Atactk [link](https://github.com/ParkerLab/atactk)- we provide special install instructions since we need to patch a source file.
* UCSC-tools [link](http://hgdownload.soe.ucsc.edu/admin/exe/)- we provide special install instructions.

Other tools already contained in CUT&RUNTools:

* Samblaster [link](https://github.com/GregoryFaust/samblaster)
* SEACR [link](https://github.com/FredHutch/SEACR)


### Atactk

This is a [python2 package](https://github.com/ParkerLab/atactk) that determines the enzyme cut frequency matrix. Originally for Tn5 transposase in ATAC-seq, the logic of the tool also applies to other digestions like CUT&RUN. 
However since original implementation is designed for ATAC-seq it needs to be patched in one small place in the code to make it suitable for CUT&RUN analysis. Otherwise, we estimate that the cut frequency is sometime off by 1 bp in its calculation.

We provide a patch for this problem. Install the patched version of the package by reading [`atactk.install.sh`](atactk.install.sh). The patches [`make_cut_matrix.patch`](make_cut_matrix.patch) and [`metrics.py.patch`](metrics.py.patch). Then install by:

```
source atactk.install.sh
```
This will use pip to install atactk to the user's home directory (~/.local/bin).


### UCSC-tools

We need two specific tools bedGraph2BigWig and fetchChromSizes, and provide a script [`ucsc-tools.install`](ucsc-tools.install) to automatically download and install them from the UCSC genome browser:
```
source ucsc-tools.install
```
This will download the two executables in the current directory.


### A note about Anaconda

If you are a conda/bioconda user, you can already obtain the pre-requisite softwares above from the bioconda package repository. On a linux system, the anaconda environment might look like the following:

```
conda info -e
# conda environments:
#
base                     *  /home/qzhu/anaconda2
bedGraphToBigWig            /home/qzhu/anaconda2/envs/bedGraphToBigWig
bedops                      /home/qzhu/anaconda2/envs/bedops
bedtools                    /home/qzhu/anaconda2/envs/bedtools
bowtie2                     /home/qzhu/anaconda2/envs/bowtie2
deeptools                   /home/qzhu/anaconda2/envs/deeptools
macs2                       /home/qzhu/anaconda2/envs/macs2
samtools                    /home/qzhu/anaconda2/envs/samtools
trimmomatic                 /home/qzhu/anaconda2/envs/trimmomatic
```

You still need to follow the custom Atactk and UCSC-tools instructions above.

### A note about other Linux systems

If you are using Ubuntu Linux 18.04 or above, obtaining and installing the pre-requisites should be in most cases pretty simple. 
You can either use the `apt-get install` command to install things like `bedops`, `bedtools`, `bowtie2`. For other packages that are not in the Ubuntu package manager, 
they can be installed by going to the official channel, downloading the source tar file, extracting it, and doing `./configure --prefix=/usr/local` followed by `make && make install`.

On my system, which is the Harvard Medical School O2 Computing Cluster, these pre-requisites have been installed by the system admin (see `module spider` on O2), so there is no need to worry about them
if you are a user from the HMS community.


## Configuration file

The configuration file tells CutRunTools where to locate the prerequisite tools. This is a [JSON](http://www.json.org/) file. A sample JSON file is provided below ([`config.json`](config.json)). Filling in the information should be pretty easy: in most cases we need to provide the path to the `bin` directory of each tool.

```json
{
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
	"extrasettings": "/home/qz64/cutrun_pipeline", 
	"input/output": {
		"fastq_directory": "/n/scratch2/qz64/Nan_18_aug23/Nan_run_19",
		"workdir": "/n/scratch2/qz64/workdir",
		"fastq_sequence_length": 42,
		"organism_build": "hg19"
	},
	"motif_finding": {
		"num_bp_from_summit": 150,
		"num_peaks": 5000,
		"total_peaks": 15000,
		"motif_scanning_pval": 0.0005,
		"num_motifs": 20
	},
	"cluster": {
		"email": "johndoe@gmail.com",
		"step_alignment": {
			"queue": "short",
			"memory": 32000,
			"time_limit": "0-12:00"
		},
		"step_process_bam": {
			"queue": "short",
			"memory": 32000,
			"time_limit": "0-12:00"
		},
		"step_motif_find": {
			"queue": "short",
			"memory": 32000,
			"time_limit": "0-12:00"
		},
		"step_footprinting": {
			"queue": "short",
			"memory": 32000,
			"time_limit": "0-12:00"
		}
	}
}
```
Pay attention to the **the first 20 lines of this file** which concern the software installation. The rest is related to an actual analysis (explained in [USAGE.md](USAGE.md)). 


### Validate prerequisites

To check if the paths are correct and if the softwares in these paths indeed work:
```
./validate.py config.json --ignore-input-output --software
```

### Common errors flagged by validate.py

#### Validation of macs2, make_cut_matrix failed
This could be due to that the versions of Python are inconsistent. The version of Python indicated in the field `pythonbin` should be consistent with the version of Python used to install MACS2, and version used to install make_cut_matrix (`makecutmatrixbin`). 
To confirm this, on our system, we enter
```bash
head -n 1 /n/app/macs2/2.1.1.20160309/bin/macs2
#!/n/app/python/2.7.12/bin/python
head -n 1 /home/qz64/.local/bin/make_cut_matrix
#!/n/app/python/2.7.12/bin/python
```
Use `/n/app/python/2.7.12/bin` for the `pythonbin` field in `config.json`. 
Similarly for Perl (`perlbin`) and meme-chip (`memebin`). Version of Perl used to install meme should be used for `perlbin`:
```bash
head -n 1 /home/qz64/meme/bin/meme-chip
#!/n/app/perl/5.24.0/bin/perl
```
Use `/n/app/perl/5.24.0/bin` as value for `perlbin` in `config.json`.

#### Validation of meme, meme-chip failed, XML/Simple.pm not found
Validate script may flag an error with MEME installation: XML/Simple.pm is not found.
You may encounter this error if the XML PERL module is installed to a custom directory (i.e. home rather than in `/usr/` or `/usr/local`). 
If this is the case, the solution is to set up the Perl library environment variables.
On our system, the custom Perl modules are installed to `/home/qz64/perl5-O2/lib`, so we enter:
```bash
PERL5LIB="/home/qz64/perl5-O2/lib/perl5":$PERL5LIB
PERL_LOCAL_LIB_ROOT="/home/qz64/perl5-O2":$PERL_LOCAL_LIB_ROOT
export PERL5LIB
export PERL_LOCAL_LIB_ROOT
```
Then try validate script again.


## Download genome assemblies

The genome sequence of a specific organism build (such as hg19, hg38) is required for motif finding. We provide a script [`assemblies.install`](assemblies.install) to download this automatically from NCBI. We specifically require repeat-**masked** version of genome sequence file. See [`assemblies.install`](assemblies.install) for more details. 
```
source assemblies.install
```

## Install CENTIPEDE

In R, run:
```
install.packages("CENTIPEDE", repos="http://R-Forge.R-project.org")
```

## Install CutRunTools

We first install a special trimmer we wrote `kseq`. This tool further trims the reads by 6 nt to get rid of the problem of possible adapter run-through. To install:
```
source make_kseq_test.sh
```

Once completed, CutRunTools is installed.

CutRunTools does not need to be installed to any special location such as /usr/bin, or /usr/local/bin. 

It can be run directly from the directory containing the CutRunTools scripts.

The main program consists of `create_scripts.py`, `validate.py`, and a set of scripts in `aligned.aug10`, and in `macs2.narrow.aug18` which perform the important motif finding and footprinting analyses.

We suggest placing the CUT&RUNTools scripts to a more permanent location (like `~/cutruntools-scripts` or `/usr/local/cutruntools-scripts`). For example:

```
#suppose this is where cutruntools scripts are downloaded to:
pwd
/home/qz64/Downloads/qzhudfci-cutruntools-1806e65e5b28
#copy it in a more permanent location like home directory
mkdir ~/cutruntools-scripts
cp -r /home/qz64/Downloads/qzhudfci-cutruntools-1806e65e5b28/* ~/cutruntools-scripts/.
```

See [USAGE.md](USAGE.md) for details. Briefly, a user writes a `config.json` configuration file for a new analysis. CutRunTools uses this to generate a set of slurm-based job-submission scripts, customized to the user's samples and environment. These job-submission scripts can be directly used by the user to start analyzing his/her Cut&Run samples.

```
#Example
#Use the config.json from ~/cutruntools-scripts as a template for a new analysis
cp ~/cutruntools-scripts/config.json config.jan4.json
#modify config.jan4.json to your samples and needs
~/cutruntools-scripts/create_scripts.py config.jan4.json
#the SLURM-based job submission scripts would now have been created. See USAGE.md how to continue.
```