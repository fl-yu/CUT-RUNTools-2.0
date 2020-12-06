# CUT&RUNTools 2.0 bulk data analysis

This package contains the pipeline for conducting a bulk CUT&RUN or CUT&Tag data analysis.
The pipeline comprises of read trimming, alignment steps, motif finding steps, and finally the motif footprinting step. 

To get started, please see [INSTALL.md](INSTALL.md) about how to set up the pipeline.

Once the package is installed, please see [USAGE.md](USAGE.md) about usage. Basic usage is provided below.


## What's New for bulk data analysis

11/19/20: spike-in sequence alignment and data normalizaion   
11/19/20: new options for CUT&RUN or CUT&Tag data analysis  
11/19/20: new functions for peaks annotation  
10/19/20: flexiable option for fragments selections (>120bp)  
05/06/20: supporting different peak calling strategies  
05/06/20: compatiable with more computational platforms  


## Basic Usage

When first starting an analysis, CUT&RUNTools requires a JSON configuration file (named `config.json`) to be written. This file specifies all that is needed to run an analysis. A sample configuration file is below. 

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
The first 10-20 lines define the software paths for various sequencing tools. For each software, please note the directory to the binary executable (usually the `bin` directory, such as `/usr/bin`, `/usr/local/bin`, or sometimes it may depend on how softwares are installed on the cluster). For more information, refer to [INSTALL.md](INSTALL.md). 

*Configuring Python and Perl*: The version of Python that is indicated in `pythonbin` field should be consistent with the version of Python used to install MACS2, and version used to install make_cut_matrix (`makecutmatrixbin`). 
Similarly for Perl (`perlbin`) and meme-chip (`memebin`). The version of Perl used to install meme should be used for `perlbin` field.

Next, check that the software prerequisites are met.
```bash
./validate.py config.json --ignore-input-output --software
```
### Start an analysis
Put sample fastq files (R1_001.fastq.gz, and R2_001.fastq.gz suffix) in a directory (in my case `/n/scratch2/qz64/Nan_18_aug23/Nan_run_19`).

Next, modify the following lines in the `config.json` file: **adapterpath** (line 8), **bt2idx** (line 17); **genome_sequence** (line 18);
section **input/output**: **fastq_directory** (line 22), **workdir** (line 23), **fastq_sequence_length** (line 24),
**organism_build** (line 25); section **cluster**: **email** (line 35).

*  `fastq_directory` is the directory containing paired-end CUT&RUN sequences (with _R1_001.fastq.gz and _R2_001.fastq.gz suffix). 
*  `workdir` is the output directory, where results are stored.
*  `organism_build` is one of supported genome assemblies: hg38, hg19, mm10, and mm9. 
*  `adapterpath` contains Illumina Truseq3-PE adapter sequences (we provide them). 
*  `genome_sequence` is the whole-genome **masked** sequence which matches with the appropriate organism build. See INSTALL.md for details.

### Create job submission scripts
```bash
./validate.py config.json
```
This script checks that your configuration file is correct and all paths are correct. You will get an empty line if the validate.py runs without errors.
```bash
./create_scripts.py config.json
```
This creates a set of Slurm job-submission scripts based on the configuration file above (named config.json) in the `workdir` directory. The scripts can be directly, easily executed.

### Easy 1-step execution

The previous step created a working directory and deposited a set of submission scripts. We can now use those scripts to start the analysis on a sample.
```bash
cd /n/scratch2/qz64/workdir
./integrated.all.steps.sh GATA1_D7_30min_chr11_R1_001.fastq.gz
```
Note the single parameter the input file. There are supposed to be two files per input sample (R1_001.fastq.gz, R2_001.fastq.gz).
**Please just use the R1 fastq as input, and do not specify both files**. The toolkit is smart enough to use the filename to automatically look for the R2 file.

What happens next is that CUT&RUNTools will sequentially run the analysis pipeline which consists of 4 major steps (for details about the steps, see option 2 below). By entering squeue -u <username> you will see 5 scripts submitted to slurm, but they will be executed sequentially one after the other via a dependency.



### Option 2: Manual four-step process

Optionally, users can run the analysis steps individually to have more control over the steps of the pipeline.

Step 1. **Read trimming, alignment.** We suppose the `workdir` is defined as `/n/scratch2/qz64/workdir`
```bash
cd /n/scratch2/qz64/workdir
sbatch ./integrated.sh CR_BCL11A_W9_r1_S17_R1_001.fastq.gz
```
The parameter is the fastq file. Even though we specify the _R1_001.fastq.gz, CUT&RUNTools actually checks that both forward and reverse fastq files are present. Always use the _R1_001 of the pair as parameter of this command.

Step 2. **BAM processing, peak calling.** It marks duplicates in bam files, and filter fragments by size. It then performs peak calling using both **MACS2 and SEACR**. Based on the results, users can choose to proceed with either MACS2 or SEACR's results.

```bash
cd aligned.aug10
sbatch ./integrated.step2.sh CR_BCL11A_W9_r1_S17_aligned_reads.bam
```

After this step, CUT&RUNTools has varied through different peak calling settings and generated multiple results for these settings in the following directories. Users need to **select only one** to go to steps 3 & 4.

| Result Directory (in `../`) | Tool | Config    | Fragments used | Use duplicates (y/n) |
| ---------------------|------|-----------|-----------|----------------------|
|macs2.narrow.aug18         | MACS2| narrowPeak| <120bp | y |
|macs2.broad.aug18          | MACS2| broadPeak | <120bp | y |
|macs2.narrow.all.frag.aug18| MACS2| narrowPeak| all    | y |
|macs2.broad.all.frag.aug18 | MACS2| broadPeak | all    | y |
|macs2.narrow.aug18.dedup   | MACS2| narrowPeak| <120bp | n |
|macs2.broad.aug18.dedup    | MACS2| broadPeak | <120bp | n |
|macs2.narrow.all.frag.aug18.dedup | MACS2| narrowPeak| all    | n |
|macs2.broad.all.frag.aug18.dedup  | MACS2| broadPeak | all    | n |
|seacr.aug12                | SEACR| stringent | <120bp | y |
|seacr.aug12.all.frag       | SEACR| stringent | all    | y |
|seacr.aug12.dedup          | SEACR| stringent | <120bp | n | 
|seacr.aug12.all.frag.dedup | SEACR| stringent | all    | n |

* Which directory to use: if **TF CUT&RUN**, I prefer **macs2.narrow.aug18** or **macs2.narrow.aug18.dedup**. If **histone CUT&RUN**, use **macs2.broad.all.frag.aug18**. If **SEACR**, use **seacr.aug12.all.frag** (histone) or **seacr.aug12** (TF) and use the **stringent** peaks within each folder.

* **Large fragment (>120bp) peaks**:
    * use the peak calling directories that end in `*.all.frag` (e.g. macs2.broad.all.frag.aug18). 

Step 3. **Motif finding.** CUT&RUNTools uses MEME-chip for de novo motif finding on sequences surrounding the peak summits.
```bash
#Use macs2.narrow.aug18 or any of the peak calling result directory in the above table
cd ../macs2.narrow.aug18
#For narrow setting, the peak file ends in .narrowPeak
#For broad setting, the peak file ends in .broadPeak
#For seacr setting, the peak file ends in .stringent.sort.bed
#Use the right peak file accordingly
sbatch ./integrate.motif.find.sh CR_BCL11A_W9_r1_S17_aligned_reads_peaks.narrowPeak
```
Similar procedure applies in other peak calling directories for broadPeaks, all fragment results, or SEACR peaks.
```bash
#SEACR
cd ../seacr.aug12.all.frag
sbatch ./integrate.motif.find.sh GATA1_D7_30min_chr11_aligned_reads_treat.stringent.sort.bed
```

Step 4. **Motif footprinting.**
```bash
cd ../macs2.narrow.aug18
#supports file type narrowPeak, broadPeak, or stringent.sort.bed (SEACR)
sbatch ./integrate.footprinting.sh CR_BCL11A_W9_r1_S17_aligned_reads_peaks.narrowPeak
```
Beautiful footprinting figures will be located in the directory `fimo.result`. Footprinting figures are created for every motif found by MEME-chip, but only the right motif (associated with TF) will have a proper looking shape. Users can scan through all the motifs' footprints.


### Outputs

Please see [USAGE.md](USAGE.md) for more details for the output files.

Briefly, CUT&RUNTools will generate the following directory structure in the output folder.

```
+ macs2.narrow.aug18
  + random.10000
    + MEME_GATA1_HDP2_30min_S13_aligned_reads_shuf
      ++ summary.tsv
  + fimo.result
    + GATA1_HDP2_30min_S13_aligned_reads_peaks
      + fimo2.DREME-10.CCWATCAG
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
      + fimo2.DREME-11.CTCCWCCC
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
      + fimo2.DREME-12.CACGTG
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
      + fimo2.DREME-13.CAKCTGB
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
      + fimo2.DREME-14.CWGATA
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
      + fimo2.DREME-1.HGATAA
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
      + fimo2.DREME-20.CCCTYCC
        ++ fimo.png
        ++ fimo.logratio.txt
        ++ fimo.bed
```

Files are denoted by "++", and directories are denoted by "+".
Here the sample name is called GATA1_HDP2_30min_S13.
The summary.tsv lists all the motif found by motif searching analysis. The files within "fimo2..." folder shows the motif footprinting figure (.png), the motif sites for the motif (.bed), and the logratio binding scores at these sites (.logratio.txt).

### Single locus cut profile

CUT&RUNTools allows users to obtain a single nucleotide resolution cut profile for a region of interest.

For more details, please see [USAGE.md](USAGE.md).

### Footprinting for a user-specified motif

CUT&RUNTools can run the motif scanning and motif footprinting step on a user-specified motif, such as a motif from the public JASPAR database. The motif should be in the MEME format. 

The script that you will want to look at is `generate.footprinting.factor.specific.centipede.py` that is in the directory `macs2.narrow.aug18` or `macs2.narrow.aug18.dedup`.

To do this analysis:

```bash
pwd
/n/scratch2/qz64/workdir
cd macs2.narrow.aug18
./generate.footprinting.factor.specific.centipede.py -b MA.00001.agataa.meme -p 0.001 -n GATA1
```

Option **-b** specifies the MEME file. Option **-p** is the motif scanning p-value (recommended 0.0005, but for this example we will use 0.001 since GATA1 motif is quite short). Option **-n** is the name you give it. 

The script will generate a custom script for GATA1 motif, called `integrate.footprinting.GATA1.centipede.sh`. With this script, then you can run it on a narrowPeak file as follows:

```bash
sbatch ./integrate.footprinting.GATA1.centipede.sh GATA1_HDP2_30min_S13_aligned_reads_peaks.narrowPeak
```

The output will be in the `fimo.GATA1.result` folder, with the cut frequency estimated, and log odds binding score estimated.


## Want to try?

Download a small CUT&RUN experiment [qzhudfci/datasets/src](https://bitbucket.org/qzhudfci/datasets/src). This is GATA1 chr11 only.
Follow along the instructions [INSTALL.md](INSTALL.md), [USAGE.md](USAGE.md).