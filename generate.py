#!/usr/bin/python
import sys
import math
import os
import numpy as np
import scipy
import json
from collections import OrderedDict

if __name__=="__main__":
	config = OrderedDict([
	("Rscriptbin", "/n/app/R/3.3.3/bin"),
	("pythonbin", "/n/app/python/2.7.12/bin"),
	("trimmomaticbin", "/n/app/trimmomatic/0.36/bin"),
	("trimmomaticjarfile", "trimmomatic-0.36.jar"),
	("bowtie2bin", "/n/app/bowtie2/2.2.9/bin"),
	("samtoolsbin", "/n/app/samtools/1.3.1/bin"),
	("adapterpath", "/home/qz64"),
	("picardbin", "/n/app/picard/2.8.0/bin"), 
	("picardjarfile", "picard-2.8.0.jar"),
	("macs2bin", "/n/app/macs2/2.1.1.20160309/bin"),
	("kseqbin", "/home/qz64"),
	("memebin", "/home/qz64/meme/bin"),
	("bedopsbin", "/home/qz64/bin"),
	("bedtoolsbin", "/n/app/bedtools/2.26.0/bin"),
	("makecutmatrixbin", "/home/qz64/.local/bin"),
	("bt2idx", "/n/groups/shared_databases/bowtie2_indexes"),
	("genome_sequence", "/home/qz64/chrom.hg19/hg19.fa"),
	("extratoolsbin", "/home/qz64"),
	("extrasettings", "/home/qz64"),
	("input/output", OrderedDict([
		("fastq_directory", "/n/scratch2/qz64/Nan_18_aug23/Nan_run_19"),
		("workdir", "/n/scratch2/qz64/workdir"),
		("fastq_sequence_length", 40),
		("organism_build", "hg19"),
	])),
	("motif_finding", OrderedDict([
		("num_bp_from_summit", 150),
		("num_peaks", 5000),
		("total_peaks", 15000),
		("motif_scanning_pval", 0.0005),
		("num_motifs", 20),
	])),
	("cluster", OrderedDict([
		("email", "bernardzhu@gmail.com"),
		("step_alignment", OrderedDict([	
			("queue", "short"),
			("memory", 32000),
			("time_limit", "0-12:00"),
		])), 
		("step_process_bam", OrderedDict([
			("queue", "short"),
			("memory", 32000),
			("time_limit", "0-12:00"),
		])), 
		("step_motif_find", OrderedDict([
			("queue", "short"),
			("memory", 32000),
			("time_limit", "0-12:00"),
		])),
		("step_footprinting", OrderedDict([
			("queue", "short"),
			("memory", 32000),
			("time_limit", "0-12:00"),
		])), 
	])),
	])
	fw = open("config.json", "w")
	json.dump(config, fw, indent=4)
	fw.close()
