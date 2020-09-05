#!/usr/bin/python
# Copyright (C) 2020 Fulong Yu
#
# CUT&RUNTools 2.0 is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools 2.0 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools 2.0 and is found in LICENSE.md.

import shutil
import sys
import os
import re
import json
import argparse

def check_program_exists(path, program):
	if not os.path.isfile(path + "/" + program):
		print program, "is not found in", path
		return False
	return True
	

if __name__=="__main__":
	parser = argparse.ArgumentParser(prog="validate")
	parser.add_argument("--software", action="store_true")
	parser.add_argument("--ignore-input-output", action="store_true")
	parser.add_argument("config", type=file)

	args = vars(parser.parse_args(sys.argv[1:]))
	#f = open(args["config"]) #config.json
	f = args["config"]
	config = json.load(f)
	f.close()

	#check input and output portion of config
	if not args["ignore_input_output"]:
		#validate the input files
		flist = set([])
		for filename in os.listdir(config["input/output"]["fastq_directory"]):
			if filename.endswith("fastq.gz"):
				flist.add(filename)
		for x in flist:
			m1 = re.match("(.*)_R1_001.fastq.gz", x)
			m2 = re.match("(.*)_R2_001.fastq.gz", x)
			if m1 is None and m2 is None:
				print "Error:", x, " no pattern of _R1_001.fastq.gz, or _R2_001.fastq.gz detected"
				sys.exit(1)
			if m1 is not None:
				x2 = m1.group(1) + "_R2_001.fastq.gz"
				if not x2 in flist:
					print x, "does not have a corresponding _R2_001.fastq.gz file"
					sys.exit(1)
			if m2 is not None:
				x2 = m2.group(1) + "_R1_001.fastq.gz"
				if not x2 in flist:
					print x, "does not have a corresponding _R1_001.fastq.gz file"
					sys.exit(1)

	#valide software path
	if not check_program_exists(config["bowtie2bin"], "bowtie2"):
		sys.exit(1)
	if not check_program_exists(config["memebin"], "meme-chip"):
		sys.exit(1)
	if not check_program_exists(config["perlbin"], "perl"):
		sys.exit(1)
	if not check_program_exists(config["javabin"], "java"):
		sys.exit(1)

	if not check_program_exists(config["bedopsbin"], "bedops"):
		sys.exit(1)
	if not check_program_exists(config["bedopsbin"], "gff2bed"):
		sys.exit(1)
	if not check_program_exists(config["bedopsbin"], "sort-bed"):
		sys.exit(1)

	if not check_program_exists(config["samtoolsbin"], "samtools"):
		sys.exit(1)
	if not check_program_exists(config["trimmomaticbin"], config["trimmomaticjarfile"]):
		sys.exit(1)
	if not check_program_exists(config["macs2bin"], "macs2"):
		sys.exit(1)
	if not check_program_exists(config["picardbin"], config["picardjarfile"]):
		sys.exit(1)
	if not check_program_exists(config["kseqbin"], "kseq_test"):
		sys.exit(1)
	if not check_program_exists(config["bedtoolsbin"], "bedtools"):
		sys.exit(1)
	if not check_program_exists(config["makecutmatrixbin"], "make_cut_matrix"):
		sys.exit(1)
	if not check_program_exists(config["Rscriptbin"], "Rscript"):
		sys.exit(1)
	if not check_program_exists(config["pythonbin"], "python"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "bedGraphToBigWig"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "fetchChromSizes"):
		sys.exit(1)
	if not check_program_exists(config["extrasettings"], "filter_below.awk"):
		sys.exit(1)
	if not check_program_exists(config["adapterpath"], "Truseq3.PE.fa"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "SEACR_1.1.sh"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "SEACR_1.1.R"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "change.bdg.py"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "get_summits_seacr.py"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "get_summits_broadPeak.py"):
		sys.exit(1)

	#check genome sequence exists, and check bowtie2 indices
	if not args["ignore_input_output"]:
		if not check_program_exists(os.path.dirname(config["genome_sequence"]), "%s.chrom.sizes" % 
		config["input/output"]["organism_build"]):
			sys.exit(1)

		org = config["input/output"]["organism_build"]
		if org=="hg38":
			org = "GRCh38"
		for ff in ["%s.1.bt2" % org, "%s.2.bt2" % org, "%s.3.bt2" % org, \
		"%s.4.bt2" % org, "%s.rev.1.bt2" % org, "%s.rev.2.bt2" % org]:	
			if not check_program_exists(config["bt2idx"], ff):
				sys.exit(1)
	
	#test software works
	if args["software"]:
		print "======================Testing Rscript...======================"
		os.system("%s/Rscript --version" % config["Rscriptbin"])
		print "======================Testing python...======================"
		os.system("%s/python --version" % config["pythonbin"])
		print "======================Testing perl...======================"
		os.system("%s/perl -version" % config["perlbin"])
		print "======================Testing java...======================"
		os.system("%s/java -version" % config["javabin"])
		print "======================Testing trimmomatic...======================"
		os.system("%s/java -jar %s/%s -version" % (config["javabin"], config["trimmomaticbin"], config["trimmomaticjarfile"]))
		print "======================Testing bowtie2...======================"
		os.system("%s/bowtie2 --version" % config["bowtie2bin"])
		print "======================Testing samtools...======================"
		os.system("%s/samtools --version" % config["samtoolsbin"])
		print "======================Testing picard...======================"
		os.system("%s/java -jar %s/%s -h" % (config["javabin"], config["picardbin"], config["picardjarfile"]))
		print "======================Testing macs2...======================"
		p_pythonbin = config["pythonbin"]
		p_pythonbin_suffix = p_pythonbin.rstrip("/").rstrip("/bin")
		p_pythoninclude = p_pythonbin_suffix + "/include"
		p_pythonlib = p_pythonbin_suffix + "/lib"

		#cmd1 = "pythonpath=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythonpath" % (p_pythonbin, p_pythonbin)
		#cmd2 = "pythoninclude=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythoninclude" % (p_pythoninclude, p_pythoninclude)
		cmd3 = "ldlibrary=`echo $LD_LIBRARY_PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset LD_LIBRARY_PATH && export LD_LIBRARY_PATH=%s:$ldlibrary" % (p_pythonlib, p_pythonlib)
		cmd4 = "pythonlib=`echo $PYTHONPATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PYTHONPATH && export PYTHONPATH=%s:$pythonlib" % (config["macs2pythonlib"], config["macs2pythonlib"])
		cmd5 = "%s/macs2 --version" % config["macs2bin"]
		all_cmd = " && ".join([cmd3, cmd4, cmd5])
		#print all_cmd
		os.system(all_cmd)
		#os.system("pythonlib=`echo $PYTHONPATH | tr : \\n | grep -v %s | paste -s -d:` && unset $PYTHONPATH && export PYTHONPATH=%s:$pythonlib && %s/macs2 --version" % (config["macs2pythonlib"], config["macs2pythonlib"], config["macs2bin"]))
		print "======================Testing kseq...======================"
		os.system("%s/kseq_test --help" % config["kseqbin"])
		print "======================Testing meme...======================"
		os.system("%s/meme -version" % config["memebin"])
		print "======================Testing meme-chip...======================"
		os.system("%s/meme-chip -version" % (config["memebin"]))
		print "======================Testing bedops...======================"
		os.system("%s/bedops --version" % config["bedopsbin"])
		print "======================Testing bedtools...======================"
		os.system("%s/bedtools --version" % config["bedtoolsbin"])
		print "======================Testing make_cut_matrix...======================"
		#cmd1 = "pythonpath=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythonpath" % (p_pythonbin, p_pythonbin)
		#cmd2 = "pythoninclude=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythoninclude" % (p_pythoninclude, p_pythoninclude)
		cmd3 = "ldlibrary=`echo $LD_LIBRARY_PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset LD_LIBRARY_PATH && export LD_LIBRARY_PATH=%s:$ldlibrary" % (p_pythonlib, p_pythonlib)
		#cmd4 = "pythonlib=`echo $PYTHONPATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PYTHONPATH && export PYTHONPATH=%s:$pythonlib" % (config["macs2pythonlib"], config["macs2pythonlib"])
		cmd5 = "%s/make_cut_matrix --version" % config["makecutmatrixbin"]
		all_cmd = " && ".join([cmd3, cmd5])
		#print all_cmd
		os.system(all_cmd)
		#os.system("%s/make_cut_matrix --version" % config["makecutmatrixbin"])

