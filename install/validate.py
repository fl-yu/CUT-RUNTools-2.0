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

# check python module
def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
		return False
		print module_name, "is not found in python"
    else:
        return True

# ./validate.py config.json --ignore-input-output --software

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
		for filename in os.listdir(config["input_output"]["fastq_directory"]):
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
	if not check_program_exists(config["software_config"]["bowtie2bin"], "bowtie2"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["memebin"], "meme-chip"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["perlbin"], "perl"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["javabin"], "java"):
		sys.exit(1)

	if not check_program_exists(config["software_config"]["bedopsbin"], "bedops"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["bedopsbin"], "gff2bed"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["bedopsbin"], "sort-bed"):
		sys.exit(1)

	if not check_program_exists(config["software_config"]["samtoolsbin"], "samtools"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["trimmomaticbin"], config["software_config"]["trimmomaticjarfile"]):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["macs2bin"], "macs2"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["picardbin"], config["software_config"]["picardjarfile"]):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["kseqbin"], "kseq_test"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["bedtoolsbin"], "bedtools"):
		sys.exit(1)
#	if not check_program_exists(config["software_config"]["makecutmatrixbin"], "make_cut_matrix"):
#		sys.exit(1)
	if not check_program_exists(config["software_config"]["Rscriptbin"], "Rscript"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["pythonbin"], "python"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["adapterpath"], "Truseq3.PE.fa"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["extratoolsbin"], "SEACR_1.1.sh"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["extratoolsbin"], "SEACR_1.1.R"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["extratoolsbin"], "change.bdg.py"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["extratoolsbin"], "get_summits_seacr.py"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["extratoolsbin"], "get_summits_broadPeak.py"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["path_parallel"], "parallel"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["path_deeptools"], "deeptools"):
		sys.exit(1)
	if not check_program_exists(config["software_config"]["path_tabix"], "tabix"):
		sys.exit(1)

	#valide python packages
	if not module_exists("umap-learn"):
    		sys.exit(1)
	if not module_exists("leidenalg"):
    		sys.exit(1)
	if not module_exists("igraph"):
    		sys.exit(1)

	#check genome sequence exists, and check bowtie2 indices
	if not args["ignore_input_output"]:
		if not check_program_exists(os.path.dirname(config["software_config"]["genome_sequence"]), "%s.chrom.sizes" % 
		config["input_output"]["genome"]):
			sys.exit(1)

		org = config["input_output"]["genome"]
		if org=="hg38":
			org = "GRCh38"
		for ff in ["%s.1.bt2" % org, "%s.2.bt2" % org, "%s.3.bt2" % org, \
		"%s.4.bt2" % org, "%s.rev.1.bt2" % org, "%s.rev.2.bt2" % org]:	
			if not check_program_exists(config["software_config"]["bt2idx"], ff):
				sys.exit(1)
	print "abc3"
	#test software works
	if args["software"]:
		print "======================Testing Rscript...======================"
		os.system("%s/Rscript --version" % config["software_config"]["Rscriptbin"])
		print "======================Testing python...======================"
		os.system("%s/python --version" % config["software_config"]["pythonbin"])
		print "======================Testing perl...======================"
		os.system("%s/perl -version" % config["software_config"]["perlbin"])
		print "======================Testing java...======================"
		os.system("%s/java -version" % config["software_config"]["javabin"])
		print "======================Testing trimmomatic...======================"
		os.system("%s/java -jar %s/%s -version" % (config["software_config"]["javabin"], config["software_config"]["trimmomaticbin"], config["software_config"]["trimmomaticjarfile"]))
		print "======================Testing bowtie2...======================"
		os.system("%s/bowtie2 --version" % config["software_config"]["bowtie2bin"])
		print "======================Testing samtools...======================"
		os.system("%s/samtools --version" % config["software_config"]["samtoolsbin"])
		print "======================Testing picard...======================"
		os.system("%s/java -jar %s/%s -h" % (config["software_config"]["javabin"], config["software_config"]["picardbin"], config["software_config"]["picardjarfile"]))
		print "======================Testing macs2...======================"
		p_pythonbin = config["software_config"]["pythonbin"]
		p_pythonbin_suffix = p_pythonbin.rstrip("/").rstrip("/bin")
		p_pythoninclude = p_pythonbin_suffix + "/include"
		p_pythonlib = p_pythonbin_suffix + "/lib"

		#cmd1 = "pythonpath=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythonpath" % (p_pythonbin, p_pythonbin)
		#cmd2 = "pythoninclude=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythoninclude" % (p_pythoninclude, p_pythoninclude)
		cmd3 = "ldlibrary=`echo $LD_LIBRARY_PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset LD_LIBRARY_PATH && export LD_LIBRARY_PATH=%s:$ldlibrary" % (p_pythonlib, p_pythonlib)
		cmd4 = "pythonlib=`echo $PYTHONPATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PYTHONPATH && export PYTHONPATH=%s:$pythonlib" % (config["software_config"]["macs2pythonlib"], config["software_config"]["macs2pythonlib"])
		cmd5 = "%s/macs2 --version" % config["software_config"]["macs2bin"]
		all_cmd = " && ".join([cmd3, cmd4, cmd5])
		#print all_cmd
		os.system(all_cmd)
		#os.system("pythonlib=`echo $PYTHONPATH | tr : \\n | grep -v %s | paste -s -d:` && unset $PYTHONPATH && export PYTHONPATH=%s:$pythonlib && %s/macs2 --version" % (config["macs2pythonlib"], config["macs2pythonlib"], config["macs2bin"]))
		print "======================Testing kseq...======================"
		os.system("%s/kseq_test --help" % config["software_config"]["kseqbin"])
		print "======================Testing meme...======================"
		os.system("%s/meme -version" % config["software_config"]["memebin"])
		print "======================Testing meme-chip...======================"
		os.system("%s/meme-chip -version" % (config["software_config"]["memebin"]))
		print "======================Testing bedops...======================"
		os.system("%s/bedops --version" % config["software_config"]["bedopsbin"])
		print "======================Testing bedtools...======================"
		os.system("%s/bedtools --version" % config["software_config"]["bedtoolsbin"])
		print "======================Testing make_cut_matrix...======================"
		#cmd1 = "pythonpath=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythonpath" % (p_pythonbin, p_pythonbin)
		#cmd2 = "pythoninclude=`echo $PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PATH && export PATH=%s:$pythoninclude" % (p_pythoninclude, p_pythoninclude)
		cmd3 = "ldlibrary=`echo $LD_LIBRARY_PATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset LD_LIBRARY_PATH && export LD_LIBRARY_PATH=%s:$ldlibrary" % (p_pythonlib, p_pythonlib)
		#cmd4 = "pythonlib=`echo $PYTHONPATH | tr \":\" \"\\n\" | grep -v %s | paste -s -d:` && unset PYTHONPATH && export PYTHONPATH=%s:$pythonlib" % (config["macs2pythonlib"], config["macs2pythonlib"])
		cmd5 = "%s/make_cut_matrix --version" % config["software_config"]["makecutmatrixbin"]
		all_cmd = " && ".join([cmd3, cmd5])
		#print all_cmd
		os.system(all_cmd)
		#os.system("%s/make_cut_matrix --version" % config["makecutmatrixbin"])

	print("[INFO] All the software seems OK!")

