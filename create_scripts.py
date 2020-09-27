#!/usr/bin/python
# Copyright (C) 2019 Qian Zhu
#
# CUT&RUNTools is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.


import shutil
import sys
import os
import re
import json

def make_executable(filename):
	st = os.stat(filename)
	os.chmod(filename, st.st_mode | 0o111)
	
def generate_integrated_sh(config, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
		
	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%d                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_alignment"]["time_limit"], config["cluster"]["step_alignment"]["queue"], 
	config["cluster"]["step_alignment"]["memory"], config["cluster"]["email"])

	path="""trimmomaticbin=%s
trimmomaticjarfile=%s
adapterpath=%s
bowtie2bin=%s
samtoolsbin=%s
javabin=%s

bt2idx=%s
kseqbin=%s

infile=$1
#expand the path of infile
relinfile=`realpath -s $infile`
dirname=`dirname $relinfile`
base=`basename $infile _R1_001.fastq.gz`
>&2 echo "Input file is $relinfile"
>&2 date

#cd to current directory
cd $dirname
workdir=`pwd`

len=`cat length`
trimdir=$workdir/trimmed
trimdir2=$workdir/trimmed3
logdir=$workdir/logs
aligndir=$workdir/aligned.aug10

for d in $trimdir $trimdir2 $logdir $aligndir; do
if [ ! -d $d ]; then
mkdir $d
fi
done

""" % (config["trimmomaticbin"], config["trimmomaticjarfile"], config["adapterpath"], 
	config["bowtie2bin"], config["samtoolsbin"], config["javabin"],
	config["bt2idx"], config["kseqbin"])

	bowtie2_org = config["input/output"]["organism_build"]
	if config["input/output"]["organism_build"]=="hg38":
		bowtie2_org = "GRCh38"

	scripts="""
#trimming paired-end
#good version
>&2 echo "Trimming file $base ..."
>&2 date
$javabin/java -jar $trimmomaticbin/$trimmomaticjarfile PE -threads 1 -phred33 $dirname/"$base"_R1_001.fastq.gz $dirname/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

>&2 echo "Second stage trimming $base ..."
>&2 date
$kseqbin/kseq_test $trimdir/"$base"_1.paired.fastq.gz $len $trimdir2/"$base"_1.paired.fastq.gz
$kseqbin/kseq_test $trimdir/"$base"_2.paired.fastq.gz $len $trimdir2/"$base"_2.paired.fastq.gz

>&2 echo "Aligning file $base ..."
>&2 date
($bowtie2bin/bowtie2 -p 2 --dovetail --phred33 -x $bt2idx/%s -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/"$base"_aligned_reads.bam

>&2 echo "Finished"
>&2 date
""" % bowtie2_org
	outp.write(header + "\n")
	outp.write(path + "\n")
	outp.write(scripts + "\n")

	if output is not None:
		outp.close()

def generate_integrated_step2_sh(config, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%d                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_process_bam"]["time_limit"], config["cluster"]["step_process_bam"]["queue"], 
	config["cluster"]["step_process_bam"]["memory"], config["cluster"]["email"])

	p_pythonbase = config["pythonbin"].rstrip("/").rstrip("/bin")

	path = """Rscriptbin=%s
pythonbin=%s
bedopsbin=%s
picardbin=%s
picardjarfile=%s
samtoolsbin=%s
macs2bin=%s
javabin=%s
extratoolsbin=%s
extrasettings=%s
chromsizedir=`dirname %s`
macs2pythonlib=%s

pythonlib=`echo $PYTHONPATH | tr : "\\n" | grep -v $macs2pythonlib | paste -s -d:`
unset PYTHONPATH
export PYTHONPATH=$macs2pythonlib:$pythonlib

pythonldlibrary=%s
ldlibrary=`echo $LD_LIBRARY_PATH | tr : "\\n" | grep -v $pythonldlibrary | paste -s -d:`
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$pythonldlibrary:$ldlibrary

""" % (config["Rscriptbin"], config["pythonbin"], config["bedopsbin"], 
	config["picardbin"], config["picardjarfile"], config["samtoolsbin"], 
	config["macs2bin"], config["javabin"],
	config["extratoolsbin"], config["extrasettings"], 
	config["genome_sequence"], config["macs2pythonlib"], p_pythonbase + "/lib")

	scripts =""">&2 echo "Input parameters are: $1"
>&2 date

#expand the path of $1
relinfile=`realpath -s $1`
dirname=`dirname $relinfile`
base=`basename $1 .bam`

#cd to current directory (aligned.aug10)
cd $dirname

workdir=`pwd`
logdir=$workdir/logs

for d in $logdir sorted dup.marked dedup; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "Filtering unmapped fragments... ""$base".bam
>&2 date
$samtoolsbin/samtools view -bh -f 3 -F 4 -F 8 $dirname/"$base".bam > sorted/"$base".step1.bam

>&2 echo "Sorting BAM... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile SortSam \
INPUT=sorted/"$base".step1.bam OUTPUT=sorted/"$base".bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
rm -rf sorted/"$base".step1.bam

>&2 echo "Marking duplicates... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dup.marked/"$base".bam VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=metrics."$base".txt

>&2 echo "Removing duplicates... ""$base".bam
>&2 date
$samtoolsbin/samtools view -bh -F 1024 dup.marked/"$base".bam > dedup/"$base".bam

for d in dup.marked.120bp dedup.120bp; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "Filtering to <120bp... ""$base".bam
>&2 date
$samtoolsbin/samtools view -h dup.marked/"$base".bam |LC_ALL=C awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dup.marked.120bp/"$base".bam
$samtoolsbin/samtools view -h dedup/"$base".bam |LC_ALL=C awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dedup.120bp/"$base".bam

>&2 echo "Creating bam index files... ""$base".bam
>&2 date
$samtoolsbin/samtools index sorted/"$base".bam
$samtoolsbin/samtools index dup.marked/"$base".bam
$samtoolsbin/samtools index dedup/"$base".bam
$samtoolsbin/samtools index dup.marked.120bp/"$base".bam
$samtoolsbin/samtools index dedup.120bp/"$base".bam

>&2 echo "Peak calling using MACS2... ""$base".bam
>&2 echo "Logs are stored in $logdir"
>&2 date
bam_file=dup.marked.120bp/"$base".bam
dir=`dirname $bam_file`
base_file=`basename $bam_file .bam`
"""

	macs2_org = "hs"
	if config["input/output"]["organism_build"]=="hg19" or \
	config["input/output"]["organism_build"]=="hg38":
		macs2_org = "hs"
	elif config["input/output"]["organism_build"]=="mm10" or \
	config["input/output"]["organism_build"]=="mm9":
		macs2_org = "mm"

	macs_script = """
outdir=$workdir/../macs2.narrow.aug18 #for macs2
outdir2=$workdir/../macs2.narrow.aug18.dedup #for macs2 dedup version

outdirbroad=$workdir/../macs2.broad.aug18 #for macs2
outdirbroad2=$workdir/../macs2.broad.aug18.dedup #for macs2 dedup version

outdirseac=$workdir/../seacr.aug12 #for seacr
outdirseac2=$workdir/../seacr.aug12.dedup #for seacr dedup version

for d in $outdir $outdir2 $outdirbroad $outdirbroad2 $outdirseac $outdirseac2; do
if [ ! -d $d ]; then
mkdir $d
fi
done

$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdir2 -q 0.01 -B --SPMR 2> $logdir/"$base_file".dedup.macs2

#broad peak calls
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> $logdir/"$base_file".broad.all.frag.macs2
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirbroad2 --broad --broad-cutoff 0.1 -B --SPMR 2> $logdir/"$base_file".broad.all.frag.dedup.macs2
$pythonbin/python $extratoolsbin/get_summits_broadPeak.py $outdirbroad/"$base_file"_peaks.broadPeak|$bedopsbin/sort-bed - > $outdirbroad/"$base_file"_summits.bed
$pythonbin/python $extratoolsbin/get_summits_broadPeak.py $outdirbroad2/"$base_file"_peaks.broadPeak|$bedopsbin/sort-bed - > $outdirbroad2/"$base_file"_summits.bed

#SEACR peak calls
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirseac -q 0.01 -B --keep-dup all
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirseac2 -q 0.01 -B
$pythonbin/python $extratoolsbin/change.bdg.py $outdirseac/"$base_file"_treat_pileup.bdg > $outdirseac/"$base_file"_treat_integer.bdg
$pythonbin/python $extratoolsbin/change.bdg.py $outdirseac2/"$base_file"_treat_pileup.bdg > $outdirseac2/"$base_file"_treat_integer.bdg
$extratoolsbin/SEACR_1.1.sh $outdirseac/"$base_file"_treat_integer.bdg 0.01 non stringent $outdirseac/"$base_file"_treat $Rscriptbin
$extratoolsbin/SEACR_1.1.sh $outdirseac2/"$base_file"_treat_integer.bdg 0.01 non stringent $outdirseac2/"$base_file"_treat $Rscriptbin
$bedopsbin/sort-bed $outdirseac/"$base_file"_treat.stringent.bed > $outdirseac/"$base_file"_treat.stringent.sort.bed
$bedopsbin/sort-bed $outdirseac2/"$base_file"_treat.stringent.bed > $outdirseac2/"$base_file"_treat.stringent.sort.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac/"$base_file"_treat.stringent.bed|$bedopsbin/sort-bed - > $outdirseac/"$base_file"_treat.stringent.sort.summits.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac2/"$base_file"_treat.stringent.bed|$bedopsbin/sort-bed - > $outdirseac2/"$base_file"_treat.stringent.sort.summits.bed
for i in _summits.bed _peaks.xls _peaks.narrowPeak _control_lambda.bdg _treat_pileup.bdg; do 
rm -rf $outdirseac/"$base_file"$i
rm -rf $outdirseac2/"$base_file"$i
done

#SEACR relaxed peak calls
$extratoolsbin/SEACR_1.1.sh $outdirseac/"$base_file"_treat_integer.bdg 0.01 non relaxed $outdirseac/"$base_file"_treat $Rscriptbin
$extratoolsbin/SEACR_1.1.sh $outdirseac2/"$base_file"_treat_integer.bdg 0.01 non relaxed $outdirseac2/"$base_file"_treat $Rscriptbin
$bedopsbin/sort-bed $outdirseac/"$base_file"_treat.relaxed.bed > $outdirseac/"$base_file"_treat.relaxed.sort.bed
$bedopsbin/sort-bed $outdirseac2/"$base_file"_treat.relaxed.bed > $outdirseac2/"$base_file"_treat.relaxed.sort.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac/"$base_file"_treat.relaxed.bed|$bedopsbin/sort-bed - > $outdirseac/"$base_file"_treat.relaxed.sort.summits.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac2/"$base_file"_treat.relaxed.bed|$bedopsbin/sort-bed - > $outdirseac2/"$base_file"_treat.relaxed.sort.summits.bed
""" % (macs2_org, macs2_org, macs2_org, macs2_org, macs2_org, macs2_org)

	scripts2 = """
cur=`pwd`
>&2 echo "Converting bedgraph to bigwig... ""$base".bam
>&2 date
cd $outdir
LC_ALL=C sort -k1,1 -k2,2n $outdir/"$base_file"_treat_pileup.bdg > $outdir/"$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$base_file".sort.bdg $chromsizedir/%s.chrom.sizes $outdir/"$base_file".sorted.bw
rm -rf "$base_file".sort.bdg

cd $outdir2
LC_ALL=C sort -k1,1 -k2,2n $outdir2/"$base_file"_treat_pileup.bdg > $outdir2/"$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig $outdir2/"$base_file".sort.bdg $chromsizedir/%s.chrom.sizes $outdir2/"$base_file".sorted.bw
rm -rf "$base_file".sort.bdg
""" % (config["input/output"]["organism_build"], config["input/output"]["organism_build"])


	scriptallfrag = """#====================================================================================================================================

#all fragments
cd $cur
bam_file=dup.marked/"$base".bam
dir=`dirname $bam_file`
base_file=`basename $bam_file .bam`

outdir=$workdir/../macs2.narrow.all.frag.aug18 #for macs2
outdir2=$workdir/../macs2.narrow.all.frag.aug18.dedup #for macs2 dedup version

outdirbroad=$workdir/../macs2.broad.all.frag.aug18 #for macs2
outdirbroad2=$workdir/../macs2.broad.all.frag.aug18.dedup #for macs2 dedup version

#SEACR peak calling
outdirseac=$workdir/../seacr.aug12.all.frag #for seacr
outdirseac2=$workdir/../seacr.aug12.all.frag.dedup #for seacr dedup version

for d in $outdir $outdir2 $outdirbroad $outdirbroad2 $outdirseac $outdirseac2; do
if [ ! -d $d ]; then
mkdir $d
fi
done

$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdir2 -q 0.01 -B --SPMR 2> $logdir/"$base_file".dedup.macs2

#broad peak calls
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --keep-dup all 2> $logdir/"$base_file".broad.all.frag.macs2
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirbroad2 --broad --broad-cutoff 0.1 -B 2> $logdir/"$base_file".broad.all.frag.dedup.macs2
$pythonbin/python $extratoolsbin/get_summits_broadPeak.py $outdirbroad/"$base_file"_peaks.broadPeak|$bedopsbin/sort-bed - > $outdirbroad/"$base_file"_summits.bed
$pythonbin/python $extratoolsbin/get_summits_broadPeak.py $outdirbroad2/"$base_file"_peaks.broadPeak|$bedopsbin/sort-bed - > $outdirbroad2/"$base_file"_summits.bed

#SEACR peak calling
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirseac -q 0.01 -B --keep-dup all
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdirseac2 -q 0.01 -B 
$pythonbin/python $extratoolsbin/change.bdg.py $outdirseac/"$base_file"_treat_pileup.bdg > $outdirseac/"$base_file"_treat_integer.bdg
$pythonbin/python $extratoolsbin/change.bdg.py $outdirseac2/"$base_file"_treat_pileup.bdg > $outdirseac2/"$base_file"_treat_integer.bdg
$extratoolsbin/SEACR_1.1.sh $outdirseac/"$base_file"_treat_integer.bdg 0.01 non stringent $outdirseac/"$base_file"_treat $Rscriptbin
$extratoolsbin/SEACR_1.1.sh $outdirseac2/"$base_file"_treat_integer.bdg 0.01 non stringent $outdirseac2/"$base_file"_treat $Rscriptbin
$bedopsbin/sort-bed $outdirseac/"$base_file"_treat.stringent.bed > $outdirseac/"$base_file"_treat.stringent.sort.bed
$bedopsbin/sort-bed $outdirseac2/"$base_file"_treat.stringent.bed > $outdirseac2/"$base_file"_treat.stringent.sort.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac/"$base_file"_treat.stringent.bed|$bedopsbin/sort-bed - > $outdirseac/"$base_file"_treat.stringent.sort.summits.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac2/"$base_file"_treat.stringent.bed|$bedopsbin/sort-bed - > $outdirseac2/"$base_file"_treat.stringent.sort.summits.bed
for i in _summits.bed _peaks.xls _peaks.narrowPeak _control_lambda.bdg _treat_pileup.bdg; do 
rm -rf $outdirseac/"$base_file"$i
rm -rf $outdirseac2/"$base_file"$i
done

$extratoolsbin/SEACR_1.1.sh $outdirseac/"$base_file"_treat_integer.bdg 0.01 non relaxed $outdirseac/"$base_file"_treat $Rscriptbin
$extratoolsbin/SEACR_1.1.sh $outdirseac2/"$base_file"_treat_integer.bdg 0.01 non relaxed $outdirseac2/"$base_file"_treat $Rscriptbin
$bedopsbin/sort-bed $outdirseac/"$base_file"_treat.relaxed.bed > $outdirseac/"$base_file"_treat.relaxed.sort.bed
$bedopsbin/sort-bed $outdirseac2/"$base_file"_treat.relaxed.bed > $outdirseac2/"$base_file"_treat.relaxed.sort.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac/"$base_file"_treat.relaxed.bed|$bedopsbin/sort-bed - > $outdirseac/"$base_file"_treat.relaxed.sort.summits.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac2/"$base_file"_treat.relaxed.bed|$bedopsbin/sort-bed - > $outdirseac2/"$base_file"_treat.relaxed.sort.summits.bed
""" % (macs2_org, macs2_org, macs2_org, macs2_org, macs2_org, macs2_org)

	scripts2allfrag = """
>&2 echo "Converting bedgraph to bigwig... ""$base".bam
>&2 date
cd $outdir
LC_ALL=C sort -k1,1 -k2,2n $outdir/"$base_file"_treat_pileup.bdg > $outdir/"$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$base_file".sort.bdg $chromsizedir/%s.chrom.sizes $outdir/"$base_file".sorted.bw
rm -rf "$base_file".sort.bdg

cd $outdir2
LC_ALL=C sort -k1,1 -k2,2n $outdir2/"$base_file"_treat_pileup.bdg > $outdir2/"$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig $outdir2/"$base_file".sort.bdg $chromsizedir/%s.chrom.sizes $outdir2/"$base_file".sorted.bw
rm -rf "$base_file".sort.bdg
""" % (config["input/output"]["organism_build"], config["input/output"]["organism_build"])


	scriptfinal = """
>&2 echo "Finished"
>&2 date
"""

	outp.write(header + "\n")
	outp.write(path + "\n")
	outp.write(scripts + "\n")
	outp.write(macs_script + "\n")
	outp.write(scripts2 + "\n")
	outp.write(scriptallfrag + "\n")
	outp.write(scripts2allfrag + "\n")
	outp.write(scriptfinal + "\n")


	if output is not None:
		outp.close()

def generate_integrated_motif_find_sh(config, output=None, peak="narrowPeak", is_seacr=False):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
	
	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%s                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_motif_find"]["time_limit"], config["cluster"]["step_motif_find"]["queue"], 
	config["cluster"]["step_motif_find"]["memory"], config["cluster"]["email"])
	
	scripts = """memebin=%s
bedopsbin=%s
bedtoolsbin=%s
pythonbin=%s
perlbin=%s
genome_sequence=%s
extrasettings=%s
blacklist=$extrasettings/%s.blacklist.bed

i=$1 #filename must end with .narrowPeak or .broadPeak or .bed (if SEACR)
>&2 echo "Input file is $i"

#expand the path for $1
relinfile=`realpath -s $i`
dirname=`dirname $relinfile`

#cd to current directory
cd $dirname

for d in blk_filtered; do
if [ ! -d $d ]; then
mkdir $d
fi
done
""" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], config["pythonbin"], \
	config["perlbin"], config["genome_sequence"], config["extrasettings"], \
	config["input/output"]["organism_build"])


	suffix = "_peaks.narrowPeak"
	summit_suffix = "_summits.bed"
	summit_padded_suffix = "_summits_padded.fa"
	if is_seacr==False and peak=="narrowPeak":
		suffix="_peaks.narrowPeak"
	elif is_seacr==False and peak=="broadPeak":
		suffix="_peaks.broadPeak"
	elif is_seacr==True:
		suffix="_treat.stringent.sort.bed"
		summit_suffix="_treat.stringent.sort.summits.bed"
		summit_padded_suffix="_treat.stringent.sort.summits_padded.fa"

	sc2 = """workdir=`pwd`
fname=`basename $i %s`
peak=$fname"%s"
summit=$fname"%s"
summitfa=$fname"%s"
""" % (suffix, suffix, summit_suffix, summit_padded_suffix)

	scripts2 = """>&2 echo "Get filtered peaks..."
cat $peak | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blk_filtered/$peak
cat $summit | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blk_filtered/$summit
"""
	scripts3 = """#motif discovery starts here
motif_dir=random.%d
msummit=$motif_dir/summits
mpadded=$motif_dir/padded
mpaddedfa=$motif_dir/padded.fa

for d in $motif_dir $msummit $mpadded $mpaddedfa; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "Get randomized %d peaks..."
cat blk_filtered/$peak | sort -t"	" -g -k8 -r | head -n %d | shuf | head -n %d | $bedopsbin/sort-bed - > $motif_dir/$peak
$bedopsbin/bedops -e 1 blk_filtered/$summit $motif_dir/$peak > $msummit/$summit
$bedopsbin/bedops --range %d -u $msummit/$summit > $mpadded/$summit

$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $mpadded/$summit -fo $mpaddedfa/$summitfa

>&2 echo "Start MEME analysis for de novo motif finding..."
meme_outdir=$motif_dir/MEME_"$fname"_shuf
$memebin/meme-chip -oc $meme_outdir -dreme-m %d -meme-nmotifs %d $mpaddedfa/$summitfa

>&2 echo "Finished"
""" % (config["motif_finding"]["num_peaks"], config["motif_finding"]["num_peaks"], 
	config["motif_finding"]["total_peaks"], 
	config["motif_finding"]["num_peaks"], config["motif_finding"]["num_bp_from_summit"], 
	config["motif_finding"]["num_motifs"], config["motif_finding"]["num_motifs"])

	outp.write(header + "\n")
	outp.write(scripts + "\n")
	outp.write(sc2 + "\n")
	outp.write(scripts2 + "\n")
	outp.write(scripts3 + "\n")

	if output is not None:
		outp.close()

def generate_integrated_footprinting_sh(config, dedup=False, output=None, is_all_frag=False, peak="narrowPeak", is_seacr=False):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	suffix = "_peaks.narrowPeak"
	if is_seacr==False and peak=="narrowPeak":
		suffix="_peaks.narrowPeak"
	elif is_seacr==False and peak=="broadPeak":
		suffix="_peaks.broadPeak"
	elif is_seacr==True:
		suffix="_treat.stringent.sort.bed"

	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%s                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_footprinting"]["time_limit"], config["cluster"]["step_footprinting"]["queue"], 
	config["cluster"]["step_footprinting"]["memory"], config["cluster"]["email"])

	scripts = """
pythonbin=%s
peak_file=$1 #a narrowPeak/broadPeak/SEACR bed file
mbase=`basename $peak_file %s`
peak=$mbase"%s"
mdiscovery=random.%d/MEME_"$mbase"_shuf

#expand the path for $peak_file
relinfile=`realpath -s $peak_file`
dirname=`dirname $relinfile`

#cd to current directory (macs2.narrow.aug10)
cd $dirname

$pythonbin/python read.meme.py $mdiscovery
""" % (config["pythonbin"], suffix, suffix, config["motif_finding"]["num_peaks"])

	p_pythonbase = config["pythonbin"].rstrip("/").rstrip("/bin")

	scripts2 = """
memebin=%s
bedopsbin=%s
bedtoolsbin=%s
genome_sequence=%s
samtoolsbin=%s
makecutmatrixbin=%s
Rscriptbin=%s
extrasettings=%s

pythonldlibrary=%s
ldlibrary=`echo $LD_LIBRARY_PATH | tr : "\\n" | grep -v $pythonldlibrary | paste -s -d:`
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$pythonldlibrary:$ldlibrary

p=%.5f
motif_dir=$mdiscovery/motifs #a directory containing a list of *.meme files
peak_filename=`basename $peak_file`
workdir=`pwd`
dir=blk_filtered
fa_dir=blk_filtered.fa

if [ ! -d $fa_dir ]; then
mkdir $fa_dir
fi

if [ ! -d $dir ] || [ ! -f $dir/$peak ] ; then
blacklist=$extrasettings/%s.blacklist.bed
cat $workdir/$dir/"$peak_filename" | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > $workdir/$dir/$peak
fi

""" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], 
	config["genome_sequence"], config["samtoolsbin"], 
	config["makecutmatrixbin"], config["Rscriptbin"],
	config["extrasettings"], p_pythonbase + "/lib", 
	config["motif_finding"]["motif_scanning_pval"], 
	config["input/output"]["organism_build"])

	scripts3 = """
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/$peak -fo $fa_dir/"$mbase".fa
$pythonbin/python fix_sequence.py $fa_dir/"$mbase".fa

outdir=fimo.result
for d in $outdir $outdir/$mbase; do
if [ ! -d $d ]; then
mkdir $d
fi
done

for m in `ls -1 $motif_dir`; do
motif=`basename $m .meme`
fimo_d=$outdir/$mbase/fimo2.$motif
if [ ! -d $fimo_d ]; then
mkdir $fimo_d
fi
$memebin/fimo --thresh $p --parse-genomic-coord -oc $fimo_d $motif_dir/"$motif".meme $fa_dir/"$mbase".fa
cur_path=`echo $PATH | tr : "\\n" | grep -v $bedopsbin | paste -s -d:`
unset PATH
export PATH=$cur_path:$bedopsbin

$bedopsbin/gff2bed < $fimo_d/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $fimo_d/fimo.bed
done
""" 

	bamdir = "dup.marked.120bp"
	if is_all_frag==True and dedup==True:
		bamdir = "dedup"
	elif is_all_frag==True:
		bamdir = "dup.marked"
	elif is_all_frag==False and dedup==True:
		bamdir = "dedup.120bp"
	else:
		bamdir = "dup.marked.120bp"

	scripts4a = """
bamfile=../aligned.aug10/%s/"$mbase".bam
""" % bamdir

	scripts4 = """
workdir=`pwd`
dir=`dirname $bamfile`
bambase=`basename $bamfile .bam`

dest=centipede.bam
outbam=$dest/"$bambase".bam
if [ ! -d $dest ]; then
mkdir $dest
fi
"""

	scripts5 = """
cd $dest
ln -s ../../aligned.aug10/%s/"$mbase".bam .
ln -s ../../aligned.aug10/%s/"$mbase".bam.bai .
cd ..
""" % (bamdir, bamdir)

	scripts6 = """
fimo_dir=$outdir/$mbase

for i in `ls -1 $fimo_dir`; do #shows a list of motifs
echo "Doing $i..."
fimo_d=$fimo_dir/$i
tmp=`echo $i|cut -d "." -f3|wc -c`
mlen=$(( tmp - 1 ))
$makecutmatrixbin/make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p 1 -f 3 -F 4 -F 8 -q 0 $outbam $fimo_d/fimo.bed > $fimo_d/fimo.cuts.freq.txt
$Rscriptbin/Rscript run_centipede_parker.R $fimo_d/fimo.cuts.freq.txt $fimo_d/fimo.bed $fimo_d/fimo.png $mlen
done
"""
	outp.write(header + "\n")
	outp.write(scripts + "\n")
	outp.write(scripts2 + "\n")
	outp.write(scripts3 + "\n")
	outp.write(scripts4a + "\n")
	outp.write(scripts4 + "\n")
	outp.write(scripts5 + "\n")
	outp.write(scripts6 + "\n")

	if output is not None:
		outp.close()

def generate_integrated_all_steps_sh(output=None):
	script = """#!/bin/bash

sample=$1
#expand the path for $sample
relsample=`realpath -s $sample`
dirname=`dirname $relsample`
base=`basename $relsample _R1_001.fastq.gz`

cd $dirname

echo "Submitting job 1..."
jid1=$(sbatch ./integrated.sh $relsample)
jid1=${jid1##* }
echo "Submitting job 2..."
jid2=$(sbatch -d afterany:$jid1 aligned.aug10/integrated.step2.sh aligned.aug10/"$base"_aligned_reads.bam)
jid2=${jid2##* }
echo "Submitting job 3..."
jid3=$(sbatch -d afterany:$jid2 macs2.narrow.aug18/integrate.motif.find.sh macs2.narrow.aug18/"$base"_aligned_reads_peaks.narrowPeak)
jid3=${jid3##* }
echo "Submitting job 4..."
jid4=$(sbatch -d afterany:$jid3 macs2.narrow.aug18/integrate.footprinting.sh macs2.narrow.aug18/"$base"_aligned_reads_peaks.narrowPeak)
jid4=${jid4##* }

echo "Submitting job 5..."
jid5=$(sbatch -d afterany:$jid4 macs2.narrow.aug18.dedup/integrate.motif.find.sh macs2.narrow.aug18.dedup/"$base"_aligned_reads_peaks.narrowPeak)
jid5=${jid5##* }
echo "Submitting job 6..."
jid6=$(sbatch -d afterany:$jid5 macs2.narrow.aug18.dedup/integrate.footprinting.sh macs2.narrow.aug18.dedup/"$base"_aligned_reads_peaks.narrowPeak)
jid6=${jid6##* }
"""
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
	outp.write(script + "\n")
	if output is not None:
		outp.close()


	


def generate_single_locus_script_sh(config, dedup=False, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	if dedup:
		samtools_flags = "-f 3 -F 4 -F 8 -F 1024"
	else:
		samtools_flags = "-f 3 -F 4 -F 8"

	script = """#!/bin/bash

region=$1
bamfile=$2
outdir=$3
chromsizedir=`dirname %s`
chromsizefile=$chromsizedir/%s.chrom.sizes
pythonbin=%s
samtoolsbin=%s
bedtoolsbin=%s
bedopsbin=%s
extratoolsbin=%s
samtoolsflags="%s"

regionname=`echo $region|sed "s/:/-/g"`
basename=`basename $bamfile .bam`
newbamfile="$basename"-"$regionname".bam
newbase=`basename $newbamfile .bam`

if [ ! -d $outdir ]; then
mkdir $outdir
fi
$samtoolsbin/samtools view -bh $samtoolsflags $bamfile "$region" > $outdir/$newbamfile
$samtoolsbin/samtools index $outdir/$newbamfile
$samtoolsbin/samtools view -b $outdir/$newbamfile|$samtoolsbin/samtools sort -O bam -n - -T tmp.test|$bedtoolsbin/bedtools bamtobed -i stdin -bedpe > $outdir/"$newbase".frag.ends.txt

$pythonbin/python check_coordinate.py $chromsizefile $outdir/"$newbase".frag.ends.txt > $outdir/"$newbase".frag.ends.checked.txt

$pythonbin/python quantify_separate.py $outdir/"$newbase".frag.ends.checked.txt $outdir/"$newbase".frag.ends.R1.bed $outdir/"$newbase".frag.ends.R2.bed
$bedopsbin/sort-bed $outdir/"$newbase".frag.ends.R1.bed > $outdir/"$newbase".frag.ends.R1.sorted.bed
$bedopsbin/sort-bed $outdir/"$newbase".frag.ends.R2.bed > $outdir/"$newbase".frag.ends.R2.sorted.bed
$bedtoolsbin/groupBy -i $outdir/"$newbase".frag.ends.R1.sorted.bed -g 1,2,3 -c 2 -o count > $outdir/"$newbase".frag.ends.R1.bdg
$bedtoolsbin/groupBy -i $outdir/"$newbase".frag.ends.R2.sorted.bed -g 1,2,3 -c 2 -o count > $outdir/"$newbase".frag.ends.R2.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$newbase".frag.ends.R1.bdg $chromsizefile $outdir/"$newbase".frag.ends.R1.bw
$extratoolsbin/bedGraphToBigWig $outdir/"$newbase".frag.ends.R2.bdg $chromsizefile $outdir/"$newbase".frag.ends.R2.bw

$pythonbin/python quantify.py $outdir/"$newbase".frag.ends.checked.txt $outdir/"$newbase".frag.ends.bed
$bedopsbin/sort-bed $outdir/"$newbase".frag.ends.bed > $outdir/"$newbase".frag.ends.sorted.bed
$bedtoolsbin/groupBy -i $outdir/"$newbase".frag.ends.sorted.bed -g 1,2,3 -c 2 -o count > $outdir/"$newbase".frag.ends.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$newbase".frag.ends.bdg $chromsizefile $outdir/"$newbase".frag.ends.bw

""" % (config["genome_sequence"], config["input/output"]["organism_build"], config["pythonbin"], config["samtoolsbin"], config["bedtoolsbin"], \
config["bedopsbin"], config["extratoolsbin"], samtools_flags)

	outp.write(script + "\n")
	if output is not None:
		outp.close()



def write_length_file(n, length):
	fw = open(n, "w")
	fw.write(str(length) + "\n")
	fw.close()


def generate_TF_footprinting_script_sh(config, dedup=False, output=None, is_all_frag=False, peak="narrowPeak", is_seacr=False):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	dedup_str = "False"
	if dedup:
		dedup_str = "True"
	else:
		dedup_str = "False"

	suffix = "_peaks.narrowPeak"
	if is_seacr==False and peak=="narrowPeak":
		suffix="_peaks.narrowPeak"
	elif is_seacr==False and peak=="broadPeak":
		suffix="_peaks.broadPeak"
	elif is_seacr==True:
		suffix="_treat.stringent.sort.bed"

	bamdir = "dup.marked.120bp"
	if is_all_frag==True and dedup==True:
		bamdir = "dedup"
	elif is_all_frag==True:
		bamdir = "dup.marked"
	elif is_all_frag==False and dedup==True:
		bamdir = "dedup.120bp"
	else:
		bamdir = "dup.marked.120bp"


	script="""#!/usr/bin/python
import shutil
import sys
import os
import re
import json
import argparse

def generate_TF_specific_footprint_script_sh(config, pval, motif_file, tf_interest, output=None, dedup=False):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
	header = \"\"\"#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%s                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
\"\"\" % (config["cluster"]["step_footprinting"]["time_limit"], config["cluster"]["step_footprinting"]["queue"], 
	config["cluster"]["step_footprinting"]["memory"], config["cluster"]["email"])

	script = \"\"\"
pythonbin=%s

peak_file=$1 #a narrowPeak file
peak_filename=`basename $peak_file`
mbase=""
summit=""
summitfa=""
if [[ "$peak_file" == *narrowPeak ]]
then
	mbase=`basename $peak_file _peaks.narrowPeak`
	summit=$mbase"_summits.bed"
	summitfa=$mbase"_summits_padded.fa"
elif [[ "$peak_file" == *broadPeak ]]
then
	mbase=`basename $peak_file _peaks.broadPeak`
	summit=$mbase"_summits.bed"
	summitfa=$mbase"_summits_padded.fa"
elif [[ "$peak_file" == *stringent.sort.bed ]]
then
	mbase=`basename $peak_file _treat.stringent.sort.bed`
	summit=$mbase"_treat.stringent.sort.summits.bed"
	summitfa=$mbase"_treat.stringent.sort.summits_padded.fa"
fi

#expand the path for $peak_file
relinfile=`realpath -s $peak_file`
dirname=`dirname $relinfile`

#cd to current directory (macs2.narrow.aug10)
cd $dirname
\"\"\" % config["pythonbin"]

	p_pythonbase = config["pythonbin"].rstrip("/").rstrip("/bin")

	script2 = \"\"\"
memebin=%s
bedopsbin=%s
bedtoolsbin=%s
genome_sequence=%s
samtoolsbin=%s
makecutmatrixbin=%s
Rscriptbin=%s
extrasettings=%s

pythonldlibrary=%s
ldlibrary=`echo $LD_LIBRARY_PATH | tr : "\\n" | grep -v $pythonldlibrary | paste -s -d:`
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$pythonldlibrary:$ldlibrary

p=%.5f
motif_file=%s
workdir=`pwd`

dir=blk_filtered
fa_dir=blk_filtered.fa
if [ ! -d $fa_dir ]; then
mkdir $fa_dir
fi
if [ ! -d $dir ]; then
mkdir $dir
fi

if [ ! -f $dir/"$peak_filename" ]; then
blacklist=$extrasettings/%s.blacklist.bed
cat $peak_file | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > $dir/$peak_filename
cat $summit | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > $dir/$summit
fi

\"\"\" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], 
	config["genome_sequence"], config["samtoolsbin"], 
	config["makecutmatrixbin"], config["Rscriptbin"],
	config["extrasettings"], p_pythonbase + "/lib", 
	pval, motif_file, 
	config["input/output"]["organism_build"])

	script3 = \"\"\"
echo "get fasta"
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/"$peak_filename" -fo $fa_dir/"$mbase".fa
echo "fix sequence"
$pythonbin/python fix_sequence.py $fa_dir/"$mbase".fa

outdir=fimo.%s.result
for d in $outdir $outdir/$mbase; do
if [ ! -d $d ]; then
mkdir $d
fi
done

motif=`basename $motif_file .meme`
fimo_d=$outdir/$mbase/fimo2.$motif
if [ ! -d $fimo_d ]; then
mkdir $fimo_d
fi
echo "get fimo"
$memebin/fimo --thresh $p --parse-genomic-coord -oc $fimo_d "$motif".meme $fa_dir/"$mbase".fa

cur_path=`echo $PATH | tr : "\n" | grep -v $bedopsbin | paste -s -d:`
unset PATH
export PATH=$cur_path:$bedopsbin

$bedopsbin/gff2bed < $fimo_d/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $fimo_d/fimo.bed
\"\"\" % (tf_interest)
	script4 = \"\"\"
"""


	scripta="""
bamfile=../aligned.aug10/%s/"$mbase".bam

workdir=`pwd`
dir=`dirname $bamfile`
bambase=`basename $bamfile .bam`

dest=centipede.bam
outbam=$dest/"$bambase".bam
if [ ! -d $dest ]; then
mkdir $dest
fi
\"\"\"
	script5 = \"\"\"
cd $dest
ln -s ../../aligned.aug10/%s/"$mbase".bam .
ln -s ../../aligned.aug10/%s/"$mbase".bam.bai .
cd ..
\"\"\"
""" % (bamdir, bamdir, bamdir)

	scriptb="""
	script7 = \"\"\"
#peakfile=blk_filtered/"$base"_peaks.narrowPeak
fimo_dir=$outdir/"$mbase"

for i in `ls -1 $fimo_dir`; do #shows a list of motifs
echo "Doing $i..."
fimo_d=$fimo_dir/$i
tmp=`cat $motif_file |grep "MOTIF"|cut -d" " -f3|wc -c`
mlen=$(( tmp - 1 ))
$makecutmatrixbin/make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p 1 -f 3 -F 4 -F 8 -q 0 $outbam $fimo_d/fimo.bed > $fimo_d/fimo.cuts.freq.txt
$Rscriptbin/Rscript run_centipede_parker.R $fimo_d/fimo.cuts.freq.txt $fimo_d/fimo.bed $fimo_d/fimo.png $mlen
done
\"\"\"
	outp.write(header + "\\n")
	outp.write(script + "\\n")
	outp.write(script2 + "\\n")
	outp.write(script3 + "\\n")
	outp.write(script4 + "\\n")
	outp.write(script5 + "\\n")
	outp.write(script7 + "\\n")

	if output is not None:
		outp.close()

if __name__=="__main__":
	f = open("../current_config.json")
	config = json.load(f)
	f.close()

	parser = argparse.ArgumentParser(description="generates a footprint script for a given TF PWM matrix", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-b", "--motif", dest="motif", help="motif file in MEME format", type=str, required=True)
	parser.add_argument("-p", "--pvalue", dest="pvalue", help="pvalue cutoff for motif scanning in FIMO", type=float, required=True)
	parser.add_argument("-n", "--name", dest="name", help="name of factor", type=str, required=True)

	args = parser.parse_args()

	pval = args.pvalue
	motif_file = args.motif
	tf_interest = args.name
""" 
	script2 = """
	outfile = "integrate.footprinting.%%s.centipede.sh" %% tf_interest
	generate_TF_specific_footprint_script_sh(config, pval, motif_file, tf_interest, output=outfile, dedup=%s)
	st = os.stat(outfile)
	os.chmod(outfile, st.st_mode | 0o111)
""" % dedup_str

	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
	outp.write(script + "\n")
	outp.write(scripta + "\n")
	outp.write(scriptb + "\n")
	outp.write(script2 + "\n")
	if output is not None:
		outp.close()
	
if __name__=="__main__":
	curpath = os.path.dirname(os.path.abspath(__file__))
	#print curpath
	f = open(sys.argv[1]) #config.json
	config = json.load(f)
	f.close()

	outdir = config["input/output"]["workdir"]
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	write_length_file(config["input/output"]["workdir"] + "/length", config["input/output"]["fastq_sequence_length"])

	if not os.path.isdir(outdir + "/aligned.aug10"):
		os.mkdir(outdir + "/aligned.aug10")

	dir_create = ["macs2.narrow.aug18", "macs2.broad.aug18", \
	"macs2.narrow.all.frag.aug18", "macs2.broad.all.frag.aug18", \
	"macs2.narrow.aug18.dedup", "macs2.broad.aug18.dedup", \
	"macs2.narrow.all.frag.aug18.dedup", "macs2.broad.all.frag.aug18.dedup", \
	"seacr.aug12", "seacr.aug12.all.frag", "seacr.aug12.dedup", "seacr.aug12.all.frag.dedup"]

	for dc in dir_create:
		if not os.path.isdir(outdir + "/" + dc):
			os.mkdir(outdir + "/" + dc)

	generate_integrated_sh(config, output=outdir+"/integrated.sh")
	generate_integrated_step2_sh(config, output=outdir+"/aligned.aug10/integrated.step2.sh")

	settings = {\
	"macs2.narrow.aug18":          {"dedup":False, "is_all_frag":False, "peak":"narrowPeak", "is_seacr":False}, \
	"macs2.broad.aug18":           {"dedup":False, "is_all_frag":False, "peak":"broadPeak", "is_seacr":False}, \
	"macs2.narrow.all.frag.aug18": {"dedup":False, "is_all_frag":True,  "peak":"narrowPeak", "is_seacr":False}, \
	"macs2.broad.all.frag.aug18":  {"dedup":False, "is_all_frag":True,  "peak":"broadPeak", "is_seacr":False}, \
	"macs2.narrow.all.frag.aug18.dedup":{"dedup":True, "is_all_frag":True, "peak":"narrowPeak", "is_seacr":False}, \
	"macs2.broad.all.frag.aug18.dedup": {"dedup":True, "is_all_frag":True, "peak":"broadPeak", "is_seacr":False}, \
	"macs2.narrow.aug18.dedup":{"dedup":True, "is_all_frag":False, "peak":"narrowPeak", "is_seacr":False}, \
	"macs2.broad.aug18.dedup": {"dedup":True, "is_all_frag":False, "peak":"broadPeak", "is_seacr":False}, \
	"seacr.aug12":       {"dedup":False, "peak":None, "is_all_frag":False, "is_seacr":True}, \
	"seacr.aug12.dedup": {"dedup":True,  "peak":None, "is_all_frag":False, "is_seacr":True}, \
	"seacr.aug12.all.frag":       {"dedup":False, "peak":None, "is_all_frag":True, "is_seacr":True}, \
	"seacr.aug12.all.frag.dedup": {"dedup":True,  "peak":None, "is_all_frag":True, "is_seacr":True} }

	for dc in dir_create:
		par = [settings[dc]["dedup"], settings[dc]["is_all_frag"], settings[dc]["peak"], settings[dc]["is_seacr"]]
		generate_integrated_motif_find_sh(config, peak=par[2], is_seacr=par[3], output=outdir+"/%s/integrate.motif.find.sh" % dc)
		generate_integrated_footprinting_sh(config, dedup=par[0], is_all_frag=par[1], peak=par[2], is_seacr=par[3], output=outdir+"/%s/integrate.footprinting.sh" % dc)
		generate_single_locus_script_sh(config, dedup=par[0], output=outdir+"/%s/get_cuts_single_locus.sh" % dc)
		generate_TF_footprinting_script_sh(config, dedup=par[0], is_all_frag=par[1], peak=par[2], is_seacr=par[3], output=outdir+"/%s/generate.footprinting.factor.specific.centipede.py" % dc)

	generate_integrated_all_steps_sh(output=outdir+"/integrated.all.steps.sh")
	make_executable(outdir+"/integrated.sh")
	make_executable(outdir+"/integrated.all.steps.sh")
	make_executable(outdir+"/aligned.aug10/integrated.step2.sh")

	for dc in dir_create:
		make_executable(outdir+"/%s/integrate.motif.find.sh" % dc)
		make_executable(outdir+"/%s/integrate.footprinting.sh" % dc)
		make_executable(outdir+"/%s/get_cuts_single_locus.sh" % dc)
		make_executable(outdir+"/%s/generate.footprinting.factor.specific.centipede.py" % dc)
		shutil.copyfile("%s/macs2.narrow.aug18/filter.py" % (curpath), outdir+"/%s/filter.py" % dc)
		shutil.copyfile("%s/macs2.narrow.aug18/fix_sequence.py" % (curpath), outdir+"/%s/fix_sequence.py" % dc)
		shutil.copyfile("%s/macs2.narrow.aug18/read.meme.py" % (curpath), outdir+"/%s/read.meme.py" % dc)
		shutil.copyfile("%s/macs2.narrow.aug18/run_centipede_parker.R" % (curpath), outdir+"/%s/run_centipede_parker.R" % dc)
		shutil.copyfile("%s/quantify.py" % curpath, outdir+"/%s/quantify.py" % dc)
		shutil.copyfile("%s/quantify_separate.py" % curpath, outdir+"/%s/quantify_separate.py" % dc)
		shutil.copyfile("%s/check_coordinate.py" % curpath, outdir+"/%s/check_coordinate.py" % dc)

	shutil.copyfile(sys.argv[1], outdir+"/current_config.json")

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

	if config["input/output"]["fastq_directory"]!=config["input/output"]["workdir"]:
		for fx in flist:
			if os.path.islink(config["input/output"]["workdir"] + "/" + fx):
				os.remove(config["input/output"]["workdir"] + "/" + fx)
			os.symlink(config["input/output"]["fastq_directory"] + "/" + fx, \
			config["input/output"]["workdir"] + "/" + fx)
