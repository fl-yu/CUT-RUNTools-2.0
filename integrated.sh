#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=32000                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bernardzhu@gmail.com   # Email to which notifications will be sent

#=================================================
#remember to change hg19 to mm10 if needed
#remember to create/modify file "length" with the length of the reads (currently 42)
#=================================================

#Usage: ./integrated.sh CR_NFYA_RNP115_r3_S5_R1_001.fastq.gz
#Note: an experiment will have both *R1_001.fastq.gz and *R2_001.fastq.gz files. Specify the *R1_001.fastq.gz file.


module load trimmomatic/0.36
module load samtools/1.3.1
module load gcc/6.2.0
module load bowtie2/2.2.9

trimmomaticbin=/n/app/trimmomatic/0.36/bin
adapterpath=/home/qz64
bowtie2bin=/n/app/bowtie2/2.2.9/bin
samtoolsbin=/n/app/samtools/1.3.1/bin
len=`cat length`

workdir=`pwd`
trimdir=$workdir/trimmed
trimdir2=$workdir/trimmed3
logdir=$workdir/logs
bt2idx=/n/groups/shared_databases/bowtie2_indexes
aligndir=$workdir/aligned.aug10

mkdir $trimdir
mkdir $trimdir2
mkdir $logdir
mkdir $aligndir

infile=$1
base=`basename $infile _R1_001.fastq.gz`

>&2 echo "Input file is $infile"
>&2 date

#trimming paired-end
#good version
>&2 echo "Trimming file $base ..."
>&2 date
java -jar $trimmomaticbin/trimmomatic-0.36.jar PE -threads 1 -phred33 $workdir/"$base"_R1_001.fastq.gz $workdir/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

>&2 echo "Second stage trimming $base ..."
>&2 date
/home/qz64/kseq_test $trimdir/"$base"_1.paired.fastq.gz $len $trimdir2/"$base"_1.paired.fastq.gz
/home/qz64/kseq_test $trimdir/"$base"_2.paired.fastq.gz $len $trimdir2/"$base"_2.paired.fastq.gz

>&2 echo "Aligning file $base ..."
>&2 date
($bowtie2bin/bowtie2 -p 2 --dovetail --phred33 -x $bt2idx/hg19 -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/"$base"_aligned_reads.bam

>&2 echo "Finished"
>&2 date
