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


#=============================================
#remember in macs2 peak calling to change hs to mm for mouse
#remember in the line "extratoolsbin/bedGraphToBigWig" to change the hg19.chrom.sizes to mm10.chrom.sizes
#=============================================

#Usage: ./integrated.step2.sh CR_NFYA_RNP115_Day4_r1_S12_aligned_reads.bam
#input is a bam file generated from step 1

module load picard/2.8.0
module load samtools
module load gcc/6.2.0 python/2.7.12
module load macs2/2.1.1.20160309

picardbin=/n/app/picard/2.8.0/bin
samtoolsbin=/n/app/samtools/1.3.1/bin
macs2bin=/n/app/macs2/2.1.1.20160309/bin
extratoolsbin=/home/qz64
extrasettings=/home/qz64

workdir=`pwd`
logdir=$workdir/logs
outdir=../macs2.narrow.aug18 #for macs2

mkdir $logdir
mkdir $outdir

>&2 echo "Input parameters are: $1"
>&2 date
base=`basename $1 .bam`
#requires the filter_below.awk file for filtering

mkdir sorted dup.marked dedup
>&2 echo "Sorting bam... ""$base".bam
>&2 date
java -jar $picardbin/picard-2.8.0.jar SortSam \
INPUT="$base".bam OUTPUT=sorted/"$base".bam SORT_ORDER=coordinate

>&2 echo "Marking duplicates... ""$base".bam
>&2 date
java -jar $picardbin/picard-2.8.0.jar MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dup.marked/"$base".bam \
METRICS_FILE=metrics."$base".txt

>&2 echo "Removing duplicates... ""$base".bam
>&2 date
java -jar $picardbin/picard-2.8.0.jar MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dedup/"$base".bam \
METRICS_FILE=metrics."$base".txt \
REMOVE_DUPLICATES=true

mkdir sorted.120bp dup.marked.120bp dedup.120bp
>&2 echo "Filtering to <120bp... ""$base".bam
>&2 date
$samtoolsbin/samtools view -h sorted/"$base".bam |awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > sorted.120bp/"$base".bam
$samtoolsbin/samtools view -h dup.marked/"$base".bam |awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dup.marked.120bp/"$base".bam
$samtoolsbin/samtools view -h dedup/"$base".bam |awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dedup.120bp/"$base".bam

>&2 echo "Creating bam index files... ""$base".bam
>&2 date
$samtoolsbin/samtools index sorted/"$base".bam
$samtoolsbin/samtools index dup.marked/"$base".bam
$samtoolsbin/samtools index dedup/"$base".bam
$samtoolsbin/samtools index sorted.120bp/"$base".bam
$samtoolsbin/samtools index dup.marked.120bp/"$base".bam
$samtoolsbin/samtools index dedup.120bp/"$base".bam

>&2 echo "Peak calling using MACS2... ""$base".bam
>&2 echo "Logs are stored in $logdir"
>&2 date
bam_file=dup.marked.120bp/"$base".bam
dir=`dirname $bam_file`
base_file=`basename $bam_file .bam`

cd $outdir
#$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g hs -f BAMPE -n $base_file --outdir $outdir -q 0.01 -B --SPMR 2> $logdir/"$base_file".macs2
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g hs -f BAMPE -n $base_file --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2

>&2 echo "Converting bedgraph to bigwig... ""$base".bam
>&2 date
sort -k1,1 -k2,2n "$base_file"_treat_pileup.bdg > "$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig "$base_file".sort.bdg $extrasettings/hg19.chrom.sizes "$base_file".sorted.bw
rm -rf "$base_file".sort.bdg

>&2 echo "Finished"
>&2 date
