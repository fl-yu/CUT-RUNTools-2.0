# ==================================================================
# step1 fastq->peak (raw data processing and quality assessment)
# ==================================================================

#len=`cat length` 
trimdir=$workdir/sc_trimmed
trimdir2=$workdir/sc_trimmed2
logdir=$workdir/sc_logs
aligndir=$workdir/sc_aligned.aug10
qcdir=$workdir/sc_qc
len=25

echo "# "
echo "#         step1 fastq->peak (raw data processing and quality assessment)"
echo "# "

mkdir -p $workdir

for d in $trimdir $trimdir2 $logdir $aligndir $qcdir; do
if [ ! -d $d ]; then
mkdir $d
fi
done

cd $fastq_directory
infiles=`echo *_R1_001.fastq.gz | sed s/_R1_001.fastq.gz//g`
file_num=`ls *_R1_001.fastq.gz | wc -l`
>&2 echo "[info] $file_num PE fastq files were detected ..."

# begin to work
cd $workdir
# ----------------------------------------------------------------------------------------------------
# choose if trim the FASTQ files, default: true
two_step_trim=true
if $two_step_trim
then
    mkdir -p $logdir/trim
    # trimming paired-end
    >&2 echo "[info] First stage trimming ..."
    >&2 date
    $path_parallel/parallel -k -j $cores "$javabin/java -jar $trimmomaticbin/$trimmomaticjarfile PE -phred33 $fastq_directory/{}_R1_001.fastq.gz $fastq_directory/{}_R2_001.fastq.gz $trimdir/{}_R1_001.fastq.gz $trimdir/{}_R1_001.unpaired.fastq.gz $trimdir/{}_R2_001.fastq.gz $trimdir/{}_R2_001.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true > $logdir/trim/{}_trim1.log 2>&1" ::: $infiles >/dev/null 2>&1

    >&2 echo "[info] Second stage trimming ..."
    >&2 date
    $path_parallel/parallel -k -j $cores "$kseqbin/kseq_test $trimdir/{}_R1_001.fastq.gz $len $trimdir2/{}_R1_001.fastq.gz" ::: $infiles >/dev/null 2>&1
    $path_parallel/parallel -k -j $cores "$kseqbin/kseq_test $trimdir/{}_R2_001.fastq.gz $len $trimdir2/{}_R2_001.fastq.gz" ::: $infiles >/dev/null 2>&1

else
    trimdir2=$fastq_directory
fi

# ----------------------------------------------------------------------------------------------------
# confirm the reference genome
if [ $genome = "hg38" ]
then
    genome_prefix=GRCh38
elif [ $genome = "hg19" ]
then
    genome_prefix=hg19
elif [ $genome = "mm10" ]
then
    genome_prefix=mm10
elif [ $genome = "mm9" ]
then
    genome_prefix=mm9
else
    echo "[info] Parameter genome should be one of hg38, hg19, mm10, mm9"
fi

>&2 echo "[info] Aligning FASTQ files ..."
>&2 date
mkdir -p $logdir/bowtie2
$path_parallel/parallel -k -j $cores "($bowtie2bin/bowtie2 --phred33 -x $bt2idx/$genome_prefix -1 $trimdir2/{}_R1_001.fastq.gz -2 $trimdir2/{}_R2_001.fastq.gz) 2> $logdir/bowtie2/{}.bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/{}_aligned_reads.bam" ::: $infiles >/dev/null 2>&1
>&2 echo "[info] Alignment finished ..."
>&2 date

cd $aligndir
for d in sorted dup.marked dup.marked.clean; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "[info] Sorting BAMs... "
>&2 date
mkdir -p $logdir/bamSort
$path_parallel/parallel -j $cores "$javabin/java -jar $picardbin/$picardjarfile SortSam INPUT=$aligndir/{}_aligned_reads.bam OUTPUT=sorted/{}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT > $logdir/bamSort/{}.sort 2>&1" ::: $infiles >/dev/null 2>&1

>&2 echo "[info] Marking duplicates... "
>&2 date
mkdir -p $logdir/bamMarkdup
$path_parallel/parallel -j $cores "$javabin/java -jar $picardbin/$picardjarfile MarkDuplicates INPUT=sorted/{}.bam OUTPUT=dup.marked/{}.bam VALIDATION_STRINGENCY=SILENT METRICS_FILE=$logdir/bamMarkdup/{}.metrics.txt > $logdir/bamMarkdup/{}.markdup 2>&1" ::: $infiles >/dev/null 2>&1
$path_parallel/parallel -j $cores "$samtoolsbin/samtools index dup.marked/{}.bam" ::: $infiles >/dev/null 2>&1

>&2 echo "[info] Filtering unmapped, low quality, unproper paired fragments... "
>&2 date
$path_parallel/parallel -j $cores "$samtoolsbin/samtools view -bh -q 30 -f 0x2 -F 1024 dup.marked/{}.bam > dup.marked.clean/{}.bam" ::: $infiles >/dev/null 2>&1
$path_parallel/parallel -j $cores "$samtoolsbin/samtools index dup.marked.clean/{}.bam" ::: $infiles >/dev/null 2>&1


>&2 echo "[info] Generating coverage .bed files used for single-cell genome track visualization ... "
>&2 date
$path_parallel/parallel -j $cores "$bedtoolsbin/bedtools genomecov -ibam dup.marked.clean/{}.bam -bg | cut -f 1-3 > dup.marked.clean/{}.bed" ::: $infiles >/dev/null 2>&1

>&2 echo "[info] Aggregation analysis of individual cells... "
>&2 date
# define the parameters
r_function_dir=$scriptdir/Rscript
bash_function_dir=$scriptdir/BASHscript
bam_dir=$workdir/sc_aligned.aug10/dup.marked.clean
sc_pseudoBulk_dir=$workdir/sc_pseudoBulk
# SEACR script should be in $extratoolsbin
mkdir -p $sc_pseudoBulk_dir 
$Rscriptbin/Rscript $r_function_dir/pseudobulk_analysis_aggregation.r \
                    $bam_dir \
                    $sc_pseudoBulk_dir \
                    $samtoolsbin \
                    $macs2bin \
                    $path_deeptools \
                    $genome \
                    $path_tabix \
                    $bash_function_dir \
                    $cores \
                    $extratoolsbin \
                    $Rscriptbin \
                    $peak_type \
                    $peak_caller

sc_pseudoBulk_dir=$workdir/sc_pseudoBulk
r_function_dir=$scriptdir/Rscript
>&2 echo "[info] Generating QC reports and figures... "
>&2 date
bamdir=$aligndir/dup.marked

if [[ $peak_caller = "macs2" ]]
then
    # ${peak_type}=narrow or broad
    peakfile=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/macs2/groups_aggregation_peaks.${peak_type}Peak
elif [[ $peak_caller = "SEACR" ]]
then
    peakfile=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/SEACR/groups_aggregation.stringent.bed
else
    echo "The \$peak_caller parameter is not set as \"macs2\" or \"SEACR\", please check it! Use \"macs2\" for the following processing! "
    peakfile=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/macs2/groups_aggregation_peaks.narrowPeak
fi
echo "[info] $peakfile will used for peak analysis .."

pseudo_bam=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/groups_aggregation_pseudo_sort.bam
# >&2 echo "[info] Removing the temporary files... "
rm -rf $qcdir/temp
mkdir -p $qcdir/temp
cd $logdir/bowtie2
####1 Barcode_name 
ls *.bowtie2 |sed 's/\.bowtie2//' | tr " " "\n" > $qcdir/temp/sc_barcode_name.temp
>&2 echo -e "[ok]Barcode_name"
####3 Overall_alignment_ratio
ls *.bowtie2 | xargs grep "overall alignment rate" | cut -d":" -f 2 | sed 's/% overall alignment rate//' > $qcdir/temp/sc_alignment_ratio.temp
>&2 echo -e "[ok]Overall_alignment_ratio"

mkdir -p $qcdir/temp/inset.temp
cd $bamdir
for i in `ls *.bam | sed 's/\.bam//'`; do
$samtoolsbin/samtools stats ${i}.bam > $qcdir/temp/inset.temp/${i}.stats
done
cd $qcdir/temp/inset.temp

####2 Total_reads (not fragment number) SN      raw total sequences:    100182
ls | xargs grep -e "raw total sequences" | cut -f 3  > $qcdir/temp/sc_total_pair.temp
>&2 echo -e "[ok]Total_reads_pair"

####4 Reads_properly_paired
ls | xargs grep -e "reads properly paired" | cut -f 3  > $qcdir/temp/sc_proper-pair.temp
>&2 echo -e "[ok]Reads_properly_paired" #SN      reads properly paired:  99948   # proper-pair bit set
####5 Reads_properly_paired_percentage 
ls | xargs grep -e "percentage of properly paired reads" | cut -f 3  > $qcdir/temp/sc_proper-pair_percentage.temp
>&2 echo -e "[ok]Reads_properly_paired_percentage" #SN      percentage of properly paired reads (%):	84.3
####6 Reads_duplicated
ls | xargs grep -e "reads duplicated" | cut -f 3  > $qcdir/temp/sc_reads_duplicated.temp
>&2 echo -e "[ok]Reads_duplicated" #reads duplicated:	8979	# PCR or optical duplicate bit set
####7 Reads_duplicated_percentage
paste $qcdir/temp/sc_reads_duplicated.temp $qcdir/temp/sc_total_pair.temp | awk '{printf "%.2f\n", $1*100/$2}' > $qcdir/temp/sc_reads_duplicated_percentage.temp
>&2 echo -e "[ok]Reads_duplicated_percentage"
####8 insert size average
ls | xargs grep -e "insert size average" | cut -f 3  > $qcdir/temp/sc_is_average.temp
>&2 echo -e "[ok]Insert_size_average"
####9 insert size standard deviation
ls | xargs grep -e "insert size standard deviation" | cut -f 3  > $qcdir/temp/sc_is_sd.temp
>&2 echo -e "[ok]Insert_size_standard_deviation"
#### 10 Reads_nuclear
cd $bamdir
for i in `ls *.bam`; do
$samtoolsbin/samtools idxstats ${i} | grep chrM | cut -f 3 >> $qcdir/temp/sc_chrM.temp
done
paste $qcdir/temp/sc_chrM.temp $qcdir/temp/sc_total_pair.temp | awk '{print $2-$1}' > $qcdir/temp/sc_nuclear.temp
>&2 echo -e "[ok]Reads_nuclear"
####11 Reads_nuclear_percentage
paste $qcdir/temp/sc_chrM.temp $qcdir/temp/sc_total_pair.temp | awk '{printf "%.2f\n", ($2-$1)*100/$2}' > $qcdir/temp/sc_nuclear_percentage.temp
>&2 echo -e "[ok]Reads_nuclear_percentage"
####12 Reads_in_peak
cd $bamdir
# for i in `ls *.bam | sed 's/\.bam//'`; do
# $bedtoolsbin/bedtools sort -i $peakfile | $bedtoolsbin/bedtools merge -i stdin | $bedtoolsbin/bedtools intersect -u -a ${i}.bam -b stdin -ubam | $samtoolsbin/samtools view -c >> $qcdir/temp/sc_rip.temp
# done
$path_parallel/parallel -k -j $cores "$bedtoolsbin/bedtools sort -i $peakfile | $bedtoolsbin/bedtools merge -i stdin | $bedtoolsbin/bedtools intersect -u -a {}.bam -b stdin -ubam | $samtoolsbin/samtools view -c " ::: `ls *.bam | sed 's/\.bam//'` > $qcdir/temp/sc_rip.temp
>&2 echo -e "[ok]Reads_in_peak"
####13 Reads_in_peak_percentage
paste $qcdir/temp/sc_rip.temp $qcdir/temp/sc_total_pair.temp | awk '{printf "%.2f\n", ($1/$2)*100}' > $qcdir/temp/sc_rip_percentage.temp
>&2 echo -e "[ok]Reads_in_peak_percentage"

####14 Reads_MAPQ30
cd $bamdir
for i in `ls *.bam`; do
samtools view -c -q 30 $i >> $qcdir/temp/sc_mapq30.temp
done
>&2 echo -e "[ok]Reads_MAPQ30"
####15 Reads_MAPQ30_percentage
paste $qcdir/temp/sc_mapq30.temp $qcdir/temp/sc_total_pair.temp | awk '{printf "%.2f\n", ($1/$2)*100}' > $qcdir/temp/sc_mapq30_percentage.temp
>&2 echo -e "[ok]Reads_MAPQ30_percentage"
####scQC.report
echo -e "Barcode_name\tTotal_reads\tOverall_alignment_ratio\tReads_properly_paired\tReads_properly_paired_percentage\tReads_duplicated\tReads_duplicated_percentage\tInsert_size_average\tInsert_size_sd\tReads_nuclear\tReads_nuclear_percentage\tReads_in_peak\tReads_in_peak_percentage\tReads_MAPQ30\tReads_MAPQ30_percentage" > $qcdir/temp/header.temp
paste $qcdir/temp/sc_barcode_name.temp $qcdir/temp/sc_total_pair.temp $qcdir/temp/sc_alignment_ratio.temp $qcdir/temp/sc_proper-pair.temp $qcdir/temp/sc_proper-pair_percentage.temp $qcdir/temp/sc_reads_duplicated.temp $qcdir/temp/sc_reads_duplicated_percentage.temp $qcdir/temp/sc_is_average.temp $qcdir/temp/sc_is_sd.temp $qcdir/temp/sc_nuclear.temp $qcdir/temp/sc_nuclear_percentage.temp $qcdir/temp/sc_rip.temp $qcdir/temp/sc_rip_percentage.temp $qcdir/temp/sc_mapq30.temp $qcdir/temp/sc_mapq30_percentage.temp| cat $qcdir/temp/header.temp - > $qcdir/temp/scQC.report

$samtoolsbin/samtools stats -@ $cores $pseudo_bam > $qcdir/temp/groups_aggregation_pseudo_bam.stats
cat $qcdir/temp/groups_aggregation_pseudo_bam.stats | grep ^IS | cut -f 2-3 > $qcdir/temp/insert_size_distribution.stat

#### QC filter for barcode after step1
num_reads_threshold=$num_reads_threshold
percentage_rip=$percentage_rip
mkdir -p $qcdir/figure
$Rscriptbin/Rscript $r_function_dir/qc_figure.r \
                    $qcdir/temp/scQC.report \
                    $qcdir/figure \
                    $qcdir/temp/insert_size_distribution.stat \
                    $num_reads_threshold \
                    $percentage_rip

>&2 echo "[info] QC figures generated ... "
>&2 echo "[info] QC reports generated... "
mkdir -p $qcdir/report
$Rscriptbin/Rscript $r_function_dir/qc_report.r \
                    $qcdir/temp/scQC.report \
                    $qcdir/report \
                    $num_reads_threshold \
                    $percentage_rip

# >&2 echo "[info] Removing the temporary files... "
rm -rf $qcdir/temp


>&2 echo "[info] Quality control criteria used to filter barcode cells"
>&2 echo "[info] Reads properly paired > $num_reads_threshold"
>&2 echo "[info] Percentage of Reads Enrichment in Peaks > $percentage_rip%"

if test -f "$qcdir/report/statistics_QCfailed.txt"; then
    echo "[info] `sed 1d $qcdir/report/statistics_QCfailed.txt | wc -l` cells fail to pass the QC thresholds, and `sed 1d $qcdir/report/statistics_QCpassed.txt | wc -l` cells retained for further analysis"
else 
    echo "[info] All cells passed QC criteria!"
fi

echo "# "
echo "#        step1 (raw data processing and quality assessment) is completed"
echo "# "
echo "-------------------------------------------------------------------------------------------------------------------------"
