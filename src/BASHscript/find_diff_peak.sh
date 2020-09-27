#!/bin/bash
bam1=$1
bam2=$2
peak1=$3
peak2=$4
summit1=$5
summit2=$6
experiment_name1=$7
experiment_name2=$8
outdir=$9
genome=${10}
extrasettings=${11}
bedopsbin=${12}
Homerbin=${13}
bedtoolsbin=${14}


genome=$genome # one of hg38 hg19 mm10 mm9
if [ $genome != "hg38" ] && [ $genome != "hg19" ] && [ $genome != "mm10" ] && [ $genome != "mm9" ]
then
    echo "Parameter genome should be one of hg38, hg19, mm10, mm9"
fi
blacklist=$extrasettings/$genome.blacklist.bed

date
echo "Filter and sort the peaks"

cat $peak1 |  egrep "^[chr0-9XY$]" | awk '$1 !~ "_"{ print $0}' |  grep -v chrM | $bedopsbin/sort-bed - | $bedtoolsbin/bedtools intersect -wa -v -a - -b $blacklist > $outdir/${experiment_name1}_peak.filter
cat $peak2 |  egrep "^[chr0-9XY$]" | awk '$1 !~ "_"{ print $0}' |  grep -v chrM | $bedopsbin/sort-bed - | $bedtoolsbin/bedtools intersect -wa -v -a - -b $blacklist > $outdir/${experiment_name2}_peak.filter

echo "Convert bam to tag"
$Homerbin/makeTagDirectory $outdir/${experiment_name1}_tag -sspe $bam1
$Homerbin/makeTagDirectory $outdir/${experiment_name2}_tag -sspe $bam2
# echo "Merge the peaks"
# cat $peak1 $peak2 | sort -k1,1 -k2,2n - | bedtools merge -d 100 -i - > $outdir/${experiment_name1}_${experiment_name2}.mergedPeak
echo "Find the differential peaks"
# find the d peaks
getDifferentialPeaks $outdir/${experiment_name1}_peak.filter $outdir/${experiment_name1}_tag $outdir/${experiment_name2}_tag > $outdir/${experiment_name1}_${experiment_name2}.diffPeak
getDifferentialPeaks $outdir/${experiment_name2}_peak.filter  $outdir/${experiment_name2}_tag $outdir/${experiment_name1}_tag > $outdir/${experiment_name2}_${experiment_name1}.diffPeak
# filter original peaks using the d peaks order
awk '$1 !~ "^#"{print $0}' $outdir/${experiment_name1}_${experiment_name2}.diffPeak | cut -f 1  | cut -d"_" -f 4 | awk 'NR==FNR{data[$1]; next}FNR in data' - $peak1 > $outdir/${experiment_name1}_${experiment_name2}.diffPeak.final
awk '$1 !~ "^#"{print $0}' $outdir/${experiment_name2}_${experiment_name1}.diffPeak | cut -f 1  | cut -d"_" -f 4 | awk 'NR==FNR{data[$1]; next}FNR in data' - $peak2 > $outdir/${experiment_name2}_${experiment_name1}.diffPeak.final
echo "Find the summits"
# filter original summits using the d peaks order
awk '$1 !~ "^#"{print $0}' $outdir/${experiment_name1}_${experiment_name2}.diffPeak | cut -f 1  | cut -d"_" -f 4 | awk 'NR==FNR{data[$1]; next}FNR in data' - $summit1 > $outdir/${experiment_name1}_${experiment_name2}.diffPeak.summit.final
awk '$1 !~ "^#"{print $0}' $outdir/${experiment_name2}_${experiment_name1}.diffPeak | cut -f 1  | cut -d"_" -f 4 | awk 'NR==FNR{data[$1]; next}FNR in data' - $summit2 > $outdir/${experiment_name2}_${experiment_name1}.diffPeak.summit.final
echo "Finish"
date

# ./find_diff_peak.sh \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/group_0/pseudo_bulk_bam/group_0_pseudo_sort.bam \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/group_1/pseudo_bulk_bam/group_1_pseudo_sort.bam \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/group_0/pseudo_bulk_bam/macs2/group_0_peaks.narrowPeak \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/group_1/pseudo_bulk_bam/macs2/group_1_peaks.narrowPeak \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/group_0/pseudo_bulk_bam/macs2/group_0_summits.bed \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/group_1/pseudo_bulk_bam/macs2/group_1_summits.bed \
# group0 \
# group1 \
# /gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/result1/sc_pseudoBulk/diffPeaks \
# hg38 \
# /gcdata/gcproj/fulong/Software/qzhudfci-cutruntools-49ddd2487e1a \
# /homes6/fulong/miniconda3/bin \
# /homes6/fulong/miniconda3/bin \
# /homes6/fulong/miniconda3/bin


# http://homer.ucsd.edu/homer/ngs/mergePeaks.html