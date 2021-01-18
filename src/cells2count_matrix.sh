# ==================================================================
# step2 cell->count_matrix (feature-by-cell matrix construction)
# ==================================================================
# define all the parameters in this step, useful to change when use this script separatedly

peak_caller=$peak_caller
bin_size=$bin_size # used when activate matrix_type=bin_by_cell
chrome_sizes_file=$chrome_sizes_file # used when activate matrix_type=bin_by_cell
genome=$genome # used when activate matrix_type=bin_by_cell
featureFile=$featureFile # if customFeature_by_cell, a bed file (such as peaks from bulk experiment) to specify the custom features
bedtoolsbin=$bedtoolsbin
cores=$cores
blacklist=`dirname $scriptdir`/blacklist/$genome.blacklist.bed
path_parallel=$path_parallel
scriptdir=$scriptdir
bash_function_dir=$scriptdir/BASHscript
. $bash_function_dir/paraCountMat.sh
matrix_type=$matrix_type # one of "peak_by_cell", "bin_by_cell", "customFeature_by_cell"

if [ "$entire_pipeline" == "TRUE" ]
then
    workdir=$workdir
    sc_countMatrix_dir=$workdir/sc_countMatrix # can make this as an input variable
    logdir=$workdir/sc_logs
    bamdir=$workdir/sc_aligned.aug10/dup.marked.clean
    sc_pseudoBulk_dir=$workdir/sc_pseudoBulk
    if [[ $peak_caller = "macs2" ]]
    then
        peakfile=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/macs2/groups_aggregation_peaks.narrowPeak
    elif [[ $peak_caller = "SEACR" ]]
    then
        peakfile=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/SEACR/groups_aggregation.stringent.bed
    else
    echo "The \$peak_caller parameter is not set as \"macs2\" or \"SEACR\", please check it! Use \"macs2\" for the following processing! "
        peakfile=$sc_pseudoBulk_dir/groups_aggregation/pseudo_bulk_data/macs2/groups_aggregation_peaks.narrowPeak
    fi
    qc_fail=$workdir/sc_qc/report/statistics_QCfailed.txt
    qc_pass=$workdir/sc_qc/report/statistics_QCpassed.txt

    echo "# "
    echo "#         step2 cell->count_matrix (feature-by-cell matrix construction)"
    echo "# "
    date
elif [ "$entire_pipeline" == "FALSE" ] && [ "$individual_step" == "cells2count_matrix" ]
then
    bamdir=$step2_bamfile_dir
    sc_countMatrix_dir=$step2_output_dir
    qc_pass=$step2_qc_pass_file
    sc_pseudoBulk_dir=$workdir/sc_pseudoBulk
    # one of "bin_by_cell", "customFeature_by_cell"
    if [ "$matrix_type" == "peak_by_cell" ]
    then
        echo "[ERROR] matrix_type can only be specified as one of bin_by_cell and customFeature_by_cell when run step2 separatedly"
        exit 1
    fi
    echo "# "
    echo "#         step2 cell->count_matrix (feature-by-cell matrix construction)"
    echo "# "
    echo "[info] Step2 will run separately... "
    date
else 
    echo "[ERROR] Please check your configuration!"
    exit 1
fi

mkdir -p $sc_countMatrix_dir
cd $bamdir
echo "[info] Filter the barcode cells"
if [ -a $qc_pass ]; then
    echo "[info] Cells will be filtered according to file: $qc_pass"
    # remember to delete the header
    bamfile=`cat $qc_pass | sed 1d | cut -f 1 | awk '{print $0".bam"}' | tr "\n" " "`
    else 
    echo "[info] The QCpass file not be found, all the cells will be used for count matrix construction"
    bamfile=`ls *.bam`
    fi

# check the parameter matrix_type
echo "[info] Check the feature type for matrix generation"
if [ "$matrix_type" != "bin_by_cell" ] && [ "$matrix_type" != "customFeature_by_cell" ] && [ "$matrix_type" != "peak_by_cell" ]; then
#if [[ "$matrix_type" != "bin_by_cell" || "$matrix_type" != "customFeature_by_cell" ]]; then
    echo "Parameter matrix_type should be specified as peak_by_cell, bin_by_cell or customFeature_by_cell ... "
    exit
    fi

# "peak_by_cell"
if [ "$matrix_type" = "peak_by_cell" ]; then
    echo "[[ $matrix_type ]] mode activate... "
    echo "[info] Peak file $peakfile will be loaded and processed "
    peakNum=`cat $peakfile | wc -l| cut -d' ' -f1`
    echo "[info] $peakNum peaks were found "
    if [ ! -f $peakfile ]; then
    echo "[info] File $peakfile not found!"
    exit
    fi
    peakfileName=`basename $peakfile`
    echo "[info] Remove peaks overlapped with ENCODE blacklist regions, mitochondrial chromosomes and other uninterest chromosomes ... "
    cat $peakfile |  egrep "^[chr0-9XY$]" | awk '$1 !~ "_"{ print $0}' |  grep -v chrM | $bedopsbin/sort-bed - | $bedtoolsbin/bedtools intersect -wa -v -a - -b $blacklist > "$sc_countMatrix_dir/$peakfileName.filter"
    # awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "+"}' "$sc_countMatrix_dir/$peakfileName.filter" > $sc_countMatrix_dir/$peakfileName".saf"
    echo "[info] Begin to count the reads across pseudo-bulk peaks..."
    peakNum2=`cat "$sc_countMatrix_dir/$peakfileName.filter" | wc -l | cut -d' ' -f1`
    echo "[info] The filtered peak num is $peakNum2..."
    echo "[info] You can find the log files in the directory $logdir..."
    date
    # cd $bamdir
    # $path_featureCount/featureCounts -p -T $cores -F SAF -a $sc_countMatrix_dir/$peakfileName".saf" -o $sc_countMatrix_dir/$matrix_type"_matrix.txt" $bamfile >> $logdir/featureCounts.log 2>&1

    countMat $sc_countMatrix_dir/$peakfileName.filter \
             $sc_countMatrix_dir \
             $sc_countMatrix_dir/$matrix_type"_matrix.txt" \
             $bamdir \
             $cores \
             $bedtoolsbin\
             $path_parallel \
             $bamfile
    fi

# "bin_by_cell"
if [ "$matrix_type" = "bin_by_cell" ] ; then
    echo "[[ $matrix_type ]] mode activate... "
    echo "[info] Begin to process the genome-wide bin files (resolution $bin_size)... "
    echo "[info] The genome-wide bin file will be generated and processed: "${sc_countMatrix_dir}/${genome}"_window_"$bin_size".bed"
    $bedtoolsbin/windowMaker -g $chrome_sizes_file -w $bin_size > $sc_countMatrix_dir/${genome}_window_${bin_size}.bed
    bin_num=`wc -l $sc_countMatrix_dir/${genome}_window_${bin_size}.bed | cut -d' ' -f1`
    echo "[info] $bin_num bins were found "
    echo "[info] Remove bins overlapped with ENCODE blacklist regions, mitochondrial chromosomes and other uninterest chromosomes ... "
    cat $sc_countMatrix_dir/$genome"_window_"$bin_size".bed" |  egrep "^[chr0-9XY$]" | awk '$1 !~ "_"{ print $0}' |  grep -v chrM | $bedopsbin/sort-bed - | $bedtoolsbin/bedtools intersect -wa -v -a - -b $blacklist > $sc_countMatrix_dir/$genome"_window_"$bin_size".bed.filter"
    bin_num2=`wc -l $sc_countMatrix_dir/${genome}_window_${bin_size}.bed.filter | cut -d' ' -f1`
    echo "[info] After filtration, the bin num is $bin_num2..."
    echo "[info] Begin to count the reads across all the bins"
    # echo "[info] You can find the log file in the directory $logdir..."
    date
    # $path_featureCount/featureCounts -p -T $cores -F SAF -a $sc_countMatrix_dir/$genome"_window_"$bin_size".saf" -o $sc_countMatrix_dir/$matrix_type"_matrix.txt" $bamfile >> $logdir/featureCounts.log 2>&1
    countMat $sc_countMatrix_dir/$genome"_window_"$bin_size".bed.filter" \
            $sc_countMatrix_dir \
            $sc_countMatrix_dir/$matrix_type"_matrix.txt" \
            $bamdir \
            $cores \
            $bedtoolsbin \
            $path_parallel \
            $bamfile 
    fi

# "customFeature_by_cell"
if [ "$matrix_type" = "customFeature_by_cell" ]; then
    echo "[[ $matrix_type ]] mode activate... "
    echo "The customFeature file $featureFile will be loaded and processed "
    echo "`cat "$featureFile" | wc -l| cut -d' ' -f1` features were found "
    if [ ! -f $featureFile ]; then
    echo "File $featureFile not found!"
    exit
    fi
    featureFileName=`basename $featureFile`
    echo "[info] Remove features overlapped with ENCODE blacklist regions, mitochondrial chromosomes and other uninterest chromosomes ... "
    cat $featureFile |  egrep "^[chr0-9XY$]" | awk '$1 !~ "_"{ print $0}' |  grep -v chrM | $bedopsbin/sort-bed - | $bedtoolsbin/bedtools intersect -wa -v -a - -b $blacklist > "$sc_countMatrix_dir/$featureFileName.filter"
    echo "[info] The filtered customFeature num is `cat "$sc_countMatrix_dir/$featureFileName.filter" | wc -l| cut -d' ' -f1`..."
    # awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "+"}' "$sc_countMatrix_dir/$featureFileName.filter" > $sc_countMatrix_dir/$featureFileName".saf"
    echo "[info] Begin to count the reads across the custom features..."
    echo "[info] You can find the log files in the directory $logdir..."
    date
    # cd $bamdir
    # $path_featureCount/featureCounts -p -T $cores -F SAF -a $sc_countMatrix_dir/$featureFileName".saf" -o $sc_countMatrix_dir/$matrix_type"_matrix.txt" $bamfile >> $logdir/featureCounts.log 2>&1
    
    countMat $sc_countMatrix_dir/$featureFileName.filter \
            $sc_countMatrix_dir \
            $sc_countMatrix_dir/$matrix_type"_matrix.txt" \
            $bamdir \
            $cores \
            $bedtoolsbin \
            $path_parallel \
            $bamfile 
    fi

echo "# "
echo "#         step 2 (feature-by-cell matrix construction) is completed"
echo "# "
echo "-------------------------------------------------------------------------------------------------------------------------"
