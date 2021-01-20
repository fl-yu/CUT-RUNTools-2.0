peak_file1=$1
peak_file2=$2
outputdir=$3
rscript_dir=`dirname "$0"`


Help()
{
   # Display Help
   echo "############################################################################" 
   echo "It is a script for Peak overlapping analysis of input peaks."
   echo
   echo "Syntax: ./PeakOverlap.sh peak1.bed peak2.bed outputdir"
   echo "options:"
   echo "peak1.bed     Peak file with BED format, e.g. cluster1.narrowpeak"
   echo "peak2.bed     Peak file with BED format, e.g. cluster2.narrowpeak"
   echo "outputdir    The directory for output, e.g. /dir/for/output"

   echo "The software BEDTools were needed, if it was not installed, please install it first https://bedtools.readthedocs.io/en/latest/content/installation.html"
   echo
   echo "############################################################################" 
}
Help


echo "-------------------------Begin to work-----------------------------------"
echo "Begin to work"
echo "[INFO] It will output how many entries in peak1.bed overlap with peak2.bed"

bedtools intersect -a $peak_file1 \
                    -b $peak_file2 \
                    -wa \
                    > $outputdir/overlap_peak1.bed
bedtools intersect -a $peak_file1 \
                    -b $peak_file2 \
                    -v \
                    > $outputdir/nooverlap_peak1.bed
bedtools intersect -a $peak_file2 \
                    -b $peak_file1 \
                    -v \
                    > $outputdir/nooverlap_peak2.bed

peak_num1=`wc -l < $outputdir/nooverlap_peak1.bed`
peak_num2=`wc -l < $outputdir/nooverlap_peak2.bed`
peak_num3=`wc -l < $outputdir/overlap_peak1.bed`
echo "[INFO] Overlapping peak number: $peak_num3"
Rscript $rscript_dir/PeakOverlap.r $peak_num1 $peak_num2 $peak_num3 $outputdir 

echo "---------------------------End here--------------------------------------"

# cd /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/src/bulk
# chmod +x ./PeakOverlap.sh
# ./PeakOverlap.sh /Users/fyu/Documents/GitHub/mpn-gwas/data/atac/panHeme.bed /Users/fyu/Documents/GitHub/mpn-gwas/data/annotations/Coding_UCSC.bed /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test


