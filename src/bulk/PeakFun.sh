
peak_file=$1
the_genome=$2
vis_num=$3
outputdir=$4
rscript_dir=`dirname "$0"`
Help()
{
   # Display Help
   echo "############################################################################"
   echo "It is a script for Gene Oncology analysis of input peaks."
   echo
   echo "Syntax: ./PeakFun.sh peak.bed genome vis_num outputdir"
   echo "options:"
   echo "peak.bed     Peak file with BED format, e.g. cluster2.narrowpeak"
   echo "genome       Genome version used, e.g. hg38"
   echo "vis_num      How many significant entries showed in the resulting barplot, e.g. 20"
   echo "outputdir    The directory for output, e.g. /dir/for/output"
   echo
  echo "############################################################################"

}
Help
echo "-------------------------Begin to work-----------------------------------"
echo "Begin to work"
Rscript $rscript_dir/PeakFun.r $the_genome $peak_file $outputdir $vis_num 

echo "---------------------------End here--------------------------------------"


# cd /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/src/bulk
# chmod +x ./PeakFun.sh
# ./PeakFun.sh /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test/FLI1_fimo.bed hg19 20 /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test

