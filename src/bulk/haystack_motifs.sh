peak_file=$1
the_genome=$2
outputdir=$3
rscript_dir=`dirname "$0"`

echo "############################################################################" 
echo "It is a script for TF motif discovery for the input peaks."
echo
echo "Syntax: ./haystack_motif.sh peak.bed genome outputdir"
echo "options:"
echo "peak1.bed     Peak file with BED format, e.g. cluster1.narrowpeak"
echo "genome        Genome version used, e.g. hg38"
echo "outputdir     The directory for output, e.g. /dir/for/output"

echo "The software haystack were needed, if it was not installed, please install it first https://github.com/pinellolab/haystack_bio"
echo
echo "############################################################################" 


echo "-------------------------Begin to work-----------------------------------"
echo "Begin to work"
echo haystack_motifs $peak_file  $the_genome --name haystack_motif --output_directory $outputdir
haystack_motifs $peak_file  $the_genome --name haystack_motif --output_directory $outputdir
# cd /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/src/bulk
# chmod +x ./haystack_motifs.sh
# ./haystack_motifs.sh /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test/FLI1_fimo.bed hg19 /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test
