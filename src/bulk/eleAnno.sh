peak_file=$1
annotation_dir=$2
outputdir=$3
# rscript_dir=`dirname "$0"`
SCRIPT=`echo "$(cd "$(dirname "$0")" && pwd -P)/$(basename "$0")"`
rscript_dir=`dirname $SCRIPT`

Help()
{
   # Display Help
   echo "############################################################################" 
   echo "It is a script for genomic annotation of input peaks."
   echo
   echo "Syntax: ./eleAnno.sh peak_file annotation_dir outputdir"
   echo "options:"
   echo "peak_file          Peak file with BED format, e.g. cluster1.narrowpeak"
   echo "annotation_dir     A directory contains predifined genomic element files, e.g. /dir/for/ele"
   echo "outputdir          The directory for output, e.g. /dir/for/output"

   echo "The software BEDTools were needed, if it was not installed, please install it first https://bedtools.readthedocs.io/en/latest/content/installation.html"
   echo
   echo "############################################################################" 
}
Help


echo "-------------------------Begin to work-----------------------------------"
echo "Begin to work"
echo "[INFO] Finding the overlap between peaks and genomic elements of 5'UTR, Promoter, Exon, Intron, 3'UTR, intragenic, intergenic, respectively"


cd $outputdir
# elementdir=/gcdata/gcproj/fulong/Data/Genomes/hg38_genomic_element
ele="utr5 promoter exon intron utr3 intragenic intergenic"
# peak_file=/gcdata/gcproj/fulong/Data/scCUT_TAG/pipeline_test_k27e3/0708-result/sc_pseudoBulk/group_1/pseudo_bulk_bam/macs2/narrowPeak_fd5_q10_sort_filter_5k2
for i in $ele
do
    echo  [ok] $i
    bedtools intersect -a $peak_file -b $annotation_dir/${i}.bed -wa > peak_overlap_${i}.bed
done

rm -rf *.temp.count
echo [INFO] Counting the peak number..
wc -l < $peak_file > peak.temp.count

echo [INFO] Counting the overlap..

for i in $ele
do
    wc -l < peak_overlap_${i}.bed >> peak_ele_overlap.temp.count
done


Rscript $rscript_dir/eleAnno.r peak_ele_overlap.temp.count peak.temp.count $outputdir 

echo "---------------------------End here--------------------------------------"

# cd /Users/fulongyu/github/CUT-RUNTools-2.0/src/bulk
# chmod +x *
# ./eleAnno.sh /Users/fulongyu/github/CUT-RUNTools-2.0/test/FLI1_fimo.bed /Users/fulongyu/github/CUT-RUNTools-2.0/test/hg38_genomic_element /Users/fulongyu/github/CUT-RUNTools-2.0/test
