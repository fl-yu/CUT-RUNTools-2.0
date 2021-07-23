# cd /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/src/bulk
# chmod +x ./get_cuts_single_locus.sh
# ./get_cuts_single_locus.sh scriptdir/bulk-config.json /cutruntools/yoursample.bam chr2:60457194-60553654 /my/singlelocus/dir

configrue_file=$1 # $scriptdir/bulk-config.json
bam_file=$2
roi=$3
outputdir=$4

echo "############################################################################" 
echo "This is a script for a single nucleotide resolution cut profile for a region of interest/the whole genome."
echo "     "
echo "Basically,  we need the configuration file (tell the code where the necessary software is), a processed BAM file of CUT&RUN/CUT&Tag, a region of interest (roi; in the format chr2:60457194-60553654), and an output directory"
echo
echo "Syntax: ./get_cuts_single_locus.sh configuration_file bam_file roi outputdir"
echo "options:"
echo "configrue_file     The configuration file, you may use one for running the bulk pipeline. Except for the software info, experiment type (cut&run or cut&tag) will be also used"
echo "bam_file               a processed BAM file from bulk pipeline"
echo "roi                    a region of interest (roi; in the format chr2:60457194-60553654). If you want to generate the whole genome profile instead of a particular genomic locus, please specify this parameter as 'whole_genome'"
echo "outputdir              The directory for output, e.g. /dir/for/output"
echo "     "
echo "Example usage:"
echo "      cd /Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/src/bulk"
echo "      chmod +x ./get_cuts_single_locus.sh"
echo "      ./get_cuts_single_locus.sh scriptdir/bulk-config.json /cutruntools/mysample.bam chr2:60457194-60553654 /my/singlelocus/dir"
echo "     "
echo "If you used CUT&RUNTools 2.0 in your work, please cite our paper (PMID:34244724)"
echo "If you run into issues and would like to report them, you can use the 'Issues' tab on our github site (https://github.com/fl-yu/CUT-RUNTools-2.0). Alternatively, you can contact authors: fyu@broadinstitute.org or guo-cheng.yuan@mssm.edu"
echo
echo "############################################################################" 

SCRIPT=`echo "$(cd "$(dirname "$0")" && pwd -P)/$(basename "$0")"`
scriptdir=`dirname $SCRIPT`
SCRIPTPATH=`echo ${scriptdir%/*/*}`
# convert configuration JSON files to bash variables
eval "$(jq -r '.software_config | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
eval "$(jq -r '.input_output | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
eval "$(jq -r '.motif_finding | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
chrom_size_file=$SCRIPTPATH/assemblies/chrom."$organism_build"/"$organism_build".chrom.sizes

echo "-------------------------Begin to work-----------------------------------"
echo "Begin to work"

sample_name=get_cut_profiles_"$roi"

# The R1 and R2 bigwigs designate strand-specific cut profiles that were created. The frag.ends.bw is one that combines cuts from both strands. The bigwigs can be displayed in any visualization tools such as IGV or UCSC genome browser.
if [ "$roi" == "whole_genome" ]
then
    echo "genomic cut profile of whole genome will be generated"
    cp $bam_file $outputdir/"$sample_name".bam
else
    $samtoolsbin/samtools view -bh $samtoolsflags $bam_file "$roi" > $outputdir/"$sample_name".bam
fi
# cd $outputdir
newbam_file=$outputdir/"$sample_name".bam
sample_temp_name="$sample_name"_temp
$samtoolsbin/samtools index $newbam_file

$samtoolsbin/samtools view -b -f 3 -F 4 -F 8 $newbam_file | 
    $samtoolsbin/samtools sort -O bam -n - | 
    $bedtoolsbin/bedtools bamtobed -i stdin > $outputdir/"$sample_temp_name".fragment.bed

awk -vshift=1 -vOFS="\t" '($6=="+"){ print $1, $2, ($2 + shift), $4, ".", $6}($6=="-"){ print $1, ($3 - shift), $3, $4, ".", $6}' $outputdir/"$sample_temp_name".fragment.bed > $outputdir/"$sample_temp_name".fragment.end.bed 

grep "-" $outputdir/"$sample_temp_name".fragment.end.bed > $outputdir/"$sample_temp_name".fragment.end.minus.bed
grep "+" $outputdir/"$sample_temp_name".fragment.end.bed > $outputdir/"$sample_temp_name".fragment.end.plus.bed

sort -k1,1 -k2,2n $outputdir/"$sample_temp_name".fragment.end.bed > $outputdir/"$sample_temp_name".fragment.end.sorted.bed
sort -k1,1 -k2,2n $outputdir/"$sample_temp_name".fragment.end.minus.bed > $outputdir/"$sample_temp_name".fragment.end.minus.sorted.bed
sort -k1,1 -k2,2n $outputdir/"$sample_temp_name".fragment.end.plus.bed > $outputdir/"$sample_temp_name".fragment.end.plus.sorted.bed

bedtools genomecov -i $outputdir/"$sample_temp_name".fragment.end.sorted.bed -g $chrom_size_file -bg  > $outputdir/"$sample_temp_name".fragment.end.sorted.bdg
bedtools genomecov -i $outputdir/"$sample_temp_name".fragment.end.minus.sorted.bed -g $chrom_size_file -bg  > $outputdir/"$sample_temp_name".fragment.end.minus.sorted.bdg
bedtools genomecov -i $outputdir/"$sample_temp_name".fragment.end.plus.sorted.bed -g $chrom_size_file -bg  > $outputdir/"$sample_temp_name".fragment.end.plus.sorted.bdg

chmod u+x $bdg2bwbin/bedGraphToBigWig
$bdg2bwbin/bedGraphToBigWig $outputdir/"$sample_temp_name".fragment.end.sorted.bdg $chrom_size_file $outputdir/"$sample_name"_combined.bw
$bdg2bwbin/bedGraphToBigWig $outputdir/"$sample_temp_name".fragment.end.minus.sorted.bdg $chrom_size_file $outputdir/"$sample_name"_minusstrand.bw
$bdg2bwbin/bedGraphToBigWig $outputdir/"$sample_temp_name".fragment.end.plus.sorted.bdg $chrom_size_file $outputdir/"$sample_name"_plusstrand.bw

rm $outputdir/*"$sample_temp_name"*
