#need sample.bed in the current directory

#url is defined as following

organism_build=$1
bedtoolsbin=$2

if [ "$organism_build" == "hg38" ]
then
    url=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFaMasked.tar.gz
elif [ "$organism_build" == "hg19" ]
then
    url=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
elif [ "$organism_build" == "mm10" ]
then
    http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFaMasked.tar.gz
elif [ "$organism_build" == "mm10" ]
then
    http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFaMasked.tar.gz
else
    echo "Only support hg38, hg19, mm10 or mm9, please check your parameter.."
fi


cur=`pwd`
mkdir -p assemblies/$organism_build
cd assemblies/$organism_build
wget $url -O chromFaMasked.tar.gz
tar -zxf chromFaMasked.tar.gz
#ls -1 maskedChroms/chr*.fa|xargs cat > hg38.fa
ls -1 chr*.fa.masked | xargs cat > ${organism_build}.fa
#should create *.fai index file
$bedtoolsbin/bedtools getfasta -fi ${organism_build}.fa -bed ../../sample.bed
#rm -rf chr*.fa.masked
cd $cur


echo "Genome assembly ${organism_build} was installed and the corresponding index file is generated.."