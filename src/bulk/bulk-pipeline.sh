#!/bin/bash
# Copyright (C) 2020 Fulong Yu
#
# CUT&RUNTools 2.0 is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools 2.0 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.


infile=$fastq_directory/$sample_name
#expand the path of infile
# relinfile=`realpath -s $infile`
dirname=`dirname $infile`
# base=`basename $infile _R1_001.fastq.gz`
base=$sample_name
>&2 echo "[info] Input file is `basename $infile`_R1_001.fastq.gz and `basename $infile`_R2_001.fastq.gz"
>&2 date

#cd to current directory
cd $dirname
workdir=$workdir
mkdir -p $workdir
len=$fastq_sequence_length
trimdir=$workdir/trimmed
trimdir2=$workdir/trimmed2
logdir=$workdir/logs
aligndir=$workdir/aligned
peakdir=$workdir/peakcalling

for d in $trimdir $trimdir2 $logdir $aligndir $peakdir; do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done

#trimming paired-end
#good version
>&2 echo "[info] Trimming file $base ..."
>&2 date
$javabin/java -jar $trimmomaticbin/$trimmomaticjarfile PE -threads 1 -phred33 $dirname/"$base"_R1_001.fastq.gz $dirname/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25 > $logdir/${base}_trim1.log 2>&1

>&2 echo "[info] Second stage trimming $base ..."
>&2 date
$kseqbin/kseq_test $trimdir/"$base"_1.paired.fastq.gz $len $trimdir2/"$base"_1.paired.fastq.gz
$kseqbin/kseq_test $trimdir/"$base"_2.paired.fastq.gz $len $trimdir2/"$base"_2.paired.fastq.gz

>&2 echo "[info] Aligning file $base to reference genome..."
>&2 date
if [ "$frag_120" == "TRUE" ]
then
    echo "[info] Bowtie2 command: --dovetail --phred33"
    echo "[info] The dovetail mode is enabled [as parameter frag_120 is on]"
    ($bowtie2bin/bowtie2 -p $cores --dovetail --phred33 -x $bt2idx/genome -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/"$base".bam
else
    echo "[info] Bowtie2 command: --very-sensitive-local --phred33 -I 10 -X 700"
    echo "[info] The dovetail mode is off [as parameter frag_120 is off]"
    ($bowtie2bin/bowtie2 -p $cores --very-sensitive-local --phred33 -I 10 -X 700 -x $bt2idx/genome -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/"$base".bam
fi

spikein_dir=$workdir/

# check the parameters
if [ "$spike_in" == "TRUE" ]
then
    mkdir -p $spikein_dir
    >&2 echo "[info] Aligning file $base to spike-in genome"
    >&2 date
    ($bowtie2bin/bowtie2 -p $cores --dovetail --phred33 -x $spike_in_bt2idx/genome -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".spikein.bowtie2 | $samtoolsbin/samtools view -bS - > $spikein_dir/"$base".bam

    $logdir/"$base".spikein.bowtie2
    total_reads=`cat $logdir/"$base".spikein.bowtie2 | grep "reads; of these:" | awk '{print $1}' - FS=' '`
    align_ratio=`cat $logdir/"$base".spikein.bowtie2 | grep "overall alignment" | awk '{print $1}' - FS=' ' | cut -f1 -d"%"`
    spikein_reads=`printf "%.0f" $(echo "$total_reads * $align_ratio"|bc)`

    >&2 echo "[info] Spikein reads number is $spikein_reads, consisting of $align_ratio % of total reads"
    >&2 echo "[info] This information could be used in spike-in normalization when generating bigwig files"
else
    >&2 echo "[info] FASTQ files won't be aligned to the spike-in genome"
fi


dirname=$aligndir
#cd to current directory (aligned)
cd $dirname

# workdir=`pwd`
# logdir=$workdir/logs

for d in $logdir sorted dup.marked dedup; do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done

>&2 echo "[info] Filtering unmapped fragments... ""$base".bam
>&2 date
$samtoolsbin/samtools view -bh -f 3 -F 4 -F 8 $dirname/"$base".bam > sorted/"$base".step1.bam

>&2 echo "[info] Sorting BAM... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile SortSam \
INPUT=sorted/"$base".step1.bam OUTPUT=sorted/"$base".bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
rm -rf sorted/"$base".step1.bam

>&2 echo "[info] Marking duplicates... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dup.marked/"$base".bam VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=metrics."$base".txt 2> $logdir/"$base".mark.dup

>&2 echo "[info] Removing duplicates... ""$base".bam
>&2 date
$samtoolsbin/samtools view -bh -F 1024 dup.marked/"$base".bam > dedup/"$base".bam


if [ "$frag_120" == "TRUE" ]
then
    for d in dup.marked.120bp dedup.120bp; do
        if [ ! -d $d ]; then
            mkdir $d
        fi
    done
    >&2 echo "[info] Filtering to <120bp... dup.marked and dedup BAMs"
    >&2 date
    $samtoolsbin/samtools view -h dup.marked/"$base".bam | LC_ALL=C awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dup.marked.120bp/"$base".bam
    $samtoolsbin/samtools view -h dedup/"$base".bam | LC_ALL=C awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dedup.120bp/"$base".bam
    dupdir=$dirname/dup.marked.120bp
    dedupdir=$dirname/dedup.120bp
else

    >&2 echo "[info] Using all the qualified fragments NOT filtering <120bp... ""$base".bam
    >&2 date
    dupdir=$dirname/dup.marked
    dedupdir=$dirname/dedup
fi

>&2 echo "[info] Creating bam index files... ""$base".bam
>&2 date
$samtoolsbin/samtools index $dupdir/"$base".bam
$samtoolsbin/samtools index $dedupdir/"$base".bam


# shifting reads if the experiment is CUT&Tag
if [ "$experiment_type" == "CUT&Tag" ]
then
    >&2 echo "[info] Reads shifting "
    >&2 date
    >&2 echo "[info] All reads aligning to the + strand were offset by +4 bp, and all reads aligning to the – strand were offset −5 bp, since the CUT&Tag approach uses Tn5 transposase which has been shown to bind as a dimer and insert two adaptors separated by 9 bp... "
    >&2 echo "[info] The shifted BAM files will be used in the following analysis "
    for d in ${dupdir}.shift ${dedupdir}.shift; do
        if [ ! -d $d ]; then
            mkdir $d
        fi
    done
    >&2 echo "[info] Shifting and indexing "
    # use --ATACshift
    $path_deeptools/alignmentSieve --numberOfProcessors $cores --ATACshift --bam $dupdir/"$base".bam -o ${dupdir}.shift/"$base".tmp.bam
    $path_deeptools/alignmentSieve --numberOfProcessors $cores --ATACshift --bam $dedupdir/"$base".bam -o ${dedupdir}.shift/"$base".tmp.bam 
    # the bam file needs to be sorted again
    $samtoolsbin/samtools sort -@ $cores -O bam -o ${dupdir}.shift/"$base".bam ${dupdir}.shift/"$base".tmp.bam
    $samtoolsbin/samtools index -@ $cores ${dupdir}.shift/"$base".bam
    rm ${dupdir}.shift/"$base".tmp.bam
    $samtoolsbin/samtools sort -@ $cores -O bam -o ${dedupdir}.shift/"$base".bam ${dedupdir}.shift/"$base".tmp.bam
    $samtoolsbin/samtools index -@ $cores ${dedupdir}.shift/"$base".bam
    rm ${dedupdir}.shift/"$base".tmp.bam
    # replace the variables
    dupdir=${dupdir}.shift
    dedupdir=${dedupdir}.shift
elif [ "$experiment_type" == "CUT&RUN" ]
then
    >&2 echo "[info] Reads shifting "
    >&2 date
    >&2 echo "[info] Your data won't be shifted as the experiment_type is specified as CUT&RUN... "
else
    >&2 echo "[info] Reads shifting "
    >&2 echo date
    >&2 echo "[warning] The parameter of experiment_type is not CUT&RUN or CUT&Tag..."
fi


>&2 echo "[info] Peak calling using MACS2... ""$base".bam
>&2 echo "[info] Logs are stored in $logdir"
>&2 date
if [ "$dup_peak_calling" == "TRUE" ]
then
    echo "[info] Peak calling with BAM file with duplications"
    bam_file=$dupdir/"$base".bam
else
    echo "[info] Peak calling with BAM file with NO duplications"
    bam_file=$dedupdir/"$base".bam
fi

dir=`dirname $bam_file`
base_file=`basename $bam_file .bam`

outdir=$workdir/peakcalling/macs2.narrow #for macs2
outdirbroad=$workdir/peakcalling/macs2.broad #for macs2
outdirseac=$workdir/peakcalling/seacr #for seacr

for d in $outdir $outdirbroad $outdirseac; do
    mkdir -p $d
done
# macs2 narrow peak calling
>&2 echo "[info] macs2 narrow peak calling"
$macs2bin/macs2 callpeak -t $bam_file -g $macs2_genome -f BAMPE -n $base_file --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2.narrow
# macs2 broad peak calling
>&2 echo "[info] macs2 broad peak calling"
$macs2bin/macs2 callpeak -t $bam_file -g $macs2_genome -f BAMPE -n $base_file --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2.broad
# macs2 broad peak calling
>&2 echo "[info] Getting broad peak summits"
$pythonbin/python $extratoolsbin/get_summits_broadPeak.py $outdirbroad/"$base_file"_peaks.broadPeak | $bedopsbin/sort-bed - > $outdirbroad/"$base_file"_summits.bed
#SEACR peak calls
>&2 echo "[info] SEACR stringent peak calling"
$macs2bin/macs2 callpeak -t $bam_file -g $macs2_genome -f BAMPE -n $base_file --outdir $outdirseac -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".seacr
$pythonbin/python $extratoolsbin/change.bdg.py $outdirseac/"$base_file"_treat_pileup.bdg > $outdirseac/"$base_file"_treat_integer.bdg
chmod +x $extratoolsbin/SEACR_1.1.sh
$extratoolsbin/SEACR_1.1.sh $outdirseac/"$base_file"_treat_integer.bdg 0.01 non stringent $outdirseac/"$base_file"_treat $Rscriptbin
$bedopsbin/sort-bed $outdirseac/"$base_file"_treat.stringent.bed > $outdirseac/"$base_file"_treat.stringent.sort.bed
$pythonbin/python $extratoolsbin/get_summits_seacr.py $outdirseac/"$base_file"_treat.stringent.bed|$bedopsbin/sort-bed - > $outdirseac/"$base_file"_treat.stringent.sort.summits.bed
for i in _summits.bed _peaks.xls _peaks.narrowPeak _control_lambda.bdg _treat_pileup.bdg; do 
    rm -rf $outdirseac/"$base_file"$i
done


>&2 echo "[info] Generating the normalized signal file with BigWig format... "
>&2 date
cd $outdir

# LC_ALL=C sort -k1,1 -k2,2n $outdir/"$base_file"_treat_pileup.bdg > $outdir/"$base_file".sort.bdg
# $extratoolsbin/bedGraphToBigWig $outdir/"$base_file".sort.bdg $chrom_size_file $outdir/"$base_file".sorted.bw # use deeptools
#     rm -rf $outdir/"$base_file".sort.bdg
$path_deeptools/bamCoverage --bam $bam_file -o $outdir/"$base_file".cpm.norm.bw \
    --binSize 10 \
    --normalizeUsing CPM \
    --effectiveGenomeSize $eGenomeSize \
    --numberOfProcessors $cores 2> $logdir/"$base".gene.bw
cp $outdir/"$base_file".cpm.norm.bw $outdirbroad
cp $outdir/"$base_file".cpm.norm.bw $outdirseac
if [ "$spike_in" == "TRUE" ]
then 
    if [ "$spikein_reads" == "0" ] || [ "$spike_in_norm" == "FALSE" ]
    then
        >&2 echo "[info] Your bigwig file won't be normalized with spike-in reads as you did not specify this parameter or the spike-in reads were 0.."
    else
        >&2 echo "[info] Your bigwig file will be normalized with spike-in reads"
        scale=$spikein_scale
        scale_factor=`printf "%.0f" $(echo "$scale / $spikein_reads"|bc)`
        >&2 echo scale_factor=$scale_factor
        $path_deeptools/bamCoverage --bam $bam_file -o $outdir/"$base_file".spikein_normalized.bw \
        --binSize 10 \
        --normalizeUsing CPM \
        --effectiveGenomeSize $eGenomeSize \
        --scaleFactor $scale_factor
        cp $outdir/"$base_file".spikein_normalized.bw $outdirbroad
        cp $outdir/"$base_file".spikein_normalized.bw $outdirseac
    fi
else
    >&2 echo "[info] Your bigwig file won't be normalized with spike-in reads"
fi


# 
# specify peak file for the motif and footprinting analysis
# 
if [[ "$peak_caller" == "macs2" ]]
then
    # using narrow by default
    peak_file=$outdir/"$base_file"_peaks.narrowPeak
    summit=$outdir/"$base_file"_summits.bed
    suffix="_peaks.narrowPeak"
    summit_suffix="_summits.bed"
    summit_padded_suffix="_summits_padded.fa"
elif [[ "$peak_caller" == "SEACR" ]]
then
    peak_file=$outdirseac/"$base_file"_treat.stringent.sort.bed
    summit=$outdirseac/"$base_file"_treat.stringent.sort.summits.bed
    suffix="_treat.stringent.sort.bed"
    summit_suffix="_treat.stringent.sort.summits.bed"
    summit_padded_suffix="_treat.stringent.sort.summits_padded.fa"
else
    echo "[info] The \$peak_caller parameter is not set as \"macs2\" or \"SEACR\", please check it! Use \"macs2\" for the following processing! "
    peak_file=$outdir/"$base_file"_peaks.narrowPeak
    summit=$outdir/"$base_file"_summits.bed
    suffix="_peaks.narrowPeak"
    summit_suffix="_summits.bed"
    summit_padded_suffix="_summits_padded.fa"
fi

#------------------------------------------------------------------------------------------------------------------------
#
# motif discovery
#
peak_file=$peak_file #filename must end with .narrowPeak or .bed (if SEACR)
>&2 echo "[info] Input file is $peak_file"

#expand the path for $sample_name
# relinfile=`realpath -s $i`
dirname=`dirname $peak_file`

#cd to current directory
cd $dirname

for d in blacklist_filtered; do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done

fname=$base_file
peak=$fname$suffix
summit=$fname$summit_suffix
summitfa=$fname$summit_padded_suffix

cat $peak | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blacklist_filtered/$peak
cat $summit | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blacklist_filtered/$summit

motif_dir=random$num_peaks
msummit=$motif_dir/summits
mpadded=$motif_dir/padded
mpaddedfa=$motif_dir/padded.fa
for d in $motif_dir $msummit $mpadded $mpaddedfa; do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done

>&2 echo "[info] Get randomized [$num_peaks] peaks from the top [$total_peaks] peaks..."
>&2 echo "[info] Filtering the blacklist regions for the selected peak files"
if [[ "$peak_caller" == "macs2" ]]
then
    cat blacklist_filtered/$peak | sort -t$'\t' -g -k8 -r | head -n $total_peaks | shuf | head -n $num_peaks | $bedopsbin/sort-bed - > $motif_dir/$peak
else
    cat blacklist_filtered/$peak | sort -t$'\t' -g -k5 -r | head -n $total_peaks | shuf | head -n $num_peaks | $bedopsbin/sort-bed - > $motif_dir/$peak
fi

$bedopsbin/bedops -e 1 blacklist_filtered/$summit $motif_dir/$peak > $msummit/$summit
$bedopsbin/bedops --range $num_bp_from_summit -u $msummit/$summit > $mpadded/$summit
>&2 echo "[info] Getting Fasta sequences"
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $mpadded/$summit -fo $mpaddedfa/$summitfa

>&2 echo "[info] Start MEME analysis for de novo motif finding ..."
>&2 echo "[info] Up to $num_motifs will be output ..."
meme_outdir=$motif_dir/MEME_"$fname"_shuf
$memebin/meme-chip -oc $meme_outdir -dreme-m $num_motifs -meme-nmotifs $num_motifs $mpaddedfa/$summitfa
>&2 echo "[info] De Novo motifs can be found: $meme_outdir ..."

#------------------------------------------------------------------------------------------------------------------------
#
# footprinting for De novo motifs
#
cd $dirname
mbase=$fname
mdiscovery=$meme_outdir

>&2 echo "[info] Loading the De Novo motifs ..."
$pythonbin/python $extrasettings/read.meme.py $mdiscovery
p=$motif_scanning_pval
>&2 echo "[info] The signficance cutoff of Fimo scaning is ${p}..."

motif_dir=$mdiscovery/motifs # a directory containing a list of *.meme files
>&2 echo "[info] Motif files can be found: $motif_dir"
peak_filename=`basename $peak_file`
workdir=`pwd`
dir=blacklist_filtered
fa_dir=blacklist_filtered.fa

if [ ! -d $fa_dir ]; then
    mkdir -p $fa_dir
fi
>&2 echo "[info] Filtering the blacklist regions for the selected peak files"
if [ ! -d $dir ] || [ ! -f $dir/$peak ] ; then
    cat $workdir/$dir/"$peak_filename" | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > $workdir/$dir/$peak
fi
>&2 echo "[info] Getting Fasta sequences"
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/$peak -fo $fa_dir/"$mbase".fa
$pythonbin/python $extrasettings/fix_sequence.py $fa_dir/"$mbase".fa

outdir=fimo.result
for d in $outdir $outdir/$mbase; do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done
>&2 echo "[info] Scaning the De Novo motifs for each peak"

for m in `ls -1 $motif_dir`; do
    motif=`basename $m .meme`
    fimo_d=$outdir/$mbase/fimo2.$motif
    if [ ! -d $fimo_d ]; then
        mkdir $fimo_d
    fi
    $memebin/fimo --thresh $p --parse-genomic-coord -oc $fimo_d $motif_dir/"$motif".meme $fa_dir/"$mbase".fa
    $bedopsbin/gff2bed < $fimo_d/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $fimo_d/fimo.bed
done
>&2 echo "[info] Output can be found: $outdir/$mbase"


####################################################################

workdir=`pwd`
dir=`dirname $bam_file`
bambase=`basename $bam_file .bam`
outbam=$bam_file

fimo_dir=$outdir/$mbase

for i in `ls -1 $fimo_dir`; do #shows a list of motifs
    echo "[info] Doing $i..."
    fimo_d=$fimo_dir/$i
    tmp=`echo $i|cut -d "." -f3 | wc -c`
    mlen=$(( tmp - 1 ))
    $makecutmatrixbin/make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p 1 -f 3 -F 4 -F 8 -q 0 $outbam $fimo_d/fimo.bed > $fimo_d/fimo.cuts.freq.txt
    $Rscriptbin/Rscript $extratoolsbin/run_centipede_parker.R $fimo_d/fimo.cuts.freq.txt $fimo_d/fimo.bed $fimo_d/fimo.pdf $mlen
done

echo "# "
echo "#        Congrats! The bulk data analysis is complete!"
echo "# "
echo "--------------------------------------------------------------------------"

