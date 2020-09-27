
TFpeakAnalysis=$TFpeakAnalysis
TA_filter_dir=/gcdata/gcproj/fulong/Software/scCutRunTools/blacklist
meme_db=/gcdata/gcproj/fulong/Software/scCutRunTools/meme_db/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

if [ "$TFpeakAnalysis" == "TRUE" ]
then 
    # -----------------------
    # motif analysis
    # revise from script integrate.motif.find
    # define the input variables of this step
    memebin=$memebin
    bedopsbin=$bedopsbin
    bedtoolsbin=$bedtoolsbin
    perlbin=$perlbin
    genome_sequence=$genome_sequence
    extrasettings=$extrasettings
    samtoolsbin=$samtoolsbin
    makecutmatrixbin=$makecutmatrixbin
    Rscriptbin=$Rscriptbin
    pythonbin=$pythonbin
    pythonldlibrary=$pythonldlibrary

    # a new variable to sepcify orignial 3 scripts in the macs2 folder (read.meme.py fix_sequence.py run_centipede_parker.R )
    ori_script_dir=$scriptdir/crtools_script
    

    # blacklist=$extrasettings/hg38.blacklist.bed
    # download mm10 from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
    # confirm the reference genome and blacklist
    genome=$genome # one of hg38 hg19 mm10 mm9
    if [ $genome != "hg38" ] && [ $genome != "hg19" ] && [ $genome != "mm10" ] && [ $genome != "mm9" ]
    then
        echo "[Warnings] Parameter genome should be one of hg38, hg19, mm10 and mm9"
    fi
    blacklist=$extrasettings/$genome.blacklist.bed
    TA_filter_file=$TA_filter_dir/${genome}_TA_repeat.bed


    
    group_name=`sed 1d $cell_anno | cut -f 2 | sort | uniq | tr "\n" " "`

    for i in $group_name
    do
        # peak=$sc_pseudoBulk_dir/group_${i}/pseudo_bulk_bam/macs2/group_${i}_peaks.narrowPeak
        peakfile=$sc_pseudoBulk_dir/group_${i}/pseudo_bulk_bam/SEACR/group_${i}.stringent.bed
        peakdir=`dirname $peakfile`
        peak=`basename $peakfile`
        $pythonbin/python $extratoolsbin/get_summits_seacr.py $peakdir/group_${i}.stringent.bed | $bedopsbin/sort-bed - > $peakdir/group_${i}_summits.bed

        # i=$1 #filename must end with .narrowPeak or .broadPeak or .bed (if SEACR)
        >&2 echo "[info] Input file is $peakfile"

        #expand the path for $1
        relinfile=`realpath -s $peak`
        # dirname=`dirname $relinfile`

        #cd to $i/SEACR ($i/macs2) directory
        cd $peakdir


        for d in blk_filtered; do
        if [ ! -d $d ]; then
        mkdir $d
        fi
        done

        workdir=`pwd`
        fname=group_${i}
        # peak=$fname"_peaks.narrowPeak"
        summit=$fname"_summits.bed"
        summitfa=$fname"_summits_padded.fa"

        >&2 echo "[info] Get filtered peaks..."
        # cat $peak | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist | $bedopsbin/bedops -n 1 - $TA_filter_file > blk_filtered/$peak
        cat $peak | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blk_filtered/$peak
        cat $summit | $bedopsbin/sort-bed - | $bedopsbin/bedops -e 1 - blk_filtered/$peak > blk_filtered/$summit

        #motif discovery starts here
        motif_dir=random.1000
        msummit=$motif_dir/summits
        mpadded=$motif_dir/padded
        mpaddedfa=$motif_dir/padded.fa

        for d in $motif_dir $msummit $mpadded $mpaddedfa; do
        if [ ! -d $d ]; then
        mkdir $d
        fi
        done

        >&2 echo "[info] Get randomized 1000 peaks..."
        # random select 1000 from top 5000 according to fold_enrichment
        cat blk_filtered/$peak | sort -t"	" -g -k8 -r | head -n 5000 | shuf | head -n 1000 | $bedopsbin/sort-bed - > $motif_dir/$peak
        $bedopsbin/bedops -e 1 blk_filtered/$summit $motif_dir/$peak > $msummit/$summit
        # 300 regions from summit 
        $bedopsbin/bedops --range 150 -u $msummit/$summit > $mpadded/$summit
        # get the sequence for new peaks
        $bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $mpadded/$summit -fo $mpaddedfa/$summitfa

        >&2 echo "[info] Start MEME analysis for de novo motif finding..."
        meme_outdir=$motif_dir/MEME_"$fname"_shuf
        # $memebin/meme-chip -oc $meme_outdir -dreme-m 10 -meme-nmotifs 10 -meme-p $cores $mpaddedfa/$summitfa
        $memebin/meme-chip -oc $meme_outdir -dreme-m 10 -meme-nmotifs 10 -db $meme_db $mpaddedfa/$summitfa
        >&2 echo "[info] Finished"
    done

    echo -e "`date` [info] Motif discovery finished \n"



    # -----------------------
    # footprint analysis
    # revise from script integrate.footprinting
    # define the input variables of this step
    echo -e "`date` [info] footprinting analysis begins\n"
    pythonldlibrary=/homes6/fulong/miniconda3/envs/py3/lib
    ldlibrary=`echo $LD_LIBRARY_PATH | tr : "\n" | grep -v $pythonldlibrary | paste -s -d:`
    unset LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$pythonldlibrary:$ldlibrary

    for i in $group_name
    do
        
        # peak_file=$sc_pseudoBulk_dir/group_${i}/pseudo_bulk_bam/macs2/group_${i}_peaks.narrowPeak
        peakfile=$sc_pseudoBulk_dir/group_${i}/pseudo_bulk_bam/SEACR/group_${i}.stringent.bed
        peakdir=`dirname $peakfile`
        peak=`basename $peakfile`
        $pythonbin/python $extratoolsbin/get_summits_seacr.py $peakdir/group_${i}.stringent.bed | $bedopsbin/sort-bed - > $peakdir/group_${i}_summits.bed

        #peak_file=$1 #a narrowPeak/broadPeak/SEACR bed file
        mbase=group_${i}
        # peak=$mbase"_peaks.narrowPeak"
        mdiscovery=random.1000/MEME_"$mbase"_shuf

        #expand the path for $peak_file
        relinfile=`realpath -s $peak_file`
        # dirname=`dirname $relinfile`

        #cd to current directory (macs2.narrow.aug10)
        cd $peakdir
        ## create motif dir and write the motif in
        echo "`date` [info] readme"
        $pythonbin/python $ori_script_dir/read.meme.py $mdiscovery

        p=0.00050
        motif_dir=$mdiscovery/motifs #a directory containing a list of *.meme files
        workdir=`pwd`
        dir=blk_filtered
        fa_dir=blk_filtered.fa

        if [ ! -d $fa_dir ]; then
        mkdir $fa_dir
        fi

        if [ ! -d $dir ] || [ ! -f $dir/$peak ] ; then
        cat $workdir/$dir/"$peak" | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist | $bedopsbin/bedops -n 1 - $TA_filter_file > $workdir/$dir/$peak
        fi


        # get peak sequence
        echo "`date` [info] get peak sequence"
        $bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/$peak -fo $fa_dir/"$mbase".fa
        $pythonbin/python $ori_script_dir/fix_sequence.py $fa_dir/"$mbase".fa

        outdir=fimo.result
        for d in $outdir $outdir/$mbase; do
        if [ ! -d $d ]; then
        mkdir $d
        fi
        done

        for m in `ls -1 $motif_dir`; do
        motif=`basename $m .meme`
        fimo_d=$outdir/$mbase/fimo2.$motif
        if [ ! -d $fimo_d ]; then
        mkdir $fimo_d
        fi
        $memebin/fimo --thresh $p --parse-genomic-coord -oc $fimo_d $motif_dir/"$motif".meme $fa_dir/"$mbase".fa
        # /homes6/fulong/miniconda3/envs/py3/bin/fimo --thresh 0.0005 --parse-genomic-coord -oc $fimo_d $motif_dir/"$motif".meme $fa_dir/"$mbase".fa

        cur_path=`echo $PATH | tr : "\n" | grep -v $bedopsbin | paste -s -d:`
        unset PATH
        export PATH=$cur_path:$bedopsbin

        $bedopsbin/gff2bed < $fimo_d/fimo.gff | awk 'BEGIN {IFS="	"; OFS="	";} {print $1,$2,$3,$4,$5,$6}' > $fimo_d/fimo.bed
        done

        #bamfile=../aligned.aug10/dup.marked.120bp/"$mbase".bam
        bamfile=../"$mbase"_pseudo_sort.bam

        workdir=`pwd`
        dir=`dirname $bamfile`
        bambase=`basename $bamfile .bam`

        dest=centipede.bam
        outbam=$dest/"$bambase".bam
        if [ ! -d $dest ]; then
        mkdir $dest
        fi

        cd $dest
        # ln -s ../../aligned.aug10/dup.marked.120bp/"$mbase".bam .
        # ln -s ../../aligned.aug10/dup.marked.120bp/"$mbase".bam.bai .
        ln -s ../../"$mbase"_pseudo_sort.bam .
        ln -s ../../"$mbase"_pseudo_sort.bam.bai .
        cd ..

        fimo_dir=$outdir/$mbase

        for i in `ls -1 $fimo_dir`; do #shows a list of motifs
        echo "[info] Doing $i..."
        fimo_d=$fimo_dir/$i
        tmp=`echo $i|cut -d "." -f3|wc -c`
        mlen=$(( tmp - 1 ))
        $makecutmatrixbin/make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p $cores -f 3 -F 4 -F 8 -q 0 $outbam $fimo_d/fimo.bed > $fimo_d/fimo.cuts.freq.txt
        $Rscriptbin/Rscript $ori_script_dir/run_centipede_parker.R $fimo_d/fimo.cuts.freq.txt $fimo_d/fimo.bed $fimo_d/fimo.png $mlen
        done

        >&2 echo "[info] Finished"
    done

    # cp /gcdata/gcproj/fulong/Data/chip_ct_comparison/0505-menin-new_data/processeddata/CUTTAG/pipeline/CUTTAG_MENIN_100k_DMSO_3d_S29/macs2.narrow.aug18/read.meme.py /gcdata/gcproj/fulong/Data/scCUT_RUN/0510_scCRT_module/processed/ctcf_v1
    # cp /gcdata/gcproj/fulong/Data/chip_ct_comparison/0505-menin-new_data/processeddata/CUTTAG/pipeline/CUTTAG_MENIN_100k_DMSO_3d_S29/macs2.narrow.aug18/fix_sequence.py /gcdata/gcproj/fulong/Data/scCUT_RUN/0510_scCRT_module/processed/ctcf_v1
    # cp /gcdata/gcproj/fulong/Data/chip_ct_comparison/0505-menin-new_data/processeddata/CUTTAG/pipeline/CUTTAG_MENIN_100k_DMSO_3d_S29/macs2.narrow.aug18/run_centipede_parker.R /gcdata/gcproj/fulong/Data/scCUT_RUN/0510_scCRT_module/processed/ctcf_v1

fi

# echo -e "`date` [info] Footprinting finished \n"

