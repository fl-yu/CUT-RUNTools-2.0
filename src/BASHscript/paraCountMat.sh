#! /bin/bash
countMat(){
    featurefile=$1
    outdir=$2
    outfile=$3
    bamdir=$4
    cores=$5
    bedtoolsbin=$6
    path_parallel=$7
    # bamfile=$8
    bamfile="${@:8}"
    mkdir -p $outdir/temp
    cut -f 1-3 $featurefile > $outdir/temp/temp.feature
    cd $bamdir
    $path_parallel/parallel -j $cores "$bedtoolsbin/bedtools coverage -a $outdir/temp/temp.feature -b {} | cut -f 4 > $outdir/temp/{}.temp" ::: $bamfile

    cd $outdir/temp
    ls -1 $outdir/temp/*.temp | split -l 150 -d - lists
    for i in lists*; do paste $(cat $i) > merge${i##lists}; done

    awk -v OFS="." '$1=$1' $outdir/temp/temp.feature > $outdir/temp/row.name.temp
    echo featureName $bamfile | awk -v OFS="\t" '$1=$1' > $outdir/temp/col.name.temp
    paste $outdir/temp/row.name.temp $outdir/temp/merge* > $outdir/temp/final.res.temp
    # echo `head -3 $outdir/temp/final.res.temp | tail -1 | tr '\t' '\n' | wc -l`
    cd $outdir
    cat $outdir/temp/col.name.temp $outdir/temp/final.res.temp > $outfile
    rm -rf $outdir/temp
}

