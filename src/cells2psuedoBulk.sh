# ========================================================================================
# step4 cells->pseudoBulk (cell-type-specific pseudo-bulk data analysis)
# ========================================================================================
# define parameters
if [ "$entire_pipeline" == "TRUE" ]
then
    Rscriptbin=$Rscriptbin
    r_function_dir=$scriptdir/Rscript
    bash_function_dir=$scriptdir/BASHscript
    bam_dir=$workdir/sc_aligned.aug10/dup.marked.clean
    cell_anno=$workdir/sc_cluster/leiden_cluster_annotation.txt
    sc_pseudoBulk_dir=$workdir/sc_pseudoBulk
    samtoolsbin=$samtoolsbin
    macs2bin=$macs2bin
    path_deeptools=$path_deeptools
    path_tabix=$path_tabix
    genome=$genome
    cores=$cores
    extratoolsbin=$extratoolsbin
    peak_type=$peak_type
    peak_caller=$peak_caller
    echo "# "
    echo "#         step4 cells2psuedoBulk (cell-type-specific pseudo-bulk data analysis)"
    echo "# "
    date
elif [ "$entire_pipeline" == "FALSE" ] && [ "$individual_step" == "cells2psuedoBulk" ]
then
    bam_dir=$step4_bamfile_dir
    cell_anno=$step4_cell_anno_file
    sc_pseudoBulk_dir=$step4_output_dir

    Rscriptbin=$Rscriptbin
    r_function_dir=$scriptdir/Rscript
    bash_function_dir=$scriptdir/BASHscript
    samtoolsbin=$samtoolsbin
    macs2bin=$macs2bin
    path_deeptools=$path_deeptools
    path_tabix=$path_tabix
    genome=$genome
    cores=$cores
    extratoolsbin=$extratoolsbin
    peak_type=$peak_type
    peak_caller=$peak_caller
    echo "# "
    echo "#         step4 cells2psuedoBulk (cell-type-specific pseudo-bulk data analysis)"
    echo "# "
    echo "[info] Step4 will run separately... "
    date
else 
    echo "[ERROR] Please check your configuration!"
    exit 1
fi



mkdir -p $sc_pseudoBulk_dir 

$Rscriptbin/Rscript $r_function_dir/pseudobulk_analysis_sub.r \
                    $cell_anno \
                    $bam_dir \
                    $sc_pseudoBulk_dir \
                    $samtoolsbin \
                    $macs2bin \
                    $path_deeptools \
                    $genome \
                    $path_tabix \
                    $bash_function_dir \
                    $cores \
                    $extratoolsbin \
                    $Rscriptbin \
                    $peak_type \
                    $peak_caller

echo "# "
echo "#         step 4 (cell-type-specific pseudo-bulk data analysis) is completed"
echo "# "
echo "-------------------------------------------------------------------------------------------------------------------------"