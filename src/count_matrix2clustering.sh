# ============================================================================================
#    step3 count_matrix->clustering analysis (dimensionality reduction and clustering analysis)
# ============================================================================================
# define all the parameters in this step, useful to change when use this script alone
Rscriptbin=$Rscriptbin
sc_countMatrix_dir=$workdir/sc_countMatrix
sc_cluster_dir=$workdir/sc_cluster
r_function_dir=$scriptdir/Rscript
matrix_type=$matrix_type
count_matrix=$sc_countMatrix_dir/$matrix_type"_matrix.txt"
experiment_name=$experiment_name
pythonbin=$pythonbin

# clustering parameters
cluster_pc=$cluster_pc
cluster_resolution=$cluster_resolution


if [ "$entire_pipeline" == "TRUE" ]
then
    Rscriptbin=$Rscriptbin
    sc_countMatrix_dir=$workdir/sc_countMatrix
    sc_cluster_dir=$workdir/sc_cluster
    r_function_dir=$scriptdir/Rscript
    matrix_type=$matrix_type
    count_matrix=$sc_countMatrix_dir/$matrix_type"_matrix.txt"
    experiment_name=$experiment_name
    pythonbin=$pythonbin

    # clustering parameters
    cluster_pc=$cluster_pc
    cluster_resolution=$cluster_resolution
    echo "# "
    echo "#         step3 count_matrix->clustering analysis (dimensionality reduction and clustering analysis)"
    echo "# "
    date
elif [ "$entire_pipeline" == "FALSE" ] && [ "$individual_step" == "count_matrix2clustering" ]
then
    sc_cluster_dir=$step3_output_dir
    count_matrix=$step3_count_matrix

    cluster_pc=$cluster_pc
    cluster_resolution=$cluster_resolution
    Rscriptbin=$Rscriptbin
    pythonbin=$pythonbin
    experiment_name=$experiment_name
    matrix_type=$matrix_type
    r_function_dir=$scriptdir/Rscript
    echo "# "
    echo "#         step3 count_matrix->clustering analysis (dimensionality reduction and clustering analysis)"
    echo "# "
    echo "[info] Step3 will run separately... "

    date
else 
    echo "[ERROR] Please check your configuration!"
    exit 1
fi



mkdir -p $sc_cluster_dir
$Rscriptbin/Rscript $r_function_dir/scmatrix_analysis_pipeline.r \
                    $count_matrix \
                    $cluster_pc \
                    $cluster_resolution \
                    $sc_cluster_dir \
                    $experiment_name \
                    $r_function_dir \
                    $pythonbin

echo "# "
echo "#         step 3 (dimensionality reduction and clustering analysis) is completed"
echo "# "
echo "-------------------------------------------------------------------------------------------------------------------------"