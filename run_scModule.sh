#!/usr/bin/bash
# Copyright (C) 2020 Fulong Yu
#
# CUT&RUNTools 2.0 is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools 2.0 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.

configrue_file=$1 # $scriptdir/sc-config.json
SCRIPT=`readlink -f $0`
SCRIPTPATH=`dirname $SCRIPT`
scriptdir=$SCRIPTPATH/src


# convert configuration JSON files to bash variables
eval "$(jq -r '.software_config | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
eval "$(jq -r '.input_output | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
eval "$(jq -r '.run_pipeline | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"


workdir=$workdir/$experiment_name

if [ "$entire_pipeline" == "TRUE" ]
then
    logdir=$workdir/sc_logs
    qcdir=$workdir/sc_qc
    sc_countMatrix_dir=$workdir/sc_countMatrix
    bamdir=$workdir/sc_aligned.aug10/dup.marked.clean
    sc_pseudoBulk_dir=$workdir/sc_pseudoBulk
    matrix_type=$matrix_type
    bin_size=$bin_size
    featureFile=$featureFile
    sc_cluster_dir=$workdir/sc_cluster
    echo "==================================== The entire pipeline will run! ======================================================================"
    >&2 echo -e "## Experiment name:               $experiment_name"
    sleep 0.1
    >&2 echo -e "## Input FASTQ folder:            $fastq_directory"
    sleep 0.1
    >&2 echo -e "## Workdir folder:                $workdir"
    sleep 0.1
    >&2 echo -e "## Log files floder:              $logdir"
    sleep 0.1
    >&2 echo -e "## QC reports and figures folder: $qcdir"
    sleep 0.1
    >&2 echo -e "## Experiment type:               $experiment_type"
    sleep 0.1
    >&2 echo -e "## Reference genome:              $genome"
    sleep 0.1
    >&2 echo -e "## Processors used:               $cores"
    sleep 0.1
    >&2 echo -e "## BAM folder of cells:           $bamdir"
    sleep 0.1
    >&2 echo -e "## Peak_caller:                   $peak_caller"
    sleep 0.1
    >&2 echo -e "## Cell filter [read in peak %]:    $percentage_rip"
    sleep 0.1
    >&2 echo -e "## Cell filter [qulified reads]:  $num_reads_threshold"
    sleep 0.1
    >&2 echo -e "## Count matrix folder:           $sc_countMatrix_dir"
    sleep 0.1
    >&2 echo -e "## Count matrix type:             $matrix_type"
    sleep 0.1
    >&2 echo -e "## Cell clustering folder:        $sc_cluster_dir"
    sleep 0.1
    >&2 echo -e "## Cell clustering [resolution]:  $cluster_resolution"
    sleep 0.1
    >&2 echo -e "## Cell clustering [PCs]:         $cluster_pc"
    sleep 0.1
    >&2 echo -e "## Pseudo-bulk data folder:       $sc_pseudoBulk_dir"
    echo -e "========================================================================================================================================="
    # step1 fastq -> peak
    . $scriptdir/fastq2peak.sh
    # step2 cells -> count_matrix
    . $scriptdir/cells2count_matrix.sh
    # step3 count_matrix -> clustering
    . $scriptdir/count_matrix2clustering.sh
    # step4 cells -> psuedoBulk
    . $scriptdir/cells2psuedoBulk.sh
fi

if [ "$entire_pipeline" != "TRUE" ]
then
        echo "==================================== Individual step will run! =========================================================="
    # run step1
    if [ "$individual_step" == "fastq2peak" ]
    then
        echo "==================================== step 1 (raw data processing) will be performed! ===================================="
        >&2 echo -e "## Experiment name:               $experiment_name"
        sleep 0.1
        >&2 echo -e "## Input FASTQ folder:            $fastq_directory"
        sleep 0.1
        >&2 echo -e "## Workdir folder:                $workdir"
        sleep 0.1
        >&2 echo -e "## Log files floder:              $logdir"
        sleep 0.1
        >&2 echo -e "## QC reports and figures folder: $qcdir"
        sleep 0.1
        >&2 echo -e "## Experiment type:               $experiment_type"
        sleep 0.1
        >&2 echo -e "## Reference genome:              $genome"
        sleep 0.1
        >&2 echo -e "## Processors used:               $cores"
        sleep 0.1
        >&2 echo -e "## BAM folder of cells:           $bamdir"
        sleep 0.1
        >&2 echo -e "## Peak_caller:                   $peak_caller"
        sleep 0.1
        >&2 echo -e "## Cell filter [read in peak]:    $percentage_rip"
        sleep 0.1
        >&2 echo -e "## Cell filter [qulified reads]:  $num_reads_threshold"
        echo "========================================================================================================================="
        . $scriptdir/fastq2peak.sh
    # run step2
    elif [ "$individual_step" == "cells2count_matrix" ]
    then

        echo "==================================== step 2 (count matrix generation) will be performed! ================================"
        >&2 echo -e "## Input individual BAM folder:  $step2_bamfile_dir"
        sleep 0.1
        >&2 echo -e "## Output count matrix folder:   $step2_output_dir"
        sleep 0.1
        >&2 echo -e "## Cells passed QC file:         $step2_qc_pass_file"
        sleep 0.1
        >&2 echo -e "## Count matrix type:            $matrix_type"
        sleep 0.1
        echo "========================================================================================================================="
        . $scriptdir/cells2count_matrix.sh
    # run step3
    elif [ "$individual_step" == "count_matrix2clustering" ]
    then
        echo "==================================== step 3 (clustering analysis) will be performed! ===================================="
        >&2 echo -e "## Input count matrix file:      $step3_count_matrix"
        sleep 0.1
        >&2 echo -e "## Output clustering folder:     $step3_output_dir"
        sleep 0.1
        >&2 echo -e "## Cell clustering [resolution]: $cluster_resolution"
        sleep 0.1
        >&2 echo -e "## Cell clustering [PCs]:        $cluster_pc"
        echo "========================================================================================================================="
        . $scriptdir/count_matrix2clustering.sh
    # run step4
    elif [ "$individual_step" == "cells2psuedoBulk" ]
    then
        echo "==================================== step 4 (cell-type-specific aggregation) will be performed! ========================="
        >&2 echo -e "## Input individual BAM folder:                   $step4_bamfile_dir"
        sleep 0.1
        >&2 echo -e "## Cell label file:                               $step4_cell_anno_file"
        sleep 0.1
        >&2 echo -e "## Output cell type-specific pseudo-bulk folder:  $step4_output_dir"
        echo "========================================================================================================================="
        . $scriptdir/cells2psuedoBulk.sh
    # error
    else
        echo "[ERROR] Please specify the individual_step as one of fastq2peak; cells2count_matrix; count_matrix2clustering; cells2psuedoBulk"
        exit 1
    fi
fi
