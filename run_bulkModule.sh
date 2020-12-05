#!/usr/bin/bash
# Copyright (C) 2020 Fulong Yu
#
# CUT&RUNTools 2.0 is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools 2.0 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.

configrue_file=$1 # $scriptdir/bulk-config.json
sample_name=$2 # $scriptdir/bulk-config.json
SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
scriptdir=$SCRIPTPATH/src/bulk
chrom_size_file=$SCRIPTPATH/assemblies/chrom."$organism_build"/"$organism_build".chrom.sizes
blacklist=$SCRIPTPATH/blacklist/"$organism_build".blacklist.bed
blacklist2=$SCRIPTPATH/blacklist/"$organism_build"_TA_repeat.bed
# convert configuration JSON files to bash variables
eval "$(jq -r '.software_config | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
eval "$(jq -r '.input_output | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"
eval "$(jq -r '.motif_finding | to_entries | .[] | .key + "=\"" + .value + "\""' < $configrue_file)"

workdir=$workdir/$experiment_name
# check the parameters
if [ "$spike_in" == "none" ]
then
    if [ "$spike_in_norm" == "TRUE" ]
    then
    >&2 echo [ERROR] Check the parameters of spike_in and spike_in_norm
    exit
    fi
fi
if [ "$organism_build" != "hg38" ] && [ "$organism_build" != "hg19" ] && [ "$organism_build" != "mm10" ] && [ "$organism_build" != "mm9" ]
    then
    echo "organism_build should be one of hg38, hg19, mm10 or mm9"
    exit 1
fi
if [ "$organism_build" == "hg38" ] || [ "$organism_build" == "hg19" ]
then
    macs2_genome=hs
else 
    macs2_genome=mm
fi


# also check T and F for spike_in_norm.....

echo "==================================== Bulk data analysis pipeline will run =============================================================="
    >&2 echo -e "## Input FASTQ folder:            $fastq_directory"
    sleep 0.1
    >&2 echo -e "## Sample name:                   $sample_name"
    sleep 0.1
    >&2 echo -e "## Workdir folder:                $workdir"
    sleep 0.1
    >&2 echo -e "## Experiment name:               $experiment_name"
    sleep 0.1
    >&2 echo -e "## Experiment type:               $experiment_type"
    sleep 0.1
    >&2 echo -e "## Reference genome:              $organism_build"
    sleep 0.1
    >&2 echo -e "## Spike-in genome:               $spike_in"
    sleep 0.1
    >&2 echo -e "## Spike-in normalization:        $spike_in_norm"
    sleep 0.1
    >&2 echo -e "## Fragment 120 filtration:       $frag_120"
    sleep 0.1
    echo -e "================================================================================================================================="

. $scriptdir/bulk-pipeline.sh







