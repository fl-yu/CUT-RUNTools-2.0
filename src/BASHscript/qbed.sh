#! /bin/bash
input_dir=$1
output_name=$2
tabix_bin=$3
output_dir=$4

cd $input_dir
bed_file=(`echo *.bed`)
bed_name=(`echo *.bed| sed s/.bed//g`)
length=${#bed_file[@]}

echo "[info] $length files will be processed into a single-cell genome track file"
for ((i=0; i<${length}; i++)); do 
    temp_name=${bed_name[$i]}
    temp_qbed_name="${temp_name}".qbed
    awk -v OFS="\t" -v temp_name="${bed_name[$i]}" -v list="$i" '{ print $0, list, "+", temp_name }' "${bed_file[$i]}" > "${temp_name}".temp.temp
done
echo "[info] Sort the qbed files"
cat *.temp.temp | sort -k1V -k2n -k3n - > $output_name.qbed
echo "[info] Generate the qbed file and compress it"

$tabix_bin/bgzip -f < $output_name.qbed > $output_dir/$output_name.qbed.gz
$tabix_bin/tabix -f -p bed $output_dir/$output_name.qbed.gz
echo "[info] Remove the temp files"
rm -rf *.temp.temp

# configure the json files https://eg.readthedocs.io/en/latest/datahub.html#example-qbed-track
# Usage:
#   sc_qbed_track /path/to/beds qbed_name /path/to/tabix /path/to/output
