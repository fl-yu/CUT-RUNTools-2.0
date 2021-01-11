## Example for bulk data analysis
### Easy 1-step execution

To walk through the pipeline quickly, we make an small example of GATA1 CUT&RUN dataset from [our paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4) only containing data from chromosome 11 (can be found in the folder `exampleData` ) to illustrate the usage of our bulk data analysis pipeline. It will take roughly less than half an hour to run the entire pipeline from raw sequencing data using a laptop with 4 cores.

Download the example data `exampleData-bulk.zip` with two fastq files suffix with R1_001.fastq.gz, and R2_001.fastq.gz in the `exampleData` folder.

To configure the exampleData_configure.sh file:   
> specify the *workdir*;  
  specify the *fastqdir* parameter as the folder of exampleData;  
  make sure you have successfully installed pre-requsite software, set the paths and validate the bulk-config.json file ;  

Then you can simply run the following command to perform the entire pipeline.

```
# Usage:
#   ./run_bulkModule.sh configure.json sample_name
chmod +x ./run_bulkModule.sh   
./run_bulkModule.sh /path/to/bulk-config.json GATA1_D7_30min_chr11
```

What happens next is that CUT&RUNTools 2.0 will sequentially run the analysis pipeline, see the [Usages](./bulk-USAGE.md) and output [Directory Structure](./bulk-DIRECTORY.md) for details. 



