## Quickstart
### Easy 1-step execution

To walk through the pipeline quickly, we make an small example dataset only containing data from chr11 (can be found in the folder exampleData). It will take roughly less than half an hour to run the entire pipeline from raw sequencing data using a laptop with 8 cores.

Uncompress the example data of two fastq files suffix with R1_001.fastq.gz, and R2_001.fastq.gz in a directory.

```
unzip exampleData-bulk.zip
```

To configure the exampleData_configure.sh file:   
> specify the *workdir*;  
  specify the *fastqdir* parameter as the folder of exampleData;  
  set the paths of pre-requsite software;  

Then you can simply run the following command to perform the entire pipeline.

```
# Usage:
#   ./run_bulkModule.sh configure.json sample_name
chmod +x ./run_bulkModule.sh   
./run_bulkModule.sh /path/to/bulk-config.json GATA1_D7_30min_chr11
```

What happens next is that CUT&RUNTools 2.0 will sequentially run the analysis pipeline (for details see the Usage below). 


See the following links for more details:

- [New updates](./bulk-news.md)
- [Installation Page](./bulk-INSTALL.md)
- [Usages](./bulk-USAGE.md)
- [Directory Structure](./bulk-DIRECTORY.md)


