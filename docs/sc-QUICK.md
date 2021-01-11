## Example for single-cell data analysis
### Easy 1-step execution

We used a published single-cell CUT&Tag dataset which detected the H3K27me3 signal of individual H1 and K562 cells ([Kaya-Okur, et al. Nat Commun 10, 1930 (2019)](https://doi.org/10.1038/s41467-019-09982-5)), as an example to illustrate the utility and functionality of our tool. To walk through the pipeline quickly, we make an small example dataset only containing data from chr11 of 200 cells (can be found in the folder `exampleData`), by randomly selecting 100 cells for the two cell types, respectively. It will take roughly half an hour to run the entire pipeline from raw sequencing data using a laptop with 8 cores.

Download and uncompress the exampleData.tar.gz

```
tar -zxvf exampleData-sc.tar.gz
```

To configure the exampleData_configure.sh file:   
> specify the *workdir*;  
  specify the *fastqdir* parameter as the folder of exampleData;  
  make sure you have successfully installed pre-requsite software, set the paths and validate the bulk-config.json file  

Then you can simply run the following command to perform the entire pipeline.

```
chmod +x ./run_scModule.sh   
./run_scModule.sh /path/to/sc-config.json
```

What happens next is that CUT&RUNTools will sequentially run the analysis pipeline, see the [Usages](./sc-USAGE.md) and output [Directory Structure](./sc-DIRECTORY.md) for details. 



