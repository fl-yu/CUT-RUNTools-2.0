## Quickstart
### Easy 1-step execution

We used a published single-cell CUT&Tag dataset which detected the H3K27me3 signal of individual H1 and K562 cells ([Kaya-Okur, et al. Nat Commun 10, 1930 (2019)](https://doi.org/10.1038/s41467-019-09982-5)), as an example to illustrate the utility and functionality of our tool. To walk through the pipeline quickly, we make an small example dataset only containing data from chr11 of 200 cells (can be found in the folder exampleData), by randomly selecting 100 cells for the two cell types, respectively. It will take roughly less than half an hour to run the entire pipeline from raw sequencing data using a laptop with 8 cores.

Uncompress the exampleData.tar.gz

```
tar -zxvf exampleData.tar.gz
```

To configure the exampleData_configure.sh file:   
> specify the *workdir*;  
  specify the *fastqdir* parameter as the folder of exampleData;  
  set the paths of pre-requsite software;  

Then you can simply run the following command to perform the entire pipeline.

```
chmod +x ./run_scModule.sh   
./run_scModule.sh /path/to/scCutRunTools/src /path/to/exampleData_configure.sh
```

What happens next is that CUT&RUNTools will sequentially run the analysis pipeline (for details see the Usage and Tutorial below). 


See the following links for more details:

- [Overview](./sc-OVERVIEW.md)
- [Installation Instructions](./sc-INSTALL.md)
- [Usage Page](./sc-USAGE.md)
- [Directory Structure](./sc-DIRECTORY.md)

