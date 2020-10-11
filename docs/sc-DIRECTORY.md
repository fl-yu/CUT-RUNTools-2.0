## The directory structure of input and output folders 

*A schematic view of the data structure for the typical single-cell data analysis*


```
+ fastq_directory (INPUT)
+ workdir         (OUTPUT)
    ++ sc_trimmed
    ++ sc_trimmed3
    ++ sc_log
	+++ trim
	+++ bowtie2
	+++ bamSort
	+++ bamMarkdup
    ++ sc_aligned.aug10
	+++ sorted
	+++ dup.marked
	+++ dup.marked.clean
    ++ sc_countMatrix
    ++ sc_cluster
    ++ sc_qc
	+++ report
	+++ figure
    ++ sc_pseudoBulk
	+++ groups_aggregation
	    ++++ single_cells
	    ++++ pseudo_bulk_data
		++++++  macs2
		++++++  SEACR
	++++ group_1
	    ++++ single_cells
	    ++++ pseudo_bulk_data
	++++ group_2
		++++++ single_cells
		++++++ pseudo_bulk_data
```

