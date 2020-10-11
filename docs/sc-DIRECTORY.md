
## Directory structure of input and output folders 

```
	+ fastq_directory
	+ workdir
	    + sc_trimmed
	    + sc_trimmed3
	    + sc_log
	    	+ trim
	    	+ bowtie2
	    	+ bamSort
	    	+ bamMarkdup
	    + sc_aligned.aug10
	    	+ sorted
	    	+ dup.marked
	    	+ dup.marked.clean
	    + sc_countMatrix
	    + sc_cluster
	    + sc_qc
	    	+ report
	    	+ figure
	    + sc_pseudoBulk
	        + groups_aggregation
	            + single_cells
	            + pseudo_bulk_data
	                +  macs2
	                +  SEACR
	        + group_1
	            + single_cells
	            + pseudo_bulk_data
	        + group_2
	        	+ single_cells
	        	+ pseudo_bulk_data
  
```

