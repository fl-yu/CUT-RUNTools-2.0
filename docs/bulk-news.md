
#### New features for bulk data analysis

In this release, we provided several new features in the following which are urgently needed from the feedback of the community. Also, we slightly revised the design of our code to be more flexible and handy for data processing and analysis.
- *supporting spike-in sequence alignment and data normalization*   
    Users can provide spike-in genome data (usually E.coli or fly genome) to align their fastq data. The statistic of aligned spike-in reads can be used for normalization of CUT&RUN or CUT&Tag 
- *options of the experiment for CUT&RUN or CUT&Tag data analysis*  
    Different enzymes are used in CUT&RUN (pAG-MNase) and CUT&Tag (Tn5) experiments. Although the overall peak calling strategies for these two types of data are basically the same, footprinting should be different. CUT&Tag reads should be shifted + 4 bp and − 5 bp for positive and negative strand respectively, to account for the 9-bp duplication created by DNA repair of the nick by Tn5 transposase and achieve base-pair resolution of TF footprint and motif-related analyses. CUT&RUN 2.0 now provides a new option for the user to select data types.
- *flexible option for fragments selections (>120bp)*  
    Reads filter step (size 120 bp) is good to find the enriched signal for TF data, while it is not suitable for histone modification data whose fragment size generally large than 150 bp. CUT&RUN 2.0 now provides a new option for the user to select whether to perform this fragment selection according to their data type.
- *Different peak calling strategies and new functions for peaks annotation*  
    As different types of data including different TFs and histone modifications can be broadly detected by CUT&RUN and CUT&Tag methods, three commonly used peak calling methods were provided in CUT&RUN 2.0. Several peak annotation functions were also provided.
- *compatible with more computational platforms*    
    To compatible with more computational platforms for bulk data analysis, the requirement of the SLURM job submission environment was removed. 

