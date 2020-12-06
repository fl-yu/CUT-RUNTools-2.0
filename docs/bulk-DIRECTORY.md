## The directory structure of input and output folders 

*A typical input and output data structure for the typical bulk data analysis*


```
+ fastq_directory (INPUT)
+ workdir         (OUTPUT)
    ++ trimmed
    ++ trimmed2
    ++ spike_in
    ++ aligned
        +++ dedup
        +++ dedup.120bp
        +++ dup.marked
        +++ dup.marked.120bp
        +++ sorted
    ++ peakcalling
        +++ macs2.narrow
            ++++ blk_filtered
            ++++ blk_filtered.fa
            ++++ centipede.bam
            ++++ random.1000
            ++++ fimo.result
                +++++ fimo2.DREME-1.CTCHGCC
                +++++ fimo2.DREME-2.RAAATA
                +++++ fimo2.DREME-3.GGCTCACD
        +++ macs2.broad
        +++ seacr
        

See the following links for more details:

- [Overview](./bulk-OVERVIEW.md)
- [Installation Instructions](./bulk-INSTALL.md)
- [Quick Start](./bulk-QUICK.md)
- [Usage Page](./bulk-USAGE.md)

