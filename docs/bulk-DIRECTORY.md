## The directory structure of input and output folders 

*A schematic view of the data structure for the typical bulk data analysis*


```
+ fastq_directory (INPUT)
+ workdir         (OUTPUT)
    ++ trimmed
    ++ trimmed3
    ++ aligned.aug10
        +++ dedup
        +++ dedup.120bp
        +++ dup.marked
        +++ dup.marked.120bp
    ++ macs2.narrow.aug18
        +++ blk_filtered
        +++ blk_filtered.fa
        +++ centipede.bam
        +++ fimo.result
            ++++ fimo2.DREME-1.CTCHGCC
            ++++ fimo2.DREME-2.RAAATA
            ++++ fimo2.DREME-3.GGCTCACD
        +++ random.1000
    ++ macs2.narrow.aug18.dedup
    ++ macs2.narrow.all.frag.aug18
    ++ macs2.narrow.all.frag.aug18.dedup
    ++ macs2.broad.aug18
    ++ macs2.broad.aug18.dedup
    ++ macs2.broad.all.frag.aug18
    ++ macs2.broad.all.frag.aug18.dedup
    ++ seacr.aug12
    ++ seacr.aug12.all.frag
    ++ seacr.aug12.all.frag.dedup
    ++ seacr.aug12.dedup
```


See the following links for more details:

- [Overview](./bulk-OVERVIEW.md)
- [Installation Instructions](./bulk-INSTALL.md)
- [Quick Start](./bulk-QUICK.md)
- [Usage Page](./bulk-USAGE.md)
- [The full tutorial of CUT&RUNTools 2.0](./2.0-TUTORIAL.md)

