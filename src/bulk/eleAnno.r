



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")
library(ChIPpeakAnno)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


assignChromosomeRegion
