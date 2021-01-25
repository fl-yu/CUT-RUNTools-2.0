eleAnno







## import the MACS output
macs <- system.file("extdata", "MACS_peaks.xls", package="ChIPpeakAnno")
macsOutput <- toGRanges(macs, format="MACS")
## annotate the peaks with precompiled ensembl annotation
data(TSS.human.GRCh38)
macs.anno <- annotatePeakInBatch(macsOutput, AnnotationData=TSS.human.GRCh38)
## add gene symbols
library(org.Hs.eg.db)
macs.anno <- addGeneIDs(annotatedPeak=macs.anno, 
                        orgAnn="org.Hs.eg.db", 
                        IDs2Add="symbol")

if(interactive()){## annotate the peaks with UCSC annotation
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    ucsc.hg38.knownGene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    macs.anno <- annotatePeakInBatch(macsOutput, 
                                     AnnotationData=ucsc.hg38.knownGene)
    macs.anno <- addGeneIDs(annotatedPeak=macs.anno, 
                            orgAnn="org.Hs.eg.db", 
                            feature_id_type="entrez_id",
                            IDs2Add="symbol")
}




if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
    aCR<-assignChromosomeRegion(gr1, nucleotideLevel=FALSE, 
                           precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                           TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
    barplot(aCR$percentage)
}

