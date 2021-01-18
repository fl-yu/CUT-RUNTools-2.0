PeakFun
PeakFun -h 


peak.bed=$1
genome=$2
vis_num=$3

Help()
{
   # Display Help
   echo "It is a script for Gene Oncology analysis of input peaks."
   echo
   echo "Syntax: ./PeakFun.sh peak.bed genome vis_num outputdir"
   echo "options:"
   echo "peak.bed     Peak file with BED format, e.g. cluster2.narrowpeak"
   echo "genome       Genome version used, e.g. hg38"
   echo "vis_num      How many significant entries showed in the resulting barplot, e.g. 20"
   echo "outputdir    The directory for output, e.g. /dir/for/output"
   echo
}



get_peak_fun() {
  echo "John"
}

# diffpeak_file <- "/Users/fulongyu/Data/temp/group1_intergenic.peak1-3"
diffpeak_file1 <- "/Users/fulongyu/Data/temp/group2_promoter.peak"
diffpeak_file2 <- "/Users/fulongyu/Data/temp/group2_intergenic.peak"
diffpeak_file3 <- "/Users/fulongyu/Data/temp/group1_promoter.peak"
diffpeak_file4 <- "/Users/fulongyu/Data/temp/group1_intergenic.peak"
out_dir <- "/Users/fulongyu/Data/temp"
out_name1 <- "H1_promoter_peak_function"
out_name2 <- "H1_intergenic_peak_function"
out_name3 <- "K562_promoter_peak_function"
out_name4 <- "K562_intergenic_peak_function"
the_genome <- "hg38"
vis_num <- 20
suppressPackageStartupMessages({
  library(rGREAT)
  library(ggplot2)
  })
diffpeak_1 <- read.csv(diffpeak_file1, header=F, sep="\t")
diffpeak_2 <- read.csv(diffpeak_file2, header=F, sep="\t")
diffpeak_3 <- read.csv(diffpeak_file3, header=F, sep="\t")
diffpeak_4 <- read.csv(diffpeak_file4, header=F, sep="\t")
# rank peak with fc and get the first 1000
diffpeak_1_1k = diffpeak[order(diffpeak_1[, 9], decreasing = T), ]
if(nrow(diffpeak_1_1k)>1000){
    diffpeak_1_1k <- diffpeak_1_1k[1:1000, ]
}
diffpeak_2_1k = diffpeak[order(diffpeak_2[, 9], decreasing = T), ]
if(nrow(diffpeak_2_1k)>1000){
    diffpeak_2_1k <- diffpeak_2_1k[1:1000, ]
}
diffpeak_3_1k = diffpeak[order(diffpeak_3[, 9], decreasing = T), ]
if(nrow(diffpeak_3_1k)>1000){
    diffpeak_3_1k <- diffpeak_3_1k[1:1000, ]
}
diffpeak_4_1k = diffpeak[order(diffpeak_4[, 9], decreasing = T), ]
if(nrow(diffpeak_4_1k)>1000){
    diffpeak_4_1k <- diffpeak_4_1k[1:1000, ]
}
job1 = submitGreatJob(diffpeak_1_1k, species=the_genome, request_interval = 10)
job2 = submitGreatJob(diffpeak_2_1k, species=the_genome, request_interval = 10)
job3 = submitGreatJob(diffpeak_3_1k, species=the_genome, request_interval = 10)
job4 = submitGreatJob(diffpeak_4_1k, species=the_genome, request_interval = 10)
tb1 = getEnrichmentTables(job1)
tb2 = getEnrichmentTables(job2)
tb3 = getEnrichmentTables(job3)
tb4 = getEnrichmentTables(job4)
#names(tb) 
tb_bp1 <- tb1$`GO Biological Process`
tb_bp2 <- tb2$`GO Biological Process`
tb_bp3 <- tb3$`GO Biological Process`
tb_bp4 <- tb4$`GO Biological Process`
# if(sum(tb_bp$Binom_Adjp_BH < 0.05)==0){
#     stop("No functions were significantly enriched!")
# } else {
#     tb_bp_sig <- tb_bp[tb_bp$Binom_Adjp_BH < 0.05, ]
#     tb_bp_sig <- tb_bp_sig[order(tb_bp_sig$Binom_Adjp_BH), ]
# }

tb_bp_sig1 <- tb_bp1[tb_bp1$Hyper_Adjp_BH < 0.05, ]
tb_bp_sig1 <- tb_bp_sig1[order(tb_bp_sig1$Hyper_Adjp_BH), ]

tb_bp_sig2 <- tb_bp2[tb_bp2$Hyper_Adjp_BH < 0.05, ]
tb_bp_sig2 <- tb_bp_sig2[order(tb_bp_sig2$Hyper_Adjp_BH), ]

tb_bp_sig3 <- tb_bp3[tb_bp3$Hyper_Adjp_BH < 0.05, ]
tb_bp_sig3 <- tb_bp_sig3[order(tb_bp_sig3$Hyper_Adjp_BH), ]

tb_bp_sig4 <- tb_bp4[tb_bp4$Hyper_Adjp_BH < 0.05, ]
tb_bp_sig4 <- tb_bp_sig4[order(tb_bp_sig4$Hyper_Adjp_BH), ]

message(nrow(tb_bp_sig1), " items of functions were enriched")
message(nrow(tb_bp_sig2), " items of functions were enriched")
message(nrow(tb_bp_sig3), " items of functions were enriched")
message(nrow(tb_bp_sig4), " items of functions were enriched")
vis_num=30
tb2_bp_top1 = tb_bp_sig1[1:vis_num, ]
tb2_bp_top2 = tb_bp_sig2[1:vis_num, ]
tb2_bp_top3 = tb_bp_sig3[1:vis_num, ]
tb2_bp_top4 = tb_bp_sig4[1:vis_num, ]

# tb2_bp_top1 <- tb2_bp_top1[order(tb2_bp_top1$Hyper_Adjp_BH, decreasing=T),]
# tb2_bp_top1$name <- factor(tb2_bp_top1$name, levels = tb2_bp_top1$name)
# tb2_bp_top2 <- tb2_bp_top2[order(tb2_bp_top2$Hyper_Adjp_BH, decreasing=T),]
# tb2_bp_top2$name <- factor(tb2_bp_top2$name, levels = tb2_bp_top2$name)
# tb2_bp_top3 <- tb2_bp_top3[order(tb2_bp_top3$Hyper_Adjp_BH, decreasing=T),]
# tb2_bp_top3$name <- factor(tb2_bp_top3$name, levels = tb2_bp_top3$name)
# tb2_bp_top4 <- tb2_bp_top4[order(tb2_bp_top4$Hyper_Adjp_BH, decreasing=T),]
# tb2_bp_top4$name <- factor(tb2_bp_top4$name, levels = tb2_bp_top4$name)
tb2_bp_top1 <- tb2_bp_top1[order(tb2_bp_top1$Hyper_Adjp_BH, decreasing=F),]
tb2_bp_top1$name <- factor(tb2_bp_top1$name, levels = tb2_bp_top1$name)
tb2_bp_top2 <- tb2_bp_top2[order(tb2_bp_top2$Hyper_Adjp_BH, decreasing=F),]
tb2_bp_top2$name <- factor(tb2_bp_top2$name, levels = tb2_bp_top2$name)
tb2_bp_top3 <- tb2_bp_top3[order(tb2_bp_top3$Hyper_Adjp_BH, decreasing=F),]
tb2_bp_top3$name <- factor(tb2_bp_top3$name, levels = tb2_bp_top3$name)
tb2_bp_top4 <- tb2_bp_top4[order(tb2_bp_top4$Hyper_Adjp_BH, decreasing=F),]
tb2_bp_top4$name <- factor(tb2_bp_top4$name, levels = tb2_bp_top4$name)
setwd("/Users/fulongyu/Data/temp/motif2")
write.table(tb2_bp_top1, paste0(out_name1, "top30.csv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(tb2_bp_top2, paste0(out_name2, "top30.csv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(tb2_bp_top3, paste0(out_name3, "top30.csv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(tb2_bp_top4, paste0(out_name4, "top30.csv"), row.names=F, col.names=F, sep="\t", quote=F)



barplot_bp = ggplot(data=tb2_bp_top, aes(x=name, y=(-log10(Hyper_Adjp_BH)))) + 
             geom_bar(stat="identity", fill="steelblue")+
             theme_bw() +
             theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
             coord_flip() + labs(y= "-log10(BH adjust p values)", x = "GO Biological Process")





