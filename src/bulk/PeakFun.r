
#!/n/app/R/3.3.3/bin/Rscript
the_genome=commandArgs(trailingOnly = T)[1]
peak_file=commandArgs(trailingOnly = T)[2]
outputdir=commandArgs(trailingOnly = T)[3]
vis_num=as.numeric(commandArgs(trailingOnly = T)[4])
# print(the_genome)
# print(peak_file)
# print(outputdir)
# print(vis_num)
# the_genome="hg38"
# peak_file="FLI1_fimo.bed"
# outputdir="/Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test"
# vis_num=30

if(the_genome != "hg19" & the_genome != "hg38" & the_genome != "mm10" & the_genome != "mm9"){
    stop("[STOP] the genome should be one of hg38, hg19, mm10 or mm9")
}

setwd(outputdir)
suppressPackageStartupMessages({
  library(rGREAT)
  library(ggplot2)
  })
peak_file <- read.csv(peak_file, header=F, sep="\t")
if(nrow(peak_file) >= 5000){
    message("[WARNING] The peak number is more than 5000, you may need to refine the peakset.")
}
job1 = submitGreatJob(peak_file, species=the_genome, request_interval = 10)
tb1 = getEnrichmentTables(job1)
#names(tb) 
tb_bp <- tb1$`GO Biological Process`
message("[INFO] Using Hyper_Adjp_BH 0.05 as the cutoff.")
if(sum(tb_bp$Hyper_Adjp_BH < 0.05)==0){
    message("[INFO] The results can be found in the ", outputdir)
    stop("[STOP] No functions were significantly enriched!")
} else {
    tb_bp_sig <- tb_bp[tb_bp$Hyper_Adjp_BH < 0.05, ]
    tb_bp_sig <- tb_bp_sig[order(tb_bp_sig$Hyper_Adjp_BH), ]
    message("[INFO]", nrow(tb_bp_sig), " items of functions were enriched")
}
write.csv(tb_bp_sig, "enriched_biological_function.csv", quote=F, row.names = FALSE)
if(nrow(tb_bp_sig) >= vis_num){
    tb_bp_top = tb_bp_sig[1:vis_num, ]
} else {
    tb_bp_top = tb_bp_sig
}
tb_bp_top$name <- factor(tb_bp_top$name, level=rev(tb_bp_top$name))
barplot_bp = ggplot(data=tb_bp_top, aes(x=name, y=(-log10(Hyper_Adjp_BH)))) + 
             geom_bar(stat="identity", fill="steelblue")+
             theme_bw() +
             theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
             coord_flip() + labs(y= "-log10(BH adjust p values)", x = "GO Biological Process") +
             geom_hline(yintercept=-log10(0.05), linetype = "dashed")
pdf(paste("Barplot of top ", vis_num, "enriched functions.pdf"))
print(barplot_bp)
dev.off()
message("[INFO] The results can be found in the ", outputdir)
