

# input
diffpeak_file <- commandArgs(trailingOnly = T)[1] # "/Users/fulongyu/Data/temp/group0_group1.diffPeak.final"
out_dir <- commandArgs(trailingOnly = T)[2] # "/Users/fulongyu/Data/temp"
out_name <- commandArgs(trailingOnly = T)[3] # "group0_group1"
the_genome <- commandArgs(trailingOnly = T)[4] # "hg38"
vis_num <- commandArgs(trailingOnly = T)[5] # "20"

suppressPackageStartupMessages({
  library(rGREAT)
  library(ggplot2)
  })

diffpeak <- read.csv(diffpeak_file, header=F, sep="\t")
# rank peak with fc and get the first 1000
diffpeak_1k = diffpeak[order(diffpeak[, 9], decreasing = T), ]
if(nrow(diffpeak_1k)>1000){
    diffpeak_1k <- diffpeak_1k[1:1000, ]
}
job = submitGreatJob(diffpeak_1k, species=the_genome)
tb = getEnrichmentTables(job)
#names(tb) 
tb_bp <- tb$`GO Biological Process`
# if(sum(tb_bp$Binom_Adjp_BH < 0.05)==0){
#     stop("No functions were significantly enriched!")
# } else {
#     tb_bp_sig <- tb_bp[tb_bp$Binom_Adjp_BH < 0.05, ]
#     tb_bp_sig <- tb_bp_sig[order(tb_bp_sig$Binom_Adjp_BH), ]
# }
if(sum(tb_bp$Hyper_Adjp_BH < 0.05)==0){
    stop("No functions were significantly enriched!")
} else {
    tb_bp_sig <- tb_bp[tb_bp$Hyper_Adjp_BH < 0.05, ]
    tb_bp_sig <- tb_bp_sig[order(tb_bp_sig$Hyper_Adjp_BH), ]
}
message(nrow(tb_bp_sig), "items of functions were enriched")
if(nrow(tb_bp_sig)>vis_num){
    message("Top", vis_num, "items will be showed in the figure")
    tb2_bp_top = tb_bp_sig[1:vis_num, ]
} else {
    message("All of them will be showed in the figure")
    tb2_bp_top = tb_bp_sig
}
tb2_bp_top$name <- factor(tb2_bp_top$name, levels = tb2_bp_top$name)
barplot_bp = ggplot(data=tb2_bp_top, aes(x=name, y=(-log10(Hyper_Adjp_BH)))) + 
             geom_bar(stat="identity", fill="steelblue")+
             theme_bw() +
             theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
             coord_flip() + labs(y= "-log10(BH adjust p values)", x = "GO Biological Process")
png(paste0(out_dir,"/", out_name, "_enriched_functions.png"))
barplot_bp
invisible(dev.off())
message("[ok]enriched function figures")

