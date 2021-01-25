
#!/n/app/R/3.3.3/bin/Rscript
peak_overlap=commandArgs(trailingOnly = T)[1]
peak_count=commandArgs(trailingOnly = T)[2]
outputdir=commandArgs(trailingOnly = T)[3]

peak_overlap <- read.table(peak_overlap, header=F)[, 1]
peak_count <- read.table(peak_count, header=F)[, 1]

peak_per <- round(peak_overlap*100/peak_count, 2)
ele_name <- factor(c("UTR-5", "Promoter", "Exon", "Intron", "UTR-3", "Intragenic", "Intergenic"), levels=rev(c("UTR-5", "Promoter", "Exon", "Intron", "UTR-3", "Intragenic", "Intergenic")))
setwd(outputdir)

library(ggplot2)
df <- data.frame(peak_per, ele_name)

peak_overlap_per <- ggplot(data=df, aes(x=peak_per, y=ele_name)) +
                    geom_bar(stat="identity", fill="steelblue")+
                    theme_minimal() + ylab("genomic elements") + xlab("percentage")

pdf("Overlap of peaks and genomic elements.pdf")
print(peak_overlap_per)
dev.off()
message("[INFO] The results can be found in the ", outputdir)

