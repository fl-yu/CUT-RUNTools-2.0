# input
qcfile <- commandArgs(trailingOnly = T)[1]
outputdir <- commandArgs(trailingOnly = T)[2]
insert_d <- commandArgs(trailingOnly = T)[3]
reads_threshold <- as.numeric(commandArgs(trailingOnly = T)[4])
percentage_threshold <- as.numeric(commandArgs(trailingOnly = T)[5])
suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(gridExtra)
  })
#### 
#### draw figures for QC metrices
#### 
a <- read.csv(qcfile, header=T, sep="\t")
#brewer.pal(4, "Dark2")
#### fig 1
p1 <- ggplot(a, aes(x="", y=Overall_alignment_ratio)) + 
    geom_violin(width=0.5, lwd=1.2, fill="#1B9E77") + 
    geom_boxplot(width=0.02, outlier.shape = NA) +
    theme_classic(base_size=12)+
    xlab("Reads overall alignment") + ylab("Percentage")
# Reads_properly_paired_percentage
p2 <- ggplot(a, aes(x="", y=Reads_properly_paired_percentage)) + 
    geom_violin(width=0.5, lwd=1.2, fill="#D95F02") + 
    geom_boxplot(width=0.02, outlier.shape = NA) +
    theme_classic(base_size=13) +
    xlab("Reads properly paired") + ylab("Percentage")
# Reads_duplicated_percentage
p3 <- ggplot(a, aes(x="", y=Reads_duplicated_percentage)) + 
    geom_violin(width=0.5, lwd=1.2, fill="#7570B3") + 
    geom_boxplot(width=0.02, outlier.shape = NA) +
    theme_classic(base_size=13) +
    xlab("Reads duplicated") + ylab("Percentage")
# Reads_nuclear_percentage
p4 <- ggplot(a, aes(x="", y=Reads_nuclear_percentage)) + 
    geom_violin(width=0.5, lwd=1.2, fill="#E7298A") + 
    geom_boxplot(width=0.02, outlier.shape = NA) +
    theme_classic(base_size=13) +
    xlab("Reads nuclear") + ylab("Percentage")
png(filename = paste0(outputdir, "/QC-statistics-1.png"), width = 1000, height = 500)
grid.arrange(p1, p2, p3, p4, ncol=4)
invisible(dev.off())
# pdf
pdf(paste0(outputdir, "/QC-statistics-1.pdf"))
grid.arrange(p1, p2, p3, p4, ncol=4)
invisible(dev.off())
message("[ok]QC-figure1")

#### fig 2
plot2 = ggplot(a, aes(x = log10(Reads_properly_paired), y = Reads_in_peak_percentage)) +
  geom_hex(bins = 50) +
  theme_bw() + scale_fill_viridis() +
  xlab("log10 Properly Paired Reads") +
  ylab("Percentage of Reads Enrichment in Peaks") +
  geom_hline(yintercept = percentage_threshold, lty = "dashed") +
  geom_vline(xintercept = log10(reads_threshold), lty = "dashed") +
  theme_classic() +
  theme(axis.title.x = element_text(size=10))
# png
png(filename = paste0(outputdir, "/QC-statistics-2.png"), width = 800, height = 500)
plot2
invisible(dev.off())
# pdf
pdf(paste0(outputdir, "/QC-statistics-2.pdf"))
plot2
invisible(dev.off())
message("[ok]QC-figure2")

#### fig 3
xx <- read.csv(insert_d, header=F, sep="\t")
# png
png(filename = paste0(outputdir, "/QC-statistics-3.png"), width = 550, height = 500)
par(mar=c(5, 4.3, 4, 2) + 0.5)
plot(xx[,1], (xx[,2])/1000000, type="l", col="indianred",lwd=4, xlab="Insert Size", ylab=expression("Fragment Counts" ~ 10^6))
invisible(dev.off())
# pdf
pdf(paste0(outputdir, "/QC-statistics-3.pdf"))
par(mar=c(5, 4.3, 4, 2) + 0.5)
plot(xx[,1], (xx[,2])/1000000, type="l", col="indianred",lwd=4, xlab="Insert Size", ylab=expression("Fragment Counts" ~ 10^6))
invisible(dev.off())
message("[ok]QC-figure3")
