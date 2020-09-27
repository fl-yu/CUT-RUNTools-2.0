qcfile <- commandArgs(trailingOnly = T)[1]
outputdir <- commandArgs(trailingOnly = T)[2]
reads_threshold <- as.numeric(commandArgs(trailingOnly = T)[3])
percentage_threshold <- as.numeric(commandArgs(trailingOnly = T)[4])
#### 
#### generating the QC reports for experiment and barcode cells
#### 
a <- read.csv(qcfile, header=T, sep="\t")
report_keys = c("Overall Experimental Statistics", "Total Barcode Cells", "Total Reads", "Statistics of Barcode Cells (median)", "Alignment Ratio", "Reads Properly Paired (%)", "Reads Duplicated (%)", "Reads Quality > MAPQ30 (%)", "Reads Nuclear (%)", "Reads in Peaks (%)", "Insert Size Average")
reads_per_cell = round(sum(a$Total_reads)/nrow(a), 0)
report_values = c("", nrow(a), sum(a$Total_reads), "", paste0(round(median(a$Overall_alignment_ratio), 1), "%"), paste(round(median(a$Reads_properly_paired), 1), paste0("(", round(median(a$Reads_properly_paired_percentage), 1),"%)")), paste(round(median(a$Reads_duplicated), 1), paste0("(", round(median(a$Reads_duplicated_percentage), 1),"%)")), paste(round(median(a$Reads_MAPQ30), 1), paste0("(", round(median(a$Reads_MAPQ30_percentage), 1),"%)")), paste(round(median(a$Reads_nuclear), 1), paste0("(", round(median(a$Reads_nuclear_percentage), 1),"%)")), paste(round(median(a$Reads_in_peak), 1), paste0("(", round(median(a$Reads_in_peak_percentage), 1),"%)")), median(a$Insert_size_average))
qc_report = data.frame(report_keys, report_values)
print(unname(qc_report))
write.table(qc_report, paste0(outputdir, "/qc_report.txt"), row.names=F, col.names=F, quote=F, sep="\t")

qc_idx = (a$Reads_in_peak_percentage > percentage_threshold) & (a$Reads_properly_paired > reads_threshold)
if(sum(!qc_idx)!=0){
    write.table(a[qc_idx, ], paste0(outputdir, "/statistics_QCpassed.txt"), row.names=F, col.names=T, quote=F, sep="\t")
    write.table(a[!qc_idx, ], paste0(outputdir, "/statistics_QCfailed.txt"), row.names=F, col.names=T, quote=F, sep="\t")
} else {
    write.table(a, paste0(outputdir, "/statistics_QCpassed.txt"), row.names=F, col.names=T, quote=F, sep="\t")
}

