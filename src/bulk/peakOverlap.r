# Further customization

#!/n/app/R/3.3.3/bin/Rscript
peak_num1=as.numeric(commandArgs(trailingOnly = T)[1])
peak_num2=as.numeric(commandArgs(trailingOnly = T)[2])
peak_num3=as.numeric(commandArgs(trailingOnly = T)[3])
outputdir=commandArgs(trailingOnly = T)[4]
# peak_num1=100
# peak_num2=100
# peak_num3=50
# outputdir="/Users/fyu/Documents/GitHub/CUT-RUNTools-2.0/test"
# if (!require(ggvenn)){
#     if (!require(devtools)) {
#         install.packages("devtools")
#     }
#     devtools::install_github("yanlinlin82/ggvenn")
# }

if (!require(VennDiagram)){
    install.packages("VennDiagram")
}

set1 <- c(paste("a", seq(1, peak_num1)), seq(1, peak_num3))
set2 <- c(paste("b", seq(1, peak_num2)), seq(1, peak_num3))
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

pdf(paste0(outputdir, "/Vennplot_for_peak_overlap.pdf"))
vennplot_peak <- display_venn(
        x=list(set1, set2),
        category.names = c("Peak set 1" , "Peak set 2"),
        # Circles
        lwd = 2,
        lty = 'blank',
        # fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
        fill = c("#999999", "#E69F00"),
        # Numbers
        cex = .9,
        fontface = "italic",
        # Set names
        cat.cex = 1,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        main = "Venn plot for the peak overlap"
        # imagetype="pdf",
        # filename = paste0(outputdir, "/Vennplot_for_peak_overlap.pdf"),
        ,output=T
)
dev.off()
message("[INFO] The results can be found in the ", outputdir)
