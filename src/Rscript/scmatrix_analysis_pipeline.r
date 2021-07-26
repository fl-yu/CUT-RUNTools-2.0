# recieve the parameters
raw_mat_file <- commandArgs(trailingOnly = T)[1]
dr_npc <- commandArgs(trailingOnly = T)[2]
cluster_resolution <- commandArgs(trailingOnly = T)[3]
output_dir <- commandArgs(trailingOnly = T)[4]
experiment_name <- commandArgs(trailingOnly = T)[5]
path_to_rscript <- commandArgs(trailingOnly = T)[6]
pythonbin <- commandArgs(trailingOnly = T)[7]

setwd(output_dir)
source(paste0(path_to_rscript, "/filter_cells.r"))
source(paste0(path_to_rscript, "/filter_features.r"))
source(paste0(path_to_rscript, "/dim_reduction.r"))
source(paste0(path_to_rscript, "/do_pca.r"))
source(paste0(path_to_rscript, "/fast_tfidf.r"))
source(paste0(path_to_rscript, "/tfidf.r"))
source(paste0(path_to_rscript, "/multiplot.r"))
message(date(), " [info] All the function scripts were loaded ...")

# load package
message(date(), " [info] Load packages")
suppressPackageStartupMessages({
    library(reticulate)
    use_python(pythonbin, required=T)
    reticulate::import("leidenalg")
    reticulate::import("igraph")
    library(leiden)
    library(data.table)
    library(Matrix)
    library(irlba)
    library(Rtsne)
    library(RANN)
    library(igraph)
    library(uwot)
    library(ggplot2)
    })

# convert character parameter to numberic
dr_npc = as.numeric(dr_npc)
cluster_resolution = as.numeric(cluster_resolution)
message(date(), " [info] The countMatrix file ", raw_mat_file, " will be loaded ...")
raw_mat <- fread(raw_mat_file, header=T, sep="\t", showProgress=T)
binary_mat <- as.matrix(raw_mat[, -1])
rownames(binary_mat) <- as.data.frame(raw_mat[, 1])[, 1]
message(date(), " [info] The countMatrix has ", nrow(binary_mat), " features and ", ncol(binary_mat), " cells ...")
message(date(), " [info] The matrix is binarized and converted to sparseMatrix format")
binary_mat <- as(binary_mat, "sparseMatrix")
binary_mat@x[binary_mat@x >0] <- 1

# filter_cell_num=5
message(date(), " [info] The rare and ubiquitous features will be discarded for further analysis ")
rare_threshold <- max(2, floor(ncol(binary_mat)/1000))
ubiq_threshold <- floor(ncol(binary_mat)*0.8)
feature_detected <- rowSums(binary_mat)
feature_idx <- (feature_detected >= rare_threshold) & (feature_detected <= ubiq_threshold)
binary_mat <- binary_mat[feature_idx, ]


# message(date(), "[info] The matrix is filtered, ", nrow(binary_mat), " features and ", ncol(binary_mat), " cells were retained")
# message(date(), "[info] The matrix is normalized by TF-IDF method")
# tfidf2 <- tfidf(binary_mat)

# cv_filter = T
# threshold_cv = 20000
# if(cv_filter){
#     if(nrow(tfidf2)>=threshold_cv){
#         message(date(), " [info] Features were further filtered based on CV to ensure high signal noise ratio and reduce calculation complexity")
#         meantfidf2 <- rowMeans(tfidf2)
#         sdtfidf2 <- apply(tfidf2, 1, sd)
#         cv2tfidf2 <- sdtfidf2/meantfidf2
#         tfidf2 = tfidf2[order(cv2tfidf2, decreasing=TRUE), ][1:threshold_cv, ]
#         message(date(), " [info] Top ", threshold_cv, " Features were retained")
#     } else {
#         message(date(), " [info] Number of features less than ", threshold_cv)
#         message(date(), " [info] All features were retained")
#     }
# }
cv_filter = F
threshold_cv = 20000
if(cv_filter){
    if(nrow(binary_mat)>=threshold_cv){
        message(date(), " [info] Features were further filtered based on CV to ensure high signal noise ratio and reduce calculation complexity")
        meanbinary_mat <- rowMeans(binary_mat)
        sdbinary_mat <- apply(binary_mat, 1, sd)
        cv2binary_mat <- sdbinary_mat/meanbinary_mat
        binary_mat <- binary_mat[order(cv2binary_mat, decreasing=TRUE), ][1:threshold_cv, ]
        message(date(), " [info] Top ", threshold_cv, " Features were retained")
    } else {
        message(date(), " [info] Number of features less than ", threshold_cv)
        message(date(), " [info] All features were retained")
    }
}


# head(apply(binary_mat, 1, quantile))
# any(is.na(binary_mat[, 1370]))
message(date(), " [info] The matrix is filtered: ", nrow(binary_mat), " features and ", ncol(binary_mat), " cells were retained")
message(date(), " [info] The matrix is normalized by TF-IDF method")
# tfidf2 <- tfidf(binary_mat, log_TF=F)
tfidf2 <- tfidf(binary_mat)
# write.table(as.matrix(tfidf2), paste0("normlized_countMat_", experiment_name, ".txt"), row.names=T, col.names=T, quote=F, sep="\t")
message(date(), " [info] Dimensionality reduction method PCA is performed")
save(tfidf2, file="tfidf2.rda")
print(dr_npc)
pca_mat <- do_pca(tfidf2, dims=dr_npc)
message(date(), " [info] Dimensionality reduction method tSNE is performed")
tsne_mat = Rtsne(pca_mat, check_duplicates = FALSE, verbose = TRUE, max_iter = 1000, is_distance = FALSE, theta = 0.5, pca = FALSE)
message(date(), " [info] Dimensionality reduction method UMAP is performed")
umap_mat <- umap(pca_mat)

colnames(pca_mat) = paste0("PCA", 1:ncol(pca_mat))
colnames(tsne_mat$Y) = c("tSNE1", "tSNE2")
colnames(umap_mat) = c("UMAP1", "UMAP2")

## leiden cluster
knn.info<- nn2(pca_mat, k=dr_npc)
## convert to adjacancy matrix
knn <- knn.info$nn.idx
adj <- matrix(0, nrow(pca_mat), nrow(pca_mat))
rownames(adj) <- colnames(adj) <- rownames(pca_mat)
for(i in seq_len(nrow(pca_mat))){
    adj[i, rownames(pca_mat)[knn[i, ]]] <- 1
    }
xx <- graph_from_adjacency_matrix(adj)
adjacency_matrix <- igraph::as_adjacency_matrix(xx)
partition <- leiden(adjacency_matrix, resolution_parameter=cluster_resolution)
message(date(), " [info] The cluster number and cell distribution showed as follow")
print(table(partition))
cluster_num <- length(table(partition))
myColor = c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F", "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF")
message(date(), " [info] Dimensionality reduction plots generating")
pdf(paste0("Dimension reduction (PCA) of ", experiment_name, ".pdf"))
plot(pca_mat[, 1], pca_mat[, 2], pch=20, main="pca", cex=0.8, cex.main=1.5, frame.plot=FALSE, col=rev(myColor)[1:cluster_num][partition])
legend("topleft", pch=c(20, 20), col=rev(myColor)[1:cluster_num], legend=paste0("Cluster", unique(partition)), bty="o",  box.col="darkgreen", cex=.8)
invisible(dev.off())
pdf(paste0("Dimension reduction (tSNE) of ", experiment_name, ".pdf"))
plot(tsne_mat$Y, pch=20, main="tsne", cex=0.8, cex.main=1.5, frame.plot=FALSE, col=rev(myColor)[1:cluster_num][partition])
legend("topleft", pch=c(20, 20), col=rev(myColor)[1:cluster_num], legend=paste0("Cluster", unique(partition)), bty="o",  box.col="darkgreen", cex=.8)
invisible(dev.off())
pdf(paste0("Dimension reduction (UMAP) of ", experiment_name, ".pdf"))
plot(umap_mat[, 1], umap_mat[, 2], pch=20, main="umap", cex=0.8,  cex.main=1.5, frame.plot=FALSE, col=rev(myColor)[1:cluster_num][partition])
legend("topleft", pch=c(20, 20), col=rev(myColor)[1:cluster_num], legend=paste0("Cluster", unique(partition)), bty="o",  box.col="darkgreen", cex=.8)
invisible(dev.off())
message(date(), " [info] Output the normalized matrix and annotation files")
write.table(as.matrix(tfidf2), paste0("normlized_countMat_", experiment_name, ".txt"), row.names=T, col.names=T, quote=F, sep="\t")
cell_anno <- data.frame(colnames(binary_mat),partition)
colnames(cell_anno) <- c("cell_name", "leiden_cluster")
write.table(cell_anno, paste0("leiden_cluster_annotation.txt"), row.names=F, col.names=T, quote=F, sep="\t")

