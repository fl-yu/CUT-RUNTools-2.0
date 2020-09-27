#` Dim reduction (tSNE/UMAP) and clustering given PCA space via Seurat
#` @param count_matrix    (sparse matrix)  Matrix to store in Seurat object
#` @param cell_embeddings (matrix)         Typically PCA coordinates of cells but could be any set of reduced dim coordinates
#` @param dims            (integer)        vector of dims to use from cell_embeddings in downstream analysis
#` @param metadata        (dataframe)      dataframe of metadata (rowonames are cell names) to add to Seurat object
#` @param reduction       (character)      reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#` @return                (Seurat object)  A Seurat object 

dim_reduction <- function(count_matrix, cell_embeddings, dims, metadata=NULL, reduction='pca') {
  if (is.null(metadata)) {
    seurat_obj <- Seurat::CreateSeuratObject(count_matrix, assay = "CutTag_all_bulk_peak", project = "scCutTag_HM")
  } else {
    seurat_obj <- Seurat::CreateSeuratObject(count_matrix, meta.data = metadata, assay = "CutTag_all_bulk_peak", project = "scCutTag_HM")
    print(2)
  }
    print(seurat_obj)
  seurat_obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', assay = "CutTag_all_bulk_peak")
  print(1)
  seurat_obj <- Seurat::L2Dim(seurat_obj, reduction='pca')
  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = reduction, dims = dims)
  seurat_obj <- Seurat::RunTSNE(seurat_obj, reduction = reduction, dims = dims)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction=reduction, nn.eps=0.25, dims=dims)
  return(seurat_obj)
}
