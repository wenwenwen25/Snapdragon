# load datasets
#test.merge <- readRDS("pbmc3k_final.rds")

# save metadata table:
test.merge$barcode <- colnames(test.merge)
test.merge$UMAP_1 <- test.merge@reductions$umap@cell.embeddings[,1]
test.merge$UMAP_2 <- test.merge@reductions$umap@cell.embeddings[,2]
write.csv(test.merge@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(test.merge, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(test.merge@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
write.csv(test.merge@reductions$umap@cell.embeddings, file='umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
write.csv(test.merge@assays$RNA@counts,file = "test.merge_counts.csv")
# load datasets
#test.merge <- readRDS("pbmc3k_final.rds")