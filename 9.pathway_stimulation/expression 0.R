TV <- subset(test.merge,orig.ident == "Tet2mut_Veh")
TD <- subset(test.merge,orig.ident == "Tet2mut_DSS")
WV <- subset(test.merge,orig.ident == "WT_Veh")
WD <- subset(test.merge,orig.ident == "WT_DSS")
TV <- Tet2mut_Veh
###IFNG######
# 获取基因表达矩阵
gene_expression <- TV@assays$RNA@counts

genes_of_interest <- c("Ifngr1", "Ifngr2", "Ptpn2", "Jak2", "Sumo1", "Pias1", "Socs1", "Socs3", "Ifng")

# 找到这些基因在基因表达矩阵中的列索引
gene_indices <- which(rownames(gene_expression) %in% genes_of_interest)

gene_expression[gene_indices, ] <- 0

TV_IFNG <- TV
# 将修改后的基因表达矩阵重新赋值给 Seurat 对象
TV_IFNG@assays$RNA@counts <- gene_expression
save(TV_IFNG,file = "TV_IFNG.RData")
setwd("~/Tet2_DSS/TOSICA")

# 假设 gene 是你要查看的基因名称
gene_index <- which(rownames(gene_expression) == "Tyk2")

gene_expression_value <- gene_expression[gene_index, ]
gene_expression_value
###IFNA######
# 获取基因表达矩阵
gene_expression <- TV@assays$RNA@counts

genes_of_interest <- c("Ptpn6", "Ifnar1", "Ifnar2", "Ptpn1", "Usp18", "Tyk2", "Stat2", "Ptpn11", 
                       "Ifnb1", "Ifna13", "Ifna4", "Ifna12", "Ifna2", "Ifna16", "Ifna9", "Ifna1", 
                       "Ifna14", "Ifna15", "Ifna5", "Ifnab", "Ifna11", "Ifna7", "Ifna6",
                       "Ifngr1", "Ifngr2", "Ptpn2", "Jak2", "Sumo1", "Pias1", "Socs1", "Socs3", "Ifng")

# 找到这些基因在基因表达矩阵中的列索引
gene_indices <- which(rownames(gene_expression) %in% genes_of_interest)

gene_expression[gene_indices,] <- 0
TV_IFN <- TV
# 将修改后的基因表达矩阵重新赋值给 Seurat 对象
TV_IFN@assays$RNA@counts <- gene_expression
save(TV_IFN,file = "TV_IFN.RData")

###IL1R######
# 获取基因表达矩阵
gene_expression <- TV@assays$RNA@counts

genes_of_interest <- c("Irak2", "Il1rn", "Ikbkb", "Chuk", "Mapk14", "Irak1", "Ecsit", "Il1b", "Il1r1",
                       "Nfkb1", "Map2k3", "Rela", "Tgfb2", "Tgfb3", "Traf6", "Ifna1", "Ifnb1", "Il1a",
                       "Jun", "Myd88", "Nfkbia", "Tgfb1", "Map2k6", "Map3k1", "Tnf", "Mapk8", "Map3k14",
                       "Irak3", "Il6", "Il1rap", "Map3k7")


# 找到这些基因在基因表达矩阵中的列索引
gene_indices <- which(rownames(gene_expression) %in% genes_of_interest)

gene_expression[gene_indices,] <- 0
TV_TNF <- TV
# 将修改后的基因表达矩阵重新赋值给 Seurat 对象
TV_TNF@assays$RNA@counts <- gene_expression
save(TV_TNF,file = "TV_TNF.RData")


###TFNA######
# 获取基因表达矩阵
gene_expression <- TV@assays$RNA@counts

genes_of_interest <- c("Mcl1", "Cd80", "Traf1", "F2rl1", "Dusp2", "Tnc", "Fosl2", "Stat5a", "Vegfa", 
                       "Efna1", "Relb", "Rela", "Cebpd", "Ptger4", "Cdkn1a", "Ptx3", "Il15ra", "Atp2b1", 
                       "Nfkbia", "Tnf", "Ier3", "Ier2", "Hes1", "Tnfaip2", "Dusp1", "Eif1", "Fosl1", 
                       "Bcl6", "Ifngr2", "Tank", "Gadd45b", "Gadd45a", "Tubb2a", "Sqstm1", "Il18", "Rhob", 
                       "Cxcl1", "Btg2", "Nfe2l2", "Tsc22d1", "Irs2", "Atf3", "Nfil3", "Btg3", "Jag1", 
                       "Ackr3", "Cxcl5", "Phlda1", "Bhlhe40", "Per1", "Plk2", "Nfkb2", "Tnfsf9", "Tnfrsf9", 
                       "Klf10", "Tgif1", "Nfkbie", "Tnfaip6", "Tnfaip3", "Ninj1", "Birc3", "Birc2", "Smad3", 
                       "Dennd5a", "Socs3", "Phlda2", "Olr1", "Snn", "Bcl2a1d", "Egr3", "Sphk1", "G0s2", "Cd83", 
                       "Ccl20", "Klf9", "Msc", "Cflar", "Ier5", "Gfpt2", "Csf2", "Csf1", "Sgk1", "Cxcl2", "Ehd1", 
                       "Fjx1", "Klf4", "Klf2", "Tlr2", "Klf6", "Map2k3", "Map3k8", "Sdc4", "Cxcl10", "Nr4a1", 
                       "Nr4a2", "Nr4a3", "Icosl", "Nfat5", "Panx1", "Cxcl11", "Plek", "Rcan1", "Ripk2", "Dnajb4", 
                       "Plpp3", "Pnrc1", "Kynu", "Ifih1", "Dram1", "Ccrl2", "Spsb1", "Rnf19b", "Ccnl1", "Tnip1", 
                       "Ppp1r15a", "B4galt5", "Pdlim5", "Litaf", "Pmepa1", "Nampt", "Clcf1", "Il23a", "Zbtb10", 
                       "Slc16a6", "Trip10", "Tnfaip8", "Tiparp", "Pfkfb3", "Zc3h12a", "Tnip2", "Yrdc", "Gpr183", 
                       "Dusp4", "Rigi", "Slc2a6", "Trib1", "Kdm6b", "Dusp5", "Areg", "Bcl3", "Bmp2", "Btg1", 
                       "Ccnd1", "Cd44", "Cd69", "Cebpb", "F3", "Ccn1", "Serpinb8", "Edn1", "Egr1", "Egr2", "Ets2", 
                       "Fos", "Fosb", "Fut4", "Gch1", "B4galt1", "Slc2a3", "Hbegf", "Icam1", "Id2", "Il12b", "Il1a", 
                       "Il1b", "Il6", "Il6st", "Il7r", "Inhba", "Irf1", "Jun", "Junb", "Ldlr", "Lif", "Marcks", "Mxd1", 
                       "Maff")


# 找到这些基因在基因表达矩阵中的列索引
gene_indices <- which(rownames(gene_expression) %in% genes_of_interest)

gene_expression[gene_indices,] <- 0
TV_TNF <- TV
# 将修改后的基因表达矩阵重新赋值给 Seurat 对象
TV_TNF@assays$RNA@counts <- gene_expression

save(TV_TNF,file = "TV_TNF.RData")


###IFNG######
# 获取基因表达矩阵
gene_expression <- TV@assays$RNA@counts

genes_of_interest <- c("Ifngr1", "Ifngr2", "Ptpn2", "Jak2", "Sumo1", "Pias1", "Socs1", "Socs3", "Ifng",
                       "Ptpn6", "Ifnar1", "Ifnar2", "Ptpn1", "Usp18", "Tyk2", "Stat2", "Ptpn11", 
                       "Ifnb1", "Ifna13", "Ifna4", "Ifna12", "Ifna2", "Ifna16", "Ifna9", "Ifna1", 
                       "Ifna14", "Ifna15", "Ifna5", "Ifnab", "Ifna11", "Ifna7", "Ifna6")

# 找到这些基因在基因表达矩阵中的列索引
gene_indices <- which(rownames(gene_expression) %in% genes_of_interest)

gene_expression[gene_indices, ] <- 0

TV_TNF <- TV
# 将修改后的基因表达矩阵重新赋值给 Seurat 对象
TV_TNF@assays$RNA@counts <- gene_expression
save(TV_TNF,file = "TV_TNF.RData")
setwd("~/Tet2_DSS/TOSICA")


setwd("~/Tet2_DSS/TOSICA/TV_TNF")
# save metadata table:
TV_TNF$barcode <- colnames(TV_TNF)
TV_TNF$UMAP_1 <- TV_TNF@reductions$umap@cell.embeddings[,1]
TV_TNF$UMAP_2 <- TV_TNF@reductions$umap@cell.embeddings[,2]
write.csv(TV_TNF@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(TV_TNF, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(TV_TNF@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
write.csv(TV_TNF@reductions$umap@cell.embeddings, file='umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
write.csv(TV_TNF@assays$RNA@counts,file = "TV_TNF_counts.csv")
write.csv(counts_matrix,file = "TV_TNF_counts.csv")
# load datasets
#TV_TNF <- readRDS("pbmc3k_final.rds")

FeaturePlot(TV_TNF,features = "Jak2")


gene_expression <- TV@assays$RNA@counts
#Adra1d	Adrb3	Adra2a	Adrb1	Adra2c	Adrb2	Adra1a	Adra1b	Adra2b
genes_of_interest <- c("Adra1d", "Adrb3", "Adra2a", "Adrb1", "Adra2c", "Adrb2", "Adra1a", "Adra1b", "Adra2b")

# 找到这些基因在基因表达矩阵中的列索引
gene_indices <- which(rownames(gene_expression) %in% genes_of_interest)

gene_expression[gene_indices, ] <- 0
gene_expression[gene_indices, ]
TV_TNF <- TV
# 将修改后的基因表达矩阵重新赋值给 Seurat 对象
TV_TNF@assays$RNA@counts <- gene_expression
save(TV_TNF,file = "TV_TNF.RData")
setwd("~/Tet2_DSS/TOSICA")

# save metadata table:
TV_TNF$barcode <- colnames(TV_TNF)
TV_TNF$UMAP_1 <- TV_TNF@reductions$umap@cell.embeddings[,1]
TV_TNF$UMAP_2 <- TV_TNF@reductions$umap@cell.embeddings[,2]
write.csv(TV_TNF@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(TV_TNF, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(TV_TNF@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
write.csv(TV_TNF@reductions$umap@cell.embeddings, file='umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
write.csv(TV_TNF@assays$RNA@counts,file = "TV_TNF_counts.csv")
write.csv(counts_matrix,file = "TV_TNF_counts.csv")
# load datasets
#TV_TNF <- readRDS("pbmc3k_final.rds")

FeaturePlot(TV_TNF,features = "Jak2")

