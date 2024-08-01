cellinfo <- HSC1@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
cellinfo <- HSC2@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
colnames(cellinfo)=c('celltype', 'sample','nGene' ,'nUMI')
sample <- as.data.frame(subset(cellinfo,select = 'sample'))
selectedResolution <- "sample"

cellsPerGroup <- split(rownames(sample), 
                       sample[,selectedResolution])
#去除extend的TF
# sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
# dim(sub_regulonAUC)

#计算平均表达
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

#scale处理\类似于热图数据的标准化
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 


regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled_select = regulonActivity_byGroup_Scaled[rownames(regulonActivity_byGroup_Scaled) %in% c("Egr1(+)","Egr1(-)","Egr2(+)","Egr3(+)","Hoxb2(+)","Fos(+)","Fos(-)",
                                                                                                                       "Thra(+)","Myb(+)","Foxd2(+)","Irf6(-)","Gata3(+)","Gata3(-)",
                                                                                                                       "Gata1(+)","Gata1(-)","Gata2(+)","Gata2(-)","Foxo3(+)","Runx1(+)","Runx1(-)",
                                                                                                                       "Hoxa10(+)","Meis1(+)","Erg(+)","Rarb(+)","Irf6(+)","Irf8(+)"),]
regulonActivity_byGroup_Scaled_select
#热图
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byGroup_Scaled_select, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=12),
                                   show_row_names = T)) 
hm
