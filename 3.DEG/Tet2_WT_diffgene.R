library(Seurat)
ident.colors <- c("#466983FF","#CE3D32FF","#F0E685FF","#5050FFFF","#BA6338FF","#6BD76BFF",
                  "#D595A7FF", "#5DB1DDFF","#749B58FF","#802268FF","#924822FF","#af2934",
                  "#ffe327","#2f4e87","#b0b9b8","#23452F","#aed4e9","#f4a69a")      
names(x = ident.colors) <- c("CMP", "MDP", "MEP", 'Macrophage', 'MP', 'T_cells', 'ErP', 'Erythroid cells', 'CLP', 'MkP', 'Myeloid cells', 'Basophil','Eosinophil', 'Pre_B', 'HSC', 'GMP', 'NP', 'Pro_B')
DimPlot(BM_big_harmony3, reduction = "umap",group.by = "cell_type")+  scale_color_manual(values = ident.colors)
DimPlot(BM_big_harmony3, reduction = "umap",group.by = "cell_type",split.by = "orig.ident",ncol = 2)+  scale_color_manual(values = ident.colors)

library(Seurat)
library(ggplot2)
library(dplyr)
# setwd("D:/KS项目/公众号文章/堆叠柱状图显示比例")

Ratio <- BM_big_harmony3@meta.data %>%group_by(orig.ident,cell_type) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)

ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cell_type))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = ident.colors)


Ratio <- BM_big_harmony3@meta.data %>%group_by(orig.ident,celltype) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)

ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = celltype))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#af2934","#ffe327","#2f4e87","#b0b9b8","#f0eedf",
                               "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
                               "#262a35","#c5942e","#a2a7ab"))


# Plot results
##########################

library(fgsea)
library(ggplot2)
library(msigdbr)
library(foreach)
library(doParallel)

library(Seurat)
BM_big_harmony3

#BM_big_harmony3 <- SCTransform(BM_big_harmony3, vars.to.regress = "percent.mt", verbose = FALSE) 
#BM_big_harmony3 <- SCTransform(BM_big_harmony3, 
                      method = "glmGamPoi", 
                      ncells = 20000, 
                      vars.to.regress = c("percent.mt","S.Score","G2M.Score"), 
                      verbose = T)
Idents(BM_big_harmony3) <- "cell_type"
deg_HSC=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="HSC")
deg_HSC$gene=rownames(deg_HSC)
deg_HSC$cluster="HSC"

deg_MDP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="MDP")
deg_MDP$gene=rownames(deg_MDP)
deg_MDP$cluster="MDP"

deg_CMP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="CMP")
deg_CMP$gene=rownames(deg_CMP)
deg_CMP$cluster="CMP"

deg_MEP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="MEP")
deg_MEP$gene=rownames(deg_MEP)
deg_MEP$cluster="MEP"

deg_Macrophage=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                           group.by = "orig.ident",subset.ident ="Macrophage")
deg_Macrophage$gene=rownames(deg_Macrophage)
deg_Macrophage$cluster="Macrophage"

deg_GMP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="GMP")
deg_GMP$gene=rownames(deg_GMP)
deg_GMP$cluster="GMP"


deg_MKP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="MkP")
deg_MKP$gene=rownames(deg_MKP)
deg_MKP$cluster="MkP"

deg_ErP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="ErP")
deg_ErP$gene=rownames(deg_ErP)
deg_ErP$cluster="ErP"


deg_Macrophage=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                           group.by = "orig.ident",subset.ident ="Macrophage")
deg_Macrophage$gene=rownames(deg_Macrophage)
deg_Macrophage$cluster="Macrophage"


deg_MP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                   group.by = "orig.ident",subset.ident ="MP")
deg_MP$gene=rownames(deg_MP)
deg_MP$cluster="MP"

deg_T_cells=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                      group.by = "orig.ident",subset.ident ="T_cells")
deg_T_cells$gene=rownames(deg_T_cells)
deg_T_cells$cluster="T_cells"

deg_Erythroid_cells=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                                group.by = "orig.ident",subset.ident ="Erythroid cells")
deg_Erythroid_cells$gene=rownames(deg_Erythroid_cells)
deg_Erythroid_cells$cluster="Erythroid cells"

deg_CLP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                    group.by = "orig.ident",subset.ident ="CLP")
deg_CLP$gene=rownames(deg_CLP)
deg_CLP$cluster="CLP"

deg_Myeloid_cells=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                  group.by = "orig.ident",subset.ident ="Myeloid cells")
deg_Myeloid_cells$gene=rownames(deg_Myeloid_cells)
deg_Myeloid_cells$cluster="Myeloid_cells"

deg_Basophil=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
             group.by = "orig.ident",subset.ident ="Basophil")
deg_Basophil$gene=rownames(deg_Basophil)
deg_Basophil$cluster="Basophil"

deg_Eosinophil=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                         group.by = "orig.ident",subset.ident ="Eosinophil")
deg_Eosinophil$gene=rownames(deg_Eosinophil)
deg_Eosinophil$cluster="Eosinophil"

deg_Pre_B=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                      group.by = "orig.ident",subset.ident ="Pre_B")
deg_Pre_B$gene=rownames(deg_Pre_B)
deg_Pre_B$cluster="Pre_B"

#deg_Pro_B=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
#                      group.by = "orig.ident",subset.ident ="Pro_B")
#deg_Pro_B$gene=rownames(deg_Pro_B)
#deg_Pro_B$cluster="Pro_B"

deg_NP=FindMarkers(BM_big_harmony3,ident.1 = "Tet2mut_Veh",ident.2 = "WT_Veh",
                           group.by = "orig.ident",subset.ident ="NP")
deg_NP$gene=rownames(deg_NP)
deg_NP$cluster="NP"


clusterdeg <- rbind(deg_HSC,deg_CMP,deg_MDP,deg_MEP,deg_Macrophage,deg_MP,
                    deg_T_cells,deg_ErP,deg_Erythroid_cells,deg_CLP,
                    deg_MKP,deg_Myeloid_cells,deg_Basophil,deg_Eosinophil,
                    deg_Pre_B,deg_GMP)
clusterdeg <- rbind(deg_Immature_B_cells,deg_MKP,deg_CMP,deg_MKP/ErP,deg_CLP,deg_GMP,
                    deg_Macrophage,deg_MEP,deg_CDP,deg_HSC,deg_NP)

color_cluster=c("#af2934","#ffe327","#2f4e87","#b0b9b8","#f0eedf",
                "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
                "#262a35","#c5942e","#a2a7ab") #颜色可以自己定义，ggsci和RColorBrewer包的颜色搭配就很NICE了
library(scRNAtoolVis)
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', 
               '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', 
               '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E',
               '#68A180', '#3A6963', '#968175')#颜色设置 
names(color_cluster)= c("HSC","MEP","MPP_1","GMP","CMP","MPP_2","Macrophage","MKP_ErP","Immature B","MDP")

ident.colors <- c("#466983FF","#CE3D32FF","#F0E685FF","#5050FFFF","#BA6338FF","#6BD76BFF",
                  "#D595A7FF", "#5DB1DDFF","#749B58FF","#802268FF","#924822FF","#af2934",
                  "#ffe327","#2f4e87","#b0b9b8","#23452F","#aed4e9","#f4a69a")      
names(x = ident.colors) <- c("CMP", "MDP", "MEP", 'Macrophage', 'MP', 'T_cells', 'ErP', 'Erythroid cells', 'CLP', 'MkP', 'Myeloid cells', 'Basophil','Eosinophil', 'Pre_B', 'HSC', 'GMP', 'NP', 'Pro_B')

library(scRNAtoolVis)
clusterdeg_1= clusterdeg[clusterdeg$p_val_adj< 0.05, ] 
p1<-jjVolcano(diffData = clusterdeg_1,
              base_size = 20,
              log2FC.cutoff = 0.5, # logfc根据自己情况设定
              pvalue.cutoff = 0.05,
              size  = 3, #设置点的大小
              #aesCol = c('blue','orange'), #设置点的颜色
              tile.col = ident.colors, #设置cluster的颜色
              #col.type = "adjustP", #设置图例显示方式
              topGeneN = 5, #设置展示topN的基因
              flip = F,
              legend.position = c(4,4),
)
p1
filtered_data <- data[!grepl(keyword,data$column_name), ]

clusterdeg_1 =clusterdeg[!grepl("MT-*",clusterdeg$gene), ]

deg_CMP_1 =deg_CMP[!grepl("MT-*",deg_CMP$gene), ]
deg_CDP_1 =deg_CDP[!grepl("MT-*",deg_CDP$gene), ]
deg_MDP_1 =deg_CDP[!grepl("MT-*",deg_MDP$gene), ]
deg_MKP_ErP_1 =deg_MKP_ErP[!grepl("MT-*",deg_MKP_ErP$gene), ]
deg_GMP_1 =deg_GMP[!grepl("MT-*",deg_GMP$gene), ]
deg_HSC_1 =deg_HSC[!grepl("MT-*",deg_HSC$gene), ]
deg_Immature_B_cells_1 =deg_Immature_B_cells[!grepl("MT-*",deg_Immature_B_cells$gene), ]
deg_MEP_1 =deg_MEP[!grepl("MT-*",deg_MEP$gene), ]
deg_Macrophage_1 =deg_Macrophage[!grepl("MT-*",deg_Macrophage$gene), ]