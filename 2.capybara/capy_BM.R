# Construction of a high-resolution reference

rm(list=ls())
#setwd("~/") # http://192.168.31.163:8787/
#setwd("~/downloads_01/Capybara_reproducibility-main")
#rm(list=ls())

setwd("~/downsample_HehangDSS")
library("devtools")
devtools::install_github("morris-lab/Capybara")
#library(splatter)
library(Capybara)
library(scater)
library(mixdist)
library(MASS)
require(pastecs)
library(Seurat)
library(ggsci)
setwd("~/hh028/Tet2_DSS_scRNAseq/BM/240204/capybara")
source("~/hh028/Tet2_DSS_scRNAseq/BM/240204/capybara/function.R")

load("HSPC_downsample_harmony_rename_final.RData")

### STOP
test1.seu<- subset(test.merge, subset=orig.ident=="WT_Veh")
test2.seu<- subset(test.merge, subset=orig.ident=="Tet2mut_Veh")
test3.seu<- subset(test.merge, subset=orig.ident=="WT_DSS")
test4.seu<- subset(test.merge, subset=orig.ident=="Tet2mut_DSS")

library(scales)
library(ggsci)
library(Seurat)
mycolor2<- pal_simpsons()(16)
mycolor2
show_col(mycolor2)

mycolor<- c("#FED439FF", "#709AE1FF", "#8A9197FF", "blue",
            "#FD7446FF",  "#197EC0FF", "#D2AF81FF", 
            "#71D0F5FF", "#370335FF", "#075149FF",
            "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
show_col(mycolor)
ref.seu<- test1.seu ## WT AS Ref
DimPlot(test1.seu, reduction = "umap", group.by = "cell_type", cols = mycolor)

DimPlot(test1.seu, reduction = "umap", group.by = "cell_type", cols = mycolor)
DimPlot(test2.seu, reduction = "umap", group.by = "cell_type", cols = mycolor)
DimPlot(test3.seu, reduction = "umap", group.by = "cell_type", cols = mycolor)
DimPlot(test4.seu, reduction = "umap", group.by = "cell_type", cols = mycolor)


count<- GetAssayData(object = ref.seu, layer = "counts") ## slot or layer 

count.test1<- GetAssayData(object = test1.seu, layer = "counts") ## slot or layer 
count.test2<- GetAssayData(object = test2.seu, layer = "counts") ## slot or layer 
count.test3<- GetAssayData(object = test3.seu, layer = "counts") ## slot or layer 
count.test4<- GetAssayData(object = test4.seu, layer = "counts") ## slot or layer 

metaDF<- ref.seu@meta.data
metaDF[1:4,]

# construct.high.res.reference
ref.list <- construct.high.res.reference(count, coldata.df = metaDF, criteria = "cell_type", cell.num.for.ref = 2000) ## adjust cell.num.for.ref
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.df <- ref.construction(ref.list[[1]], ref.list[[2]], "cell.type")


single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 10, save.to.path = "./", save.to.filename = "reference_WT", unix.par = TRUE)

# test1 to test4
single.round.QP.analysis(ref.df, count.test1, n.cores = 10, save.to.path = "./", save.to.filename = "WT_querry", unix.par = TRUE) ## 3 minutes
single.round.QP.analysis(ref.df, count.test2, n.cores = 10, save.to.path = "./", save.to.filename = "TET2_querry", unix.par = TRUE) ## 3 minutes
single.round.QP.analysis(ref.df, count.test3, n.cores = 10, save.to.path = "./", save.to.filename = "WTDSS_querry", unix.par = TRUE) ## 3 minutes
single.round.QP.analysis(ref.df, count.test4, n.cores = 10, save.to.path = "./", save.to.filename = "TET2DSS_querry", unix.par = TRUE) ## around 10 minutes


# Read in background and testing identity scores
background.mtx <- read.csv("./Hehang_reference_WT_scale.csv", header = T, row.names = 1, stringsAsFactors = F)


mtx.test <- read.csv("WT_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./TET2_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./WTDSS_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./TET2DSS_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
dim(mtx.test)
metaDF<- test1.seu@meta.data #CHANGE META DATA ACCORDINGLY
metaDF<- test2.seu@meta.data #CHANGE META DATA ACCORDINGLY
metaDF<- test3.seu@meta.data #CHANGE META DATA ACCORDINGLY
metaDF<- test4.seu@meta.data #CHANGE META DATA ACCORDINGLY
dim(mtx.test)
mtx.test[1:4,]
#View(mtx.test)

metaDF_new<- cbind(metaDF,mtx.test)
#metaDF_new[1:4,]
#test1.seu@reductions$umap$UMAP_1

test1.seu@reductions$umap@cell.embeddings[1:4,]

metaDF_new_final<- cbind(metaDF_new,test4.seu@reductions$umap@cell.embeddings)
metaDF_new_final[1:4,]
dim(metaDF_new_final)

library(RColorBrewer) 
##Run good
ggplot(metaDF_new_final, aes(x=UMAP_1, y=UMAP_2, color=frxn_cell.type_CLP))+
  geom_point(size=1, alpha=0.9, shape=20)+
  scale_color_gradient(low="gray",high = "#990000")+
  theme_classic()
?shuf()
devtools::install_github("caleblareau/BuenColors")
library(BuenColors)
pdf("./cell_type_prob_heatmap/WT_CLP.pdf", width = 8, height = 8)
ggplot(shuf(metaDF_new_final), aes(x=UMAP_1, y=UMAP_2, color=frxn_cell.type_HSC1))+
  geom_point(size=1, alpha=0.9, shape=20)+
  scale_color_gradient(low="gray",high = "#990000")+
  theme_classic()+
  theme(legend.title = element_blank())+
  ggtitle("WT_HSC1_prob")
dev.off()



ggplot(shuf(metaDF_new_final), aes(x=frxn_cell.type_CLP))+
  stat_ecdf(geom = "step")+
  theme_classic()



ggplot(metaDF_new_final, aes(x=frxn_cell.type_HSC1))+
  stat_ecdf(geom = "step")+
  xlim(0,1)+
  theme_classic()
ggplot(metaDF_new_final, aes(x=frxn_cell.type_HSC2))+
  stat_ecdf(geom = "step")+
  xlim(0,1)+
  theme_classic()
ggplot(metaDF_new_final, aes(x=frxn_cell.type_CLP))+
  stat_ecdf(geom = "step")+
  xlim(0,1)+
  theme_classic()


ggplot(metaDF_new_final, aes(x=frxn_cell.type_ErP))+
  stat_ecdf(geom = "step")+
  xlim(0,1)+
  theme_classic()

ggplot(metaDF_new_final, aes(x=frxn_cell.type_MEP))+
  stat_ecdf(geom = "step")+
  xlim(0,1)+
  theme_classic()

FeaturePlot()+
  scale_color_aaas()



## GOOD one


metaDF_new_final[1:4,]
ggplot(metaDF_new_final, aes(x=frxn_cell.type_Pro_Mast, color=cell_type))+
  stat_ecdf()+
  scale_color_manual(values=mycolor)+
  xlim(0,1)+
  theme_classic()

=====================
  =====================
  =====================
  
  
  +
  scale_color_gradient(low="gray",high = "#990000")+
  theme_classic()+
  theme(legend.title = element_blank())+
  ggtitle("WT_CLP_prob")




reorder(data$Pathway,data$Count,decreasing=TRUE)
table(metaDF_new_final$cell_type)
"#370335FF"
metaDF_new_final<- metaDF_new_final[order(metaDF_new_final$frxn_cell.type_CLP),]

## run GOOD
pdf("./cell_type_prob_heatmap/WT_CLP.pdf", width = 8, height = 8)
ggplot(metaDF_new_final, aes(x=UMAP_1, y=UMAP_2, color=frxn_cell.type_CLP))+
  geom_point(size=2, alpha=0.7, shape=20)+
  scale_color_gradient(low="#71D0F5FF",high = "#C80813FF")+
  theme_classic()+
  theme(legend.title = element_blank())+
  ggtitle("WT_CLP_prob")
dev.off()







#library(RColorBrewer)  
ggplot(metaDF_new_final, aes(x=UMAP_1, y=UMAP_2, color=frxn_cell.type_CLP))+
  geom_point(size=2, alpha=0.7, shape=20)+
  #scale_color_gradient(low="gray",high = "#990000")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme_classic()

install.packages("ggthemes")
install.packages("ggtech")
library(ggthemes)
## Run good
ggplot(metaDF_new_final, aes(x=UMAP_1, y=UMAP_2, color=frxn_cell.type_HSC1))+
  geom_point(size=1)+
  #scale_color_gradient(low="gray",high = "#990000")+
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  #scale_color_viridis_c(begin = 0.15, end = 0.85, option = "C") +
  #geom_jitter(width=0.25, height=0.25)+
  scale_color_viridis_c(option = "C") +
  ggtitle("WT_on_Reference")+
  theme_classic()
#theme_economist()
#theme_solid()

# 安装并加载所需的包
install.packages("umap")
install.packages("ggplot2")
library(umap)
library(ggplot2)

# 读取数据，假设数据包含细胞ID、UMAP坐标和细胞类型等列
# 假设您的数据框名为 'data'，列名为 'cell_id', 'UMAP1', 'UMAP2', 'cell_type'
# 替换 'your_data.csv' 为您的数据文件路径

####################
# 创建UMAP对象
umap_obj <- umap(data[, c("UMAP1", "UMAP2")])

# 将UMAP结果添加到原始数据框
data_umap <- cbind(data, umap_obj$layout)

# 绘制UMAP图并高亮显示特定细胞类型
ggplot(metaDF_new_final, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 1) +
  labs(title = "HSC1 in Seurat",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal()+
  theme_classic()

# 绘制所有细胞的UMAP图
p <- ggplot(test2.seu, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "lightgrey", size = 1) +
  labs(title = "UMAP Plot with Highlighted HSC1",
       x = "UMAP1", y = "UMAP2")+
  theme_classic()
p
# 高亮显示特定细胞类型
highlighted_cells <- subset(test2.seu, cell_type == "MEP")
p1 <- p + geom_point(data = highlighted_cells, aes(color = cell_type), size = 1)
p1
# 显示图形
print(p)


# 绘制所有细胞的UMAP图
p <- ggplot(metaDF_new_final, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "lightgrey", size = 1) +
  labs(title = "UMAP Plot with Highlighted Pro_Mast",
       x = "UMAP1", y = "UMAP2")+
  theme_classic()
p
# 高亮显示特定细胞类型
highlighted_cells <- subset(metaDF_new_final, cell_type == "Pro_Mast")
p1 <- p + geom_point(data = highlighted_cells, aes(color = cell_type), size = 1)
p1
# 显示图形
print(p)


#FeaturePlot(test1.seu, reduction = "umap", features="Egr1", order=TRUE)
table(metaDF_new_final$cell_type)
#write.csv(metaDF_new_final, file="TET2DSS_metaDF_new_final.csv")
WT <- subset(test.merge,orig.ident == "WT_Veh")
FeaturePlot(WT,reduction = "umap",features = "Nr4a1")
FeaturePlot(WT,reduction = "umap",features = "Cd34")
FeaturePlot(WT,reduction = "umap",features = "Elane")
FeaturePlot(WT,reduction = "umap",features = "Itga2b")
FeaturePlot(WT,reduction = "umap",features = "Alox5")

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(viridis)
library(scCustomize)
#viridis_plasma_light_high
#viridis_magma_dark_high
#viridis_magma_light_high
#viridis_inferno_dark_high
#viridis_inferno_light_high
#viridis_dark_high
#viridis_light_high
gene<-c("Nr4a1","Cd34","Elane","Itga2b","Alox5")
gene<-c("Nr4a1")

i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = WT, 
                                  colors_use = viridis_inferno_light_high, 
                                  features = gene[i])+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
library(patchwork)
p<-wrap_plots(plots, ncol = 3);p
wrap_plots(plots)
ggsave(p,file="featureplot2.pdf",width = 8,height = 6)
scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

background.mtx[1:4,]
col.sub <- ncol(background.mtx) - 2
ncol(background.mtx)
col.sub

# Conduct reference randomization to get empirical p-value matrix
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])
save(ref.perc.list, file= "./ref.perc.list.RData")




# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))


# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.list[[2]], perc.ls = perc.list)
# Classificationn
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
rownames(classification) <- classification$barcode



multi.classification.list <- multi.id.curate.qp(binary.counts = bin.count, classification = classification, qp.matrix = mtx.test)
# Reassign variables
actual.multi <- multi.classification.list[[1]]
new.classification <- multi.classification.list[[2]]
new.classification[1:4,]

## VISUALIZATION


dim(metaDF)
metaDF[1:4,]
metaDF$call<- new.classification$call

actual.multi$qp.score[1:4]
actual.multi[1:4,]
#View(actual.multi)
dim(actual.multi)
score.df <- transition.score(actual.multi)
score.df[1:4,]
score.df
#write.csv(score.df, file="score_TET2dss.csv")

rslt <- table(metaDF$call, metaDF$cell_type)
rslt <- as.data.frame(apply(rslt, 2, function(x) round(x * 100/sum(x), digits = 3)))
rownames(rslt) <- paste0("Capy.", rownames(rslt))
colnames(rslt) <- paste0("Actual.", colnames(rslt))
rslt$capy <- rownames(rslt)
rslt.stk <- reshape2::melt(rslt)


rslt.stk[1:4,]


rslt.stk[1:4,]
rslt.stk
write.csv(rslt.stk, file="./for_aTET2_stacking.csv")
table(metaDF_new_final$cell_type)
#rslt.stk$capy 

rslt.stk$capy <- factor(rslt.stk$capy, levels = c("Capy.HSC1", 
                                                  "Capy.CLP",
                                                  "Capy.MEP",
                                                  "Capy.HSC2", 
                                                  "Capy.MDP", 
                                                  "Capy.NMP",
                                                  "Capy.Pro_NE",
                                                  "Capy.Pro_Mono",
                                                  "Capy.Pro_Mast",
                                                  "Capy.Pro_DC",
                                                  "Capy.Pro_Mk",
                                                  "Capy.ErP",
                                                  "Capy.Pro_B"), ordered = T)


rslt.stk$variable <- factor(rslt.stk$variable, levels = c("Actual.HSC1", 
                                                          "Actual.CLP",
                                                          "Actual.MEP",
                                                          "Actual.HSC2", 
                                                          "Actual.MDP", 
                                                          "Actual.NMP",
                                                          "Actual.Pro_NE",
                                                          "Actual.Pro_Mono",
                                                          "Actual.Pro_Mast",
                                                          "Actual.Pro_DC",
                                                          "Actual.Pro_Mk",
                                                          "Actual.ErP",
                                                          "Actual.Pro_B"), ordered = T)






ggplot(rslt.stk, aes(x = variable, y = capy, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(begin = 0.15, end = 0.85, option = "A") +
  labs(x = "Actual Annotation", y = "Capybara Annotation") +
  ggtitle("Simulation Classification Result") +
  theme(legend.position="right",
        axis.text.x = element_text(face = "bold.italic", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold.italic", size = 12),
        axis.title = element_text(face = "bold.italic", size = 14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        title = element_text(face = "bold.italic", size = 14),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_blank())

ggplot(rslt.stk, aes(x = variable, y = capy, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "A") +
  labs(x = "Actual Annotation", y = "Capybara Annotation") +
  ggtitle("Simulation Classification Result") +
  theme(legend.position="right",
        axis.text.x = element_text(face = "bold.italic", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold.italic", size = 12),
        axis.title = element_text(face = "bold.italic", size = 14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        title = element_text(face = "bold.italic", size = 14),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_blank())


# #Visualized in tilePlot
# pdf("./WT.pdf", width = 12, height = 12)     ##这里，需要改下面的编号
# ggplot(rslt.stk, aes(x = variable, y = capy, fill = value)) +
#   geom_tile(
#     color = "#FAFAFA", # 设置小格子边框颜色
#     lwd = 0.8, # 设置小格子的边框宽度
#     linetype = 1) + # 设置小格子的边框线型+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1))+
#   scale_fill_viridis_c(option = "C")+
#   ggtitle("WT_transition plot")   
# dev.off()

pdf("./WTdss_transitionMatrix.pdf", width = 12, height = 12)
ggplot(rslt.stk, aes(x = variable, y = capy, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "A") +
  labs(x = "Seurat Cell Type", y = "Capybara Transition Potential") +
  ggtitle("Simulation Classification Result") +
  theme(legend.position="right",
        axis.text.x = element_text(face = "bold.italic", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold.italic", size = 12),
        axis.title = element_text(face = "bold.italic", size = 14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        title = element_text(face = "bold.italic", size = 14),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_blank())
dev.off()

## WTDSS


df1<- read.csv("score_WT.csv")
df2<- read.csv("score_TET2.csv")
df3<- read.csv("score_WTdss.csv")
df4<- read.csv("score_TET2dss.csv")

df_score<- cbind(df1, df2, df3, df4)
row.names(df_score)<- df_score$X
df_score[1:4,]
write.csv(df_score, file="df_score.csv")

colnames(df_score)
df_score<- df_score[,c(3,6,9,12)]
colnames(df_score)<- c("WT_Veh", "TET2_Veh", "WT_DSS", "TET2_DSS")


library(pheatmap)
pheatmap(df_score, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "column",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("Navy", "white", "firebrick3"))(60))

df_score
write.csv(df_score, file="df_score_orig.csv")
pheatmap(df_score, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "row",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("Navy", "white", "firebrick3"))(60))




df<- read.csv("df_score_norm.csv")
df
row.names(df)<- df$X
df<- df[,-1]
df

pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "column",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("Navy", "white", "firebrick3"))(60))

df1<- log2(df)
df1

pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("Navy", "white", "firebrick3"))(60))


pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("#276419", "white", "#8E0152"))(60))


pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("#276419", "white", "#8E0152"))(60))



pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("#276419", "white", "#8E0152"))(60))


pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "gray",
         colorRampPalette(c("#5E4FA2", "gray", "#9E0142"))(60))

pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "white",
         colorRampPalette(c("#5E4FA2", "gray", "#9E0142"))(60))

pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "grey",
         colorRampPalette(c("#5E4FA2", "white", "#9E0142"))(60))


pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "grey",
         colorRampPalette(c("blue", "white", "red"))(60))

pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "grey",
         colorRampPalette(c("#5E4FA2", "white", "#9E0142"))(60))

pheatmap(df1, cellwidth = 30, cellheight = 30,
         show_rownames = T,
         cluster_rows = F,
         cluster_col =F,
         scale = "none",
         #breaks =bk, #标尺赋值
         display_numbers = T,
         na_col = "white",
         border_color = "grey",
         colorRampPalette(c("#3C5488B2", "white", "#E64B35B2"))(30))


library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal

library("scales")
show_col(mypal)


library(RColorBrewer)
brewer.pal(9, "Set1")
brewer.pal.info
brewer.pal(11,"PiYG")
brewer.pal(11, "Spectral")
library(scales)
show_col(brewer.pal(11, "Spectral"))

df<- read.csv("for_aWT_stacking.csv")
df<- df[,-1]
df
#df_matrix<- table(df$variable, df$capy)
#df_matrix

ggplot(df, aes(x = variable, y = capy, fill = value)) +
  geom_tile(
    color = "#FAFAFA", # 设置小格子边框颜色
    lwd = 0.8, # 设置小格子的边框宽度
    linetype = 1) + # 设置小格子的边框线型+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_viridis_c(option = "C")

ggplot(df, aes(x = variable, y = capy, fill = value)) +
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_viridis_c(option = "C")


ggplot(df, aes(x = variable, y = capy, fill = value)) +
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_viridis_c(option = "C")

df


df$capy <- factor(df$capy, levels = c("Capy.HSC1", 
                                      "Capy.CLP",
                                      "Capy.MEP",
                                      "Capy.HSC2", 
                                      "Capy.MDP", 
                                      "Capy.NMP",
                                      "Capy.Pro_NE",
                                      "Capy.Pro_Mono",
                                      "Capy.Pro_Mast",
                                      "Capy.Pro_DC",
                                      "Capy.Pro_Mk",
                                      "Capy.ErP",
                                      "Capy.Pro_B"), ordered = T)


df$variable <- factor(df$variable, levels = c("Actual.HSC1", 
                                              "Actual.CLP",
                                              "Actual.MEP",
                                              "Actual.HSC2", 
                                              "Actual.MDP", 
                                              "Actual.NMP",
                                              "Actual.Pro_NE",
                                              "Actual.Pro_Mono",
                                              "Actual.Pro_Mast",
                                              "Actual.Pro_DC",
                                              "Actual.Pro_Mk",
                                              "Actual.ErP",
                                              "Actual.Pro_B"), ordered = T)

ggplot(df, aes(x = variable, y = capy, fill = value)) +
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_viridis_c(option = "A")

df

metaDF_new<- metaDF_new_final
#metaDF_new<- read.csv(file="./plots_WT/WT_metaDF_new.csv")
metaDF_new[1:4,]
metaDF_new$total<- metaDF_new$frxn_cell.type_HSC1+metaDF_new$frxn_cell.type_HSC2
metaDF_new$bias<- metaDF_new$frxn_cell.type_HSC1-metaDF_new$frxn_cell.type_HSC2
metaDF_new[1:4,]

mtx.test <- read.csv(".~/capybara_HehangDSS/qp_files/WT_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./TET2_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./WTDSS_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./TET2DSS_querry_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
dim(mtx.test)
metaDF<- test1.seu@meta.data #CHANGE META DATA ACCORDINGLY
metaDF<- test2.seu@meta.data #CHANGE META DATA ACCORDINGLY
metaDF<- test3.seu@meta.data #CHANGE META DATA ACCORDINGLY
metaDF<- test4.seu@meta.data #CHANGE META DATA ACCORDINGLY
dim(mtx.test)
mtx.test[1:4,]
#View(mtx.test)

#metaDF_new<- cbind(metaDF,mtx.test)
metaDF_new_final<- cbind(metaDF_new,test1.seu@reductions$umap@cell.embeddings)
bias_df<- metaDF_new_final
bias_df<- bias_df[bias_df$total>0.5,]
bias_df[1:4,]
dim(bias_df)
bias_df$bias<- bias_df$bias/bias_df$total
bias_df$bias_Z<- scale(bias_df$bias)
bias_df[1:4,]
write.csv(bias_df,file = "WT_HSC_GMP_bais.csv")
library(RColorBrewer)
library(ggplot2)
#brewer.pal.info
pdf("./plots_TET2DSS/TET2DSS_HSC1_HSC2_biasZ.pdf", width = 10, height = 10) 
ggplot()+
  geom_point(data=metaDF_new_final, aes(x=UMAP_1, y=UMAP_2),size=1, alpha=0.7, shape=20, color="grey")+
  geom_point(data=bias_df, aes(x=UMAP_1, y=UMAP_2, color=bias_Z),size=1, alpha=0.7, shape=20)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")))+
  theme_classic()+
  ggtitle("HSC1 vs. HSC2 bias in WT")
dev.off()
