HM <- subset(scRNA_2,group == "HM")
LM <- subset(scRNA_2,group == "LM")

write.csv(t(as.matrix(HM@assays$RNA@counts)),file = "HM_sce_exp.csv")
write.csv(t(as.matrix(LM@assays$RNA@counts)),file = "LM_sce_exp.csv")
# cellInfo <- human_data@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
# colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
# head(cellInfo)
# write.csv(cellInfo, file = "cellInfo.csv")



#############################################################################
#                                在R中可视化 
#############################################################################
install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR")
devtools::load_all("~/hh028/Tet2_DSS_scRNAseq/BM/pySCENIC/20231212/SCopeLoomR-master")
setwd("/home/shpc_100828/Pyscenic/")
#加载分析包
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
#可视化相关包，多加载点没毛病
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
setwd("~/hh028/Tet2_DSS_scRNAseq/BM/240204/pyscenic")
##读取pyscenic第三步分析的文件sce_SCENIC.loom
sce_SCENIC <- open_loom("sce_SCENIC.loom")
exprMat <- get_dgem(sce_SCENIC)#从sce_SCENIC文件提取表达矩阵
exprMat_log <- log2(exprMat+1) # log处理
write.csv(exprMat,file = "exprMat.csv")
write.csv(exprMat_log,file = "exprMat_log.csv")
#这里的表达矩阵其实就是我们在pyscenic分析第一步的输入矩阵，可见这些文件都是在一起的

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
#提取第二步分析的regulons,column.attr.name填Regulons，具体按照实际情况提示选择
#也就是我们所输入的基因和找到的转录因子组成的表达文件
regulons <- regulonsToGeneLists(regulons_incidMat)#将上一步矩阵文件转化为list
class(regulons)


#提取pyscenic第三步分析中AUC结果
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
write.csv(regulons, file = "regulonAUC.csv")
#以上就是一些主要文件了、够后续分析和可视化
#####################################################################################
##==============================加载seurat对象、RSS分析=======================================
#在可视化之前，我们再做一个分析，计算RSS值，计算regulon特异性评分
#BM_data <- readRDS("~/Pyscenic/test.merge.rds")
test.merge<- HSPC_sub3_harmony
cellinfo <- test.merge@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
Idents(test.merge) <- "cell_type"
GMP <- subset(scRNA_,idents = "GMP")
TD <- subset(test.merge, orig.ident == "Tet2mut_DSS")
TV <- subset(test.merge, orig.ident == "Tet2mut_Veh")
WD <- subset(test.merge, orig.ident == "WT_DSS")
WV <- subset(test.merge, orig.ident == "WT_Veh")
cellinfo <- test.merge@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
HSC2 <-subset(test.merge,idents = "HSC2")
cellinfo <- HSC1@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
cellinfo <- HSC2@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
cellinfo <- TD@meta.data[,c('cell_type','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息

colnames(cellinfo)=c('celltype', 'sample','nGene' ,'nUMI')
######计算细胞特异性TF
#在实际数据分析应用中，我认为比较靠谱的应用在于，细胞分了亚群，例如macrophage，有不同的特征
#我们可以查看不同亚群特异性的TF，有助于了解亚群的功能！！！！
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
sample <- as.data.frame(subset(cellinfo,select = 'sample'))
selectedResolution <- "celltype"
selectedResolution <- "sample"
sub_regulonAUC <- regulonAUC
#cellsPerGroup_TD <-cellsPerGroup$Tet2mut_DSS
#sub_regulonAUC <- t(sub_regulonAUC)
#sub_regulonAUC_TD = sub_regulonAUC[colnames(sub_regulonAUC) %in% cellsPerGroup_TD,]
#sub_regulonAUC_TD <-subset(sub_regulonAUC,cell_type)
rownames(sub_regulonAUC)
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),#从aucellresults获取AUC矩阵
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)#去除含有NA的行
#可视化细胞特异性TF
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 1.5,
    cluster_columns = FALSE,
    order_rows = TRUE,
    trh=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')

rssPlot

#提取数据，可以自己可视化dotplot，或者热图
rss_data <- rssPlot$df
devtools::install_github("XiaoLuo-boy/ggheatmap")
library(ggheatmap)
library(ggpubr)
library(reshape2)
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
col_ann <- data.frame(group= c(rep("HSC1",1),
                               rep("HSC2",1),
                               rep("MEP",1),
                               rep("NMP",1),
                               rep("Pro_B",1),
                               rep("MDP",1),
                               rep("Pro_DC",1),
                               rep("ErP",1),
                               rep("Pro_Mast",1),
                               rep("CLP",1),
                               rep("Pro_Mono",1),
                               rep("Pro_NE",1),
                               rep("Pro_MK",1)))#列注释
rownames(col_ann) <- colnames(rss_data)
mycolor<- c("#FED439FF", "#709AE1FF", "#8A9197FF", "blue",
            "#FD7446FF",  "#197EC0FF", "#D2AF81FF", 
            "#71D0F5FF", "#370335FF", "#075149FF",
            "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
names(x = mycolor) <- c("CLP", "ErP", "HSC1", 'HSC2', 'MDP', 'MEP', 'NMP', 'Pro_B', 'Pro_DC', 'Pro_Mast', 'Pro_Mk', 'Pro_Mono','Pro_NE')
col <- list(group=mycolor)

text_columns <- sample(colnames(rss_data),0)#不显示列名
typeof(rss_data)
as.matrix(rss)
p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = F,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p

rownames(rss) <- rss[,1]
rss <- rss[,-1]
colnames(rss)
col_ann <- data.frame(group= c(rep("HSC",1),
                               rep("CMP",1),
                               rep("GMP",1),
                               rep("MEP",1),
                               rep("MP",1),
                               rep("MDP",1),
                               rep("MkP",1),
                               rep("NP",1),
                               rep("CLP",1),
                               rep("ErP",1),
                               rep("Pro_B",1)))#列注释
rownames(col_ann) <- colnames(rss)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC","#466983FF",
              "#5050FFFF","#CE3D32FF","#F0E685FF","#5DB1DDFF","#BA6338FF")
names(groupcol) <- c("MEP","GMP" ,  "CMP" ,  "Pro_B", "MDP"  , "CLP"  , "ErP",   "MP" ,   "HSC",   "MkP",   "NP")
col <- list(group=groupcol)

text_columns <- sample(colnames(rss),0)#不显示列名

p <- ggheatmap(rss,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = F,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p





ggsave("a.pdf",p,limitsize = F)
######既然可以分析细胞中的特异性TF，那么也可以分析样本中特异性的TF
#我们可以提取某一个细胞的表达矩阵，做pyscenic，后期再进行RSS分析的时候
#可以使用样本，就可以看出哪些TF是这个样本细胞中特异性的TF，这样生物学意义也就有了
##############################################################################

##==============================TF_AUC与seurat结合===========================
#普通展示
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)
next_regulonAUC <- regulonAUC[,match(colnames(test.merge),colnames(regulonAUC))]
dim(next_regulonAUC)#将AUC结果于seurat对象结合

regulon_AUC <- regulonAUC@NAMES
test.merge@meta.data = cbind(test.merge@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))
HSC1 <- subset(test.merge,idents = "HSC1")
View(test.merge@meta.data)
Idents(HSC1) <- "orig.ident"
#自己选定感兴趣的或者比较重要的转录因子，这里我是随机的
TF_plot <- c("Egr1(+)","Egr2(+)","Egr3(+)","Hoxb2(+)",
             "Thra(+)","Myb(+)","Foxd2(+)","Irf6(-)","Gata3(-)",
             "Hoxa10(+)","Meis1(+)","Erg(+)","Rarb(+)","Irf6(+)","Irf8(+)")

TF_plot <- c("Rarb(+)","Rarb(-)","Nfkb2(+)","Nfkb2(-)","Klf10(+)","Klf10(-)","Pbx2(+)","Pbx2(-)",
             "Junb(+)","Junb(-)","Sox4(+)","Sox4(-)","Irf1(+)","Irf1(-)","Egr1(+)","Egr1(-)",
             "Sp4(+)","Sp4(-)","Ikzf1(+)","Ikzf1(-)","Hoxc9(+)","Hoxc9(-)","Stat1(+)","Stat1(-)",
             "Nfatc2(+)","Nfatc2(-)","Rela(+)","Rela(-)","Relb(+)","Relb(-)",
             "Fosb(+)","Fosb(-)","Ikzf2(+)","Ikzf2(-)","Gli3(+)","Gli3(-)")
DotPlot(TV, features = TF_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

####

DotPlot(HSC1, features = TF_plot, group.by = 'orig.ident')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  theme(legend.direction = "horizontal", 
        legend.position = "bottom")+
  labs(x=NULL,y=NULL)

###
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(Seurat)
pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)
FeaturePlot(test.merge, features =c("Egr3(+)","Erg(+)","Cebpe(+)","Gata2(+)","Cebpa(+)"),split.by = "orig.ident",cols = pal, order = T)
FeaturePlot(test.merge, features =c("Egr3","Erg","Cebpe","Gata2","Cebpa"),split.by = "orig.ident",cols = pal, order = T)

#############################################################################


cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

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


################################################------------------------
#------rank可视化rss-----------------------------------------------------

B_rss <- as.data.frame(rss_WD)#rss特异性TF结果
#需要作图的细胞类型
celltype <- c("WD_HSC1","WD_HSC2","WD_MEP","WD_Pro_Mast","WD_Pro_NE")
rssRanklist <- list()
library(ggrepel)
for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), color="black", size=3, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', 
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[i]] <- p
}


library(cowplot)

plot_grid(rssRanklist[[1]],rssRanklist[[2]],rssRanklist[[3]],
          rssRanklist[[4]],rssRanklist[[5]])



#=======================================================================================
#                              pySCENIC的差异分析及其他思路
#=======================================================================================
#加载R包
library(limma)
library(SCENIC)
library(AUCell)
library(data.table)


anaAUC <- regulonAUC
anaAUC <- anaAUC[onlyNonDuplicatedExtended(rownames(anaAUC)),]
ana.cellinfo<-test.merge@meta.data#细胞信息

#选取需要分析的细胞，这里我们以T cell为例子
Pro_NEcell.cellinfo <- subset(ana.cellinfo, ana.cellinfo$cell_type=='Pro_NE')
targets<-data.table(FileName=rownames(Pro_NEcell.cellinfo),Target=Pro_NEcell.cellinfo$orig.ident)#提取组合分组

#接下来其实就是常规的差异分析了
lev<-unique(targets$Target)
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #样本矩阵
colnames(design) <- lev #design矩阵重命名

eset=getAUC(anaAUC) #获取所有细胞TF的AUC值
eset=eset[,targets$FileName]#我们只选取要分析的T细胞
eset<-t(scale(t(eset))) #scale

#比较分析limma
cont.wt <- makeContrasts("WT_DSS-WT_Veh",levels=design) #GM vs BM
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)#差异分析结果

tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")#选择关键的三列即可
write.table(tT,file="Tcell_TF_GM_BM.csv",sep="\t",quote=F)#保存结果

#接下来就可以按照不同的阈值筛选显著的TF了，其实这样的差异分析对于结果的解释和
#应用上更有说服力，这样我们可以直观的看到某些转录因子在不同细胞状态下的差异
#可视化的话可以用热图、火山图等展示。
setwd("~/hh028/Tet2_DSS_scRNAseq/BM/240204/pyscenic/diff_sig/WD_vs_WV")
logFoldChange<-0.5
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]#显著筛选
diffSig=na.omit(diffSig)
rownames(diffSig)<-gsub("\\(","",rownames(diffSig))
rownames(diffSig)<-gsub("\\)","",rownames(diffSig))
rownames(diffSig)<-gsub("\\+","",rownames(diffSig))
write.csv(diffSig, file = "Pro_NE_WD_WV_diffSig_TF.csv")
rownames(diffSig)


#我们用热图展示以下差异结果
#先计算每个样本细胞的平均表达值
regulonActivity_byCellType <- sapply(split(rownames(HSC2cell.cellinfo), HSC2cell.cellinfo$orig.ident),
                                     function(cells) rowMeans(getAUC(anaAUC)[,cells]))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType
#这里有一个基础知识，很多文章展示TF不显示+号，用gsub函数就可以去除
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\(","",rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\)","",rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\+","",rownames(regulonActivity_byCellType_Scaled))

#选择差异显著的TF进行可视化
####TD vs WD HSC1
selTF <- c("Tgif1(+)","Junb(+)","Egr1(-)","Ltf(-)","Jun(+)","Gfi1(-)","Rara(-)",
           "Fosl1(-)","Fosb(+)","Egr2(+)","E2f1(+)", "Egr1(+)","Gata1(-)","Irf9(+)",
           "Meis1(-)","Myb(-)","Ets1(+)","Gata2(+)","Foxp1(+)","Fos(+)","Jund(+)",
           "Pbx2(+)","Cebpa(+)","Irf1(-)","Irf6(+)","Irf8(+)","Sox4(+)","Foxo1(+)",
           "Irf4(-)","Hoxa2(+)","Crem(+)","Pou3f1(+)","Sox4(-)")

####TV vs WV HSC1 
selTF <- c("Egr2(+)","Ikzf2(+)","Junb(+)","Rnf114(-)","Tcf7l2(+)","Rela(+)","Klf10(+)","Irf8(-)","Fos(+)","Foxd2(+)","Stat1(+)","Hoxb8(-)","Nfkb2(+)",
           "Meis1(+)","Hes1(-)","Rarb(+)","Irf1(+)","Jun(+)","Pou6f1(+)","Cebpb(-)","Stat3(+)","Irf1(-)","Erf(-)","Cebpd(-)","Elf2(+)","Ikzf3(-)",
           "Nr3c1(+)","Rfxank(+)","Irf7(+)","Sox4(-)","Maz(-)","Fli1(-)","Zscan22(+)")


####WD vs WV HSC1
selTF <- c("Atf4(+)","Fosl2(+)","Myc(+)","Hes1(+)","Myb(-)","Myb(+)","Fos(+)","Atf3(+)","Irf1(+)","Tal1(-)","Irf3(+)","Fosl2(-)","Xbp1(+)","Ikzf2(+)",
           "Irf8(-)","Mafk(+)","Ikzf3(-)","Cux1(-)","Nfkb2(+)","Gata2(-)","Fli1(-)","Stat3(+)","Pou6f1(+)","Rnf114(-)","Cebpa(+)","Irf8(+)","Foxo1(+)","Pou2f1(+)",
           "Bclaf1(+)","Bptf(-)","Etv5(+)","Etv1(+)","Egr3(+)","Erf(+)","Hmbox1(-)","Tcf7l2(+)","Hoxb2(-)","Runx2(+)","Cebpb(+)","Zfhx3(-)","Eomes(-)","Erg(+)",
           "Cebpd(+)","Thra(-)","Egr2(+)","Gata1(-)","Nfic(-)")

####TD vs TV HSC1
selTF <- c("Atf4(+)","Hes1(+)","Fosl2(+)","Myc(+)","Irf1(+)","Nfkb2(+)","Gata1(-)","Mafk(+)","Rarb(+)","Egr3(+)","Myb(-)","Fosl2(-)","Tal1(-)","Klf10(+)",
           "Fos(+)","Cebpa(+)","Pou2f1(+)","Xbp1(+)","Elk3(+)","Gata2(-)","Myb(+)","Irf3(+)","Rela(+)","Stat1(+)","Atf3(+)","Bclaf1(+)","Tgif1(+)","Cebpd(+)",
           "Foxp1(+)","Ikzf1(-)","Elk1(+)","Egr1(+)","Foxo1(-)","Creb5(-)","Fiz1(-)","Irf1(-)","Hoxa3(-)","Cebpb(+)","Srebf1(+)","Fosb(+)","Pbx1(-)","Yy1(+)",
           "Rfx1(-)")

sel_regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[selTF,]
#'#B53E2B','#F3B1A0',"black","#9FA3A8"
annotation_col = data.frame(group = c("WT_Veh","Tet2mut_Veh","WT_DSS","Tet2mut_DSS"))
rownames(annotation_col)<-colnames(regulonActivity_byCellType_Scaled)
ann_colors =  list(group = c("WT_Veh" = '#B53E2B', "Tet2mut_Veh" = '#F3B1A0',"WT_DSS" = 'black', "Tet2mut_DSS" = '#9FA3A8'))

pheatmap::pheatmap(sel_regulonActivity_byCellType_Scaled,
                   color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100),
                   breaks=seq(-2, 2, length.out = 100),
                   border_color=NA,
                   cluster_rows = F,
                   cluster_cols=F,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   scale="row",
                   main = "T cell TF")

#=======================================================================================
#                                   pySCENIC 转录因子与靶基因
#=======================================================================================
#我们还可以可视化感兴趣的TF调控的靶基因
#TF与靶基因的关系在pyscenic分析得到的第二个文件

sce_regulons <- read.csv("sce.regulons.csv")
sce_regulons <- sce_regulons[-2, ]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ]
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

#举例子我这里关注FOXP3和TGIF1这两个TF

FOXP3 <- subset(sce_regulons, TF=='FOXP3')
FOXP3 <- FOXP3[which(FOXP3$AUC>0.1),]
FOXP3 <- FOXP3[, c("TF","TargetGenes")]
FOXP3$TargetGenes <-gsub("\\[","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\]","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\(","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\)","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\'","",FOXP3$TargetGenes)
library(stringr)
split_FOXP3<-str_split(FOXP3$TargetGenes,",")
FOXP31 <- as.data.frame(split_FOXP3[[1]])
FOXP32<- as.data.frame(split_FOXP3[[2]])
FOXP33<-as.data.frame(split_FOXP3[[3]])
FOXP32<-as.data.frame(split_FOXP3[[4]])
FOXP35<- as.data.frame(split_FOXP3[[5]])

names(FOXP31) <- 'TF'
names(FOXP32) <- 'TF'
names(FOXP33) <- 'TF'
names(FOXP34) <- 'TF'
names(FOXP35) <- 'TF'

FOXP3 <- rbind(FOXP31,FOXP32,FOXP33,FOXP34,FOXP35)

FOXP3_target <- FOXP3[seq(1,nrow(FOXP3),2), ]
FOXP3_score <- FOXP3[seq(0,nrow(FOXP3),2), ]

FOXP3_gene <- data.frame(FOXP3_target,FOXP3_score)
FOXP3_gene <- FOXP3_gene[!duplicated(FOXP3_gene$FOXP3_target), ]
FOXP3_gene$gene <- 'FOXP3'
colnames(FOXP3_gene) <- c("target","score",'tf')

#同理得到TGIF1及其靶基因，此处省略1万字
#two thousand years later
#TGIF1_gene

TF_target <- rbind(FOXP3_gene,TGIF1_gene)
TF_target$score <- as.numeric(TF_target$score)
###接下来就是网络图了

#节点数据
paths <- c("FOXP3", "TGIF1")#列重命名
nodelist <- list()
for (i in 1:length(paths)){
  node <- subset(TF_target, tf == paths[i])#提取数据
  
  nodes <- data.frame(name = unique(union(node$tf, node$target)))#整理为datafram
  nodes$value <- c(sum(node$score)/10, node$score)#加上values
  
  nodelist[[i]] <- nodes
}  #提取每个大节点数据


nodes <- rbind(nodelist[[1]],nodelist[[2]])#将三个节点文件合并
nodes$cluster <- c(rep("FOXP3",1),rep("FOXP3_gene",20),
                   rep("TGIF1",1),rep("TGIF1_gene",6))#分组，为了后续颜色设置

edges <- TF_target[c("tf","target","score")]#边缘文件
edges$class <- edges$tf

library(ggraph)
library(tidygraph)
layout_cir <- tbl_graph(nodes = nodes, edges = edges)#构建ggraph作图文件
#作图
ggraph(layout_cir,layout='linear',circular = TRUE) +#选择circle
  geom_node_point(aes(size=value,colour = cluster))+#节点，大小用我们赋的值表示，颜色用分组
  geom_node_text(aes(x = 1.03 * x,
                     y = 1.03 * y,
                     label=name,
                     color=cluster,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                 hjust='outward') +#文字设置。x，y是为了调整位置。angle是为了调整角度，以后其他所有网络图angle都用此公式，文字会向外发散排列
  geom_edge_arc(aes(colour=class))+#连线为曲线
  theme_void()+#theme主题
  theme(legend.position = "none")+
  scale_colour_manual(values =c('#407972',
                                '#961E28',
                                '#D46724',
                                '#0f8096'))+#节点颜色
  scale_edge_colour_manual(values = c('#961E28',
                                      '#D46724',
                                      '#0f8096'))+#连线颜色
  scale_size_continuous(range = c(2,8))+#点的大小范围设置
  coord_cartesian(xlim=c(-1.5,1.5),ylim = c(-1.5,1.5))#设置坐标位置，防止图溢出作图边界显示不全




