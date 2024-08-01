pesudotime <- rbind(WT_Veh_pseudotime,Tet2_Veh_pseudotime,WT_DSS_pseudotime,Tet2_DSS_pseudotime)
colnames(pesudotime) <- "pesudotime"
select <- subset(test.merge, idents = c("HSC1","HSC2","MEP","Pro_Mast","Pro_NE"))
Idents(select) <- "orig.ident"
Tet2mut_DSS <-subset(select,idents = "Tet2mut_DSS")
Tet2mut_Veh <- subset(select,idents = "Tet2mut_Veh")
WT_DSS <-subset(select,idents = "WT_DSS")
WT_Veh <- subset(select,idents = "WT_Veh")
View(Tet2mut_DSS@meta.data)
View(Tet2mut_Veh@meta.data)
View(WT_DSS@meta.data)
exp <- select@assays$RNA@counts
exp <- as.matrix(exp)
exp <- t(exp)
library(tidyverse)
col_vec <- c("Gse1","Zfp36l2","Marcksl1","Mdga1","Nr4a1","Dusp6","Egr3",
             "Jund","Fos","Fosb","Ndrg1","Socs3","Sox4","Vezf1",
             "Osbpl1a","Egr2","Gata2","Junb","Ikzf2")

exp_s <- exp[,col_vec]
ident <- select@meta.data[,c("orig.ident","cell_type")]
data <- cbind(ident,exp_s)
data$barcode <- rownames(data)
pesudotime$barcode <- rownames(pesudotime)
data <-merge(data,pesudotime,,by = "barcode",all=T) 


#######TF
data <- cbind(select@meta.data,pesudotime)
#library(plyr)
#data <- bind.fill(data,pesudotime)

NE = data %>% filter(cell_type == c("HSC1","HSC2","Pro_NE"))
Mast = data %>% filter(cell_type == c("HSC1","MEP","Pro_Mast"))
Mast <- subset(data, cell_type %in% c("HSC1",  "MEP", "Pro_Mast"))
NE <- subset(data, cell_type %in% c("HSC1", "HSC2", "Pro_NE"))

head(re2)

#选择需要的列即可，我这里的origin.ident就是分组
Mast_s <-Mast %>% select("orig.ident","pesudotime",
                      "Atf1","Atf4", "Cd34", "Cebpb","Cebpd","Cebpe","Churc1","Egr1",
                      "Fos","Foxo4","Gata2","Hes1","Hoxa10","Hoxb5","Ikzf2","Il1b",
                      "Il1r1","Irf7","Junb","Meis1","Myb","Runx2","Rxrb",
                      "Smad3","Sp4","Srebf2","Srf","Stat1","Tfap4","Zbtb33","Ifitm1",
                      "Ifitm3","Ly6a")

NE_s <- NE %>% select("orig.ident","pesudotime","Nfe2l2(+)","Meis1(+)","E2f2(+)","Egr1(+)","Egr1(-)","Egr2(+)","Erg(+)","Cebpa(+)","Cebpe(+)","Myb(+)","Gata2(-)","Cebpb(-)",
                    "Pbx1(+)","Ikzf2(+)","Gfi1(+)","Cebpd(+)","Hoxa10(+)","Hoxa10(-)","Gata3(-)","Klf4(-)","Runx1(-)","Spi1(-)")
Mast_s <- Mast %>% select("orig.ident","pesudotime","Nfe2l2(+)","Meis1(+)","E2f2(+)","Egr1(+)","Egr1(-)","Egr2(+)","Erg(+)","Cebpa(+)","Cebpe(+)","Myb(+)","Gata2(-)","Cebpb(-)",
                      "Pbx1(+)","Ikzf2(+)","Gfi1(+)","Cebpd(+)","Hoxa10(+)","Hoxa10(-)","Gata3(-)","Klf4(-)","Runx1(-)","Spi1(-)")

NE_s <-NE %>% select("orig.ident","pesudotime",
                         "Gse1","Zfp36l2","Marcksl1","Mdga1","Nr4a1","Dusp6","Egr3",
                         "Jund","Fos","Fosb","Ndrg1","Socs3","Sox4","Vezf1",
                         "Osbpl1a","Egr2","Gata2","Junb","Ikzf2")


#使用分屏，应该就是文献种的办法
#首先将data宽数据转化为长数据
library(reshape)
NE_long_m<-melt(NE_s, id.vars = c("orig.ident", "pesudotime"), #需保留的不参与聚合的变量列名
                  measure.vars = 3:21,#选择需要转化的列
                  variable.name = c('gene'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名
colnames(NE_long_m) <- c("orig.ident","pesudotime","gene","value")

p=ggplot(NE_long_m, aes(x=pesudotime, y=value, color=orig.ident))+
  geom_smooth(aes(fill=orig.ident))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative activity') +
  facet_wrap(~gene, scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c('#B53E2B','#F3B1A0',"black","#9FA3A8"))+#修改颜色
  scale_fill_manual(name=NULL, values = c('#B53E2B','#F3B1A0',"black","#9FA3A8"))#修改颜色
p
p + theme(text = element_text(face = "bold", color="black"), 
          plot.title = element_text(size = 16), 
          plot.subtitle = element_text(size = 10), 
          #legend.position = "none",
          axis.text = element_text(size = 10, color="black"),
          axis.title = element_text(size = 10),
          axis.line = element_line(colour = "black", size = 1),
          axis.ticks.length=unit(.35, "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
library(fgsea)
library(nichenetr)

pathway_file <- "/disk2/cai039/Tet2_DSS/20231127/monocle/m5.go.v2023.2.Mm.symbols.gmt"
Egr1_target <- c("Gse1","Zfp36l2","Marcksl1","Gata2","Mdga1","Egr1",
                 "Nr4a1","Dusp6","Egr3","Jund","Fosb","Ndrg1","Ikzf2",
                 "Socs3","Sox4","Fos","Junb","Vezf1","Osbpl1a","Egr2")
Egr1_target <- as.data.frame(Egr1_target)
Egr1_target_gene <- Egr1_target$Egr1_target
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
View(pathway_names)

#选择需要的通路查看基因（自己关注的或者感兴趣的通路）
#这里我随便找两个举例子
IFNA <- pathways[["BIOCARTA_IFNA_PATHWAY"]]
IFNG <- pathways[["BIOCARTA_IFNG_PATHWAY"]]
IL1R <- pathways[["BIOCARTA_IL1R_PATHWAY"]]
Inflamation <- pathways[["BIOCARTA_INFLAM_PATHWAY"]]
Neu <- pathways[["BIOCARTA_NEUTROPHIL_PATHWAY"]]
TNFR2 <- pathways[["BIOCARTA_TNFR2_PATHWAY"]]
TNFR1 <- pathways[["BIOCARTA_TNFR1_PATHWAY"]]
TGFB <- pathways[["BIOCARTA_TGFB_PATHWAY"]]


#我用的是鼠的示例数据，所以需要将人的代谢基因转化为鼠的。如果是人的
#就不必操作了
IFNA <- as.data.frame(IFNA)
#IL2$gene = IL2$IL2 %>% convert_human_to_mouse_symbols()
#去除NA
IFNA <- na.omit(IFNA)
IFNA_gene <- IFNA$IFNA


IFNG <- as.data.frame(IFNG)
#IL6$gene = IL6$IL6 %>% convert_human_to_mouse_symbols()
#去除NA
IFNG <- na.omit(IFNG)
IFNG_gene <- IFNG$IFNG

IL1R <- as.data.frame(IL1R)
#IL2$gene = IL2$IL2 %>% convert_human_to_mouse_symbols()
#去除NA
IL1R <- na.omit(IL1R)
IL1R_gene <- IL1R$IL1R

Inflamation <- as.data.frame(Inflamation)
#IL6$gene = IL6$IL6 %>% convert_human_to_mouse_symbols()
#去除NA
Inflamation <- na.omit(Inflamation)
Inflamation_gene <- Inflamation$Inflamation


Neu <- as.data.frame(Neu)
#IL2$gene = IL2$IL2 %>% convert_human_to_mouse_symbols()
#去除NA
Neu <- na.omit(Neu)
Neu_gene <- Neu$Neu

TNFR1 <- as.data.frame(TNFR1)
#IL2$gene = IL2$IL2 %>% convert_human_to_mouse_symbols()
#去除NA
TNFR1 <- na.omit(TNFR1)
TNFR1_gene <- TNFR1$TNFR1

TNFR2 <- as.data.frame(TNFR2)
#IL2$gene = IL2$IL2 %>% convert_human_to_mouse_symbols()
#去除NA
TNFR2 <- na.omit(TNFR2)
TNFR2_gene <- TNFR2$TNFR2

TGFB <- as.data.frame(TGFB)
#IL2$gene = IL2$IL2 %>% convert_human_to_mouse_symbols()
#去除NA
TGFB <- na.omit(TGFB)
TGFB_gene <- TGFB$TGFB

library(Seurat)
#计算代谢通路评分，添加到seurat对象
DefaultAssay(select) <- "SCT"




select <- AddModuleScore(select,
                       features=list('IFNA' = IFNA_gene,
                                     'IFNG'=IFNG_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('IFNA', 'IFNG'))

select <- AddModuleScore(select,
                         features=list('Egr1_target' = Egr1_target_gene),
                         pool = rownames(select), k=F, nbin=24,
                         name = c('Egr1_target'))


select <- AddModuleScore(select,
                       features=list('IFNA' = IFNA_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('IFNA'))

select <- AddModuleScore(select,
                       features=list('IFNG' = IFNG_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('IFNG'))

select <- AddModuleScore(select,
                       features=list('IL1R' = IL1R_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('IL1R'))

select <- AddModuleScore(select,
                       features=list('Inflamation' = Inflamation_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('Inflamation'))

select <- AddModuleScore(select,
                       features=list('Neu' = Neu_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('Neu'))

select <- AddModuleScore(select,
                       features=list('TNFR2' = TNFR2_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('TNFR2'))

select <- AddModuleScore(select,
                       features=list('TNFR1' = TNFR1_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('TNFR1'))

select <- AddModuleScore(select,
                       features=list('TGFB' = TGFB_gene),
                       pool = rownames(select), k=F, nbin=24,
                       name = c('TGFB'))
data1 <- select@meta.data
data1 <-cbind(data1,pesudotime) 
Mast <- subset(data1, cell_type %in% c("HSC1",  "MEP", "Pro_Mast"))
NE <- subset(data1, cell_type %in% c("HSC1", "HSC2", "Pro_NE"))

#提取作图数据

colnames(NE)
data1 <- data1[, c("orig.ident", "Pseudotime", "IFNA1")]
data1 <- data1[, c("orig.ident", "pseudotime", "IFNG2")]
data1 <- data1[, c("orig.ident", "pseudotime", "IL1R1")]
#data1 <- data1[, c("orig.ident", "pseudotime", "IL21")]
data1 <- data1[, c("orig.ident", "pseudotime", "TGFB1")]
NE1 <- NE[, c("orig.ident", "pesudotime", "Egr1_target1")]
Mast1 <- Mast[, c("orig.ident", "pesudotime", "Egr1_target1")]
data1_long   <-melt(NE1, id.vars = c("orig.ident", "pseudotime"), #需保留的不参与聚合的变量列名
                    measure.vars = 3:6,#选择需要转化的列
                    variable.name = c('pathway'),#聚合变量的新列名
                    value.name = 'value')#聚合值的新列名

data1_long   <-melt(Mast1, id.vars = c("orig.ident", "pesudotime"), #需保留的不参与聚合的变量列名
                    measure.vars = 3,#选择需要转化的列
                    variable.name = c('pathway'),#聚合变量的新列名
                    value.name = 'value')#聚合值的新列名


data1_new <- data1_long
colnames(data1_new) <- c("orig.ident","pesudotime","pathway","value" )
levels(data1_new$pathway) <- c('BIOCARTA_IFNA_PATHWAY',
                               'BIOCARTA_IFNG_PATHWAY',
                               'BIOCARTA_IL1R_PATHWAY',
                               'BIOCARTA_INFLAM_PATHWAY',
                               "BIOCARTA_NEUTROPHIL_PATHWAY",
                               "BIOCARTA_TNFR1_PATHWAY",
                               "BIOCARTA_TGFB_PATHWAY")

levels(data1_new$pathway) <- c('BIOCARTA_TGFB1_PATHWAY')

levels(data1_new$pathway) <- c('BIOCARTA_NEUTROPHIL_PATHWAY')

levels(data1_new$pathway) <- c('Egr1_target genes')


my_comparisons <- list(c('Tet2mut_DSS','Tet2mut_Veh'),c('Tet2mut_DSS','WT_DSS'),
                       c("WT_DSS","WT_Veh"),c("Tet2_Veh","WT_Veh"))
library(ggpubr)

library(ggsignif)
#作图和上述一样
p =ggplot(data1_new, aes(x=pesudotime, y=value, color=orig.ident))+
  geom_smooth(aes(fill=orig.ident))+ #平滑的填充
  xlab('pesudotime') + 
  ylab('Relative expression') +
  facet_wrap(~levels(pathway), scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 12),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c('#B53E2B','#F3B1A0',"black","#9FA3A8"))+ #修改颜色
  scale_fill_manual(name=NULL, values = c('#B53E2B','#F3B1A0',"black","#9FA3A8")) #修改颜色
p
p + theme(text = element_text(face = "bold", color="black"), 
          plot.title = element_text(size = 16), 
          plot.subtitle = element_text(size = 10), 
          #legend.position = "none",
          axis.text = element_text(size = 10, color="black"),
          axis.title = element_text(size = 10),
          axis.line = element_line(colour = "black", size = 1),
          axis.ticks.length=unit(.35, "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)

FeaturePlot(scRNA_,features = c("Ly6c2","Ms4a6c","Mpo","Ctsg") ,cols = pal, order = T,split.by = "orig.ident")


