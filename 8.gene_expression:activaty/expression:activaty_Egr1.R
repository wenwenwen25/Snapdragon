DefaultAssay(test.merge)="RNA" 
#在化学/氧化应激条件下，细胞中NRF2可逃离KEAP1的抑制，
#通过NRF2上Neh1域的一个核转位信号进入胞核，
#启动抗氧化反应元件(ARE)调节的靶基因的转录。
#NRF2的靶基因有：醌氧化还原酶-1(NQO1)，血红素加氧酶(HO-1)和谷氨酸-半胱氨酸连接酶(GCL)等。
#以上NRF2的靶基因能够保护细胞应对各种有毒物质引起的应激。
Idents(test.merge)<- "cell_type"
FeaturePlot(test.merge,features = "Adra2a",split.by = "orig.ident")
features <-c("Adra2a")  
expr = as.matrix(test.merge@assays$SCT@data[features,])
colnames(expr) <- "Adra2a"
posi= match(rownames(expr),rownames(test.merge@meta.data))
ts = cbind(expr,test.merge@meta.data)
ts = cbind(expr,test.merge@meta.data[posi,c("orig.ident","cell_type")])
library(ggpubr)
# 这里是画了单个基因的图，把基因名替换掉即可绘制其他基因的图
## 在这里我们指定了2组比较："compact" 和 "midsize"； "minivan"和"suv"
comparisons <- list(c("WT_Veh", "WT_DSS"), c("WT_Veh", "Tet2mut_Veh"),
                    c("WT_DSS", "Tet2mut_DSS"),c("Tet2mut_Veh", "Tet2mut_DSS"))
my_comparisons <- list( c("WT_Veh", "WT_DSS"), c("WT_Veh", "Tet2mut_Veh"), 
                        c("WT_DSS","Tet2mut_DSS"),c("Tet2mut_Veh", "Tet2mut_DSS") )
gene = "Adra2a"
expr_single = ts[,gene]
# 筛选特定细胞类型的数据
specific_cells <- ts[ts$cell_type %in% c("Pro_NE"), ]
expr_single = specific_cells[,gene]

specific_cells$orig.ident<- factor(specific_cells$orig.ident,levels = c('WT_Veh',"Tet2mut_Veh",'WT_DSS',"Tet2mut_DSS"))


ggplot(ts,aes(x = orig.ident, y = expr_single,color=orig.ident))+
  geom_violin(alpha=0.5,width=0.8,size=0.5)+ 
  scale_color_manual(values = c("#9FA3A8",'#F3B1A0',"black",'#B53E2B'))+
  xlab("")+
  ylab("Expression Level")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),
        axis.text.x = element_text(angle=45,size=10,vjust = 1,hjust =1,color = "black"),
        axis.text.y = element_text(size =10))+
  labs(fill = gene) +
  ggtitle(gene) + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_signif(comparisons = list(
    c("WT_Veh", "WT_DSS"),c("WT_Veh", "Tet2mut_Veh"),c("WT_DSS", "Tet2mut_DSS"),c("Tet2mut_Veh", "Tet2mut_DSS")
  ),y_position = c(4.4, 4, 4,4.8), map_signif_level = F, textsize=3)



pt + geom_signif(comparisons = list(
  c("WT_Veh", "WT_DSS"), 
  c("WT_Veh", "Tet2mut_Veh"), 
  c("WT_DSS", "Tet2mut_DSS"), 
  c("Tet2mut_Veh", "Tet2mut_DSS")
  ), map_signif_level = TRUE, textsize=3)

pt <- ggplot(ts,aes(x = cell_type, y = expr_single,fill=orig.ident))+
  geom_violin(scale = "width",alpha=0.8,width=0.6,size=0.8)+ 
  scale_fill_manual(values = c("#9FA3A8","black",'#F3B1A0','#B53E2B'))+
  geom_signif(comparisons = comparisons,
              test = "wilcox.test",#指定使用的检验方法
              textsize = 4#指定标记中文字的大小
  )+
  xlab("")+
  ylab("Expression Level")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),
        axis.text.x = element_text(angle=45,size=10,vjust = 1,hjust =1,color = "black"),
        axis.text.y = element_text(size =10))+
  labs(fill = gene) +
  ggtitle(gene) + 
  theme(plot.title = element_text(hjust = 0.5))
pt

pt + geom_signif(comparisons = list(
  c("WT_Veh", "WT_DSS"), 
  c("WT_Veh", "Tet2mut_Veh"), 
  c("WT_DSS", "Tet2mut_DSS"), 
  c("Tet2mut_Veh", "Tet2mut_DSS")
), map_signif_level = TRUE, textsize=3)


posi= rownames(test.merge@meta.data)
ts = test.merge@meta.data[posi,c("orig.ident","cell_type","Egr1(+)")]


gene = "Egr1(+)"
expr_single = ts[,gene]
# 筛选特定细胞类型的数据
specific_cells <- ts[ts$cell_type %in% c("HSC2"), ]
expr_single = specific_cells[,gene]

expr_single$orig.ident<- factor(expr_single$orig.ident,levels = c('WT_Veh',"Tet2mut_Veh",'WT_DSS',"Tet2mut_DSS"))
pt <- ggplot(ts,aes(x = cell_type, y = expr_single,color=orig.ident))+
  geom_boxplot(alpha=0.8,width=0.6,size=0.8, outlier.shape = NA)+ 
  scale_color_manual(values = c("#9FA3A8","black",'#F3B1A0','#B53E2B'))+
  xlab("")+
  ylab("Expression Level")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),
        axis.text.x = element_text(angle=45,size=10,vjust = 1,hjust =1,color = "black"),
        axis.text.y = element_text(size =10))+
  labs(fill = gene) +
  ggtitle(gene) + 
  theme(plot.title = element_text(hjust = 0.5))
pt

pt + stat_compare_means(comparisons = my_comparisons)


p2 <- ggboxplot(ts,x="cell_type", y= "Egr1",color = "orig.ident",
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                #outlier.shape=NA, #不显示outlier
                legend = "right")+ 
  scale_fill_manual(values = c("#9FA3A8","black",'#F3B1A0','#B53E2B'))

p2 <- ggboxplot(ts,x="cell_type", y="Egr1",color = "orig.ident")
p2
p2 + stat_compare_means(aes(group = group))

