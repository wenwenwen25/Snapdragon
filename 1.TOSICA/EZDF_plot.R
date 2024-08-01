setwd("~/TOSICA_GRN_GOBP")

df1= read.csv("TOSICA_PreMatrix_WT.csv")
df1[1:4,]
df1= df1[,-1]
colnames(df1)
colnames(df1)<- c("CLP",
                  "ErP",
                  "HSC1",
                  "HSC2",
                  "MDP",
                  "MEP",
                  "NMP",
                  "Pro_B",
                  "Pro_DC",
                  "Pro_Mast",
                  "Pro_Mk",
                  "Pro_Mono",
                  "Pro_NE" )
df1[1:4,]
dim(df1)



df2= read.csv("TOSICA_PreMeta_WT.csv")
df2[1:4,]
df2= df2[,-1]
df2[1:4,]
dim(df2)
table(df2$cell_type)


library(ggplot2)
df = cbind(df2,df1)

mycolor<- c("#FED439FF", "#709AE1FF", "#8A9197FF", "blue",
            "#FD7446FF",  "#197EC0FF", "#D2AF81FF", 
            "#71D0F5FF", "#370335FF", "#075149FF",
            "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
#show_col(mycolor)


ggplot(data=df,aes(x=HSC2, color=cell_type))+
  stat_ecdf()+
  scale_color_manual(values=mycolor)+
  xlim(0,1)+
  theme_classic() ##Good one


ggplot(data=df,aes(x=HSC1, color=cell_type))+
  stat_ecdf()+
  scale_color_manual(values=mycolor)+
  xlim(0,1)+
  theme_classic() ##Good one

library(RColorBrewer)  
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=HSC2))+
  geom_point(size=1, alpha=0.7, shape=20)+
  #scale_color_gradient(low="gray",high = "#990000")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme_classic()



#downsample to 4500
load("/disk1/cai001/HehangDSS_Capybara/HSPC_downsample_harmony_rename_final.RData")
test.merge = test.down
DimPlot(test.merge, reduction = "umap", group.by = "cell_type", split.by= "orig.ident", ncol = 2, cols =mycolor)
table(test.merge@meta.data$orig.ident)

test1.seu<- subset(test.merge, subset=orig.ident=="WT_Veh")
test2.seu<- subset(test.merge, subset=orig.ident=="Tet2mut_Veh")
test3.seu<- subset(test.merge, subset=orig.ident=="WT_DSS")
test4.seu<- subset(test.merge, subset=orig.ident=="Tet2mut_DSS")

metaDF<- test1.seu@meta.data
metaDF[1:4,]
metaDF$barcode.1<- row.names(metaDF)
metaDF$barcode.1<- substr(metaDF$barcode.1,1,24)
metaDF$barcode.1<- paste("WT_1:",metaDF$barcode.1,sep = '')
metaDF[1:4,]

df[1:4,]
metaDF_new<- merge(metaDF,df, by="barcode.1")
metaDF_new[1:4,]

library(RColorBrewer)  

pdf("./plots_WT/WT_HSC2_prob_background.pdf", width = 10, height = 10) 
ggplot(metaDF_new, aes(x=UMAP_1, y=UMAP_2, color=HSC2))+
  geom_point(size=1, alpha=0.7, shape=20)+
  #scale_color_gradient(low="gray",high = "#990000")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme_classic()
dev.off()

pdf("./plots_WT/WT_HSC2_ECDF.pdf", width = 10, height = 10) 
ggplot(data=df,aes(x=HSC2, color=cell_type))+
  stat_ecdf()+
  scale_color_manual(values=mycolor)+
  xlim(0,1)+
  theme_classic() ##Good one
dev.off()






