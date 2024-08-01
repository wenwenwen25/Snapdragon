d1= read.csv("~/TOSICA/pathway/IFN/TOSICA_PreMatrix_TET2_IFN.csv")
d2= read.csv("~/TOSICA/pathway/TNF/TOSICA_PreMatrix_TET2_TNF.csv")
d3= read.csv("~/TOSICA/pathway/ADRB/TOSICA_PreMatrix_TET2_ADRB.csv")
d4= read.csv("~/TOSICA/pathway/IL1R/TOSICA_PreMatrix_TET2_IL1R.csv")
d5= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMatrix_WT.csv")
d6= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMatrix_TET2.csv")
d7= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMatrix_WTDSS.csv")
d8= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMatrix_TET2DSS.csv")
d9= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/ADRB_IL1R/TOSICA_PreMatrix_TET2_ADRB_IL1R.csv")
WT <- subset(test.down,orig.ident == "WT_DSS")
TD <- subset(test.down,orig.ident == "Tet2mut_DSS")

WT_d <- WT@meta.data
TD_d <- TD@meta.data

TD_d$barcode.1 <- rownames(TD_d)
d1 <- d1[,-1]
d2 <- d2[,-1]
d3 <- d3[,-1]
d4 <- d4[,-1]
d5 <- d5[,-1]
d6 <- d6[,-1]
d7 <- d7[,-1]
d8 <- d8[,-1]
d9 <- d9[,-1]
colnames(d9)<- c("CLP",
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

df1= read.csv("~/TOSICA/pathway/IFN/TOSICA_PreMeta_TET2_IFN.csv")
df2= read.csv("~/TOSICA/pathway/TNF/TOSICA_PreMeta_TET2_TNF.csv")
df3= read.csv("~/TOSICA/pathway/ADRB/TOSICA_PreMeta_TET2_ADRB.csv")
df4= read.csv("~/TOSICA/pathway/IL1R/TOSICA_PreMeta_TET2_IL1R.csv")
df5= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMeta_WT.csv")
df6= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMeta_TET2.csv")
df7= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMeta_WTDSS.csv")
df8= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/TOSICA_PreMeta_TET2DSS.csv")
df9= read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/ADRB_IL1R/TOSICA_PreMeta_TET2_ADRB_IL1R.csv")
a1 <- cbind(d1,df1)
a2 <- cbind(d2,df2)
a3 <- cbind(d3,df3)
a4 <- cbind(d4,df4)

a5 <- cbind(d5,df5)
a6 <- cbind(d6,df6)
a7 <- cbind(d7,df7)
a8 <- cbind(d8,df8)
a9 <- cbind(d9,df9)

modified_string <- sub("(.*)_4", "Tet_DSS:\\1", original_string)
metaDF_new<- a9
metaDF_new$total<- metaDF_new$HSC2+metaDF_new$Pro_NE
metaDF_new$bias<- metaDF_new$HSC2-metaDF_new$Pro_NE
metaDF_new[1:4,]

metaDF_new_final <- metaDF_new
bias_df<- metaDF_new_final
bias_df<- bias_df[bias_df$total>0.5,]
bias_df[1:4,]
dim(bias_df)
bias_df$bias_n<- bias_df$bias/bias_df$total

# 找到 "bias" 列的最小值和最大值
min_bias <- min(bias_df$bias_n)
max_bias <- max(bias_df$bias_n)

# 归一化 "bias" 列
normalized_bias <- (bias_df$bias - min_bias) / (max_bias - min_bias)

# 将归一化后的值替换到原始数据框中
bias_df$bias_Z <- normalized_bias


bias_df[1:4,]

# 计算每0.1范围内的百分比
calculate_percentage <- function(data, lower_bound, upper_bound) {
  subset_data <- subset(data, bias_Z >= lower_bound & bias_Z < upper_bound)
  percentage <- nrow(subset_data) / nrow(data) * 100
  return(percentage)
}
percentage <- calculate_percentage(bias_df, 0.4, 0.6)
print(paste("Percentage in the range 0.4-0.6:", round(percentage, 2), "%"))

library(ggplot2)
library(RColorBrewer)
ggplot()+
  geom_point(data=metaDF_new_final, aes(x=UMAP_1, y=UMAP_2),size=1, alpha=0.7, shape=20, color="grey")+
  geom_point(data=bias_df, aes(x=UMAP_1, y=UMAP_2, color=bias_Z),size=1, alpha=0.7, shape=20)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")))+
  theme_classic()+
  ggtitle("HSC vs. MEP bias in Tet2mutDSS")



# 定义数据框
data <- data.frame(
  Sample = c("WT_Veh", "Tetmut_Veh", "Tetmut_Veh_IFN_KO", "Tetmut_Veh_TNF_KO", "Tetmut_Veh_ADRB_KO", "Tetmut_Veh_IL1R_KO","Tetmut_Veh_IL1R_ADRB_KO"),
  HSC_vs_MEP = c(3.22, 7.25, 6.62, 7.38, 6.61, 7.26,7.26),
  GMP_vs_Pro_NE = c(8.79, 11.89, 11.28, 10.24, 11.32, 11.61, 11.6),
  GMP_vs_MEP = c(6.19, 4.47, 5.52, 7.98, 5.56, 5.73, 5.81)
)

# 打印数据框
print(data)

# 设置图形布局为 2x2
par(mfrow = c(1, 3))

# 定义颜色
colors <- c("cornflowerblue", "seagreen3", "salmon", "gold", "mediumorchid", "tomato")

# 绘制 HSC vs MEP 的条形图
barplot(data$HSC_vs_MEP, 
        names.arg = data$Sample, 
        main = "HSC vs MEP", 
        xlab = "orig.ident", 
        ylab = "Precentage", 
        col = colors, 
        las = 2,  # x轴标签旋转90度
        ylim = c(0, max(data$HSC_vs_MEP) + 1))

# 绘制 GMP vs Pro_NE 的条形图
barplot(data$GMP_vs_Pro_NE, 
        names.arg = data$Sample, 
        main = "GMP vs Pro_NE", 
        xlab = "orig.ident", 
        ylab = "Precentage", 
        col = colors, 
        las = 2,  # x轴标签旋转90度
        ylim = c(0, max(data$GMP_vs_Pro_NE) + 1))

# 绘制 GMP vs MEP 的条形图
barplot(data$GMP_vs_MEP, 
        names.arg = data$Sample, 
        main = "GMP vs MEP", 
        xlab = "orig.ident", 
        ylab = "Precentage", 
        col = colors, 
        las = 2,  # x轴标签旋转90度
        ylim = c(0, max(data$GMP_vs_MEP) + 1))

# 恢复默认单图布局
par(mfrow = c(1, 1))
