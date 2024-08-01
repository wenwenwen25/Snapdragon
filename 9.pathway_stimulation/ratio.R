WT_TOSICA <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/TOSICA_PreMeta_WT.csv", row.names=1)
Tet2_TOSICA <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/TOSICA_PreMeta_TET2.csv", row.names=1)
WTDSS_TOSICA <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/TOSICA_PreMeta_WTDSS.csv", row.names=1)
Tet2DSS_TOSICA <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/TOSICA_PreMeta_TET2DSS.csv", row.names=1)
Tet2_IFN <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/Tet2_IFN_TOSICA_prediction_meta_table.csv", row.names=1)
WTDSS_TNF <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/WTDSS_TNFA_TOSICA_prediction_meta_table.csv", row.names=1)
Tet2DSS_IL1R <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/Tet2mutDSS_IL1R_TOSICA_prediction_meta_table.csv", row.names=1)
Tet2_ADRB_IL1R <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/ADRB_IL1R/TOSICA_PreMeta_TET2_ADRB_IL1R.csv",row.names = 1)
library(Seurat)
library(ggplot2)
library(dplyr)
# setwd("D:/KS项目/公众号文章/堆叠柱状图显示比例")

Ratio <- Tet2_ADRB_IL1R %>%group_by(orig.ident,Prediction) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
mycolor<- c("#FED439FF", "#709AE1FF", "#8A9197FF", "blue",
            "#FD7446FF",  "#197EC0FF", "#D2AF81FF", 
            "#71D0F5FF", "#370335FF", "#075149FF",
            "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
names(x = mycolor) <- c("CLP", "ErP", "HSC1", 'HSC2', 'MDP', 'MEP', 'NMP', 'Pro_B', 'Pro_DC', 'Pro_Mast', 'Pro_Mk', 'Pro_Mono','Pro_NE')
ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = Prediction))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = mycolor)


WTDSS_TNF <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/WTDSS_TNFA_TOSICA_prediction_meta_table.csv", row.names=1)
Tet2DSS_ADRB <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/Tet2mutDSS_IL1R_TOSICA_prediction_meta_table.csv", row.names=1)
library(Seurat)
library(ggplot2)
library(dplyr)
# setwd("D:/KS项目/公众号文章/堆叠柱状图显示比例")

Ratio <- Tet2DSS_IL1R %>%group_by(orig.ident,Prediction) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
mycolor<- c("#FED439FF", "#709AE1FF", "#8A9197FF", "blue",
            "#FD7446FF",  "#197EC0FF", "#D2AF81FF", 
            "#71D0F5FF", "#370335FF", "#075149FF",
            "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
names(x = mycolor) <- c("CLP", "ErP", "HSC1", 'HSC2', 'MDP', 'MEP', 'NMP', 'Pro_B', 'Pro_DC', 'Pro_Mast', 'Pro_Mk', 'Pro_Mono','Pro_NE')
ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = Prediction))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = mycolor)

WT <- read.csv("~/hh028/Tet2_DSS_scRNAseq/BM/240204/TOSICA_PreMeta_WT.csv", row.names=1)
library(Seurat)
library(ggplot2)
library(dplyr)
# setwd("D:/KS项目/公众号文章/堆叠柱状图显示比例")

Ratio <- WT %>%group_by(orig.ident,Prediction) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
mycolor<- c("#FED439FF", "#709AE1FF", "#8A9197FF", "blue",
            "#FD7446FF",  "#197EC0FF", "#D2AF81FF", 
            "#71D0F5FF", "#370335FF", "#075149FF",
            "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
names(x = mycolor) <- c("CLP", "ErP", "HSC1", 'HSC2', 'MDP', 'MEP', 'NMP', 'Pro_B', 'Pro_DC', 'Pro_Mast', 'Pro_Mk', 'Pro_Mono','Pro_NE')
ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = Prediction))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = mycolor)
