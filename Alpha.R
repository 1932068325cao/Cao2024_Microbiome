#安装包
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("ggsignif")
install.packages("vegan")
install.packages("ggprism")
install.packages("dplyr")
install.packages("RColorBrewer")
#加载包
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)

##############α-多样性指数的计算###############
##导入数据，所需是数据行名为样本名、列名为OTUxxx的数据表，输入表中必须全为整数
df <- read.delim("zotu_table_norm.txt",header = T, row.names = 1, check.names = F)
#使用vegan包计算多样性指数
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)#spe.rich =sobs
###将以上多样性指数统计成表格
index <- as.data.frame(cbind(Shannon, Simpson, Richness))
tdf <- t(df)#转置表格
tdf<-ceiling(as.data.frame(t(df)))
#计算obs，chao，ace指数
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]#统一行名
#将obs，chao，ace指数与前面指数计算结果进行合并
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
#计算Pielou及覆盖度
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)
#导出表格
write.table(cbind(sample=c(rownames(index)),index),'diversity.index.txt', row.names = F, sep = '\t', quote = F)
#读入文件，添加分组信息
df1 <- read.table('diversity.index_n.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

#使用 ggpubr（ggplot2 的扩展包），作图展示 wilcox 比较的结果

p<-ggboxplot(data = df1, x = 'group', y = 'Ace', color = 'black',fill="group", legend="none",palette=c( "#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","gray",'#2baeb5',"#4DBBD5FF","#00A087FF","#DC0000FF" )
             , add.params=list(size=1.2))+
  stat_compare_means(label.x = 2.5,label.y = 600)
p

ggsave('Ace_n_16s.pdf',p, height = 5.5, width = 5.5)

