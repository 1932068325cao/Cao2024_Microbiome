#加载包
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)
library(patchwork)
#读取 OTUs 丰度表
df <- read.table('zotu_table_norm.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)
#使用vegan包计算多样性指数
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)#spe.rich =sobs
###将以上多样性指数统计成表格
diversity <- as.data.frame(cbind(Shannon, Simpson, Richness))
#导入经纬度
geo <- read.table("geo.txt", row.names = 1,sep = '\t',header = T)
aa <- cbind(geo, diversity)
#按照分组提取数据
b<-subset(aa,aa$group=="N")
c<-subset(aa,aa$group=="M")
#计算松材线虫的纬度多样性
summary(lm(b$Richness ~ b$Latitude))
summary(lm(b$Shannon ~ b$Latitude))
p1 <- ggplot(b,aes(x = Latitude, y = Richness))+
  geom_point(size = 3)+
  geom_smooth(method = "lm",formula = 'y ~ x',alpha = 0.2)
p1
ggsave('松材线虫Richness纬度多样性.pdf',p1,width=5.5,height=5.5)


p2 <- ggplot(b,aes(x = Latitude, y = Shannon))+
  geom_point(size = 3)+
  geom_smooth(method = "lm",formula = 'y ~ x',alpha = 0.2)
ggsave('松材线虫Shannon纬度多样性.pdf',p2,width=5.5,height=5.5)
#计算松树纬度多样性
summary(lm(c$Richness ~ c$Latitude))
summary(lm(c$Shannon ~ c$Latitude))
p3 <- ggplot(c,aes(x = Latitude, y = Richness))+
  geom_point(size = 3)+
  geom_smooth(method = "lm",formula = 'y ~ x',alpha = 0.2)
ggsave('松树Richness纬度多样性.pdf',p3,width=5.5,height=5.5)

p4 <- ggplot(c,aes(x = Latitude, y = Shannon))+
  geom_point(size = 3)+
  geom_smooth(method = "lm",formula = 'y ~ x',alpha = 0.2)
ggsave('松树Shannon纬度多样性.pdf',p4,width=5.5,height=5.5)
