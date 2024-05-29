#setwd("D:/R/R Script/23.4.2测序样品结果绘图/Mantel test")
library(vegan)
install.packages('geosphere')
library(geosphere)
#读取上述数据集
df <- read.csv('genus_zj_wood.csv', header= TRUE,sep=',')

##计算距离
#根据物种丰度数据，计算样方间的 Bray-curtis 距离
abund <- df[ ,4:ncol(df)]
dist.abund <- vegdist(abund, method = 'bray',na.rm = TRUE)
#根据经纬度，计算样方间实际的地理距离
geo <- data.frame(df$Longitude, df$Latitude)
d.geo <- distm(geo, fun = distHaversine)       #library(geosphere)
dist.geo <- as.dist(d.geo)
#物种丰度和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
abund_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 999, na.rm = TRUE)
abund_geo

library(ggplot2)
#将上文获得的距离测度，转化为数据框，一一对应起来
aa <- as.vector(dist.abund)
gg <- as.vector(dist.geo)
mat <- data.frame(aa, gg)
#基于物种丰度的距离与样方间地理距离之间的相关性散点图，上文已知二者无相关性
mm <- ggplot(mat, aes(y = aa, x = gg/1000)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2,formula = 'y ~ x') + 
  labs(x = "Physical separation (km)", y = "Bray-Curtis Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))

mm
ggsave('geo-abundance_zj_wood.pdf',mm,width = 5.5,height = 5.5)
