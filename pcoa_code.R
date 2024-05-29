rm(list=ls())#clear Global Environment
#setwd('D:\\桌面\\公众号运营\\自动回复及关键词\\R可视化\\PCoA')#设置工作路径

#安装所需R???
#install.packages("vegan")
#install.packages("ggplot2")
#加载???
library(vegan)#计算距离时需要的???
library(ggplot2)#绘图???

#读取数据，一般所需是数据行名为样本名、列名为OTUxxx的数据表
otu_raw <- read.table(file="zotu_table_freqs_norm_N.csv",sep=",",header=T,check.names=FALSE ,row.names=1)
#由于排序分析函数所需数据格式原因，需要对数据进行转置
otu <- t(otu_raw)

#计算bray_curtis距离
otu.distance <- vegdist(otu)
#pcoa分析#PCoA 函数中本身也提供了校正参数（add），避免负特征值的产生
#Cailliez 校正
pcoa <- cmdscale (otu.distance,k = (nrow(otu) - 1),eig=TRUE)

pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释???

###绘图###
#pcl2原来是matrix,转化为data.frame
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$samples <- row.names(pc12)
head(pc12)

#绘图
p <- ggplot(pc12,aes(x=V1, y=V2))+#指定数据、X轴、Y???
  geom_point(size=3)+#绘制点图并设定大小???
  theme_bw()#主题
p

#读入分组文件
group <- read.table("group_n.txt", sep='\t', header=T)
#修改列名
colnames(group) <- c("samples","group")
#将绘图数据和分组合并
df <- merge(pc12,group,by="samples")
head(df)
#library(ggsci)
#绘图
color=c("#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","gray",'#2baeb5',"#4DBBD5FF","#00A087FF","#DC0000FF" 
)#颜色变量
p1 <-ggplot(data=df,aes(x=V1,y=V2,
                       color=group))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=3)+#绘制散点图的相关参数，大小，点的形状，点的颜色，点的大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+#图中虚线
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=2)+#添加数据点的标签
 # guides(color=guide_legend(title=NULL))+#去除图例标题
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+#将x、y轴标题改为贡献度
  #stat_ellipse(data=df,
               #geom = "polygon",level=0.95,
               #linetype = 2,linewidth=0.5,
               #aes(fill=group),
               #alpha=0.2,
               #show.legend = T)+
  scale_color_manual(values = color) +#点的颜色设置
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=12),#修改X轴标题文???
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文???
        axis.text.y=element_text(size=10),#修改x轴刻度标签文???
        axis.text.x=element_text(size=10),#修改y轴刻度标签文???
        panel.grid=element_blank())#隐藏网格???
p1


#执行 PERMANOVA，详??? ?adonis
##使用 vegan 包的置换多元方差分析（PERMANOVA），比较用药前后两组中肠道菌群结构是否存在不???

#执行 PERMANOVA，详??? ?adonis
#样本间的成对距离选择使用物种多样性中最常使用的 Bray-curtis 距离，并通过 999 次置换估??? p ???
adonis_drug <- adonis2(otu~group, df, distance = 'bray', permutations = 999)
adonis_drug
#在图中添加上文的 PERMANOVA 的结果（P 值）
p1 <- p1 +labs(title=paste0("adonis R2: ",round(adonis_drug$R2,2), "; P: ", adonis_drug$`Pr(>F)`))
  
p1
ggsave('PCOA-otu-N.pdf',width=5.5,height =5.5 )
ggsave('PCOA.png',width=5.5,height = 5.5)
#计算不同分组内部的离散程度，如果P>0.05则说明，adonis的结果是可信的
dispersion <- betadisper(otu.distance, group$group)
permutest(dispersion)

