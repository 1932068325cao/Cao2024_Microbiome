##计算微生物丰度间的相关系数
library(Hmisc)
# 读取otu-sample矩阵，行为sample，列为otu
otu <- read.table('ALL_freqs_M.txt', sep='\t',head=T, row.names=1)
otu <-t(otu)
#计算相关性
sp.cor<-rcorr(t(otu),type="spearman")
#提取r、p值矩阵；
occor.r<-sp.cor$r
occor.p <- sp.cor$P
#使用Benjamini-Hochberg("FDR-BH")法进行多重检验校正；
p <- p.adjust(occor.p, method="BH")
# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.7] = 0
diag(occor.r) <- 0
#将occor.r保存为csv文件
write.table(occor.r,file="相关性计算结果.txt",sep="\t")
#根据上述筛选的 r 值和 p 值保留数据
z <- occor.r * occor.p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
z[abs(z)>0]=1
z
adjacency_unweight <- z

##igraph 包计算网络模块
library(igraph)
#输入数据示例，邻接矩阵
#这是一个微生物互作网络，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
head(adjacency_unweight)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph    #igraph 的邻接列表

#计算节点度
V(igraph)$degree <- degree(igraph)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块

write.table(nodes_list, 'nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#上述的邻接矩阵类型的网络文件
adjacency_unweight 

#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, 'zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.csv(zi_pi,"zipi结果.csv")
ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.8, linewidth = 6,shape=17) +
  scale_y_continuous(limits=c(-2,3))+
  scale_color_manual(values = c("#8491B4FF","#91D1C2FF","#F39B7FFF", "#4DBBD5FF"), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62,linetype=2,linewidth=1) +
  geom_hline(yintercept = 2.5,linetype=2,linewidth=1)+  
  theme_bw()+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
ggsave('zipi_m.pdf',width=5.5,height=5.5)



