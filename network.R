##计算微生物丰度间的相关系数
library(Hmisc)

#以属水平丰度为例，“genus_table.txt” 是一个属水平的微生物丰度表
genus <- read.table('zotu_table_freqs_norm_N.CSV', row.name = 1, header=T,check.names = FALSE,sep=',')
#genus<-t(genus)
#可选事先过滤一些低丰度或低频的类群
genus <- genus[which(rowSums(genus) >= 0.005), ]    #例如只保留相对丰度总和高于 0.005 的属

genus1 <- genus
genus1[genus1>0] <- 1
genus <- genus[which(rowSums(genus1) >= 10), ]    #例如只保留在 5 个及以上样本中出现的属

genus<-t(genus)
write.table(genus,'ALL_freqs_N.txt', col.names = NA, sep = '\t', quote = FALSE)
#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'OTU_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('genus_taxonomy.txt',check.names = FALSE, row.name=1,stringsAsFactors = FALSE)
tax <- tax[as.character(V(g)$name), ]

V(g)$kingdom <- tax$kingdom
V(g)$kingdom
V(g)$phylum <- tax$phylum
V(g)$class <- tax$class
V(g)$order <- tax$order
V(g)$family <- tax$family
V(g)$genus <- tax$genus
g
plot(g)
##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.CSV', col.names = NA, sep = ',', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)

write.table(edge_list, 'network.edge_list.CSV', sep = ',', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(g)),
  kingdom = V(g)$kingdom,
  phylum = V(g)$phylum,
  class = V(g)$class,
  order = V(g)$order,
  family = V(g)$family,
  genus = V(g)$genus
)
head(node_list)

write.table(node_list, 'network.node_list.CSV', sep = ',', row.names = FALSE, quote = FALSE)
#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'network.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')
