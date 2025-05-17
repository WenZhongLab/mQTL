
# co-expression clustering analysis, corresponding to figure 1b-d

# da_forcluster: rows-metabolite, columns-samples 

# ------------------UMAP-----------------------

library(umap)
library(tidyverse)

custom.config_cluster <- umap.defaults
custom.config_cluster$random_state = 14 

u1 <- umap( t(scale(t(da_forcluster))) , custom.config_cluster) 
umapdata_527 <- data.frame(u1$layout) %>% 
  rownames_to_column(var = "metabolite") %>% 
  left_join(met527_anno[,c("metabolite","class","subclass")], by="metabolite")



# ---------------------- cluster KNN -------------------

# KNN
mat <- scale(t(da_forcluster))

kn <- 20
knn.info <- RANN::nn2( t(mat ), k=kn)

# get adjacency matrix
knn <- knn.info$nn.idx
adj_527 <- matrix(0, ncol(mat), ncol(mat))
rownames(adj_527) <- colnames(adj_527) <- colnames(mat)
for(i in seq_len(ncol(mat))) {
  adj_527[i,colnames(mat)[knn[i,]]] <- 1
}

# ---------------------- cluster SNN ---------------------

# SNN
sn <- 5
snn.info <- dbscan::sNN(adj_527, k = sn) 

snn.info.class <- dbscan::adjacencylist(snn.info)

# get adjacency matrix
adj_snn_527 <- matrix(0, ncol(mat), ncol(mat))
rownames(adj_snn_527) <- colnames(adj_snn_527) <- colnames(mat)
for(i in seq_len(ncol(mat))) {
  adj_snn_527[i,colnames(mat)[snn.info.class[[i]]]] <- 1
}


# Create graphs from adjacency matrices of SNN
g_snn_527 <- igraph::graph.adjacency(adj_snn_527, mode="undirected")
# remove self loops
g_snn_527 <- igraph::simplify(g_snn_527) 


# -------------------- cluster: louvain find community -----------------


# community detection
set.seed(12)
km_527 <- igraph::cluster_louvain(g_snn_527)

# community membership
community_527 <- km_527$membership
names(community_527) <- km_527$names

primary_cluster <- data.frame(community_527) %>% 
  rownames_to_column(var="metabolite")
colnames(primary_cluster)[1] <- "primary_cluster"


# -------------------cluster: hierarchical clustering----------------------

# calculate mean abundance of each cluster for hierarchical clustering
da_cluster <- left_join(met527, primary_cluster)


daforcluster <- Rmisc::summarySE(da_cluster[da_cluster$iid %in% colnames(da_forcluster),c("iid","num_log2","primary_cluster")], 
                                 measurevar="num_log2", groupvars=c("primary_cluster","iid"), na.rm = T)

daforcluster_wide <- data.frame(pivot_wider(daforcluster[,c("iid","num_log2","primary_cluster")],
                                            names_from ="iid", values_from = "num_log2" ))

rownames(daforcluster_wide) <- paste("cluster", daforcluster_wide$primary_cluster, sep = "")
daforcluster_wide <- daforcluster_wide[,-c(1)]

# hierarchical clustering
out.dist_527 <- dist(t(scale(t(daforcluster_wide))),method="euclidean") 

out.hclust_527 <- hclust(out.dist_527, method="ward.D2")


# cut tree
cou.id_527 <- cutree(out.hclust_527, k = 8)

hc_res_527 <- data.frame(cou.id_527)
hc_res_527$primary_cluster <- unlist(lapply(rownames(hc_res_527),
                                            function(x){strsplit(x, split = "cluster")[[1]][2]}))
colnames(hc_res_527)[1] <- "hcluster"
hc_res_527 <- data.frame(hc_res_527) %>% 
  mutate(primary_cluster = as.numeric(primary_cluster),
         hcluster = as.numeric(hcluster))


met527_cluster_result <- unique(da_cluster[,c("metabolite","primary_cluster")]) %>% 
  left_join(  hc_res_527, by = "primary_cluster")



# ------------ visualization: UMAP -------------------

hcluster_color <- c("#CC0C00", "#5C88DA", "#8ce600", "#FFCD00", "#7C878E", "#87cefa", "#00AF66", "#b553f4", "#e69995")
hcluster_color_dot <- c("#a30a00", "#2e64ca", "#70b800", "#cca400", "#4f565b", "#3db0f7", "#007042", "#690ba3", "#d85e58")

# prepare data for visualization
da_forUMAPvisu <- umapdata_527
da_forUMAPvisu <- left_join(da_forUMAPvisu, met527_cluster_result)
da_forUMAPvisu$primary_cluster <- factor(da_forUMAPvisu$primary_cluster)
da_forUMAPvisu$class <- factor(da_forUMAPvisu$class, levels = class9)
da_forUMAPvisu$hcluster <- factor(da_forUMAPvisu$hcluster)
da_forUMAPvisu <- left_join(da_forUMAPvisu, met527_anno[,c("metabolite",'metabolite_name')])


# figure 1: color by cluster
p1 <- ggplot(da_forUMAPvisu, aes(X1, X2)) +
  geom_point(aes(col =  hcluster), size = 2.5)+
  scale_color_manual(values = hcluster_color) +
  labs(x="UMAP1", y= "UMAP2", title = element_blank(), color = "Clusters") +
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position = "right",
        axis.title.x =element_text(size=15, colour = "black"), 
        axis.title.y =element_text(size=15, colour = "black"),
        axis.text.x =element_text(size=15, colour = "black"), 
        axis.text.y =element_text(size=15, colour = "black"),
        legend.text=element_text(size=15, colour = "black"), 
        legend.title=element_text(size=15, colour = "black"))



# figure 2: color by class

library(ggnewscale)
library(ggrepel)

expid_less <- c("lc_192", "SM_100", "gc_257", "gc_98",  "TG_46",  "lc_75",   "LPC_91","gc_53")

p2 <- ggplot(da_forUMAPvisu, aes(X1, X2)) +
  geom_point(aes(col =  class), size = 2.5)+
  scale_color_manual(values = color9_c) + 
  
  # add example metabolites
  new_scale_color()+
  scale_color_manual(values = hcluster_color_dot) +
  geom_point(data = subset(da_forUMAPvisu, metabolite %in% expid_less),
             aes(col = hcluster), alpha=1, size = 5, shape = 17) +
  geom_text_repel(data = subset(da_forUMAPvisu, metabolite %in% expid_less),
                  aes(color=hcluster, label = metabolite_name), size = 5) +
  
  # add ellipse
  new_scale_color()+
  scale_color_manual(values = hcluster_color) +
  stat_ellipse(aes(color = hcluster), type = "t", linewidth=0.3) +
  
  
  labs(x="UMAP1", y= "UMAP2", title = element_blank(), color = "Clusters") +
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position = "none",
        axis.title.x =element_text(size=15, colour = "black"), 
        axis.title.y =element_text(size=15, colour = "black"),
        axis.text.x =element_text(size=15, colour = "black"), 
        axis.text.y =element_text(size=15, colour = "black"),
        legend.text=element_text(size=15, colour = "black"), 
        legend.title=element_text(size=15, colour = "black"))





# ----------- visualization: sankey plot ------------

da <- met527_cluster_result %>% 
  left_join(  met527_anno[,c("metabolite","metabolite_name","class","subclass")],  by ="metabolite")
da$connection <- paste(da$hcluster, da$class, sep = "-")
da[which(da$class == "Lipid"),"connection"] <- paste(da[which(da$class == "Lipid"),]$connection, 
                                                     da[which(da$class == "Lipid"),]$subclass, sep = "-")
da$connection <- gsub(" ", "-", da$connection)
da <- data.frame(table(da$connection))
da$Var1 <- as.character(da$Var1)

da$cluster <- as.character(unlist(lapply(da$Var1, function(x){strsplit(x, split = "-")[[1]][1]})))
da$class <- as.character(unlist(lapply(da$Var1, function(x){substr(x,3,nchar(x))})))
da$color <- as.character(unlist(lapply(da$Var1, function(x){strsplit(x, split = "-")[[1]][2]})))

da[which(da$color == "Lipid" & da$Freq < 4),4] <- "Lipid-Others"
da$Var1 <- paste(da$cluster, da$class, sep = "-")

for(i in 1:8){
  da[which(da$cluster == i & da$class == "Lipid-Others"),"Freq"] <- sum(da[which(da$cluster == i & da$class == "Lipid-Others"),]$Freq)
}

da <- unique(da)

da_edge <- da[,c(4,3,2,5)]
colnames(da_edge) <- c("source","target","value","color")

da_node <- unique(da_edge[,c(1,4)])
temp <- data.frame(scale_fill_palette_class9)
temp$name <- as.character(unlist(lapply(rownames(temp), function(x){strsplit(x, split = " ")[[1]][1]})))
  
colnames(temp)[1] <- "color"
da_node <- left_join(da_node, temp, by=c("color"="name"))
da_node <- da_node[,c(1,3)]
colnames(da_node) <- c("name","color")
temp <- data.frame("color" = hcluster_color[1:8],
                   name = as.character(1:8))


da_node <- rbind(da_node, temp[,c(2,1)])

da_node <- da_node[c(6,7,5,1,10,4,9,3,2,8,11,13,12,14,15:22),]

da_edge$target <- paste("Cluster-",da_edge$target, sep="")
da_node[which(da_node$name %in% 1:8),1] <- paste("Cluster-",da_node[which(da_node$name %in% 1:8),1], sep="")
da_node$node_group <- gsub(" ", "-", da_node$name)
da_node[da_node=="Lipid-Steroids-and-steroid-derivatives"] <- "Lipid-Steroids"

my_color <- 'd3.scaleOrdinal() .domain(["Amino-acid","Carbohydrate","Cofactors-and-vitamins","Energy","Lipid-Others","Lipid-Fatty-Acyls","Lipid-Glycerolipids","Nucleotide","Others","Peptide","Xenobiotics","Lipid-Steroids","Lipid-Glycerophospholipids","Lipid-Sphingolipids","Cluster-1","Cluster-2","Cluster-3","Cluster-4","Cluster-5","Cluster-6","Cluster-7","Cluster-8"]) .range(["#45aee8","#f57070","#add547","#ad85d5","#ffa510","#ffa510","#ffa510","#41b7ac","#C3DEE0","#f6c9aa","#EB99C2","#ffa510","#ffa510","#ffa510","#CC0C00","#5C88DA","#8ce600","#FFCD00","#7C878E","#87cefa","#00AF66","#b553f4"])'


da_edge$IDsource <- match(da_edge$source, da_node$name)-1
da_edge$IDtarget <- match(da_edge$target, da_node$name)-1

p3 <- networkD3::sankeyNetwork(Links = da_edge, Nodes = da_node, Source = "IDsource",iterations = 0,
                                 fontFamily = "Arial", 
                                 Target = "IDtarget", nodePadding= 2,Value = "value", 
                                 NodeID = "name",NodeGroup = "node_group",height = 400, width=375,
                                 fontSize = 14, unit = "Letter(s)",  colourScale = my_color)


# --------------- enrichment analysis -----------------

term75 <- met527_anno[,c("metabolite","customed_term")]
terms_id <- as.character(data.frame(table(term75$customed_term))$Var1)

# enrichment analysis: fisher's exact analysis
go_custom_527 <- data.frame(matrix(ncol = 7,nrow = 1))
colnames(go_custom_527) <- c("cluster", "cluster_metN", "term",  "total", "hit","pvalue", "OR")

n=0
for(i in 1:max(met527_cluster_result$hcluster)){
  da <- met527_cluster_result[met527_cluster_result$hcluster == i,]
  for(term_x in terms_id){
    n=n+1
    term_met <- term75[term75$customed_term == term_x,1]
    overlap_n <- length(intersect(term_met, da$metabolite))
    
    datax <- matrix(c(overlap_n, nrow(da)-overlap_n, 
                      length(term_met)-overlap_n, 527-nrow(da)-length(term_met)+overlap_n ), nrow = 2)
    pp <- fisher.test(datax, alternative = "greater")
    go_custom_527[n,] <- c(i, nrow(da), term_x, length(term_met), overlap_n, 
                           pp$p.value, pp$estimate[["odds ratio"]])
  }
  
}

for(i in c(1,2,4,5,6,7)){
  go_custom_527[,i] <- as.numeric(go_custom_527[,i] )
}


go_custom_527$padj_BH <- p.adjust(go_custom_527$pvalue, method = "BH")
go_custom_527 <- go_custom_527[go_custom_527$hit > 0,]
go_custom_527 <- go_custom_527[order(go_custom_527$pvalue, decreasing = F),]


# ----------- visualization: enrichment analysis --------------

# get first 10 significant terms for each cluster
for(i in 1:max(met527_cluster_result$hcluster)){
  h1 <- go_custom_527[which(go_custom_527$cluster == i),]
  row.names(h1) <- 1:nrow(h1)
  h1 <- na.omit(h1[c(1:10),])
  if(i == 1){
    go_custom_select_527 <- h1
  }else{
    go_custom_select_527 <- rbind(go_custom_select_527, h1)
  }
  
}


# result filtering
go_custom_select_527 <- go_custom_select_527[which(go_custom_select_527$total != 1),]
go_custom_select_527 <- go_custom_select_527[which(go_custom_select_527$hit != 1),]


# prepare data for visualization
go_custom_select_527$lable <- go_custom_select_527$term
rownames(go_custom_select_527) <- 1:49
# select top significant terms for visualization
go_custom_select_527 <- go_custom_select_527[c(1,3,4,9,11,12,15:17,21,23,25,31:33,38,39,40,41,44:46),]

rownames(go_custom_select_527) <- 1:nrow(go_custom_select_527)
go_custom_select_527$xlab <- factor(as.numeric(rownames(go_custom_select_527)), levels = nrow(go_custom_select_527):1 )
go_custom_select_527$logp <- -log10(go_custom_select_527$pvalue)
go_custom_select_527$logq_bh <- -log10(go_custom_select_527$padj_BH)

go_custom_select_527$cluster <- factor(go_custom_select_527$cluster)


p4 <- ggplot(go_custom_select_527,
               aes(y=logp,x=xlab,fill=cluster)) + 
  geom_bar(stat="identity",position = "dodge") +
  geom_text(aes(label = lable, y=0.25),hjust = 0, size = 5, color = "black", alpha = 1)+
  scale_fill_manual(values = hcluster_color) +
  geom_vline(xintercept = c(3.5 , 5.5, 7.5, 10.5, 13.5, 16.5,19.5), linetype = "dashed")+
  coord_flip() + 
  scale_y_continuous(expand = c(0,0))+ 
  labs(fill="Cluster", x=element_blank(), y="-Log10(p-value)")+
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.title.x =element_text(size=15, colour = "black"), 
        axis.title.y =element_text(size=15, colour = "black"),
        axis.text.x =element_text(size=15, colour = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())




