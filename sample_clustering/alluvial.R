# Are samples from the same individual in the same group? ####
str <- strsplit(rownames(metsyn.mat), ":") %>% unlist()
ind <- str[seq(1,length(str),2)]
visit <- str[seq(2,length(str),2)]

group.uniq <- unique(metsyn.mat$kmeans)
visit.uniq <- unique(metadata$visit)

visit.list <- vector(mode="list", length(visit.uniq))
for (n in 1:length(visit.uniq)){
  visit.list[[n]] <- group.uniq
}
names(visit.list) <- paste0("visit", visit.uniq)

#all combinations across the visits
all.comb <- expand.grid(visit.list)  

#group across visits, per individual
df.ind <- data.frame(subject=unique(ind))
df.ind <- cbind(df.ind, matrix(NA, nrow(df.ind), length(visit.uniq)))
for (n in 1:nrow(df.ind)){
  for (m in 1:length(visit.uniq)){
    idx <- which(ind==df.ind$subject[n] & visit==as.character(m))
    if (length(idx)>0){
      df.ind[n,m+1] <- metsyn.mat$kmeans[idx]
    }
  }
}
colnames(df.ind)[2:ncol(df.ind)] <- paste0("visit", 1:length(visit.uniq))


# Alluvial plot
df.plot <- data.frame()
for (n in 1:6){
  x <- paste0("visit",n) %>% rep(nrow(df.ind)) #ith visit
  node <- df.ind[,n+1] #status at ith visit
  
  if (n<6){ #next visit and next status
    next_x <- paste0("visit",n+1) %>% rep(nrow(df.ind)) #ith visit
    next_node <- df.ind[,n+2]
  } else {
    next_x <- rep(NA, nrow(df.ind))
    next_node <- rep(NA, nrow(df.ind))
  }
  
  df.n <- data.frame(x=x, node=node, next_x=next_x, next_node=next_node)
  df.plot <- rbind(df.plot, df.n)
}
df.plot <- df.plot[which(!is.na(df.plot$node)),]


df.plot$x <- factor(df.plot$x, levels=paste0("visit", 1:6))
df.plot$next_x <- factor(df.plot$next_x, levels=paste0("visit", 1:6))
ggplot(df.plot, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = node,
               label=node)) +
  geom_alluvial(flow.alpha = .6, width=0.2) +
  geom_alluvial_text(size = 5, color = "white") + 
  scale_fill_manual(values=c("#55ddffff", "#c83737ff")) + 
  theme_minimal() +  
  theme(text=element_text(size=20, family="Arial"), plot.title = element_text(hjust = 0.5)) + 
  xlab("") + ylab("Individual")
