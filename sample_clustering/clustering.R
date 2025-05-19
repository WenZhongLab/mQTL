#Metabolic syndrome variables
metsyn.var <- c("HDL", "BMI", "SBP", "TG", "Gluc")

metsyn.mat <- metadata[,metsyn.var] %>% as.data.frame() #all visits
rownames(metsyn.mat) <- metadata$id
metsyn.mat <- zoo::na.aggregate(metsyn.mat, FUN=mean)
rownames(metsyn.mat) <- metadata$id

#remove sex effect
idx.f <- which(metadata$Gender=="f")
idx.m <- which(metadata$Gender=="m")
for (n in 1:ncol(metsyn.mat)){
  df <- data.frame(y=metsyn.mat[,n], sex=metadata$Gender)
  lmFit <- lm(y~0 + sex, df)
  
  metsyn.mat[idx.f,n] <- metsyn.mat[idx.f,n]-coef(lmFit)["sexf"]
  metsyn.mat[idx.m,n] <- metsyn.mat[idx.m,n]-coef(lmFit)["sexm"]
  metsyn.mat[,n]  <- metsyn.mat[,n] / sd(metsyn.mat[,n] )
}


#Identify different groups based on metabolic syndrome parameters ####
# K-means
n.cluster <- 2
color.cluster <- c(rgb(85/256,221/256,255/256), rgb(200/256,55/256,55/256))

cl.kmeans <- kmeans(metsyn.mat[,metsyn.var], n.cluster, iter.max=500)$cluster

#assign 1 to High risk
if (median(metadata$BMI[cl.kmeans==1]) > median(metadata$BMI[cl.kmeans==2])){
  cl.kmeans <- ifelse(cl.kmeans==1, 1, 0)
} else {
  cl.kmeans <- ifelse(cl.kmeans==2, 1, 0)
}

metsyn.mat$kmeans <- as.character(cl.kmeans) #store groups

#Bar plot of other health variables
health.var <- c("ALAT", "GGT", "Urate", "TNT", "WBC")
health.mat <- metadata[,health.var] %>% as.data.frame()
health.mat <- zoo::na.aggregate(health.mat, FUN=mean)

#remove sex effect
idx.f <- which(metadata$Gender=="f")
idx.m <- which(metadata$Gender=="m")
for (n in 1:ncol(health.mat)){
  df <- data.frame(y=health.mat[,n], sex=metadata$Gender)
  lmFit <- lm(y~0 + sex, df)
  
  health.mat[idx.f,n] <- health.mat[idx.f,n]-coef(lmFit)["sexf"]
  health.mat[idx.m,n] <- health.mat[idx.m,n]-coef(lmFit)["sexm"]
  health.mat[,n]  <- health.mat[,n] /sd(health.mat[,n] )
}

#add clustering information
health.mat$kmeans <- metsyn.mat$kmeans

df.plot <- data.frame()
df.plot.full <- data.frame()
cluster.uniq <- unique(metsyn.mat$kmeans)
for (n in 1:n.cluster){
  for (m in 1:length(health.var)){
    
    values <- health.mat[health.mat$kmeans==cluster.uniq[n], health.var[m]] %>% na.omit()
    df.plot.full <- rbind(df.plot.full, data.frame(values=values, group=rep(cluster.uniq[n], length(values)),
                                                   var=rep(health.var[m], length(values))))
    
    #summary
    values.mean <- mean(values)
    values.ths <- quantile(values, c(0.05,.5,.95), na.rm = T)
    values.sd <- sd(values, na.rm = T)
    df.plot <- rbind(df.plot, list(group=cluster.uniq[n], var=health.var[m], 
                                   ths5=values.ths[1], ths50=values.ths[2], ths95=values.ths[3], 
                                   sd=values.sd, mean=values.mean))
  }
}
df.plot$group <- as.character(df.plot$group)
df.plot.full$group <- as.character(df.plot.full$group)
df.plot.full$group <- factor(df.plot.full$group, levels=cluster.uniq)
df.plot$group <- factor(df.plot$group, levels=cluster.uniq)


comb.group <- combn(cluster.uniq,2)
comparisons <- vector(mode="list", ncol(comb.group))
for (n in 1:ncol(comb.group)){
  comparisons[[n]] <- c(comb.group[1,n], comb.group[2,n])
}

df.plot.full %>%
  ggplot(aes(x=group, y=values, fill=group)) +
  geom_boxplot() + xlab("") + ylab("") +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_jitterdodge(0.2), alpha=0.5) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), comparisons=comparisons, method="t.test") +
  scale_y_continuous(expand = expansion(mult = .1)) +
  facet_wrap(~var, scales="free") + theme_minimal() +
  theme(text=element_text(size=20, family="Arial"))

