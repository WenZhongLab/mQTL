# Seasonal pattern, sex- and BMI- association analysis and visualization
# corresponding to figure 2c-h

# --------------Mix effect model + KR-----------------
met527name <- unique(met527$metabolite)

pv_SABS_month <- data.frame(matrix(ncol = 15, nrow = 527))
colnames(pv_SABS_month) <- c("metabolite", 
                       "month", "month_sin_beta", "month_sin_se", "month_cos_beta", "month_cos_se",
                       "sex", "sex_betaf", "sex_sef",
                       "BMI", "BMI_beta", "BMI_se", 
                       "age", "age_beta", "age_se")
pv_SABS_month$metabolite <- met527_anno$metabolite
rownames(pv_SABS_month) <- pv_SABS_month$metabolite




for(x in met527name){
  da <- met527_timelabel[met527_timelabel$metabolite == x,]
  da$Gender <- factor(da$Gender, levels = c("m", "f") )
 
  # full model  
  a <- lme4::lmer(num_log2 ~ sin((2*pi*lable12ture)/12) + cos((2*pi*lable12ture)/12) +
                    Gender + BMI + Age_at_Visit  +  (1|sample), da)
  
  # get beta and SE
  pv_SABS_month[x,c("month_sin_beta", "month_cos_beta", "sex_betaf", "BMI_beta", "age_beta" )] <- a@beta[2:length(a@beta)]
  pv_SABS_month[x,c("month_sin_se", "month_cos_se", "sex_sef", "BMI_se", "age_se" )] <- 
    summary(a)$coefficients[2:length(a@beta),2]
  
  # season
  b <- lme4::lmer(num_log2 ~  Gender + BMI + Age_at_Visit  + (1|sample), da)
  kr <- pbkrtest::KRmodcomp(a, b)
  pv_SABS_month[x,"month"] <- kr$test[[5]][1]

  # gender
  fm1 <- update(a, .~. -Gender)
  kr <- pbkrtest::KRmodcomp(a, fm1)
  pv_SABS_month[x,"sex"] <- kr$test[[5]][1]

  # BMI
  fm2 <- update(a, .~. -BMI)
  kr2 <- pbkrtest::KRmodcomp(a, fm2)
  pv_SABS_month[x,"BMI"] <- kr2$test[[5]][1]
  
  # age
  fm3 <- update(a, .~. -Age_at_Visit)
  kr3 <- pbkrtest::KRmodcomp(a, fm3)
  pv_SABS_month[x,"age"] <- kr3$test[[5]][1]
  
}

# Multiple Testing Corrections
pv_SABS_month$month_q <- p.adjust(pv_SABS_month$month, method = "BH")
pv_SABS_month$sex_q <- p.adjust(pv_SABS_month$sex, method = "BH")
pv_SABS_month$BMI_q <- p.adjust(pv_SABS_month$BMI, method = "BH")
pv_SABS_month$age_q <- p.adjust(pv_SABS_month$age, method = "BH")


table(pv_SABS_month$month_q < 0.05)
table(pv_SABS_month$sex_q < 0.05)
table(pv_SABS_month$BMI_q < 0.05)



# ------------- visualization: amplitude -----------------
da <- pv_SABS_month
da$amplitude <- sqrt(da$month_sin_beta^2 + da$month_cos_beta^2)
da$log10p <- -log10(da$month_q)
da <- left_join(da, met527_anno[,c("metabolite","metabolite_name","class")])

metid <- c("lc_1","gc_170","gc_257","lc_12","lc_146","lc_192","gc_55")
p1 <- ggplot(da, aes(amplitude, -log10(month_q)))+
  geom_point(aes( color = class, size=log10p), alpha=0.5)+
  scale_color_manual( values = scale_fill_palette_class9)+
  geom_hline(yintercept = c(-log10(0.05)) , linetype = "dashed")+
  
  # example metabolite
  geom_point(da = da[da$metabolite %in% metid,], 
             mapping=aes(color=class, size=log10p))+
  geom_text_repel(da = da[da$metabolite %in%  metid,], 
                  mapping=aes(label=metabolite_name, color=class),size=4.25)+
  
  labs(x="Amplitude", y= "-Log10(adjusted p-value)", title = element_blank(), 
       color = "Class",size="-Log10(adj.p)") +
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position = "none",
        axis.title.x =element_text(size=12, colour = "black"), 
        axis.title.y =element_text(size=12, colour = "black"),
        axis.text.x =element_text(size=12, colour = "black"), 
        axis.text.y =element_text(size=12, colour = "black"))




# ------------- visualization: heatmap ---------------------------

da <- left_join(met527, visit_label)
da <- Rmisc::summarySE(da, measurevar="num_log2", groupvars=c("metabolite","month"),na.rm = T)

da_plot <- as.data.frame(pivot_wider(da[,c("metabolite","month","num_log2")], names_from = "month", values_from = "num_log2"))
rownames(da_plot) <- da_plot$metabolite
da_plot <- da_plot[,-c(1)]
 
a <- pv_SABS_month[which(pv_SABS_month$month_q < 0.05),]

anno <- unique(met527_anno[which(met527_anno$metabolite %in% a$metabolite),c(1,3)])
rownames(anno) <- anno$metabolite


anno <- data.frame(
  row.names = anno$metabolite,
  "class"=anno$class
)

ann_colors = list(
  class=scale_fill_palette_class9
)

heatmap_met121 <- pheatmap::pheatmap( t(scale(t(da_plot[a$metabolite,]))), treeheight_row=10,
                            annotation_row = anno, annotation_colors = ann_colors,
                            fontsize = 5, 
                            angle_col = "0",
                            cluster_cols = F, scale = "none",show_rownames=F,legend = F,
                            border_color = NA,clustering_method="ward.D2",cutree_rows=4,
                            color = colorRampPalette(colors = c("black","white","yellow3"))(100))





# -------------- seasonal pattern recognition ---------------
cluster <- heatmap_met121$tree_row

cut_527 <- data.frame(cutree(cluster,4))
 
cut_527$metabolite <- row.names(cut_527)
colnames(cut_527)[1] <- "cluster"


# -------------- visualization: seasonal pattern ----------------

da <- left_join(met527, visit_label)
da <- Rmisc::summarySE(da, measurevar="num_log2", groupvars=c("metabolite","month"),na.rm = T)

da_plot <- as.data.frame(pivot_wider(da[,c("metabolite","month","num_log2")], names_from = "month", values_from = "num_log2"))
rownames(da_plot) <- da_plot$metabolite
da_plot <- da_plot[,-c(1)]

da_plot$metabolite <- row.names(da_plot)
da_plot <- na.omit(left_join(da_plot, cut_527))
rownames(da_plot) <- da_plot$metabolite

da_for_example <- data.frame(t(scale(t(da_plot[a$metabolite,1:12]))))
colnames(da_for_example) <- c(1:12)
da_for_example$metabolite <- rownames(da_for_example)

p2=list()
index <- c(2,3,1,4)
for(i in 1:4){
  da <- da_for_example[which(da_for_example$metabolite %in% cut_527[which(cut_527$cluster == i),]$metabolite),1:13]
  da <- reshape2::melt(da)
  colnames(da)[2] <-"lable12ture"
  da$lable12ture <- as.numeric(da$lable12ture)
  da <- left_join(da, met527_anno[,c("metabolite","class")])
  temp <- Rmisc::summarySE(da, measurevar="value", groupvars=c("lable12ture"),na.rm = T)
  
  
  p2[[index[i]]] <- ggplot(da, aes(x=lable12ture, y=value))+
    geom_jitter(aes(color = class), width = 0.2)+
    geom_point(data=temp, mapping = aes(lable12ture,value), color = "grey30",size=2.5,shape=16)+
    geom_smooth(method="lm",  formula = y ~ sin((2*pi*x)/12) + cos((2*pi*x)/12), color = "grey10",se=F )+
    geom_errorbar(da=temp, mapping=aes(ymin=value-sd,ymax=value+sd),width=0.3, color="grey30")+ 
    scale_color_manual( values = scale_fill_palette_class9)+
    
    labs(x="Month", y= "Concentration (Scaled)", title = paste("Cluster M",index[i], sep =""),color="Sex") + 
    scale_x_continuous(breaks=unique(da$lable12ture), labels = unique(da$lable12ture))+
    
    theme_classic() +
    theme(panel.grid=element_blank(),
          panel.background = element_blank(),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position = "none",
          axis.title.x =element_text(size=12, color = "black"), 
          axis.title.y =element_text(size=12, color = "black"),
          axis.text.x =element_text(size=12, color = "black"),  
          axis.text.y =element_text(size=12, color = "black"),  
          legend.text=element_text(size=12, color = "black"), 
          legend.title=element_text(size=12, color = "black"))
   
}
 


# ----------------- visualization: seasonal pattern and enrichment analysis --------------

# enrichment analysis of metabolites belonging to each cluster was performed by using the ORA module of HUBMet

# result filtering
c1_enrich <- c1_enrich %>%  
  filter(!Pathway.Subject %in% c("Disease","Drug Action")) %>% 
  mutate(Cluster = "M2")

c2_enrich <- c2_enrich%>%  
  filter(!Pathway.Subject %in% c("Disease","Drug Action")) %>% 
  mutate(Cluster = "M3")

c3_enrich <- c3_enrich%>%  
  filter(!Pathway.Subject %in% c("Disease","Drug Action")) %>% 
  mutate(Cluster = "M1")

c4_enrich <- c4_enrich%>%  
  filter(!Pathway.Subject %in% c("Disease","Drug Action")) %>% 
  mutate(Cluster = "M4")


# prepare data for visualization
c_smpdb <- do.call(rbind, list(
  c1_enrich[2:5,],
  c2_enrich[1:4,],
  c3_enrich[c(1,2,4,5),],
  c4_enrich[1:4,]
))

c_smpdb$cluster <- c(rep(2,4),rep(3,4),rep(1,4),rep(4,4))
c_smpdb <- c_smpdb[c(9:12,1:8,13:16),]

c_smpdb$xlab <- as.factor(nrow(c_smpdb):1)
c_smpdb$logp <- -log10(c_smpdb$pvalue)
c_smpdb$logadjp <- -log10(c_smpdb$p.adjust)

color_season4 <- c( "#67D5B5", "#ee7785", "#c89ec6", "#84e1ed")[4:1]

p3 <- ggplot(c_smpdb,aes(y=logadjp,x=xlab,fill=as.factor(cluster))) + 
  geom_bar(stat="identity",position = "dodge") +
  geom_text(aes(label = Pathway.Name, y=0.25),hjust = 0, size = 4.2, color = "black", alpha = 1)+
  scale_fill_manual(values =color_season4) + 
  coord_flip() + 
  geom_vline(xintercept = c(4.5 , 8.5, 12.5), linetype = "dashed")+
  scale_y_continuous(expand = c(0,0))+ 
  labs(fill="Cluster", x="SMPDB term", y="-Log10(adjusted p-value)")+
  
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.title.x =element_text(size=12, colour = "black"), 
        axis.title.y =element_text(size=12, colour = "black"),
        axis.text.x =element_text(size=12, colour = "black"), 
        axis.text.y =element_blank(),
        legend.text=element_text(size=12, colour = "black"), 
        legend.title=element_text(size=12, colour = "black"))
