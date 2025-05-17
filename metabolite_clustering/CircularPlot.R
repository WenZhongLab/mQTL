# ---------------- figure 1a circle plot -------------------

library(tidyverse)

# class 9
color9_c <- c("#ffa510", "#45aee8", "#EB99C2", "#f6c9aa", "#f57070", "#add547", "#41b7ac", "#ad85d5", "#C3DEE0")

color9dot_c <- c("#d98700", "#0a6a9e", "#de589b", "#ee995f", "#cf2a2a", "#7aa300", "#008080", "#793eb4", "#A8B8BD")


class9 <- c("Lipid", "Amino acid", "Xenobiotics", "Peptide", "Carbohydrate",
            "Cofactors and vitamins", "Nucleotide", "Energy", "Others")


# met527 includes the abundance of each metabolite in all samples, as well the annotation of each metabolite (class, subclass)

# calculate mean, sd, se for subclass
met527_class <- Rmisc::summarySE(met527[,c("num_log2","class","subclass_simple")], measurevar="num_log2", 
                                 groupvars=c("class", "subclass_simple"),na.rm = T)

# calculate mean for each metabolite
met527_meancv <- Rmisc::summarySE(met527[,c("num_log2","metabolite")], measurevar="num_log2", groupvars=c("metabolite"),na.rm = T) %>% 
  left_join( met527_anno, by = "metabolite")


# only visualize main sub classes
da <- met527_class[met527_class$subclass_simple != "Others",]

rownames(da) <- 1:nrow(da)
da <- da[c(2,3,4,7,8,5,1,6,9:41),]
rownames(da) <- 1:nrow(da)
da["NA", ] <- NA

da$subclass_simple <- factor(da$subclass_simple, levels = da$subclass_simple[1:nrow(da)] )

da$subclass_number <- 1:nrow(da)

# data for visualization, angle of text
angle <-  90 - 360 * ( as.numeric(da$subclass_number) -0.5 ) /length(da$subclass_number)
da$hjust <- ifelse( angle < -90, 1, 0)
da$angle <- ifelse(angle < -90, angle+180, angle)


# prepare data for dotplot
dapoint <- met527_meancv[met527_meancv$subclass_simple != "Others",]
dapoint$subclass_simple <- factor(dapoint$subclass_simple, levels = da$subclass_simple)


coordy <- tibble('coordylocation' = seq(from = 0, to = max(dapoint$num_log2), 5),
                 'coordytext' = as.character(round(coordylocation, 2)),
                 'x' = 0.15)

griddata <- expand.grid('locationx' = seq(0,25,5), 'locationy' = coordy$coordylocation)



da$label <- as.character(da$subclass_simple)


# ----------- visualization -----------

f1a_circle_met527 <- ggplot()+
  # barplot
  geom_bar(da=da, aes(x=subclass_simple, y = num_log2,fill = class),stat="identity",position="dodge")+
  geom_errorbar(da= da, aes(x=subclass_simple, ymax=num_log2+sd,ymin=num_log2-sd),
                width=0.35,linewidth = 0.75, color = "gray40")+
  # dotplot
  geom_jitter(da=dapoint, aes(x = subclass_simple, y = num_log2, colour = class), 
              alpha = 0.9, position = position_jitter(0.2)) + 
  # color theme
  scale_fill_manual( values = color9_c, limits = class9)+
  scale_color_manual( values = color9dot_c, limits = class9)+
  
  # y-axis
  geom_vline(xintercept = 0.5, color = "gray40")+
  geom_segment(data = griddata, 
               aes(x = 0.4, xend = 0.5, y = locationy, yend = locationy),
               colour = "gray40", alpha=1, size=0.75) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="gray30", size=3.5 , angle=0, fontface="bold") + 
  labs(x=element_blank(), title = element_blank(), y= element_blank(), fill = "Sub classes") +
  scale_y_continuous(expand = c(0,0), limits = c(-5,25))+  
  
  # label for subclass
  geom_text(da=da, aes(x=subclass_simple, y=1, label=label, hjust="outward", angle = angle), 
            color="black", size=3.5,fontface="bold",
            inherit.aes = FALSE )+
  
  # theme
  coord_polar(start = 0)+
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        plot.title=element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        axis.title.x =element_blank(), 
        axis.title.y =element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y =element_blank(),
        legend.text=element_blank(), 
        legend.title=element_blank()) +
  theme(plot.margin=unit(rep(0,4),'cm'))
