
# ---- canonical correspondence analysis ----------------
# corresponding to figure 3a 

library(vegan)
library(tidyverse)
library(ggnewscale)
library(ggrepel)
library(ggsci)

# da1: rows-sample, columns-metabolite
# da2: rows-samples, columns-clinical variables
# sampleid187: sample id of these without missing

da2$Gender <- factor(da2$Gender, levels = c("m","f"))

set.seed(12)
met527_cli_cca_scale_new <- vegan::cca( X=da1[sampleid187 ,],  Y=da2 , scale = T)

CCA_res_new <- summary(met527_cli_cca_scale_new)

# ---------------- visualization: CCA --------------------

# prepare data for visualization

met_cca_new <- data.frame(CCA_res_new$species)*50 %>% 
  rownames_to_column(var="metabolite") %>% 
  left_join( met527_anno, by = "metabolite")
colnames(met_cca_new)[8] <- 'name_for_annotation'

sample_cca_new <- data.frame(CCA_res_new$sites) %>% 
  rownames_to_column(var="iid") %>% 
  left_join(sampleinfo[sampleid187,], by = "iid")

cli_cca_new <- data.frame(CCA_res_new$biplot)*5 %>% 
  rownames_to_column(var="Short.name") %>% 
  left_join( cliname, by = "Short.name") %>% 
  mutate(name_label = as.character(Short.name))

# ------------- visualization -------------
p <- ggplot()+
  # sample
  scale_color_manual( values = scale_color_palette_sex, breaks = c("f", "m"), labels = c("Female sample", "Male sample"))+
  geom_point(data = sample_cca_new, aes(x=CCA1, y=CCA2, color = Gender), alpha = 1, shape = 17, size = 2)+
  
  # metabolites
  new_scale_color()+
  geom_point(data = met_cca_new, aes(x=CCA1, y=CCA2), color = "gray20", alpha=0.5,shape = 15)+
  
  geom_point(data = met_cca_new[(met_cca_new$CCA1 > 2.25 | met_cca_new$CCA1 < -1.5 | met_cca_new$CCA2 > 2 | met_cca_new$CCA2 < -1.3 ) , ], aes(x=CCA1, y=CCA2), color = "gray20", alpha=0.5, size = 2,shape = 15)+
  geom_text_repel(data = met_cca_new[(met_cca_new$CCA1 > 2.25 | met_cca_new$CCA1 < -1.5 | met_cca_new$CCA2 > 2 | met_cca_new$CCA2 < -1.3 ) , ],
                  aes(x=CCA1, y=CCA2, label = name_for_annotation), min.segment.length= 0.25, 
                  color = "gray10", size = 4,fontface="plain" ) +
  # clinical measures
  new_scale_color()+ #c2c1d3
  scale_color_manual( values = append(append(c("#fc1658","#fc1658"),pal_d3("category20")(11)[-c(1,4,6,8,10)]),c("#a0bb65" )),
                      name = "Clinical measure groups")+
  geom_point(data = cli_cca_new, aes(x=CCA1, y=CCA2, color = Biomarker.group), alpha = 0.5) +
  geom_point(data = cli_cca_new[(cli_cca_new$CCA1 > 1 | cli_cca_new$CCA1 < -1.2 |
                                     cli_cca_new$CCA2 > 1 | cli_cca_new$CCA2 < -1 ) , ],
             aes(x=CCA1, y=CCA2, color = Biomarker.group), alpha = 1) +
  geom_text_repel(data = cli_cca_new[(cli_cca_new$CCA1 > 1 | cli_cca_new$CCA1 < -1.2 | 
                                          cli_cca_new$CCA2 > 1 | cli_cca_new$CCA2 < -1 ) , ], 
                  aes(x=CCA1, y=CCA2, label = name_label, color = Biomarker.group), size = 4) +
  
  # add arrow
  geom_segment(da = cli_cca_new[cli_cca_new$Short.name %in% c("BMI", "Gender_female"),], aes(x=0,y=0,xend=CCA1,yend=CCA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "#ef5c7a",lineend = "round", linewidth = 1) +
  
  scale_y_continuous( position = "right")+ 
  
  labs(x="CCA1", y= "CCA2", title = element_blank(), color = "Biomarker groups") +
  theme_bw() + theme(panel.grid=element_blank(),panel.background = element_blank(),
                     panel.border = element_rect(linewidth = 0.75, colour = "black"),) + 
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        legend.position = "none",
        axis.title.x =element_text(size=12, colour = "black"), 
        axis.title.y =element_text(size=12, colour = "black"),
        axis.text.x =element_text(size=12, colour = "black"), 
        axis.text.y =element_text(size=12, colour = "black"),
        legend.text=element_text(size=12, colour = "black"), 
        legend.title=element_text(size=12, colour = "black"))

