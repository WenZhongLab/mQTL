# Netowrk of the 200 top metabolite-proiein associations ####
#Color palette
color_palette <- data.frame(group=c("Genetics", "Clinical measurements", "Lifestyle", "Visit", "Sex", "Age",
                                    "Body composition", "Lipid profile", "Erythrocytes","Leukocytes", "Platelets", "Glucose homeostasis", 
                                    "Blood pressure", "Urate", "Kidney", "Liver", "Heart", "Acute phase",
                                    "Housing_change", "NSAID_painmed", "Antibiotics_med", "Smoking"),
                            color=c("#35b779",	"#0073C2FF",	"#EFC000FF",	"#7f7f7f",	"#A30059",	"#8c564b",	
                                    "#4a6fe3", "#aec7e8",	"#FB9A99", "#FDE8DC",	"#dbdb8dff" ,	"#E377C2FF",
                                    "#ef9708",	"#968d6aff",	"#98df8a",	"#aa87deff",	"#f0b98d",	"#1CE6FF",
                                    "gold", "khaki", "goldenrod", "sandybrown"))


color_palette <- rbind(color_palette, 
                       data.frame(group=c("enzyme", "cytokine", "immunity", "hormone"),
                                  color=c("#ff80b2ff", "#ffab00ff", "#6f7c91ff", "#2ad4ffff")))


# Functions to make networks
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

layout_in_circles <- function(g, group=1, inner.scale=1, outer.scale=1) {
  layout <- lapply(split(V(g), group), function(x) {
    layout_in_circle(induced_subgraph(g,x))
  })
  layout <- Map(`*`, layout, seq_along(layout))
  layout[[1]] <- layout[[1]]*inner.scale #make the internal circle bigger (labels won't fit)
  layout[[2]] <- layout[[2]]*outer.scale #make the outer circle bigger (labels won't fit)
  
  lab.locs.inner <- radian.rescale(x=1:nrow(layout[[1]]), direction=-1, start=0)
  lab.locs.outer <- radian.rescale(x=1:nrow(layout[[2]]), direction=-1, start=0)
  lab.locs <- c(lab.locs.inner,lab.locs.outer)
  
  x <- matrix(0, nrow=vcount(g), ncol=2)
  split(x, group) <- layout
  return(list(layout=x, lab.locs=lab.locs))
}

# Create metabolites and proteins groups ####

#full graph
lmm_pm$metabolite.id <- anno_metabolite$name_for_annotation[match(lmm_pm$metabolite, anno_metabolite$metabolite)]
lmm_pm$protein.id <- protein_anno$symbol_id[match(lmm_pm$protein, protein_anno$Assay_now)]

full_link <- lmm_pm[, c("protein.id", "metabolite.id")]
full_node <- data.frame(node=unique(c(full_link$protein.id, full_link$metabolite.id)))
lmm_net_full <- graph_from_data_frame(d=full_link, vertices=full_node, directed=F) %>% igraph::simplify()

#comm.full <- cluster_fast_greedy(lmm_net_full)



# Metabolite-protein pairs from lmm analysis ####
Niter <- 200 # n. interactions to select
protein_metabolite <- lmm_pm[1:Niter,c("protein", "metabolite")] %>% unique()
protein_metabolite$metabolite.id <- anno_metabolite$name_for_annotation[match(protein_metabolite$metabolite, anno_metabolite$metabolite)]
protein_metabolite$protein.id <- protein_anno$symbol_id[match(protein_metabolite$protein, protein_anno$Assay_now)]
protein_metabolite <- protein_metabolite[!is.na(protein_metabolite$metabolite.id),]
protein_metabolite$direction <- ""
for (n in 1:nrow(protein_metabolite)){
  idx <- which(lmm_pm$protein==protein_metabolite$protein[n] & lmm_pm$metabolite==protein_metabolite$metabolite[n])
  stopifnot(length(idx)==1)
  b <- lmm_pm$estimate[idx]
  if (b>=0){
    protein_metabolite$direction[n] <- "Positive"
  } else {
    protein_metabolite$direction[n] <- "Negative"
  }
}


# LMM network ####

#interactions
net <- protein_metabolite[,c("protein.id", "metabolite.id", "direction")]

#add if edge is in MR too
net$MR <- F
for (n in 1:nrow(net)){
  if (any(MR.table$protein.id==net$protein.id[n] & MR.table$metabolite.id==net$metabolite.id[n])){
    net$MR[n] <- T
  }
}

#proteins groups
proteins_enzyme <- protein_class$Symbol[protein_class$Enzymes_HPA==T] %>% setdiff(NA) %>% intersect(net$protein)
proteins_growth <- protein_class$Symbol[grep("Growth|Hormone", protein_class$Molecular_function_HPA)] %>% intersect(net$protein)
proteins_cytokine <- protein_class$Symbol[grep("Cytokine", protein_class$Molecular_function_HPA)] %>% 
  union(protein_class$Symbol[grep("Chemotaxis", protein_class$Biological_process_HPA)]) %>% 
  intersect(net$protein)
proteins_immunity <- protein_class$Symbol[grep("Immunity|virus|Antiviral|B-cell", protein_class$Biological_process_HPA)] %>% 
  intersect(net$protein)


#create annotation
df <- data.frame(vertex=c(unique(net$protein.id), unique(net$metabolite.id)))
df$group <- "1"

#metabolites
df$group <- clinical_associations$Biomarker_group[match(df$vertex, clinical_associations$id)]
for (n in 1:nrow(df)){
  if (df$vertex[n] %in% net$metabolite.id){
    idx <- which(clinical_associations$name_for_annotation_new == df$vertex[n])
    if (length(idx)>0){
     idx.min <- idx[which.min(clinical_associations$p.value[idx])] 
      
     df$group[n] <- clinical_associations$Biomarker_group[idx.min]
     if (df$group[n]=="Lifestyle"){
       df$group[n] <- clinical_associations$clinical[idx.min]
     }
    }
  }
}

#proteins
df$group[df$vertex %in% proteins_enzyme] <- "enzyme"
df$group[df$vertex %in% proteins_cytokine] <- "cytokine"
df$group[df$vertex %in% proteins_growth] <- "hormone"
df$group[df$vertex %in% proteins_immunity] <- "immunity"

df$group[is.na(df$group)] <- "none"

#create network with endowed annotation
lmm_net <- graph_from_data_frame(net, vertices=df)

#Assign vertex to color groups
V(lmm_net)$color <- color_palette$color[match(V(lmm_net)$group, color_palette$group)]
V(lmm_net)$color[is.na(V(lmm_net)$color)] <- "#e7e7e7ff"

#color edges
net$direction[net$direction=="Negative"] <- rgb(255/256,102/256,102/256,100/256)
net$direction[net$direction=="Positive"] <- rgb(0/256,204/256,51/256,100/256)
E(lmm_net)$color <- net$direction
E(lmm_net)$lty <- ifelse(net$MR, 1,5)
E(lmm_net)$width <- ifelse(net$MR, 2,1)


#change order based on group
df$group <- gsub("Smoking", "AAASmoking", df$group) #to have lifestyle categories close
df$group <- gsub("Antibiotics_med", "AAAntibiotics_med", df$group) #to have lifestyle categories close
df$group <- gsub("NSAID_painmed", "AANSAID_painmed", df$group) #to have lifestyle categories close
df$group <- gsub("Housing_change", "AAHousing_change", df$group) #to have lifestyle categories close
group.uniq <- unique(df$group) %>% sort()
vertex.order <- character(0)
for (n in 1:length(group.uniq)){
  vertex.order <- c(vertex.order, df$vertex[which(df$group == group.uniq[n])])
}
df$group  <- gsub("AA", "", df$group)
group.uniq <- gsub("AA", "", group.uniq)
lmm_net <- permute(lmm_net, match(V(lmm_net)$name, vertex.order))

#create double circle layout
group <- V(lmm_net)$name %in% protein_metabolite$metabolite.id
layout <- layout_in_circles(lmm_net, group=group)
shape <- ifelse( names(V(lmm_net)) %in% net$protein.id, "circle", "square")
plot(lmm_net, layout=layout$layout, vertex.label.dist=2,
     vertex.size=4, vertex.label.cex=1, vertex.label.degree=layout$lab.locs,
     vertex.lab.color="black", edge.arrow.size=0, vertex.shape=shape,
     vertex.label="", alpha=0.1)

x <- layout$layout[,1]*.54
y <- layout$layout[,2]*.54

#create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(layout$layout[,1]/layout$layout[,2]))*(180/pi) < 0,  
               90 + atan(-(layout$layout[,1]/layout$layout[,2]))*(180/pi), 
               270 + atan(-layout$layout[,1]/layout$layout[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(lmm_net)$name[i], 
       adj=as.numeric(x[i]<0), pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T, family="Arial")
}

# MR network ####

#interactions
protein_metabolite <- lmm_pm[,c("protein", "metabolite")] 
protein_metabolite$direction <- ""
for (n in 1:nrow(protein_metabolite)){
  idx <- which(lmm_pm$protein==protein_metabolite$protein[n] & lmm_pm$metabolite==protein_metabolite$metabolite[n])
  stopifnot(length(idx)==1)
  b <- lmm_pm$estimate[idx]
  if (b>=0){
    protein_metabolite$direction[n] <- "Positive"
  } else {
    protein_metabolite$direction[n] <- "Negative"
  }
}


#add if edge is in MR too
protein_metabolite$MR <- F
for (n in 1:nrow(protein_metabolite)){
  if (any(MR.table$exposure==protein_metabolite$protein[n] & MR.table$outcome==protein_metabolite$metabolite[n])){
    protein_metabolite$MR[n] <- T
  }
}

protein_metabolite <- protein_metabolite[which(protein_metabolite$MR==T), ]
net <- protein_metabolite[,c("protein", "metabolite", "direction")]
net$metabolite <- anno_metabolite$name_for_annotation[match(net$metabolite, anno_metabolite$metabolite)]

#net <- MR.table[c("exposure", "id.outcome", "lmm_beta")] %>% unique()
#colnames(net) <- c("protein", "metabolite", "direction")

#proteins groups
proteins_enzyme <- protein_class$Symbol[protein_class$Enzymes_HPA==T] %>% setdiff(NA) %>% intersect(net$protein)
proteins_growth <- protein_class$Symbol[grep("Growth|Hormone", protein_class$Molecular_function_HPA)] %>% intersect(net$protein)
proteins_cytokine <- protein_class$Symbol[grep("Cytokine", protein_class$Molecular_function_HPA)] %>% 
  union(protein_class$Symbol[grep("Chemotaxis", protein_class$Biological_process_HPA)]) %>% 
  intersect(net$protein)
proteins_immunity <- protein_class$Symbol[grep("Immunity|virus|Antiviral|B-cell", protein_class$Biological_process_HPA)] %>% 
  intersect(net$protein)


#create annotation
df <- data.frame(vertex=c(unique(net$protein), unique(net$metabolite)))
df$group <- "1"

#metabolites
df$group <- clinical_associations$Biomarker_group[match(df$vertex, clinical_associations$id)]
for (n in 1:nrow(df)){
  if (df$vertex[n] %in% net$metabolite){
    idx <- which(clinical_associations$name_for_annotation_new == df$vertex[n])
    if (length(idx)>0){
      idx.min <- idx[which.min(clinical_associations$p.value[idx])] 
      
      df$group[n] <- clinical_associations$Biomarker_group[idx.min]
      if (df$group[n]=="Lifestyle"){
        df$group[n] <- clinical_associations$clinical[idx.min]
      }
    }
  }
}

#proteins
df$group[df$vertex %in% proteins_enzyme] <- "enzyme"
df$group[df$vertex %in% proteins_cytokine] <- "cytokine"
df$group[df$vertex %in% proteins_growth] <- "hormone"
df$group[df$vertex %in% proteins_immunity] <- "immunity"

df$group[is.na(df$group)] <- "none"

#create network with endowed annotation
MR_net <- graph_from_data_frame(net, vertices=df)

#Assign vertex to color groups
V(MR_net)$color <- color_palette$color[match(V(MR_net)$group, color_palette$group)]
V(MR_net)$color[is.na(V(MR_net)$color)] <- "#e7e7e7ff"

#color edges
net$direction[net$direction=="Negative"] <- rgb(255/256,102/256,102/256,100/256)
net$direction[net$direction=="Positive"] <- rgb(0/256,204/256,51/256,100/256)
E(MR_net)$color <- net$direction


#change order based on group
df$group <- gsub("Antibiotics_med", "AAAntibiotics_med", df$group) #to have lifestyle categories close
df$group <- gsub("NSAID_painmed", "AANSAID_painmed", df$group) #to have lifestyle categories close
df$group <- gsub("Housing_change", "AAHousing_change", df$group) #to have lifestyle categories close
group.uniq <- unique(df$group) %>% sort()
vertex.order <- character(0)
for (n in 1:length(group.uniq)){
  vertex.order <- c(vertex.order, df$vertex[which(df$group == group.uniq[n])])
}
df$group  <- gsub("AA", "", df$group)
group.uniq <- gsub("AA", "", group.uniq)
MR_net <- permute(MR_net, match(V(MR_net)$name, vertex.order))

#create double circle layout
group <- V(MR_net)$name %in% net$metabolite
layout <- layout_in_circles(MR_net, group=group)
shape <- ifelse( names(V(MR_net)) %in% net$protein, "circle", "square")
plot(MR_net, layout=layout$layout, vertex.label.dist=2,
     vertex.size=4, vertex.label.cex=1, vertex.label.degree=layout$lab.locs,
     vertex.lab.color="black", edge.arrow.size=0, vertex.shape=shape,
     vertex.label="", alpha=0.1)

x <- layout$layout[,1]*.54
y <- layout$layout[,2]*.54

#create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(layout$layout[,1]/layout$layout[,2]))*(180/pi) < 0,  
               90 + atan(-(layout$layout[,1]/layout$layout[,2]))*(180/pi), 
               270 + atan(-layout$layout[,1]/layout$layout[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(MR_net)$name[i], 
       adj=as.numeric(x[i]<0), pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T, family="Arial")
}
