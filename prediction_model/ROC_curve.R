
# plot ROC curve based on true class and probability (group = 1)

ROCplot <- function(filename, color_need = "red"){
  library(pROC)
  library(tidyverse)
  data <- read.csv(filename, header = F)
  y_true <- c(t(data[1,]))
  y_pred <- c(t(data[2,]))
  rocobj <- plot.roc( y_true,y_pred,
                      main = "Confidence intervals", 
                      percent=TRUE,
                      ci = TRUE,                 
                      print.auc = TRUE)          
  ciobj <- ci.se(rocobj, conf.level=0.95,                        
                 specificities = seq(0, 100, 5)) 
  
  rocres = pROC::roc(y_true, y_pred)
  roxrange <- ci(rocres)
  str = c("AUC: ",round(as.numeric(substr(rocres$auc, 1,7)),2)," (", 
          round(roxrange[1],2),"-", round(roxrange[3],2),")"  )
  note = paste(str, collapse="")
  #roc
  ciobj_df <- data.frame(ciobj)
  ciobj_df$x <- 1 - as.numeric( substr( rownames(ciobj_df), 1, nchar(rownames(ciobj_df))-1 ) )/100
  colnames(ciobj_df)[1:3] <- c("y1", "y2","y3")
  ciobj_df$y1 <- ciobj_df$y1/100
  ciobj_df$y2 <- ciobj_df$y3/100
  ciobj_df$y3 <- ciobj_df$y3/100
  
  p <- ggroc(rocres, col = color_need,legacy.axes = TRUE, size = 1.5) +
    geom_ribbon(aes(x = x,ymax = y3, ymin = y1 ),data = ciobj_df, fill = "gray", alpha = 0.5) +
    scale_y_continuous(expand=c(0,0), breaks = c(0,0.25,0.5,0.75,1),labels = c(0,"",0.5,"",1))+
    scale_x_continuous(expand=c(0,0), breaks = c(0,0.25,0.5,0.75,1),labels = c(0,"",0.5,"",1))+
    geom_abline(intercept=0,slope=1,linetype = "dashed" )+
    labs(x = "False Positive Rate", y = "True Positive Rate") +
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey",
                                          size = 0.25 ),
          panel.border = element_blank(), 
          axis.line = element_line(size=0.5, colour = "black"),
          legend.position = "none",
          axis.title.x =element_text(size=20,color="black"), 
          axis.title.y =element_text(size=20,color="black"),
          axis.text.x =element_text(size=20,color="black"), 
          axis.text.y =element_text(size=20,color="black"),
          plot.margin=unit(rep(0.5,4),'cm'))
  return(p)
  
  
}



# run the function

ROCplot("./result/tier1protein_100balance_pro13scale_forROC.csv", "#68775e")
ggsave("tier1protein_ROC.pdf", width = 4.5, height = 5.25)

