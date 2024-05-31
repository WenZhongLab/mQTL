## setup
library(lme4) 
library(lmerTest) 
library(broom.mixed) 
library(broom)
library(skimr)

## Load data


print("############load data ###################")

df_p_m = read_delim("combined_qPCR_proteins_metabolites.txt")  



#### perform linear mixed model regression #######
print("#### perform linear mixed model regression #######")


args=commandArgs(T)

m = print(args)


df_res = data.frame()

linear <- function(x) {
  ### x: protein 
  x = df_p_m[,x] %>% as.data.frame() 
  .var = paste0(colnames(x), ":", m) 
  x = x[,1] 
  
  ### y: metabolite
  y = df_p_m[,m] %>% as.data.frame() 
  y = y[,1]
  
  ##### model: metabolite ~ protein + sex + age + BMI + (1 | visit) + (1 | sample)
  lme = lmerTest::lmer(y ~ x + Sex + Age_at_Visit + BMI + (1 | visit) + (1 | sample), data = df_p_m) ## linear mixed model
  
  res = tidy(lme)
  
  res[2, 3] = .var 
  res[2, 3:8]  
  
  
}

resList = lapply(2:795, FUN = linear)
res = do.call(rbind, resList) # combine all the results

df_res = rbind(df_res, res)


print("############ output data  ###################")

write.table(df_res, file = paste0("results/", m, ".lmm.res.txt"), col.names = T, row.names = F, sep = '\t', quote = F)
