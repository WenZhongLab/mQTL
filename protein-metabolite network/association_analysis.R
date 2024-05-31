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

###### For each metabolite, the above code was looped and then all the results were combined to calculate the corrected p-values ###############



### load results

list <- list.files(path="C:/Users/jinwa42/OneDrive - Linköpings universitet/02_project/2022_09_03_mQTL/mQTL_wfr/data/04.m_p_lmm/results",pattern = "lmm.res.txt")

df_lmm_res <- data.frame()

for(i in list){
  path <- i
  tmp <- read_delim(file = paste0("C:/Users/jinwa42/OneDrive - Linköpings universitet/02_project/2022_09_03_mQTL/mQTL_wfr/data/04.m_p_lmm/results/",path), delim = "\t", col_names = T) %>% mutate(metabolite = i)
    
  df_lmm_res <- rbind(df_lmm_res,tmp)
}


### add p_adj
df_lmm_res2 = df_lmm_res %>%
  rename(lmm_p = `p.value`) %>%
  mutate(lmm_p_adj = p.adjust(lmm_p, method = "fdr")) 

df_lmm_res2 %>%
  write.table(here::here("output", DOCNAME, "m_p_lmm_res.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
