# Tidyverse
library(tidyverse)

# data.table
library(data.table)

library(broom)
# plot
library(patchwork)
library(ggpubr)
library(cowplot)
library(ggtern)
library(plotly)
theme_set(theme_cowplot())

###### load umap results
df_umap = read_delim(here::here("data", DOCNAME, "UMAP.LMM.nolipid_20240517.txt")) %>%
  mutate(iid = paste0(subject, "_", visit)) %>%
  select(iid, UMAP1, UMAP2) %>%
  pivot_longer(cols = -iid, names_to = 'umap', values_to = 'num') 

### load annotation information
df_anno = read_tsv(here::here("data/Qced", "QCed_metabolites.txt")) %>%
  dplyr::select(fid, iid, sample, IID = geno_id_2, visit) %>%
  distinct()

### load clinical data
df_cli = read_tsv(here::here("data/Qced", "QCed_wellness_clinical_alt.txt")) %>%
  filter(visit < 5) %>%
  filter(!is.na(Biomarker_group)) %>%
  dplyr::select(iid, clinical, num) %>%
  pivot_wider(id_cols = "iid", names_from = "clinical", values_from = "num")

### load lifestyle data
df_lifestyle = read_delim(here::here("output", "11.lifestyle_v3", "wellness_lifestyle.txt"))

df_VD = df_VD = read_tsv(here::here("data", "11.lifestyle_v3", "VitD")) %>%
  mutate(fid = paste0(subject_id, "_", Visit)) %>%
  select(fid, VitD)

### load mQTL genotype
df_mQTL =  read_tsv(here::here("output","03.GWAS_v7_conditional_analysis",  "coeff_mQTL_2.17-09_independent_condi_geneName.txt")) %>%
  filter(OBS_CT > 50) %>% ### remove lipid metabolites
  select(SNP) %>%
  distinct()

df_mQTL_geno = read_delim(here::here("data", "03.GWAS_v7", "coeff_snp_5e-08_geno.txt"), delim = " ") %>%
  mutate(ref = str_split(SNP, ':', simplify = T)[,3],
           alt = str_split(SNP, ':', simplify = T)[,4]) %>%
  mutate(geno = case_when(
    Genotype == paste0(ref,ref) ~ 0,
    Genotype == paste0(ref,alt) ~ 1,
    Genotype == paste0(alt,ref) ~ 1,
    Genotype == paste0(alt,alt) ~ 2)) %>%
  select(IID, SNP, geno) %>%
  inner_join(df_mQTL,  by = "SNP")

#### load pQTL  data
df_pQTL = read_tsv(here::here("data", DOCNAME, "qPCR_pQTL_all_indep_5e-08_condi.txt")) %>%
  filter(label == "median") %>%
  filter(P < 6e-11) %>%
  select(SNP) %>%
  distinct()

### pQTL genotype data
df_pQTL_geno = read_delim(here::here("data", "05.one_sample_MR_v3", "pQTL_geno.txt"), col_names = T, delim = " ") %>%
  mutate(ref = str_split(SNP, ':', simplify = T)[,3],
           alt = str_split(SNP, ':', simplify = T)[,4]) %>%
  mutate(geno = case_when(
    genotype == paste0(ref,ref) ~ 0,
    genotype == paste0(ref,alt) ~ 1,
    genotype == paste0(alt,ref) ~ 1,
    genotype == paste0(alt,alt) ~ 2)) %>%
  select(IID, SNP, geno) %>%
  inner_join(df_pQTL,  by = "SNP")

### all genotype data
df_snp = rbind(df_mQTL_geno, df_pQTL_geno) %>% distinct() %>%
  mutate(SNP = paste0("chr",SNP)) %>%
  select(IID, SNP, geno) %>%
  pivot_wider(id_cols = "IID", names_from = "SNP", values_from = "geno")


######################################## linear model #######################################

variables = c("UMAP1", "UMAP2")

res_all = data.frame()
res_var_all = data.frame()

for (i in variables){
  
  
  df_all = df_umap %>%
    filter(umap == i) %>%
    inner_join(df_anno, by = "iid") %>%
    inner_join(df_cli, by = "iid") %>%
    inner_join(df_lifestyle, by = "iid") %>%
    inner_join(df_VD, by = "fid") %>%
    inner_join(df_snp, by = "IID")   
    
  ### data for lm  
  dat = df_all %>%
    select(-fid, -iid, -sample, -IID, -umap)
  
  
  ### linear model
  lme = lm(num ~ ., data = dat) ## multiple linear model
  
  
  ### statistical results
  res = tidy(lme) %>%
    filter(term != "(Intercept)") %>%
    mutate(umap = i) %>%
    dplyr::select(umap, term:p.value)
    
  
  ### variance
  res_var = tidy(anova(lme)) %>%   ###variance analysis
    mutate(sum = sum(sumsq)) %>%
    mutate(R2 = sumsq/sum) %>%
    mutate(umap = i) %>%
    select(umap, term, R2) %>%
    filter(term != "Residuals")
  
  
  res_all = rbind(res_all, res)
  res_var_all = rbind(res_var_all, res_var)
    
}

res_all_2 = res_all %>% mutate(id = paste0(umap, "_", term)) %>% mutate(p_adj = p.adjust(`p.value`, method = "BH"))
res_var_all_2 = res_var_all %>% mutate(id = paste0(umap, "_", term)) %>% select(id, R2)

df_lm_res = res_all_2 %>%
  left_join(res_var_all_2, by = "id") %>%
  select(-id)

df_lm_res %>%
   write.table(file = here::here("output", DOCNAME, "umap_lm_lifestyle_clinical_genotype_results_variance.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
