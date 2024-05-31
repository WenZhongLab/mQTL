# Tidyverse
library(tidyverse)
# data.table
library(data.table)
library(broom)


##### metabolite table
df_m = read_tsv(here::here("data/Qced", "QCed_metabolites.txt")) %>%
  dplyr::select(iid, sample, IID = geno_id_2, num_log2, metabolite, visit)

#### clinical table
df_clinical = read_tsv(here::here("data/Qced", "QCed_wellness_clinical.txt"))

df_cli = read_tsv(here::here("data/Qced", "QCed_wellness_clinical_alt.txt")) %>%
  filter(visit < 5) %>%
  filter(!is.na(Biomarker_group)) %>%
  dplyr::select(iid, clinical, num) %>%
  pivot_wider(id_cols = "iid", names_from = "clinical", values_from = "num")

### GWAS results
df_condi =  read_tsv(here::here("output", "03.GWAS_v7_conditional_analysis", "coeff_mQTL_2.17-09_independent_condi_geneName.txt")) 

#### snp list
df_condi_snp = df_condi %>% select(SNP) %>% distinct()

#### genotype table
df_geno = read_delim(here::here("data","03.GWAS_v7", "coeff_snp_5e-08_geno.txt"),delim = " ") %>%
    mutate(ref = str_split(SNP, ':', simplify = T)[,3],
           alt = str_split(SNP, ':', simplify = T)[,4]) %>%
  mutate(geno = case_when(
    Genotype == paste0(ref,ref) ~ 0,
    Genotype == paste0(ref,alt) ~ 1,
    Genotype == paste0(alt,ref) ~ 1,
    Genotype == paste0(alt,alt) ~ 2)) %>%
  select(IID, SNP, geno) %>%
  inner_join(df_condi_snp,  by = "SNP")
  
##### lifestyle table
df_lifestyle = read_delim(here::here("output", DOCNAME, "wellness_lifestyle.txt"))

df_VD = read_tsv(here::here("data", DOCNAME, "VitD")) %>%
  mutate(fid = paste0(subject_id, "_", Visit)) %>%
  inner_join(df_clinical, by = "fid") %>%
  select(iid, VitD)


##### linear model

m = df_m %>%
  dplyr::select(metabolite) %>%
  distinct() %>%
  pull

#m = m[1:10]

res_all = data.frame()
res_var_all = data.frame()

for (i in m){
  
  ### merge data
  df_snp = df_condi %>%
    dplyr::select(metabolite, SNP) %>%
    distinct() %>%
    filter(metabolite == i) %>%
    select(SNP) %>%
    left_join(df_geno, by = "SNP") %>%
    mutate(SNP = paste0("chr",SNP)) %>%
    pivot_wider(id_cols = "IID", names_from = "SNP", values_from = "geno")
  
  
  df_all = df_m %>%
    filter(metabolite == i) %>%
    left_join(df_lifestyle, by = "iid") %>%
    left_join(df_cli, by = "iid") %>%
    left_join(df_snp, by = "IID")  %>%
    left_join(df_VD, by = "iid")
    
  ### data for lm  
  dat = df_all %>%
    dplyr::select(-iid, -sample, -IID,  -metabolite)
    
  ### linear model
  lme = lm(num_log2 ~ ., data = dat) ## multiple linear model
  
  ### statistical results
  res = tidy(lme) %>%
    filter(term != "(Intercept)") %>%
    mutate(metabolite = i) %>%
    dplyr::select(metabolite, term:p.value)
    
  
  ### variance
  res_var = tidy(anova(lme)) %>%   ### variance 
    mutate(sum = sum(sumsq)) %>%
    mutate(R2 = sumsq/sum) %>%
    mutate(metabolite = i) %>%
    select(metabolite, term, R2) %>%
    filter(term != "Residuals")
  
  
  res_all = rbind(res_all, res)
  res_var_all = rbind(res_var_all, res_var)
    
}


###### results for beta, se, P 
res_all_2 = res_all %>% mutate(id = paste0(metabolite, "_", term)) %>% mutate(p_adj = p.adjust(`p.value`, method = "BH"))
#### results for R2
res_var_all_2 = res_var_all %>% mutate(id = paste0(metabolite, "_", term)) %>% select(id, R2)

##### merge
df_lm_res = res_all_2 %>%
  left_join(res_var_all_2, by = "id") %>%
  select(-id)

df_lm_res %>%
write.table(file = here::here("output", DOCNAME, "lm_lifestyle_clinical_genotype_results.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

