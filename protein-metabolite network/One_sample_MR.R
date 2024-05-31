###### One sample MR
##### by wangjing 2023.11.09

## Setup

## Load packages

# Tidyverse
library(tidyverse)

# data.table
library(data.table)

## MR
library(AER)


## workdir

setwd("/proj/sens2022015/jingw/wellness_20221101/OneSample_MR_v3/olink_qPCR")


args=commandArgs(T)


p = print(args)


### Load anno data

df_anno_m =  read_tsv("../QCed_metabolites_name.txt") %>%
  select(outcome = metabolite, id.outcome = name_for_annotation_new)

df_anno_p =  read_tsv("protein_gene_anno.txt") %>%
  select(exposure = Assay_now, id.exposure = protein_name)

metabolite = df_anno_m %>% dplyr::select(outcome) %>% pull

### load  cov data
df_cov = read_tsv("../HPA_101_cov.txt") %>%
  select(-FID)

### indep_condi_snp_list
df_indep_condi = read_delim("qPCR_pQTL_all_indep_5e-08_condi.txt")  %>% select (label, Assay, SNP) %>% distinct()


############### label loop ############

  
print(paste0("############### start  loop ---- ",p, " #################"))
Sys.time()

### output file
output = data.frame()

### protein
df_p = read_tsv("qPCR_794_median.txt") %>%
  select(IID, !!sym(p))%>%
  distinct()



### independent & conditinal snp
df_snp = df_indep_condi %>%
  filter(label == l) %>%
  filter(Assay == p) %>%
  dplyr::select(Assay_now = Assay, SNP) %>%
  distinct()

snp_number = nrow(df_snp)

if (snp_number > 0) {
  
  ############### get genotype ##############
  df_geno = read_delim("median.clumping.snp.p5e-8.geno.txt", delim = " ", col_names = F) %>%
    rename(SNP = X1, IID = X2, genotype = X3, geno = X4) %>%
    inner_join(df_snp, by = "SNP") %>%
    mutate(SNP = str_replace_all(SNP, ":", "_")) %>%
    mutate(SNP = paste0("chr", SNP)) %>%
    dplyr::select(IID, SNP, geno) %>%
    pivot_wider(id_cols = "IID", names_from = "SNP", values_from = "geno")
  
  ################  metabolite loop ##############
  for (m in metabolite){
    
    df_m = read_tsv(paste0("../mQTL.coeff.txt")) %>%
      dplyr::select(IID, !!sym(m))%>%
      distinct()
    
    ############### data for MR ##############
    df_all = df_p %>%
      full_join(df_m, by = "IID") %>%
      full_join(df_cov, by = "IID") %>%
      full_join(df_geno, by = "IID") %>%
      distinct()
    
    ########### genorate formula ###############
    tmp = colnames(df_all)
    
    exposure = tmp[2] ## e.g. CCL15
    
    outcome = tmp[3] ## e.g. PG_81
    
    snp = str_c(tmp[-1:-7], collapse = " + ") ## e.g. chr17_36001607_C_CAGGGCAG + chr8_7038210_A_G
    
    cov = c("Sex + Age_at_Visit + BMI")
    
    part_front = paste(exposure, cov, sep =" + ") ## e.g. CCL15 + Sex + Age_at_Visit + BMI
    
    part_back = paste(snp, cov, sep =" + ") ## e.g. chr17_36001607_C_CAGGGCAG + chr8_7038210_A_G + Sex + Age_at_Visit + BMI
    
    merge = paste(part_front, part_back, sep = " | ")
    
    ############ one-sampe MR analysis ############
    
    ### e.g. PG_81 ~ CCL15 + Sex + Age_at_Visit + BMI | chr17_36001607_C_CAGGGCAG + chr8_7038210_A_G + Sex + Age_at_Visit + BMI
    
    fit = ivreg(formula(paste(outcome, merge, sep = " ~ ")), data = df_all)
    
    
    ############### results ############
    df_res = summary(fit,diagnostics = TRUE)
    
    
    
    df_beta = df_res[["coefficients"]] %>% as.data.frame() %>%
      rownames_to_column(., var = "exposure") %>%
      filter(exposure == p) %>% ##### protein name
      rename(beta = Estimate, se = `Std. Error`, t_value = `t value`, MR_P = `Pr(>|t|)`) %>%
      mutate(low_ci = beta - 1.96*se,
             up_ci = beta + 1.96*se) %>%
      select(exposure, beta, se, low_ci, up_ci, t_value, MR_P)
    
    #df_beta
    
    
    df_diagnostics = df_res[["diagnostics"]] %>% t() %>% as.data.frame() %>%
      rownames_to_column(., var = "tmp") %>%
      mutate(exposure = p)  %>% ##### protein name
      filter(tmp == "p-value") %>%
      dplyr::select(exposure, Weak_instruments_P = `Weak instruments`, Wu_Hausman_P = Sargan, Sargan_P = Sargan)
    
    #df_diagnostics
    
    
    df_r2 = df_res[["r.squared"]] %>% as.data.frame() %>% rename(r2 = ".") %>%
      mutate(exposure = p)  %>% ##### protein name
      dplyr::select(exposure, r2)
    
    #df_r2
    
    
    df_r2_adj =  df_res[["adj.r.squared"]] %>% as.data.frame() %>% rename(r2_adj = ".") %>%
      mutate(exposure = p)  %>% ##### protein name
      dplyr::select(exposure, r2_adj)
    
    #df_r2_adj
    
    df_wald_test = df_res[["waldtest"]] %>% t() %>% as.data.frame() %>%
      mutate(exposure = p)  %>% ##### protein name
      dplyr::select(exposure, Wald_test_statistic = V1, Wald_test_P = V2)
    
    #df_wald_test
    
    
    ############### combine MR results ############
    df_res_all = df_beta %>%
      mutate(label = "soft_threshold") %>%
      mutate(outcome = m) %>% ### outcome name
      left_join(df_anno_m, by = "outcome") %>%
      left_join(df_anno_p, by = "exposure") %>%
      left_join(df_diagnostics, by = "exposure") %>%
      left_join(df_r2, by = "exposure") %>%
      left_join(df_r2_adj, by = "exposure") %>%
      left_join(df_wald_test, by = "exposure") %>%
      select(label, exposure, outcome, id.exposure, id.outcome, everything())
    
    
    print("################### head(df_res_all) #####################")
    
    head(df_res_all)
    
    
    
    ############### merge all the results ############
    
    output = rbind(output, df_res_all)
    
    
  }
  
  write.table(output, file = paste0("results/",p,"_one_sample_mr_result_qPCR_soft.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
  
  
} else {
  
  print(paste0("############## protein: ", p, " no siginficant SNP ##############"))
}
Sys.time()


##### load all of results and calculate the adjusted P
df_soft = read_tsv(here::here("data", DOCNAME, "all_one_sample_mr_result_qPCR_soft.txt"))%>%
  mutate(MR_P_BHadj_qPCR = p.adjust(MR_P_qPCR, method = "BH")) %>% ##### BH adj
  select(p_m_id,exposure:MR_P_qPCR, MR_P_BHadj_qPCR, MR_Weak_instruments_P_qPCR:MR_Wald_test_P_qPCR,lmm_beta:lmm_p_adj) %>%
  arrange(lmm_p_adj)


df_soft %>%
  write.table(file = here::here("output", DOCNAME, "one_sample_mr.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
