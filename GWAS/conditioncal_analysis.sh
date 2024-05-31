
``` bash
#!/bin/bash -l

#SBATCH -A sens2022015
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J coeff_lc_2_condi
#SBATCH -e coeff_lc_2_condi.err
#SBATCH -o coeff_lc_2_condi.out


module load bioinfo-tools
module load plink
module load plink2

cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/conditional_analysis

i="lc_2"


echo '################### coeff gwas for protein `$i`  ######################'


echo Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)


awk '{if($1 == "lc_2"){print $2}}' coeff_leader_SNP.txt > results/${i}_leader_SNP

awk '{if($1 == "lc_2"){print $2}}' coeff_all_SNP.txt > results/${i}_all_SNP


## re-calculated the P value, using leader SNP as convariate.


plink2 --bfile ../../HPA_101_qc_maf0.05 --extract results/${i}_all_SNP --condition-list results/${i}_leader_SNP --ci 0.95 --linear hide-covar --pheno ../mQTL.coeff.txt --pheno-name ${i} --covar ../mQTL.cov_log2.txt  --covar-name Sex, Age_b, BMI_b --covar-variance-standardize  --out results/condition

cd results

##### only the association with conditional P < 0.01 were considered to be independent mQTLs. 

awk 'BEGIN {FS="\t"}{OFS="\t"} {if($17 != "NA" && strtonum($17) < 0.01) {print "lc_2",$3}}' condition.${i}.glm.linear > coeff_${i}_codition_SNP.txt

chmod +x coeff_${i}_codition_SNP.txt


work_end_time=$(date +%s)
echo "#######-- cv gwas End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))

```
