## example of GWAS analysis


``` bash
#!/bin/bash -l

#SBATCH -A sens2022015
#SBATCH -p node
#SBATCH -n 2
#SBATCH -t 100:00:00
#SBATCH -J lc_2_gwas
#SBATCH -e lc_2_gwas.err
#SBATCH -o lc_2_gwas.out

module load bioinfo-tools plink plink2


cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7


i="lc_2"


echo "################### coeff gwas for metabolite $i  ######################"


echo Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)



plink2 --bfile ../HPA_101_qc_maf0.05  --ci 0.95 --linear hide-covar --pheno mQTL.coeff.txt --pheno-name ${i} --covar mQTL.cov_log2.txt  --covar-name Sex, Age_b, BMI_b   --covar-variance-standardize  --out coeff/1_all_results/coef


echo '################### GWAS filter: P < 1e-05 ######################'

awk 'BEGIN {FS="\t"}{OFS="\t"} {if($17 != "NA" && strtonum($17) < 0.00001) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$16,$17}}' coeff/1_all_results/coeff.${i}.glm.linear > coeff/2_p1e-05/coeff.${i}.glm.linear.p1e-5.txt


echo "################### clumping  ######################"


sed -i 's/ID/SNP/g'  coeff/2_p1e-05/coeff.${i}.glm.linear.p1e-5.txt

plink --bfile ../HPA_101_qc_maf0.05 --clump  coeff/2_p1e-05/coeff.${i}.glm.linear.p1e-5.txt --clump-p1 0.00001 --clump-p2 0.01  --clump-r2 0.1 --clump-kb 500 --out coeff/3_clumping/coeff.${i}.glm.linear.p1e-5

awk 'BEGIN {OFS="\t"} NR==FNR{a[$3]=$0;next} NR>FNR{if($3 in a) print $0}' coeff/3_clumping/coeff.${i}.glm.linear.p1e-5.clumped coeff/2_p1e-05/coeff.${i}.glm.linear.p1e-5.txt > coeff/3_clumping/coeff.${i}.glm.linear.p1e-5.clumping.txt


work_end_time=$(date +%s)
work_end_time=$(date +%s)
echo "#######-- coeff gwas End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))
echo "#######-- coeff gwas Elapsed Time:$elapsed_time sec --#######"
```
