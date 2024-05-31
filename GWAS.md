# GWAS pipeline

pwd: `/proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7`

.
├── coeff (gwas result)
│   ├── 1_all_results
│   ├── 1_pdf
│   │   ├── novel_mQTL
│   │   └── validated_mQTL
│   ├── 2_p1e-05
│   └── 3_clumping
├── coeff_adjGAB (gwas result)
│   ├── 1_all_results
│   ├── 2_p1e-05
│   └── 3_clumping
├── conditional_analysis (coeff)
│   ├── results
│   └── sbatch
├── extract_anno (coeff and coeff_adjGAB)
├── extract_geno (coeff)
└── sbatch (gwas.sh)

## GWAS 

test tow models:

- coeff ~ SNP + sex _age_b + BMI_b

- coeff_adjGAB ~ SNP

``` bash
for i in `cat ../mQTL.name.txt`;do cp zzzz_gwas.sh ${i}_gwas.sh;sed -i "s/zzzz/$i/g" ${i}_gwas.sh;sbatch ${i}_gwas.sh;done
```

### zzzz_gwas.sh

``` bash
#!/bin/bash -l

#SBATCH -A sens2022015
#SBATCH -p node
#SBATCH -n 2
#SBATCH -t 100:00:00
#SBATCH -J zzzz_gwas
#SBATCH -e zzzz_gwas.err
#SBATCH -o zzzz_gwas.out


module load bioinfo-tools plink plink2


cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7


i="zzzz"


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




echo "################### coeff_adjGAB gwas for metabolite $i  ######################"


echo Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)


cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7

plink2 --bfile ../HPA_101_qc_maf0.05  --ci 0.95 --linear allow-no-covars  --pheno mQTL.coeff_adjGAB.txt --pheno-name ${i}  --out coeff_adjGAB/1_all_results/coeff_adjGAB


echo "################### GWAS filter: P < 1e-05 ######################"

awk 'BEGIN {FS="\t"}{OFS="\t"} {if($17 != "NA" && strtonum($17) < 0.00001) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$16,$17}}' coeff_adjGAB/1_all_results/coeff_adjGAB.${i}.glm.linear > coeff_adjGAB/2_p1e-05/coeff_adjGAB.${i}


echo "################### clumping  ######################"


sed -i 's/ID/SNP/g'  coeff_adjGAB/2_p1e-05/coeff_adjGAB.${i}.glm.linear.p1e-5.txt

plink --bfile ../HPA_101_qc_maf0.05 --clump  coeff_adjGAB/2_p1e-05/coeff_adjGAB.${i}.glm.linear.p1e-5.txt --clump-p1 0.00001 --clump-p2 0.01  --clump-r2 0.1 --clump-kb 500 --out coeff_adjGAB/3_clumping/coeff_adjGAB.${i}.glm.line

awk 'BEGIN {OFS="\t"} NR==FNR{a[$3]=$0;next} NR>FNR{if($3 in a) print $0}' coeff_adjGAB/3_clumping/coeff_adjGAB.${i}.glm.linear.p1e-5.clumped coeff_adjGAB/2_p1e-05/coeff_adjGAB.${i}.glm.linear.p1e-5.txt > coeff_adjGAB/3_clumping/



work_end_time=$(date +%s)
echo "#######-- cv gwas End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))
echo "#######-- cv gwas Elapsed Time:$elapsed_time sec --#######"



echo '###################  process finished :)  ######################'

```



## extract anno

pwd: `/proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/extract_anno`

###  extract_anno.sh

``` bash
#!/bin/bash -l

awk 'BEGIN {FS="\t"}{OFS="\t"} NR==FNR{a[$1]=$0;next} NR>FNR{if($1 in a) print $0}' all.p1e-5.SNP.txt /proj/sens2022015/jingw/wellness_20221101/WGS_QC/anno_yh/HPA_101_qc_maf0.05.vep.txt > gwas_v7_all_p1e-5_SNP_anno.txt
```

<img src="C:\Users\jinwa42\AppData\Roaming\Typora\typora-user-images\image-20231107231139978.png" alt="image-20231107231139978" style="zoom:67%;" />





## full GWAS table



- After analyzing the GWAS results, we finally only select the coeff model

``` bash
#!/bin/bash -l

#SBATCH -A sens2022015
#SBATCH -p node
#SBATCH -n 2
#SBATCH -t 100:00:00
#SBATCH -J full_table
#SBATCH -e full_table.err
#SBATCH -o full_table.out

#cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/coeff/1_pdf/validated_mQTL

#for i in `cat validate_metabolite_list`;do cp ../../1_all_results/coeff.${i}.glm.linear .;done

#for i in `cat validate_metabolite_list`;do tar -czvf coeff.${i}.glm.linear.tar.gz coeff.${i}.glm.linear;done


cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/coeff/1_pdf/novel_mQTL

#for i in `cat  novel_metabolite_list` ;do cp ../../1_all_results/coeff.${i}.glm.linear .;done

for i in `cat  novel_metabolite_list` ; do tar -czvf coeff.${i}.glm.linear.tar.gz coeff.${i}.glm.linear;done
```





## conditional analysis



pwd: `/proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/conditional_analysis`

.
├── coeff_all_SNP.txt
├── coeff_leader_SNP.txt
├── coeff_name_conditional.txt

├── results/
└── sbatch/

``` bash 
cd sbatch

for i in `cat ../coeff_name_conditional.txt`; do cp coeff_zzzz_condi.sh coeff_${i}_condi.sh;sed -i "s/zzzz/$i/g" coeff_${i}_condi.sh; sbatch coeff_${i}_condi.sh; done

for i in `ls *.txt`;do sed -i '1d' ${i};done

cat *.txt > coeff_codition_SNP.txt

cp coeff_codition_SNP.txt  /proj/nobackup/sens2022015/wharf/jingw/jingw-sens2022015

to: C:\Users\jinwa42\OneDrive - Linköpings universitet\02_project\2022_09_03_mQTL\mQTL_wfr\output\03.GWAS_v7_conditional_analysis

next step: 03.GWAS_v7_conditional_analysis.Rmd
```

### coeff_zzzz_condi.sh

``` bash
#!/bin/bash -l

#SBATCH -A sens2022015
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J coeff_zzzz_condi
#SBATCH -e coeff_zzzz_condi.err
#SBATCH -o coeff_zzzz_condi.out


module load bioinfo-tools
module load plink
module load plink2

cd /proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/conditional_analysis

i="zzzz"


echo '################### coeff gwas for protein `$i`  ######################'


echo Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)


awk '{if($1 == "zzzz"){print $2}}' coeff_leader_SNP.txt > results/${i}_leader_SNP

awk '{if($1 == "zzzz"){print $2}}' coeff_all_SNP.txt > results/${i}_all_SNP


## re-calculated the P value, using leader SNP as convariate.


plink2 --bfile ../../HPA_101_qc_maf0.05 --extract results/${i}_all_SNP --condition-list results/${i}_leader_SNP --ci 0.95 --linear hide-covar --pheno ../mQTL.coeff.txt --pheno-name ${i} --covar ../mQTL.cov_log2.txt  --covar-name Sex, Age_b, BMI_b --covar-variance-standardize  --out results/condition

cd results

##### only the association with conditional P < 0.01 were considered to be independent mQTLs. 

awk 'BEGIN {FS="\t"}{OFS="\t"} {if($17 != "NA" && strtonum($17) < 0.01) {print "zzzz",$3}}' condition.${i}.glm.linear > coeff_${i}_codition_SNP.txt

chmod +x coeff_${i}_codition_SNP.txt


work_end_time=$(date +%s)
echo "#######-- cv gwas End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))

```



## extract geno



pwd: `/proj/sens2022015/jingw/wellness_20221101/GWAS/mQTL_v7/extract_geno`



### coeff_extract_geno.sh

```bash
#!/bin/bash -l

#SBATCH -A sens2022015
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J coeff_extract_geno
#SBATCH -e coeff_extract_geno.err
#SBATCH -o coeff_extract_geno.out

module load bioinfo-tools plink bcftools

#### extract all the snps with P <= 5e-08

plink --bfile ../../HPA_101_qc_maf0.05 --extract coeff_snp_5e-08 --recode vcf-iid  --out coeff_snp_5e-08

bcftools query -f '[%ID\t%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' coeff_snp_5e-08.vcf > coeff_snp_5e-08_geno

awk '{if ($7=="1/1"){print $1,$6,$5$5} else if ($7=="0/1") {print $1,$6,$4$5} else {print $1,$6,$4$4}}' coeff_snp_5e-08_geno > coeff_snp_5e-08_geno.txt

```



