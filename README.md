# Burden analysis of rare predicted loss of function variants in COMMD related genes for Parkinson's disease.
Bernabe Bustos, Ph.D. (Krainc Lab, Department of Neurology, Northwestern University Feinberg School of Medicine)

Genetic analysis for the paper "Commander complex regulates lysosomal function and is implicated in Parkinson’s disease risk" published in the journal Science https://www.science.org/doi/10.1126/science.adq6650.

### This pipeline starts with already QCed and processed AMP-PD and UK Biobank datasets as described in the methods section of our publication. This includes the creation of the phenotype file that includes the covariates: sex, age at onset for cases and recruitment for controls, and principal components 1 to 10.

### Exected QCs from AMP-PD and UK Biobank

#### For AMP-PD genomes:
##### 1. All variants must have genotype depth (DP) above 10, genotype quality (GQ) above 20.
##### 2. Left aligned indels, and split multi allelic variants.
##### 3. Removed samples with IBD pi_hat scores > 0.187.
##### 4. Removed samples with low coverage (< 90%).
##### 5. Removed samples with more than 6 standard deviations from the 1000 genomes European population in PCA analysis.

#### For UK Biobank exomes:
##### 1. Excluded variants that falied the 90% readings with DP≥10 filter provided by the UK Biobank team. (file OQFE.90pct10dp_qc_variants.txt).
##### 2. Removed samples with IBD pi_hat scores > 0.187.
##### 3. Removed samples with low coverage (< 90%).
##### 4. Kept samples with reported genetic European ancestry.
##### 5. Removed control samples with neurological/neurodegenerative diseases (see paper methods).


### **Variant annotation and selection for AMP-PD genomes and UK Biobank exomes data**

```

# Annotation with ANNOVAR (version from Jun 8 2020)

# Code used for AMP-PD run in a Linux bash console. It applies to UK Biobank as well

./table_annovar.pl chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.vcf.gz \
/projects/b1049/genetics_programs/annovar_2022/annovar/humandb/ -buildver hg38 \
--thread 20 \
-out chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.annotated \
-remove -protocol refGene,genomicSuperDups,gnomad211_exome,gnomad30_genome \
-operation g,r,f,f \
-nastring . \
-vcfinput \
-otherinfo

# Selection of rare predicted loss of function damaging variants

cut -f1-9,11,12,29,47 chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.annotated.hg38_multianno.txt | awk '$6 == "exonic" || $6 == "exonic;splicing" || $6 == "splicing" {print}' | awk '$9 == "frameshift" || $9 == "stopgain" || $9 == "stoploss" {print $(NF-7)}' > AMPPD.MAF01.LOF.txt

# Extraction of selected variants from VCF file

bcftools view -i 'I=@AMPPD.MAF01.LOF.txt' chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.vcf.gz -Oz -o chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.LOF.vcf.gz 

# Index VCF file

tabix -p vcf chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.LOF.vcf.gz 

```
## **Burden analysis in RVTest version 20171009**
### **Gene-wise burden**
```
# AMP-PD

./rvtest --inVcf chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.LOF.vcf.gz  \
--pheno AMPPD_2024_pheno_file.txt \
--pheno-name PHENO \
--covar AMPPD_2024_pheno_file.txt \
--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,SITE \
--hide-covar \
--gene VPS16,VPS11,VPS18,VPS41,VPS33A,VPS39,COMMD1,COMMD2,COMMD3,COMMD3-BMI1,COMMD4,COMMD5,COMMD6,COMMD7,COMMD8,COMMD9,COMMD10,CCDC22,CCDC93,DENND10,VPS35L,VPS26C,VPS29 \
--geneFile refFlat_hg38.txt \
--out AMPPD_2024.LOF.COMMD_genes \
--burden cmcWald \
--numThread 20

# UK Biobank

./rvtest --inVcf ukb23158_c1-X_b0_v1.QC_Samples_Vars_MAF01.LOF.vcf.gz \
--peopleIncludeFile clinically_accurate.PD_CC.IDs.txt \
--pheno UKBB_PD_WES_470k.phe \
--pheno-name PHENO \
--covar UKBB_PD_WES_470k.phe \
--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--hide-covar \
--geneFile refFlat_hg38.txt \
--gene VPS16,VPS11,VPS18,VPS41,VPS33A,VPS39,COMMD1,COMMD2,COMMD3,COMMD3-BMI1,COMMD4,COMMD5,COMMD6,COMMD7,COMMD8,COMMD9,COMMD10,CCDC22,CCDC93,DENND10,VPS35L,VPS26C,VPS29 \
--out UKBB.LOF.COMMD_genes \
--burden cmcWald
```
### **Gene-set burden**
```
# Prepare set file for 5 gene sets:

1. COMMD_ring = COMMD1, COMMD2, COMMD3, COMMD3-BMI1, COMMD4, COMMD5, COMMD6, COMMD7, COMMD8, COMMD9, COMMD10.

2. HOPS_complex = VPS16, VPS11, VPS18, VPS41, VPS33A, VPS39.

3. Commander_complex = COMMD1, COMMD2, COMMD3, COMMD3-BMI1, COMMD4, COMMD5, COMMD6, COMMD7, COMMD8, COMMD9, COMMD10, CCDC22, CCDC93, DENND10, VPS35L, VPS26C, VPS29.

4. CCC_complex = COMMD1, COMMD2, COMMD3, COMMD3-BMI1, COMMD4, COMMD5, COMMD6, COMMD7, COMMD8, COMMD9, COMMD10, CCDC22, CCDC93.

5. Retriever_complex = VPS35L, VPS26C, VPS29.


# Extraction of gene coordinates (chr:start-end) from gene-wise results: UKBB.LOF.COMMD_genes

head COMMD_ring.txt
COMMD1
COMMD2
COMMD3
...

# Generate the file 'All_sets.txt' that contains the gene coordinates of all genes for each gene set.

for i in COMMD_ring HOPS_complex Commander_complex CCC_complex Retriever_complex
do

grep -F -wf $i UKBB.LOF.COMMD_genes.CMCWald.assoc | cut -f2 | sort | uniq | sed -z 's/\n/,/g;s/,$/\n/' >> All_sets.set
 
done

# Create the file 'set_names.txt' with the gene sets names.

head set_names.txt
COMMD_ring
HOPS_complex
Commander_complex
CCC_complex
Retriever_complex

# Append the gene set names to the coordinates.

paste set_names.txt All_sets.set > All_sets_wNames.set

# Example of a resulting 'All_sets_wNames.set' file.

COMMD_ring      10:22316382-22331485,10:22316387-22320305,11:36272291-36289424,11:36272291-36289424,11:36272291-36289424,11:36272291-36289424,13:75525213-75537872,13:75525218-75537872,13:75525218-75549439,13:75525218-75537872,13:75525218-75537872,15:75336062-75340273,15:75336062-75343226,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,15:75336062-75340273,20:32702700-32743467,20:32702700-32743467,2:61888699-62141185,2:61888699-62141185,2:61905641-62136070,2:61888390-62141185,3:149738471-149752489,3:149738469-149752499,4:47450788-47463702,4:47450788-47463730,5:116085277-116190938,5:116085461-116293290,5:116085024-116190940,5:116085024-116293287,5:116085024-116190940,8:144850165-144853067,8:144850175-144853048,8:144850175-144853048,8:144850175-144853556

HOPS_complex    11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89416-103571,11_KN196481v1_fix:89611-103571,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119067817-119081972,11:119068012-119081972,11:119067817-119081972,11:119067817-119081972,12:122229563-122266494,12:122229563-122266494,12:122229563-122266494,12:122229563-122266494,12:122252861-122266494,15:40894449-40903975,15:42158710-42208304,15:42158700-42208304,20:2840696-2866732,20:2840744-2866732,7:38722973-38909191,7:38723942-38909200

# Gene set burden analysis

# AMP-PD

./rvtest --inVcf chr1-X.OnlyPASS_DP10_GQ20_CR90.MultiSplit.SNPsInDels.RawID.bcftools.PD_CC_Nov2022.Plink_rareMAF01.LOF.vcf.gz  \
--pheno AMPPD_2024_pheno_file.txt \
--pheno-name PHENO \
--covar AMPPD_2024_pheno_file.txt \
--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,SITE \
--hide-covar \
--setFile All_sets_wNames.set \
--out AMPPD_2024.LOF.COMMD_gene_sets \
--burden cmcWald \
--numThread 20

# UK Biobank

./rvtest --inVcf ukb23158_c1-X_b0_v1.QC_Samples_Vars_MAF01.LOF.vcf.gz \
--peopleIncludeFile clinically_accurate.PD_CC.IDs.txt \
--pheno UKBB_PD_WES_470k.phe \
--pheno-name PHENO \
--covar UKBB_PD_WES_470k.phe \
--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--hide-covar \
--setFile All_sets_wNames.set \
--out UKBB.LOF.COMMD_gene_sets \
--burden cmcWald
```

### Result

```

# Extract relevant columns from output file, select rows with the gene-set association (first for each gene set), and sort for significance.

cut -f1,3- AMPPD_2024.LOF.COMMD_gene_sets.CMCWald.assoc | awk '!seen[$1]++ {print}' | sort -k8g

Range   N_INFORMATIVE   NumVar  NumPolyVar      NonRefSite      Beta    SE      Pvalue
COMMD_ring      5659    27      18      23      1.93614 0.656567        0.00318928
CCC_complex     5659    30      21      26      1.63563 0.596798        0.0061314
Commander_complex       5659    42      27      35      0.914293        0.440866        0.0380925
HOPS_complex    5659    20      11      67      0.229454        0.290973        0.430361
Retriever_complex       5659    10      7       10      -0.194984       0.699447        0.780423


```
### **Meta-analysis of gene set analysis for AMP-PD and UK Biobank**
```
# Analysis performed in R version 4.2.3

# Install and load required packages

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)
}

if (!require("stringr", quietly = TRUE)) {
  install.packages("stringr")
  library(stringr)
} else {
  library(stringr)
}

if (!require("metafor", quietly = TRUE)) {
  install.packages("metafor")
  library(metafor)
} else {
  library(metafor)
}

# Transform effect sizes to OR and confidence intervals.

# Read the data file (adjust the file name and separator if needed)
data <- read.table("input_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Calculate the odds ratio and the lower and upper 95% confidence limits
data$OR      <- exp(data$Beta)
data$Lower95 <- exp(data$Beta - 1.96 * data$SE)
data$Upper95 <- exp(data$Beta + 1.96 * data$SE)

# Select the desired columns to output
result <- data[, c("GeneSet", "OR", "Lower95", "Upper95", "Pvalue")]


# Define the data.

#data <- results

# Or define the values manually

data <- data.frame(
  Geneset = c(
    "COMMD_ring_AMPPD", "COMMD_ring_UKBB",
    "CCC_complex_AMPPD", "CCC_complex_UKBB",
    "Commander_complex_AMPPD", "Commander_complex_UKBB",
    "HOPS_complex_AMPPD", "HOPS_complex_UKBB",
    "Retriever_complex_AMPPD", "Retriever_complex_UKBB"
  ),
  OR = c(3.54, 2.13, 3.51, 2.33, 2.25, 2.09, 1.00, 0.53, 0.98, 1.59),
  CI95_L = c(1.59, 1.33, 1.64, 1.53, 1.19, 1.42, 0.64, 0.22, 0.29, 0.58),
  CI95_U = c(7.88, 3.40, 7.52, 3.53, 4.24, 3.07, 1.57, 1.29, 3.35, 4.37),
  Pvalue = c(2.02e-03, 1.59e-03, 1.22e-03, 7.85e-05, 1.24e-02, 1.92e-04, 9.83e-01, 1.63e-01, 9.74e-01, 3.67e-01),
  stringsAsFactors = FALSE
)

# Parse data and calculate effect sizes

data <- data %>%
  mutate(
    gene_set = str_replace(Geneset, "_(AMPPD|UKBB)$", ""),
    study = ifelse(str_detect(Geneset, "AMPPD$"), "AMP-PD", "UK Biobank"),
    logOR = log(OR),
    SE = (log(CI95_U) - log(OR)) / 1.96  # SE derived from 95% CI
  )

# Perform meta-analysis for each gene set

results_list <- list()
gene_sets <- unique(data$gene_set)

for (gs in gene_sets) {
  subset_dat <- subset(data, gene_set == gs)
  
  if (nrow(subset_dat) < 2) {
    warning(paste("Not enough studies for gene set:", gs, "- skipping meta-analysis."))
    next
  }
  
  res <- tryCatch(
    rma(yi = subset_dat$logOR, sei = subset_dat$SE, method = "FE"),
    error = function(e) {
      warning(paste("Meta-analysis failed for gene set:", gs))
      return(NULL)
    }
  )
  
  if (!is.null(res)) {
    results_list[[gs]] <- list(
      OR = round(exp(res$b), 4),
      CI.Lower = round(exp(res$ci.lb), 4),
      CI.Upper = round(exp(res$ci.ub), 4),
      p.value = signif(res$pval, 4),
      I2 = round(res$I2, 2)
    )
  }
}


# Save meta-Analysis results to a table

if (length(results_list) > 0) {
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    as.data.frame(t(unlist(x)))
  }))
  rownames(results_df) <- names(results_list)
  
  # Reorder columns for clarity
  results_df <- results_df[, c("OR", "CI.Lower", "CI.Upper", "p.value", "I2")]
  
  # Save the results to a TSV file
  write.table(results_df, file = "meta_analysis_results.tsv", sep = "\t", 
              quote = FALSE, col.names = NA)
  
  # Print the results to the console
  cat("\nMeta-Analysis Results:\n")
  print(results_df)
} else {
  cat("No meta-analysis results to compile.\n")
}

```
### Gene-set meta-analysis results

```
                      OR CI.Lower CI.Upper   p.value    I2
COMMD_ring        2.4241   1.6188   3.6299 1.721e-05 13.35
CCC_complex       2.5594   1.7772   3.6858 4.418e-07  0.00
Commander_complex 2.1319   1.5346   2.9616 6.370e-06  0.00
HOPS_complex      0.8782   0.5873   1.3132 5.269e-01 35.76
Retriever_complex 1.3079   0.5991   2.8555 5.004e-01  0.00
```
