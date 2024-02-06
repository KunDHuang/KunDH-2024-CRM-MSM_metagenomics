## First time setup

While no pre-installation is required for running workflows in this repository, you first need to download the repository to your local computer:
~~~bash
git clone https://github.com/KunDHuang/KunDH-2024-CRM-MSM_metagenomics.git
~~~ 

NOTE: As scripts provided in this repository depend on many third-party tools, we highly recommend [conda](https://conda.io/projects/conda/en/latest/index.html) environment for managing required dependencies which are listed in each of sections in the tutorials.

## General usages

This repository provides two types of utilities: 1) Python scripts and 2) Importable R functions

1) Python scripts are codes encapulated for executing specific analysis.
Example: [visualizing cumulative distribution function.](../docs/cumulative_distribution_function.md)
~~~bash
$ cumulative_distribution_function.py --input_table <example_reads_stats.tsv> --output_figure <nr_QC_reads_pairs.svg> --value_header <nr_QC_reads_pairs> --palette_map <reads_stats_color_map.tsv>
~~~  

2) Importable R functions are R codes wrapped up for solving specific problems and should be reusable by just importing their scripts. 

Example: [estimate PERMANOVA significance.](../docs/beta_diversity_analysis.md)

~~~R
>source(file = "path_to_the_package/KunDH-2023-CRM-MSM_metagenomics/scripts/functions/beta_diversity_funcs.R")
>est_permanova(mat = matrix, 
               md = metadata, 
               variable = "condom_use", 
               covariables = c("Antibiotics_6mo", "HIV_status", "inflammatory_bowel_disease", "BMI_kg_m2_WHO", "diet"),
               nper = 999, 
               to_rm = c("no_receptive_anal_intercourse"),
               by_method = "margin")

                           Df SumOfSqs      R2      F Pr(>F)   
condom_use                  4   1.2161 0.08194 1.5789  0.008 **
Antibiotics_6mo             2   0.4869 0.03281 1.2643  0.160   
HIV_status                  1   0.3686 0.02484 1.9146  0.030 * 
inflammatory_bowel_disease  1   0.2990 0.02015 1.5529  0.066 . 
BMI_kg_m2_WHO               5   1.8376 0.12382 1.9087  0.002 **
diet                        3   0.8579 0.05781 1.4853  0.036 * 
Residual                   49   9.4347 0.63571                 
Total                      65  14.8412 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
~~~
