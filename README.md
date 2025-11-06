# 1. Overview
This repository contains the following files, each serving a specific purpose:

- Main functions.R: It includes main functions used in our paper, along with functions for other competitors such as Trans-Lasso. The specific role of each function can be found in the corresponding inline annotations.
- Adaptive rule.R : This file includes the RIW-TL function with an adaptive selection rule for tuning parameters.
- Example.R: It provides an example of running the above functions and demonstrates how to use our method.
- gene_id.csv: This file contains the covariate information used in the real data analysis.
- target-list.RData: The processed target datasets in our real analysis.
- source-list.RData: The processed source datasets in our real analysis.
- data processing.R: It includes the preprocessing steps for raw data.
- data analysis.R: This file records the codes used in analyzing the processed GTEx data (if the GTEx data is available).
- Real data analysis.md: This document provides the code, explanations, and relevant output results for the real data analysis. The output (in PNG format) is stored in the "resources" folder.

# 2. Setup Instructions
- R version: Please ensure that your R version is 4.3.3 or later.
- Packages: Our numerical study relies on the following packages: glmnet, ncvreg, ks, mvtnorm, bayesdistreg, Matrix, and stringi. These can be installed using standard methods, such as "install.packages('package name')" in R.

# 3. Data Acquisition
The raw data can be available at https://pan.baidu.com/s/1pO1rn9xL6p4sTeWr_bBcBw?pwd=ija2. However, this dataset may be relatively large, which could result in longer download time. A more convenient approach is to directly use the provided processed dataset, i.e., target-list.RData and source-list.RData.

# 4. Workflow Guide
## 4.1 Guide for simulation study
- Step 1: Execute the script Functions.R to load all functions.
- Step 2: Run Example.R to obtain parameter estimates and sample utilization rates for each method. Futhermore, For different experimental settings (e.g., varying $n_0$ or information level $d$), one can adjust the corresponding parameters as indicated in the inline annotations.

## 4.2 Guide for real data analysis
- Step 1: Open the script Real data processing.R and load the data from data_all.RData (if available). If not, proceed directly to Step 3 using the processed data provided.
- Step 2: Execute Real data processing.R to generate the processed data, referred to as target.list (target data) and source.list (source data).
- Step 3: Run Real data analysis.R to obtain the results of the real data analysis using 100 replicates.
