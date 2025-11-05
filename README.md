# 1. Overview
This repository contains the following files with their serve as below.

- Main functions.R: It includes all functions used in our paper, along with the functions for other competitors such as the Trans-Lasso. The specific role of each function can be found in the corresponding inline anotations.
- Adaptive rule.R : The RIW-TL using adaptive selection rule for the tuning parameters.
- Example.R: It provides an example of running the above functions and specifically how to use our method.
- gene_id.csv: This file contains the covariate information that we focus on in real analysis.
- target-list.RData: Since the collected GTEx data (raw data) is too large to upload, we upload the processed target data and source data (below) in our real analysis.
- source-list.RData: The source data in our real analysis.
- data processing.R: It includes the preprocessing steps for raw data.
- data analysis.R: This file records the codes used in analyzing the processed GTEx data.
- Real data analysis.md: This document provides the code, explanations, and relevant output results for the real data analysis.

# 2. Setup Instructions
- R version: It is sufficient to ensure that the R version is R 4.3.3 or updated ones.
- Packages: Our numerical study relied on the following packages: glmnet, ncvreg, ks, mvtnorm, bayesdistreg, Matrix and stringi. All these packages can be installed using conventional methods, such as input "iinstall.packages("package name") in R.

# 3. Data Acquisition
The analyzed real data has been provided in the file "data_all.RData". However, the data size is too large to upload. To solve the problem, we have uploaded the processed data (i.e., target-list.RData and source-list.RData). If one is desirable to view the raw data, please contact us with email: a1261110123@163.com

# 4. Workflow Guide
## 4.1 Guide for simulation study
Step 1: Execute the script Functions.R for comprehensiveness;
Step 2: Execute the script Example.R, so that it can output the parameter estimates and sample utilization rates for each method. Futhermore, for various experients such as the varying $n_0$ or informative level $d$, one can adjust the corrrepsonding setting. Details can be refer to the inline annotation.

## 4.2 Guide for real data analysis
Step 1: Open the scipt "Real data processing.R" and load the data in "data_all.RData" (If the data is available; if not, one can jump to the Step 3, wherein the processed data have been uploaded);
Step 2: Execute the script "Real data processing.R" and then output the processed data, which we call "target.list" for the target data, and "source.list" for the source data.
Step 3: Execute the script "Real data analysis.R" and one can get the results for real data analysis. Using the inner plot function, one can get the corresponding results on prediction errors.
