# Logistic Regression Analyses: NABCR & SHHS Projects

## Overview
This repository contains logistic regression analyses on biomedical datasets, focusing on cancer survival and heart disease risk prediction. These projects demonstrate how to model and interpret key predictors related to health outcomes using rigorous statistical methods in R.

## Projects

### NABCR Logistic Regression Analysis

Examines the association between receptor level and survival status (death), adjusting for cancer stage, with stage-specific effect estimation and interaction testing.

### SHHS Logistic Regression Analysis

Performs a comprehensive analysis on the Scottish Heart Health Survey dataset to identify predictors of coronary heart disease (CHD) through univariable and multivariable logistic regression, confounding and interaction assessments, and model validation.

---

## Prerequisites

### Software

- R (version 4.0 or later recommended)

### Required R Packages

- For NABCR project:  
  `car`, `lmtest`, `broom`, `dplyr`, `ggplot2`

- For SHHS project:  
  `dplyr`, `car`, `ResourceSelection`, `lmtest`, `broom`, `sandwich`

Install missing packages in R with:

```r
install.packages(c("car", "lmtest", "broom", "dplyr", "ggplot2", "ResourceSelection", "sandwich"))
```

## Project Details

### 1. NABCR Logistic Regression Analysis

**Data Requirements**

CSV file named `NABCR.csv` with variables:

- `level`: receptor level (coded 1 or 2)  
- `stage`: cancer stage (coded 1 to 3)  
- `dead`: survival status (0 = alive, 1 = dead)  

The script converts numeric codes to factor variables with descriptive labels.

**Analysis Workflow**

- Load required packages and dataset  
- Prepare data (convert variables to factors)  
- Descriptive statistics and frequency tables  
- Logistic regression models:  
  - Model 1: Unadjusted effect of receptor level on death  
  - Model 2: Adjusted for cancer stage (confounding assessment)  
  - Model 3: Interaction model between receptor level and stage  
- Stage-specific logistic regression analyses  
- Visualization of odds ratios with 95% confidence intervals by cancer stage  
- Save model summaries and plots  

**Output Files**

- `output/tables/stage_specific_ORs.csv` — Odds ratios by cancer stage  
- `output/tables/model1_summary.csv` — Model 1 summary  
- `output/tables/model2_summary.csv` — Model 2 summary  
- `output/tables/model3_summary.csv` — Model 3 summary  
- `output/plots/stage_specific_ORs_plot.png` — Odds ratio plot by stage  

...


## 2. SHHS Logistic Regression Analysis

### Data Requirements

CSV file with CHD outcome and predictors including:

- Coronary heart disease (CHD) status  
- Total cholesterol  
- BMI  
- Smoking status  
- Age  
- Systolic blood pressure  
- Physical activity  

### Analysis Workflow

- Load dataset and prepare variables (categorical conversion, missing data check)  
- Descriptive statistics for predictors  
- Univariable logistic regression models  
- Multivariable logistic regression with confounder assessment  
- Confounding and interaction analyses  
- Functional form assessment of continuous predictors  
- Final model selection and validation (Hosmer-Lemeshow test)  
- Compare logistic, log-binomial, and Poisson regression models for CHD risk  
- Assess smoking effect using propensity score weighting with IPTW  
- Save tables and diagnostic plots  

### Output Files

- Summary tables and model outputs saved under `output/tables/`  
- Diagnostic and exploratory plots saved under `output/plots/`  
- Model comparison and IPTW analysis outputs included  

### How to Use

1. Place your input CSV files (`NABCR.csv` for NABCR project, and 'SHHS.csv' for SHHS project) in accessible directories.  
2. Update file paths in the corresponding R scripts to correctly point to your input data and desired output folders.  
3. Install all required R packages if missing.  
4. Run the R scripts for each project.  
5. Review generated output files (tables and plots) under the `output/` directory.  

### Directory Structure

```
output/
├── tables/ # CSV summary tables and logistic regression model results
└── plots/ # Diagnostic and exploratory plots (e.g., odds ratio plots, IPTW weight histograms)
```


If you need further assistance or code examples for either project, feel free to ask!
