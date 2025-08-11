# Logistic Regression Analyses: NABCR & SHHS Projects

## Overview  
This repository contains logistic regression analyses on biomedical datasets, focusing on cancer survival and heart disease risk prediction. These projects demonstrate how to model and interpret key predictors related to health outcomes using rigorous statistical methods in R.

---

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

