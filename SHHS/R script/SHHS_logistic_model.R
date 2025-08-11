#########################################
# SHHS Logistic Regression Analysis
########################################

# 1. Load Required Libraries -----------------------------------
library(dplyr)
library(car)                  # For Anova()
library(ResourceSelection)    # For hoslem.test
library(lmtest)               # For likelihood ratio test (LR.test)

# 2. Load Data -------------------------------------------------

# Read CSV data file into R
mydata1 = read.csv(file = "................................./data/SHHS.csv", header = TRUE)

# Check structure of the dataset to understand variable types
str(mydata1)


# 3. Data Preparation -------------------------------------------

# Convert categorical variables to factors for modeling
mydata1$smoking_c = as.factor(mydata1$smoking)
mydata1$activity_c = as.factor(mydata1$activity)

# Check for missing values across dataset
total_missing = sum(is.na(mydata1))
cat("Total missing values in dataset:", total_missing, "\n")


# 4. Descriptive Statistics -------------------------------------

# Frequency tables for categorical variables to check distribution
table(mydata1$smoking_c)
table(mydata1$activity_c)
table(mydata1$chd)

# Calculate descriptive statistics (mean, sd, missing) for continuous variables
# Useful for understanding data distribution and quality
descriptive_stats <- mydata1 %>%
  summarise(
    n_age = sum(!is.na(age)), missing_age = sum(is.na(age)),
    n_totchol = sum(!is.na(totchol)), missing_totchol = sum(is.na(totchol)),
    n_bmi = sum(!is.na(bmi)), missing_bmi = sum(is.na(bmi)),
    n_systol = sum(!is.na(systol)), missing_systol = sum(is.na(systol)),
    mean_age = mean(age, na.rm = TRUE),
    mean_totchol = mean(totchol, na.rm = TRUE),
    mean_bmi = mean(bmi, na.rm = TRUE),
    mean_systol = mean(systol, na.rm = TRUE),
    sd_age = sd(age, na.rm = TRUE),
    sd_totchol = sd(totchol, na.rm = TRUE),
    sd_bmi = sd(bmi, na.rm = TRUE),
    sd_systol = sd(systol, na.rm = TRUE)
  )

# Calculate 95% Confidence Intervals (CI) for means of continuous variables
ci_continuous <- mydata1 %>%
  summarise(
    ci_age = mean(age, na.rm = TRUE) + c(-1, 1) * qt(0.975, df = sum(!is.na(age)) - 1) * sd(age, na.rm = TRUE) / sqrt(sum(!is.na(age))),
    ci_totchol = mean(totchol, na.rm = TRUE) + c(-1, 1) * qt(0.975, df = sum(!is.na(totchol)) - 1) * sd(totchol, na.rm = TRUE) / sqrt(sum(!is.na(totchol))),
    ci_bmi = mean(bmi, na.rm = TRUE) + c(-1, 1) * qt(0.975, df = sum(!is.na(bmi)) - 1) * sd(bmi, na.rm = TRUE) / sqrt(sum(!is.na(bmi))),
    ci_systol = mean(systol, na.rm = TRUE) + c(-1, 1) * qt(0.975, df = sum(!is.na(systol)) - 1) * sd(systol, na.rm = TRUE) / sqrt(sum(!is.na(systol)))
  )

# Output summary statistics and confidence intervals
list(
  descriptive_stats = descriptive_stats,
  confidence_intervals = ci_continuous
)

# Save the results ---
output_dir <- "................................./SHHS/output"
tables_dir <- file.path(output_dir, "tables")

if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(tables_dir)) dir.create(tables_dir)

desc_stats_df <- as.data.frame(descriptive_stats)
ci_df <- data.frame(
  Variable = c("age", "totchol", "bmi", "systol"),
  Lower_CI = c(ci_continuous$ci_age[1], ci_continuous$ci_totchol[1], ci_continuous$ci_bmi[1], ci_continuous$ci_systol[1]),
  Upper_CI = c(ci_continuous$ci_age[2], ci_continuous$ci_totchol[2], ci_continuous$ci_bmi[2], ci_continuous$ci_systol[2])
)

write.csv(desc_stats_df, file = file.path(tables_dir, "descriptive_stats_continuous.csv"), row.names = FALSE)
write.csv(ci_df, file = file.path(tables_dir, "confidence_intervals_continuous.csv"), row.names = FALSE)

cat("Descriptive statistics and confidence intervals saved to:", tables_dir, "\n")

#------------------------------------------------------------------------
# BUILDING A PARSIMONIOUS LOGISTIC REGRESSION MODEL FOR CHD PREDICTION
#------------------------------------------------------------------------

# Objective:
# To develop a parsimonious and interpretable logistic regression model predicting
# coronary heart disease (CHD) by selecting key predictors and assessing confounding,
# functional forms, and interactions.

# Modeling steps:
# 1) Fit univariable logistic regression models to evaluate each predictor separately.
# 2) Build a multivariable logistic regression model using significant predictors.
# 3) Assess confounding by examining changes in coefficients when adding/removing variables.
# 4) Reconsider additional variables excluded initially to check for confounding effects.
# 5) Explore functional forms of continuous variables (e.g., categorize or use splines).
# 6) Test for potential interaction effects between predictors.
# 7) Validate the final model using goodness-of-fit tests (e.g., Hosmer-Lemeshow test).


## STEP 1: Univariable Logistic Regression Models
# Fit separate logistic regression models for each predictor individually.
# This step screens variables for possible inclusion based on their unadjusted associations with CHD.

model1 <- glm(chd ~ age, family = binomial(logit), data = mydata1)
model2 <- glm(chd ~ totchol, family = binomial(logit), data = mydata1)
model3 <- glm(chd ~ bmi, family = binomial(logit), data = mydata1)
model4 <- glm(chd ~ systol, family = binomial(logit), data = mydata1)
model5 <- glm(chd ~ smoking_c, family = binomial(logit), data = mydata1)
model6 <- glm(chd ~ activity_c, family = binomial(logit), data = mydata1)

# Calculate odds ratios (OR) and 95% confidence intervals (CI) for each model's coefficients.
model1OR <- exp(cbind(OR = coef(model1), confint.default(model1)))
model2OR <- exp(cbind(OR = coef(model2), confint.default(model2)))
model3OR <- exp(cbind(OR = coef(model3), confint.default(model3)))
model4OR <- exp(cbind(OR = coef(model4), confint.default(model4)))
model5OR <- exp(cbind(OR = coef(model5), confint.default(model5)))
model6OR <- exp(cbind(OR = coef(model6), confint.default(model6)))

# Perform likelihood ratio tests (LRT) to assess the significance of each predictor.
model1g <- Anova(model1, type = 3)
model2g <- Anova(model2, type = 3)
model3g <- Anova(model3, type = 3)
model4g <- Anova(model4, type = 3)
model5g <- Anova(model5, type = 3)
model6g <- Anova(model6, type = 3)

# Extract coefficients, standard errors, ORs, confidence limits, and p-values for a summary table.
coef_vals <- sapply(list(model1, model2, model3, model4, model5, model6), 
                    function(m) if(length(coef(m)) >= 2) coef(m)[2] else NA)

se_vals <- sapply(list(model1, model2, model3, model4, model5, model6), 
                  function(m) if(length(coef(m)) >= 2) sqrt(vcov(m)[2,2]) else NA)

or_vals <- c(
  if(length(model1OR) >= 2) model1OR[2,1] else NA,
  if(length(model2OR) >= 2) model2OR[2,1] else NA,
  if(length(model3OR) >= 2) model3OR[2,1] else NA,
  if(length(model4OR) >= 2) model4OR[2,1] else NA,
  if(length(model5OR) >= 2) model5OR[2,1] else NA,
  if(length(model6OR) >= 2) model6OR[2,1] else NA
)

lcl_vals <- c(
  if(length(model1OR) >= 2) model1OR[2,2] else NA,
  if(length(model2OR) >= 2) model2OR[2,2] else NA,
  if(length(model3OR) >= 2) model3OR[2,2] else NA,
  if(length(model4OR) >= 2) model4OR[2,2] else NA,
  if(length(model5OR) >= 2) model5OR[2,2] else NA,
  if(length(model6OR) >= 2) model6OR[2,2] else NA
)

ucl_vals <- c(
  if(length(model1OR) >= 2) model1OR[2,3] else NA,
  if(length(model2OR) >= 2) model2OR[2,3] else NA,
  if(length(model3OR) >= 2) model3OR[2,3] else NA,
  if(length(model4OR) >= 2) model4OR[2,3] else NA,
  if(length(model5OR) >= 2) model5OR[2,3] else NA,
  if(length(model6OR) >= 2) model6OR[2,3] else NA
)

lrt_vals <- c(
  if(length(model1g$`LR Chisq`) > 0) model1g$`LR Chisq` else NA,
  if(length(model2g$`LR Chisq`) > 0) model2g$`LR Chisq` else NA,
  if(length(model3g$`LR Chisq`) > 0) model3g$`LR Chisq` else NA,
  if(length(model4g$`LR Chisq`) > 0) model4g$`LR Chisq` else NA,
  if(length(model5g$`LR Chisq`) > 0) model5g$`LR Chisq` else NA,
  if(length(model6g$`LR Chisq`) > 0) model6g$`LR Chisq` else NA
)

pval_vals <- c(
  if(length(model1g$`Pr(>Chisq)`) > 0) model1g$`Pr(>Chisq)` else NA,
  if(length(model2g$`Pr(>Chisq)`) > 0) model2g$`Pr(>Chisq)` else NA,
  if(length(model3g$`Pr(>Chisq)`) > 0) model3g$`Pr(>Chisq)` else NA,
  if(length(model4g$`Pr(>Chisq)`) > 0) model4g$`Pr(>Chisq)` else NA,
  if(length(model5g$`Pr(>Chisq)`) > 0) model5g$`Pr(>Chisq)` else NA,
  if(length(model6g$`Pr(>Chisq)`) > 0) model6g$`Pr(>Chisq)` else NA
)

# Combine all results into a single data frame for easy comparison.
STEP1 <- data.frame(
  Coefficient = coef_vals,
  Std_Error = se_vals,
  OR = or_vals,
  Lower_CI = lcl_vals,
  Upper_CI = ucl_vals,
  LRT_ChiSq = lrt_vals,
  p_value = pval_vals
)

# Display summary of univariate models.
print(STEP1)


## STEP 2: Multivariable Logistic Regression Model with Selected Variables
# Fit a multivariable model including predictors significant in Step 1.
# This estimates adjusted effects controlling for other variables.

step2model <- glm(chd ~ age + totchol + bmi + systol + smoking_c, 
                  family = binomial(logit), data = mydata1)

# Show model summary with coefficients and statistical significance.
summary(step2model)

# Perform Type 3 ANOVA to test the overall significance of each predictor.
Anova(step2model, type = 3)


## STEP 3: Assess Confounding by Removing Variables One at a Time
# Fit reduced models excluding one variable at a time (e.g., age or bmi)
# and examine how coefficients of remaining variables change to assess confounding.

# Step 3: Assess Confounding by Removing Variables One at a Time

# Get all smoking_c coefficient names (handles factor with multiple levels)
smoking_coef_names <- grep("^smoking_c", names(coef(step2model)), value = TRUE)

# Model excluding age
step3model1 <- glm(chd ~ totchol + bmi + systol + smoking_c, family = binomial(logit), data = mydata1)

# Calculate percent change safely for numeric variables
Delta_totchol <- ((coef(step2model)["totchol"] - coef(step3model1)["totchol"]) / abs(coef(step3model1)["totchol"])) * 100
Delta_bmi     <- ((coef(step2model)["bmi"] - coef(step3model1)["bmi"]) / abs(coef(step3model1)["bmi"])) * 100
Delta_systol  <- ((coef(step2model)["systol"] - coef(step3model1)["systol"]) / abs(coef(step3model1)["systol"])) * 100

# Calculate percent changes for all smoking_c coefficients
Delta_smoking_c <- sapply(smoking_coef_names, function(coef_name) {
  ((coef(step2model)[coef_name] - coef(step3model1)[coef_name]) / abs(coef(step3model1)[coef_name])) * 100
})

confounding1 <- data.frame(
  Variable = c("totchol", "bmi", "systol", smoking_coef_names),
  Coef_with_Age = c(coef(step2model)["totchol"], coef(step2model)["bmi"], coef(step2model)["systol"], coef(step2model)[smoking_coef_names]),
  Coef_without_Age = c(coef(step3model1)["totchol"], coef(step3model1)["bmi"], coef(step3model1)["systol"], coef(step3model1)[smoking_coef_names]),
  Percent_Change = c(Delta_totchol, Delta_bmi, Delta_systol, Delta_smoking_c)
)

print(confounding1)


# Model excluding bmi
step3model2 <- glm(chd ~ age + totchol + systol + smoking_c, family = binomial(logit), data = mydata1)

Delta_age <- ((coef(step2model)["age"] - coef(step3model2)["age"]) / abs(coef(step3model2)["age"])) * 100
Delta_totchol <- ((coef(step2model)["totchol"] - coef(step3model2)["totchol"]) / abs(coef(step3model2)["totchol"])) * 100
Delta_systol <- ((coef(step2model)["systol"] - coef(step3model2)["systol"]) / abs(coef(step3model2)["systol"])) * 100

Delta_smoking_c <- sapply(smoking_coef_names, function(coef_name) {
  ((coef(step2model)[coef_name] - coef(step3model2)[coef_name]) / abs(coef(step3model2)[coef_name])) * 100
})

confounding2 <- data.frame(
  Variable = c("age", "totchol", "systol", smoking_coef_names),
  Coef_with_BMI = c(coef(step2model)["age"], coef(step2model)["totchol"], coef(step2model)["systol"], coef(step2model)[smoking_coef_names]),
  Coef_without_BMI = c(coef(step3model2)["age"], coef(step3model2)["totchol"], coef(step3model2)["systol"], coef(step3model2)[smoking_coef_names]),
  Percent_Change = c(Delta_age, Delta_totchol, Delta_systol, Delta_smoking_c)
)

print(confounding2)

# Combine confounding results
colnames(confounding2) <- colnames(confounding1)
STEP3 <- rbind(
  cbind(Model = "Excluding Age", confounding1),
  cbind(Model = "Excluding BMI", confounding2)
)

## STEP 4: Reconsider Variables Excluded Earlier
# Add previously excluded variables one at a time to the multivariable model
# to assess their contribution.

test1 <- glm(chd ~ totchol + bmi + systol + smoking_c + activity_c, 
             family = binomial(logit), data = mydata1)

# Calculate Wald Chi-square statistics and p-values for each variable.
wald_results <- data.frame(Variable = character(), Wald_ChiSq = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (i in 2:length(coef(test1))) {
  wald_stat <- (coef(test1)[i] / sqrt(vcov(test1)[i, i]))^2
  p_val <- pchisq(wald_stat, df = 1, lower.tail = FALSE)
  wald_results <- rbind(wald_results, data.frame(Variable = names(coef(test1))[i], Wald_ChiSq = wald_stat, p_value = p_val))
}

print(wald_results)


## STEP 5: Assess Functional Form of Continuous Variables
# Categorize continuous variables (Total Cholesterol, BMI, Systolic BP) into quartiles
# to check for possible non-linear relationships with CHD.

# --- Total Cholesterol ---
mydata1$TOTCHOL4 <- cut(mydata1$totchol,
                        breaks = quantile(mydata1$totchol, probs = seq(0, 1, 0.25), na.rm = TRUE),
                        include.lowest = TRUE, labels = FALSE)
mydata1$TOTCHOL4_f <- as.factor(mydata1$TOTCHOL4)

# --- BMI ---
mydata1$BMI4 <- cut(mydata1$bmi,
                    breaks = quantile(mydata1$bmi, probs = seq(0, 1, 0.25), na.rm = TRUE),
                    include.lowest = TRUE, labels = FALSE)
mydata1$BMI4_f <- as.factor(mydata1$BMI4)

# --- Systolic Blood Pressure ---
mydata1$SYSTOL4 <- cut(mydata1$systol,
                       breaks = quantile(mydata1$systol, probs = seq(0, 1, 0.25), na.rm = TRUE),
                       include.lowest = TRUE, labels = FALSE)
mydata1$SYSTOL4_f <- as.factor(mydata1$SYSTOL4)

# Fit logistic model using quartile categories for each variable separately 
# while adjusting for smoking_c (example)

# Model 1: Total Cholesterol quartiles + smoking
model_totchol <- glm(chd ~ TOTCHOL4_f + smoking_c, family = binomial(logit), data = mydata1)
summary(model_totchol)

# Model 2: BMI quartiles + smoking
model_bmi <- glm(chd ~ BMI4_f + smoking_c, family = binomial(logit), data = mydata1)
summary(model_bmi)

# Model 3: Systolic BP quartiles + smoking
model_systol <- glm(chd ~ SYSTOL4_f + smoking_c, family = binomial(logit), data = mydata1)
summary(model_systol)


# Prepare data frames for plotting log-odds by quartile for each variable

plotdata_totchol <- data.frame(
  Quartile = levels(mydata1$TOTCHOL4_f),
  LogOdds = c(0, coef(model_totchol)[2:4])  # Reference quartile log-odds = 0
)

plotdata_bmi <- data.frame(
  Quartile = levels(mydata1$BMI4_f),
  LogOdds = c(0, coef(model_bmi)[2:4])
)

plotdata_systol <- data.frame(
  Quartile = levels(mydata1$SYSTOL4_f),
  LogOdds = c(0, coef(model_systol)[2:4])
)


# Plot log-odds vs quartiles for each variable

par(mfrow = c(1, 3))  # 3 plots in 1 row

plot(as.numeric(plotdata_totchol$Quartile), plotdata_totchol$LogOdds,
     xlab = 'Total Cholesterol Quartile',
     ylab = 'Log Odds',
     main = 'Log Odds by Total Cholesterol Quartile',
     type = "b", pch = 19, lty = 1)

plot(as.numeric(plotdata_bmi$Quartile), plotdata_bmi$LogOdds,
     xlab = 'BMI Quartile',
     ylab = 'Log Odds',
     main = 'Log Odds by BMI Quartile',
     type = "b", pch = 19, lty = 1)

plot(as.numeric(plotdata_systol$Quartile), plotdata_systol$LogOdds,
     xlab = 'Systolic BP Quartile',
     ylab = 'Log Odds',
     main = 'Log Odds by Systolic BP Quartile',
     type = "b", pch = 19, lty = 1)

par(mfrow = c(1, 1))  # reset plotting layout

## STEP 6: Test Interaction Effects
# Fit models with pairwise interactions to assess if effects of predictors depend on each other.

maineff <- glm(chd ~ totchol + bmi + systol + smoking_c, family = binomial(logit), data = mydata1)

int1 <- glm(chd ~ totchol + bmi + systol + smoking_c + totchol * bmi, family = binomial(logit), data = mydata1)
int2 <- glm(chd ~ totchol + bmi + systol + smoking_c + totchol * systol, family = binomial(logit), data = mydata1)
int3 <- glm(chd ~ totchol + bmi + systol + smoking_c + totchol * smoking_c, family = binomial(logit), data = mydata1)
int4 <- glm(chd ~ totchol + bmi + systol + smoking_c + bmi * systol, family = binomial(logit), data = mydata1)
int5 <- glm(chd ~ totchol + bmi + systol + smoking_c + bmi * smoking_c, family = binomial(logit), data = mydata1)
int6 <- glm(chd ~ totchol + bmi + systol + smoking_c + systol * smoking_c, family = binomial(logit), data = mydata1)

# Calculate likelihood ratio test statistics and p-values comparing main effects model to interaction models.
lrt_stats <- function(main_model, interaction_model) {
  dev_diff <- main_model$deviance - interaction_model$deviance
  df_diff <- main_model$df.residual - interaction_model$df.residual
  p_val <- pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
  return(c(Deviance_Diff = dev_diff, DF_Diff = df_diff, p_value = p_val))
}

lrt_results <- rbind(
  int1 = lrt_stats(maineff, int1),
  int2 = lrt_stats(maineff, int2),
  int3 = lrt_stats(maineff, int3),
  int4 = lrt_stats(maineff, int4),
  int5 = lrt_stats(maineff, int5),
  int6 = lrt_stats(maineff, int6)
)

STEP6 <- data.frame(
  Interaction = c("totchol*bmi", "totchol*systol", "totchol*smoking_c", "bmi*systol", "bmi*smoking_c", "systol*smoking_c"),
  Deviance_Change = lrt_results[, "Deviance_Diff"],
  DF_Change = lrt_results[, "DF_Diff"],
  p_value = lrt_results[, "p_value"]
)

print(STEP6)


## STEP 7: Final Model Selection and Model Fit Assessment
# Build a preliminary final model including main effects and significant interactions.
intmodel1 <- glm(chd ~ totchol + bmi + systol + smoking_c + totchol*smoking_c + bmi*systol + bmi*smoking_c, 
                 family = binomial(logit), data = mydata1)

summary(intmodel1)

# Likelihood ratio test comparing main effects model to interaction model.
lrt_g <- maineff$deviance - intmodel1$deviance
lrt_d <- maineff$df.residual - intmodel1$df.residual
lrt_p <- pchisq(lrt_g, df = lrt_d, lower.tail = FALSE)

cat("LRT Chi-square:", lrt_g, "\nDegrees of freedom:", lrt_d, "\np-value:", lrt_p, "\n")

# Consider a simplified model by removing some interactions and compare.
intmodel2 <- glm(chd ~ totchol + bmi + systol + smoking_c + bmi*smoking_c, family = binomial(logit), data = mydata1)
summary(intmodel2)

lrt_g2 <- maineff$deviance - intmodel2$deviance
lrt_d2 <- maineff$df.residual - intmodel2$df.residual
lrt_p2 <- pchisq(lrt_g2, df = lrt_d2, lower.tail = FALSE)

cat("LRT Chi-square (simplified):", lrt_g2, "\nDegrees of freedom:", lrt_d2, "\np-value:", lrt_p2, "\n")

# Compare the two interaction models to assess parsimony.
comp_g <- intmodel2$deviance - intmodel1$deviance
comp_d <- intmodel2$df.residual - intmodel1$df.residual
comp_p <- pchisq(comp_g, df = comp_d, lower.tail = FALSE)

cat("Comparison Chi-square:", comp_g, "\nDegrees of freedom:", comp_d, "\np-value:", comp_p, "\n")

# Perform Hosmer-Lemeshow goodness-of-fit test on the final selected model.
hoslem.test(mydata1$chd, fitted(intmodel2), g = 10)


## Step 8. Save Results 

# Make sure output directories exist
output_dir <- "../output"
tables_dir <- file.path(output_dir, "tables")

if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(tables_dir)) dir.create(tables_dir)

# 1: Univariable logistic regression results
write.csv(STEP1, file = file.path(tables_dir, "STEP1_univariable_results.csv"), row.names = FALSE)

# 2: Multivariable logistic regression coefficients summary (extracting coefficients table)
step2_summary <- summary(step2model)
step2_coefs <- as.data.frame(coef(step2_summary))
step2_coefs$Variable <- rownames(step2_coefs)
rownames(step2_coefs) <- NULL
write.csv(step2_coefs, file = file.path(tables_dir, "STEP2_multivariable_coefficients.csv"), row.names = FALSE)

# 3: Confounding assessment tables
write.csv(confounding1, file = file.path(tables_dir, "STEP3_confounding_excluding_age.csv"), row.names = FALSE)
write.csv(confounding2, file = file.path(tables_dir, "STEP3_confounding_excluding_bmi.csv"), row.names = FALSE)
write.csv(STEP3, file = file.path(tables_dir, "STEP3_confounding_combined.csv"), row.names = FALSE)

# 4: Wald test results for adding back excluded variables
write.csv(wald_results, file = file.path(tables_dir, "STEP4_wald_tests.csv"), row.names = FALSE)

# 5: Functional form assessment - log odds by quartiles for each variable
write.csv(plotdata_totchol, file = file.path(tables_dir, "STEP5_logodds_totchol_quartiles.csv"), row.names = FALSE)
write.csv(plotdata_bmi, file = file.path(tables_dir, "STEP5_logodds_bmi_quartiles.csv"), row.names = FALSE)
write.csv(plotdata_systol, file = file.path(tables_dir, "STEP5_logodds_systol_quartiles.csv"), row.names = FALSE)

# 6: Interaction effects - Likelihood ratio test results
write.csv(STEP6, file = file.path(tables_dir, "STEP6_interaction_tests.csv"), row.names = FALSE)

# 7: Final models coefficient summaries
intmodel1_summary <- summary(intmodel1)
intmodel1_coefs <- as.data.frame(coef(intmodel1_summary))
intmodel1_coefs$Variable <- rownames(intmodel1_coefs)
rownames(intmodel1_coefs) <- NULL
write.csv(intmodel1_coefs, file = file.path(tables_dir, "STEP7_final_model_full_interaction_coefficients.csv"), row.names = FALSE)

intmodel2_summary <- summary(intmodel2)
intmodel2_coefs <- as.data.frame(coef(intmodel2_summary))
intmodel2_coefs$Variable <- rownames(intmodel2_coefs)
rownames(intmodel2_coefs) <- NULL
write.csv(intmodel2_coefs, file = file.path(tables_dir, "STEP7_final_model_simplified_interaction_coefficients.csv"), row.names = FALSE)

# Optionally, save a small summary table for model comparison stats
model_comparison <- data.frame(
  Comparison = c("Main vs Full Interaction", "Main vs Simplified Interaction", "Full vs Simplified Interaction"),
  Chi_Square = c(lrt_g, lrt_g2, comp_g),
  DF = c(lrt_d, lrt_d2, comp_d),
  p_value = c(lrt_p, lrt_p2, comp_p)
)
write.csv(model_comparison, file = file.path(tables_dir, "STEP7_model_comparisons.csv"), row.names = FALSE)

plots_dir <- file.path(output_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Save the Total Cholesterol log-odds plot
png(filename = file.path(plots_dir, "logodds_totchol_quartiles.png"), width = 600, height = 400)
plot(as.numeric(plotdata_totchol$Quartile), plotdata_totchol$LogOdds,
     xlab = 'Total Cholesterol Quartile',
     ylab = 'Log Odds',
     main = 'Log Odds by Total Cholesterol Quartile',
     type = "b", pch = 19, lty = 1)
dev.off()

# Save the BMI log-odds plot
png(filename = file.path(plots_dir, "logodds_bmi_quartiles.png"), width = 600, height = 400)
plot(as.numeric(plotdata_bmi$Quartile), plotdata_bmi$LogOdds,
     xlab = 'BMI Quartile',
     ylab = 'Log Odds',
     main = 'Log Odds by BMI Quartile',
     type = "b", pch = 19, lty = 1)
dev.off()

# Save the Systolic BP log-odds plot
png(filename = file.path(plots_dir, "logodds_systol_quartiles.png"), width = 600, height = 400)
plot(as.numeric(plotdata_systol$Quartile), plotdata_systol$LogOdds,
     xlab = 'Systolic BP Quartile',
     ylab = 'Log Odds',
     main = 'Log Odds by Systolic BP Quartile',
     type = "b", pch = 19, lty = 1)
dev.off()



#----------------------------------------------------------------------------------
# COMPARISON OF LOGISTIC, LOG-BINOMIAL, AND POISSON REGRESSION MODELS FOR CHD RISK
#-----------------------------------------------------------------------------------

# Objective:
# To estimate and compare the association between predictors (total cholesterol, BMI, smoking status)
# and coronary heart disease (CHD) using three regression models:
# 1) Logistic regression (estimates odds ratios, OR)
# 2) Log-binomial regression (directly estimates relative risks, RR)
# 3) Poisson regression with robust standard errors (alternative method to estimate RR)
# The comparison includes effect estimates and confidence intervals from each model.

library(sandwich)  # for robust standard errors
library(lmtest)    # for coefficient testing with robust SEs

# Load and prepare data
mydata2 <- read.csv("................................./data/SHHS.csv", header = TRUE)
mydata2$smoking_c <- as.factor(mydata2$smoking)

# Fit logistic regression model (logit link)
logit_model <- glm(chd ~ totchol + bmi + smoking_c, family = binomial(link = "logit"), data = mydata2)

# Fit log-binomial regression model (log link)
logbin_model <- glm(chd ~ totchol + bmi + smoking_c, family = binomial(link = "log"), data = mydata2)

# Fit Poisson regression model with log link and robust SEs for RR estimation
poisson_model <- glm(chd ~ totchol + bmi + smoking_c, family = poisson(link = "log"), data = mydata2)
cov_poisson_robust <- vcovHC(poisson_model, type = "HC0")
poisson_robust_se <- sqrt(diag(cov_poisson_robust))
poisson_coefs_robust <- coeftest(poisson_model, vcov = cov_poisson_robust)

# Helper function to extract exponentiated estimates and 95% CI
extract_estimates <- function(model) {
  est <- coef(model)
  ci <- confint.default(model)
  data.frame(
    Estimate = round(exp(est), 3),
    Lower_CI = round(exp(ci[,1]), 3),
    Upper_CI = round(exp(ci[,2]), 3)
  )
}

# Extract and label results for each model
logit_results <- extract_estimates(logit_model)
names(logit_results) <- paste0("Logit_", names(logit_results))

logbin_results <- extract_estimates(logbin_model)
names(logbin_results) <- paste0("LogBin_", names(logbin_results))

poisson_est <- coef(poisson_model)
poisson_est_exp <- exp(poisson_est)
poisson_lower_ci <- exp(poisson_est - 1.96 * poisson_robust_se)
poisson_upper_ci <- exp(poisson_est + 1.96 * poisson_robust_se)

poisson_results <- data.frame(
  Poisson_Estimate = round(poisson_est_exp, 3),
  Poisson_SE = round(poisson_robust_se, 3),
  Poisson_Lower_CI = round(poisson_lower_ci, 3),
  Poisson_Upper_CI = round(poisson_upper_ci, 3)
)

# Combine all model results into one table for comparison
results_table <- cbind(
  Variable = names(coef(logit_model)),
  logit_results,
  logbin_results,
  poisson_results
)

# Print the combined results
print(results_table)

# Define output directory path (the folder where you want to save tables)
output_tables_dir <- file.path(output_dir, "tables")  

# Create the directory if it doesn't exist
if (!dir.exists(output_tables_dir)) {
  dir.create(output_tables_dir, recursive = TRUE)
}

# Save the results table as CSV
write.csv(results_table, file.path(output_tables_dir, "CHD_Model_Comparison_Results.csv"), row.names = FALSE)

# Confirmation message
cat("Model comparison results saved to", file.path(output_tables_dir, "CHD_Model_Comparison_Results.csv"), "\n")



#-------------------------------------------------------------------------------------
# ASSESSING CURRENT SMOKING AS A RISK FACTOR FOR CHD USING PROPENSITY SCORE WEIGHTING
#-------------------------------------------------------------------------------------

# Objective:
# To evaluate whether being a current smoker increases the risk of coronary heart disease (CHD),
# controlling for potential confounders (age, total cholesterol, BMI, systolic blood pressure, and activity level)
# using propensity score weighting (IPTW) without matching.

library(broom)   # for augment() function
library(car)     # for Anova()

# Load and inspect data
mydata3 <- read.csv("................................./data/SHHS.csv", header = TRUE)
str(mydata3)

# Define current smoker variable: 1 = current smoker, 0 = non-current smoker
table(mydata3$smoking)
mydata3$current_smoker <- ifelse(mydata3$smoking == 3, 1, 0)
mydata3$activity <- as.factor(mydata3$activity)  # Factorize activity for modeling

# Fit propensity score model predicting current smoking from confounders
ps_model <- glm(current_smoker ~ age + totchol + bmi + systol + activity,
                family = binomial(link = "logit"), data = mydata3)
summary(ps_model)

# Extract propensity scores (probability of being a current smoker)
ps_model.results <- augment(ps_model, data = mydata3)
ps_model.results$p <- exp(ps_model.results$.fitted) / (1 + exp(ps_model.results$.fitted))

# Calculate inverse probability of treatment weights (IPTW)
ps_model.results$wgt <- ifelse(ps_model.results$current_smoker == 1,
                               1 / ps_model.results$p,
                               1 / (1 - ps_model.results$p))

# Summarize IPTW weights distribution and save histogram
summary(ps_model.results$wgt)

# Define output folders (update the path as per your system)
output_dir <- "................................./SHHS/output"
tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")

# Create folders if they do not exist
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(tables_dir)) dir.create(tables_dir)
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# --- Save Tables ---

# Save combined logistic, log-binomial, and poisson regression results table
write.csv(results_table, file.path(tables_dir, "CHD_Model_Comparison_Results.csv"), row.names = FALSE)

# Save unweighted logistic regression model results (from propensity score section)
write.csv(unweighted_results, file.path(tables_dir, "CHD_unweighted_model_results.csv"), row.names = FALSE)

# Save weighted logistic regression model results (from propensity score section)
write.csv(weighted_results, file.path(tables_dir, "CHD_weighted_model_results.csv"), row.names = FALSE)

# --- Save Plot ---

# Save IPTW weights histogram plot to output/plots folder
png(filename = file.path(plots_dir, "IPTW_weights_histogram.png"), width = 800, height = 600)
hist(ps_model.results$wgt, breaks = 50, main = "Histogram of IPTW Weights", xlab = "Weights")
dev.off()

# Informative messages
cat("Tables saved to:", tables_dir, "\n")
cat("Plots saved to:", plots_dir, "\n")
