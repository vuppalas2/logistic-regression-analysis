#########################################
# NABCR Logistic Regression Analysis
########################################

# 1. Load Required Libraries -----------------------------------
library(car)        # For Anova on glm models
library(lmtest)     # For likelihood ratio tests (lrtest)
library(broom)      # For tidy model outputs
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting

# 2. Load Data -------------------------------------------------
nabcr_data <- read.csv(".........../data/NABCR.csv")

# Inspect data structure and summary
str(nabcr_data)
summary(nabcr_data)

# 3. Data Preparation -------------------------------------------
# Convert variables to factors with meaningful labels
nabcr_data <- nabcr_data %>%
  mutate(
    level_c = factor(level, levels = c(1, 2), labels = c("Low", "High")),
    stage_c = factor(stage, levels = c(1, 2, 3), labels = c("Stage I", "Stage II", "Stage III")),
    dead_c = factor(dead, levels = c(0, 1), labels = c("Alive", "Dead"))
  )

# 4. Descriptive Statistics -------------------------------------
# Frequencies of each variable
cat("Frequency of Receptor Level:\n")
print(table(nabcr_data$level_c))

cat("\nFrequency of Survival Status:\n")
print(table(nabcr_data$dead_c))

cat("\nFrequency of Cancer Stage:\n")
print(table(nabcr_data$stage_c))

# Cross-tabulations to see distributions and relationships
cat("\nReceptor Level vs Survival Status:\n")
print(table(nabcr_data$level_c, nabcr_data$dead_c))

cat("\nCancer Stage vs Survival Status:\n")
print(table(nabcr_data$stage_c, nabcr_data$dead_c))

cat("\nReceptor Level vs Cancer Stage:\n")
print(table(nabcr_data$level_c, nabcr_data$stage_c))

# Optional: Stratified receptor level vs survival by cancer stage
for(stage_level in levels(nabcr_data$stage_c)) {
  cat(sprintf("\nReceptor Level vs Survival Status for %s:\n", stage_level))
  print(table(
    subset(nabcr_data, stage_c == stage_level)$level_c,
    subset(nabcr_data, stage_c == stage_level)$dead_c
  ))
}

# 5. Logistic Regression Models ---------------------------------

# Model 1: Unadjusted effect of receptor level on death
model1 <- glm(dead ~ level_c, family = binomial(logit), data = nabcr_data)
summary(model1)

# Extract Odds Ratio (OR) and 95% Confidence Interval (CI)
OR1 <- exp(coef(model1)[2])
CI1 <- exp(confint(model1)[2])
cat("Model 1 (Unadjusted) OR:", OR1, "\n95% CI:", CI1, "\n\n")

# Model 2: Adjusted for cancer stage to assess confounding
model2 <- glm(dead ~ level_c + stage_c, family = binomial(logit), data = nabcr_data)
summary(model2)

OR2 <- exp(coef(model2)[2])
CI2 <- exp(confint(model2)[2])
cat("Model 2 (Adjusted) OR:", OR2, "\n95% CI:", CI2, "\n\n")

# Likelihood ratio test comparing Model 1 and Model 2
lrtest_res <- lrtest(model2, model1)
print(lrtest_res)

# Calculate percent change in OR to quantify confounding
confounding_pct <- ((OR1 - OR2) / OR2) * 100
cat("Percent change in OR due to confounding by stage:", round(confounding_pct, 2), "%\n\n")

# Model 3: Add interaction term to test effect modification
model3 <- glm(dead ~ level_c * stage_c, family = binomial(logit), data = nabcr_data)
summary(model3)

# Likelihood ratio test comparing Model 2 and Model 3
lrtest_interaction <- lrtest(model3, model2)
print(lrtest_interaction)

# 6. Stage-Specific Logistic Regression ------------------------

stage_ORs <- data.frame(Stage = character(),
                        OR = numeric(),
                        Lower_CI = numeric(),
                        Upper_CI = numeric(),
                        stringsAsFactors = FALSE)

for (stage_level in levels(nabcr_data$stage_c)) {
  cat("Running logistic regression for", stage_level, "...\n")
  data_subset <- filter(nabcr_data, stage_c == stage_level)
  mod <- glm(dead ~ level_c, family = binomial(logit), data = data_subset)
  
  OR <- exp(coef(mod)[2])
  CI <- exp(confint(mod)[2, ])
  
  cat("OR:", OR, "\n95% CI:", CI, "\n\n")
  
  stage_ORs <- rbind(stage_ORs, data.frame(
    Stage = stage_level,
    OR = OR,
    Lower_CI = CI[1],
    Upper_CI = CI[2]
  ))
}

# 7. Plot Stage-specific Odds Ratios ---------------------------
ggplot(stage_ORs, aes(x = Stage, y = OR)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Adjusted Odds Ratios of Death by Receptor Level Across Cancer Stages",
       x = "Cancer Stage",
       y = "Odds Ratio (95% CI)") +
  theme_minimal()

# 8. Save Results ------------------------------------------------

output_dir <- "../output"
tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")

# Create folders if not exist
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(tables_dir)) dir.create(tables_dir)
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Then save files with these paths:
write.csv(stage_ORs, file.path(tables_dir, "stage_specific_ORs.csv"), row.names = FALSE)
write.csv(broom::tidy(model1, exponentiate = TRUE, conf.int = TRUE),
          file.path(tables_dir, "model1_summary.csv"), row.names = FALSE)
write.csv(broom::tidy(model2, exponentiate = TRUE, conf.int = TRUE),
          file.path(tables_dir, "model2_summary.csv"), row.names = FALSE)
write.csv(broom::tidy(model3, exponentiate = TRUE, conf.int = TRUE),
          file.path(tables_dir, "model3_summary.csv"), row.names = FALSE)

ggsave(file.path(plots_dir, "stage_specific_ORs_plot.png"), width = 6, height = 4)


