# Comparative Antihypertensive Efficacy of Ramipril, Lisinopril, and Moexipril
# Rongrong Qian | 10 Sep 2025

# Load required packages
library(readr)       # read CSV data
library(dplyr)       # data manipulation
library(tidyr)       # data reshaping
library(ggplot2)     # plotting
library(patchwork)   # combine plots
library(broom)       # tidy model outputs
library(car)         # type-II ANOVA, Levene’s test
library(emmeans)     # estimated marginal means

# -------------------------------------------------------
# 1. Read and prepare data
# -------------------------------------------------------
# Please change the file path to your local path
bp <- read_csv("/Users/rongrongqian/Desktop/660/1/bloodpressurestudy.csv", 
               show_col_types = FALSE)

# Calculate absolute and percent SBP change
bp <- bp %>%
  mutate(
    SBP_Change = before - after,
    SBP_Change_Pct = (before - after) / before * 100
  )

# -------------------------------------------------------
# 2. Descriptive statistics by medication
# -------------------------------------------------------
# Summary statistics table
summary_stats <- bp %>%
  group_by(drug) %>%
  summarise(
    N = n(),
    Baseline_Mean = mean(before),
    Baseline_SD = sd(before),
    Post_Mean = mean(after),
    Post_SD = sd(after),
    Reduction_Mean = mean(SBP_Change),
    Reduction_SD = sd(SBP_Change),
    Reduction_Pct_Mean = mean(SBP_Change_Pct),
    .groups = "drop"
  )

# Format and display the summary statistics table
summary_stats_formatted <- summary_stats %>%
  mutate(across(where(is.numeric) & !N, ~ round(., 1))) %>%
  rename(
    Medication = drug,
    "Baseline Mean" = Baseline_Mean,
    "Baseline SD" = Baseline_SD,
    "Post-treatment Mean" = Post_Mean,
    "Post-treatment SD" = Post_SD,
    "Mean Reduction" = Reduction_Mean,
    "Reduction SD" = Reduction_SD,
    "Mean % Reduction" = Reduction_Pct_Mean
  )

cat("\n========== Descriptive Statistics by Medication ==========\n")
summary_stats_formatted


# -------------------------------------------------------
# 3. EDA
# -------------------------------------------------------
# Convert data from 'wide' to 'long' format for plotting
bp_long <- bp %>%
  pivot_longer(
    cols = c(before, after),
    names_to = "time",
    values_to = "sbp"
  ) %>%
  # Ensure 'time' is an ordered factor for correct line plotting
  mutate(time = factor(time, levels = c("before", "after")))

# (Fig 1) Individual lines by drug
ggplot(data = bp_long, aes(x = time, y = sbp, 
                           group = subject, color = drug)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  # Create separate panels for each medication
  facet_wrap(~ drug) +
  labs(
    title = "Figure 1: Individual SBP Change Before and After Treatment",
    subtitle = "Each line represents a single patient",
    x = "Measurement Time",
    y = "Systolic Blood Pressure (mmHg)",
    color = "Medication"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

# (Fig 2) Jitter plots of absolute and percent reduction
# Plot (a): Absolute SBP Change
plot_a <- ggplot(data = bp, aes(x = drug, y = SBP_Change, color = drug)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "errorbar", 
               aes(x = drug, y = SBP_Change, color = drug),
               width = 0.3, color = "black") +
  labs(
    title = "Figure 2: Jitter Plots of SBP Change by Medication",
    subtitle = "(a) Absolute SBP Change",
    x = NULL,
    y = "Absolute SBP Change (mmHg)",
    color = "Medication"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Plot (b): Percent SBP Change
plot_b <- ggplot(data = bp, aes(x = drug, y = SBP_Change_Pct, color = drug)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "errorbar", 
               aes(x = drug, y = SBP_Change_Pct, color = drug),
               width = 0.3, color = "black") +
  labs(
    subtitle = "(b) Percent SBP Change",
    x = NULL,
    y = "Percent Reduction (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Combine the two plots side-by-side
plot_a + plot_b

# (Fig 3) Scatter plots with regression lines for baseline SBP vs. SBP change
ggplot(data = bp, 
       aes(x = before, y = SBP_Change, color = drug, fill = drug)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(
    title = "Figure 3: Relationship Between Baseline SBP and SBP Reduction",
    subtitle = "Scatter plot with regression line",
    x = "Baseline SBP (mmHg)",
    y = "Absolute SBP Change (mmHg)",
    color = "Medication",
    fill = "Medication"
  ) +
  theme_minimal()

# (Fig 4) plot a scatter plot without obs 6
graph_a <- ggplot(data = bp, 
                  aes(x = before, y = after, color = drug, fill = drug)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(
    title = "Figure 4: Relationship Between Baseline SBP and Post-treatment SBP",
    subtitle = "(a) Scatter plot with regression line",
    x = "Baseline SBP (mmHg)",
    y = "Post-treatment SBP (mmHg)",
    color = "Medication", 
    fill = "Medication"
  ) +
  theme_minimal()

# plot a scatter plot without obs 6
graph_b <- ggplot(data = bp %>% filter(subject != 6), 
                  aes(x = before, y = after, color = drug, fill = drug)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(
    subtitle = "(b) Scatter plot leave observation 6 out",
    x = "Baseline SBP (mmHg)",
    y = "Post-treatment SBP (mmHg)",
    color = "Medication",   
    fill = "Medication"
  ) +
  theme_minimal()

# Combine the two plots side-by-side with shared legend
(graph_a | graph_b) + plot_layout(guides = "collect") & theme(legend.position = "right")


# -------------------------------------------------------
# 4. ANCOVA model
# -------------------------------------------------------
# Ensure drug is a factor before releveling
bp$drug <- factor(bp$drug)

# Set Lisinopril as the reference group
bp$drug <- relevel(bp$drug, ref = "Lisinopril")

# Interaction model: check homogeneity of slopes
fit_int <- lm(after ~ before * drug, data = bp)
cat("\n========== Interaction Model (after ~ before * drug) ==========\n")
summary(fit_int)
cat("\nType-II ANOVA (Interaction Model):\n")
print(car::Anova(fit_int, type = 2))

# ANCOVA model without interaction
fit_ancova <- lm(after ~ before + drug, data = bp)
cat("\n========== ANCOVA (after ~ before + drug) ==========\n")
summary(fit_ancova)


# Residual diagnostics
# Levene's test for homogeneity of variance
cat("\n========== Residual Diagnostics ==========\n")
levene_test_result <- leveneTest(residuals(fit_ancova) ~ bp$drug)
cat("\nLevene's Test for Homogeneity of Variance (on residuals by drug):\n")
print(levene_test_result)

# Shapiro-Wilk test for normality of residuals
shapiro_test_result <- shapiro.test(residuals(fit_ancova))
cat("\nShapiro-Wilk Test for Normality of Residuals:\n")
print(shapiro_test_result)

# Cook's distance for influential points
cooks_d <- cooks.distance(fit_ancova)
infl_idx <- which(cooks_d > (4 / nrow(bp)))
cat("\nCook's Distance (threshold = 4/n):\n")
cd_tab <- data.frame(
  Row = seq_along(cooks_d),
  Subject = if ("subject" %in% names(bp)) bp$subject else seq_len(nrow(bp)),
  CooksD = round(cooks_d, 4),
  Influential = ifelse(seq_along(cooks_d) %in% infl_idx, "Yes", "No")
)
print(cd_tab)
if (length(infl_idx)) {
  cat("\nPotentially influential rows:", paste(infl_idx, collapse = ", "), "\n")
} else {
  cat("\nNo influential points flagged by threshold.\n")
}

# Model summary stats
sa <- summary(fit_ancova)
cat("\n========== Model Summary Stats (ANCOVA) ==========\n")
cat(sprintf("R-squared      : %.3f\n", sa$r.squared))
cat(sprintf("Adj R-squared  : %.3f\n", sa$adj.r.squared))
cat(sprintf("F-statistic    : %.2f on %d and %d DF (p = %.4g)\n",
            sa$fstatistic["value"], sa$fstatistic["numdf"], sa$fstatistic["dendf"],
            pf(sa$fstatistic["value"], sa$fstatistic["numdf"], sa$fstatistic["dendf"], lower.tail = FALSE)))
cat(sprintf("RMSE (sigma)   : %.3f\n", sigma(fit_ancova)))

# Model coefficients with 95% CI
cat("\n========== Coefficients with 95% CI (ANCOVA) ==========\n")
coef_tab <- broom::tidy(fit_ancova)
coef_tab$CI_Lower <- round(coef_ci[,1], 4)
coef_tab$CI_Upper <- round(coef_ci[,2], 4)
print(coef_tab)

# Type-II ANOVA table
anova_ancova <- Anova(fit_ancova, type = 2)
cat("\nType-II ANOVA (ANCOVA Model):\n")
anova_ancova

# Estimated marginal means with 95% CI and pairwise comparisons
cat("\n========== Estimated Marginal Means (by drug) ==========\n")
emm <- emmeans(fit_ancova, ~ drug)
summary(emm, infer = TRUE)

# Pairwise comparisons with Tukey adjustment
cat("\nPairwise Comparisons (Tukey-adjusted):\n")
pairwise_results <- emmeans(fit_ancova, pairwise ~ drug, adjust = "tukey")
summary(pairwise_results$contrasts, infer = TRUE)


# -------------------------------------------------------
# 5. Randomization check using permutation test
# -------------------------------------------------------
cat("\n========== Randomization Check via Permutation Test ==========\n")
# Rank-based permutation test for baseline SBP differences
rankdata <- bp %>%
  mutate(BP.pre.rank = rank(before)) %>%
  select(drug, before, BP.pre.rank)

# Observed F-statistic from ANOVA on ranks
raw.results <- lm(BP.pre.rank ~ drug, data = rankdata)
Fobs <- summary(raw.results)$f[1]
Fobs

# Permutation procedure
set.seed(123)
reps <- replicate(9999, {
  perm.data <- data.frame(
    drug = rankdata$drug,
    permuted.rank = sample(rankdata$BP.pre.rank)
  )
  summary(lm(permuted.rank ~ drug, data = perm.data))$f[1]
})

# Plot null distribution
ggplot(data.frame(reps), aes(x = reps)) +
  geom_histogram() +
  geom_vline(xintercept = Fobs, color = "red", linetype = "longdash") +
  labs(title = "Null Distribution of Permuted F-statistic",
       x = "Simulated F values under H0")

# p-value
pvalue <- mean(c(Fobs, reps) >= Fobs)
cat("The observed TS is F =", Fobs,"; p-value =", pvalue, "\n")