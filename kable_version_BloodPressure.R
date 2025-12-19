# Comparative Antihypertensive Efficacy of Ramipril, Lisinopril, and Moexipril
# Rongrong Qian | 10 Sep 2025

# Load required packages
library(readr)       # read CSV data
library(dplyr)       # data manipulation
library(tidyr)       # data reshaping
library(ggplot2)     # plotting
library(knitr)       # tables
library(kableExtra)  # table formatting
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

knitr::kable(
  summary_stats_formatted,
  caption = "Table 1: Descriptive Summary of SBP (mmHg) by Medication",
  booktabs = TRUE,
  align = 'lccrrrrr'
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

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

# ggplot for Individual SBP Change
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

# Jitter plots of absolute and percent reduction
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


# Scatter plots with regression lines for baseline SBP vs. SBP change
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

# plot a scatter plot without obs 6
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
anova_int <- Anova(fit_int, type = 2)
knitr::kable(anova_int, caption = "Table 2: ANOVA for Interaction Model") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# ANCOVA model without interaction
fit_ancova <- lm(after ~ before + drug, data = bp)

# Residual diagnostics
# Levene's test for homogeneity of variance
levene_test_result <- leveneTest(residuals(fit_ancova) ~ bp$drug)
knitr::kable(as.data.frame(levene_test_result),
             caption = "Table 3: Levene's Test for Homogeneity of Variance",
             digits = 3, booktabs = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# Shapiro-Wilk test for normality of residuals
shapiro_test_result <- shapiro.test(residuals(fit_ancova))
knitr::kable(
  tibble(Statistic = c("W", "p-value"),
         Value = c(round(shapiro_test_result$statistic, 3),
                   round(shapiro_test_result$p.value, 3))),
  caption = "Table 4: Shapiro-Wilk Normality Test"
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# Cook's distance for influential points
cooks_d <- cooks.distance(fit_ancova)
influential_points <- which(cooks_d > (4 / nrow(bp)))
knitr::kable(
  tibble(
    Subject = names(cooks_d),
    Cook_s_Distance = round(cooks_d, 3),
    Influential = ifelse(names(cooks_d) %in% names(influential_points),
                         "Yes", "No")
  ),
  caption = "Table 5: Cook's Distance for Influential Observations",
  booktabs = TRUE
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# Model coefficients with 95% CI
coef_ci <- confint(fit_ancova)
tidy(fit_ancova) %>%
  mutate(
    CI_Lower = round(coef_ci[, 1], 3),
    CI_Upper = round(coef_ci[, 2], 3)
  ) %>%
  knitr::kable(
    caption = "Table 6: ANCOVA Coefficients with 95% CI",
    digits = 3, booktabs = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# Model summary stats and ANOVA
summary_ancova <- summary(fit_ancova)
knitr::kable(
  tibble(
    Statistic = c("R-squared", "Adj R-squared", "F-statistic", "RMSE"),
    Value = c(round(summary_ancova$r.squared, 3),
              round(summary_ancova$adj.r.squared, 3),
              round(summary_ancova$fstatistic[1], 2),
              round(summary_ancova$sigma, 2))
  ),
  caption = "Table 7: ANCOVA Model Summary"
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

anova_ancova <- Anova(fit_ancova, type = 2)
knitr::kable(anova_ancova,
             caption = "Table 8: Type-II ANOVA for ANCOVA Model",
             digits = 3, booktabs = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# Estimated marginal means with 95% CI and pairwise comparisons
emm <- emmeans(fit_ancova, ~ drug)

knitr::kable(
  summary(emm, infer = TRUE),   # infer = TRUE gives 95% CI by default
  caption = "Table 9: Adjusted Means of Post-treatment SBP with 95% CI",
  digits = 2, booktabs = TRUE
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

# Pairwise comparisons with Tukey adjustment
pairwise_results <- emmeans(fit_ancova, pairwise ~ drug, adjust = "tukey")

knitr::kable(
  tidy(pairwise_results$contrasts),
  caption = "Table 10: Pairwise Comparisons of Adjusted Means",
  digits = 3, booktabs = TRUE
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")


# -------------------------------------------------------
# 5. Randomization check using permutation test
# -------------------------------------------------------
# Rank-based permutation test for baseline SBP differences
rankdata <- bp %>%
  mutate(BP.pre.rank = rank(before)) %>%
  select(drug, before, BP.pre.rank)

# Observed F-statistic from ANOVA on ranks
raw.results <- lm(BP.pre.rank ~ drug, data = rankdata)
Fobs <- summary(raw.results)$f[1]

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

# Permutation p-value
pvalue <- mean(c(Fobs, reps) >= Fobs)

# Display results
knitr::kable(
  tibble(
    Statistic = c("Observed F-statistic", "Permutation p-value"),
    Value = c(round(Fobs, 3), round(pvalue, 3))
  ),
  caption = "Table 11: Permutation Test for Baseline SBP Differences",
  booktabs = TRUE
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")