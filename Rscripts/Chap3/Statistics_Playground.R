## playground for statistics
# Load necessary libraries
# If you don't have them, install via install.packages(c("tidyverse", "lme4", "lmerTest", "emmeans"))
library(tidyverse)  # For data manipulation and plotting
library(lme4)       # For Linear Mixed Models
library(lmerTest)   # To get p-values for Mixed Models
library(emmeans)    # For pairwise comparisons (post-hoc tests)

# ==========================================
# 1. GENERATE DUMMY DATA (Replace this with your actual data)
# ==========================================
set.seed(123)
n_samples <- 18
methods <- c("Raw", "Method_X", "Method_Y")

# Create a simulated dataset
data_sim <- expand.grid(SampleID = factor(1:n_samples), Method = methods) %>%
  mutate(
    # Simulate Column A (Significant Species Count)
    # Raw has most, Methods X and Y might reduce it slightly
    Col_A_SigSpecies = case_when(
      Method == "Raw" ~ rpois(n(), 150),
      Method == "Method_X" ~ rpois(n(), 145), # X preserves well
      Method == "Method_Y" ~ rpois(n(), 100)  # Y removes too much
    ),
    
    # Simulate Column B (Contaminant Intersection)
    # Raw has contaminants, X and Y reduce them
    Col_B_Contaminants = case_when(
      Method == "Raw" ~ rpois(n(), 20),
      Method == "Method_X" ~ rpois(n(), 2),   # X works well
      Method == "Method_Y" ~ rpois(n(), 1)    # Y works well but maybe too aggressive
    )
  ) %>%
  # Calculate the Proportion of Contaminants remaining
  mutate(Prop_Contamination = Col_B_Contaminants / Col_A_SigSpecies)

print("Head of the dataset:")
print(head(data_sim))

# ==========================================
# 2. VISUALIZATION
# ==========================================
# We use boxplots with connecting lines to show the paired nature of the data

plot_metric <- function(data, y_var, title) {
  ggplot(data, aes(x = Method, y = .data[[y_var]], fill = Method)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(alpha = 0.6) +
    # Add lines connecting the same sample across methods (Critical for paired view)
    geom_line(aes(group = SampleID), alpha = 0.2, color = "gray30") +
    theme_minimal() +
    labs(title = title, y = "Count/Ratio") +
    theme(legend.position = "none")
}

p1 <- plot_metric(data_sim, "Col_A_SigSpecies", "Preservation: Total Significant Species (A)")
p2 <- plot_metric(data_sim, "Col_B_Contaminants", "Efficacy: Remaining Contaminants (B)")
p3 <- plot_metric(data_sim, "Prop_Contamination", "Relative: Proportion of Contaminants (B/A)")

# Display plots
print(p1)
print(p2)
print(p3)

# ==========================================
# 3. STATISTICAL MODELING (Linear Mixed Models)
# ==========================================

# We use (1|SampleID) to account for the fact that the 18 samples are repeated measures.

# --- Model 1: Comparing Contaminant Counts (Col B) ---
# We want to see which method yields the lowest B.
# Note: counts often require log transformation or Poisson GLMM. 
# Here we use log(x+1) to handle zeros and normalize residuals.
mod_contam <- lmer(log(Col_B_Contaminants + 1) ~ Method + (1|SampleID), data = data_sim)

cat("\n=== ANOVA Results for Contaminant Removal ===\n")
print(anova(mod_contam))

# Post-hoc Pairwise comparisons
cat("\n=== Pairwise Comparisons (Contaminants) ===\n")
emm_contam <- emmeans(mod_contam, specs = pairwise ~ Method)
print(emm_contam$contrasts)

# --- Model 2: Comparing Species Preservation (Col A) ---
# We want to ensure the method didn't significantly lower A compared to Raw.
mod_preservation <- lmer(Col_A_SigSpecies ~ Method + (1|SampleID), data = data_sim)

cat("\n=== Pairwise Comparisons (Species Preservation) ===\n")
emm_pres <- emmeans(mod_preservation, specs = pairwise ~ Method)
print(emm_pres$contrasts)

# ==========================================
# 4. AUTOMATED INTERPRETATION FUNCTION
# ==========================================
interpret_results <- function(model_emmeans, label) {
  df <- as.data.frame(model_emmeans$contrasts)
  sig <- df %>% filter(p.value < 0.05)
  if(nrow(sig) > 0) {
    print(paste("For", label, "the following differences are significant:"))
    print(sig[,c("contrast", "p.value")])
  } else {
    print(paste("No significant differences found for", label))
  }
}

interpret_results(emm_contam, "Contaminant Removal (Lower is better)")
interpret_results(emm_pres, "Species Preservation (Higher is usually better)")
