#############################################
###------  Linear Mixed Model with Fixed and random effect ------------
## Fixed: Decontamination method
## Radom: SampleID
# Load necessary libraries
library(tidyverse)  # Data manipulation and plotting
library(lme4)       # Linear Mixed Models (for repeated measures)
library(lmerTest)   # P-values for Mixed Models
library(emmeans)    # Pairwise comparisons

# ==========================================
# 1. DATA PREPARATION
# ==========================================
# --- Load data
df_origin <- read.table("results/Chap3/survey_overlap_NCT/df_origin_concat.tsv", header = T, sep = "\t", check.names = F) %>% 
  mutate(Method = "origin")
df_decontam <- read.table("results/Chap3/survey_overlap_NCT/df_decontam_concat.tsv", header = T, sep = "\t", check.names = F) %>% 
  mutate(Method = "decontam")
df_Fisher <- read.table("results/Chap3/survey_overlap_NCT/df_Fisher_concat.tsv", header = T, sep = "\t", check.names = F) %>% 
  mutate(Method = "Fisher")
df_restrictive <- read.table("results/Chap3/survey_overlap_NCT/df_restrictive_concat.tsv", header = T, sep = "\t", check.names = F) %>% 
  mutate(Method = "restrictive")
df_SCRuB <- read.table("results/Chap3/survey_overlap_NCT/df_SCRuB_concat.tsv", header = T, sep = "\t", check.names = F) %>% 
  mutate(Method = "SCRuB")

df_all <- rbind(df_origin, df_decontam, df_Fisher, df_restrictive, df_SCRuB)
## save table for later usage
write.table(df_all, "results/Chap3/survey_overlap_NCT/df_all_concat.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

# --- Generate Dummy Data (Replace this block with your actual data loading) ---
# set.seed(42)
# samples <- 18
# method_list <- c("Method_X", "Method_Y", "Method_Z")
# 
# data_sim <- expand.grid(SampleID = factor(1:samples), Method = method_list) %>%
#   mutate(
#     # Column A: Original Total Species (Baseline richness varies by sample)
#     Col_A_Original = rep(rpois(samples, lambda = 500), times = length(method_list)),
#     
#     # Column B: Species Kept (Varies by method aggressiveness)
#     Col_B_Kept = case_when(
#       Method == "Method_X" ~ floor(Col_A_Original * 0.95), # Keeps almost everything
#       Method == "Method_Y" ~ floor(Col_A_Original * 0.60), # Aggressive cleaning
#       Method == "Method_Z" ~ floor(Col_A_Original * 0.85)  # Balanced
#     ),
#     
#     # Column C: Overlap with Contaminants (The "Noise")
#     Col_C_Contam = case_when(
#       Method == "Method_X" ~ floor(Col_B_Kept * 0.15), # High contamination left
#       Method == "Method_Y" ~ floor(Col_B_Kept * 0.01), # Very clean, but lost data in B
#       Method == "Method_Z" ~ floor(Col_B_Kept * 0.02)  # Low contamination
#     )
#   )
# ---------------------------------------------------------------------------

# ==========================================
# 2. CALCULATE THE SINGLE METRIC
# ==========================================

# We handle potential division by zero if a method removes ALL species (B=0)
df_analysis <- df_all %>%
  select(SampleID, Method, observed, B = NumPrev_Abd, C = NumInNCT_both) %>% 
  mutate(
    # 1. Calculate "Clean Species" count
    Clean_Count = B - C,
    
    # 2. Calculate Yield: (Clean Species / Original Total)
    Yield = Clean_Count / observed,
    
    # 3. Calculate Purity: (Clean Species / Total Kept)
    # If Col_B_Kept is 0, Purity is 0
    Purity = ifelse(B > 0, Clean_Count / B, 0),
    
    # 4. COMPOSITE METRIC: Penalized Clean Yield
    # This rewards keeping good data AND having a clean final table
    Composite_Score = Yield * Purity
  )

# Inspect the calculated metric
head(df_analysis)

# ==========================================
# 3. VISUALIZATION
# ==========================================

plt <- ggplot(df_analysis, aes(x = Method, y = Composite_Score, fill = Method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  # Connect lines for the same sample to show paired changes
  geom_line(aes(group = SampleID), color = "gray", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Comparison of Decontamination Methods",
    subtitle = "Metric: Penalized Clean Yield (Higher is Better)",
    y = "Composite Score\n(Yield x Purity)"
  ) +
  theme(legend.position = "none")
ggsave("results/Chap3/img/Both_compare.png", plt)
# ==========================================
# 4. STATISTICAL MODELING (Linear Mixed Model)
# ==========================================

# Why LMM? 
# 1. Fixed Effect: 'Method' (What we want to compare)
# 2. Random Effect: 'SampleID' (Accounting for biological variation between samples)
model <- lmer(Composite_Score ~ Method + (1|SampleID), data = df_analysis)

# Check assumptions (Residuals should be roughly normal)
qqnorm(resid(model))
qqline(resid(model))

# Print ANOVA table
cat("\n=== ANOVA Table ===\n")
print(anova(model))

# ==========================================
# 5. POST-HOC COMPARISONS
# ==========================================

# Perform pairwise comparisons with Tukey adjustment for multiple testing
emm <- emmeans(model, specs = pairwise ~ Method)

cat("\n=== Pairwise Differences ===\n")
print(emm$contrasts)

# Interpretation Helper
emm_df <- as.data.frame(emm$emmeans)
best_method <- emm_df %>% arrange(desc(emmean)) %>% slice(1) %>% pull(Method)
cat(paste("\nBased on the model, the method performing best (highest mean score) is:", best_method, "\n"))
