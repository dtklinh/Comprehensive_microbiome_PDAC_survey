## inspect

library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)
library(lme4)       # For alpha diversity mixed models
library(lmerTest)   # To get p-values for lme4
library(vegan)      # For beta diversity (adonis2)

rm(list = ls())
pseq <- readRDS("./data/Chap2/px_J_U_NCT_v02.rds") %>% 
  ps_get() %>% 
  tax_filter(min_total_abundance = 3, min_prevalence = 2) %>% 
  #WrenchWrapper(grp = "HOOD_BENCH")
  rarefy_even_depth()
  
metadata <- data.frame(sample_data(pseq))
## Alpha
# Calculate alpha diversity metrics
alpha_div <- estimate_richness(pseq, measures = c("Observed", "Shannon"))
alpha_div$SampleID <- rownames(alpha_div)

# Merge with metadata
df_alpha <- merge(metadata, alpha_div, by = "row.names")

# Run Linear Mixed Model (LMM)
# Formula: Diversity ~ Fixed_Effect + (1 | Random_Effect)
model_shannon <- lmer(Observed ~ HOOD_BENCH + (1 | person), data = df_alpha)
#model_shannon <- lmer(Observed ~ person + (1|HOOD_BENCH) , data = df_alpha)

# Check results
summary(model_shannon)
anova(model_shannon) # P-value for Condition after accounting for Technician

# 1. Calculate Distance Matrix (e.g., Bray-Curtis)
dist_matrix <- phyloseq::distance(pseq, method = "bray")

# 2. Define the permutation structure (Blocking by Technician)
# This is the "proper" way to handle batch as a nuisance factor
perm_ctrl_person <- how(blocks = metadata$person, nperm = 99999)
perm_ctrl_cnd <- how(blocks = metadata$HOOD_BENCH, nperm = 99999)

## check the permutation to see if it perform as expected
check(metadata, control = perm_ctrl)
metadata[shuffle(nrow(metadata), control = perm_ctrl_cnd),]

# 3. Run PERMANOVA
# We test 'Condition', and permutations are restricted within each 'Technician'
adonis_res <- adonis2(dist_matrix ~ HOOD_BENCH, 
                      data = metadata, 
                      permutations = 9999
                      #permutations = perm_ctrl_person
                      )

print(adonis_res)

ord <- ordinate(pseq, method = "PCoA", distance = dist_bc_batch_corrected)
plot_ordination(pseq, ord, color = "HOOD_BENCH") +
  geom_point(size = 3) +
  theme_classic() +
  stat_ellipse()

### remove batch effect
dist_bc <- phyloseq::distance(pseq, method = "bray")
db <- capscale(dist_bc ~ person, data = metadata)
# Extract residual distance matrix
dist_bc_batch_corrected <- residuals(db, type = "distance")

