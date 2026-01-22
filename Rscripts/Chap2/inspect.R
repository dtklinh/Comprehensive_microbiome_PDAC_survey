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

### remove effect of batch for ord plot

#ord <- ordinate(pseq, method = "PCoA", distance = dist_bc_batch_corrected)
ord <- ordinate(pseq, method = "PCoA", distance = dist_matrix)
plot_ordination(pseq, ord, color = "HOOD_BENCH", shape = "person") +
  geom_point(size = 3) +
  #theme_classic() +
  theme_bw() +
  stat_ellipse() +
  ggtitle("PCoA: Condition vs Technician Batch Effect")

### remove batch effect
dist_bc <- phyloseq::distance(pseq, method = "bray")
db <- capscale(dist_bc ~ person, data = metadata)
# Extract residual distance matrix
dist_bc_batch_corrected <- residuals(db, type = "distance")

## Covariate approach
# Technician is added first to "soak up" its variance before testing Condition
adonis_covar <- adonis2(dist_matrix ~  HOOD_BENCH + person, 
                        data = metadata, 
                        by = "terms")

###------------ limma for visualization
# Required libraries
library(phyloseq)
library(limma)
library(microbiome) # For easy CLR transformation
library(patchwork)  # To show plots side-by-side

# 1. Prepare Data (Use CLR transformation)
# Linear models like limma work best on log-like, continuous data
#ps_clr <- microbiome::transform(pseq, "compositional")
rm(list = ls())
trans <- "clr" ## compositional, clr
dist_mx <- "euclidean"   ## euclidean, bray
ps <- readRDS("./data/Chap2/px_J_U_NCT_v02.rds") %>% 
  ps_get() %>% 
  tax_filter(min_total_abundance = 2, min_prevalence = 2) %>% 
  microbiome::transform(transform = trans)
#ps_comp <- microbiome::transform(pseq, transform = "compositional")

# Extract components
otu_table_clr <- as(otu_table(ps_comp), "matrix")
metadata <- data.frame(sample_data(ps_comp))

# 2. Batch Correction
# IMPORTANT: Include 'design' to protect the effect of 'Condition' 
# so you don't accidentally remove your biological signal!
design <- model.matrix(~HOOD_BENCH, data = metadata)
otu_corrected <- removeBatchEffect(otu_table_clr, 
                                   batch = metadata$person, 
                                   design = design)

# 3. Create a New Phyloseq Object with Corrected Data
ps_corrected <- ps_comp
otu_table(ps_corrected) <- otu_table(otu_corrected, taxa_are_rows = TRUE)

# 4. Generate Ordinations
# For CLR data, Euclidean distance is equivalent to Aitchison distance
ord_before <- ordinate(ps_comp, method = "PCoA", distance = "bray")
ord_after  <- ordinate(ps_corrected, method = "PCoA", distance = "bray")

# 5. Plot Comparison
p1 <- plot_ordination(ps_comp, ord_before, color = "HOOD_BENCH" 
                      #shape = "person"
                      ) +
  geom_point(size = 3) +
  theme_bw() +
  stat_ellipse() +
  ggtitle("Before Batch Correction (Technician Effect Visible)")

p2 <- plot_ordination(ps_corrected, ord_after, color = "HOOD_BENCH" 
                      #shape = "person"
                      ) +
  geom_point(size = 3) +
  theme_bw() +
  stat_ellipse() +
  ggtitle("After Batch Correction (Technician Effect Removed)")

# Show side-by-side
p1 + p2 + plot_layout(guides = "collect")


#### residual approach
otu_tab <- microbiome::abundances(ps)
metadata <- microbiome::meta(ps)

# Function to get residuals from Technician effect
get_residuals <- function(taxa_vector, batch_factor) {
  # We fit: Abundance ~ Technician
  fit <- lm(taxa_vector ~ batch_factor)
  # We return the residuals
  return(residuals(fit))
}

otu_residuals <- apply(otu_tab, 1, function(x) get_residuals(x, metadata$person))
otu_residuals <- t(otu_residuals)

ps_resid <- ps
otu_table(ps_resid) <- otu_table(otu_residuals, taxa_are_rows = TRUE)

# Before (Original CLR)
ord_before <- ordinate(ps, method = "PCoA", distance = dist_mx)
p1 <- plot_ordination(ps, ord_before, color = "HOOD_BENCH", shape = "person") +
  theme_bw() + ggtitle("Original (With TA as Batch Effect)")

# After (Residuals)
ord_after <- ordinate(ps_resid, method = "PCoA", distance = dist_mx)
p2 <- plot_ordination(ps_resid, ord_after, color = "HOOD_BENCH", shape = "person") +
  theme_bw() + ggtitle("Residuals (Technician Regressed Out)")

library(patchwork)
p1 / p2
