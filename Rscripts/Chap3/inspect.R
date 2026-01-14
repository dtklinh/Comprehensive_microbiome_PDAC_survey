## inspect SCRuB

Ts <- c("Peptonophilus_genitalis", "Duncaniella_sp._B8", "Lachnoclostridium_phocaeense")

tmp1 <- pseq_no %>% taxa_sums()
tmp1[names(tmp1) %in% Ts]

tmp2 <- pseq_restrictive %>% taxa_sums()
tmp2[names(tmp2) %in% Ts]

###---------------------

pseq_no_WN <- pseq_no %>% 
  append_AN_NR(df_additionalInfo) %>% 
  WrenchWrapper(grp = "Sex")
otu_table <- abundances(pseq_no_WN)

tmp_prev <- names(prevalence(pseq_no_WN)[prevalence(pseq_no_WN) >=0.75])
taxa_per_sample <- apply(otu_table, 2, function(x){
  intersect(rownames(otu_table)[x > 0], tmp_prev)
})

survey_NCT(pseq_no_WN, lst_NCT, by = "both", thres_prev = 0.15, thres_abd = 0.01)


### 
# --- Bio-statistical Analysis Function ---

analyze_enrichment <- function(n1, N1, n2, N2) {
  
  # 1. Setup the Contingency Table
  # We compare counts of (Positive vs Negative) in (True vs Background)
  # Rows: Group (True, Background)
  # Cols: E. coli Status (Present, Absent)
  
  true_neg <- N1 - n1
  back_neg <- N2 - n2
  
  contingency_table <- matrix(
    c(n1, true_neg,    # True Samples row
      n2, back_neg),   # Background Samples row
    nrow = 2,
    byrow = TRUE,
    dimnames = list("Group" = c("True_Sample", "Background"),
                    "E.coli" = c("Present", "Absent"))
  )
  
  print("--- Contingency Table ---")
  print(contingency_table)
  
  # 2. Run Fisher's Exact Test
  # alternative = "greater" tests if True Sample > Background
  fisher_res <- fisher.test(contingency_table, alternative = "greater")
  
  # 3. Output Results
  p_val <- fisher_res$p.value
  odds_ratio <- fisher_res$estimate
  
  print("--- Statistical Results ---")
  cat(sprintf("Prevalence in True Samples: %.2f%%\n", (n1/N1)*100))
  cat(sprintf("Prevalence in Background:   %.2f%%\n", (n2/N2)*100))
  cat(sprintf("Odds Ratio: %.2f\n", odds_ratio))
  cat(sprintf("P-value: %.5f\n", p_val))
  
  # 4. Bio-statistical Interpretation
  cat("\n--- Conclusion ---\n")
  if(p_val < 0.05) {
    cat("RESULT: SIGNIFICANT.\n")
    cat("E. coli is significantly enriched in your true samples compared to background.\n")
    cat("It is likely a 'real' species in this context (Signal > Noise).")
  } else {
    cat("RESULT: NOT SIGNIFICANT.\n")
    cat("The presence of E. coli in true samples is not statistically different from the background.\n")
    cat("We cannot rule out that it is simply background contamination.")
  }
}

# --- EXAMPLE USAGE ---

# Example data:
# True samples: 100 collected, 15 have E. coli
# Background:   100 collected, 2 have E. coli
analyze_enrichment(n1 = 15, N1 = 100, n2 = 2, N2 = 100)

## test function
pseq_t <- pseq_true_Nj
pseq_n <- pseq_NCT
prefix <- paste0(c("pcr", "seq"), collapse = "_")
NCT_batch = list("NCT_type" = c("pcr", "seq"))
N1 <- nsamples(pseq_t)
pseq_n <- pseq_n %>% ps_filter(.data[[names(NCT_batch)]] %in% NCT_batch[[1]])
N2 <- nsamples(pseq_n)
df_t_prev <- pseq_t %>%
  ps_get() %>%
  prevalence(count = T, detection = 1) %>% 
  tibble(Tax = names(.) ,Prev_true = .)
new_name <- paste0(prefix, "Prev_nct")
df_n_prev <- pseq_n %>% 
  ps_get() %>%
  prevalence(count = T, detection = 1) %>%
  { tibble(Tax = names(.) , !!paste0(prefix, "Prev_nct") := .) }
df_t_n_prev <- df_t_prev %>% 
  left_join(., df_n_prev, by = "Tax") %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  rowwise() %>% 
  mutate(!!paste0(prefix, "Binom") := binom.test(c_across(2), N1, c_across(3)/N2, alternative = "greater")$p.value,
         !!paste0(prefix, "Fisher") := fisher.test(matrix(c(c_across(2), N1 - c_across(2), 
                                                            c_across(3), N2 - c_across(3)), nrow = 2, byrow = TRUE), alternative = "greater")$p.value) %>% 
  ungroup()
###-----------------------
### JSD testing
pseq_true <- readRDS("data/Chap3/pseq_bulk_Fisher_v02.rds") %>% 
  ps_
pseq_true2 <- readRDS("data/Chap3/pseq_FFPE_Biom.rds")
pseq <- merge_phyloseq(pseq_true, pseq_true2)
dis <- distance(pseq, method = "jsd")
#### wrench with only one group
df_additionalInfo <- readxl::read_xlsx("./meta/Mice_meta.xlsx")
pseq_raw <- append_AN_NR(pseq_raw, df_additionalInfo)
cnt_table <- pseq_raw %>% otu_table() %>% as.matrix()
group <- pseq_raw %>% sample_data() %>% pull("Sex")
#group <- c(rep("A", 18))
w <- wrench(cnt_table, condition = group)
norm_counts <- Wrench::normalize(w)
# Update your phyloseq object
ps_wrench <- PhyloObjct
otu_table(ps_wrench) <- otu_table(norm_counts, taxa_are_rows = TRUE)

### rarefaction curve with different colors
group <- meta(pseq_raw_both)$true.control
cols <- c("true" = "blue", "NCT" = "red")
xxx <- rarecurve(t(abundances(pseq_raw_both)), step = 100, col = cols[group], label = FALSE, xlim = c(0, 80000), tidy = TRUE, se = FALSE)
xxx <- xxx %>%
  mutate(Group = group)
xx2 <- legend(
  "bottomright",
  legend = names(cols),
  col = cols,
  lty = 1,
  bty = "n"
)

###### t paired-test for phyloseq whose reads count is not whole number
library(lme4)       # Linear Mixed Models (for repeated measures)
library(lmerTest)   # P-values for Mixed Models
library(emmeans)    # Pairwise comparisons
rm(list = ls())
source("./add_code/Functions.R")
df_additionalInfo <- readxl::read_xlsx("./meta/Mice_meta.xlsx")
pseq_raw <- readRDS("./data/Chap3/pseq_Proj5_postFilter_v04.rds") %>% 
  ps_filter(ffpe.bulk == "bulk") %>% 
  ps_filter(true.control == "true") %>% 
  append_AN_NR(df_additional = df_additionalInfo) %>% 
  WrenchWrapper(grp = "Sex") 
pseq_restrictive <- readRDS("./data/Chap3/pseq_bulk_restrictive.rds") %>% 
  append_AN_NR(df_additionalInfo) %>% 
  WrenchWrapper(grp = "Sex")
pseq_decontam <- readRDS("./data/Chap3/pseq_bulk_decontam_p0.5.rds") %>% 
  append_AN_NR(df_additionalInfo) %>% 
  WrenchWrapper(grp = "Sex")
pseq_SCRuB <- readRDS("./data/Chap3/pseq_bulk_SCRuB.rds") %>% 
  append_AN_NR(df_additionalInfo) %>% 
  WrenchWrapper(grp = "Sex")
pseq_Nj <- readRDS("./data/Chap3/pseq_bulk_Fisher_v02.rds") %>% 
  append_AN_NR(df_additionalInfo) %>% 
  WrenchWrapper(grp = "Sex")

## meta
df_meta <- meta(pseq_raw) %>% 
  rownames_to_column(var = "SampleID") %>% 
  dplyr::select(SampleID, AN_NR)
## observed species
df_raw <- pseq_raw %>% 
  otu_table() %>% 
  specnumber(MARGIN = 2) %>% 
  tibble(SampleID = names(.), Raw = .)
df_restric <- pseq_restrictive %>% 
  otu_table() %>% 
  specnumber(MARGIN = 2) %>% 
  tibble(SampleID = names(.), Restrictive = .)
df_Decontam <- pseq_decontam %>% 
  otu_table() %>% 
  specnumber(MARGIN = 2) %>% 
  tibble(SampleID = names(.), decontam = .)
df_SCRuB <- pseq_SCRuB %>% 
  otu_table() %>% 
  specnumber(MARGIN = 2) %>% 
  tibble(SampleID = names(.), SCRuB = .)
df_Nj <- pseq_Nj %>% 
  otu_table() %>% 
  specnumber(MARGIN = 2) %>% 
  tibble(SampleID = names(.), Nj = .)
df_merge <- df_meta %>% 
  merge(., df_raw, by = "SampleID") %>% 
  merge(., df_restric, by = "SampleID") %>% 
  merge(., df_Decontam, by = "SampleID") %>% 
  merge(., df_SCRuB, by = "SampleID") %>% 
  merge(., df_Nj, by = "SampleID") %>% 
  select(-SampleID)
df_merge_loner <- df_merge %>% 
  pivot_longer(cols = -AN_NR, names_to = "Method", values_to = "ObsSpec")

model <- lmer(ObsSpec ~ Method + (1|AN_NR), data = df_merge_loner)
# Check assumptions (Residuals should be roughly normal)
qqnorm(resid(model))
qqline(resid(model))

print(anova(model))
emm <- emmeans(model, specs = pairwise ~ Method)
contr <- contrast(emm, method = "pairwise")

tmp <- contr %>% 
  as_tibble() %>%
  transmute(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".* - ", "", contrast),
    estimate,
    se = SE,
    statistic = t.ratio,
    df,
    p = p.value
    ##method = "emmeans model-based t-test"
  ) 
tmp <- tmp %>% 
  add_significance() %>% 
  mutate(y.position = c(445, 481.6667, 518.3333, 555, 591.6667, 628.3333, 665, 701.6667, 738.3333, 775))

## barplot
plt <- ggplot(df_merge_loner, aes(x = Method, y = ObsSpec, fill = Method)) +
  geom_violin() +
  #geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_boxplot(width=0.2) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  # Connect lines for the same sample to show paired changes
  geom_line(aes(group = AN_NR), color = "gray", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Comparison of Decontamination Methods",
    subtitle = "Metric: Penalized Clean Yield (Higher is Better)",
    y = "Composite Score\n(Yield x Purity)"
  ) +
  theme(legend.position = "none")
plt + stat_pvalue_manual(tmp, label = "p.signif", inherit.aes = FALSE, tip.length = 0.01)

ob <- t_test(formula = ObsSpec ~ Method, data = df_merge_loner , paired = T) %>% 
  #adjust_pvalue(method = "BH") %>%  
  #add_significance("p.adj") %>% 
  add_xy_position()
## wrap to a function
AlphaPlot_Violin_LMM <- function(df, SampleID, strata, val){
  ## df: in pivot longer manner
  ## calc p value in manner of paired
  ## get y.position for ggplot
  y_pos <- t_test(formula = as.formula(sprintf("%s ~ %s", val, strata)), data = df) %>% add_xy_position() %>% pull(y.position)
  model <- lmer(as.formula(sprintf("%s ~ %s + (1|%s)", val, strata, SampleID)), data = df)
  emm <- emmeans(model, specs = as.formula(sprintf("pairwise ~ %s", strata)))
  contr <- contrast(emm, method = "pairwise")
  
  p_add <- contr %>% 
    as_tibble() %>%
    transmute(
      group1 = sub(" - .*", "", contrast),
      group2 = sub(".* - ", "", contrast),
      estimate,
      se = SE,
      statistic = t.ratio,
      df,
      p = p.value
      ##method = "emmeans model-based t-test"
    ) %>% 
    add_significance() %>% 
    mutate(y.position = y_pos)
  plt <- ggplot(df, aes(x = .data[[strata]], y = .data[[val]], fill = .data[[strata]])) +
    geom_violin() +
    #geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_boxplot(width=0.2) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
    # Connect lines for the same sample to show paired changes
    geom_line(aes(group = AN_NR), color = "gray", alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Comparison of Decontamination Methods",
      subtitle = "Metric: Penalized Clean Yield (Higher is Better)",
      y = "Composite Score\n(Yield x Purity)"
    ) +
    theme(legend.position = "none")
  plt <- plt + stat_pvalue_manual(p_add, label = "p.signif", inherit.aes = FALSE, tip.length = 0.01)
  return(plt)
}
AlphaPlot_Violin_LMM(df_merge_loner, SampleID = "AN_NR", strata = "Method", val = "ObsSpec")
extract_df <- function(lst_pseq, diversity = "obs"){
  ## diversity: obs, shannon, simpson, invsimpson
  ## extract meta
  df_meta <- meta(lst_pseq[[1]]) %>% 
    rownames_to_column(var = "SampleID") %>% 
    dplyr::select(SampleID, AN_NR)
  df_full <- df_meta
  for(nm in names(lst_pseq)){
    print(nm)
    if(diversity == "obs"){
      div <- lst_pseq[[nm]] %>% 
        otu_table() %>% 
        specnumber(MARGIN = 2)
      df_tmp <- tibble(SampleID = names(div), !!nm := unname(div))
    }else{
      div <- lst_pseq[[nm]] %>% 
        otu_table() %>% 
        vegan::diversity(index = diversity, MARGIN = 2)
      df_tmp <- tibble(SampleID = names(div), !!nm := unname(div))
    }
    
    df_full <- df_full %>% 
      merge(., df_tmp, by = "SampleID")
    print(df_full)
  }
  df_full_longer <- df_full %>% 
    select(-SampleID) %>% 
    pivot_longer(cols = -AN_NR, names_to = "Method", values_to = "Val")
  return(df_full_longer)
}
ls <- list(raw = pseq_raw, decontam = pseq_decontam, restrictive = pseq_restrictive, SCRuB = pseq_SCRuB, Nj = pseq_Nj)
extract_df(ls)

###### Beta diversity with PERMANOVA and paired (constraint permutation)
pseq <- pseq_merge
taxa_rank <- "genus"
m_group <- "Decon_type"
unconstrained_aitchison_pca <- pseq %>% 
  tax_agg(taxa_rank) %>% 
  tax_transform("clr") %>% 
  ord_calc()
pca_plot <- unconstrained_aitchison_pca %>% 
  ord_plot(
    plot_taxa = 1:4, colour = m_group, size = 1.25, tax_vec_length = 0.325, 
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 1), auto_caption = 8
  )
customised_plot <- pca_plot +
  stat_ellipse(aes(linetype = .data[[m_group]], colour = .data[[m_group]]), linewidth = 0.3) + # linewidth not size, since ggplot 3.4.0
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 1, clip = "off") # makes rotated labels align correctly
## calculate p value
aitchison_dists <- pseq %>%
  tax_transform("identity", rank = taxa_rank) %>%
  dist_calc("aitchison")
aitchison_perm <- aitchison_dists %>%
  dist_permanova(
    seed = 210488, # for set.seed to ensure reproducibility of random process
    n_processes = 2, n_perms = 9999, # you should use at least 999!
    variables = m_group
    #strata = meta(aitchison_dists)[["AN_NR"]]
  )

p_val <- perm_get(aitchison_perm) %>% as.data.frame() %>% pull(`Pr(>F)`) %>% .[1]
p <- case_when(
  p_val > 0.05 ~ paste("p =", round(p_val,4), "n.s.", sep = " "),
  p_val < 0.05 &  p_val > 0.01 ~ paste("p =", round(p_val,4), "*", sep = " "),
  p_val <= 0.01 & p_val > 0.001  ~ paste("p =", round(p_val,4), "**", sep = " "),
  p_val <= 0.001 ~ paste("p =",round(p_val,4), "***", sep = " "),
)
# annotation
annotations <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = p,
  hjustvar = c(-0.2) ,
  vjustvar = c(1.5))
res_plot <- customised_plot +
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,
                                 label=annotateText), size = 4.5, inherit.aes = FALSE)

ad <- aitchison_dists@dist
res <- pairwise.adonis2(x = ad ~ Decon_type, data = meta(aitchison_dists), strata = "AN_NR")

### another Venn diagram
# Install packages if you haven't already
# install.packages("VennDiagram")
# install.packages("gridExtra")

library(VennDiagram)
library(grid)

# 1. Define the data based on your image
raw_species <- c(
  "Sphingomonas sp. AAP5", "Sphingomonas alpina", "Escherichia coli", 
  "Sphingomonas sp. QA11", "Streptococcus suis", "Alistipes shahii", 
  "Sphingomonas sp. So64.6b", "Alistipes finegoldii", "Mammaliicoccus sciuri"
)

both_species <- c(
  "Mammaliicoccus lentus", "Streptococcus respiraculi", "Alistipes senegalensis", 
  "Staphylococcus xylosus", "Ligilactobacillus murinus", "Lactobacillus johnsonii"
)

restrictive_species <- c(
  "Faecalibaculum rodentium", "Dialister succinatiphilus", "Streptococcus constellatus", 
  "Dialister massiliensis", "Megasphaera stantonii", "Kineothrix sp. MB12-C1", 
  "Limosilactobacillus reuteri", "Lactobacillus taiwanensis", "Subdoligranulum variabile"
)

# 2. Create the Venn Diagram Object
# Note: category.names are "Raw" and "restrictive"
venn_plot <- draw.pairwise.venn(
  area1 = length(raw_species) + length(both_species),
  area2 = length(restrictive_species) + length(both_species),
  cross.area = length(both_species),
  category = c("Raw", "restrictive"),
  fill = c("lightblue", "pink"),
  alpha = c(0.5, 0.5),
  lty = "solid",
  cex = 2,               # Size of the numbers (9, 6, 9)
  cat.cex = 1.5,         # Size of category titles
  cat.pos = c(-20, 20),  # Position of titles
  ext.text = FALSE
)

# 3. Create a layout with a plot on the left and legend on the right
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.5, 0.5), "npc"))))

# Draw the Venn on the left
pushViewport(viewport(layout.pos.col = 1))
grid.draw(venn_plot)
popViewport()

# Draw the Legend on the right
pushViewport(viewport(layout.pos.col = 2, x = 0.1, just = "left"))

# Function to draw species list with a color header
draw_list <- function(title, species, color, y_start) {
  grid.text(title, x = 0, y = y_start, just = "left", 
            gp = gpar(col = color, fontface = "bold", fontsize = 12))
  for (i in seq_along(species)) {
    grid.text(species[i], x = 0.05, y = y_start - (i * 0.03), 
              just = "left", gp = gpar(col = "black", fontsize = 10, fontitalic = TRUE))
  }
}

# Add the three sections of the legend
draw_list("Raw Only", raw_species, "blue", 0.9)
draw_list("Common (Intersection)", both_species, "purple", 0.55)
draw_list("Restrictive Only", restrictive_species, "red", 0.3)

popViewport()


## for 3 sets
# install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# 1. Define Example Data (7 groups)
set_a_only <- c("Species A1", "Species A2")
set_b_only <- c("Species B1", "Species B2")
set_c_only <- c("Species C1", "Species C2")

int_ab <- c("Species AB1")       # Only in A and B
int_bc <- c("Species BC1")       # Only in B and C
int_ac <- c("Species AC1")       # Only in A and C
int_abc <- c("Species ABC1")     # In all three

# 2. Calculate totals for the Venn function
# VennDiagram needs the total count per set and total count per intersection
n1 <- length(set_a_only) + length(int_ab) + length(int_ac) + length(int_abc)
n2 <- length(set_b_only) + length(int_ab) + length(int_bc) + length(int_abc)
n3 <- length(set_c_only) + length(int_ac) + length(int_bc) + length(int_abc)
n12 <- length(int_ab) + length(int_abc)
n23 <- length(int_bc) + length(int_abc)
n13 <- length(int_ac) + length(int_abc)
n123 <- length(int_abc)

# 3. Create the Venn Diagram Object
venn_plot <- draw.triple.venn(
  area1 = n1, area2 = n2, area3 = n3,
  n12 = n12, n23 = n23, n13 = n13, n123 = n123,
  category = c("Group A", "Group B", "Group C"),
  fill = c("skyblue", "pink", "palegreen"),
  alpha = 0.5,
  lty = "blank",
  cex = 2,            # Size of numbers
  cat.cex = 1.5,      # Size of Group Titles
  cat.col = c("blue", "red", "darkgreen")
)

# 4. Setup Layout (Venn on Left, Legend on Right)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.4, 0.6), "npc"))))

# Draw Venn
pushViewport(viewport(layout.pos.col = 1))
grid.draw(venn_plot)
popViewport()

# 5. Create Legend on the Right
pushViewport(viewport(layout.pos.col = 2, x = 0.1, just = "left"))

draw_sub_legend <- function(title, species, color, y_pos) {
  grid.text(title, x = 0, y = y_pos, just = "left", 
            gp = gpar(col = color, fontface = "bold", fontsize = 10))
  if(length(species) > 0) {
    species_text <- paste(species, collapse = ", ")
    # Wrap text if it's too long
    wrapped_text <- paste(strwrap(species_text, width = 50), collapse = "\n")
    grid.text(wrapped_text, x = 0.02, y = y_pos - 0.05, just = c("left", "top"), 
              gp = gpar(col = "black", fontsize = 9, fontitalic = TRUE))
  }
}

# Positioning the 7 categories
draw_sub_legend("Unique to A", set_a_only, "blue", 0.95)
draw_sub_legend("Unique to B", set_b_only, "red", 0.82)
draw_sub_legend("Unique to C", set_c_only, "darkgreen", 0.69)
draw_sub_legend("A & B Intersection", int_ab, "purple", 0.56)
draw_sub_legend("B & C Intersection", int_bc, "brown", 0.43)
draw_sub_legend("A & C Intersection", int_ac, "darkcyan", 0.30)
draw_sub_legend("Common to All (A,B,C)", int_abc, "black", 0.17)

popViewport()
