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
