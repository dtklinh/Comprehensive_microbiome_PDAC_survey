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
