## a function to calculate different statistical test
## pseq_t: phyloseq of true sample
## pseq_n: phyloseq of negative control sample
## NCT_batch: a list: key is names in column of sample_data, and values are the one we wannna keep.
## prefix: prefix to add to column name of the results, if NULL, it will concatenate values in NCT_batch.
stat_test <- function(pseq_t, pseq_n, NCT_batch = list("NCT_type" = c("pcr", "seq")), prefix = NULL){
  if(!is.null(NCT_batch) && is.null(prefix)){
    prefix = paste0(NCT_batch[[1]], collapse = "_")
  }
  if(length(NCT_batch)!=1){
    warning("NCT_batch must be a list, with one key and values!")
    return(NULL)
  }
  N1 <- nsamples(pseq_t)
  pseq_n <- pseq_n %>% ps_filter(.data[[names(NCT_batch)]] %in% NCT_batch[[1]])
  N2 <- nsamples(pseq_n)
  df_t_prev <- pseq_t %>%
    ps_get() %>%
    prevalence(count = T, detection = 1) %>% 
    tibble(Tax = names(.) ,Prev_true = .)
  df_n_prev <- pseq_n %>% 
    ps_get() %>%
    prevalence(count = T, detection = 1) %>% 
    { tibble(Tax = names(.) , !!paste0(prefix, "Prev_nct") := .) }
  #tibble(Tax = names(.) ,!!paste0(prefix, "Prev_nct", collapse = "_") := .)
  df_t_n_prev <- df_t_prev %>% 
    left_join(., df_n_prev, by = "Tax") %>% 
    mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
    rowwise() %>% 
    mutate(!!paste0(prefix, "Binom") := binom.test(c_across(2), N1, c_across(3)/N2, alternative = "greater")$p.value,
           !!paste0(prefix, "Fisher") := fisher.test(matrix(c(c_across(2), N1 - c_across(2), 
                                                              c_across(3), N2 - c_across(3)), nrow = 2, byrow = TRUE), alternative = "greater")$p.value) %>% 
    ungroup()
  return(df_t_n_prev)
}

## write a function to plot a curve, in order to determine a cutoff value for high prevalence in NCT samples.
HighPrevalence_Data <- function(True.Sample, NCT.Sample){
  num_taxa <- ntaxa(True.Sample)
  # subset of taxa which in in overlap
  idx_overlap <- intersect(taxa_names(True.Sample) %>% unique(),
                           taxa_names(NCT.Sample) %>% unique())
  NCT.Sample.Subset <- prune_taxa(idx_overlap, NCT.Sample)
  prev.nct <- tibble(
    prev = prevalence(NCT.Sample.Subset, detection  = 0, sort = TRUE, count = FALSE),
    OTU =as.factor(names(prevalence(NCT.Sample.Subset , detection = 0, sort = TRUE, count = FALSE)))
  )
  x <- prev.nct$prev
  factorx <- factor(cut(x, breaks=nclass.Sturges(x)))
  xout <- as.data.frame(table(factorx)) %>% map_df(rev)
  xout <- mutate(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
  xout <- xout %>% mutate(rel_tax=cumFreq/num_taxa)
}