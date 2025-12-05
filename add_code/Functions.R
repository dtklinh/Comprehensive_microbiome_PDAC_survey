### Addiontally code for chap 3

## filter taxa according to their prevalence and get the top
get_most_prev <- function(pseq, prev_threshold = 0.5, N = 25){
  res <- pseq %>% 
    prevalence() %>% 
    tibble(taxa = names(.), value = .) %>% 
    filter(value >= prev_threshold) %>% 
    arrange(desc(value)) %>% 
    head(N)
  return(res)
}

## append meta to have more info about mouse
append_AN_NR <- function(pseq, df_additional){
  df_meta <- meta(pseq)
  if(!"AN_NR" %in%  colnames(df_meta)){
    print("There must be AN_NR column!")
    return(NA)
  }
  ## check if info is already in the main data frame, warning if already exits
  col_names_fea <- setdiff(colnames(df_additional), "AN_NR")
  new_col_names_fea <- setdiff(col_names_fea, colnames(df_meta))
  if(length(new_col_names_fea)==0){
    warning("All new feature already exits, nothing to be added")
    return(pseq)
  }
  if(length(new_col_names_fea) < length(col_names_fea)){
    warning("partially info already there, add the rest!")
    df_additional <- df_additional %>% 
      select(c("AN_NR", new_col_names_fea))
  }
  
  df_meta_new <- df_meta %>% 
    rownames_to_column(var = "SampleName") %>% 
    left_join(., df_additional, by = "AN_NR") %>% 
    column_to_rownames(var = "SampleName")
  sample_data(pseq) <- df_meta_new
  return(pseq)
}

## Wrench
WrenchWrapper <- function(PhyloObjct, grp, roundUp = F){
  cnt_table <- PhyloObjct %>% otu_table()
  group <- PhyloObjct %>% sample_data() %>% pull(grp)
  w <- wrench(cnt_table, condition = group)
  
  # deseq.obj <- DESeqDataSetFromMatrix(cnt_table %>% as.data.frame(), DataFrame(group), ~group)
  # DESeq2::sizeFactors(deseq.obj) <- w$nf
  # cnt_table_normalized <- DESeq2::counts(deseq.obj, normalized=TRUE)
  
  norm_factors <- w$nf
  norm_counts <- sweep(cnt_table, 2, norm_factors, FUN = '/')
  if(roundUp){norm_counts <- norm_counts %>% round()}
  if(!is.null(phy_tree(PhyloObjct, errorIfNULL = F))){
    return(phyloseq(otu_table(norm_counts, taxa_are_rows = T), tax_table(PhyloObjct %>% tax_table()), sample_data(PhyloObjct %>% sample_data()),
                    phy_tree(PhyloObjct)))
  } else{
    return(phyloseq(otu_table(norm_counts, taxa_are_rows = T), tax_table(PhyloObjct %>% tax_table()), sample_data(PhyloObjct %>% sample_data())))
  }
}
## Chap 3, extract info
## pseq with normalization and additional info added
survey_NCT <- function(pseq, lst_NCT, by = "abund", thres_abd = 0.01, thres_prev = 0.5){
  ## by abd: keep only taxa with abundant >=0.01, and inspect its overlap with taxa in NCT
  ## by prev: each sample, keep track of species whose prev is bigger than a threshold (in 18 samples). Inspect
  ## them with the list if NCT
  df_final <- pseq %>% 
    microbiome::alpha(index = c("observed", "diversity_shannon")) %>% 
    rownames_to_column(var = "SampleID")
  df_read <- pseq %>% 
    sample_sums() %>% 
    tibble(SampleID = names(.), reads = .)
  df_final <- df_final %>% 
    left_join(., df_read, by = "SampleID")
  otu_table <- abundances(pseq) 
  taxa_per_sample <- NULL
  if(!by %in% c("abund", "prev", "both")){
    warning("by must be either abund or prev, or both!")
    return(NA)
  }
  if(by == "abund"){
    ## Abundance: 
    taxa_per_sample <- apply(otu_table, 2, function(x) {
      rownames(otu_table)[x/sum(x) > thres_abd] 
    })
    sample_taxa_df <- data.frame(
      SampleID = names(taxa_per_sample),
      Taxa_List = I(taxa_per_sample)  # I() preserves the list inside the column
    ) %>% 
      rowwise() %>% 
      mutate(NumAbund_1Per = length(Taxa_List),
             NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
      ungroup() %>% 
      select(c("SampleID", "NumAbund_1Per", "NumInNCT"))
  } else if(by == "prev"){
    ## Prevalent
    tmp_prev <- names(prevalence(pseq)[prevalence(pseq) >=thres_prev])
    taxa_per_sample <- apply(otu_table, 2, function(x){
      intersect(rownames(otu_table)[x > 0], tmp_prev)
    })
    sample_taxa_df <- data.frame(
      SampleID = names(taxa_per_sample),
      Taxa_List = I(taxa_per_sample)  # I() preserves the list inside the column
    ) %>% 
      rowwise() %>% 
      mutate(NumPrev_50Per = length(Taxa_List),
             NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
      ungroup() %>% 
      select(c("SampleID", "NumPrev_50Per", "NumInNCT"))
  }
  else if(by == "both"){
    tmp_prev <- names(prevalence(pseq)[prevalence(pseq) >=thres_prev])
    taxa_per_sample <- apply(otu_table, 2, function(x){
      intersect(rownames(otu_table)[x/sum(x) > thres_abd], tmp_prev)
    })
    sample_taxa_df <- data.frame(
      SampleID = names(taxa_per_sample),
      Taxa_List = I(taxa_per_sample)  # I() preserves the list inside the column
    ) %>% 
      rowwise() %>% 
      mutate(NumPrev_Abd = length(Taxa_List),
             NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
      ungroup() %>% 
      select(c("SampleID", "NumPrev_Abd", "NumInNCT"))
  }
  
  ##--- cont
  # sample_taxa_df <- data.frame(
  #   SampleID = names(taxa_per_sample),
  #   Taxa_List = I(taxa_per_sample)  # I() preserves the list inside the column
  # ) %>% 
  #   rowwise() %>% 
  #   mutate(NumAbund_1Per = length(Taxa_List),
  #          NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
  #   ungroup() %>% 
  #   select(c("SampleID", "NumAbund_1Per", "NumInNCT"))
  df_final <- df_final %>% 
    left_join(., sample_taxa_df, by ="SampleID")
  return(df_final)
}

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
###-----------------------
stat_test_all <- function(pseq_t, pseq_n, ctrl_types = "NCT_type"){
  lst_ctrl <- meta(pseq_n)[[ctrl_types]] %>% unique()
  lst_ctrl <- lst_ctrl[!is.na(lst_ctrl)]
  if(length(lst_ctrl)==0){
    warning("Sth wrong with control types in NCT samples!")
    return(NULL)
  }
  N1 <- nsamples(pseq_t)
  df_res <- pseq_t %>%
    ps_get() %>%
    prevalence(count = T, detection = 1) %>% 
    {tibble(Tax = names(.) , !!sprintf("Prev_true(N=%d)", N1) := .)}
  #{tibble(Tax = names(.) , !!paste0("Prev_true(N=",nsamples(pseq_t), ")") := .)}
  
  #ctrl <- lst_ctrl[1]
  
  
  for(ctrl in lst_ctrl){
    ps <- pseq_n %>% 
      ps_get() %>% 
      ps_filter(!!sym(ctrl_types) == ctrl)
    N2 <- nsamples(ps)
    df_tmp <- ps %>% 
      prevalence(count = T, detection = 1) %>% 
      {tibble(Tax = names(.) , !!sprintf("Prev_%s(N=%d)", ctrl, N2) := .)}
    df_res <- df_res %>% 
      left_join(., df_tmp, by = "Tax") %>% 
      mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
      rowwise() %>% 
      mutate(!!sprintf("%s_Binom", ctrl) := binom.test(x = c_across(contains("Prev_true")), n = N1, p = c_across(contains(sprintf("Prev_%s", ctrl)))/N2, alternative = "greater")$p.value,
             !!sprintf("%s_Fisher", ctrl) :=fisher.test(matrix(c(c_across(contains("Prev_true")), N1 - c_across(contains("Prev_true")), c_across(contains(sprintf("Prev_%s", ctrl))), N2- c_across(contains(sprintf("Prev_%s", ctrl)))), nrow = 2, byrow = TRUE), alternative = "greater")$p.value)
  }
  return(df_res)
}
