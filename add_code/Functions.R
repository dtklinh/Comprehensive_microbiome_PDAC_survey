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
survey_NCT <- function(pseq, lst_NCT, by = "abund", thres = 0.01){
  ## keep only taxa with abundant >=0.01, and inspect its overlap with taxa in NCT
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
  if(!by %in% c("abund", "prev")){
    warning("by must be either abund or prev!")
    return(NA)
  }
  if(by == "abund"){
    taxa_per_sample <- apply(otu_table, 2, function(x) {
      rownames(otu_table)[x/sum(x) > 0.01] 
    })
  }
  
  ##--- cont
  sample_taxa_df <- data.frame(
    SampleID = names(taxa_per_sample),
    Taxa_List = I(taxa_per_sample)  # I() preserves the list inside the column
  ) %>% 
    rowwise() %>% 
    mutate(NumAbund_1Per = length(Taxa_List),
           NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
    ungroup() %>% 
    select(c("SampleID", "NumAbund_1Per", "NumInNCT"))
  df_final <- df_final %>% 
    left_join(., sample_taxa_df, by ="SampleID")
  return(df_final)
}

## write to file of survey
