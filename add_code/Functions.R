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
survey_NCT <- function(pseq, lst_NCT){
  df_final <- pseq %>% 
    microbiome::alpha(index = c("observed", "diversity_shannon")) %>% 
    rownames_to_column(var = "SampleID")
  df_read <- pseq %>% 
    sample_sums() %>% 
    tibble(SampleID = names(.), reads = .)
  df_final <- df_final %>% 
    left_join(., df_read, by = "SampleID")
  otu_table <- abundances(pseq) 
  taxa_per_sample <- apply(otu_table, 2, function(x) {
    rownames(otu_table)[x/sum(x) > 0.01] 
  })
}
