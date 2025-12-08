### Addiontally code for chap 3
library(rstatix)
library(ggpubr)

Kolors <- c("#9d547c","#56ca63","#a357d6","cornflowerblue","#419d2a","sandybrown","red3","peachpuff","cyan","paleturquoise3","mistyrose","mediumpurple","mediumseagreen","mediumorchid","moccasin","orange4","olivedrab","midnightblue","papayawhip","palevioletred4","brown1","greenyellow","orchid","navy","darkred","navajowhite1","mistyrose1","grey85","#525fd6","red2","#8cbe3a","#c944aa","indianred3","#5ba557","#9e66cb","#c1b735","#6d82ec","grey25","#e69728","#6654b0","lightsalmon3","lightcyan1","khaki1","seagreen1","plum1","lightsteelblue1","palevioletred3","mintcream","magenta3","#799330","#da7fdf","#3c782c","#e44586","blue4","#63c996","#dc3f53","#49cbc8","#cf3f29","#4fabda","#da6c2b","#598bd1","#b78c24","#8d4191","#a0b971","slategray1","sienna","plum1","lightyellow1","lightskyblue3","linen","limegreen","cornsilk1","mediumaquamarine","gray14","gold3","darkviolet","#b2386a","#479d71","#ae4341","#2ba198","#e07557","#5361a3","#dda353","#aa98df","#5b6114","#dc89bf","#327243","slateblue1","#e57b94","#277257","#9b62a0","#bbab59","#98495a","#526229","#d8827d","#857624","gray40","#9a4a22","#7c7d46","mediumslateblue","lemonchiffon1","#e3a073","#9e6b33", "gray74","slateblue1","rosybrown3", "lawngreen","gainsboro","dodgerblue3","deeppink3","firebrick3", "orchid2", "olivedrab1", "ivory3", "darkseagreen", "bisque2", "darkgoldenrod2", "blue2", "skyblue", "seashell2", "turquoise", "tan1", "seagreen2", "palevioletred3", "linen", "steelblue4","ghostwhite","dodgerblue1","deeppink1","firebrick1", "limegreen", "purple3", "khaki3", "snow3", "darkslategray","darkorchid","lavender", "magenta2", "palegreen", "salmon", "maroon", "cyan2","#671408","#FAEBD7","#7FFFD4","#F0FFFF","#A52A2A","burlywood","cadetblue","#7FFF00","chocolate","cornsilk","slateblue1","#FF7F50","red1","#008B8B","darkgoldenrod1","darkolivegreen","darkorange4","white","hotpink","honeydew1","goldenrod2","darkgreen","oldlace","darkslategray3","navajowhite3","orchid4","gray25","#F0924D")

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
             !!sprintf("%s_Fisher", ctrl) :=fisher.test(matrix(c(c_across(contains("Prev_true")), N1 - c_across(contains("Prev_true")), c_across(contains(sprintf("Prev_%s", ctrl))), N2- c_across(contains(sprintf("Prev_%s", ctrl)))), nrow = 2, byrow = TRUE), alternative = "greater")$p.value) %>% 
      ungroup()
  }
  return(df_res)
}
##- 
stat_test_specific <- function(pseq_t, pseq_n, prefix = "rnd"){
  N1 <- nsamples(pseq_t)
  N2 <- nsamples(pseq_n)
  if(N1==0 || N2==0){ warning("Zero entry phyloseq object!"); return(NULL)}
  df1 <- pseq_t %>%
    ps_get() %>%
    prevalence(count = T, detection = 1) %>% 
    {tibble(Tax = names(.) , !!sprintf("Prev_true(N=%d)", N1) := .)}
  df2 <- pseq_n %>% 
    get_ps() %>% 
    prevalence(count = T, detection = 1) %>% 
    {tibble(Tax = names(.) , !!sprintf("%s_Prev_NCT(N=%d)", prefix, N2) := .)}
  df_res <- df1 %>% 
    left_join(., df2, by = "Tax") %>% 
    mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
    rowwise() %>% 
    mutate(!!sprintf("%s_Binom", prefix) := binom.test(x = c_across(2), n = N1, p = c_across(3)/N2, alternative = "greater")$p.value,
           !!sprintf("%s_Fisher", prefix) := fisher.test(matrix(c(c_across(2), N1 -c_across(2), c_across(3), N2 - c_across(3)) ,nrow = 2, byrow = T), alternative = "greater")$p.value) %>% 
    ungroup()
  return(df_res)
}
####------ microbiome analysis
AlphaPlot_Violin <- function(PhyloObj, index = "Observed", strata = "treatment", y_label = "Observed species", add_legend = FALSE, gg_title = NULL){
  #group.colors <- c(ctrl = "dodgerblue3", Gem = "firebrick2")
  condition_names <- sample_data(PhyloObj)[[strata]] %>% table() %>% names()
  group.colors <- Kolors[1:length(condition_names)]
  names(group.colors) <- condition_names
  #my.labels <- c("Ctrl", "Gem")
  my.labels <- condition_names
  tmp2 <- PhyloObj %>% estimate_richness()
  new_names <- rownames(tmp2) %>% 
    gsub("^X", "", .) %>% 
    gsub("\\.", "-", .)
  rownames(tmp2) <- new_names
  rich_meta <- merge(PhyloObj %>% sample_data(), tmp2, by = "row.names")
  ##-- only using statistics test
  ob <- rich_meta %>% t_test(as.formula(paste0(index, " ~ ", strata))) %>% adjust_pvalue(method = "BH") %>%  add_significance("p.adj") %>% add_xy_position()
  
  p1 <- ggplot (rich_meta, aes(x = !!sym(strata), y = !!sym(index), fill=!!sym(strata)))+ 
    #geom_boxplot()+
    geom_violin() +
    geom_boxplot(width=0.2) +
    geom_point (position=position_jitterdodge( jitter.width = 0.05))+
    theme_gray() + 
    scale_x_discrete(labels= my.labels)+
    #xlab(strata)+
    ylab(y_label)+
    # facet_grid(.~tp)+
    #scale_x_discrete(labels= my.labels)+
    scale_fill_manual(values=group.colors, labels = my.labels)+
    ggtitle(gg_title) +
    #ggtitle("KPC tumor vs. Healthy pancreas - Observed species") +
    theme(axis.text.y = element_text (size=11),
          axis.title = element_text(size=12, face="bold"))
  
  if(!add_legend){    
    p1 <- p1 + theme(legend.position = "none",
                     # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                     legend.key.size = unit(4,"mm"),
                     #axis.text.x = element_blank(),
                     plot.title = element_text(size = 12))
    
  } else {
    p1 <- p1 + theme(legend.text = element_text(size = 12),
                     #legend.title = element_blank(),
                     # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                     legend.key.size = unit(4,"mm"),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12))
  }
  p1 <- p1 + stat_pvalue_manual(ob, label = "p.adj", inherit.aes = FALSE, tip.length = 0.01)
  return(p1)
}
##--
AlphaPlotWrapper_Violin <- function(PhyloObj, strata = NULL, roundUp = TRUE){
  ## round up/down otu_table
  if(roundUp){
    otu_table(PhyloObj) <- PhyloObj %>% otu_table() %>% round()
  } else{
    otu_table(PhyloObj) <- PhyloObj %>% otu_table() %>% ceiling()
  }
  plt.1 <- AlphaPlot_Violin(PhyloObj, index = "Observed", strata = strata, y_label = "Observed species", add_legend = F)
  plt.2 <- AlphaPlot_Violin(PhyloObj, index = "Shannon", strata = strata, y_label = "Shannon index", add_legend = F)
  plt.3 <- AlphaPlot_Violin(PhyloObj, index = "InvSimpson", strata = strata, y_label = "Inv Simpson index", add_legend = F)
  return(list("Observed" = plt.1, "Shannon" = plt.2, "InvSimpson" = plt.3))
}