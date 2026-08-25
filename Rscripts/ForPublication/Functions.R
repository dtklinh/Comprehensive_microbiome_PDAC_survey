## collection of all functions
Kolors <- c("red1", "dodgerblue3", "seagreen1", "#277257", "#9d547c","#56ca63", "navy", "#a357d6","cornflowerblue","#419d2a","sandybrown","red3","peachpuff","cyan","paleturquoise3","mistyrose","mediumpurple","mediumseagreen","mediumorchid","moccasin","orange4","olivedrab","midnightblue","papayawhip","palevioletred4","brown1","greenyellow","orchid","darkred","navajowhite1","mistyrose1","grey85","#525fd6","red2","#8cbe3a","#c944aa","indianred3","#5ba557","#9e66cb","#c1b735","#6d82ec","grey25","#e69728","#6654b0","lightsalmon3","lightcyan1","khaki1","plum1","lightsteelblue1","palevioletred3","mintcream","magenta3","#799330","#da7fdf","#3c782c","#e44586","blue4","#63c996","#dc3f53","#49cbc8","#cf3f29","#4fabda","#da6c2b","#598bd1","#b78c24","#8d4191","#a0b971","slategray1","sienna","plum1","lightyellow1","lightskyblue3","linen","limegreen","cornsilk1","mediumaquamarine","gray14","gold3","darkviolet","#b2386a","#479d71","#ae4341","#2ba198","#e07557","#5361a3","#dda353","#aa98df","#5b6114","#dc89bf","#327243","slateblue1","#e57b94","#9b62a0","#bbab59","#98495a","#526229","#d8827d","#857624","gray40","#9a4a22","#7c7d46","mediumslateblue","lemonchiffon1","#e3a073","#9e6b33", "gray74","slateblue1","rosybrown3", "lawngreen","gainsboro","deeppink3","firebrick3", "orchid2", "olivedrab1", "ivory3", "darkseagreen", "bisque2", "darkgoldenrod2", "blue2", "skyblue", "seashell2", "turquoise", "tan1", "seagreen2", "palevioletred3", "linen", "steelblue4","ghostwhite","dodgerblue1","deeppink1","firebrick1", "limegreen", "purple3", "khaki3", "snow3", "darkslategray","darkorchid","lavender", "magenta2", "palegreen", "salmon", "maroon", "cyan2","#671408","#FAEBD7","#7FFFD4","#F0FFFF","#A52A2A","burlywood","cadetblue","#7FFF00","chocolate","cornsilk","slateblue1","#FF7F50","#008B8B","darkgoldenrod1","darkolivegreen","darkorange4","white","hotpink","honeydew1","goldenrod2","darkgreen","oldlace","darkslategray3","navajowhite3","orchid4","gray25","#F0924D")
Kolors_Chap2 <- c("#9f5590", "#365c4d", "#d3edea", "#edf4fc")
##----------------------------------------------------------------------------
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
#####-------------------------------------------------------
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
      dplyr::select(c("AN_NR", new_col_names_fea))
  }
  
  df_meta_new <- df_meta %>%
    rownames_to_column(var = "SampleName") %>%
    left_join(., df_additional, by = "AN_NR") %>%
    column_to_rownames(var = "SampleName")
  sample_data(pseq) <- df_meta_new
  return(pseq)
}
##---------------------------------------------------------------------
CompositionalPlot_microViz <- function(PhyloseqObj, target_rank='order', facet = NULL, num = 10, tax_order = sum, my_title=""){
  PhyloseqObj %>%
    comp_barplot(
      tax_level = target_rank, n_taxa = num,
      tax_order = tax_order, other_name = "Other",
      #taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
      palette = distinct_palette(n = num, add = "grey90"),
      merge_other = FALSE, bar_outline_colour = "darkgrey"
    ) +
    coord_flip() +
    ggtitle(paste0(my_title))+
    facet_wrap(facet, nrow = 1, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 14),
          plot.title = element_text(size=16, face="bold"))
}
##---------------------------------------------------------------------------------
#### create comp_barplot function wrapping 3 properties
CompositionalPlot_Viz_v02 <- function(PhyloseqObj,
                                      grp = "Control_type",
                                      supper_rank='phylum',
                                      smaller_rank = "family",
                                      target_taxa_name = "Sphingomonadaceae",
                                      supper_rank_list = NULL,
                                      facet = NULL, each_num = 4, tax_order_by = sum,
                                      my_title=""){
  hueRank <- supper_rank
  hueRankPlural <- switch (hueRank,
                           "phylum" = "phyla",
                           "class" = "classes",
                           "order" = "orders",
                           "family" = "families",
                           "genus" = "genera"
  )
  shadeRank <- smaller_rank
  my_order <- supper_rank_list
  
  # Sort phyloseq at lower, and then higher ranks
  pseq2 <- PhyloseqObj %>%
    #ps_filter(gender == "male") %>%
    tax_sort(by = tax_order_by, at = shadeRank) %>%
    tax_sort(by = tax_order_by, at = hueRank) %>%
    tax_agg(rank = shadeRank)
  
  ## sort taxa
  tax_tbl <- tax_tibble(pseq2) %>%
    dplyr::arrange(match(.data[[hueRank]], my_order)) %>%
    remove_rownames() %>%
    column_to_rownames(var = "FeatureID")
  pseq2 <- pseq2 %>%
    tax_reorder(tax_order = rownames(tax_tbl))
  
  # Specify number of hues and shades desired
  nHues <- length(supper_rank_list) # "Other" phyla will be shades of grey
  nShades <- each_num # "Other" families will be the lightest shade of each hue
  
  hierarchicalPalInfo <- data.frame(
    hue = as.vector(tt_get(pseq2)[, hueRank]),
    shade = as.vector(tt_get(pseq2)[, shadeRank]),
    counts = taxa_sums(otu_get(pseq2))
  )
  
  hierarchicalPalInfo <- hierarchicalPalInfo %>%
    dplyr::mutate(
      hue = forcats::fct_other(
        f = hue, keep = unique(hue)[seq_len(nHues)],
        other_level = paste("Other", hueRankPlural)
      ),
      nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
    ) %>%
    dplyr::group_by(hue) %>%
    dplyr::mutate(
      shade = forcats::fct_other(
        f = shade, keep = unique(shade)[seq_len(nShades - 1)],
        other_level = "Other"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
      Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
    )
  
  hierarchicalPalMatrix <- matrix(
    data = sapply(
      X = seq(from = 30, to = 75, length.out = nShades),
      FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
    ),
    byrow = TRUE, ncol = nHues
  )
  hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))
  
  hierarchicalPal <- hierarchicalPalMatrix %>%
    as.vector() %>%
    setNames(unique(hierarchicalPalInfo$Taxa))
  
  col_nm <- sprintf("%s: %s", hueRank, shadeRank)
  plt <- pseq2 %>%
    ps_get() %>%
    #tax_mutate("Phylum: Family" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
    tax_mutate(!!col_nm := hierarchicalPalInfo$Taxa, .keep = "none") %>%
    tax_transform("compositional") %>%
    ps_arrange(desc(!!sym(target_taxa_name)), .target = "otu_table") %>%
    comp_barplot(
      tax_level = col_nm, n_taxa = length(hierarchicalPal),
      #tax_order = "asis",
      sample_order = "asis", counts_warn = FALSE,
      tax_order = "asis",
      palette = hierarchicalPal, bar_width = 0.975
    ) +
    coord_flip() +
    #ggtitle(paste0(my_title)) +
    #theme(legend.text = element_text(family = "mono")) # for text alignment
    facet_wrap(grp, nrow = 1, scales = "free") +
    labs(x = NULL, y = NULL,
         title = "Composition of NCT Samples",
         subtitle = my_title) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 14),
          plot.title = element_text(size=16, face="bold"))
  return(plt)
}

##--------------------------------------------------------------------------
## Heatmap for composition
HeatMap_MicroViz <- function(ps, grp = "Control_type",
                             rnk = "species",
                             tax_trans = "compositional",
                             tax_order_by = sum,
                             n_taxa = 20,
                             my_title = ""){
  #cols <- distinct_palette(n = 3, add = NA)
  cols_tmp <- Kolors[1:length(unique(meta(pseq)[["Control_type"]]))]
  names(cols_tmp) <- unique(meta(pseq)[["Control_type"]])
  p1 <- ps %>%
    ps_filter(!is.na(!!sym(grp))) %>%
    ps_arrange(!!sym(grp)) %>%
    tax_transform(tax_trans, rank = rnk) %>%
    comp_heatmap(grid_col = NA,
                 cluster_rows = FALSE, row_title = NULL,
                 taxa = tax_top(., n= n_taxa, by = tax_order_by),
                 #colors = heat_palette(sym = FALSE), name = "Abund.",
                 tax_anno = taxAnnotation(
                   Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(2.5, "cm")),
                   Abund. = anno_tax_box(size = grid::unit(2.5, "cm"))
                 ),
                 sample_anno = sampleAnnotation(
                   ControlType = anno_sample(grp),
                   col = list(ControlType = cols_tmp), 
                   border = FALSE
                   #Kolors[1:length(unique(meta(pseq)[["Control_type"]]))]
                   #TissueType = anno_sample("TissueType"),
                   #ControlType = anno_sample_cat(grp, col = Kolors[1:length(unique(meta(pseq)[["Control_type"]]))])
                   #col = list(ControlType = c("#9d547c","#56ca63","#a357d6")), border = F
                   #Early_Late = anno_sample("Stage")
                 ),
                 # group_rows = FALSE,
                 # group_cols = TRUE,
                 sample_seriation = "Identity",
                 row_names_gp = grid::gpar(fontsize = 12, fontfamily = "ArialMT"),
                 #width = grid::unit(9, "cm"),
                 rect_gp = grid::gpar(col = NA)
    )
  p1 <- p1 %>%
    ComplexHeatmap::draw(
      annotation_legend_list = attr(p1, "AnnoLegends"),
      column_title = my_title)
  return(p1)
}

## Alpha with glm
AlphaPlot_Violin_glm <- function(PhyloObj, index = "Observed", main_var = "", 
                                 covariate_var = "", x_label, y_label, add_legend = FALSE, 
                                 gg_title = NULL, gg_subtitle = NULL, m_Kolors = Kolors, m_family = gaussian, is_nb = F){
  condition_names <- meta(PhyloObj)[[main_var]] %>% table() %>% names()
  group.colors <- m_Kolors[1:length(condition_names)]
  #print(group.colors)
  names(group.colors) <- condition_names
  #my.labels <- c("Ctrl", "Gem")
  my.labels <- condition_names
  ## calculate alpha diversity
  
  tmp2 <- PhyloObj %>% estimate_richness()
  ##------------------------------
  new_names <- rownames(tmp2) %>% 
    gsub("^X", "", .) %>% 
    gsub("\\.", "-", .)
  rownames(tmp2) <- new_names
  rich_meta <- merge(PhyloObj %>% sample_data(), tmp2, by = "row.names")
  
  ##-------------------------
  y_pos <- t_test(formula = as.formula(sprintf("%s ~ %s", index, main_var)), data = rich_meta) %>% 
    add_xy_position() %>% 
    dplyr::select(group1, group2, y.position)
  if(!is_nb){
    model <- glm(formula = reformulate(termlabels = c(main_var, covariate_var), response = index), data = rich_meta, family = m_family)
  }else{
    model <- MASS::glm.nb(formula = reformulate(termlabels = c(main_var, covariate_var), response = index), data = rich_meta)
  }
  emm <- emmeans(model, specs = as.formula(sprintf("pairwise ~ %s", main_var)), adjust = "BH")
  contr <- contrast(emm, method = "pairwise")
  
  p_add <- contr %>% 
    as_tibble() %>%
    transmute(
      group1 = sub(" - .*", "", contrast),
      group2 = sub(".* - ", "", contrast),
      estimate,
      se = SE,
      statistic = if ("t.ratio" %in% names(.)) t.ratio else z.ratio,
      df,
      p = p.value
      ##method = "emmeans model-based t-test"
    ) %>% 
    add_significance() %>% 
    mutate(p_show = case_when(
      p < 0.001 ~ "<0.001***",
      p < 0.01 ~ sprintf("%.3f**", round(p,3)),
      p < 0.05 ~ sprintf("%.3f*", round(p,3)),
      TRUE ~ sprintf("%.3f", round(p,3))
    )) %>% 
    left_join(., y_pos, by = c("group1", "group2"))
  print(p_add)
  
  plt <- ggplot(rich_meta, aes(x = .data[[main_var]], y = .data[[index]], fill = .data[[main_var]])) +
    geom_violin() +
    #geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_boxplot(width=0.2) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
    # Connect lines for the same sample to show paired changes
    #geom_line(aes(group = .data[[SampleID]]), color = "gray", alpha = 0.5) +
    theme_minimal() +
    scale_fill_manual(values=group.colors) +
    labs(
      title = gg_title, #"Comparison of Decontamination Methods",
      subtitle = gg_subtitle, #sprintf("Alpha Diversity: %s", metric),
      y = y_label,
      x = x_label
    ) +
    theme_classic(base_family = "Arial") +
    theme(legend.text = element_text(family = "Arial", size = 12),
          #legend.title = element_blank(),
          # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
          legend.key.size = unit(4,"mm"),
          axis.text.x = element_blank(),
          plot.title = element_text(size = 12)
          #legend.position = "none"
    )
  #theme(legend.position = "none")
  plt <- plt + stat_pvalue_manual(p_add, label = "p_show", inherit.aes = FALSE, tip.length = 0.01)
  return(list(plt = plt, const = contr))
}

### Alpha with univariate analyses
AlphaPlot_Violin <- function(PhyloObj, index = "Observed", 
                             strata = "Control_type", y_label = "Observed species", 
                             add_legend = FALSE, gg_title = NULL, m_Kolors = Kolors, x_label = "Control type"){
  #group.colors <- c(ctrl = "dodgerblue3", Gem = "firebrick2")
  condition_names <- sample_data(PhyloObj)[[strata]] %>% table() %>% names()
  group.colors <- m_Kolors[1:length(condition_names)]
  print(group.colors)
  names(group.colors) <- condition_names
  #my.labels <- c("Ctrl", "Gem")
  my.labels <- condition_names
  ## calculate alpha diversity
  
  tmp2 <- PhyloObj %>% estimate_richness()
  ##------------------------------
  new_names <- rownames(tmp2) %>% 
    gsub("^X", "", .) %>% 
    gsub("\\.", "-", .)
  rownames(tmp2) <- new_names
  rich_meta <- merge(PhyloObj %>% sample_data(), tmp2, by = "row.names")
  ##-- only using statistics test
  ob <- rich_meta %>% 
    t_test(as.formula(paste0(index, " ~ ", strata))) %>% 
    adjust_pvalue(method = "BH") %>%  
    add_significance("p.adj") %>% 
    add_xy_position() %>% 
    mutate(p_show = case_when(
      p.adj < 0.001 ~ "<0.001***",
      p.adj < 0.01 ~ sprintf("%.3f**", round(p.adj,3)),
      p.adj < 0.05 ~ sprintf("%.3f*", round(p.adj,3)),
      TRUE ~ sprintf("%.3f(ns)", round(p.adj,3))
    ))
  #mutate(p_show = sprintf("%.3f(%s)", round(p.adj,3), p.adj.signif))
  
  p1 <- ggplot (rich_meta, aes(x = !!sym(strata), y = !!sym(index), fill=!!sym(strata)))+ 
    geom_violin() + #geom_boxplot()
    geom_boxplot(width=0.3) +
    geom_point (position=position_jitterdodge( jitter.width = 0.05), size = 0.3, col = alpha('grey', 0.2))+
    theme_minimal() + # theme_gray(), theme_bw()
    scale_x_discrete(labels= my.labels)+
    #xlab(strata)+
    ylab(y_label)+
    # facet_grid(.~tp)+
    #scale_x_discrete(labels= my.labels)+
    scale_fill_discrete(name = x_label) + 
    scale_fill_manual(values=group.colors, labels = my.labels)+ # my.labels
    xlab(x_label) + 
    ggtitle(gg_title) +
    #theme(panel.grid =  element_blank()) + 
    #ggtitle("KPC tumor vs. Healthy pancreas - Observed species") +
    theme(axis.text.y = element_text (size=12),
          axis.title = element_text(size=12, face="bold"))
  
  if(!add_legend){    
    p1 <- p1 + theme(legend.position = "none",
                     # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                     legend.key.size = unit(4,"mm"),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12))
    
  } else {
    p1 <- p1 + theme(legend.text = element_text(size = 12),
                     #legend.title = element_blank(),
                     # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                     legend.key.size = unit(4,"mm"),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12))
  }
  p1 <- p1 + stat_pvalue_manual(ob, label = "p_show", inherit.aes = FALSE, tip.length = 0.01, step.increase = 0.1)
  return(p1)
}
## AlphaPlot_Violin version 2 with glm to calculate OR and p-value

##--
AlphaPlotWrapper_Violin <- function(PhyloObj, strata = "Control_type", roundUp = TRUE, m_Kolors = Kolors, ...){
  ## round up/down otu_table
  if(roundUp){
    otu_table(PhyloObj) <- PhyloObj %>% otu_table() %>% round()
  } else{
    otu_table(PhyloObj) <- PhyloObj %>% otu_table() %>% ceiling()
  }
  plt.1 <- AlphaPlot_Violin(PhyloObj, index = "Observed", strata = strata, y_label = "Observed species", add_legend = T, m_Kolors = m_Kolors, ...)
  plt.2 <- AlphaPlot_Violin(PhyloObj, index = "Shannon", strata = strata, y_label = "Shannon index", add_legend = T, m_Kolors = m_Kolors, ...)
  plt.3 <- AlphaPlot_Violin(PhyloObj, index = "InvSimpson", strata = strata, y_label = "Inv Simpson index", add_legend = T, m_Kolors = m_Kolors, ...)
  return(list("Observed" = plt.1, "Shannon" = plt.2, "InvSimpson" = plt.3))
}

#### BETA
#### without any blocked permutation
BetaPlot <- function(PhyloObjct, strata_f , n_per = 9999, title_method = "Original", dis_method = "bray", ordination_method = "PCoA"){
  bray_dist = phyloseq::distance(PhyloObjct, method=dis_method)
  ordination = ordinate(PhyloObjct, method=ordination_method, distance=bray_dist)
  
  pcoa1 <- paste("PCoA 1 [", round(ordination[[3]]$Relative_eig[1], digits = 3)*100, "%]", sep = "")
  pcoa2 <- paste("PCoA 2 [", round(ordination[[3]]$Relative_eig[2], digits = 3)*100, "%]", sep = "")
  
  ##p.adonis <- adonis2(bray_dist ~ sample_data(PhyloObjct)$Treatment)
  p.adonis <- adonis2(as.formula(paste0("bray_dist ~ ", strata_f)), data = PhyloObjct %>% sample_data() %>% as_tibble(), permutations = n_per)
  
  p <- case_when(
    p.adonis$`Pr(>F)`[1] > 0.05 ~ paste("p =", round(p.adonis$`Pr(>F)`[1],4), "n.s.", sep = " "),
    p.adonis$`Pr(>F)`[1] < 0.05 &  p.adonis$`Pr(>F)`[1] > 0.01 ~ paste("p =", round(p.adonis$`Pr(>F)`[1],4), "*", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.01 & p.adonis$`Pr(>F)`[1] > 0.001  ~ paste("p =", round(p.adonis$`Pr(>F)`[1],4), "**", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.001 ~ paste("p =",round(p.adonis$`Pr(>F)`[1],4), "***", sep = " "),
  )
  
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = p,
    hjustvar = c(-0.2) ,
    vjustvar = c(1.5))
  
  p1 <- plot_ordination(PhyloObjct, ordination, color = strata_f) +
    #geom_point(aes(colour = !!sym(strata_f), shape = !!sym(type)), size = 3) +
    geom_point(aes(colour = !!sym(strata_f)), size = 3) +
    #geom_point(aes(colour = .data[[strata_f]], shape = .data[[type]]), size = 3) +
    #geom_point(aes(colour = sym(strata_f), shape = sym(type)), size = 3)
    ##geom_text_repel(aes(label = id), size = 4) +
    theme(aspect.ratio=1) +
    #theme_bw()+
    theme_minimal() +
    scale_color_brewer(palette = "Set1")+
    stat_ellipse() +
    xlab(pcoa1)+
    ylab(pcoa2)+
    #theme(panel.grid =  element_blank())+
    ggtitle(paste0(title_method, "")) +
    theme (axis.text=element_text(size=14),
           axis.title=element_text(size=16,face="bold"),
           legend.text = element_text(size = 12),
           legend.title = element_blank())+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4.5, inherit.aes = FALSE)
  return(p1)
}

#### pairwise comparisons with beta plot
BetaPlot_pairwise <- function(PhyloObjct, grp = "Replica", n_per = 9999, 
                              title_method = "Original", dis_method = "bray", ordination_method = "PCoA"){
  
  df_meta <- meta(PhyloObjct) 
  Y <- abundances(PhyloObjct) %>% t()
  #Y <- decostand(Y, method = "hellinger")
  my_dist = vegan::vegdist(Y, method = dis_method)
  df_meta[[grp]] <- factor(df_meta[[grp]])
  #bray_dist = phyloseq::distance(PhyloObjct, method=dis_method)
  
  p.adonis <- pairwise.adonis2(x = reformulate(termlabels = c(grp), response = "my_dist"), data = df_meta, permutations = n_per)
  
  return(p.adonis)
}

### ordinate plot by manual
BetaPlot_v02 <- function(PhyloObjct, strata_f = NULL, color_f = NULL, shape_f = NULL, 
                         ellipse_f = NULL, n_per = 9999, title_method = "Original", dis_method = "bray", ordination_method = "PCoA"){
  bray_dist = phyloseq::distance(PhyloObjct, method=dis_method)
  ordination = ordinate(PhyloObjct, method=ordination_method, distance=bray_dist)
  
  pcoa1 <- paste("PCoA 1 [", round(ordination[[3]]$Relative_eig[1], digits = 3)*100, "%]", sep = "")
  pcoa2 <- paste("PCoA 2 [", round(ordination[[3]]$Relative_eig[2], digits = 3)*100, "%]", sep = "")
  
  ##p.adonis <- adonis2(bray_dist ~ sample_data(PhyloObjct)$Treatment)
  p.adonis <- adonis2(as.formula(paste0("bray_dist ~ ", strata_f)), data = PhyloObjct %>% sample_data() %>% as_tibble(), permutations = n_per)
  
  p <- case_when(
    p.adonis$`Pr(>F)`[1] > 0.05 ~ paste("p =", round(p.adonis$`Pr(>F)`[1],4), "n.s.", sep = " "),
    p.adonis$`Pr(>F)`[1] < 0.05 &  p.adonis$`Pr(>F)`[1] > 0.01 ~ paste("p =", round(p.adonis$`Pr(>F)`[1],4), "*", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.01 & p.adonis$`Pr(>F)`[1] > 0.001  ~ paste("p =", round(p.adonis$`Pr(>F)`[1],4), "**", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.001 ~ paste("p =",round(p.adonis$`Pr(>F)`[1],4), "***", sep = " "),
  )
  
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = p,
    hjustvar = c(-0.2) ,
    vjustvar = c(1.5))
  
  p1 <- plot_ordination(PhyloObjct, ordination, color = color_f, shape = shape_f) +
    #geom_point(aes(colour = !!sym(strata_f), shape = !!sym(type)), size = 3) +
    geom_point(aes(colour = !!sym(color_f)), size = 3) +
    #geom_point(aes(colour = .data[[strata_f]], shape = .data[[type]]), size = 3) +
    #geom_point(aes(colour = sym(strata_f), shape = sym(type)), size = 3)
    ##geom_text_repel(aes(label = id), size = 4) +
    theme(aspect.ratio=1) +
    #theme_bw()+
    theme_minimal() +
    #scale_color_brewer(palette = "Set1")+
    scale_color_manual(values = c("LP1" = "#9f5590", "LP2" = "#365c4d")) +
    #stat_ellipse() +
    stat_ellipse(
      aes(group = .data[[ellipse_f]]),
      geom = "polygon",
      alpha = 0.0,
      color = "black"
      #inherit.aes = FALSE
    ) +
    xlab(pcoa1)+
    ylab(pcoa2)+
    #theme(panel.grid =  element_blank())+
    ggtitle(paste0(title_method, "")) +
    theme (axis.text=element_text(size=14),
           axis.title=element_text(size=16,face="bold"),
           legend.text = element_text(size = 12),
           legend.title = element_blank())+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4.5, inherit.aes = FALSE)
  return(p1)
}

## alpha plot function with LMM with fixed and random variables
AlphaPlot_Violin_LMM <- function(df, SampleID, strata, val, y_lab = "Obsered species", m_title = "", m_subtitle = "", ...){
  ## df: in pivot longer manner
  ## calc p value in manner of paired
  ## get y.position for ggplot
  y_pos <- t_test(formula = as.formula(sprintf("%s ~ %s", val, strata)), data = df, ...) %>% add_xy_position() %>% dplyr::select(group1, group2, y.position)
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
    left_join(., y_pos, by = c("group1", "group2"))
  ##mutate(y.position = y_pos)
  plt <- ggplot(df, aes(x = .data[[strata]], y = .data[[val]], fill = .data[[strata]])) +
    geom_violin() +
    #geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_boxplot(width=0.2) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
    # Connect lines for the same sample to show paired changes
    geom_line(aes(group = .data[[SampleID]]), color = "gray", alpha = 0.5) +
    theme_minimal() +
    labs(
      title = m_title, #"Comparison of Decontamination Methods",
      subtitle = m_subtitle, #sprintf("Alpha Diversity: %s", metric),
      y = y_lab
    ) +
    theme(legend.position = "none")
  plt <- plt + stat_pvalue_manual(p_add, label = "p.signif", inherit.aes = FALSE, tip.length = 0.01)
  return(list(plt = plt, const = contr))
}
##-----------------------------------------------------------------------------------------------------------
## alpha plot with Wilcon_test
AlphaPlot_Violin_Wilcox <- function(df, SampleID, strata, val, y_label = "Obsered species", title = "Paired Wilcox test", subtitle = "The higher the better", ...){
  ## df: in pivot longer manner
  ## calc p value in manner of paired
  ## get y.position for ggplot
  p_add <- wilcox_test(formula = as.formula(sprintf("%s ~ %s", val, strata)), data = df, paired = TRUE, p.adjust.method = "BH", ...) %>%
    add_significance() %>%  
    add_xy_position() %>% ##%>% pull(y.position)
    mutate(p_show = case_when(
      p.adj < 0.001 ~ "<0.001***",
      p.adj < 0.01 ~ sprintf("%.3f**", round(p.adj,3)),
      p.adj < 0.05 ~ sprintf("%.3f*", round(p.adj,3)),
      TRUE ~ sprintf("%.3f", round(p.adj,3))
    )) 
  
  plt <- ggplot(df, aes(x = .data[[strata]], y = .data[[val]], fill = .data[[strata]])) +
    geom_violin() +
    #geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_boxplot(width=0.2) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
    # Connect lines for the same sample to show paired changes
    geom_line(aes(group = .data[[SampleID]]), color = "gray", alpha = 0.5) +
    theme_minimal() +
    labs(
      title = title, #"Paired Wilcox test",
      subtitle = subtitle,
      y = y_label
    ) +
    theme(legend.position = "none")
  plt <- plt + stat_pvalue_manual(p_add, label = "p_show", inherit.aes = FALSE, tip.length = 0.01)
  return(plt)
}
get_most_prev <- function(pseq, prev_threshold = 0.5, N = 25){
  res <- pseq %>% 
    prevalence() %>% 
    tibble(taxa = names(.), value = .) %>% 
    filter(value >= prev_threshold) %>% 
    arrange(desc(value)) %>% 
    head(N)
  return(res)
}

## Alpha plot with pairwise_wilcox_test, for all pairs, but shown for Nj-*
AlphaPlot_Violin_Wilcox_2 <- function(df, SampleID, strata, val, y_label = "Obsered species", title = "Paired Wilcox test", subtitle = "The higher the better", my_alternative = "greater", comparisons = list(c("Nj", "ND"))){
  ## df: in pivot longer manner
  ## calc p value in manner of paired
  ## get y.position for ggplot
  p_add <- pairwise_wilcox_test(formula = as.formula(sprintf("%s ~ %s", val, strata)), data = df, paired = TRUE, p.adjust.method = "BH", alternative = my_alternative) %>%
    dplyr::filter(group1 == "Nj") %>% 
    #add_significance() %>%  
    add_xy_position() %>% ##%>% pull(y.position)
    mutate(p_show = case_when(
      p.adj < 0.001 ~ "<0.001***",
      p.adj < 0.01 ~ sprintf("%.3f**", round(p.adj,3)),
      p.adj < 0.05 ~ sprintf("%.3f*", round(p.adj,3)),
      TRUE ~ sprintf("%.3f", round(p.adj,3))
    )) 
  
  plt <- ggplot(df, aes(x = .data[[strata]], y = .data[[val]], fill = .data[[strata]])) +
    geom_violin() +
    #geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_boxplot(width=0.2) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
    # Connect lines for the same sample to show paired changes
    geom_line(aes(group = .data[[SampleID]]), color = "gray", alpha = 0.5) +
    theme_minimal() +
    labs(
      title = title, #"Paired Wilcox test",
      subtitle = subtitle,
      y = y_label
    ) +
    theme(legend.position = "none")
  plt <- plt + stat_pvalue_manual(p_add, label = "p_show", inherit.aes = FALSE, tip.length = 0.01)
  return(plt)
}
get_most_prev <- function(pseq, prev_threshold = 0.5, N = 25){
  res <- pseq %>% 
    prevalence() %>% 
    tibble(taxa = names(.), value = .) %>% 
    filter(value >= prev_threshold) %>% 
    arrange(desc(value)) %>% 
    head(N)
  return(res)
}

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
      mutate(NumRemain = length(Taxa_List),
             NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
      ungroup() %>% 
      dplyr::select(SampleID, NumRemain, NumInNCT)
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
      mutate(NumRemain = length(Taxa_List),
             NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
      ungroup() %>% 
      dplyr::select(SampleID, NumRemain, NumInNCT)
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
      mutate(NumRemain = length(Taxa_List),
             NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
      ungroup() %>% 
      dplyr::select(SampleID, NumRemain, NumInNCT)
  }
  
  df_final <- df_final %>% 
    left_join(., sample_taxa_df, by ="SampleID")
  return(df_final)
}
process_survey <- function(psLst, lstNCT, m_by = "prev", thres_abund = 0.01, thres_prev = 0.5){
  df_all <- tibble(SampleID = character(), observed = numeric(), 
                   NumRemain = numeric(), NumInNCT = numeric(), Method = character())
  for(i in 1:length(psLst)){
    xxx <- survey_NCT(psLst[[i]], lst_NCT = lstNCT, by = m_by, 
                      thres_abd = thres_abund, thres_prev = thres_prev) %>% 
      mutate(Method = names(psLst[i])) %>% 
      dplyr::select(SampleID, observed, NumRemain, NumInNCT, Method)
    df_all <- df_all %>% 
      bind_rows(tibble(xxx))
  }
  df_all <- df_all %>% 
    mutate(Clean_Count = NumRemain - NumInNCT,
           Yield = Clean_Count/observed,
           Purity = ifelse(NumRemain > 0, Clean_Count / NumRemain, 0),
           Composite_Score = Yield * Purity)
  return(df_all)
}
#### create square table, where N_i,j tell us if Method i is less than Method j
pairwise_WT_Method <- function(df_longer, metric = "Composite_Score", strata = "Method", m_alternative = "less", m_order = c("Nj", "Restrictive", "SCRuB", "Decontam", "ND")){
  ## df: df_longer format --> AN_NR  Methods    Val
  #m_order <- c("Nj", "Restrictive", "SCRuB", "Decontam", "Raw")
  all_groups <- m_order
  
  df_longer[[strata]] <- factor(df_longer[[strata]], levels = m_order)
  res_upper <- pairwise_wilcox_test(data = df_longer, formula = as.formula(sprintf("%s ~ %s", metric, strata)), paired = TRUE, alternative = m_alternative, p.adjust.method = "BH") %>% 
    mutate(val = sprintf("%s(%s)", p.adj, p.adj.signif))
  
  df_longer[[strata]] <- factor(df_longer[[strata]], levels = rev(m_order))
  res_lower <- pairwise_wilcox_test(data = df_longer, formula = as.formula(sprintf("%s ~ %s", metric, strata)), paired = TRUE, alternative = m_alternative, p.adjust.method = "BH") %>% 
    mutate(val = sprintf("%s(%s)", p.adj, p.adj.signif))
  grid <- expand.grid(
    group1 = all_groups,
    group2 = all_groups,
    stringsAsFactors = FALSE
  )
  full <- grid %>%
    left_join(res_upper %>% dplyr::select(group1, group2, p_less = val), 
              by = c("group1", "group2")) %>%
    left_join(res_lower %>% dplyr::select(group1, group2, p_greater = val), 
              by = c("group1", "group2"))
  full <- full %>%
    mutate(value = case_when(
      !is.na(p_less) ~ as.character(p_less),        # i < j
      !is.na(p_greater) ~ as.character(p_greater),  # i > j
      TRUE ~ NA_character_
    ))
  # Pivot to wide format
  mat_df <- full %>%
    dplyr::select(group1, group2, value) %>%
    pivot_wider(names_from = group2, values_from = value)
  # convert to matrix
  final_mat <- mat_df %>%
    column_to_rownames("group1") %>%
    as.matrix()
  return(final_mat)
}