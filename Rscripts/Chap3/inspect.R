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
