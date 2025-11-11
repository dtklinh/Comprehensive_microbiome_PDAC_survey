## inspect SCRuB

Ts <- c("Peptonophilus_genitalis", "Duncaniella_sp._B8", "Lachnoclostridium_phocaeense")

tmp1 <- pseq_no %>% taxa_sums()
tmp1[names(tmp1) %in% Ts]

tmp2 <- pseq_restrictive %>% taxa_sums()
tmp2[names(tmp2) %in% Ts]
