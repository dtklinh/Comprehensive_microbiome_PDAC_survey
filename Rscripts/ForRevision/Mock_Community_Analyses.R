### Work with Mock Community
library(phyloseq)
library(microViz)
library(microbiome)
library(dplyr)
library(openxlsx)

rm(list = ls())

## supporting functions
rename_samples <- function(ps, trans = "identity"){
  ps <- ps %>% 
    ps_mutate(ID = sample_names(.)) %>% 
    #ps_mutate(Replica = m_Replica, Method = m_Method) %>% 
    #ps_mutate(Group = sprintf("%s_%s", Method, Replica)) %>% 
    ps_mutate(NewID = sprintf("%s_%s", ID, Method)) %>% ps_arrange(Dilution) %>% 
    microbiome::transform(transform = trans)
  sample_names(ps) <- meta(ps)[["NewID"]]
  return(ps)
}

#pseq_original = readRDS("./MC_data/phyloseq/MC.trivialFilter.rds")
pseq_ND = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_ND.rds") %>% 
  #ps_filter(!is.NCT) %>% 
  ps_mutate(Method = "ND") %>%
  ps_mutate(Dilution = c(1,2,4,8,16,32,64,128, 256,512,1024,2048)) %>% 
  rename_samples()
pseq_restrictive = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_restrictive.rds") %>% 
  ps_mutate(Method = "Restrictive") %>% 
  ps_mutate(Dilution = c(1,2,4,8,16,32,64,128, 256,512,1024,2048)) %>% 
  rename_samples()
pseq_decontam = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_decontam.rds") %>% 
  ps_mutate(Method = "Decontam") %>% 
  ps_mutate(Dilution = c(1,2,4,8,16,32,64,128, 256,512,1024,2048)) %>% 
  rename_samples()
pseq_SCRuB = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_scrub.rds") %>% 
  ps_mutate(Method = "SCRuB") %>% 
  ps_mutate(Dilution = c(1,2,4,8,16,32,64,128, 256,512,1024,2048)) %>% 
  rename_samples()
pseq_Nj = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_Nj.rds") %>% 
  ps_mutate(Method = "Nj") %>% 
  ps_mutate(Dilution = c(1,2,4,8,16,32,64,128, 256,512,1024,2048)) %>% 
  rename_samples()

pseq_all = merge_phyloseq(pseq_ND, pseq_restrictive, pseq_decontam, pseq_SCRuB, pseq_Nj)

## barplot
my_order <- c("Nj", "Restrictive", "SCRuB", "Decontam", "ND")
sample_data(pseq_all)$Method <- factor(sample_data(pseq_all)$Method, levels = my_order) 
sample_data(pseq_all)$Dilution <- factor(sample_data(pseq_all)$Dilution, levels = c(1,2,4,8, 16,32,64, 128, 256,512,1024,2048)) 
plt <- pseq_all %>% 
  comp_barplot(
    tax_level = "species", n_taxa = 20, bar_width = 0.9, sample_order = "asis", # don't sort
    x = "Dilution" # x argument is available since microViz 0.9.7
  ) +
  facet_wrap(
    facets = vars(Method),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Sample ID", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + # axis starts exactly at 0
  theme_bw() + # slightly clearer axes for facets
  theme(panel.spacing.x = unit(6, "mm")) # space for non-overlapping axis numbers
ggsave("./results/MC_barplot_v02.png", plot = plt, width = 9, height = 6, dpi = 300)
plt

## count: recall - precision, total species

df_theori_abund = read.xlsx("./MC_data/TheoreticalAbund.xlsx")
pseq_ND = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_ND.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_restrictive = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_restrictive.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_decontam = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_decontam.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_SCRuB = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_SCRuB.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_Nj = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_Nj.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))

calc_metrics = function(ps_meth, ps_ND, dilution =1, df_true){
  
  ps_est = ps_meth %>% ps_filter(Dilution ==dilution)
  ps_org = ps_ND %>% ps_filter(Dilution == dilution)
  
  pred_pos = ps_est %>% taxa_names()
  pred_neg = setdiff(ps_org %>% taxa_names(), ps_est %>% taxa_names())
  actual_pos = df_theori_abund$species
  actual_neg = setdiff(ps_org %>% taxa_names(), actual_pos)
  
  TP = intersect(pred_pos, actual_pos) %>% length()
  FP = intersect(pred_pos, actual_neg) %>% length()
  FN = intersect(pred_neg, actual_pos) %>% length()
  TN = intersect(pred_neg, actual_neg) %>% length()
  num_species = pred_pos %>% length()
  
  Precision = TP/(TP+FP)
  Recall = TP/(TP+FN)
  F1 = 2*TP/(2*TP+FP+FN)
  MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  return(list(TP = TP, FP = FP, FN = FN, TN = TN, Precision = Precision, Recall = Recall, F1=F1, MCC = MCC, num_species = num_species))
}

Nj_tab = data.frame(Dilution = numeric(0),
                    TP = numeric(0), FP = numeric(0),
                    FN = numeric(0), TN = numeric(0),
                    Precision = numeric(0),
                    Recall = numeric(0),
                    F1 = numeric(0),
                    MCC = numeric(0),
                    Num_species = numeric(0))
lst_dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_Nj, pseq_ND, d, df_theori_abund)
  Nj_tab[nrow(Nj_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
Nj_tab = Nj_tab %>% mutate(Method = "Nj")

restrictive_tab = data.frame(Dilution = numeric(0),
                             TP = numeric(0), FP = numeric(0),
                             FN = numeric(0), TN = numeric(0),
                             Precision = numeric(0),
                             Recall = numeric(0),
                             F1 = numeric(0),
                             MCC = numeric(0),
                             Num_species = numeric(0))
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_restrictive, pseq_ND, d, df_theori_abund)
  restrictive_tab[nrow(restrictive_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
restrictive_tab = restrictive_tab %>% mutate(Method = "Restrictive")
##-----------------
decontam_tab = data.frame(Dilution = numeric(0),
                          TP = numeric(0), FP = numeric(0),
                          FN = numeric(0), TN = numeric(0),
                          Precision = numeric(0),
                          Recall = numeric(0),
                          F1 = numeric(0),
                          MCC = numeric(0),
                          Num_species = numeric(0))
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_decontam, pseq_ND, d, df_theori_abund)
  decontam_tab[nrow(decontam_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
decontam_tab = decontam_tab %>% mutate(Method = "Decontam")
##------------------------------------
SCRuB_tab = data.frame(Dilution = numeric(0),
                       TP = numeric(0), FP = numeric(0),
                       FN = numeric(0), TN = numeric(0),
                       Precision = numeric(0),
                       Recall = numeric(0),
                       F1 = numeric(0),
                       MCC = numeric(0),
                       Num_species = numeric(0))
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_SCRuB, pseq_ND, d, df_theori_abund)
  SCRuB_tab[nrow(SCRuB_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
SCRuB_tab = SCRuB_tab %>% mutate(Method = "SCRuB")
df = rbind(restrictive_tab, decontam_tab, SCRuB_tab, Nj_tab)
write.xlsx(df, "./results/Diluted_MC_Precision_Recall_v02.xlsx")

##-------------------------
removed_lst = c("Methanobrevibacter smithii", "Salmonella enterica", "Enterococcus faecalis", "Clostridium perfringens")
df_theori_abund = df_theori_abund %>% 
  dplyr::filter(!species %in% removed_lst)
pseq_ND = pseq_ND %>% 
  prune_taxa(!taxa_names(.) %in% removed_lst, .)
pseq_restrictive = pseq_restrictive %>% 
  prune_taxa(!taxa_names(.) %in% removed_lst, .)
pseq_decontam = pseq_decontam %>% 
  prune_taxa(!taxa_names(.) %in% removed_lst, .)
pseq_SCRuB = pseq_SCRuB %>% 
  prune_taxa(!taxa_names(.) %in% removed_lst, .)
pseq_Nj = pseq_Nj %>% 
  prune_taxa(!taxa_names(.) %in% removed_lst, .)

restrictive_tab = data.frame(Dilution = numeric(0),
                             TP = numeric(0), FP = numeric(0),
                             FN = numeric(0), TN = numeric(0),
                             Precision = numeric(0),
                             Recall = numeric(0),
                             F1 = numeric(0),
                             MCC = numeric(0),
                             Num_species = numeric(0))
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_restrictive, pseq_ND, d, df_theori_abund)
  restrictive_tab[nrow(restrictive_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
restrictive_tab = restrictive_tab %>% mutate(Method = "Restrictive")

decontam_tab = data.frame(Dilution = numeric(0),
                          TP = numeric(0), FP = numeric(0),
                          FN = numeric(0), TN = numeric(0),
                          Precision = numeric(0),
                          Recall = numeric(0),
                          F1 = numeric(0),
                          MCC = numeric(0),
                          Num_species = numeric(0))
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_decontam, pseq_ND, d, df_theori_abund)
  decontam_tab[nrow(decontam_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
decontam_tab = decontam_tab %>% mutate(Method = "Decontam")

SCRuB_tab = data.frame(Dilution = numeric(0),
                       TP = numeric(0), FP = numeric(0),
                       FN = numeric(0), TN = numeric(0),
                       Precision = numeric(0),
                       Recall = numeric(0),
                       F1 = numeric(0),
                       MCC = numeric(0),
                       Num_species = numeric(0))
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_SCRuB, pseq_ND, d, df_theori_abund)
  SCRuB_tab[nrow(SCRuB_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
SCRuB_tab = SCRuB_tab %>% mutate(Method = "SCRuB")

Nj_tab = data.frame(Dilution = numeric(0),
                    TP = numeric(0), FP = numeric(0),
                    FN = numeric(0), TN = numeric(0),
                    Precision = numeric(0),
                    Recall = numeric(0),
                    F1 = numeric(0),
                    MCC = numeric(0),
                    Num_species = numeric(0))
#lst_dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
for(d in lst_dilution){
  print(sprintf("Dilution: %i" , d))
  res = calc_metrics(pseq_Nj, pseq_ND, d, df_theori_abund)
  Nj_tab[nrow(Nj_tab) +1, ] = c(d, res$TP, res$FP, res$FN, res$TN, res$Precision, res$Recall, res$F1, res$MCC, res$num_species)
}
Nj_tab = Nj_tab %>% mutate(Method = "Nj")

df = rbind(restrictive_tab, decontam_tab, SCRuB_tab, Nj_tab)
write.xlsx(df, "./results/Diluted_MC_Precision_Recall_v03_Filter1PerAbund.xlsx")

### check species prevalence for each method
df_theori_abund = read.xlsx("./MC_data/TheoreticalAbund.xlsx")
pseq_ND = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_ND.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_restrictive = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_restrictive.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_decontam = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_decontam.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_SCRuB = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_SCRuB.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
pseq_Nj = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_Nj.rds") %>% 
  ps_mutate(Dilution = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))

tb_ND = pseq_ND %>% 
  abundances() %>% as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% df_theori_abund$species) %>% 
  rownames_to_column(var = "Taxa")
tb_restrictive = pseq_restrictive %>% 
  abundances() %>% as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% df_theori_abund$species) %>% 
  rownames_to_column(var = "Taxa")
tb_decontam = pseq_decontam %>% 
  abundances() %>% as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% df_theori_abund$species) %>% 
  rownames_to_column(var = "Taxa")
tb_SCRuB = pseq_SCRuB %>% 
  abundances() %>% as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% df_theori_abund$species) %>% 
  rownames_to_column(var = "Taxa")
tb_Nj = pseq_Nj %>% 
  abundances() %>% as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% df_theori_abund$species) %>% 
  rownames_to_column(var = "Taxa")
tbs = list(ND = tb_ND, Restrictive = tb_restrictive, Decontam = tb_decontam, SCRuB = tb_SCRuB, Nj = tb_Nj)
write.xlsx(x = tbs, "./results/MC_abundance_taxa_methods.xlsx", overwrite = TRUE)
