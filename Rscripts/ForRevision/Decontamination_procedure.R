## decontamination procedures
library(phyloseq)
library(microViz)
library(microbiome)
library(decontam)
library(SCRuB)
library(openxlsx)
library(tidyverse)

## supportive functions:
source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "Decontam_Support_functions.R"))


## create phyloseq object from count table
count_tab <- read_csv("./MC_data/MoreDilution/count_table_truesample/count_table_AS1750Cov70.csv") %>% 
  column_to_rownames(var = "taxID") %>% 
  select(where(~ sum(.x, na.rm = TRUE) >= 100))
  

df_meta = read.xlsx("./MC_data/metadata/mock_community.xlsx", sheet = 2) %>% 
  mutate(dilution = c(1,2,4,8,16,32,64, rep(NA, 9)),
         is.NCT = c(rep(FALSE, 7), rep(TRUE, 9))) %>% 
  column_to_rownames(var = "barcode") %>% 
  dplyr::select(dilution, is.NCT, Control_type)

df_meta = tibble(barcode = paste0("barcode", c("01", "02", "03", "04", "05",
                                               "18", "19", "20", "21", "22",
                                               "23", "24")),
                 ID = paste0("Sample_", 1:12)) %>% 
  column_to_rownames(var = "barcode")
  
tax_tab = readRDS("../../data/Tax_table_v2.1.rds")

ps = phyloseq(otu_table(count_tab,taxa_are_rows = T), sample_data(df_meta), tax_table(tax_tab))
ps = ps %>% 
  subset_taxa(superkingdom = "Bacteria") %>% 
  microViz::tax_fix(sep= "__") %>% 
  microViz::tax_agg(rank = "species")

new_names = c("MC1:1", "MC1:2", "MC1:4", "MC1:8", "MC1:16", "MC1:32", "MC1:64",
              paste0("NCT_B_", 1:3), paste0("NCT_P_", 1:3), paste0("NCT_L_", 1:3))
sample_names(ps) <- new_names
  
## v0.1 saved
saveRDS(ps, "./MC_data/phyloseq/MC.original.rds")


## decontamination
## remove noise, trivial noise both in true and nct samples
ps_true <- ps %>% 
  microViz::ps_filter(!is.NCT)
ps_NCT <- ps %>% 
  ps_filter(is.NCT)

ps_true <- ps_true %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 1e-6)
ps_NCT <- ps_NCT %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 5e-7)
pseq_merge <- merge_phyloseq(ps_true, ps_NCT)

## save v0.2
saveRDS(pseq_merge, "./MC_data/phyloseq/MC.trivialFilter.rds")

##------------
## decontamination applied to pseq_origin_v0.2.rds
rm(list = ls())
pseq <- readRDS("./MC_data/phyloseq/MC.trivialFilter.rds")

pseq_true <- pseq %>% 
  microViz::ps_filter(!is.NCT)
pseq_NCT <- pseq %>% 
  ps_filter(is.NCT)

## Method: RESTRICTIVE
pseq <- readRDS("data/Chap3_Addition/pseq_origin_v0.2.rds")
pseq_true_restrictive <- pseq_true %>% 
  prune_taxa(!taxa_names(.) %in% taxa_names(pseq_NCT), .)
saveRDS(pseq_true_restrictive, "./MC_data/phyloseq/pseq_restrictive.rds")

## decontam
library(decontam)
pseq <- pseq %>% 
  ps_mutate(is.neg = ifelse(is.NCT == "TRUE", TRUE, FALSE))

dfcontam_prev <- pseq %>% 
  decontam::isContaminant(method="prevalence", neg="is.neg", threshold=0.5)
taxa_Conts <- dfcontam_prev[dfcontam_prev$contaminant,] %>% rownames()
pseq_true_decontam <- pseq_true %>% 
  prune_taxa(!taxa_names(.) %in% taxa_Conts, .)
saveRDS(pseq_true_decontam, "./MC_data/phyloseq/pseq_decontam.rds")

## SCRuB
s_data <- pseq %>% 
  abundances() %>% t()
s_metadata <- pseq %>% 
  meta() %>% 
  mutate(is_control = if_else(is.NCT == "TRUE", T, F)) %>% 
  dplyr::select(is_control, sample_type = Control_type) %>% 
  mutate(sample_type = replace_na(sample_type, "true_sample"))
## assume that we have NCT control types: buffer, pcr, seq
s_out <- SCRuB(s_data, s_metadata, c("C1", "C2", "C3"))
otu_out <- s_out$decontaminated_samples %>% 
  as.data.frame() %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

pseq_true_scrub <- phyloseq(otu_out, sample_data(pseq_true), tax_table(pseq_true))
pseq_true_scrub <- pseq_true_scrub %>% 
  prune_taxa(taxa_sums(.) >0, .) %>% 
  tax_filter(min_prevalence = 1,
             prev_detection_threshold = 2,
             min_total_abundance = 1e-6)
saveRDS(pseq_true_scrub, "./MC_data/phyloseq/pseq_SCRuB.rds")

#### Nejman

## determine a threshold for high prevalence filtering.
## return an array, how many percent of taxa in true sample I remove if I choose that threshold 
rm(list = ls())
## merge phyloseq object
##merge_phyloseq(ps1, ps2, ps3)
pseq <- readRDS("data/Chap3_Addition/pseq_origin_v0.2.rds")
pseq_true <- pseq %>% 
  ps_filter(true.control == "TRUE")
pseq_NCT <- pseq %>% 
  ps_filter(true.control == "NCT")


xout <- HighPrevalence_Data(pseq_true, pseq_NCT)
plt <- ggplot(xout, aes(x = factorx, y = rel_tax)) + 
  geom_col() +
  labs(
    x = "prevalence in NCT samples",
    y = "proportion of taxa in true samples"
  )
plot(xout$factorx, xout$rel_tax, xlab="prevalence in NCT samples", ylab="proportion of tax in true samples")
ggsave("./Manuscript/Figures_for_Manuscript_files/saved_png/3_Supplementary_Nj_curve.svg", plot = plt, width = 9, height = 8, dpi = 300, limitsize = F)

## from the curve/barplot, 0.45 seems to be a reasonable cutoff.
lst_NCT_highPrev <- pseq_NCT %>% 
  prevalence()
pseq_true_Nj <- pseq_true %>% 
  prune_taxa(!taxa_names(.) %in% names(lst_NCT_highPrev[lst_NCT_highPrev>=0.45]), .)
  
## Use binomial testing directly
## There are 2 bactches of NCT: not buffer, and not pcr

df1 <- stat_test(pseq_t = pseq_true_Nj, pseq_n = pseq_NCT, NCT_batch = list("Control_type" = c("C1", "C2")))
df2 <- stat_test(pseq_t = pseq_true_Nj, pseq_n = pseq_NCT, NCT_batch = list("Control_type" = c("C1", "C3")))
df3 <- stat_test(pseq_t = pseq_true_Nj, pseq_n = pseq_NCT, NCT_batch = list("Control_type" = c("C1", "C2", "C3")), prefix = "All")
df_all <- df1 %>% 
  left_join(., df2, by = c("Tax", "Prev_true")) %>% 
  left_join(., df3, by = c("Tax", "Prev_true"))
df_all_filter <- df_all %>% 
  rowwise() %>% 
  filter(max(c_across(contains("fisher", ignore.case = TRUE))) <= 0.05) %>% 
  ungroup()
pseq_true_Nj <- pseq_true_Nj %>% 
  prune_taxa(taxa_names(.) %in% df_all_filter$Tax, .)
saveRDS(pseq_true_Nj, "./MC_data/phyloseq/pseq_Nejman.rds")

##################
## clean FFPE sample with Nejman procedure
pseq <- readRDS("data/Chap3/pseq_Proj5_postFilter_v04.rds")
pseq_ffpe_true <- pseq %>% 
  ps_filter(ffpe.bulk == "FFPE") %>% 
  ps_filter(true.control == "true")
pseq_ffpe_nct <- pseq %>% 
  ps_filter(ffpe.bulk == "FFPE") %>% 
  ps_filter(true.control == "NCT")

xout <- HighPrevalence_Data(pseq_ffpe_true, pseq_ffpe_nct)
plot(xout$factorx, xout$rel_tax, xlab="prevalence in NCT samples", ylab="proportion of tax in true samples")
lst_NCT_highPrev <- pseq_ffpe_nct %>% 
  prevalence()
pseq_ffpe_true_Nj <- pseq_ffpe_true %>% 
  prune_taxa(!taxa_names(.) %in% names(lst_NCT_highPrev[lst_NCT_highPrev>=0.45]), .)

df1 <- stat_test(pseq_t = pseq_ffpe_true_Nj, pseq_n = pseq_ffpe_nct, NCT_batch = list("NCT_type" = c("paraffin", "buffer")))
df2 <- stat_test(pseq_t = pseq_ffpe_true_Nj, pseq_n = pseq_ffpe_nct, NCT_batch = list("NCT_type" = c("paraffin", "pcr")))
df3 <- stat_test(pseq_t = pseq_ffpe_true_Nj, pseq_n = pseq_ffpe_nct, NCT_batch = list("NCT_type" = c("paraffin", "pcr", "buffer")), prefix = "All")
df_all <- df1 %>% 
  left_join(., df2, by = c("Tax", "Prev_true")) %>% 
  left_join(., df3, by = c("Tax", "Prev_true"))
df_all_filter <- df_all %>% 
  rowwise() %>% 
  filter(max(c_across(contains("fisher", ignore.case = TRUE))) <= 0.1) %>% 
  ungroup()

#####---------------------------------------------------------
count_tab_true <- read_csv("./MC_data/MoreDilution/count_table_truesample/count_table_AS1250Cov70.csv") %>% 
  column_to_rownames(var = "taxID") %>% 
  select(where(~ sum(.x, na.rm = TRUE) >= 100))

ct_NCT2 = read_csv("./MC_data/MoreDilution/NCT/run02/count_table_AS1000Cov50.csv") %>% 
  column_to_rownames(var = "taxID") %>% 
  select(where(~ sum(.x, na.rm = TRUE) >= 100))
ct_NCT3 = read_csv("./MC_data/MoreDilution/NCT/run03/count_table_AS1000Cov50.csv") %>% 
  column_to_rownames(var = "taxID") %>% 
  select(where(~ sum(.x, na.rm = TRUE) >= 100))
meta_true = read.xlsx("./MC_data/MoreDilution/mock_community.xlsx", sheet = "Run04")
meta_nct2 = read.xlsx("./MC_data/MoreDilution/mock_community.xlsx", sheet = "Run02") %>% 
  mutate(Control_type = c(rep("Buffer", 3), rep("PCR", 3), "Seq"))
meta_nct3 = read.xlsx("./MC_data/MoreDilution/mock_community.xlsx", sheet = "Run03") %>% 
  mutate(Control_type = c(rep("Buffer", 2), rep("PCR", 2), rep("Seq", 2)))

tax_tab = readRDS("../../data/Tax_table_v2.1.rds")

meta_true = meta_true %>% 
  column_to_rownames(var = "barcode") %>% 
  mutate(is.neg = FALSE,
         Run = "Run04") %>% 
  dplyr::select(Dilution = pseudonym, is.neg, Run)
ps_true = phyloseq(otu_table(count_tab_true, taxa_are_rows = T), sample_data(meta_true), tax_tab)
sample_names(ps_true) = sprintf("TrueSample_%s", 1:12)
saveRDS(ps_true, "./MC_data/MoreDilution/phyloseq/pseq_true_v0.rds")

meta_nct2 = meta_nct2 %>% 
  mutate(ID = c("Buffer_1", "Buffer_2", "Buffer_3", "PCR_1", "PCR_2", "PCR_3", "Seq_1"),
         Run = "Run02",
         is.neg = rep(TRUE, 7)) %>% 
  column_to_rownames(var = "barcode") %>% 
  dplyr::select(Control_type, is.neg, Run, ID)

ps_NCT2 = phyloseq(otu_table(ct_NCT2, taxa_are_rows = TRUE), sample_data(meta_nct2), tax_tab)
sample_names(ps_NCT2) = meta(ps_NCT2)$ID
saveRDS(ps_NCT2, "./MC_data/MoreDilution/phyloseq/pseq_NCT2_v0.rds")

meta_nct3 = meta_nct3 %>% 
  mutate(ID = c("Buffer_4", "Buffer_5", "PCR_4", "PCR_5", "Seq_2", "Seq_3"),
         Run = "Run03",
         is.neg = rep(TRUE, 6)) %>% 
  column_to_rownames(var = "barcode") %>% 
  dplyr::select(Control_type, is.neg, Run, ID)
ps_NCT3 = phyloseq(otu_table(ct_NCT3, taxa_are_rows = TRUE), sample_data(meta_nct3), tax_tab)
sample_names(ps_NCT3) = meta(ps_NCT3)$ID
saveRDS(ps_NCT3, "./MC_data/MoreDilution/phyloseq/pseq_NCT3_v0.rds")

###---------------------------------------------------

ps_true = readRDS("./MC_data/MoreDilution/phyloseq/pseq_true_v0.rds")
ps_true = ps_true %>% 
  subset_taxa(!is.na(superkingdom) & superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "__") %>% 
  microViz::tax_agg(rank = "species") %>% 
  ps_get() %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 1e-6)
saveRDS(ps_true, "./MC_data/MoreDilution/phyloseq/pseq_true_ND.rds")

ps_NCT2 = readRDS("./MC_data/MoreDilution/phyloseq/pseq_NCT2_v0.rds")
ps_NCT2 = ps_NCT2 %>% 
  subset_taxa(!is.na(superkingdom) & superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "__") %>% 
  microViz::tax_agg(rank = "species") %>% 
  ps_get() %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 5e-7)

ps_NCT3 = readRDS("./MC_data/MoreDilution/phyloseq/pseq_NCT3_v0.rds")
ps_NCT3 = ps_NCT3 %>% 
  subset_taxa(!is.na(superkingdom) & superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "__") %>% 
  microViz::tax_agg(rank = "species") %>% 
  ps_get() %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 5e-7)

## restrictive
ps_true_restrictive = ps_true %>% 
  prune_taxa(!taxa_names(.) %in% taxa_names(ps_NCT2), .) %>%
  prune_taxa(!taxa_names(.) %in% taxa_names(ps_NCT3), .)
saveRDS(ps_true_restrictive, "./MC_data/MoreDilution/phyloseq/pseq_true_restrictive.rds")

## decontam
library(decontam)
pseq <- merge_phyloseq(ps_true, ps_NCT2, ps_NCT3)

dfcontam_prev <- pseq %>% 
  decontam::isContaminant(method="prevalence", neg="is.neg", threshold=0.5)
taxa_Conts <- dfcontam_prev[dfcontam_prev$contaminant,] %>% rownames()
ps_true_decontam <- ps_true %>% 
  prune_taxa(!taxa_names(.) %in% taxa_Conts, .)
saveRDS(ps_true_decontam, "./MC_data/MoreDilution/phyloseq/pseq_true_decontam.rds")

## SCRuB
s_data <- pseq %>% 
  abundances() %>% t()
s_metadata <- pseq %>% 
  meta() %>% 
  #mutate(is_control = if_else(is.NCT == "TRUE", T, F)) %>% 
  dplyr::select(is.neg, sample_type = Control_type) %>% 
  mutate(sample_type = replace_na(sample_type, "true_sample"))
## assume that we have NCT control types: buffer, pcr, seq
s_out <- SCRuB(s_data, s_metadata, c("Buffer", "PCR", "Seq"))
otu_out <- s_out$decontaminated_samples %>% 
  as.data.frame() %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

ps_true_scrub <- phyloseq(otu_out, sample_data(ps_true), tax_table(ps_true))
ps_true_scrub = ps_true_scrub %>% 
  prune_taxa(taxa_sums(.) >0, .)

saveRDS(ps_true_scrub, "./MC_data/MoreDilution/phyloseq/pseq_true_scrub.rds")

## Nejman

ps_NCT = merge_phyloseq(ps_NCT2, ps_NCT3)
xout <- HighPrevalence_Data(ps_true, ps_NCT)
plt <- ggplot(xout, aes(x = factorx, y = rel_tax)) + 
  geom_col() +
  labs(
    x = "prevalence in NCT samples",
    y = "proportion of taxa in true samples"
  )
plot(xout$factorx, xout$rel_tax, xlab="prevalence in NCT samples", ylab="proportion of tax in true samples")
ggsave("./Manuscript/Figures_for_Manuscript_files/saved_png/3_Supplementary_Nj_curve.svg", plot = plt, width = 9, height = 8, dpi = 300, limitsize = F)

## from the curve/barplot, 0.45 seems to be a reasonable cutoff.
lst_NCT_highPrev <- ps_NCT %>% 
  prevalence()
ps_true_Nj <- ps_true %>% 
  prune_taxa(!taxa_names(.) %in% names(lst_NCT_highPrev[lst_NCT_highPrev>=0.45]), .)

## Use binomial testing directly
## There are 2 bactches of NCT: not buffer, and not pcr

df1 <- stat_test(pseq_t = ps_true_Nj, pseq_n = ps_NCT, NCT_batch = list("Control_type" = c("Buffer", "Seq")))
df2 <- stat_test(pseq_t = ps_true_Nj, pseq_n = ps_NCT, NCT_batch = list("Control_type" = c("PCR", "Seq")))
df3 <- stat_test(pseq_t = ps_true_Nj, pseq_n = ps_NCT, NCT_batch = list("Control_type" = c("Buffer", "PCR", "Seq")), prefix = "All")
df_all <- df1 %>% 
  left_join(., df2, by = c("Tax", "Prev_true")) %>% 
  left_join(., df3, by = c("Tax", "Prev_true"))
df_all_filter <- df_all %>% 
  rowwise() %>% 
  filter(max(c_across(contains("fisher", ignore.case = TRUE))) <= 0.05) %>% 
  ungroup()
ps_true_Nj <- ps_true_Nj %>% 
  prune_taxa(taxa_names(.) %in% df_all_filter$Tax, .)
saveRDS(ps_true_Nj, "./MC_data/MoreDilution/phyloseq/pseq_true_Nj.rds")
