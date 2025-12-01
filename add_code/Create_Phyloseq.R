## create phyloseq onject from metapont
library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)

## True tumor
df_meta <- read_csv("data/Chap3_Addition/2025_11_07_16S_KPC_TT_UW_run01/2025_11_07_16S_KPC_TT_UW_run01_metadata.csv")
df_count <- read_csv("data/Chap3_Addition/2025_11_07_16S_KPC_TT_UW_run01/adj_count_table_AS1500Cov70.csv")

df_count <- df_count %>% 
  column_to_rownames(var = "taxID") %>% 
  rename_with(~paste0("2025_11_07_16S_KPC_TT_UW_run01_", .x))

tax_tab <- readRDS("./data/Tax_table_v2.1.rds")

df_meta <- df_meta %>% 
  column_to_rownames(var = "sample_id")

ps <- phyloseq(otu_table(df_count, taxa_are_rows = T), sample_data(df_meta), tax_tab)

## NCT
nct_meta <- read_csv("data/Chap3_Addition/2025_11_13_16S_KPC_TT_NTC_UW_run02/2025_11_13_16S_KPC_TT_NTC_UW_run02_metadata.csv")
nct_count <- read_csv("data/Chap3_Addition/2025_11_13_16S_KPC_TT_NTC_UW_run02/count_table_AS1250Cov70.csv")
nct_count <- nct_count %>% 
  select(-c(barcode10, barcode12)) %>% 
  column_to_rownames(var = "taxID") %>% 
  rename_with(~paste0("2025_11_13_16S_KPC_TT_NTC_UW_run02_", .x))
nct_meta <- nct_meta %>% 
  column_to_rownames(var = "sample_id")
ps_nct <- phyloseq(otu_table(nct_count, taxa_are_rows = T), sample_data(nct_meta), tax_tab)
ps_all <- merge_phyloseq(ps, ps_nct)
saveRDS(ps_all, "data/Chap3_Addition/pseq_origin_v0.0.rds")
