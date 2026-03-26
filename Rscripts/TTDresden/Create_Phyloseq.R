## create phyloseq onject from metapont
library(phyloseq)
library(microViz)
library(microbiome)
library(readxl)
library(tidyverse)

## True tumor
df_meta <- rprojroot::find_rstudio_root_file() %>% 
  file.path("meta/TTDresden.xlsx") %>% 
  readxl::read_xlsx()

df_count <- rprojroot::find_rstudio_root_file() %>% 
  file.path("data/TTDresden/count_table_AS1500Cov70.csv") %>% 
  read_csv()

df_count <- df_count %>% 
  column_to_rownames(var = "taxID")

tax_tab <- rprojroot::find_rstudio_root_file() %>% 
  file.path("data/Tax_table_v2.1.rds") %>% 
  readRDS()

df_meta <- df_meta %>% 
  mutate(SampleID = sprintf("%s_%s", PDAC_ID, TissueType),
         TRUE.CTRL = "TrueSample") %>% 
  column_to_rownames(var = "Barcode")

ps <- phyloseq(otu_table(df_count, taxa_are_rows = T), sample_data(df_meta), tax_tab)
ps <- ps %>% 
  subset_taxa(superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "_") %>% 
  microViz::tax_agg(rank = "species") %>% 
  ps_get() %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 1e-6)
new_sample_names <- meta(ps)[["SampleID"]]
sample_names(ps) <- new_sample_names
ps <- ps %>% 
  ps_select(-SampleID)

## NCT
nct_meta <- rprojroot::find_rstudio_root_file() %>% 
  file.path("meta/TTDresden_NCT.xlsx") %>% 
  read_xlsx()
nct_count <- rprojroot::find_rstudio_root_file() %>% 
  file.path("data/TTDresden_NCT/count_table_AS1500Cov70.csv") %>% 
  read_csv()

nct_count <- nct_count %>% 
  column_to_rownames(var = "taxID")
  
nct_meta <- nct_meta %>% 
  column_to_rownames(var = "Barcode") %>% 
  mutate(TRUE.CTRL = "CTRL")

ps_nct <- phyloseq(otu_table(nct_count, taxa_are_rows = T), sample_data(nct_meta), tax_tab)

ps_nct <- ps_nct %>% 
  subset_taxa(superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "_") %>% 
  microViz::tax_agg(rank = "species") %>% 
  ps_get() %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 5e-7)

new_sample_names <- meta(ps_nct)[["SampleID"]]
sample_names(ps_nct) <- new_sample_names
ps_nct <- ps_nct %>% 
  ps_select(-SampleID)

ps_all <- merge_phyloseq(ps, ps_nct)

## rename taxa
old_taxa_names <- taxa_names(ps_all)
new_taxa_names <- sapply(old_taxa_names, function(x){gsub("\\s+", "_",x)})
taxa_names(ps_all) <- new_taxa_names
ps_all
saveRDS(ps_all, "./Rscripts/TTDresden/pseq.rds")
