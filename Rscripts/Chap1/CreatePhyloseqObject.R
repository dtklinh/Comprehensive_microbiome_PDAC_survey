### create phyloseq object for Chap 1
library(phyloseq)
library(microViz)
library(microbiome)
library(readxl)
library(tidyverse)

tax_tab <- readRDS("./data/Tax_table_v2.1.rds")
df_count <- read_csv("./data/Chap1/16S_ReClassified_2026April/count_table_AS1000Cov50.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "taxID")
df_meta <- read_xlsx("./meta/Chap1_NCT_16S_metadata.xlsx")
df_meta$Seq_date <- as.Date(df_meta$Seq_date, "%Y.%m.%d")

## edit meta
df_meta <- df_meta %>% 
  mutate(ID = sprintf("%s_%s", Run, Barcode)) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "ID")

## check if sample name in count and metadata agree
table(df_meta$ID %in% colnames(df_count)) 

ps <- phyloseq(otu_table(df_count, taxa_are_rows = T), sample_data(df_meta), tax_table(tax_tab))

## agg to species rank
ps <- ps %>% 
  subset_taxa(superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "_") %>% 
  microViz::tax_agg(rank = "species")

## rename row names in tax table
old_names <- taxa_names(ps)
new_names <- gsub(" ", "__", old_names)

taxa_names(ps) <- new_names

## save to RDS
saveRDS(ps, "./data/Chap1/NCT_2026April_v0.rds")

