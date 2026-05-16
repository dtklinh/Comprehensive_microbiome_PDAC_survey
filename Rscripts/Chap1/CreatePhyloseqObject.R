### create phyloseq object for Chap 1
library(phyloseq)
library(microViz)
library(microbiome)
library(readxl)
library(tidyverse)

tax_tab <- readRDS("./data/Tax_table_v2.1.rds")
df_count <- read_csv("./data/Chap1/16S_ReClassified_2026April/count_table_AS1000Cov50.csv")
df_meta <- read_xlsx("./meta/Chap1_NCT_16S_metadata.xlsx")
df_meta$Seq_date <- as.Date(df_meta$Seq_date, "%Y.%m.%d")

## edit meta
