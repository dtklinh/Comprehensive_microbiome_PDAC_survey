#### tmp process, will be deleted

library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)

## Chapter 1: check sample names

pseq_nct <- readRDS("./Rscripts/ForPublication/data/longitudinal_contaminant_survey/NCT_2026April_v0.rds")

#dirname(rstudioapi::getSourceEditorContext()$path) %>% print()
df_meta <- meta(pseq_nct)
