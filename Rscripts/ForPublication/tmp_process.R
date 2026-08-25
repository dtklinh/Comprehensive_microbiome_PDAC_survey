#### tmp process, will be deleted

library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)

## Chapter 1: check sample names

pseq_nct <- readRDS("./Rscripts/ForPublication/data/longitudinal_contaminant_survey/NCT_2026April_v0.rds")

#dirname(rstudioapi::getSourceEditorContext()$path) %>% print()
df_meta <- meta(pseq_nct)

## KPC mice
ps <- rprojroot::find_rstudio_root_file() %>% 
  file.path("data/Chap3", "pseq_Proj5_postFilter_v04.rds") %>% 
  readRDS()
df_meta <- meta(ps)

## FF replica
ff_replica <- rprojroot::find_rstudio_root_file() %>%
  file.path("data/Chap3_Addition/pseq_origin_v1.0.rds") %>% 
  readRDS()
meta_ff_rep <- meta(ff_replica)
