#### tmp process, will be deleted

library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)

## Chapter 1: check sample names

pseq_nct <- readRDS("./Rscripts/ForPublication/data/longitudinal_contaminant_survey/NCT_2026April_v0.rds")

## Chap1: get metadata
df_meta <- meta(pseq)
df_meta <- df_meta %>% 
  rownames_to_column(var = "ID") %>% 
  dplyr::select(ID, Run, Barcode, LP = Person, Control_type, Sequencing_date = Seq_date) 

library(openxlsx)
write.xlsx(df_meta, "Chap1.xlsx")
  

#dirname(rstudioapi::getSourceEditorContext()$path) %>% print()
df_meta <- meta(pseq_nct)

## KPC mice
ps <- rprojroot::find_rstudio_root_file() %>% 
  file.path("data/Chap3", "pseq_Proj5_postFilter_v04.rds") %>% 
  readRDS()
df_meta <- meta(ps)
df_meta_FF_R1 <- df_meta %>% 
  dplyr::filter(ffpe.bulk == "bulk") %>% 
  rownames_to_column(var = "ID")
write.xlsx(df_meta_FF_R1, "KPC_FF_R1.xlsx")
df_meta_FFPE <- df_meta %>% 
  dplyr::filter(ffpe.bulk == "FFPE") %>% 
  rownames_to_column(var = "ID")
write.xlsx(df_meta_FFPE, "KPC_FFPE.xlsx")

## FF replica
ff_replica <- rprojroot::find_rstudio_root_file() %>%
  file.path("data/Chap3_Addition/pseq_origin_v1.0.rds") %>% 
  readRDS()
meta_ff_rep <- meta(ff_replica) %>% 
  rownames_to_column(var = "ID")
write.xlsx(meta_ff_rep, "KPC_FF_R2.xlsx")
