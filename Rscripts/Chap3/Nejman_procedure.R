## redo the Nejman approach for selected 18 samples
library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)

source("add_code/Functions.R")
pseq <- readRDS("data/Chap3/pseq_Proj5_postFilter_v04.rds")
pseq_true <- pseq %>% 
  ps_filter(ffpe.bulk == "bulk") %>% 
  ps_filter(true.control == "true")
pseq_NCT <- pseq %>% 
  ps_filter(true.control == "NCT") %>% 
  ps_filter(ffpe.bulk == "bulk")
df1 <- stat_test(pseq_t = pseq_true, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("buffer")))
#df2 <- stat_test(pseq_t = pseq_true, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("pcr", "seq")))
df3 <- stat_test(pseq_t = pseq_true, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("buffer", "pcr", "seq")), prefix = "All")

df <- df1 %>% 
  left_join(., df3, by = c("Tax", "Prev_true")) %>% 
  rowwise() %>% 
  filter(max(c_across(contains("Fisher", ignore.case = TRUE))) <= 0.05) %>% 
  ungroup()

pseq_true_Nj <- pseq_true %>% 
  prune_taxa(taxa_names(.) %in% df$Tax , .) %>% 
  prune_samples(sample_sums(.) >0 ,.)
saveRDS(pseq_true_Nj, "data/Chap3/pseq_bulk_Fisher_v02.rds")
