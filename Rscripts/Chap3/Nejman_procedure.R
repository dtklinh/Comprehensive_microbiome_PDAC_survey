## redo the Nejman approach for selected 18 samples
library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)

pseq <- readRDS("data/Chap3/pseq_Proj5_postFilter_v04.rds")
pseq_true <- pseq %>% 
  ps_filter(ffpe.bulk == "bulk") %>% 
  ps_filter(true.control == "true")
pseq_NCT <- pseq %>% 
  ps_filter(true.control == "NCT") %>% 
  ps_filter(ffpe.bulk == "bulk")
df1 <- stat_test(pseq_t = pseq_true, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("buffer")))
df2 <- stat_test(pseq_t = pseq_true, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("pcr", "seq")))
df3 <- stat_test(pseq_t = pseq_true, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("buffer", "pcr", "seq")), prefix = "All")
