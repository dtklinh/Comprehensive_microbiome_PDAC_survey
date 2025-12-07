## compare old vs new 10 samples
library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)
library(corncob)
library(skimr)
library(patchwork)

pseq_new <- readRDS("data/Chap3_Addition/pseq_Nejman.rds") %>% 
  ps_mutate(batch = "new")
lst_names <- taxa_names(pseq_new)
lst_names <- gsub(" ", "_", lst_names, )
taxa_names(pseq_new) <- lst_names

pseq_old <- readRDS("data/Chap3/pseq_bulk_Fisher_v02.rds") %>% 
  ps_filter(AN_NR %in% meta(pseq_new)$AN_NR) %>% 
  prune_taxa(taxa_sums(.) >0, .) %>% 
  ps_mutate(batch = "old")
pseq <- merge_phyloseq(pseq_new, pseq_old) %>% 
  ps_mutate(NumReads = sample_sums(.))
## remove AN866 due to low reads count

## Alpha, beta