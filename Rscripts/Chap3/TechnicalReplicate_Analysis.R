## compare old vs new 10 samples
library(phyloseq)
library(microViz)
library(microbiome)
library(tidyverse)
library(corncob)
library(skimr)
library(patchwork)

source("add_code/Functions.R")
pseq_new <- readRDS("data/Chap3_Addition/pseq_SCRuB.rds") %>% 
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

pseq <- pseq %>% 
  ps_filter(AN_NR != "AN866")

## Alpha
lst <- AlphaPlotWrapper_Violin(PhyloObj = pseq, strata = "batch", roundUp = F, m_paired = T)
## Beta

beta_plot_microViz(pseq, m_group = "batch")

## examine the overlaped species w.r.t 
abd_thres <- 0.001
lst_new <- pseq_new %>% 
  microViz::tax_filter(min_prevalence = 1, prev_detection_threshold = 2, min_total_abundance = abd_thres) %>% 
  taxa_names()
lst_old <- pseq_old %>% 
  microViz::tax_filter(min_prevalence = 1, prev_detection_threshold = 2, min_total_abundance = abd_thres) %>% 
  taxa_names()
