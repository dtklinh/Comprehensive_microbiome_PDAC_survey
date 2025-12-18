## decontamination procedures
library(phyloseq)
library(microViz)
library(microbiome)
library(decontam)
library(SCRuB)
library(tidyverse)

pseq <- readRDS("data/Chap3_Addition/pseq_origin_v0.0.rds")

pseq_2 <- pseq %>% 
  subset_taxa(superkingdom == "Bacteria") %>% 
  microViz::tax_fix(sep= "_") %>% 
  microViz::tax_agg(rank = "species")
  
## v0.1 saved
saveRDS(pseq_2, "data/Chap3_Addition/pseq_origin_v0.1.rds")
rm(list = ls())
pseq <- readRDS("data/Chap3_Addition/pseq_origin_v0.1.rds")

## remove noise, trivial noise both in true and nct samples
pseq_true <- pseq %>% 
  microViz::ps_filter(true.control == "TRUE")
pseq_NCT <- pseq %>% 
  ps_filter(true.control == "NCT")

pseq_true <- pseq_true %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 1e-6)
pseq_NCT <- pseq_NCT %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 5e-7)
pseq_merge <- merge_phyloseq(pseq_true, pseq_NCT)

## save v0.2
saveRDS(pseq_merge, "data/Chap3_Addition/pseq_origin_v0.2.rds")
rm(list = ls())
pseq <- readRDS("data/Chap3_Addition/pseq_origin_v0.2.rds")

pseq_true <- pseq %>% 
  microViz::ps_filter(true.control == "TRUE")
pseq_NCT <- pseq %>% 
  ps_filter(true.control == "NCT")

## restrictive
pseq <- readRDS("data/Chap3_Addition/pseq_origin_v1.0.rds")
pseq_true_restrictive <- pseq_true %>% 
  prune_taxa(!taxa_names(.) %in% taxa_names(pseq_NCT), .)
saveRDS(pseq_true_restrictive, "data/Chap3_Addition/pseq_restrictive.rds")
## decontam
library(decontam)
pseq <- pseq %>% 
  ps_mutate(is.neg = ifelse(true.control == "TRUE", FALSE, TRUE))

dfcontam_prev <- pseq %>% 
  decontam::isContaminant(method="prevalence", neg="is.neg", threshold=0.5)
taxa_Conts <- dfcontam_prev[dfcontam_prev$contaminant,] %>% rownames()
pseq_true_decontam <- pseq_true %>% 
  prune_taxa(!taxa_names(.) %in% taxa_Conts, .)
saveRDS(pseq_true_decontam, "data/Chap3_Addition/pseq_decontam.rds")

## SCRuB
s_data <- pseq %>% 
  abundances() %>% t()
s_metadata <- pseq %>% 
  meta() %>% 
  mutate(is_control = if_else(true.control == "TRUE", F, T)) %>% 
  select(is_control, sample_type = NCT_type) %>% 
  mutate(sample_type = replace_na(sample_type, "true_sample"))
s_out <- SCRuB(s_data, s_metadata, c("buffer", "pcr", "seq"))
otu_out <- s_out$decontaminated_samples %>% 
  as.data.frame() %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

pseq_true_scrub <- phyloseq(otu_out, sample_data(pseq_true), tax_table(pseq_true))
pseq_true_scrub <- pseq_true_scrub %>% 
  prune_taxa(taxa_sums(.) >0, .) %>% 
  tax_filter(min_prevalence = 1,
             prev_detection_threshold = 2,
             min_total_abundance = 1e-6)
saveRDS(pseq_true_scrub, "data/Chap3_Addition/pseq_SCRuB.rds")

#### Nejman

## determine a threshold for high prevalence filtering.
## return an array, how many percent of taxa in true sample I remove if I choose that threshold 
rm(list = ls())
## merge phyloseq object
##merge_phyloseq(ps1, ps2, ps3)
pseq <- readRDS("data/Chap3_Addition/pseq_origin_v1.0.rds")
pseq_true <- pseq %>% 
  ps_filter(true.control == "TRUE")
pseq_NCT <- pseq %>% 
  ps_filter(true.control == "NCT")
HighPrevalence_Data <- function(True.Sample, NCT.Sample){
  num_taxa <- ntaxa(True.Sample)
  # subset of taxa which in in overlap
  idx_overlap <- intersect(taxa_names(True.Sample) %>% unique(),
                           taxa_names(NCT.Sample) %>% unique())
  NCT.Sample.Subset <- prune_taxa(idx_overlap, NCT.Sample)
  prev.nct <- tibble(
    prev = prevalence(NCT.Sample.Subset, detection  = 0, sort = TRUE, count = FALSE),
    OTU =as.factor(names(prevalence(NCT.Sample.Subset , detection = 0, sort = TRUE, count = FALSE)))
  )
  x <- prev.nct$prev
  factorx <- factor(cut(x, breaks=nclass.Sturges(x)))
  xout <- as.data.frame(table(factorx)) %>% map_df(rev)
  xout <- mutate(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
  xout <- xout %>% mutate(rel_tax=cumFreq/num_taxa)
}
xout <- HighPrevalence_Data(pseq_true, pseq_NCT)
plot(xout$factorx, xout$rel_tax, xlab="prevalence in NCT samples", ylab="proportion of tax in true samples")

lst_NCT_highPrev <- pseq_NCT %>% 
  prevalence()
pseq_true_Nj <- pseq_true %>% 
  prune_taxa(!taxa_names(.) %in% names(lst_NCT_highPrev[lst_NCT_highPrev>=0.45]), .)
  
## Use binomial testing directly
## There are 2 bactches of NCT: not buffer, and not pcr

df1 <- stat_test(pseq_t = pseq_true_Nj, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("pcr", "seq")))
df2 <- stat_test(pseq_t = pseq_true_Nj, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("buffer", "seq")))
df3 <- stat_test(pseq_t = pseq_true_Nj, pseq_n = pseq_NCT, NCT_batch = list("NCT_type" = c("pcr", "buffer", "seq")), prefix = "All")
df_all <- df1 %>% 
  left_join(., df2, by = c("Tax", "Prev_true")) %>% 
  left_join(., df3, by = c("Tax", "Prev_true"))
df_all_filter <- df_all %>% 
  rowwise() %>% 
  filter(max(c_across(contains("fisher", ignore.case = TRUE))) <= 0.05) %>% 
  ungroup()
pseq_true_Nj <- pseq_true_Nj %>% 
  prune_taxa(taxa_names(.) %in% df_all_filter$Tax, .)
saveRDS(pseq_true_Nj, "data/Chap3_Addition/pseq_Nejman.rds")

##################
## clean FFPE sample with Nejman procedure
pseq <- readRDS("data/Chap3/pseq_Proj5_postFilter_v04.rds")
pseq_ffpe_true <- pseq %>% 
  ps_filter(ffpe.bulk == "FFPE") %>% 
  ps_filter(true.control == "true")
pseq_ffpe_nct <- pseq %>% 
  ps_filter(ffpe.bulk == "FFPE") %>% 
  ps_filter(true.control == "NCT")

xout <- HighPrevalence_Data(pseq_ffpe_true, pseq_ffpe_nct)
plot(xout$factorx, xout$rel_tax, xlab="prevalence in NCT samples", ylab="proportion of tax in true samples")
lst_NCT_highPrev <- pseq_ffpe_nct %>% 
  prevalence()
pseq_ffpe_true_Nj <- pseq_ffpe_true %>% 
  prune_taxa(!taxa_names(.) %in% names(lst_NCT_highPrev[lst_NCT_highPrev>=0.45]), .)

df1 <- stat_test(pseq_t = pseq_ffpe_true_Nj, pseq_n = pseq_ffpe_nct, NCT_batch = list("NCT_type" = c("paraffin", "buffer")))
df2 <- stat_test(pseq_t = pseq_ffpe_true_Nj, pseq_n = pseq_ffpe_nct, NCT_batch = list("NCT_type" = c("paraffin", "pcr")))
df3 <- stat_test(pseq_t = pseq_ffpe_true_Nj, pseq_n = pseq_ffpe_nct, NCT_batch = list("NCT_type" = c("paraffin", "pcr", "buffer")), prefix = "All")
df_all <- df1 %>% 
  left_join(., df2, by = c("Tax", "Prev_true")) %>% 
  left_join(., df3, by = c("Tax", "Prev_true"))
df_all_filter <- df_all %>% 
  rowwise() %>% 
  filter(max(c_across(contains("fisher", ignore.case = TRUE))) <= 0.1) %>% 
  ungroup()
