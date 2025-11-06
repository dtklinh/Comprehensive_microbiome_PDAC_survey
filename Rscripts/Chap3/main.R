### Chapter 3 in manusctipt, mouse sample including FF and FFPE
library(phyloseq)
library(microbiome)
library(microViz)
library(vegan)
library(Wrench)
library(VennDiagram)
library(ANCOMBC)
library(tidyverse)
library(dplyr)

pseq_org <-         readRDS("./data/Chap3/pseq_Proj5_postFilter_v04.rds")
pseq_no <- pseq_org %>% 
  ps_filter(ffpe.bulk== "bulk") %>% 
  ps_filter(true.control == "true")
pseq_decontam <-    readRDS("./data/Chap3/pseq_bulk_decontam_p0.5.rds")
pseq_Fisher <-      readRDS("./data/Chap3/pseq_bulk_Fisher.rds")
pseq_restrictive <- readRDS("./data/Chap3/pseq_bulk_restrictive.rds")
pseq_SCRuB <-       readRDS("./data/Chap3/pseq_bulk_SCRuB.rds")

## load serial NCT
pseq_NCT <- readRDS("./data/Chap1/NCTs_v7_wrench.rds")

## investigate most prevalenace species in each methods

prev_thres <- 0.5
N <- 100
lst_SCRuB <-      get_most_prev(pseq_SCRuB, prev_threshold = prev_thres, N)
lst_Fischer <-    get_most_prev(pseq_Fisher, prev_thres, N)
lst_decontam <-   get_most_prev(pseq_decontam, prev_thres, N)
lst_restrictive <- get_most_prev(pseq_restrictive, prev_thres, N)
lst_NCT <-        get_most_prev(pseq_NCT, 0.1, 1000)

# sets <- list(
#   SCRuB = lst_SCRuB$taxa,
#   Fischer = lst_Fischer$taxa,
#   Decontam = lst_decontam$taxa,
#   Restrictive = lst_restrictive$taxa,
#   NCT = lst_NCT$taxa
# )
# sets2 <- list(
#   #NoDecon = taxa_names(pseq_no),
#   SCRuB = taxa_names(pseq_SCRuB),
#   Fischer = taxa_names(pseq_Fisher),
#   Decontam = taxa_names(pseq_decontam),
#   Restrictive = taxa_names(pseq_restrictive),
#   NCT = lst_NCT$taxa
# )
# venn.plot <- venn.diagram(
#   sets2,
#   filename = NULL,
#   #fill = c("red", "blue", "green", "pink", "yellow", "gray"),
#   alpha = 0.35,
#   cex = 1.5,
#   cat.cex = 1.2,
#   margin = 0.1
# )
# grid.draw(venn.plot)
# ggsave("./results/Chap3/VennDia_DecontamMethod_withNCT.png", plot = venn.plot, bg = "transparent")

## venn diagram with proportional area
#install.packages("eulerr")
# library(eulerr)
# 
# # Example set sizes:
# fit <- euler(sets)
# 
# plot(fit,
#      quantities = TRUE,
#      alpha = 0.5,
#      fills = c("skyblue", "pink", "lightgreen", "lightyellow", "grey"),
#      edges = TRUE)

###--------------------------------
### modify metadata for phyloseq
df_additionalInfo <- readxl::read_xlsx("./meta/Mice_meta.xlsx")

##------------------------------
## Work with decontam
pseq_decontam <- pseq_decontam %>% 
  append_AN_NR(df_additionalInfo)
pseq_decontam_WN <- pseq_decontam %>% 
  WrenchWrapper(grp = "Sex")
df_decontam_WN <- pseq_decontam_WN %>% 
  microbiome::alpha(index = c("observed", "diversity_shannon")) %>% 
  rownames_to_column(var = "SampleID")
df_read <- pseq_decontam_WN %>% 
  sample_sums() %>% 
  tibble(SampleID = names(.), reads = .)
df_decontam_WN <- df_decontam_WN %>% 
  left_join(., df_read, by = "SampleID")

otu_table <- abundances(pseq_decontam_WN)
# taxa_df <- tibble(
#   Sample = rep(colnames(otu_table), times = apply(otu_table, 2, function(x) sum(x > 0))),
#   Taxa = unlist(lapply(1:ncol(otu_table), function(i) rownames(otu_table)[otu_table[,i] > 0]))
# )
taxa_per_sample <- apply(otu_table, 2, function(x) {
  rownames(otu_table)[x/sum(x) > 0.01] 
  #tmp[tmp %in% lst_NCT] %>% length()
})
sample_taxa_df <- data.frame(
  SampleID = names(taxa_per_sample),
  Taxa_List = I(taxa_per_sample)  # I() preserves the list inside the column
) %>% 
  rowwise() %>% 
  mutate(NumAbund_1Per = length(Taxa_List),
    NumInNCT = sum(Taxa_List %in% lst_NCT$taxa)) %>% 
  ungroup() %>% 
  select(c("SampleID", "NumAbund_1Per", "NumInNCT"))
df_decontam_WN <- df_decontam_WN %>% 
  left_join(., sample_taxa_df, by ="SampleID")
