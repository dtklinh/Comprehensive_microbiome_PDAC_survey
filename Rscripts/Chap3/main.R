### Chapter 3 in manusctipt, mouse sample including FF and FFPE
library(phyloseq)
library(microbiome)
library(microViz)
library(vegan)
library(Wrench)
library(VennDiagram)
library(tidyverse)
library(dplyr)

pseq_org <-         readRDS("./data/Chap3/pseq_Proj5_postFilter_v04.rds")
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
lst_NCT <-        get_most_prev(pseq_NCT, prev_thres, N)

sets <- list(
  SCRuB = lst_SCRuB$taxa,
  Fischer = lst_Fischer$taxa,
  Decontam = lst_decontam$taxa,
  Restrictive = lst_restrictive$taxa,
  NCT = lst_NCT$taxa
)

venn.plot <- venn.diagram(
  sets,
  filename = NULL,
  fill = c("red", "blue", "green", "pink", "gray"),
  alpha = 0.35,
  cex = 1.5,
  cat.cex = 1.2,
  margin = 0.1
)
grid.draw(venn.plot)
ggsave("./results/Chap3/VennDia_DecontamMethod_withNCT.png", plot = venn.plot, bg = "transparent")

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

