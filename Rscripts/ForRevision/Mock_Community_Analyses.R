### Work with Mock Community
library(phyloseq)
library(microViz)
library(microbiome)

rm(list = ls())

## supporting functions
rename_samples <- function(ps, trans = "identity"){
  ps <- ps %>% 
    ps_mutate(ID = sample_names(.)) %>% 
    #ps_mutate(Replica = m_Replica, Method = m_Method) %>% 
    #ps_mutate(Group = sprintf("%s_%s", Method, Replica)) %>% 
    ps_mutate(NewID = sprintf("%s_%s", ID, Method)) %>% ps_arrange(dilution) %>% 
    microbiome::transform(transform = trans)
  sample_names(ps) <- meta(ps)[["NewID"]]
  return(ps)
}

pseq_original = readRDS("./MC_data/phyloseq/MC.trivialFilter.rds")
pseq_ND = readRDS("./MC_data/phyloseq/MC.trivialFilter.rds") %>% 
  ps_filter(!is.NCT) %>% 
  ps_mutate(Method = "ND") %>% 
  rename_samples()
pseq_restrictive = readRDS("./MC_data/phyloseq/pseq_restrictive.rds") %>% 
  ps_mutate(Method = "Restrictive") %>% 
  rename_samples()
pseq_decontam = readRDS("./MC_data/phyloseq/pseq_decontam.rds") %>% 
  ps_mutate(Method = "Decontam") %>% 
  rename_samples()
pseq_SCRuB = readRDS("./MC_data/phyloseq/pseq_SCRuB.rds") %>% 
  ps_mutate(Method = "SCRuB") %>% 
  rename_samples()
pseq_Nj = readRDS("./MC_data/phyloseq/pseq_Nejman.rds") %>% 
  ps_mutate(Method = "Nj") %>% 
  rename_samples()

pseq_all = merge_phyloseq(pseq_ND, pseq_restrictive, pseq_decontam, pseq_SCRuB, pseq_Nj)

## barplot
my_order <- c("Nj", "Restrictive", "SCRuB", "Decontam", "ND")
sample_data(pseq_all)$Method <- factor(sample_data(pseq_all)$Method, levels = my_order) 
sample_data(pseq_all)$dilution <- factor(sample_data(pseq_all)$dilution, levels = c(1,2,4,8, 16,32,64)) 
plt <- pseq_all %>% 
  comp_barplot(
    tax_level = "species", n_taxa = 20, bar_width = 0.9, sample_order = "asis", # don't sort
    x = "dilution" # x argument is available since microViz 0.9.7
  ) +
  facet_wrap(
    facets = vars(Method),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Sample ID", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + # axis starts exactly at 0
  theme_bw() + # slightly clearer axes for facets
  theme(panel.spacing.x = unit(6, "mm")) # space for non-overlapping axis numbers
ggsave("./results/MC_barplot.svg", plot = plt, width = 9, height = 6, dpi = 300)
plt
