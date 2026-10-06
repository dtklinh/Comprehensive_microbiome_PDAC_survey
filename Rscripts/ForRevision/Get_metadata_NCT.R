library(phyloseq)
library(microbiome)
library(microViz)
library(openxlsx)

## remove old object
rm(list = setdiff(ls(), c(lsf.str(), "Kolors", "Kolors_Chap2")))
## load objects

LP_anonymus = read.xlsx("../ForPublication/Qiita_metadata/LP_anonymous.xlsx")

pseq <- rprojroot::find_rstudio_root_file() %>% 
  file.path("data/Chap1/NCT_2026April_v0.rds") %>% 
  readRDS()
N <- 1 # 500 ## accept all NCT samples.
pseq <- pseq %>% 
  #ps_filter(Control_type %in% c("Buffer", "Paraffin", "PCR", "Sequencing")) %>%
  ps_filter(Control_type %in% c("Buffer", "Paraffin", "PCR")) %>% 
  prune_samples(sample_sums(.) >=N, .) %>% 
  ps_get() %>% 
  microViz::tax_filter(min_prevalence = 1,
                       prev_detection_threshold = 2,
                       min_total_abundance = 1e-6)
## add year and quarter to metadata
lib_depth = sample_sums(pseq)
pseq <- pseq %>% 
  ps_mutate(Year = as.character(lubridate::year(Seq_date))) %>% 
  ps_mutate(Quarter = as.character(lubridate::quarter(Seq_date))) %>% 
  ps_mutate(Lib_depth = lib_depth) %>% 
  ps_select(Run, Barcode, LP = Person, Control_type, Seq_date, Seq_type, Host, Year, Quarter, Lib_depth) %>% 
  ps_arrange(Seq_date)

tmp2 <- pseq %>% 
  estimate_richness() %>% 
  dplyr::select(Observed, Shannon)
##------------------------------
new_names <- rownames(tmp2) %>% 
  gsub("^X", "", .) %>% 
  gsub("\\.", "-", .)
rownames(tmp2) <- new_names
rich_meta <- merge(pseq %>% sample_data(), tmp2, by = "row.names")

rich_meta = dplyr::left_join(x = rich_meta, y = LP_anonymus, by = "LP")

rich_meta %>% 
  dplyr::select(Control_type, LP = LP_anonymous, Seq_date, Year, Quarter, ReadCount = Lib_depth, Observed, Shannon) %>% 
  write.xlsx("NCT_rich_meta.xlsx")
