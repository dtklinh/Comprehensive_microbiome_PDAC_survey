## Metadata preparation of metadata
library(phyloseq)
library(microViz)
library(openxlsx)
library(dplyr)

##------------------------------------

##load LP anonymous
df_LP_ananymous <- openxlsx::read.xlsx("./Rscripts/ForPublication/Qiita_metadata/LP_anonymous.xlsx")

Chap1_nct <- openxlsx::read.xlsx("./Rscripts/ForPublication/Qiita_metadata/Chap1.xlsx", detectDates = TRUE) %>% 
  dplyr::select(-Barcode) %>% 
  rename(Batch_run = Run) %>% 
  left_join(., y = df_LP_ananymous, by = "LP") %>% 
  dplyr::select(-LP) %>% 
  rename(LP = LP_anonymous)

Chap2 <- openxlsx::read.xlsx("./Rscripts/ForPublication/Qiita_metadata/Chap2_J_U.xlsx", detectDates = TRUE)
Chap2$ID <- sub("\\.fastq\\.gz$", "", Chap2$ID)
Chap2$Sequencing_date <- c(rep("2024-06-04", 8), rep("2024-06-06", 8)) %>% as.Date(format = "%Y-%m-%d")
Chap2 <- Chap2 %>% 
  dplyr::select(-Barcode)

KPC_FF_R1 <- openxlsx::read.xlsx("./Rscripts/ForPublication/Qiita_metadata/KPC_FF_R1.xlsx", detectDates = TRUE) %>% 
  mutate(Sequencing_date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>% 
  left_join(., y = df_LP_ananymous, by = c("person" = "LP")) %>% 
  dplyr::select(-c(species, ffpe.bulk, barcode,person, year, month, day)) %>% 
  dplyr::rename(true_control = true.control, LP = LP_anonymous) %>% 
  mutate(FFPE_FF = "Fresh_frozen", Replicate = 1)

KPC_FF_R2 <- openxlsx::read.xlsx("./Rscripts/ForPublication/Qiita_metadata/KPC_FF_R2.xlsx", detectDates = TRUE) %>% 
  mutate(Sequencing_date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>% 
  left_join(., y = df_LP_ananymous, by = c("person" = "LP")) %>% 
  dplyr::select(-c(species, ffpe.bulk, barcode, pseudonym, person, is.tumor, year, month, day)) %>% 
  dplyr::rename(true_control = true.control, LP = LP_anonymous) %>% 
  mutate(FFPE_FF = "Fresh_frozen", Replicate = 2)


KPC_FFPE <- openxlsx::read.xlsx("./Rscripts/ForPublication/Qiita_metadata/KPC_FFPE.xlsx", detectDates = TRUE) %>% 
  mutate(Sequencing_date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>% 
  left_join(., y = df_LP_ananymous, by = c("person" = "LP")) %>% 
  dplyr::select(-c(species, barcode, pseudonym, person, is.tumor, year, month, day)) %>% 
  dplyr::rename(FFPE_FF = ffpe.bulk, true_control = true.control, LP = LP_anonymous) %>% 
  mutate(Replicate = NA) %>% 
  dplyr::select(ID, AN_NR, true_control, NCT_type, Sequencing_date, LP, FFPE_FF, Replicate)

## Concatenate
Intratumor_study <- rbind(KPC_FF_R1, KPC_FF_R2, KPC_FFPE)

## Concat Chap1 and Chap2
Chap1_2 <- Chap1_nct %>% 
  left_join(., y = Chap2[, c("ID", "Environment")], by = "ID")

# save to file
openxlsx::write.xlsx(Chap1_2, "./Rscripts/ForPublication/Qiita_metadata/Only_involved_NCT.xlsx")
openxlsx::write.xlsx(Intratumor_study, "./Rscripts/ForPublication/Qiita_metadata/IntratumoralStudy.xlsx")

## merge everything
colnames()