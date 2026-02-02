## playground of upset plot
## source: https://krassowski.github.io/complex-upset/articles/Examples_R.html
library(ggplot2)
library(phyloseq)
library(microViz)
library(microbiome)
library(ComplexUpset)
library(tidyverse)

if(!require(ggplot2movies)) install.packages('ggplot2movies')
movies = ggplot2movies::movies
genres = c('Action', 'Animation', 'Comedy', 'Drama', 'Documentary', 'Romance')

upset(
  movies,
  genres,
  annotations = list(
    'Length'=ggplot(mapping=aes(x=intersection, y=length)) + geom_boxplot(),
    'Rating'=ggplot(mapping=aes(x=intersection, y=rating))
    # if you do not want to install ggbeeswarm, you can use geom_jitter
    + ggbeeswarm::geom_quasirandom(aes(color=log10(votes)))
    + geom_violin(width=1.1, alpha=0.5)
  ),
  queries=list(
    upset_query(
      intersect=c('Drama', 'Comedy'),
      color='red',
      fill='red',
      only_components=c('intersections_matrix', 'Intersection size')
    ),
    upset_query(
      set='Drama',
      fill='blue'
    ),
    upset_query(
      intersect=c('Romance', 'Drama'),
      fill='yellow',
      only_components=c('Length')
    )
  ),
  min_size=10,
  width_ratio=0.1
)

set_size(8, 5)

upset(
  movies,
  genres,
  annotations = list(
    'MPAA Rating'=(
      ggplot(mapping=aes(fill=mpaa))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values=c(
        'R'='#E41A1C', 'PG'='#377EB8',
        'PG-13'='#4DAF4A', 'NC-17'='#FF7F00'
      ))
      + ylab('MPAA Rating')
    )
  ),
  width_ratio=0.1
)


#### with phyloseq###############################
source("./add_code/Functions.R")
df_additionalInfo <- readxl::read_xlsx("./meta/Mice_meta.xlsx")
pseq_raw <- readRDS("./data/Chap3/pseq_Proj5_postFilter_v04.rds") %>% 
  ps_filter(ffpe.bulk == "bulk") %>% 
  ps_filter(true.control == "true") %>% 
  append_AN_NR(df_additional = df_additionalInfo) #%>% 
  #WrenchWrapper(grp = "Sex") 
pseq_restrictive <- readRDS("./data/Chap3/pseq_bulk_restrictive.rds") %>% 
  append_AN_NR(df_additionalInfo) #%>% 
  #WrenchWrapper(grp = "Sex")
pseq_decontam <- readRDS("./data/Chap3/pseq_bulk_decontam_p0.5.rds") %>% 
  append_AN_NR(df_additionalInfo) #%>% 
  #WrenchWrapper(grp = "Sex")
pseq_SCRuB <- readRDS("./data/Chap3/pseq_bulk_SCRuB.rds") %>% 
  append_AN_NR(df_additionalInfo) #%>% 
  #WrenchWrapper(grp = "Sex")
pseq_Nj <- readRDS("./data/Chap3/pseq_bulk_Fisher_v02.rds") %>% 
  append_AN_NR(df_additionalInfo) #%>% 
  #WrenchWrapper(grp = "Sex")

df_raw <- data.frame(OTU = taxa_names(pseq_raw), Raw = 1)
df_res <- data.frame(OTU = taxa_names(pseq_restrictive), Restrictive = 1)
df_decontam <- data.frame(OTU = taxa_names(pseq_decontam), Decontam = 1)
df_SCRuB <- data.frame(OTU = taxa_names(pseq_SCRuB), SCRuB = 1)
df_Nj <- data.frame(OTU = taxa_names(pseq_Nj), Nj = 1)

df <- df_raw %>% 
  left_join(., df_res, by = "OTU") %>% 
  left_join(., df_decontam, by = "OTU") %>% 
  left_join(., df_SCRuB, by = "OTU") %>% 
  left_join(., df_Nj, by = "OTU") %>% 
  select(-Raw) %>% 
  dplyr::mutate_if(is.numeric, coalesce, 0)

df_raw_tax <- tax_tibble(pseq_raw)
df <- df %>% 
  left_join(., df_raw_tax, by = c("OTU" = "FeatureID"))

Medthods <- c("Restrictive", "Decontam", "SCRuB", "Nj")

upset(
  df,
  Medthods,
  annotations = list(
    'order'=(
      ggplot(mapping=aes(fill=order))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values= Kolors
                            #c('R'='#E41A1C', 'PG'='#377EB8', 'PG-13'='#4DAF4A', 'NC-17'='#FF7F00')
                          )
      + ylab('order')
    )
  ),
  width_ratio=0.1
)

### add abundance info
df_raw_abund <- taxa_sums(pseq_raw) %>% 
  tibble(OTU = names(.), Count = .)

df <- df %>% 
  left_join(., y = df_raw_abund, by = "OTU")
