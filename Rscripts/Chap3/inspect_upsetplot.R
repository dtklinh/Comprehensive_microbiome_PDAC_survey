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

# df_raw <- data.frame(OTU = taxa_names(pseq_raw), Raw = 1)
# df_res <- data.frame(OTU = taxa_names(pseq_restrictive), Restrictive = 1)
# df_decontam <- data.frame(OTU = taxa_names(pseq_decontam), Decontam = 1)
# df_SCRuB <- data.frame(OTU = taxa_names(pseq_SCRuB), SCRuB = 1)
# df_Nj <- data.frame(OTU = taxa_names(pseq_Nj), Nj = 1)

df <- data.frame(OTU = taxa_names(pseq_raw), Raw = 1) %>% 
  left_join(., data.frame(OTU = taxa_names(pseq_restrictive), Restrictive = 1), by = "OTU") %>% 
  left_join(., data.frame(OTU = taxa_names(pseq_decontam), Decontam = 1), by = "OTU") %>% 
  left_join(., data.frame(OTU = taxa_names(pseq_SCRuB), SCRuB = 1), by = "OTU") %>% 
  left_join(., data.frame(OTU = taxa_names(pseq_Nj), Nj = 1), by = "OTU") %>% 
  select(-Raw) %>% 
  dplyr::mutate_if(is.numeric, coalesce, 0)

df_raw_tax <- tax_tibble(pseq_raw)
df_raw_abund <- taxa_sums(pseq_raw) %>% tibble::enframe(name = "OTU", value = "Reads")

df <- df %>% 
  left_join(., df_raw_tax, by = c("OTU" = "FeatureID")) %>% 
  left_join(., df_raw_abund, by = "OTU")

Medthods <- c("Restrictive", "Decontam", "SCRuB", "Nj", "NCT")

upset(
  df,
  Medthods,
  annotations = list(
    "reads" = (
      ggplot(mapping=aes(y = Reads, fill=order))
      + geom_bar(stat = "summary", fun = "sum", position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values= Kolors
                          #c('R'='#E41A1C', 'PG'='#377EB8', 'PG-13'='#4DAF4A', 'NC-17'='#FF7F00')
      )
      + ylab('reads')
      + theme(legend.position = "none")
    ),
    'order' = list(
      aes = aes(x=intersection, fill=order),
      geom = list(
        geom_bar(stat='count', position='fill', na.rm=TRUE),
        geom_text(
          aes(label=ifelse(order == 'Sphingomonadales', 'Sphing..', NA)),
          stat='count',
          position=position_fill(vjust = .5),
          na.rm=TRUE
        ),
        scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide='none'),
        scale_fill_manual(values = Kolors)
        
      )
    )
    # 'order'=(
    #   ggplot(mapping=aes(fill=order))
    #   + geom_bar(stat='count', position='fill')+
    #     geom_text(aes(y = order, label = ifelse(order == "Sphingomonadales", "Sphingomonadales", NA)))
    #   + scale_y_continuous(labels=scales::percent_format())
    #   + scale_fill_manual(values= Kolors
    #                         #c('R'='#E41A1C', 'PG'='#377EB8', 'PG-13'='#4DAF4A', 'NC-17'='#FF7F00')
    #                       )
    #   + ylab('order')
    #   + theme(legend.position = "none")
    # ),
    
  ),
  width_ratio=0.1
)

###-------------------
## integrate NCT
pseq_NCT <- readRDS("./data/Chap1/NCTs_v7_wrench.rds")
lst_NCT <- get_most_prev(pseq_NCT, 0.20, 1000)

df <- df %>% 
  left_join(., data.frame(OTU = lst_NCT$taxa, NCT = 1), by = "OTU") %>% 
  dplyr::mutate_if(is.numeric, coalesce, 0)
