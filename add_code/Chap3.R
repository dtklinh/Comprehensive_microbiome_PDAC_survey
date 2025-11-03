### Addiontally code for chap 3

get_most_prev <- function(pseq, prev_threshold = 0.5, N = 25){
  res <- pseq %>% 
    prevalence() %>% 
    tibble(taxa = names(.), value = .) %>% 
    filter(value >= prev_threshold) %>% 
    arrange(desc(value)) %>% 
    head(N)
  return(res)
}
