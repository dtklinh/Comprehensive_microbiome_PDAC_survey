## playground of upset plot
## source: https://krassowski.github.io/complex-upset/articles/Examples_R.html
library(ggplot2)
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
