## playground of upset plot
library(ggplot2)
library(ComplexUpset)

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
