#' TripAdvisor hotel review ratings
#'
#' A length-500 TripAdvisor review ratings on the scale 1 to 5.
#'
#' @source TripAdvisor Data Set: \url{http://times.cs.uiuc.edu/~wang296/Data/}
"data.rating"


#' Document-term matrix for adjectives in TripAdvisor hotel reviews
#'
#' A 500-by-200 document-term matrix for 200 adjectives appearing in 500 TripAdvisor reviews.
#' The document-term matrix is in sparse format.
#' @seealso \code{\link{data.rating}}, \code{\link{data.hc}}.
"data.dtm"


#' Hierarchical clustering tree for adjectives in TripAdvisor data set
#'
#' An \code{hclust} tree for the 200 adjectives appearing in the TripAdvisor reviews.
#' The tree was generated with 100-dimensional word embeddings pre-trained by GloVe
#' (Pennington et al., 2014) on Gigaword5 and Wikipedia2014 corpora for the adjectives.
#'
#' @source Embeddings available at \url{http://nlp.stanford.edu/data/glove.6B.zip}
#'
#' @references Pennington, J., Socher, R., and Manning, C. D. (2014).
#' Glove: Global vectors for word representation.
#' \emph{In Empirical Methods in Natural Language Processing (EMNLP)}, pages 1532–1543.
"data.hc"
