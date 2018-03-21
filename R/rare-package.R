#' Model path for tree-based lasso framework for selecting rare features
#'
#' The package fits the linear model with tree-based lasso regularization proposed in
#' Yan and Bien (2018) using alternating direction method of multipliers
#' (ADMM). The ADMM algorithm is proposed in Algorithm 1 of the same paper.
#' The package also provides tools for tuning regularization
#' parameters, making predictions from the fitted model and visualizing recovered
#' groups of the covariates in a dendrogram.
#'
#' Its main functions are \code{\link{rarefit}}, \code{\link{rarefit.cv}},
#' \code{\link{rarefit.predict}}, \code{\link{group.recover}} and
#' \code{\link{group.plot}}.
#'
#' @author Xiaohan Yan \email{xy257@@cornell.edu}, Jacob Bien
#' @references Yan, X. and Bien, J. (2018) \emph{Rare Feature Selection in High Dimensions}, \url{https://arxiv.org/abs/1803.06675}.
#' @name rare-package
#' @docType package
#' @useDynLib rare
#' @importFrom Rcpp sourceCpp
#' @importFrom glmnet glmnet
#' @import Matrix
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot
#' @importFrom stats as.dendrogram
NULL
