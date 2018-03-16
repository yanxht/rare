#' Recover aggregated groups of leaf indices
#'
#' The function finds aggregated groups of leaf indices by traversing non-zero
#' \code{gamma} elements and finding descendant leaves at each \code{gamma} element. In our problem,
#' \code{gamma} are latent variables corresponding to tree nodes. The order
#' of the traversal is post-order, i.e., a node is visited after its descendants.
#'
#' @param gamma Length-\code{nnodes} latent variable coefficients. Note that \code{\link{rarefit}}
#' returns \code{NA} as \code{gamma} value when \code{alpha} is zero,
#' in which case our problem becomes the lasso on \code{beta}.
#' @param A \code{nvars}-by-\code{nnodes} binary matrix encoding ancestor-descendant relationships
#' between leaves and nodes in the tree.
#' @param postorder Length-\code{nnodes} integer vector encoding post-order traversal of the tree
#' nodes such that \code{seq(nnodes)[postorder]} ensures a node appear after its descendants.
#' Default is \code{seq(nnodes)}, which gives post-order when \code{A} is generated using \code{\link{tree.matrix}}
#' for an \code{hclust} tree.
#' @return Returns a list of recovered groups of leaf indices.
#'
#' @examples
#' \dontrun{
#' # See vignette for more details.
#' set.seed(100)
#' ts <- sample(1:length(data.rating), 400) # Train set indices
#' # Fit the model on train set
#' ourfit <- rarefit(y = data.rating[ts], X = data.dtm[ts, ], hc = data.hc, lam.min.ratio = 1e-6,
#'                   nlam = 20, nalpha = 10, rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Cross validation
#' ourfit.cv <- rarefit.cv(ourfit, y = data.rating[ts], X = data.dtm[ts, ],
#'                         rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Group recovered at optimal beta and gamma
#' ibest.lambda <- ourfit.cv$ibest[1]
#' ibest.alpha <- ourfit.cv$ibest[2]
#' gamma.opt <- ourfit$gamma[[ibest.alpha]][, ibest.lambda] # works if ibest.alpha > 1
#' groups.opt <- group.recover(gamma.opt, ourfit$A)
#' }
#'
#' @export
group.recover <- function(gamma, A, postorder = seq(ncol(A))) {
  if (any(is.na(gamma))) stop("Missing gamma value. It is likely due to the input
                              gamma corresponding to alpha = 0. See ?rarefit for details.")
  nz <- postorder[postorder %in% which(gamma != 0)]
  groups <- list()
  for (i in seq(length(nz))) {
    groups[[i]] <- setdiff(which(A[, nz[i]] == 1), unlist(groups))
  }
  # Remove all empty elements in groups
  groups[unlist(lapply(groups, function(x) length(x) > 0))]
}


#' Visualize groups by coloring branches and leaves of an hclust tree
#'
#' The function plots an \code{hclust} tree with branches and leaves colored
#' based on group membership. The groups span the covariate indices \{1, ..., \code{nvars}\}.
#' Covariates from the same group share equal coefficient (\code{beta}), and sibling
#' groups have different coefficients. The function determines groups based on
#' the sparsity in \code{gamma}. In an \code{hclust} tree with \code{beta[i]} on the
#' \code{i}th leaf, the branch and leaf are colored in blue, red or gray according to \code{beta[i]}
#' being positive, negative or zero, respectively. The larger the magnitude of \code{beta[i]} is,
#' the darker the color will be. So branches and leaves from the same group will have the
#' same color.
#'
#' @param beta Length-\code{nvars} vector of covariate coefficient.
#' @param gamma Length-\code{nnodes} vector of latent variable coefficient. Note that \code{\link{rarefit}}
#' returns \code{NA} as \code{gamma} value when \code{alpha} is zero,
#' in which case our problem becomes the lasso on \code{beta}.
#' @param A \code{nvars}-by-\code{nnodes} binary matrix encoding ancestor-descendant relationships
#' between leaves and nodes in the tree.
#' @param hc An \code{hclust} tree of \code{nvars} leaves where each leaf corresponds to a covariate.
#' @param nbreaks Number of breaks in binning \code{beta} elements (positive part and negative part
#' are done separately). Each bin is associated with a color based on the magnitude and
#' positivity/negativity of \code{beta} elements in the bin.
#'
#' @examples
#' \dontrun{
#' # See vignette for more details.
#' set.seed(100)
#' ts <- sample(1:length(data.rating), 400) # Train set indices
#' # Fit the model on train set
#' ourfit <- rarefit(y = data.rating[ts], X = data.dtm[ts, ], hc = data.hc, lam.min.ratio = 1e-6,
#'                   nlam = 20, nalpha = 10, rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Cross validation
#' ourfit.cv <- rarefit.cv(ourfit, y = data.rating[ts], X = data.dtm[ts, ],
#'                         rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Visualize the groups at optimal beta and gamma
#' ibest.lambda <- ourfit.cv$ibest[1]
#' ibest.alpha <- ourfit.cv$ibest[2]
#' beta.opt <- ourfit$beta[[ibest.alpha]][, ibest.lambda]
#' gamma.opt <- ourfit$gamma[[ibest.alpha]][, ibest.lambda] # works if ibest.alpha > 1
#' # Visualize the groups at optimal beta and gamma
#' group.plot(beta.opt, gamma.opt, ourfit$A, data.hc)
#' }
#'
#' @export
group.plot <- function(beta, gamma, A, hc, nbreaks = 20) {
  group <- group.recover(gamma, A) # recover groups from gamma
  jj <- sapply(seq(length(group)), function(x) rep(x, len=length(group[[x]])))
  # the ith row gives membership of the ith leaf in group: membership[i, ] %*% seq(length(group))
  membership <- sparseMatrix(i = unlist(group), j = unlist(jj), x = rep(1, len=length(unlist(group))))
  groupvec <- as.vector(membership %*% seq(length(group)))

  # convert hc to dendrogram
  hcd = as.dendrogram(hc)
  beta_grp <- beta[unlist(lapply(group, function(x) x[1]))] # clustered beta_grp values
  beta_grp_neg <- beta_grp[beta_grp < 0]
  beta_grp_pos <- beta_grp[beta_grp > 0]
  # make breaks for color
  breaks_neg <- seq(min(beta_grp_neg), max(beta_grp_neg), len = nbreaks)
  breaks_pos <- seq(min(beta_grp_pos), max(beta_grp_pos), len = nbreaks)
  # make color palettes for negative and positive parts respectively
  # (blue-gray95-red = positive-0-negative)
  colfunc1 <- colorRampPalette(c("lightgray", "blue"))
  colfunc2 <- colorRampPalette(c("red", "lightgray"))
  col0 <- "gray95"
  colpos <- colfunc1(nbreaks + 1)[-1]
  colneg <- colfunc2(nbreaks + 1)[-21]
  # color labels and branches
  col_cluster <- rep(NA, len = length(beta_grp))
  col_cluster[beta_grp == 0] <- col0
  col_cluster[beta_grp < 0 & beta_grp != 0] <- colneg[findInterval(beta_grp_neg, breaks_neg)]
  col_cluster[beta_grp > 0 & beta_grp != 0] <- colpos[findInterval(beta_grp_pos, breaks_pos)]
  dend_col <- dendextend::color_labels(hcd, labels = hc$labels[hc$order], col = col_cluster[groupvec][hc$order])
  dend_col <- dendextend::color_branches(dend_col, col = col_cluster[unique(groupvec[hc$order])], clusters = groupvec[hc$order])
  plot(dend_col)
}
