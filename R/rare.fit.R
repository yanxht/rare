#' Find all descendant leaves of a node in an hclust tree
#'
#' The function recursively finds all leaves that are descendants of a
#' node in an \code{hclust} tree.
#'
#' @param ind Index of the tree node. For an \code{hclust} tree
#' of \code{p} leaves, -\code{j} denotes the \code{j}th leaf and \code{k}
#' denotes the interior node formed at the \code{k}th merging in constructing
#' the tree. The range of \code{ind} is
#' \{-1, ..., -p, 1,..., p-1\} where \code{p-1} is the number of interior nodes.
#' @param merge A (\code{p-1})-by-2 matrix that encodes the order of
#' mergings in constructing the tree. \code{merge} uses the same notation for
#' nodes and mergings in an \code{hclust} object.
#' See \code{\link[stats]{hclust}} for details.
#'
#' @return Returns a sequence of indices for descendant leaves
#' in the leaf set \{1, ..., p\}. Unlike the notation used in
#' \code{ind}, we use positive integers to denote leaves here.
#'
#' @examples
#' \dontrun{
#' hc <- hclust(dist(USArrests), "ave")
#' # Descendant leaves of the 10th leaf (should be iteself)
#' find.leaves(-10, hc$merge)
#'
#' # Descendant leaves of the 10th interior node
#' find.leaves(10, hc$merge)
#'
#' # Descendant leaves of the root (should be all leaves)
#' ind_root <- nrow(hc$merge)
#' all.equal(find.leaves(ind_root, hc$merge), hc$order)
#' }
#'
#' @export
find.leaves <- function(ind, merge) {
  n_interior <- nrow(merge)
  stopifnot(ind >= -n_interior - 1 & ind <= n_interior)
  c(find_leaves(ind, merge)) # return as vector
}


#' Generate matrix A encoding ancestor-descendant relationships
#' in an hclust tree
#'
#' The function generates the binary matrix \code{A} defined in Yan
#' and Bien (2018). The matrix encodes ancestor-descendant relationships between leaves
#' and tree nodes in an \code{hclust} tree.
#'
#' @param hc An \code{hclust} object.
#'
#' @return Returns a \code{nvars}-by-\code{nnodes} binary matrix \code{A}
#' where \code{nvars} is the number of leaves (we associate covariate with leaf),
#' and \code{nnodes} is the number of tree nodes (including both leaves and interior nodes).
#' For an \code{hclust} tree, \code{nnodes} = \code{2*nvars-1}. \code{A[i,j]} is 1 if the \code{i}th leaf
#' is a descendant of the \code{j}th node in the tree, and 0 otherwise. \emph{By default, we
#' let the first \code{nvars} columns correspond to leaves and the remaining
#' \code{nvars-1} columns correspond to interior nodes.}
#' \code{A} is in sparse matrix format (inherit from class
#' \code{\link[Matrix]{sparseMatrix}} as in package \code{Matrix}).
#'
#' @seealso \code{\link{find.leaves}} for finding descendant leaves of a node.
#'
#' @examples
#' \dontrun{
#' # For a perfect binary tree of depth 2 below
#' #
#' #      3
#' #      /\
#' #    1    2
#' #   /\    /\
#' # -1 -2 -3 -4
#' #
#' # A can expressed as the following:
#' A_true <- cbind(diag(4),
#'                 as.matrix(c(1, 1, 0, 0)),
#'                 as.matrix(c(0, 0, 1, 1)),
#'                 as.matrix(c(1, 1, 1, 1)))
#' # Now use tree.matrix to generate A
#' tree0 <- list()
#' tree0$merge <- matrix(c(-1, -2, -3, -4, 1, 2),
#'                       ncol = 2, byrow = TRUE)
#' tree0$labels <- c("leaf1", "leaf2", "leaf3", "leaf4")
#' A <- tree.matrix(tree0)
#' all(A_true == as.matrix(A))
#'
#' # Another example
#' hc <- hclust(dist(USArrests), "ave")
#' A <- tree.matrix(hc)
#' }
#'
#' @references Yan, X. and Bien, J. (2018) \emph{Rare Feature Selection in High Dimensions}, \url{https://arxiv.org/abs/1803.06675}.
#'
#' @export
tree.matrix <- function(hc) {
  p <- nrow(hc$merge) + 1 # number of leaves
  n_interior <- nrow(hc$merge) # number of interior nodes
  A_i <- c(as.list(seq(p)), sapply(seq(n_interior), function(x) find.leaves(x, hc$merge)))
  A_j <- sapply(seq(length(A_i)), function(x) rep(x, len=length(A_i[[x]])))
  A <- sparseMatrix(i = unlist(A_i), j = unlist(A_j), x = rep(1, len=length(unlist(A_i))))
  A
}


#' Fit the rare feature selection model
#'
#' Fit the rare feature selection model proposed in Yan and Bien (2018):
#' \deqn{min_{\beta, \gamma} 0.5 * ||y - X\beta - \beta_01_n||_2^2 +
#' \lambda * (\alpha * ||\gamma_{-root}||_1 + (1-\alpha) * ||\beta||_1)}
#' using an alternating direction method of multipliers (ADMM) algorithm
#' described in Algorithm 1 of the same paper.
#' The regularization path is computed over a two-dimensional grid of
#' regularization parameters: \code{lambda} and \code{alpha}. Of the two,
#' \code{lambda} controls the overall amount of regularization, and \code{alpha}
#' controls the tradeoff between sparsity and fusion of \eqn{\beta} (larger \code{alpha}
#' induces more fusion in \eqn{\beta}).
#'
#' The function splits model fitting path by \code{alpha}. At each \code{alpha} value,
#' the model is fit on the entire sequence of \code{lambda} with warm start. We recommend
#' including an intercept (by setting \code{intercept=T}) unless the input data have been
#' centered.
#'
#' @param y Length-\code{nobs} response variable.
#' @param X \code{nobs}-by-\code{nvars} input matrix:
#' each row is an observation vector and each column stores a count covariate.
#' @param A \code{nvars}-by-\code{nnodes} binary matrix encoding ancestor-descendant relationships
#' between leaves and tree nodes, where \code{nnodes} is the total number of tree nodes.
#' \code{A[i,j]} is 1 if the \code{i}th leaf is a descendant of the \code{j}th
#' node in the tree, and 0 otherwise. \code{A} should be in sparse matrix format
#' (inherit from class \code{\link[Matrix]{sparseMatrix}} as in package \code{Matrix}).
#' When \code{A} is \code{NULL}, the function will learn \code{A} from \code{hc}.
#' @param Q \code{(nvars+nnodes)}-by-\code{nnodes} matrix with columns forming an orthonormal
#' basis for the null space of \eqn{[I_nvars:-A]}. When \code{Q} is \code{NULL}, the function will learn
#' \code{Q} using the singular value decomposition.
#' @param hc An \code{hclust} tree of \code{nvars} leaves where each leaf corresponds
#' to a covariate. If the tree is not an \code{hclust} object, user needs to provide the matrix \code{A} instead.
#' @param intercept Whether intercept be fitted (default = TRUE) or set to zero (FALSE).
#' @param lambda A user-supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on
#' \code{nlam} and \code{lam.min.ratio}.
#' @param alpha A user-supplied \code{alpha} sequence. If letting the program
#' compute its own \code{alpha} sequence, a length-\code{nalpha} sequence of
#' equally-spaced \code{alpha} values between 0 and 1 will be used. In practice,
#' user may want to provide a more fine \code{alpha} sequence to tune
#' the model to its best performance (e.g., \code{alpha = c(1-exp(seq(0, log(1e-2), len = nalpha - 1)), 1)}).
#' @param nlam Number of \code{lambda} values (default = 50).
#' @param lam.min.ratio Smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max} (i.e., the smallest value for which all coefficients are
#' zero). The default value is \code{1e-4}.
#' @param nalpha Number of \code{alpha} values (default = 10).
#' @param rho Penalty parameter for the quadratic penalty in the ADMM algorithm.
#' The default value is \code{1e-2}.
#' @param eps1 Convergence threshold in terms of the absolute tolerance level
#' for the ADMMM algorithm. The default value is \code{1e-6}.
#' @param eps2 Convergence threshold in terms of the relative tolerance level
#' for the ADMM algorithm. The default value is \code{1e-5}.
#' @param maxite Maximum number of passes over the data for every pair of
#' (\code{lambda}, \code{alpha}). The default value is \code{1e6}.
#'
#' @return Returns regression coefficients for \code{beta} and \code{gamma} and
#' intercept \code{beta0}. We use a \emph{matrix-nested-within-list} structure to store the coefficients: each list
#' item corresponds to an \code{alpha} value; matrix (or vector) in that list item stores
#' coefficients at various \code{lambda} values by columns (or entries).
#'
#' \item{beta0}{Length-\code{nalpha} list with each item storing
#' intercept across various \code{lambda} in a vector: \code{beta0[[j]][i]}
#' is intercept fitted at (\code{lambda[i]}, \code{alpha[j]}).
#' If \code{intercept = FALSE}, \code{beta0} is \code{NULL}.}
#' \item{beta}{Length-\code{nalpha} list with each item storing
#' \code{beta} coefficient at various \code{lambda} in columns of a \code{nvars}-by-\code{nlam} matrix:
#' \code{beta[[j]][, i]} is \code{beta} coeffcient fitted at (\code{lambda[i]}, \code{alpha[j]}).}
#' \item{gamma}{Length-\code{nalpha} list with each item storing
#' \code{gamma} coefficient at various \code{lambda} in columns of a \code{nnodes}-by-\code{nlam} matrix:
#' \code{gamma[[j]][, i]} is \code{gamma} coeffcient vector fitted at (\code{lambda[i]}, \code{alpha[j]}).
#' If \code{alpha[j] = 0}, the problem becomes the lasso on \code{beta} and is solved
#' with \code{\link[glmnet]{glmnet}} on \code{beta}, in which case \code{gamma[[j]] = NA}.}
#' \item{lambda}{Sequence of \code{lambda} values used in model fit.}
#' \item{alpha}{Sequence of \code{alpha} values used in model fit.}
#' \item{A}{Binary matrix encoding ancestor-descendant relationship between leaves and nodes in the tree.}
#' \item{Q}{Matrix with columns forming an orthonormal basis for the null space of \eqn{[I_nvars:-A]}.}
#' \item{intercept}{Whether an intercept is included in model fit.}
#'
#' @examples
#' \dontrun{
#' # See vignette for more details.
#' set.seed(100)
#' ts <- sample(1:length(data.rating), 400) # Train set indices
#' # Fit the model on train set
#' ourfit <- rarefit(y = data.rating[ts], X = data.dtm[ts, ], hc = data.hc, lam.min.ratio = 1e-6,
#'                   nlam = 20, nalpha = 10, rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' }
#'
#' @seealso \code{\link{rarefit.cv}}, \code{\link{rarefit.predict}}
#'
#' @references Yan, X. and Bien, J. (2018) \emph{Rare Feature Selection in High Dimensions}, \url{https://arxiv.org/abs/1803.06675}.
#'
#' @export
rarefit <- function(y, X, A = NULL, Q = NULL, hc, intercept = T, lambda = NULL, alpha = NULL,
                    nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
                    rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5, maxite = 1e6) {
  if (is.null(A)) {
    if (class(hc) != "hclust") stop("The tree needs to be an hclust object.")
    A <- tree.matrix(hc)
  } else {
    if (class(A) != "dgCMatrix") stop("A is not in sparse matrix format.")
  }
  n <- nrow(X)
  p <- ncol(X)
  nnodes <- ncol(A)

  X_use <- X <- as.matrix(X) # X to-be-used for model fit
  y_use <- as.vector(y) # y to-be-used for model fit
  # Center X and y and intercept is to be included
  if (intercept) {
    X_use <- X - matrix(rep(1, times = n), ncol = 1) %*% matrix(colMeans(X), nrow = 1)
    y_use <- y - mean(y)
  }

  # Generate alpha sequence
  if (is.null(alpha)) {
    alpha <- seq(0, 1, len = nalpha)
  } else {
    if (min(alpha) < 0 || max(alpha) > 1) stop("alpha range is out of [0, 1].")
    #alpha <- sort(alpha)
    nalpha <- length(alpha)
  }
  # Generate lambda sequence
  if(is.null(lambda)) {
    lambda <- max(abs(t(X_use) %*% y_use))/n * exp(seq(0, log(lam.min.ratio), len = nlam))
  } else {
    if (min(lambda) < 0) stop("lambda cannot be negative.")
    #lambda <- sort(lambda)
    nlam <- length(lambda)
  }

  # SVD of (I:-A)
  if (is.null(Q)) Q <- svdA(A)
  # SVD of X
  E <- svdX(X_use, rho)
  # Two implications up to this point: 1. Q and E will be stored in memory and be
  # passed in the solver as arguments (which can be problematic when n or p is large).
  # 2. ADMM iterates at a fixed rho. If varying rho is adopted, E should be updated as well.

  # Iterate solution paths along alpha
  beta0 <- beta <- gamma <- list()
  for (i in seq(nalpha)) {
    if (alpha[i] == 0) {
      # lasso on beta
      ret <- glmnet(X_use, y_use, family = "gaussian", lambda = lambda, standardize = F,
                    intercept = F, thresh = min(eps1, eps2), maxit = maxite)
      beta[[i]] <- as.matrix(ret$beta)
      gamma[[i]] <- NA
    } else if (alpha[i] == 1) {
      # lasso on non-root gamma
      ret <- glmnet(X_use %*% A, y_use, family = "gaussian", lambda = lambda, standardize = F,
                    intercept = F, penalty.factor = c(rep(1, nnodes-1), 0),
                    thresh = min(eps1, eps2), maxit = maxite)
      beta[[i]] <- as.matrix(A %*% ret$beta)
      gamma[[i]] <- as.matrix(ret$beta)
    } else {
      # general case when 0 < alpha < 1
      ret <- our_solver(X_use, as.matrix(y_use), Q, E, lambda, alpha[i], rho, eps1, eps2, maxite)
      beta[[i]] <- ret$beta
      gamma[[i]] <- ret$gamma
    }
    cat(sprintf("Finished model fits for alpha[%s].\n", i))
  }
  # Take care of intercept
  if (intercept) {
    beta0 <- lapply(beta, function(b) (sum(y) - c(matrix(colSums(X), nrow = 1) %*% b))/n)
  }
  list(beta0 = beta0, beta = beta, gamma = gamma, lambda = lambda, alpha = alpha, A = A, Q = Q, intercept = intercept)
}
