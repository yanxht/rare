#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]] 
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector threshold(const NumericVector& b1, const NumericVector& b2, const NumericVector& b3, const NumericVector& b,
                        const NumericVector& r1, const NumericVector& r2, const NumericVector& r, const NumericVector& v1,
                        const NumericVector& v2, const NumericVector& v3, const NumericVector& u1, const NumericVector& u2, double eps1, double eps2) {
  // Compute threshold values for primal residual and dual residual in ADMM.
  int p = b.size(), t_size = b.size();
  NumericVector thres(3);
  thres[0] = std::sqrt((double) 3.0 * p + 2.0 * t_size) * eps1; // absolute error bound
  thres[1] = std::max(std::sqrt(std::inner_product(b1.begin(), b1.end(), b1.begin(), 0.0)
                                + std::inner_product(b2.begin(), b2.end(), b2.end(), 0.0)
                                + std::inner_product(b3.begin(), b3.end(), b3.begin(), 0.0)
                                + std::inner_product(r1.begin(), r1.end(), r1.begin(), 0.0)
                                + std::inner_product(r2.begin(), r2.end(), r2.begin(), 0.0)),
                      std::sqrt(3.0 * std::inner_product(b.begin(), b.end(), b.begin(), 0.0)
                                + 2.0 * std::inner_product(r.begin(), r.end(), r.begin(), 0.0))) * eps2;
  thres[2] = std::sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), 0.0)
                       + std::inner_product(v2.begin(), v2.end(), v2.begin(), 0.0)
                       + std::inner_product(v3.begin(), v3.end(), v3.begin(), 0.0)
                       + std::inner_product(u1.begin(), u1.end(), u1.begin(), 0.0)
                       + std::inner_product(u2.begin(), u2.end(), u2.begin(), 0.0)) * eps2;
  return thres;
}

// [[Rcpp::export]]
double pri_Resid(const NumericVector& b1, const NumericVector& b2, const NumericVector& b3, const NumericVector& b,
                 const NumericVector& r1, const NumericVector& r2, const NumericVector& r) {
  // Compute primal residual in ADMM.
  NumericVector diff_b1 = b1 - b;
  NumericVector diff_b2 = b2 - b;
  NumericVector diff_b3 = b3 - b;
  NumericVector diff_r1 = r1 - r;
  NumericVector diff_r2 = r2 - r;
  double ret = std::sqrt(std::inner_product(diff_b1.begin(), diff_b1.end(), diff_b1.begin(), 0.0)
                         + std::inner_product(diff_b2.begin(), diff_b2.end(), diff_b2.begin(), 0.0)
                         + std::inner_product(diff_b3.begin(), diff_b3.end(), diff_b3.begin(), 0.0)
                         + std::inner_product(diff_r1.begin(), diff_r1.end(), diff_r1.begin(), 0.0)
                         + std::inner_product(diff_r2.begin(), diff_r2.end(), diff_r2.begin(), 0.0));
  return ret;
}

// [[Rcpp::export]]
double dual_Resid(const NumericVector& b0, const NumericVector& b, const NumericVector& r0, const NumericVector& r, double rho) {
  NumericVector diff_b = b - b0;
  NumericVector diff_r = r - r0;
  double ret = std::sqrt(3.0 * std::inner_product(diff_b.begin(), diff_b.end(), diff_b.begin(), 0.0)
                         + 2.0 * std::inner_product(diff_r.begin(), diff_r.end(), diff_r.begin(), 0.0)) * rho;
  return ret;
}

// [[Rcpp::export]]
double objval(const arma::mat& X, const arma::mat& y, const NumericVector& b, const NumericVector& r, double lam, double alpha, int n, int t_size) {
  arma::mat b_vec(as<arma::vec>(b));
  arma::mat r_vec(as<arma::vec>(r));
  arma::mat resid = y - X * b_vec;
  double ret = dot(resid, resid)/(2.0 * n) + lam * alpha * accu(abs(r_vec.rows(0, t_size - 2))) + lam * (1.0 - alpha) * accu(abs(b_vec));
  return ret;
}

// [[Rcpp::export]]
List our_solver(const arma::mat& X, const arma::mat& y, const arma::mat& Q, const arma::mat& E, 
                const NumericVector& lam, double alpha, double rho, double eps1, double eps2, int maxite) {
  // Solve our optimization problem at the given alpha value and for a sequence of lambda values
  arma::mat Xty0(X.t() * y);
  NumericVector Xty = wrap(Xty0);
  // Iterate all lambda at the given alpha
  int n = X.n_rows;
  int p = X.n_cols;
  int t_size = Q.n_cols;
  int nlam = lam.size();
  IntegerVector ites(nlam);
  NumericMatrix bma(p, nlam), rma(t_size, nlam);
  NumericVector b(p), v1(p), v2(p), v3(p), r(t_size), u1(t_size), u2(t_size);
  for (int k = 0; k < nlam; ++k) {
    int ite = 0;
    double pri_resid, dual_resid;
    NumericVector thres, b1, b2, b3(p), b_mean0, r1, r2(t_size), r_mean0, temp1, temp2(p + t_size), b3r2;
    do {
      checkUserInterrupt();
      temp1 = Xty + n * rho * b - n * v1;
      b1 = wrap(E * as<arma::vec>(temp1));
      b2 = (b - v2/rho) * pmax((1.0 - lam[k] * (1.0 - alpha)/abs(rho * b - v2)), 0.0);
      r1 = (r - u1/rho) * pmax((1.0 - lam[k] * alpha/abs(rho * r - u1)), 0.0);
      r1[t_size - 1] = r[t_size - 1] - u1[t_size - 1]/rho;
      for (int i = 0; i < (p + t_size); ++i) {
        if (i < p) {
          temp2[i] = b[i] - v3[i]/rho;
        } else {
          temp2[i] = r[i - p] - u2[i - p]/rho;
        }
      }
      b3r2 = wrap(Q * (Q.t() * as<arma::vec>(temp2)));
      for (int j = 0; j < (p + t_size); ++j) {
        if (j < p) {
          b3[j] = b3r2[j];
        } else {
          r2[j - p] = b3r2[j];
        }
      }
      b_mean0 = clone(b);
      b = (b1 + b2 + b3)/3.0;
      v1 += rho * (b1 - b);
      v2 += rho * (b2 - b);
      v3 += rho * (b3 - b);
      r_mean0 = clone(r);
      r = (r1 + r2)/2.0;
      u1 += rho * (r1 - r);
      u2 += rho * (r2 - r);
      // Update ite, threshold values and residual values for comparison
      ite += 1;
      thres = threshold(b1, b2, b3, b, r1, r2, r, v1, v2, v3, u1, u2, eps1, eps2);
      pri_resid = pri_Resid(b1, b2, b3, b, r1, r2, r);
      dual_resid = dual_Resid(b_mean0, b, r_mean0, r, rho);
      //if (ite % 100 == 0) {
      //  Rcout << "Iteration " << ite << "at lam[" << k << "]) with objval = " << objval(X, y, b, r, lam[k], alpha, n, t_size)
      //        << ", pri_resid = " << pri_resid << " and dual_resid = " << dual_resid << "." << std::endl;
      //}
    } while (ite < maxite && !(pri_resid <= thres[0] + thres[1] && dual_resid <= thres[0] + thres[2]) && !(pri_resid <= 1.0e-7 && dual_resid <= 1.0e-7));
    bma(_, k) = b2; // returning b2 for exact sparsity in beta
    rma(_, k) = r1; // returning r1 for exact sparsity in gamma
    ites(k) = ite;
  }
  List ret;
  ret["beta"] = bma;
  ret["gamma"] = rma;
  ret["ites"] = ites;
  return ret;
}

// [[Rcpp::export]]
arma::mat svdA (const arma::sp_mat& A) {
  // Compute the null space of (I:-A)
  int p = A.n_rows;
  arma::mat M(join_rows(arma::speye(p, p), -A));
  arma::mat Q = null(M);
  return Q;
}

// [[Rcpp::export]]
arma::mat svdX (const arma::mat& X, double rho) {
  // Compute matrix E in updating beta 1
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd_econ(U, s, V, X);
  arma::mat E;
  if (n < p) {
    E = V * diagmat(pow(square(s) + 1.0 * n * rho, -1)) * V.t() + (arma::eye(p, p) - V * V.t()) / (1.0 * n * rho);
  } else {
    E = V * diagmat(pow(square(s) + 1.0 * n * rho, -1)) * V.t();
  }
  return E;
}

// [[Rcpp::export]]
arma::Col<int> find_leaves (int ind, const arma::Mat<int>& merge) {
  arma::Col<int> v;
  if (ind < 0) {
    // leaf
    v.set_size(1);
    v[0] = -ind;
  } else {
    // interior node; recursively find its descendant leaves
    arma::Col<int> left = find_leaves(merge(ind-1, 0), merge);
    arma::Col<int> right = find_leaves(merge(ind-1, 1), merge);
    v = join_cols(left, right);
  }
  return v;
}
