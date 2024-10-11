#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;



//' covariance_cal
//'
//' This function computes three matrices, \eqn{\Sigma_1},\eqn{\Sigma_2}, and \eqn{V}.  \eqn{V} is the asymptotic variance, defined as \eqn{\Sigma_1^{-1} \Sigma_2 \Sigma_1^{-1}}.
//'
//' @param a A vector of regression coefficients.
//' @param b  A vector of biomarker values.
//' @param c The \code{rho01} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param d The \code{rho02} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{1-Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param e A transposed matrix where each row represents a covariate, including an intercept as the first row. The rows are ordered to correspond with the regression coefficients in \code{a}.
//' @param f The \code{i - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param g The \code{l - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @return A list containing:
//' \itemize{
//'   \item \code{sigma1}: The matrix \eqn{\Sigma_1}, used in calculating the asymptotic variance.
//'   \item \code{sigma2}: The matrix \eqn{\Sigma_2}, also used in calculating the asymptotic variance.
//'   \item \code{v}: The asymptotic variance \eqn{V}.
//' }
// [[Rcpp::export]]
Rcpp::List covariance_cal(arma::vec a, arma::vec b, arma::vec c,
                            arma::vec d, arma::mat e,
                            arma::ivec f, arma::ivec g) {
  // Define the inputs
  arma::vec beta = a;
  arma::vec M = b;
  arma::vec rho01= c;
  arma::vec rho02 = d;
  arma::mat XS = e;
  arma::ivec dati = f;  
  arma::ivec datj = g;  

  int n0 = M.size();
  int nij = rho01.size();
  int nb = beta.size();

  arma::mat sigma1(nb, nb, arma::fill::zeros);
  arma::mat sigma2(nb, nb, arma::fill::zeros);
  arma::cube f_store(n0, n0, nb, arma::fill::zeros);
  arma::vec poly(nb, arma::fill::zeros);

  // First loop
  for (int i = 0; i < nij; ++i) {
    poly = XS.col(i);
    double temp = exp(arma::dot(beta, poly));
    sigma1 += (poly * poly.t()) * (rho02[i] + rho01[i]) * temp / ((1.0 + temp) * (1.0 + temp));
    arma::vec f_storeij = poly / (1.0 + temp) * (rho01[i] - rho02[i] * temp);
    f_store.subcube(dati[i], datj[i], 0, dati[i], datj[i], nb - 1) = f_storeij;
  }

  // Second loop
  for (int i = 0; i < nij; ++i) {
    int i2 = dati[i];
    int j2 = datj[i];
    for (int jj = 0; jj < n0; ++jj) {
      bool keep = (i2 != j2) && (i2 != jj) && (j2 != jj);
      if (keep) {
        arma::vec vleft = f_store.subcube(i2, j2, 0, i2, j2, nb - 1) + f_store.subcube(j2, i2, 0, j2, i2, nb - 1);
        arma::vec vright = f_store.subcube(i2, jj, 0, i2, jj, nb - 1) + f_store.subcube(jj, i2, 0, jj, i2, nb - 1);
        sigma2 += vleft * vright.t();

        arma::vec vleft1 = f_store.subcube(j2, i2, 0, j2, i2, nb - 1) + f_store.subcube(i2, j2, 0, i2, j2, nb - 1);
        arma::vec vright1 = f_store.subcube(j2, jj, 0, j2, jj, nb - 1) + f_store.subcube(jj, j2, 0, jj, j2, nb - 1);
        sigma2 += vleft1 * vright1.t();
      }
    }
  }
  
  arma::mat sigma1_inv = arma::inv(sigma1);
  arma::mat v = sigma1_inv * sigma2 * sigma1_inv;
  return Rcpp::List::create(Named("sigma1") = sigma1, Named("sigma2") = sigma2, Named("v") = v);
}

