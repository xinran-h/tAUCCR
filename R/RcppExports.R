# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' covariance_cal
#'
#' This function computes three matrices, \eqn{\Sigma_1},\eqn{\Sigma_2}, and \eqn{V}.  \eqn{V} is the asymptotic variance, defined as \eqn{\Sigma_1^{-1} \Sigma_2 \Sigma_1^{-1}}.
#'
#' @param a A vector of regression coefficients.
#' @param b  A vector of biomarker values.
#' @param c The \code{rho01} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param d The \code{rho02} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{1-Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param e A transposed matrix where each row represents a covariate, including an intercept as the first row. The rows are ordered to correspond with the regression coefficients in \code{a}.
#' @param f The \code{i - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param g The \code{l - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @return A list containing:
#' \itemize{
#'   \item \code{sigma1}: The matrix \eqn{\Sigma_1}, used in calculating the asymptotic variance.
#'   \item \code{sigma2}: The matrix \eqn{\Sigma_2}, also used in calculating the asymptotic variance.
#'   \item \code{v}: The asymptotic variance \eqn{V}.
#' }
covariance_cal <- function(a, b, c, d, e, f, g) {
    .Call(`_tAUCCR_covariance_cal`, a, b, c, d, e, f, g)
}

