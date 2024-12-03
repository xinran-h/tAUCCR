####  All functions ####


#' data_crossingdc
#' 
#' This is to reshape data for each cause of failure, the result reshaped data is ready for parameter estimation.
#' @import dplyr
#' @import rlang
#' @import tidyr
#' @param n   Total sample size
#' @param vt  Vector of observed time
#' @param vc Vector of censoring indicator
#' @param vm vector of biomarker values
#' @param vxs Covariate matrix
#' @param eta Censoring indicator for each event type
#' @param c0 A constant value used for the bandwidth h, h = c0*n^(-1/3)
#' @param cont_vars A vector of continuous variables in vxs.
#' @param discrete_vars A vector of discrete variables in vxs.
#' @returns A dataframe where each row corresponds to a case-control pair at each observed event time.
#' The dataframe has the following components:
#' \itemize{
#'   \item \code{i}: Index for cases.
#'   \item \code{yi}: Observed event time.
#'   \item \code{Covariates with "i" suffix}: Covariate columns, with names based on the original covariate name and an "i" suffix for case values.
#'   \item \code{l}: Index for controls.
#'   \item \code{Iil}: Indicator for comparison of baseline biomarker values (1 if case biomarker > control biomarker, 0 otherwise).
#'   \item \code{dif}: Euclidean norm among continuous covariates, when available.
#'   \item \code{rho01} and \code{rho02}: Numeric vectors representing concordant and discordant event values, respectively, when continuous covariates are present.
#'   \item \code{rhoweight}: Kernel weight vector, calculated with continuous covariates.
#' }
#' @examples
#' \dontrun{
#' data_crossingdc(n, Y0,C0,M0,X,1,2, c("age","bmi"),c("gender","DRUG"))
#'
#' } 
#' @export
data_crossingdc = function(n, vt,vc,vm,vxs,eta,c0, cont_vars,discrete_vars)
{
  nx<-ncol(vxs) 
  ni<-length(vc) 
  index<-which(vc == eta)
  dat<-cbind(l=1:ni,vt,vm,vc,vxs) 
  dati<-dat[index,]%>%dplyr::rename(i=l)
  
  names(dati)[1:4]<-c("i","yi","mi","event_indicator_i")
  names(dati)[5:(nx+4)]<-paste0(colnames(vxs),"i") 
  names(dat)[1:4]<-c("l","yl","ml","event_indicator_l") 
  names(dat)[5:(nx+4)]<-paste0(colnames(vxs),"l")
  
  
  if (!is.null(cont_vars)) {
    h <- c0 * n^(-1 / 3)
    
    squared_diffs <- vector("list", length(cont_vars))
    names(squared_diffs) <- cont_vars
    for (var in cont_vars) {
      case_var <- paste0(var, "i")
      control_var <- paste0(var, "l")
      squared_diffs[[var]] <- paste0("(", case_var, " - ", control_var, ")^2")
    }
    sum_of_squares_expr <- paste(squared_diffs, collapse = " + ")
    final_expr <- paste0("sqrt(", sum_of_squares_expr, ")")
  }
 
  if (!is.null(discrete_vars)) {
    conditions <- lapply(discrete_vars, function(var) {
      sym_case <- rlang::sym(paste0(var, "i"))
      sym_control <- rlang::sym(paste0(var, "l"))
      rlang::expr(!!sym_case == !!sym_control)
    })
  }
  
  
  
  tableij <- tidyr::crossing(dati, dat, .name_repair = "universal") %>%
    dplyr::filter(yi < yl) %>%
    dplyr::mutate(Iil = as.numeric(mi > ml))
  if (!is.null(discrete_vars)) {
    tableij <- tableij %>% dplyr::filter(!!!conditions)
  }
  tableij <- tableij %>% dplyr::group_by(i) %>%
    dplyr::arrange(i, ml)
  
  if (!is.null(cont_vars)) {
    tableij <- tableij %>%
      dplyr::mutate(
        dif = eval(parse(text = final_expr)),
        kh = 0.75 * (1 - (dif / h)^2) / h * as.numeric(abs(dif / h) < 1),
        rho01 = Iil * kh,
        rho02 = (1 - Iil) * kh,
        rhoweight = rho01 + rho02
      )
  }

  select_columns <- c("i", "l", "yi", "Iil")
  if (!is.null(cont_vars)) {
    select_columns <- c(select_columns, paste0(cont_vars, "i"),"rho01", "rho02", "rhoweight", "dif")
  }
  if (!is.null(discrete_vars)) {
    select_columns <- c(select_columns, paste0(discrete_vars, "i"))
  }
  
  tableij <- tableij %>%
    dplyr::select(dplyr::any_of(select_columns))
  
  
  return(data.frame(tableij))
  
}




#' auc_pred
#' 
#' This is to estimate the the time-dependent AUC with confidence interval for each cause of failure.

#' @param beta.hat   Estimated regression coefficients
#' @param V  The asymptotic variance returned from the function covariance_cal
#' @param t0 A vector of time
#' @param cov.val A vector where each element represents a covariate value, ordered to correspond to the elements in \code{beta.hat} for each covariate
#' @param nf A numeric value of 3 or 7 for the number of polynomials. When nf = 7, the vector of polynomials = c(t^(-2),t^(-1),t^(-.5),log(t),t^(.5),t,t^2). When nf = 3, the vector of polynomials = c(log(t),t, t^2)
#' @returns A matrix of confidence intervals for the time-dependent AUC for each cause of failure.
#' The matrix has three rows and the number of columns equal to the length of \code{t0}. The first row contains the estimated AUC, the second row contains the lower bound of the confidence interval, and the third row contains the upper bound of the confidence interval.
#' 
#' @examples
#' \dontrun{
#' 
#' auc_pred(c(1,0.2,0.1,0.3,0.3,0.5),matrix(rep(1,36),nrow=6),c(1,2),c(1,0), nf=3)
#'
#' } 
#' @export
auc_pred <- function(beta.hat,V,t0,cov.val, 
                        nf=7)
{
  n0 <- length(beta.hat)-nf-1
  se.store.1<- matrix(0,ncol=length(t0),nrow=3)
  for(i in 1:length(t0)){
    if (nf == 3) {
      poly <- c(1, log(t0[i]), t0[i], t0[i]^2, cov.val)
    }else if (nf == 7){
      poly <- c(1,t0[i]^(-2),t0[i]^(-1),t0[i]^(-.5),log(t0[i]),t0[i]^(.5),t0[i],t0[i]^2,cov.val)
    }
    
    temp <- c(beta.hat %*% poly)
    lower <- temp- 1.96*sqrt(max(0,c(poly %*% V %*% poly)) )
    upper <- temp+ 1.96*sqrt(max(0,c(poly %*% V %*% poly)) )
    
    se.store.1[1,i] <- 1/(1+exp(-temp))
    se.store.1[2,i] <- 1/(1+exp(-lower))
    se.store.1[3,i] <-  1/(1+exp(-upper))
  }
  return(se.store.1)
}


#' Sum a Numeric Vector Using a C Function
#'
#' This wrapper provides an R interface to the `covariance_cal` C++ function.
#' @name covariance_cal_wrapper
#' @param a A vector of regression coefficients.
#' @param b  A vector of biomarker values.
#' @param c The \code{rho01} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param d The \code{rho02} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{1-Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param e A transposed matrix where each row represents a covariate, including an intercept as the first row. The rows are ordered to correspond with the regression coefficients in \code{a}.
#' @param f The \code{i - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param g The \code{l - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @return A list containing:
#'\itemize{
#'  \item \code{sigma1}: The matrix \eqn{\Sigma_1}, used in calculating the asymptotic variance.
#'  \item \code{sigma2}: The matrix \eqn{\Sigma_2}, also used in calculating the asymptotic variance.
#'  \item \code{v}: The asymptotic variance \eqn{V}.
#' }
#' @export
covariance_cal_wrapper <- function(a,b,c,d,e,f,g) {
  # Call the C function
  covariance_cal(a, b, c, d, e, f, g)
}
