#' Estimate the average response value of the input data given a treatment regime
#'
#' @description This function estimates the average response value of the input
#' data given a 'DTR.KernSmooth' / 'DTR.Boots.KernSmooth' model object or an
#' estimated optimal treatment regime vector
#' @usage obj_value(X, y, a, object, beta, prob = 0.5)
#' @param X Input matrix, of dimension n_obs x n_vars; each row is an observation vector.
#' @param y Response variable to be maximized on average if every subject follows
#' the treatment recommended by the optimal regime.
#' @param a Received treatments for n_obs subjects. Must be bivariate, and labeled as \{0,1\}.
#' @param object Fitted "DTR.KernSmooth" or "DTR.Boots.KernSmooth" model object.
#' @param beta The treatment regime vector. Cannot be missing if "object" is not provided.
#' @param prob The propensity score for n_obs subjects, i.e., P(a=1|X). If \code{NULL},
#' it would be estimated by logistic regression a~X.
#' @details \code{object} and \code{beta} cannot be both missing. If the input
#' data (X, y, a) is missing but \code{object} is provided, the function will
#' return the optimal value of the input object.
#' @return The estimated average response value if all n_obs subjects follows the
#' treatment recommendations according to the fitted model or the estimated
#' treatment regime.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{DTR.KernSmooth}}, \code{\link{DTR.Boots.KernSmooth}}
#' @references Wu, Y. and Wang, L. (2021),
#' \emph{Resampling-based Confidence Intervals for Model-free Robust Inference
#' on Optimal Treatment Regimes, Biometrics, 77: 465– 476}, \doi{10.1111/biom.13337}.
#' @examples
#' n <- 5e3; p <- 3
#' beta <- c(0.2,1,-0.5,-0.8)
#' beta1 <- c(1,-0.5,-0.5,0.5)
#'
#' X <- matrix(rnorm(n*p),n)
#' a <- rbinom(n,1,0.5)
#' mean1 <- exp(cbind(1,X) %*% beta1)
#' mean2 <- 8/(1 + exp(-cbind(1,X) %*% beta)) - 4
#' y <- mean1 + a * mean2 + rnorm(n)
#'
#' smooth_model <- DTR.KernSmooth(X, y, a)
#' boots_smooth_model <- DTR.Boots.KernSmooth(X, y, a)
#'
#' newX <- matrix(rnorm(n*p),n)
#' newa <- rbinom(n,1,0.5)
#' newmean1 <- exp(cbind(1,newX) %*% beta1)
#' newmean2 <- 8/(1 + exp(-cbind(1,newX) %*% beta)) - 4
#' newy <- newmean1 + newa * newmean2 + rnorm(n)
#'
#' obj_value(newX, newy, newa, smooth_model)
#' obj_value(newX, newy, newa, boots_smooth_model)
#' obj_value(newX, newy, newa, beta = smooth_model$beta_smooth)
#'
obj_value<-function(X, y, a, object, beta, prob = 0.5){
  if(missing(object)) {
    if(missing(beta)){
      stop("Please supply either a fitted 'DTR.KernSmooth' / 'DTR.Boots.KernSmooth' model object,
           or a fitted optimal treatment regime vector beta")
    }
    if(missing(X)|missing(y)|missing(a)){
      stop("Please supply the input data (X, y, a) to estimate the optimal value")
    }
    if(!all(a%in%c(0,1))){
      stop("The treatment 'a' should be labeled as 0 or 1")
    }
    if(is.null(prob)){
      prob_model <- glm(a ~ X, family = binomial)
      prob <- prob_model$fitted.values
    }else{
      if(any(prob<=0)|any(prob>=1)){
        stop("The propensity score 'prob' should be between 0 and 1")
      }
      nobs <- length(y)
      prob <- rep(prob,length=nobs)
    }
    res <- obj_value_C(cbind(1,X), y, a, beta, prob)
  } else{
    if(missing(X)|missing(y)|missing(a)){
      if(attr(object, "class")=="DTR.KernSmooth"){
        res <- object$value_smooth
      } else if(attr(object, "class")=="DTR.Boots.KernSmooth"){
        res <- object$smooth_est$value_smooth
      } else{
        stop("Please supply a valid 'DTR.KernSmooth' or 'DTR.Boots.KernSmooth' object")
      }
      warning("The input data (X, y, a) is not complete. Return the optimal value of the input object")
    }else{
      if(!all(a%in%c(0,1))){
        stop("The treatment a should be labeled as 0 or 1")
      }
      if(attr(object, "class")=="DTR.KernSmooth"){
        beta <- object$beta_smooth
      } else if(attr(object, "class")=="DTR.Boots.KernSmooth"){
        beta <- object$smooth_est$beta_smooth
      } else{
        stop("Please supply a valid 'DTR.KernSmooth' or 'DTR.Boots.KernSmooth' object")
      }
      if(is.null(prob)){
        prob_model <- glm(a ~ X, family = binomial)
        prob <- prob_model$fitted.values
      }else{
        nobs <- length(y)
        prob <- rep(prob,length=nobs)
      }
      res <- obj_value_C(cbind(1,X), y, a, beta, prob)
    }
  }
  return(res)
}


#' Predict the optimal treatments given a 'DTR.KernSmooth' object
#'
#' @description This function predicts the optimal treatments for new subjects
#' from a fitted DTR.KernSmooth model.
#' @param object Fitted "DTR.KernSmooth" model object.
#' @param newX Matrix of new values for X at which predictions are to be made.
#' @param ... Not used. Other arguments to predict.
#' @details All the predicted optimal treatments are labeled as \{0,1\}.
#' @return A vector of predicted optimal treatments for the new subjects given
#' the fitted DTR.KernSmooth model.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{predict.DTR.Boots.KernSmooth}}, \code{\link{DTR.KernSmooth}},
#' \code{\link{DTR.Boots.KernSmooth}}
#' @references Wu, Y. and Wang, L. (2021),
#' \emph{Resampling-based Confidence Intervals for Model-free Robust Inference
#' on Optimal Treatment Regimes, Biometrics, 77: 465– 476}, \doi{10.1111/biom.13337}.
#' @examples
#' n <- 5e3; p <- 3
#' beta <- c(0.2,1,-0.5,-0.8)
#' beta1 <- c(1,-0.5,-0.5,0.5)
#'
#' X <- matrix(rnorm(n*p),n)
#' a <- rbinom(n,1,0.5)
#' mean1 <- exp(cbind(1,X) %*% beta1)
#' mean2 <- 8/(1 + exp(-cbind(1,X) %*% beta)) - 4
#' y <- mean1 + a * mean2 + rnorm(n)
#'
#' smooth_model <- DTR.KernSmooth(X, y, a)
#' newX <- matrix(rnorm(n*p),n)
#' predict(smooth_model, newX)
#'
#' @method predict DTR.KernSmooth
#'
predict.DTR.KernSmooth<-function(object, newX, ...){
  if (missing(newX)) {
    stop("Please supply a value for 'newX'")
  }
  if(attr(object, "class")=="DTR.KernSmooth"){
    opt.trt <- (cbind(1,newX)%*%object$beta_smooth>0)*1
  } else{
    stop("Please supply a valid 'DTR.KernSmooth' object")
  }
  return(as.numeric(opt.trt))
}

#' Predict the optimal treatment given a 'DTR.Boots.KernSmooth' object
#'
#' @description This function predicts the optimal treatments for new subjects
#' from a fitted DTR.Boots.KernSmooth model.
#' @param object Fitted "DTR.Boots.KernSmooth" model object.
#' @param newX Matrix of new values for X at which predictions are to be made.
#' @param ... Not used. Other arguments to predict.
#' @details All the predicted optimal treatments are labeled as \{0,1\}.
#' @return A vector of predicted optimal treatments for the new subjects given
#' the fitted DTR.Boots.KernSmooth model.
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{predict.DTR.KernSmooth}}, \code{\link{DTR.KernSmooth}},
#' \code{\link{DTR.Boots.KernSmooth}}
#' @references Wu, Y. and Wang, L. (2021),
#' \emph{Resampling-based Confidence Intervals for Model-free Robust Inference
#' on Optimal Treatment Regimes, Biometrics, 77: 465– 476}, \doi{10.1111/biom.13337}.
#' @examples
#' n <- 5e3; p <- 3
#' beta <- c(0.2,1,-0.5,-0.8)
#' beta1 <- c(1,-0.5,-0.5,0.5)
#'
#' X <- matrix(rnorm(n*p),n)
#' a <- rbinom(n,1,0.5)
#' mean1 <- exp(cbind(1,X) %*% beta1)
#' mean2 <- 8/(1 + exp(-cbind(1,X) %*% beta)) - 4
#' y <- mean1 + a * mean2 + rnorm(n)
#'
#' boots_smooth_model <- DTR.Boots.KernSmooth(X, y, a)
#' newX <- matrix(rnorm(n*p),n)
#' predict(boots_smooth_model, newX)
#'
#' @method predict DTR.Boots.KernSmooth
#'
predict.DTR.Boots.KernSmooth<-function(object, newX, ...){
  if (missing(newX)) {
    stop("Please supply a value for 'newX'")
  }
  if(attr(object, "class")=="DTR.Boots.KernSmooth"){
    opt.trt <- (cbind(1,newX)%*%object$smooth_est$beta_smooth>0)*1
  } else{
    stop("Please supply a valid 'DTR.Boots.KernSmooth' object")
  }
  return(as.numeric(opt.trt))
}

#' Estimate the optimal treatment regime among all linear regimes with smoothed
#' estimation methods
#'
#' @description This function estimates the optimal treatment regime among all
#' linear regimes with smoothed estimation methods, and outputs a 'DTR.KernSmooth'
#' model object
#' @usage DTR.KernSmooth(X, y, a, prob = 0.5, kernel = "normal",
#'                phi0 = 1, gamma = 2, err_tol = 1e-8,
#'                iter_tol = 100)
#' @param X Input matrix, of dimension n_obs x n_vars; each row is an observation vector.
#' @param y Response variable to be maximized on average if every subject follows
#' the treatment recommended by the optimal regime.
#' @param a Received treatments for n_obs subjects. Must be bivariate, and labeled as \{0,1\}.
#' @param prob The propensity score for n_obs subjects, i.e., P(a=1|X). If \code{NULL},
#' it would be estimated by logistic regression a~X.
#' @param kernel The kernel function to be used in smoothed estimation. Should be
#' one of "normal" and "poly". See more details in "Details".
#' @param phi0 The initial step size to be used in the Proximal Algorithm.
#' @param gamma The multiplier of the step sizes to be used in the Proximal
#' Algorithm. Must be 0 < gamma < 1.
#' @param err_tol The desired accuracy in the estimation.
#' @param iter_tol The maximum number of iterations in the estimation algorithm.
#' @details This function estimates the optimal linear treatment regime to maximizes
#' the average outcome among the population if every individual follows the treatment
#' recommended by this treatment regime.\cr
#' Assume the propensity score \eqn{\pi(\bm{x})=P(A=1|\bm{x})}  can be modeled as
#' \eqn{\pi(\bm{x},\bm{\xi})} where \eqn{\bm{\xi}} is a finite-dimensional parameter
#' (e.g., via logistic regression). Let \eqn{\widehat{\bm{\xi}}} be an estimate
#' of \eqn{\bm{\xi}}. Hence, our goal is to estimate \eqn{\bm{\beta}} which maximizes:
#' \deqn{V_n(\bm{\beta})=n^{-1}\sum_{i=1}^n \frac{\left[A_i I\left(\bm{x}_i^T\bm{\beta}>0\right)+(1-A_i)I\left(\bm{x}_i^T\bm{\beta}\leq 0\right)\right]Y_i}
#' {A_i\pi(\bm{x}_i, \widehat{\bm{\xi}}) +(1-A_i)\left[1-\pi(\bm{x}_i, \widehat{\bm{\xi}})\right]}.}
#' For the identifability, we normalize the estimator such that the second element
#' has magnitude 1, i.e., \eqn{|\widehat{\beta}_2|=1}.\cr
#' To alleviates the computational challenge due to the nonsmooth indicator function,
#' and derive asymptotic distribution of the estimators, we consider to use a smoothed
#' function \eqn{K(\cdot)} to approximate the indicator function \eqn{I(\cdot)}.
#' That is, we will estimate \eqn{\bm{\beta}} which maximizes:
#' \deqn{n^{-1}\sum_{i=1}^n \frac{\left[A_i K\left(\frac{\bm{x}_i^T\bm{\beta}}{h_n}\right)+(1-A_i)\left\{1-K\left(\frac{\bm{x}_i^T\bm{\beta}}{h_n}\right)\right\}\right]Y_i}
#' {A_i\pi(\bm{x}_i, \widehat{\bm{\xi}}) +(1-A_i)\left[1-\pi(\bm{x}_i, \widehat{\bm{\xi}})\right]}.}
#' In this function, we provide two options for the smoothed kernel functions:
#' \itemize{
#'  \item{"normal" }{ The c.d.f of N(0,1) distribution.}
#'  \item{"poly" }{   A polynomial function \eqn{K(v) =\left[0.5 + \frac{105}{64}\{\frac{v}{5}-\frac{5}{3}(\frac{v}{5})^3 +\frac{7}{5}(\frac{v}{5})^5 - \frac{3}{7}(\frac{v}{5})^7\}\right]I( -5\leq v \leq 5)+I(v>5)}.}
#'  }
#' The bandwidth is set as \eqn{h_n=0.9n^{-0.2} \min\{std (\bm{x}_i^T\bm{\beta}),IQR(\bm{x}_i^T\bm{\beta})/1.34\}}.\cr
#' To solve the non-convexity problem of the optimization, we employ a proximal
#' gradient descent algorithm for estimation. See more details in the reference.
#' @return An object of class "DTR.KernSmooth", which is a list containing at
#' least the following components:
#'  \item{X}{The input matrix used.}
#'  \item{y}{The response variable used.}
#'  \item{a}{The treatment vector received by each subject.}
#'  \item{prob}{The propensity score vector for each subject.}
#'  \item{kernel}{The kernel function used in smoothed estimation.}
#'  \item{beta_smooth}{The estimated optimal treatment regime vector.}
#'  \item{opt_treatment}{The predicted optimal treatments for the input data
#'  given the estimated optimal regime.}
#'  \item{value_smooth}{The estimated optimal average response value among all
#'  linear treatment regimes.}
#'  \item{converge}{Logical. \code{TRUE} if the estimation algorithm converges,
#'  and \code{FALSE} if not.}
#'  \item{iter_num}{The number of iterations used for the algorithm convergence.}
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{predict.DTR.KernSmooth}}, \code{\link{obj_value}},
#' \code{\link{DTR.Boots.KernSmooth}}
#' @references Wu, Y. and Wang, L. (2021),
#' \emph{Resampling-based Confidence Intervals for Model-free Robust Inference
#' on Optimal Treatment Regimes, Biometrics, 77: 465– 476}, \doi{10.1111/biom.13337}.\cr
#' Nesterov, Y. (2007).
#' \emph{Gradient methods for minimizing composite objective function. Core
#' discussion papers, Université catholique de Louvain, Center for Operations
#' Research and Econometrics (CORE)}.
#' @examples
#' n <- 5e3; p <- 3
#' beta <- c(0.2,1,-0.5,-0.8)
#' beta1 <- c(1,-0.5,-0.5,0.5)
#'
#' X <- matrix(rnorm(n*p),n)
#' a <- rbinom(n,1,0.5)
#' mean1 <- exp(cbind(1,X) %*% beta1)
#' mean2 <- 8/(1 + exp(-cbind(1,X) %*% beta)) - 4
#' y <- mean1 + a * mean2 + rnorm(n)
#'
#' smooth_model <- DTR.KernSmooth(X, y, a)
#' smooth_model$beta_smooth
#' smooth_model$value_smooth
#'
DTR.KernSmooth<-function(X, y, a, prob = 0.5, kernel = "normal", phi0 = 1, gamma = 2,
                         err_tol = 1e-8, iter_tol = 100){
  if(missing(X)|missing(y)|missing(a)){
    stop("Please supply the data (X, y, a) to estimate the optimal value")
  }
  if(!all(a%in%c(0,1))){
    stop("The treatment a should be labeled as 0 or 1")
  }
  if(!kernel %in% c("normal","poly")){
    stop("'kernel' should be one of 'normal' and 'poly'")
  }
  nvars <- dim(X)[2]
  initial <- rep(0,nvars+1)

  if(is.null(prob)){
    prob_model <- glm(a ~ X, family = binomial)
    prob <- prob_model$fitted.values
  }else{
    if(any(prob<=0)|any(prob>=1)){
      stop("The propensity score 'prob' should be between 0 and 1")
    }
    nobs <- dim(X)[1]
    prob <- rep(prob,length=nobs)
  }
  res <- Smooth_C(cbind(1,X), y, a, initial, prob, kernel,
                  phi0, gamma, err_tol, iter_tol)
  attr(res, "class") <- "DTR.KernSmooth"
  return(res)
}

#' Make inference about the linear treatment regime vector and the optimal value
#'
#' @description This function estimates the optimal treatment regime among all
#' linear regimes with smoothed estimation methods, and construct element-wise
#' 100(1-alpha)\% confidence intervals for the optimal linear treatment regime
#' vector, and the 100(1-alpha)\% confidence interval for the optimal value if
#' the population follows treatments recommended by the optimal linear regime.
#' It outputs a 'DTR.Boots.KernSmooth' model object
#' @usage DTR.Boots.KernSmooth(X, y, a, prob = 0.5,
#'                      B = 500, alpha = 0.05,
#'                      kernel = "normal",
#'                      phi0 = 1, gamma = 2,
#'                      err_tol = 1e-8,
#'                      iter_tol = 100)
#' @param X Input matrix, of dimension n_obs x n_vars; each row is an observation vector.
#' @param y Response variable to be maximized on average if every subject follows
#' the treatment recommended by the optimal regime.
#' @param a Received treatments for n_obs subjects. Must be bivariate, and labeled as \{0,1\}.
#' @param prob The propensity score for n_obs subjects, i.e., P(a=1|X). If \code{NULL},
#' it would be estimated by logistic regression a~X.
#' @param B The number of repetitions in the inference procedure by weighted bootstrap.
#' @param alpha The confidence level of the confidence interval.
#' @param kernel The kernel function to be used in smoothed estimation. Should be
#' one of "normal" and "poly". See more details in the "Details" section of
#' \code{\link{DTR.KernSmooth}}.
#' @param phi0 The initial step size to be used in the Proximal Algorithm.
#' @param gamma The multiplier of the step sizes to be used in the Proximal
#' Algorithm. Must be 0 < gamma < 1.
#' @param err_tol The desired accuracy in the estimation.
#' @param iter_tol The maximum number of iterations in the estimation algorithm.
#' @details This function constructs confidence intervals for the optimal linear
#' treatment regime vector by wild bootstrap procedures. The bootstrapped estimate
#' of the smoothed robust estimator is defined as the vector \eqn{\widehat{\bm{\beta}}^*}
#' that maximizes
#' \deqn{n^{-1}\sum_{i=1}^n \frac{\left[A_i K\left(\frac{\bm{x}_i^T\bm{\beta}}{h_n}\right)+(1-A_i)\left\{1-K\left(\frac{\bm{x}_i^T\bm{\beta}}{h_n}\right)\right\}\right]r_iY_i}{A_i\pi(\bm{x}_i, \widehat{\bm{\xi}}) +(1-A_i)\left[1-\pi(\bm{x}_i, \widehat{\bm{\xi}})\right]},}
#' where \eqn{r_i}'s are i.i.d N(1,1). \cr
#' Let \eqn{\xi_j^{\circ(\alpha/2)}} and \eqn{\xi_j^{\circ(1-\alpha/2)}} be the \eqn{(\alpha/2)}-th
#' and \eqn{(1-\alpha/2)}-th quantile of the bootstrap distribution of
#' \eqn{(nh_n)^{1/2}(\widetilde{\bm{\beta}}_j^*-\widetilde{\bm{\beta}}_j)}, respectively,
#' where \eqn{\alpha} is a small positive number. We can estimate \eqn{\xi_j^{\circ(\alpha/2)}}
#' and \eqn{\xi_j^{\circ(1-\alpha/2)}} from a large number of bootstrap samples. An asymptotic
#' \eqn{100(1-\alpha)\%} bootstrap confidence interval for \eqn{\beta_{0j}}, is given by
#' \deqn{\left\{\widetilde{\bm{\beta}}_j-(nh_n)^{-1/2}\xi_j^{\circ(1-\alpha/2)}, \widetilde{\bm{\beta}}_j-(nh_n)^{-1/2}\xi_j^{\circ(\alpha/2)}\right\}.}
#' To construct confidence intervals for the optimal value \eqn{V(\bm{\beta}_0)}, we define
#' \deqn{V_n^*(\widetilde{\bm{\beta}}) = n^{-1}\sum_{i=1}^n \frac{\left[A_i I\left(\bm{x}_i^T\widetilde{\bm{\beta}}>0\right)+(1-A_i)I\left(\bm{x}_i^T\widetilde{\bm{\beta}}\leq 0 \right) \right]r_iY_i}{A_i\pi(\bm{x}_i, \widehat{\bm{\xi}}) +(1-A_i)\left[1-\pi(\bm{x}_i, \widehat{\bm{\xi}})\right]},}
#' where \eqn{r_i}'s are i.i.d N(1,1). Let \eqn{d^{\circ(\alpha/2)}} and \eqn{d^{\circ(1-\alpha/2)}}
#' be the \eqn{(\alpha/2)}-th and \eqn{(1-\alpha/2)}-th quantile of the bootstrap
#' distribution of \eqn{n^{1/2}\{V_n^*(\widehat{\bm{\beta}})-V_n(\widehat{\bm{\beta}})\}},
#' respectively. An asymptotic \eqn{100(1-\alpha)\%} bootstrap confidence interval for
#' \eqn{V(\bm{\beta}_0)} is
#' \deqn{\left\{V_n(\widehat{\bm{\beta}})-n^{-1/2}d^{\circ(1-\alpha/2)}, V_n(\widehat{\bm{\beta}})-n^{-1/2}d^{\circ(\alpha/2)}\right\}.}
#' See more details in the reference.
#' @return An object of class "DTR.Boots.KernSmooth", which is a list containing
#' the following components:
#'  \item{alpha}{The confidence level of the confidence interval.}
#'  \item{B}{The number of repetitions in the inference procedure by weighted
#'  bootstrap.}
#'  \item{smooth_est}{The fitted "DTR.KernSmooth" object based on the input data
#'  and parameters.}
#'  \item{Beta_CI}{The 100(1-alpha)\% confidence intervals for each element of
#'  the optimal treatment regime vector.}
#'  \item{value_CI}{The 100(1-alpha)\% confidence interval for the optimal average
#'  response value among all linear treatment regimes.}
#' @author Yunan Wu and Lan Wang\cr Maintainer:
#' Yunan Wu <yunan.wu@@utdallas.edu>
#' @seealso \code{\link{predict.DTR.Boots.KernSmooth}}, \code{\link{obj_value}},
#' \code{\link{DTR.KernSmooth}}
#' @references Wu, Y. and Wang, L. (2021),
#' \emph{Resampling-based Confidence Intervals for Model-free Robust Inference
#' on Optimal Treatment Regimes, Biometrics, 77: 465– 476}, \doi{10.1111/biom.13337}.
#' @examples
#' n <- 5e3; p <- 3
#' beta <- c(0.2,1,-0.5,-0.8)
#' beta1 <- c(1,-0.5,-0.5,0.5)
#'
#' X <- matrix(rnorm(n*p),n)
#' a <- rbinom(n,1,0.5)
#' mean1 <- exp(cbind(1,X) %*% beta1)
#' mean2 <- 8/(1 + exp(-cbind(1,X) %*% beta)) - 4
#' y <- mean1 + a * mean2 + rnorm(n)
#'
#' boots_smooth_model <- DTR.Boots.KernSmooth(X, y, a)
#' boots_smooth_model$Beta_CI
#' boots_smooth_model$value_CI
#'

DTR.Boots.KernSmooth<-function(X, y, a, prob = 0.5, B = 500, alpha = 0.05, kernel = "normal",
                               phi0 = 1, gamma = 2, err_tol = 1e-8, iter_tol = 100){
  if(missing(X)|missing(y)|missing(a)){
    stop("Please supply the data (X, y, a) to estimate the optimal value")
  }
  if(!all(a%in%c(0,1))){
    stop("The treatment a should be labeled as 0 or 1")
  }
  if(!kernel %in% c("normal","poly")){
    stop("'kernel' should be one of 'normal' and 'poly'")
  }
  nvars <- dim(X)[2]
  nobs <- dim(X)[1]
  initial <- rep(0,nvars+1)

  if(is.null(prob)){
    prob_model <- glm(a ~ X, family = binomial)
    prob <- prob_model$fitted.values
  }else{
    if(any(prob<=0)|any(prob>=1)){
      stop("The propensity score 'prob' should be between 0 and 1")
    }
    prob <- rep(prob,length=nobs)
  }
  weights <- matrix(rnorm(nobs*B,1,1),nobs,B)
  res <- Boots_C(cbind(1,X), y, a, initial, prob, kernel, weights,
                 alpha, phi0, gamma, err_tol, iter_tol)
  attr(res, "class") <- "DTR.Boots.KernSmooth"
  return(res)
}
