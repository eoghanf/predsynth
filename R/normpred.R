#------------------------------------------------------------------------------
#
#                                 predsynth
#
#                                    by
#
#                              Eoghan Flanagan
#
#------------------------------------------------------------------------------
#
#
# An implementation/adaptation of the methods presented in the paper "Dynamic
# Bayesian Predictive Sythesis in Time Series Forecasting by Kenichiro McAlinn &
# Mike West. This implementation assumes that the past and future predictive
# distributions are Normal.
#
#
#------------------------------------------------------------------------------

#' An implementation of Bayesian Predictive Synthesis for Time Series
#' Forecasting
#'
#' \code{predsynth} synthesises predictive distributions to produce future
#' predictions from a number of different models or agents
#'
#' @param past.preds a matrix containing the mean of the historic normal
#' forecast density. The matrix should have n rows, for each timepoint, and
#' p columns, for each of the p predictors
#'
#' @param past.vars an nxp matrix containing the variance of the historic
#' normal forecast density
#'
#' @param past.y an 1xn matrix or n-length vector containing the outcome of the
#' historic timeseries
#'
#' @return future.preds is the mean of the synthesised predictive distribution
#' @return future.vars is the variance of the synthesised predictive distribution
#
normpred <- function(past.preds, past.vars, past.y, future.pred, future.var,
                     W.disc.factor=0.95, v=0.1, MC.steps=2000) {

  # The covars argument must be a matrix of real values

  if (!(is.matrix(past.preds))){
    stop ("past.preds must be a matrix of predictions")
  }

  if (!(is.numeric(past.preds))){
    stop ("past.preds must be a matrix of predictions")
  }

  # define n and p as the number of rows, columns of the predictions matrix as
  # specified above.

  n <- nrow(past.preds)
  p <- ncol(past.preds)

  if (!(is.matrix(past.vars))){
    stop ("past.vars must be a matrix of variances")
  }

  if (!(is.numeric(past.vars))){
    stop ("past.vars must be a matrix of positive variances")
  }

  if (!((nrow(past.vars) == n) & (ncol(past.vars) == p))){
    stop("past.vars must have the same shape as past.preds")
  }

  if (!((sum(past.vars > 0)) == (n*p))){
    stop("All elements of past.vars must be positive")
  }


  if (!(is.vector(past.y))){
    if (!(is.matrix(past.y))){
      stop("past.y must be a vector of observations or 1xn matrix")
    }
    if (ncol(past.y) != 1){
      stop("past.y must be a vector of observations or 1xn matrix")
    }
  }

  if (!(is.numeric(past.y))){
    stop ("past.y must be a vector or matrix of real values")
  }

  past.y<-matrix(past.y, ncol=1)

  if (!(is.vector(future.pred))){
    if (!(is.matrix(future.pred))){
      stop("future.pred must be a vector of predictions or px1 matrix")
    }
    if (!((ncol(future.pred) = 1) & (nrow(future.pred) =p))){
      stop("future.pred must be a vector of predictions or px1 matrix")
    }
  }

  if (!(is.vector(future.var))){
    if (!(is.matrix(future.var))){
      stop("future.var must be a vector of variances or px1 matrix")
    }
    if (!((ncol(future.var) = 1) & (nrow(future.var) =p))){
      stop("future.var must be a vector of variances or px1 matrix")
    }
  }

  # The past predictions and past prediction variances define the h distribution
  # completely in the Gaussian case

  h.means <- past.preds
  h.vars <- past.vars

  # Initialize the feature matrix F, including the dummy 1 feature in the
  # first row

  F <- matrix(0, nrow=n, ncol=p+1)
  F[,1] <- 1

  # sample the predictive densities to provide the initial 'latent states'
  # This is a n x (p+1) matrix
  F[,-1] <- sample.func(h.means, h.vars)

  y.out<-matrix(0, ncol=20, nrow=1)

  for (ss in (seq(MC.steps))){

    # Use FFBS to sample the latent state coefficients. This returns
    # a n x (p+1) matrix which is a draw from the distribution theta
    output <- bayesdlm::dlmffbs(y=past.y, X=F, W.disc.factor=W.disc.factor,
                               num.samples=1, v=v, prior.var=100)

    # Update the latent state distributions with a "Kalman filter" type
    # update. The forecast error c(t) is a n x 1 vector
    theta <- output[[1]]
    e <- matrix(past.y - rowSums (F * theta), nrow=n, ncol=p, byrow=TRUE)

    g <- rep(v, n) + rowSums((F[,-1] ^2) * h.vars)

    # Calculate the Kalman gain K(t)
    # referred to as b_t in McAlinn, West
    K <- h.vars * F[,-1] / g

    # Kalman update of previous plus forecast error times gain
    h.means <- h.means + (K * e)
    h.vars <- h.vars - ((K^2) * g)

    # Resample the latent states and recompose the F vector
    F <- matrix(0, nrow=n, ncol=p+1)
    F[,1] <- 1
    F[,-1] <- sample.func(h.means, h.vars)

    if (ss>1000 & ss%%50==0){

      y.out[(ss-1000)/50] <- sum(theta[n,]*F[n,])
    }

  }




  return(as.vector(y.out))
}


sample.func <- function(d.mean, d.vars){
  # returns a samples from the normal distribution with means given by the
  # matrix d.mean and variances given by the matrix d.vars. Note this is NOT a
  # multi-variate normal sampler.

  n<-nrow(d.vars)
  p<-ncol(d.vars)

  ans <- matrix(rnorm(n*p, d.mean, d.vars), nrow=n, ncol=p)

  return(ans)
}

