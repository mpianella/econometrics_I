## ---------------------------
##
## Script name: PS2_pianella.R
##
## Purpose of script: to solve the second problem set for the PhD course econometrics I
##
## Author: Matteo Pianella
##
## Date Created: 30th Jan 24
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# ex.3.e --------
n <- c(50, 100, 250, 1000)

#' Compute the cumulative density when lambda equals to the MLE estimate
#' 
#' @param {obs} {a vector of observations}
#' @param {x_tilde} {a scalar that defines value the argument in the cumulative distribution function}
comp_cfd_expon_1 <- function(obs, x_tilde) {
    mean_obs <- mean(obs, na.rm = TRUE) 
    lambda_mle <- 1/mean_obs # define lambda mle as the inverse of the sample mean
    value_cdf <- 1 - exp(-lambda_mle * x_tilde)
    return(value_cdf)
}

# sample n observations from the exponential distribution with lambda = 1
lambda <- 1
obs <- list()
obs <- lapply(n, function(x) rexp(x, rate = lambda))



x_tilda <- c(1/2, 1, 2, 3)
n_iter <- 