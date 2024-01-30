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

#' Compute the cumulative density using an exponential cdf and lambda equals to 
#' the MLE estimate
#' 
#' @param {obs} {a vector of observations}
#' @param {x_tilde} {a scalar that defines value the argument in the cumulative distribution function}
comp_cfd_expon_1 <- function(obs, x_tilde){
    mean_obs <- mean(obs, na.rm = TRUE) 
    lambda_mle <- 1/mean_obs # define lambda mle as the inverse of the sample mean
    value_cdf <- 1 - exp(-lambda_mle * x_tilde)
    return(value_cdf)
}

#' Compute the cumulative density using an exponential cdf and
#'  lambda equals the mean of the first two moments of MLE
comp_cdf_expon_2 <- function(obs, x_tilde){
    mean_obs <- mean(obs, na.rm = TRUE)
    sd_obs <- sd(obs, na.rm = TRUE)
    mean_first_two_moments <- ((1/mean_obs)/2) + ((1/sd_obs)/2)
    value_cdf <- 1 - exp(-mean_first_two_moments * x_tilde)
    return(value_cdf)
}

#' Compute the cumulative density using a normal cdf and mean and variance 
#' are the MLE estimate in the normal case. 
#' @param {n} {number of observations}
comp_cdf_normal <- function(n, obs, x_tilde){
    mu_MLE <- mean(obs, na.rm = TRUE)
    sd_obs_MLE <- 1/n * sum((obs - x_tilde)^2)
    value_cdf <- pnorm(x_tilde, mean = mu_MLE, sd = sd_obs_MLE)
    return(value_cdf)
}


# sample n observations from the exponential distribution with lambda = 1
lambda <- 1
obs <- list()
obs <- lapply(n, function(x) rexp(x, rate = lambda))



x_tilde <- c(1/2, 1, 2, 3)
n_iter <- 