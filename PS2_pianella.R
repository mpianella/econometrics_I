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
# define functions and arguments

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
comp_cdf_normal <- function(obs, n, x_tilde){
    mu_MLE <- mean(obs, na.rm = TRUE)
    sd_obs_MLE <- 1/n * sum((obs - x_tilde)^2)
    value_cdf <- pnorm(x_tilde, mean = mu_MLE, sd = sd_obs_MLE)
    return(value_cdf)
}

# sample n observations from the exponential distribution with lambda = 1
lambda <- 1
n <- c(50, 100, 250, 1000)
obs <- list()
obs <- lapply(n, function(x) rexp(x, rate = lambda))

# define remaining parameter for the functions
x_tilde <- c(1/2, 1, 2, 3)
n_iter <- 500

# generate data for the case of the normal distribution 
data_normal <- data.frame()
set.seed(314)
for (i in 1:length(x_tilde)){
    x_t <- x_tilde[i]
    data <- data.frame()
    for (j in 1:n_iter) {
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cdf_normal(obs = x, n = length(obs), x_tilde = x_t))
        d <- cbind(y, n_vec)
        data <- rbind(data, d)
    }
    data$x_tilde <- x_tilde[i]
    data_normal <- rbind(data_normal, data)
}

# generate data for the case of the first exponential distribution 
data_expon_1 <- data.frame()
set.seed(3141)
for (i in 1:length(x_tilde)){
    x_t <- x_tilde[i]
    data <- data.frame()
    for (j in 1:n_iter) {
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cfd_expon_1(obs = x, x_tilde = x_t))
        d <- cbind(y, n_vec)
        data <- rbind(data, d)
    }
    data$x_tilde <- x_tilde[i]
    data_expon_1 <- rbind(data_normal, data)
}

# generate data for the case of the second exponential distribition
data_expon_2 <- data.frame()
set.seed(31415)
for (i in 1:length(x_tilde)){
    x_t <- x_tilde[i]
    data <- data.frame()
    for (j in 1:n_iter) {
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cfd_expon_2(obs = x, x_tilde = x_t))
        d <- cbind(y, n_vec)
        data <- rbind(data, d)
    }
    data$x_tilde <- x_tilde[i]
    data_expon_2 <- rbind(data_normal, data)
}

# generate data for the true distribution 
data_true <- data.frame(cdf_true = sapply(x_tilde, pexp), x_tilde = x_tilde)

