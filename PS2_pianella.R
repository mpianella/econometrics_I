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

library(dplyr)
library(ggplot2)
library(gmm)
library(stats)

# please set your directory to the folder econometrics_I

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
#' @param {n_obs} {number of observations}
comp_cdf_normal <- function(obs, n_obs, x_tilde){
    mu_MLE <- mean(obs, na.rm = TRUE)
    sd_obs_MLE <- sqrt(1/n_obs * sum((obs - x_tilde)^2))
    value_cdf <- pnorm(x_tilde, mean = mu_MLE, sd = sd_obs_MLE)
    return(value_cdf)
}

#' Compute the cumulative density the share of the sample below x_tilde as estimator.
comp_less_x_tilde <- function(obs, x_tilde){
    value_cdf <- sum(obs < x_tilde)/length(obs)
    return(value_cdf)
}

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
        lambda <- 1 # sample n observations from the exponential distribution with lambda = 1
        n_obs <- c(50, 100, 250, 1000)
        obs <- lapply(n_obs, function(x) rexp(x, rate = lambda))
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cdf_normal(obs = x, n_obs = length(obs), x_tilde = x_t))
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
        lambda <- 1 # sample n observations from the exponential distribution with lambda = 1
        n_obs <- c(50, 100, 250, 1000)
        obs <- lapply(n_obs, function(x) rexp(x, rate = lambda))
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cfd_expon_1(obs = x, x_tilde = x_t))
        d <- cbind(y, n_vec)
        data <- rbind(data, d)
    }
    data$x_tilde <- x_tilde[i]
    data_expon_1 <- rbind(data_expon_1, data)
}

# generate data for the case of the second exponential distribition
data_expon_2 <- data.frame()
set.seed(31415)
for (i in 1:length(x_tilde)){
    x_t <- x_tilde[i]
    data <- data.frame()
    for (j in 1:n_iter) {
        lambda <- 1 # sample n observations from the exponential distribution with lambda = 1
        n_obs <- c(50, 100, 250, 1000)
        obs <- lapply(n_obs, function(x) rexp(x, rate = lambda))
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cdf_expon_2(obs = x, x_tilde = x_t))
        d <- cbind(y, n_vec)
        data <- rbind(data, d)
    }
    data$x_tilde <- x_tilde[i]
    data_expon_2 <- rbind(data_expon_2, data)
}

# generate data for the case of the less than x_tilde estimator
data_less_x_tilde <- data.frame()
set.seed(314159)
for (i in 1:length(x_tilde)){
    x_t <- x_tilde[i]
    data <- data.frame()
    for (j in 1:n_iter) {
        lambda <- 1 # sample n observations from the exponential distribution with lambda = 1
        n_obs <- c(50, 100, 250, 1000)
        obs <- lapply(n_obs, function(x) rexp(x, rate = lambda))
        n_vec <- sapply(obs, length)
        y <- sapply(obs, function(x) comp_cdf_less_x_tilde(obs = x, x_tilde = x_t))
        d <- cbind(y, n_vec)
        data <- rbind(data, d)
    }
    data$x_tilde <- x_tilde[i]
    data_less_x_tilde <- rbind(data_less_x_tilde, data)
}



# generate data for the true distribution 
data_true <- data.frame(cdf_true = sapply(x_tilde, pexp), x_tilde = x_tilde)

# data preprocessing for plot 
data1 <- data.frame(data_expon_1, distrib = "exponential1")
data2 <- data.frame(data_expon_2, distrib = "exponential2")
data3 <- data.frame(data_normal, distrib = "normal")
data4 <- data.frame(data_less_x_tilde, distrib = "less_x_tilde")

data_all <- rbind(data1, data2, data3, data4) %>% 
    left_join(data_true, by = "x_tilde") %>% 
    rename(estim_cum_prob = y, n_obs = n_vec) 

data_plot_ex_3 <- data_all %>% 
    group_by(distrib, x_tilde, n_obs) %>% 
    dplyr::summarise(emp_avg_bias =  mean(estim_cum_prob - cdf_true), emp_avg_var = var(estim_cum_prob), 
              MSE = emp_avg_bias^2 + emp_avg_var) 
    
# plot of average bias
plot_ex_3_e1 <- ggplot(data_plot_ex_3, aes(x = as.factor(n_obs), y = emp_avg_bias, color = as.factor(distrib))) +
    geom_point(size = 4) +
    scale_color_viridis_d() +
    theme_bw() + 
    labs(title = "Empirical average bias by x_tilde and by estimator", color = "Distribution", 
         y = "Emp avg bias", x = "Number of observations") +
    facet_wrap(~ x_tilde)
plot_ex_3_e1
ggsave("./png/ex_3_e1.png")


plot_ex_3_e2 <- ggplot(data_plot_ex_3, aes(x = as.factor(n_obs), y = emp_avg_var, color = as.factor(distrib))) +
    geom_point(size = 3) +
    scale_color_viridis_d() +
    labs(title = "Empirical average variance by x_tilde and by estimator", color = "Distribution", 
         y = "Emp avg variance", x = "Number of observations") +
    theme_bw() + 
    facet_wrap(~ x_tilde)
plot_ex_3_e2
ggsave("./png/ex_3_e2.png")


plot_ex_3_e3 <- ggplot(data_plot_ex_3, aes(x = as.factor(n_obs), y = MSE, color = as.factor(distrib))) +
    geom_point(size = 3) +
    scale_color_viridis_d() + 
    labs(title = "MSE by x_tilde and by estimator", color = "Distribution", 
         y = "MSE", x = "Number of observations") +
    theme_bw() + 
    facet_wrap(~ x_tilde)
plot_ex_3_e3
ggsave("./png/ex_3_e3.png")

#ex.4.e ---------
## please set your working directory to the econometrics_I folder

gmmdata <- read.csv(file.path("csv", "gmmdata.csv"))

# prepare data for the analysis
data_ex4 <- gmmdata %>% 
    mutate(r_bar = mean(r), r_lead1 = dplyr::lead(r), c_lead1 = dplyr::lead(c))

#' Write the g function to be used for computation
#'
#'@param {beta} {description}
#'@param {gamma} {description}
#'@param {x} {data}
#'
#'@returns {it returns a t-1 x 2 matrix}
g <- function(tet, x){
    beta <- tet[1]
    gamma <- tet[2]
    c <- x$c 
    r <- x$r
    c_lead1 <- x$c_lead1 
    r_lead1 <- x$r_lead1
    r_bar  <- x$r_bar
    g_b <- c(c^(gamma) - (beta * r_lead1 * (c_lead1 ^ (gamma - 1))))[-length(c)] # dropping the last observation because NA
    g_c <- c(c^(gamma-1) * (r_lead1 - r_bar) )[-length(c)] # dropping last because NA
    f <- cbind(g_b, g_c)
    return(f)
}

#' Writing the Jacobian of g
Dg <- function(tet, x){
    beta <- tet[1]
    gamma <- tet[2]
    c <- head(x$c, -1)
    r <- head(x$r, -1)
    c_lead1 <- head(x$c_lead1, -1) 
    r_lead1 <- head(x$r_lead1, -1)
    r_bar  <- head(x$r_bar, -1)
    jacobian <- matrix(c(mean(- r_lead1 * (c_lead1 ^(gamma -1)), na.rm = TRUE), 0, 
                         mean(c^(gamma) * log(c) - (beta * r_lead1 * c_lead1^(gamma-1)) * log(c_lead1), na.rm = TRUE),
                         mean(c^(gamma-1) *(r_lead1 - r_bar) *log(c), na.rm = TRUE)), nrow = 2, ncol =2)
    return(jacobian)
}

gmm(g = g, gradv = Dg, x = data_ex4, t0 = c(1,1), tol = 1e-117)

