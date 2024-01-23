## ---------------------------
##
## Script name: PS1_pianella.R
##
## Purpose of script: to provide solutions for the progamming part of the first
##                    problem set in econometrics I
##                  
## Author: Matteo Pianella
##
## Date Created: 23rd January 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(dplyr)
library(ggplot2)

# ex 4.a -------
## I run 819 simulations using a binomial distribution
s_bar <- .556
n_cities_between5500_and_6500 <- 102
set.seed(314)
n_sim <- 819
n_below_6000 <- rbinom(n_sim, n_cities_between5500_and_6500, s_bar)
s_c <- n_below_6000/102 # calculate s_c for each city 
sd_s_c <- sqrt((sum((s_c - s_bar)^2))/n_sim)
print(c(sd_s_c, min(s_c), max(s_c))) # my sd is consistently lower than what they have in Table A1 (0.118)

# ex 4.b.1 ------
## I use a binomial distribution to run 100 simulation for each possible scenario of N
N <- c(50, 100, 250, 500, 1000)
p <- 1/2 # this is the probability of selecting a female in a Bernoulli trial 
k <- 40  # number of students in each class
beta <- 1/3
n_sim <- 100
d_beta_hat <- list()
sigma_beta <- list()
get_sd <- function(N, p, k, beta, n_sim){
    for(i in 1:length(N)){
        d <- c()
        for(j in 1: n_sim){
            n_female_p_class <- rbinom(N[i], size = k, prob = p) # generate a list of values for each possible scenario of N
            fract_fem_p_class <- n_female_p_class/k
            epsilon <- rnorm(N[i])  # generate a list of values for epsilon
            y <- beta * fract_fem_p_class + epsilon
            data_lm <- data_frame(y, epsilon)
            model <- lm(y ~ 0 + fract_fem_p_class)
            x <- coef(model)
            d <- rbind(d, x)
        }
        # d_beta_hat[[i]] <- data.frame(beta_hat = d, size = rep(N[i], nrow(d)), row.names = NULL)
        sigma_beta[[i]] <- sd(d)
    }
    # # calculate standard deviations
    # sigma_beta <- do.call(rbind, lapply(d_beta_hat, function(x) sd(x[[1]])))
    # return(sigma_beta)
    return(sigma_beta)
    
}
set.seed(3141)
sigma_beta <- get_sd(N = N ,p=p,k=k,beta=beta,n_sim=n_sim)
d_plot_4b1 <- data_frame(sigma_beta = unlist(sigma_beta), size = N)

p <- ggplot(d_plot_4b1, aes(x = N, y = sigma_beta)) +
    geom_point() +
    labs(title = "Number of classes and standard deviation of beta coefficient", x = "N", y = "Sigma beta")
p

# ex 4.b.2 ----
# set up coefficients again
N <- c(50, 100, 250, 500, 1000)
p <- 1/2 # this is the probability of selecting a female in a Bernoulli trial 
k <- 40  # number of students in each class
beta <- 1/3
n_sim <- 1500
d_beta_hat <- list()
sigma_beta <- list()
k <- c(10, 20, 60)
data <- data_frame()
for(w in 1:length(k)){
    sigma_beta <- get_sd(N = N, p=p, k = k[w], beta = beta, n_sim = n_sim2)
    d <- data_frame(sigma_beta = unlist(sigma_beta), size = N, class_size = k[w])
    data <- rbind(data, d)
}

d <- data_frame(d_plot_4b1, class_size = 40)
data_plot_4b2 <- rbind(d, data)

p2 <- ggplot(data_plot_4b2, aes(x = size, y = sigma_beta, color = class_size)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Number of classes and standard deviation of beta coefficient", x = "N", y = "Sigma beta")
p2


