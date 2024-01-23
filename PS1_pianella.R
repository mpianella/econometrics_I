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

# ex 4.b.i ------
## I use a binomial distribution to run 100 simulation for each possible scenario of N
N <- c(50, 100, 250, 500, 1000)
p <- 1/2 # this is the probability of selecting a female in a Bernoulli trial 
k <- 40 # number of students in each class
beta <- 1/3
d_beta_hat <- list()
for(i in 1:length(N)){
    d <- c()
    for(j in 1: 100){
        n_female_p_class <- rbinom(N[i], size = k, prob = p) # generate a list of values for each possible scenario of N
        fract_fem_p_class <- n_female_p_class/k
        epsilon <- rnorm(N[i])  # generate a list of values for epsilon
        y <- beta * fract_fem_p_class + epsilon
        data_lm <- data_frame(y, epsilon)
        model <- lm(y ~ 0 + fract_fem_p_class)
        x <- coef(model)
        d <- rbind(d, x)
    }
    d_beta_hat[[i]] <- data.frame(beta_hat = d, size = rep(N[i], nrow(d)), row.names = NULL)
}







