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
library(scales)



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
p <- 1/2 # this is the probability of selecting a male in a Bernoulli trial 
k <- 40  # number of students in each class
beta <- 1/3
n_sim <- 100
get_sd <- function(N, p, k, beta, n_sim){
    sigma_beta <- list()
    for(i in 1:length(N)){
        d <- c()
        for(j in 1: n_sim){
            n_male_p_class <- rbinom(N[i], size = k, prob = p) # generate a vector of number of males for each class
            fract_male_p_class <- n_male_p_class/k
            epsilon <- rnorm(N[i])  # generate a list of values for epsilon
            y <- beta * fract_male_p_class + epsilon
            data_lm <- data_frame(y, epsilon)
            model <- lm(y ~ 0 + fract_male_p_class)
            x <- coef(model)
            d <- rbind(d, x)
        }
        sigma_beta[[i]] <- sd(d)
    }
    return(sigma_beta)
}
set.seed(3141)
sigma_beta <- get_sd(N = N ,p=p,k=k,beta=beta,n_sim=n_sim)
d_plot_4b1 <- data_frame(sigma_beta = unlist(sigma_beta), size = N)

p4b1 <- ggplot(d_plot_4b1, aes(x = N, y = sigma_beta)) +
    geom_point() +
    labs(title = "Number of classes and standard deviation of beta coefficient", x = "N", y = "Sigma beta")
p4b1

# ex 4.b.2 ----
# set up coefficients 
N <- c(50, 100, 250, 500, 1000)
p <- 1/2 # this is the probability of selecting a male in a Bernoulli trial 
k <- c(10, 20, 60) # vector with number of students for each class
beta <- 1/3
n_sim <- 100
sigma_beta <- list()
data <- data_frame()
for(w in 1:length(k)){
    sigma_beta <- get_sd(N = N, p=p, k = k[w], beta = beta, n_sim = n_sim)
    d <- data_frame(sigma_beta = unlist(sigma_beta), size = N, class_size = k[w])
    data <- rbind(data, d)
}
d <- data_frame(d_plot_4b1, class_size = 40)
data_plot_4b2 <- rbind(d, data)

p4b2 <- ggplot(data_plot_4b2, aes(x = size, y = sigma_beta, color = factor(class_size))) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd")) +  
    labs(title = "Number of classes and standard deviation of beta coefficient", x = "N", y = "Sigma beta")
p4b2


#ex 4.b.4 -------
# calculate standard deviation for the estimate of her initial sample
N <- c(200)
k <- c(30)
p <- 1/2
beta <- 1/3
n_sim <- 100
sigma_beta <- list()
set.seed(31415)
sigma_beta <- get_sd(N = N, p=p, k = k, beta = beta, n_sim = n_sim)
print(sigma_beta)
ex_ante_data <- data_frame(sigma_beta = unlist(sigma_beta), N = 200, funding = "ex_ante")


# calculate the standard deviation when she uses funding in urban area
# modify the function to return beta_hat instead of sigma_beta
get_beta_hat <- function(N, p, k, beta, n_sim){
    beta_hat <- list()
    for(i in 1:length(N)){
        d <- c()
        for(j in 1: n_sim){
            n_male_p_class <- rbinom(N[i], size = k, prob = p) # generate a vector of number of males for each class
            fract_male_p_class <- n_male_p_class/k
            epsilon <- rnorm(N[i])  # generate a list of values for epsilon
            y <- beta * fract_male_p_class + epsilon
            data_lm <- data_frame(y, epsilon)
            model <- lm(y ~ 0 + fract_male_p_class)
            x <- coef(model)
            d <- rbind(d, x)
        }
        beta_hat <- d
    }
    return(beta_hat)
}

N <- c(200, 100)
k <- c(30, 40)
p <- 1/2
beta <- 1/3
n_sim <- 100
data <- data_frame()
for(w in 1:length(k)){
    beta_hat <- get_beta_hat(N = N[w], p=p, k = k[w], beta = beta, n_sim = n_sim)
    d <- data_frame(beta_hat = unlist(beta_hat), size = N[w], class_size = k[w])
    data <- rbind(data, d)
}
set.seed(314159)
sigma_beta <- sd(data$beta_hat)
urban_data <- data_frame(sigma_beta = sigma_beta, N = 300, funding = "urban")

# run the same exercise for rural

N <- c(200, 50)
k <- c(30, 15)
p <- 1/2
beta <- 1/3
n_sim <- 100
data <- data_frame()
for(w in 1:length(k)){
    beta_hat <- get_beta_hat(N = N[w], p=p, k = k[w], beta = beta, n_sim = n_sim)
    d <- data_frame(beta_hat = unlist(beta_hat), size = N[w], class_size = k[w])
    data <- rbind(data, d)
}
set.seed(3141593)
sigma_beta <- sd(data$beta_hat)
rural_data <- data_frame(sigma_beta = sigma_beta, N = 250, funding = "rural")

d_plot_4b4 <- rbind(ex_ante_data, urban_data, rural_data)

p4b4_1 <- ggplot(d_plot_4b4, aes(x = factor(N), y = sigma_beta, fill = factor(funding))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#9467bd")) +  
    labs(title = "Number of classes and standard deviation of beta coefficient", x = "N", y = "Sigma beta")
p4b4_1

# ex.4.c ----
N <- c(10, 20, 60)
p <- 1/2
data_4c <- data_frame()
for (i in 1:length(N)) {
    perc_25 <- floor((N[i] -1 )* 0.25)
    perc_75 <- ceiling((N[i] - 1) * 0.75) + 1
    perc_less_25_higher_75 <- sum(dbinom(c(0:perc_25, perc_75:N[i]), N[i], p))
    d <- data_frame(perc_less_25_higher_75 = perc_less_25_higher_75, size_class = N[i])
    data_4c <- rbind(data_4c, d) 
}
data_4c
p4c <- ggplot(data_4c, aes(x = factor(size_class), y = perc_less_25_higher_75, fill = factor(size_class)), log10 = "y") +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(perc_less_25_higher_75, 6)), 
              position = position_dodge(width = 0.9),    # Adjust position to align with bars
              vjust = -0.2,   # Adjust vertical position to be just above the bars
              size = 6) +
    scale_fill_manual(values = c("#008080", "#FF7F50", "#FFC300")) +  
    coord_trans(y= "sqrt")+  # Log scale for y-axis with original values as labels
    labs(title = "% of students with <25% or >75% female teachers by number of attended classes", 
         x = "Class size", y = "% of students with >25% or <75% female teachers",
         fill = "Class size") +
    theme_bw()  
p4c


