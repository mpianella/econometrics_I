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

# exercise 4(a) -------
## I run a simulation for 1000 cities using a binomial distribution

s_bar <- .556
n_cities_below_6000 <- 102
set.seed(314)
n_below_6000 <- rbinom(1000, n_cities_below_6000, s_bar)
s_c <- n_below_6000/102 # calculate s_c for each city 
sd_s_c <- sqrt((sum((s_c - s_bar)^2))/n_cities_below_6000)
print(sd_s_c)


