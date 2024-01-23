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
## I run a simulation of 819 cities using a binomial distribution

s_bar <- .556
n_cities_between5500_and_6500 <- 102
set.seed(314)
n_sim <- 819
n_below_6000 <- rbinom(n_sim, n_cities_between5500_and_6500, s_bar)
s_c <- n_below_6000/102 # calculate s_c for each city 
sd_s_c <- sqrt((sum((s_c - s_bar)^2))/n_sim)
print(sd_s_c) # my sd is consistently lower than what they have in Table A1 (0.118)




