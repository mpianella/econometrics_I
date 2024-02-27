## ---------------------------
##
## Script name: ps5_pianella.R
##
## Purpose of script: solution to the fifth problem set for econometrics I 
##
## Author: Matteo Pianella
##
## Date Created: 27th February
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(readxl)
library(dplyr)
library(GenSA)

# Please select the econometrics_I folder as your working directory

# ex.4 --------
# read data 
data <- read_xlsx(file.path("xlsx", "Nerlove1963.xlsx"))

# prepare data
data <- data %>% 
    mutate(log_Q = log(output), log_PL = log(Plabor), log_PK = log(Pcapital), 
           log_PF = log(Pfuel), log_TC = log(Cost)) 

# chose a range for gamma
gamma_range <- c(as.numeric(sort(data$log_Q)[15]), as.numeric(sort(data$log_Q, decreasing = TRUE)[15])) 


## estimating via global search 
expected_sqrt_error <- function(par, data){
    gamma <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    beta3 <- par[4]
    beta4 <- par[5]
    n <- nrow(data)
    
    log_Q <- data$log_Q
    log_PK <- data$log_PK
    log_PF <- data$log_PF
    log_PL <- data$log_PL
    log_TC <- data$log_TC
    
    x_3 <- log_PL + log_PK + log_PF
    x_4 <- log_Q/(1 + exp(-(log_Q - gamma)))
    return((1/n) * sum((log_TC - beta1 + beta2 * log_Q + beta3 * x_3 + beta4 * x_4)^2))
}

model1 <- optim(par = c(0,0,0,0,0), fn = expected_sqrt_error, data = data)
gamma_unconcentrated <- model1$par[1]

# concentrated search for gamma
lower_bound <- c(gamma_range[1], -Inf, -Inf, -Inf, -Inf)
upper_bound <- c(gamma_range[2], Inf, Inf, Inf, Inf)

model2 <- optim(par = c(mean(gamma_range), 0,0,0,0), fn = expected_sqrt_error, data = data, 
      method = "L-BFGS-B", lower = lower_bound, upper = upper_bound)
gamma_concentrated <- model2$par[1]


# calculate standard error for each parameter
log_PL <- data$log_PL
log_PF <- data$log_PF
log_PK <- data$log_PK
log_Q  <- data$log_Q
vec1 <- c(rep(1, nrow(data)))
vec2 <- c(log_Q)
vec3 <- c(log_PL + log_PF + log_PK)
vec4_unconcentrated <- (log_Q)/(1+exp(-(log_Q - gamma_unconcentrated)))
vec4_concentrated  <- (log_Q)/(1+exp(-(log_Q - gamma_concentrated)))
vec5_unconcentrated <- (-log_Q * (exp(-(log_Q - gamma_unconcentrated))) * gamma_unconcentrated)/(1 + exp(-(log_Q - gamma_unconcentrated)))
vec5_concentrated <- (-log_Q * (exp(-(log_Q - gamma_concentrated))) * gamma_concentrated)/(1 + exp(-(log_Q - gamma_concentrated)))

m_hat_concentr <- cbind(vec1, vec2, vec3, vec4_concentrated, vec5_concentrated)
m_hat_unconcentr <- cbind(vec1, vec2, vec3, vec4_unconcentrated, vec5_unconcentrated)


# compute Q_hat matrix 
get_Q_hat <- function(m_theta, n){
    dot_product <- apply(m_theta , 1, function(row) sum(row * row))
    
}
dot_product <- apply(m_hat_concentr , 1, function(row) sum(row * row))
Q_hat_concentrated   <- apply(m_hat_unconcentr, 1, function(row) sum(sum * row))