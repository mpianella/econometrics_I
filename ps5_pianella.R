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
library(AER)
library(stargazer)
library(sandwich)
library(lmtest)

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
beta1_unconcentrated <- model1$par[2]
beta2_unconcentrated <- model1$par[3]
beta3_unconcentrated <- model1$par[4]
beta4_unconcentrated <- model1$par[5]

# concentrated search for gamma
lower_bound <- c(gamma_range[1], -Inf, -Inf, -Inf, -Inf)
upper_bound <- c(gamma_range[2], Inf, Inf, Inf, Inf)

model2 <- optim(par = c(mean(gamma_range), 0,0,0,0), fn = expected_sqrt_error, data = data, 
      method = "L-BFGS-B", lower = lower_bound, upper = upper_bound)
gamma_concentrated <- model2$par[1]
beta1_concentrated <- model2$par[2]
beta2_concentrated <- model2$par[3]
beta3_concentrated <- model2$par[4]
beta4_concentrated <- model2$par[5]



# calculate standard error for each parameter
log_PL <- data$log_PL
log_PF <- data$log_PF
log_PK <- data$log_PK
log_Q  <- data$log_Q
log_TC <- data$log_TC
vec1 <- c(rep(1, nrow(data)))
vec2 <- c(log_Q)
vec3 <- c(log_PL + log_PF + log_PK)
vec4_unconcentrated <- (log_Q)/(1+exp(-(log_Q - gamma_unconcentrated)))
vec4_concentrated  <- (log_Q)/(1+exp(-(log_Q - gamma_concentrated)))
vec5_unconcentrated <- (-log_Q * (exp(-(log_Q - gamma_unconcentrated))) * gamma_unconcentrated)/(1 + exp(-(log_Q - gamma_unconcentrated)))^2
vec5_concentrated <- (-log_Q * (exp(-(log_Q - gamma_concentrated))) * gamma_concentrated)/(1 + exp(-(log_Q - gamma_concentrated)))^2

m_hat_concentr <- cbind(vec1, vec2, vec3, vec4_concentrated, vec5_concentrated)
m_hat_unconcentr <- cbind(vec1, vec2, vec3, vec4_unconcentrated, vec5_unconcentrated)


# compute Q_hat matrix 

# Iterate over each row of M to calculate the product and sum the results

get_Q_hat <- function(m_theta){
    sum_matrix <- matrix(0, 5, 5) # 5x5 matrix filled with 0s
    # Iterate over each row of M to calculate the product and sum the results
    for (i in 1:nrow(m_theta)) {
        row_vector <- matrix(m_theta[i, ], nrow = 1) # Extract the row as a 1x5 matrix
        product <- t(row_vector) %*% row_vector # Calculate the product
        sum_matrix <- sum_matrix + product # Sum the product matrices
    }
    return(sum_matrix/nrow(m_theta))
}
Q_hat_concentrated   <- get_Q_hat(m_hat_concentr)
Q_hat_unconcentrated <- get_Q_hat(m_hat_unconcentr)

# compute the residual for each regression
e_hat_unconcentr <- log_TC - beta1_unconcentrated + beta2_unconcentrated * log_Q + beta3_unconcentrated * (log_PL + log_PK + log_PF) + beta4_unconcentrated * log_Q/(1 + exp(-(log_Q - gamma_unconcentrated)))
e_hat_concentr <- log_TC - beta1_concentrated + beta2_concentrated * log_Q + beta3_concentrated * (log_PL + log_PK + log_PF) + beta4_concentrated * log_Q/(1 + exp(-(log_Q - gamma_concentrated)))

# compute omega_hat
get_Omega_hat <- function(m_theta, e_hat){
    sum_matrix <- matrix(0, 5, 5) # 5x5 matrix filled with 0s
    # Iterate over each row of M to calculate the product and sum the results
    for (i in 1:nrow(m_theta)) {
        row_vector <- matrix(m_theta[i, ], nrow = 1) # Extract the row as a 1x5 matrix
        product <- (t(row_vector) %*% row_vector) * (e_hat[i]^2) # Calculate the product
        sum_matrix <- sum_matrix + product # Sum the product matrices
    }
    return(sum_matrix/nrow(m_theta))
}

# compute the Omega_hat for each regression
Omega_hat_unconcentr <- get_Omega_hat(m_hat_unconcentr, e_hat = e_hat_unconcentr)
Omega_hat_concentr <- get_Omega_hat(m_hat_concentr, e_hat = e_hat_concentr)

# compute the covariance matrix 
V_hat_unconcentrated <- solve(Q_hat_unconcentrated) %*% Omega_hat_unconcentr %*% solve(Q_hat_unconcentrated)
V_hat_concentrated <- solve(Q_hat_concentrated) %*% Omega_hat_concentr %*% solve(Q_hat_concentrated)

# standard errors 
standard_error_unconc <- sqrt(diag(V_hat_unconcentrated))
standard_error_conc <- sqrt(diag(V_hat_concentrated))

# ex.9 -----
d_AJR2001 <- read_xlsx(file.path("xlsx", "AJR2001.xlsx"))

# data preparation

data <- d_AJR2001

# OLS estimate for regression 12.86
model_1 <- lm(loggdp ~ risk, data = data)
summary(model_1)
stargazer(model_1, title = "Table 1", type = "text", out = "./txt/ps5_table1.txt")

# First stage regresssion (eq 12.87)
first_stage <- lm(risk ~ logmort0, data = data)
data$risk_hat <- predict(first_stage)
summary(first_stage)
stargazer(first_stage, title = "Table 2", type = "text", out = "./txt/ps5_table2.txt")


# Second stage regression (eq. 12.88)
second_stage <- ivreg(loggdp ~ risk_hat | logmort0, data = data)
summary(second_stage)
stargazer(second_stage, title = "Table 3", type = "text", out = "./txt/ps5_table3.txt")


# calculate heteroskedastic-robust standard errors
coeftest(model_1, vcov = vcovHC(model_1, type = "HC1"))
coeftest(first_stage, vcov = vcovHC(first_stage, type = "HC1"))
coeftest(second_stage, vcov = vcovHC(second_stage, type = "HC1"))

## calculating the regression using Indirect Least Square
# calculate gamma hat
model_indirect_gamma <- lm(risk ~ logmort0, data = data)
gamma_hat <- model_indirect_gamma$coefficients[2]
# calculate lambda hat
model_indirect_lambda <- lm(loggdp ~ logmort0, data = data)
lambda_hat <- model_indirect_lambda$coefficients[2]

# calculate ILS estimator
ils_estimator <- (1/gamma_hat) * lambda_hat
print(ils_estimator)


