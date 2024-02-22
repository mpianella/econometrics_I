## ---------------------------
##
## Script name: PS4_pianella.R
##
## Purpose of script:
##
## Author: Matteo Pianella
##
## Date Created: 22nd Feb 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(readxl)
library(cols)
library(dplyr)
library(distrMod)
library(gmm)
library(car)
library(boot)

# please set repository on the econometrics_I folder

data <- read_xlsx(file.path("xlsx", "Nerlove1963.xlsx")) 

# estimate the Cobb-Douglas
cd_model <- lm(log(Cost) ~ log(output) + log(Plabor) + log(Pcapital) + log(Pfuel), data = data)

# report coefficients and standard error
summary_cd <- summary(cd_model)
print(summary_cd)

# estimate constrained least square
y <- data$Cost
x <- data %>% select(-Cost) %>% as.matrix()
R <- c(0,1,1,1)
cd_model_2 <- cls(log(y), log(x), R, 1)
summary_cd_2 <- summary(cd_model_2)
print(summary_cd_2)

# estimate with efficient minimum distance

# Wald test for the above restriction

linearHypothesis(cd_model, "log(Plabor) + log(Pcapital) + log(Pfuel) =  1",test = "Chisq")

# ex.5 -------
data1 <- read_xlsx(file.path("xlsx", "cps09mar.xlsx"))

data_model <- data1 %>% filter(race == 2)
model_1 <- lm(log(earnings) ~ as.factor(education), data = data_model) 
model_1_summ <- summary(model_1)
print(model_1_summ)

# ex.8 -------
data <- read_xlsx(file.path("xlsx", "Nerlove1963.xlsx"))

# estimating standard error by asymptotic
model <- lm(cd_model <- lm(log(Cost) ~ log(output) + log(Plabor) + log(Pcapital) + 
                               log(Pfuel), data = data))
vcov_matrix <- vcov(model)
asymptotic_se <- sqrt(sum(vcov_matrix[c("log(Plabor)", "log(Pcapital)", "log(Pfuel)"), c("log(Plabor)", "log(Pcapital)", "log(Pfuel)")]))  # Standard errors of the coefficients
print(asymptotic_se)

# estimating standard error by bootstrapping 
# Define the function to obtain the regression coefficients
theta_fn <- function(data, indices) {
    d <- data[indices, ]  # Resample the data
    model <- lm(log(Cost) ~ log(output) + log(Plabor) + log(Pcapital) + 
                  log(Pfuel), data = d)
    return(sum(coef(model)[c("log(Plabor)", "log(Pcapital)", "log(Pfuel)")]))  # Return the coefficients
}

# Apply the bootstrap
set.seed(31415)  
boot_results <- boot(data = data, statistic = theta_fn, R = 1000)
bootstrap_se <- sd(boot_results$t)
print(bootstrap_se)

# Calculate CI using the percentile method
percentile_ci <- boot.ci(boot_results, type = "perc", conf = 0.95)

# Calculate CI using the BCa method
bca_ci <- boot.ci(boot_results, type = "bca", conf = 0.95)

# Print the CIs
print(percentile_ci)
print(bca_ci)

# estimation using jackknife
jackknife_theta <- numeric(nrow(data))
n <- nrow(data)
for (i in 1:n) {
    model <- lm(log(Cost) ~ log(output) + log(Plabor) + log(Pcapital) + 
                    log(Pfuel), data = data[-i, ])
    jackknife_theta[i] <- sum(coef(model)[c("log(Plabor)", "log(Pcapital)", "log(Pfuel)")])
}

mean_theta <- mean(jackknife_theta)
jackknife_se <- sqrt((n-1)/n * sum((jackknife_theta - mean_theta)^2))
print(jackknife_se)

# ex.9 ----
data1 <- read_xlsx(file.path("xlsx", "cps09mar.xlsx")) %>% 
    filter(hisp ==1, female ==0, region == 2, race == 1, marital == 7) 

# create a variable for experience
data1 <- data1 %>% mutate(experience = age - education - 6)
model <- lm(log(earnings) ~ education + experience + I((experience^2)/100), data = data1)

# Extract coefficients
beta1 <- summary(model)$coefficients['education', 'Estimate']
beta2 <- summary(model)$coefficients['experience', 'Estimate']
beta3 <- summary(model)$coefficients['I((experience^2)/100)', 'Estimate']

# Calculate theta for experience = 10
theta <- beta1 / (beta2 + 2 * beta3 * 0.1)

## Calculate the asymptotic standard error using the delta method
# Calculate the gradient of theta w.r.t. beta1, beta2, and beta3
gradient <- c(1 / (beta2 + 2 * beta3 * 0.1), 
              -beta1 / (beta2 + 2 * beta3 * 0.1)^2,
              -0.2 * beta1 / (beta2 + 2 * beta3 * 0.1)^2)

# Covariance matrix of the coefficients
cov_matrix <- vcov(model)

# Subset the covariance matrix to include only beta1, beta2, and beta3
cov_matrix_subset <- cov_matrix[c("education", "experience", "I((experience^2)/100)"), 
                                c("education", "experience", "I((experience^2)/100)")]

# Calculate the variance of theta using the delta method
var_theta <- t(gradient) %*% cov_matrix_subset %*% gradient

# The standard error is the square root of the variance
se_theta_asymptotic <- sqrt(var_theta)

cat("Asymptotic Standard Error of Theta: ", se_theta_asymptotic, "\n")

## Calculate the standard error using jackknife
theta_jackknife <- numeric(nrow(data1))

for (i in 1:n) {
    jackknife_sample <- data1[-i,]  # Exclude the i-th observation
    jackknife_model <- lm(log(earnings) ~ education + experience + I((experience^2)/100), data = jackknife_sample)
    
    beta1_j <- summary(jackknife_model)$coefficients['education', 'Estimate']
    beta2_j <- summary(jackknife_model)$coefficients['experience', 'Estimate']
    beta3_j <- summary(jackknife_model)$coefficients['I((experience^2)/100)', 'Estimate']
    
    theta_jackknife[i] <- beta1_j / (beta2_j + 2 * beta3_j * 0.1)
}

# Calculate the jackknife estimate and standard error
theta_mean <- mean(theta_jackknife)
se_theta_jackknife <- sqrt((n-1)/n * sum((theta_jackknife - theta_mean)^2))

cat("Jackknife Standard Error of Theta: ", se_theta_jackknife, "\n")

## Bootstrap to calculate standard error
set.seed(3141) 
n <- nrow(data1)
n_bootstraps <- 1000
theta_bootstraps <- numeric(n_bootstraps)

for(i in 1:n_bootstraps) {
    sample_indices <- sample(1:n, replace = TRUE)
    bootstrap_sample <- data1[sample_indices, ]
    bootstrap_model <- lm(log(earnings) ~ education + experience + I((experience^2)/100), data = bootstrap_sample)
    
    beta1_b <- summary(bootstrap_model)$coefficients['education', 'Estimate']
    beta2_b <- summary(bootstrap_model)$coefficients['experience', 'Estimate']
    beta3_b <- summary(bootstrap_model)$coefficients['I((experience^2)/100)', 'Estimate']
    
    theta_b <- beta1_b / (beta2_b + 2 * beta3_b * 0.1)
    theta_bootstraps[i] <- theta_b
}

# Calculate standard error of theta from bootstrap
se_theta_bootstrap <- sd(theta_bootstraps)

# Print results
cat("Estimated Theta at Experience = 10: ", theta, "\n")
cat("Bootstrap Standard Error of Theta: ", se_theta_bootstrap, "\n")

## Calculate confidence intervals using the BS percentile method
# Calculate the bias-correction factor z0
p <- mean(theta_bootstraps < theta)
z0 <- qnorm(p)

# Calculate adjusted percentiles for the confidence interval
alpha1 <- pnorm(2 * z0 + qnorm(0.025))
alpha2 <- pnorm(2 * z0 + qnorm(0.975))

# Find the corresponding percentiles in the bootstrap distribution
ci_low <- quantile(theta_bootstraps, alpha1)
ci_high <- quantile(theta_bootstraps, alpha2)

cat("BCa 95% Confidence Interval for Theta: [", ci_low, ", ", ci_high, "]\n")
