## ---------------------------
##
## Script name: PS3_pianella.R
##
## Purpose of script: my solution to the third problem set for the PhD course 
##                    econometrics I 
##
## Author: Matteo Pianella
##
## Date Created: 13th February 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(ggplot2)
library(GGally)
library(MASS)


# ex.1.b --------

# set the parameters 
n_sim <- 500 
sigma_sqr <- 1
beta <- 5

# run the simulation and calculate the correlatin
set.seed(314)
u_0 <- rnorm(n = n_sim, mean = 0, sd = sigma_sqr)
u_1 <- rnorm(n = n_sim, mean = 0, sd = 1)
u_2 <- rnorm(n = n_sim, mean = 0, sd = 1)
x_1 <- u_0 + u_2
x_2 <- u_1 + u_2
epsilon <- rnorm(n = n_sim, mean = 0, sd = 1)
y <- x_1 + beta * x_2 + epsilon

print(cor(x_1, x_2))

# predict the value of y given x_1 and calculate mean squared error
d1 <- data.frame(y, x_1)
model1 <- lm(y ~ x_1, d1)
y_hat1 <- model1$fitted.values
mse1 <- sum((y - y_hat1)^2)/n_sim
print(mse1)

#estimate the second model and calculate the mean squared error
model2 <- lm(y ~ x_1 + x_2, d1)
beta_hat2_intc <- model2$coefficients["(Intercept)"]
beta_hat2_x1 <- model2$coefficients["x_1"]
beta_hat2_x2 <- model2$coefficients["x_2"]
x_2_bar <- mean(x_2)
y_hat2 <- beta_hat2_intc + beta_hat2_x1 * x_1 + beta_hat2_x2 * x_2_bar
mse2 <- sum((y - y_hat2)^2)/n_sim
print(mse2)

# estimate x_2 from x_1 
d2 <- data.frame(x_1, x_2)
model_x_2 <- lm(x_2 ~ x_1, data = d2)
x_2_hat <- model_x_2$fitted.values
y_hat3 <- beta_hat2_intc + beta_hat2_x1 * x_1 + beta_hat2_x2 * x_2_hat
mse3 <- sum((y - y_hat3)^2)/n_sim 
print(mse3)

d_plot <- data.frame(y_hat1, y_hat2, y_hat3)
ggpairs(d_plot)

# generate fitted valus for x_2 when fitting polinomial of x_1
d3 <- data.frame(x_1, x_2)
model_x_21 <- lm(x_2 ~ poly(x_1, 20), d3)
x_3 <- model_x_21$fitted.values


# estimate y and calculate mse 
d4 <- data.frame(y, x_1, x_3)
model4 <- lm(y ~ x_1 + x_3, d4)
y_hat4 <- model4$fitted.values
mse4 <- sum((y - y_hat4)^2)/n_sim 
print(mse4)

d_plot1 <- data.frame(y_hat1, y_hat2, y_hat3, y_hat4)
ggpairs(d_plot1)

# # creating a proxi for x2 using info on correlation between x1 and x2 and the info on mean and sd
# r <- cor(x_1, x_2)
# sigma_x2 <- sd(x_2)
# mu_x2 <- mean(x_2)
# print(r)
# set.seed(31415)
# e <- rnorm(n_sim, mean = 0, sd = 1)
# x_3_no_mean <- r * x_1 + sqrt(1 - r^2) * e
# z_x_3 <- (x_3 - mean(x_3))/sd(x_3) #standardize x_3
# x_3 <- sigma_x2 * z_x_3 + mu_x2 # transform to same mean and sd as x_2
# print(cor(x_1, x_3))


# # generate estimates for y using x_1 and the proxi for x_2 and calculate mse
# d3 <- data.frame(y, x_1, x_3)
# model4 <- lm(y ~ x_1 + x_3, d3)
# y_hat4 <- model4$fitted.values
# mse4 <- sum((y - y_hat4)^2)/n_sim 
# print(mse4)
# d_plot1 <- data.frame(y_hat1, y_hat2, y_hat3, y_hat4)
# ggpairs(d_plot1)

# calculate mse when both variables are available
d4 <- data.frame(y, x_1, x_2)
model5 <- lm(y ~ x_1 + x_2, d4)
y_hat5 <- model5$fitted.values
mse5 <- sum((y - y_hat5)^2)/n_sim
print(mse5)

# ex.2 -----
n_sim <- 500
set.seed(3141)
earnings <- rnorm(n = n_sim, mean = 19, sd = 1)
capital_gains <- rnorm(n = n_sim, mean = 1, sd = 1)
u <- rnorm(n = n_sim, mean = 1, sd = 0)
e <- rnorm(n = n_sim, mean = 1, sd = 0)
occup_status <- earnings + u 
child_outcomes <- earnings - capital_gains + e
income <- earnings + capital_gains

# fraction of income that comes from labour marker earnings
print(mean(earnings/income))

# regress child outcomes on earnings and capital gains and show that controlling for 
# occupational status does not affect the other coefficients 
d1 <- data.frame(child_outcomes, earnings, capital_gains, occup_status, income)
model1 <- lm(child_outcomes ~ earnings + capital_gains, data = d1)
model2 <- lm(child_outcomes ~ earnings + capital_gains + occup_status, data = d1)
beta_hat_earnings <- model1$coefficients["earnings"]
beta_hat_capital_gains <- model1$coefficients["capital_gains"]
beta_hat_occup_status <- model2$coefficients["occup_status"]
print(c(beta_hat_earnings, beta_hat_capital_gains))
print(c(beta_hat_earnings, beta_hat_capital_gains, beta_hat_occup_status)) # NA for occup_status because of multicollinearity

# regress child outcomes on income and compare to previous model
model3 <- lm(child_outcomes ~ income, data = d1)
beta_hat_income <- model3$coefficients["income"]
print(beta_hat_income)

# correlation between occupational status and income 
print(cor(occup_status, income))

# regress child outcome on income controlling for occupational status
model4 <- lm(child_outcomes ~ income + occup_status, data = d1)
beta_hat_income1 <- model4$coefficients["income"]
beta_hat_occup_status1 <- model4$coefficients["occup_status"]
print(c(beta_hat_income1, beta_hat_occup_status1))
