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
y_hat2 <- beta_hat2_intc + beta_hat2_x1 * x_2 + beta_hat2_x2 * x_2_bar
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

