## ---------------------------
##
## Script name: PS6_pianella.R
##
## Purpose of script: solution to the sixth problem set in econometrics I 
##
## Author: Matteo Pianella
##
## Date Created: 12th March 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# Prepare data
data <- read_dta("AJR2001.dta")
data$logm2 <- data$logmort0ˆ2
x <- as.matrix(cbind(1,data$risk))
z <- as.matrix(cbind(1,data$logmort0,data$logm2))
y <- as.matrix(data$loggdp)

## a) -------
# Use 2SLS to get fitted values
iv.a <- ivreg(loggdp ~ risk| logmort0 + logm2, data=data)
# Initial matrices
Omega.hat <- matrix(0, ncol = 3, nrow = 3)
Q.hat <- matrix(0, ncol = 2, nrow = 3)
n <- nrow(y)
for (i in 1:n) {
    e.tilde <- iv.a$residuals[i]
    g.tilde <- as.matrix(z[i,] * e.tilde)
    Omega.hat <- Omega.hat + 1/n * g.tilde %*% t(g.tilde)
    Q.hat <- Q.hat + 1/n * z[i,] %*% t(x[i,])
}
w.hat <- solve(Omega.hat)
beta_gmm <- solve(t(x)%*%z%*%w.hat%*%t(z)%*%x) %*% t(x)%*%z%*%w.hat%*%t(z)%*%y

## b) -------
# Calculate J statistic
g.bar <- 1/n*(t(z)%*%y - t(z)%*%x%*%beta_gmm)
J <- n*t(g.bar)%*%w.hat%*%g.bar

## c) -------
# 2SLS
iv.a$coefficients

beta_gmm