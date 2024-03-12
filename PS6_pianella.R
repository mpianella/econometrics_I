## ---------------------------
##
## Script name: PS6_pianella.R
##
## Purpose of script: solution to the sixth problem set in econometrics I 
##
## Author: Matteo Pianella
##
## Date Created: 12th March 2024

library(readxl)
library(AER)
library(dplyr)


# exercise 3 ------
# Prepare data
data <- read_xlsx(file.path( "xlsx", "AJR2001.xlsx"))
data$logm2 <-  data$logmort0^2
x <- as.matrix(cbind(1,data$risk))
z <- as.matrix(cbind(1,data$logmort0,data$logm2))
y <- as.matrix(data$loggdp)

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
cat("(a) Efficient GMM estimate: ", beta_gmm)

# Calculate J statistic
g.bar <- 1/n*(t(z)%*%y - t(z)%*%x%*%beta_gmm)
J <- n*t(g.bar)%*%w.hat%*%g.bar
cat("(b) J statistic for overindentification: ", J)

# 2SLS
two_sls_coeff <- iv.a$coefficients
cat("(c) 2SLS estimate: ", two_sls_coeff, "\n",
     "beta_gmm: ", beta_gmm)
beta_gmm

# exercise 6 -----
# Prepare data
data <- read_xlsx(file.path("xlsx", "FRED-QD.xlsx"))
data <- data %>% 
    mutate(y = (pnfix/lag(pnfix,n = 1) - 1)) %>% # Compute growth rate
    mutate(L1.y = lag(y, n = 1)) %>% # Compute kth lag of y
    mutate(L2.y = lag(L1.y, n = 1)) %>% 
    mutate(L3.y = lag(L2.y, n = 1)) %>% 
    mutate(L4.y = lag(L3.y, n = 1)) 
data <- na.omit(data) # Drop missing values (e.g. L1 is missing in t = 1)

# Estimate model via OLS
ols <- lm(y ~ L1.y + L2.y + L3.y + L4.y, data)
# Adj. heteroskedasticity robust variance estimator
v.hc1 <- vcovHC(ols, type="HC1")
se.hc1 <- sqrt(diag(v.hc1))
coeftest(ols, vcov = v.hc1)

# Newey-West variance estimator with M = 0
v.nw.0 <- NeweyWest(ols, lag = 0, prewhite = F, adjust = T)
se.hac.0 <- sqrt(diag(v.nw.0))

# This should equal to HC1 Se's
se.hc1

se.hac.0

# Newey-West variance estimator with M = 5
v.nw.5 <- NeweyWest(ols, lag = 5, prewhite = F, adjust = T)
se.hac.5 <- sqrt(diag(v.nw.5))
# Compare results
se.hc1

se.hac.5

# Get coefficients from OLS
alpha <- coef(ols)
# Manually set IRF_k = 0 for k < 1 and IRF_1 = 1
IRF <- c(0,0,0,1,rep(NA,10))

for (j in 5:14) {
    IRF[j] <- alpha[2]*IRF[j-1] + alpha[3]*IRF[j-2] + alpha[4]*IRF[j-3] + alpha[5]*IRF[j-4]
}
IRF <- IRF[5:14]
# Plot results
j <- seq(1:10)
plot(j,IRF,type='b')



