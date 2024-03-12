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
library(ggplot2)


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
# (a) prepare data
data <- read_xlsx(file.path("xlsx", "FRED-QD.xlsx"))
data <- data %>% 
    mutate(y = (pnfix/lag(pnfix,n = 1) - 1)) %>% # Compute growth rate
    mutate(lag_1_y = lag(y, n = 1)) %>% # Compute lags of y
    mutate(lag_2_y = lag(lag_1_y, n = 1)) %>% 
    mutate(lag_3_y = lag(lag_2_y, n = 1)) %>% 
    mutate(lag_4_y = lag(lag_3_y, n = 1)) 
data <- na.omit(data) 

# Estimate model via OLS
ols <- lm(y ~ lag_1_y + lag_2_y + lag_3_y + lag_4_y, data)
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
se_5 <- sqrt(diag(v.nw.5))
# Compare results
se.hc1

cat("(c) Newey-West standard erros , with M = 5: ", se_5)

# Get coefficients from OLS
alpha <- coef(ols)
# Manually set IRF_k = 0 for k < 1 and IRF_1 = 1
IRF <- c(0,0,0,1,rep(NA,10))

for (j in 5:14) {
    IRF[j] <- alpha[2]*IRF[j-1] + alpha[3]*IRF[j-2] + alpha[4]*IRF[j-3] + alpha[5]*IRF[j-4]
}

IRF <- IRF[5:14]
j <- seq(1:10)
data_plot <- data.frame(IRF, j) 

# Plot results

ggplot(data = data_plot, aes(y = IRF, x = j)) +
    geom_line() +
    theme_bw() + 
    theme(
        # Adjusting size of axis text (tick labels)
        axis.text.x = element_text(size = 16), # Change 12 to your preferred size for x-axis text
        axis.text.y = element_text(size = 16), # Change 12 to your preferred size for y-axis text
        # If you also want to adjust title sizes, you can uncomment the following:
        axis.title.x = element_text(size = 16), # Change 14 to your preferred size for x-axis title
        axis.title.y = element_text(size = 16) # Change 14 to your preferred size for y-axis title
    )



