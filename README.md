# Boosting multivariate distributional regression models 
Component-wise gradient boosting algorithm for modeling multivariate distributional regression models within the framework of generalized additive models for location, scale and shape, which enables the simultaneous modeling of all distribution parameters of a parametric distribution of a multivariate response conditional on
explanatory variables.


# Example 
```
require('mvtnorm')
require('gamboostLSS')

source('bivariateGaussian.R') # load family

# data generation
n = 5000
set.seed(1)

x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- runif(n)
x5 <- runif(n)
x6 <- runif(n)
x <- data.frame(x1,x2,x3,x4,x5,x6)

mu1 <- -2 * x[,1] + 1 * x[,2] - 1 * x[,3] + 2 * x[,4] 
mu2 <- -1 * x[,4] + 1 * x[,5] - 1.5* x[,6]

sigma1 <- exp( 1 + -0.5 * x[,1] + 0.25 * x[,3] )
sigma2 <- exp( 1 + 0.5 * x[,2] - 0.25 *x[,4])

eta_rho <- -1 * x[,1] + 1* x[,3] 
rho <- (eta_rho) / sqrt(1+ eta_rho^2) 

y <- NULL
for(i in 1:n) {
  y <- rbind(y,rmvnorm(1, mean = c(mu1[i],mu2[i]), sigma = matrix(c(sigma1[i]^2, rho[i]*sigma1[i]*sigma2[i], rho[i]*sigma1[i]*sigma2[i], sigma2[i]^2), 2, 2)))
}

colnames(y) <- c('y1', 'y2')
dat <- cbind(y,x)

plot(dat$y1, dat$y2)

# --- model
model1 <- glmboostLSS(formula = cbind(y1,y2)~., data = dat, control = boost_control(trace = FALSE, mstop = 10, nu = 0.1), 
                      families = GaussianBV(mu1 = NULL, mu2 = NULL, sigma1 = NULL , sigma2 = NULL,  rho = NULL), method = 'noncyclic')

coef(model1[500])

```
