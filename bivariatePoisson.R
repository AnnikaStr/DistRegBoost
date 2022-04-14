

biPoisLambda1 <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL){
  
  y_check <- function(y) {
    if ((is.matrix(y) && NCOL(y)!=2))
      stop("response should be a two-column matrix (y1 and y2) for this bivariate family")
  }
  
  # neg log-likelihood
  loss <- function(lambda2, lambda3, y, f=f){
    y2 <- y[,2]
    y1 <- y[,1]
    
    if(length(lambda2) ==1) lambda2 = rep(lambda2,length(y1))
    if(length(lambda3) ==1) lambda3 = rep(lambda3,length(y1))
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (lambda3[j] / (lambda2[j]*exp(f[j])))^a}))}) 
    
    loss_vec <- sapply(1:length(y1), function(i) exp(f[i]) + lambda2[i] + lambda3[i] - y1[i] * f[i] + lfactorial(y1[i]) - y2[i] * log(lambda2[i]) + lfactorial(y2[i]) - log(B[i]))
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(lambda2 = lambda2, lambda3 = lambda3, y = y, f = f))
  }
  
  # ngradient is the negative derivate
  ngradient <- function(y, f, w = 1){
    y2 <- y[,2]
    y1 <- y[,1]
    
    if(length(lambda2) ==1) lambda2 = rep(lambda2,length(y1))
    if(length(lambda3) ==1) lambda3 = rep(lambda3,length(y1))
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (lambda3[j] / (lambda2[j]*exp(f[j])))^a}))})
    B1 <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (lambda3[j] / lambda2[j] )^a  * (-a) * exp(-a*f[j])}))})
    
    sapply(1:length(y1), function(i){- exp(f[i]) + y1[i] + B1[i]/B[i]} )
  }
  
  offset <- function(y, w) {
    if (!is.null(lambda1)) {
      RET <- log(lambda1)
    }
    else {
      lambda3 = 1
      RET <- log(max( 0.1, weighted.mean(y[,1], w = w, na.rm = TRUE) - lambda3 ))  # log(weighted.mean(y[,1], w = w, na.rm = TRUE)- cov(y[,1],y[,2]))
      
    }
    return(RET)
  }
  

  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), offset = offset,
                  name = "bivariate Poisson-distribution: lambda1 (log link)")
}




biPoisLambda2 <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL){
  
  
  ## neg log-likelihood
  loss <- function(lambda1, lambda3, y, f=f){
    y2 <- y[,2]
    y1 <- y[,1]
    
    if(length(lambda1) ==1) lambda1 = rep(lambda1,length(y1))
    if(length(lambda3) ==1) lambda3 = rep(lambda3,length(y1))
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (lambda3[j] / (lambda1[j]*exp(f[j])))^a}))}) 
    
    sapply(1:length(y1), function(i) exp(f[i]) + lambda1[i] + lambda3[i] - y1[i] * log(lambda1[i])  + lfactorial(y1[i]) - y2[i] * f[i] + lfactorial(y2[i]) - log(B[i]))
    
    
  }
  
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(lambda1 = lambda1, lambda3 = lambda3, y = y, f = f))
  }
  
  # ngradient is the negative derivate
  ngradient <- function(y, f, w = 1){
    y2 <- y[,2]
    y1 <- y[,1]
    
    if(length(lambda1) ==1) lambda1 = rep(lambda1,length(y1))
    if(length(lambda3) ==1) lambda3 = rep(lambda3,length(y1))
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (lambda3[j] / (lambda1[j]*exp(f[j])))^a}))}) 
    B1 <- sapply(1:length(x), function(j){sum(sapply(0:x[j],function(a){choose(y1[j],a) * choose(y2[j],a) * factorial(a) *  (lambda3[j]/lambda1[j])^a * (-a) * exp(-a*f[j])}))})
    
    sapply(1:length(y1), function(i){- exp(f[i]) + y2[i] + B1[i]/B[i]})
    
  }
  

  offset <- function(y, w) {
    if (!is.null(lambda2)) {
      RET <- log(lambda2)
    }
    else {
      lambda3 = 1
      RET <- log(max( 0.1, weighted.mean(y[,2], w = w, na.rm = TRUE) - lambda3 ))
      
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), offset = offset,
                  name = "bivariate Poisson-distribution: lambda2 (log link)")
}




biPoisLambda3 <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL){
  
  
  loss <- function(lambda1, lambda2, y, f=f){
    y2 <- y[,2]
    y1 <- y[,1]
    
    if(length(lambda1) ==1) lambda1 = rep(lambda1,length(y1))
    if(length(lambda2) ==1) lambda2 = rep(lambda2,length(y1))
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (exp(f[j]) / (lambda1[j]*lambda2[j]))^a}))}) 
    
    sapply(1:length(y1), function(i) exp(f[i]) + lambda1[i] + lambda2[i] - y1[i] * log(lambda1[i]) + lfactorial(y1[i]) - y2[i] * log(lambda2[i]) + lfactorial(y2[i]) - log(B[i]))
    
  }
  
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(lambda1 = lambda1, lambda2 = lambda2, y = y, f = f))
  }
  
  # ngradient is the negative derivate w.r.t. lambda3
  ngradient <- function(y, f, w = 1){
    y2 <- y[,2]
    y1 <- y[,1]
    
    if(length(lambda1) ==1) lambda1 = rep(lambda1,length(y1))
    if(length(lambda2) ==1) lambda2 = rep(lambda2,length(y1))
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (exp(f[j]) / (lambda1[j]*lambda2[j]))^a}))})
    B1 <- sapply(1:length(x), function(j){sum(sapply(0:x[j],function(a){choose(y1[j],a) * choose(y2[j],a) * factorial(a) *(a * exp(a*f[j]))  *  (lambda1[j]*lambda2[j])^(-a)}))})
    
    
    sapply(1:length(y1), function(i){- exp(f[i]) + B1[i]/B[i]})
    
  }
  
  
  offset <- function(y, w) {
    if (!is.null(lambda3)) {
      RET <- log(lambda3)
    }
    else {
      RET = log(1)
      
    }
    return(RET)
  }
  
 
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), offset = offset,
                  name = "bivariate Poisson-distribution: lambda3 (log link)")
}



## families object for new distribution
PoissonBV <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL){
  Families(lambda1 = biPoisLambda1(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3),
           lambda2 = biPoisLambda2(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3),
           lambda3 = biPoisLambda3(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3),        
           name = "PoissonBV")
}




