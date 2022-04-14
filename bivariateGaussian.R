

GaussianBVMu1 <- function (mu1 = NULL, sigma1 = NULL ,mu2 = NULL, sigma2 = NULL, rho = NULL) 
{
  
  
  y_check <- function(y) {
    if ((is.matrix(y) && NCOL(y)!=2))
      stop("response should be a two-column matrix (y1 and y2) for this bivariate family")
  }
  

  
  loss <- function(mu2 = mu2, sigma2, sigma1, y = y, f = f, rho){ 
    y2 <- y[,2]
    y1 <- y[,1]
    
    
    -(-log(2*pi) - log(sigma1) - log(sigma2) - .5*log(1 - rho^2) - 1/(2*(1 - rho^2)) *((y1 - f)^2/sigma1^2  +  (y2 - mu2)^2/sigma2^2  
                                                                                       - 2*rho*((y1 - f)*(y2 - mu2)/(sigma1*sigma2))))
  } 
  
  ngradient <- function(y, f, w = 1, scale.grad = scale.grad) {
    y2 <- y[,2]
    y1 <- y[,1]
    ngr <- 1/((1 - rho^2)*sigma1) * ( (y1 - f)/sigma1 - rho * (y2 - mu2)/sigma2 )

    return(ngr)
  }
  
  
  
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma1 = sigma1, sigma2 = sigma2,  mu2 = mu2, rho = rho))
  }
  
  
  
  offset <- function(y, w) {
    if (!is.null(mu1)) {
      RET <- mu1
    }
    else {
      RET <- weighted.mean(y[,1], w = w, na.rm = TRUE)
    }
    return(RET)
  }
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) f, 
                 offset = offset, name = "Bivariate Normal distribution: mu1 (id link)")
}


GaussianBVMu2 <- function (mu1 = NULL, sigma1 = NULL ,mu2 = NULL, sigma2 = NULL, rho = NULL) 
{
  loss <- function(mu1 = mu1, sigma2, sigma1, y, f = f, rho){ 
    y1 <- y[,1]
    y2 <- y[,2]
    
    -(-log(2*pi) - log(sigma1) - log(sigma2) - .5*log(1 - rho^2) - 1/(2*(1 - rho^2)) *((y1 - mu1)^2/sigma1^2  +  (y2 - f)^2/sigma2^2  
                                                                                       - 2*rho*((y1 - mu1)*(y2 - f)/(sigma1*sigma2))))
    
  } 
  
  ngradient <- function(y, f, w = 1) {
    y1 <- y[,1]
    y2 <- y[,2]
    ngr <-  1/((1 - rho^2)*sigma2) * ( (y2 - f)/sigma2 - rho * (y1 - mu1)/sigma1 )

    return(ngr)
  }
  
  
  
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma1 = sigma1, sigma2 = sigma2,  mu1 = mu1, rho = rho))
  }
  
  
  
  offset <- function(y, w) {
    if (!is.null(mu2)) {
      RET <- mu2
    }
    else {
      RET <- weighted.mean(y[,2], w = w, na.rm = TRUE)
    }
    return(RET)
  }
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) f, 
                 offset = offset, name = "Bivariate Normal distribution: mu2 (id link)")
}


GaussianBVSigma2 <- function (mu1 = NULL, sigma1 = NULL ,mu2 = NULL, sigma2 = NULL, rho = NULL , scale.grad = TRUE) 
{

  
  loss <- function(mu2, sigma1,  mu1, y,  f = f, rho){ 
    y1 <- y[,1]
    y2 <- y[,2]
    
    -( -log(2*pi) - log(sigma1) - f - .5*log(1 - rho^2) - 1/(2*(1 - rho^2)) *((y1 - mu1)^2/sigma1^2  +  (y2 - mu2)^2/exp(2*f)  
                                                                              - 2*rho*((y1 - mu1)*(y2 - mu2)/(sigma1*exp(f)))))
  } 
  
  ngradient <- function(y, f, w = 1) {
    y1 <- y[,1]
    y2 <- y[,2]
    ngr <-   -1 + (1 / (1 - rho^2)) * ((y2 - mu2) / exp(f) ) * ( (y2 - mu2)/exp(f) - rho * (y1 - mu1)/sigma1 )

    return(ngr)
  }
  
  
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu1 = mu1, sigma1 = sigma1,  mu2 = mu2, rho = rho))
  }
  
  
  
  offset <- function(y, w) {
    if (!is.null(sigma2)) {
      RET <- log(sigma2)
    }
    else {
      RET <- log(sd(y[w!=0,2],na.rm = TRUE))
    }
    return(RET)
  }
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), 
                 offset = offset, name = "Bivariate Normal distribution: sigma2 (log link)")
}


GaussianBVSigma1 <- function (mu1 = NULL, sigma1 = NULL ,mu2 = NULL, sigma2 = NULL, rho = NULL , scale.grad = TRUE) 
{

  loss <- function(mu2, sigma2,  mu1, y,  f = f, rho){ 
    y1 <- y[,1]
    y2 <- y[,2]
    
    -(-log(2*pi) - f - log(sigma2) - .5*log(1 - rho^2) - 1/(2*(1 - rho^2)) *((y1 - mu1)^2/exp(2*f)  +  (y2 - mu2)^2/sigma2^2  
                                                                             - 2*rho*((y1 - mu1)*(y2 - mu2)/(exp(f)*sigma2))))
  } 
  
  ngradient <- function(y, f, w = 1) {
    y1 <- y[,1]
    y2 <- y[,2]
    ngr <-  -1 + (1 / (1 - rho^2)) * ((y1 - mu1) / exp(f) ) * ( (y1 - mu1)/exp(f) - rho * (y2 - mu2)/sigma2 ) 

    return(ngr)
  }
  
  
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu1 = mu1, sigma2 = sigma2,  mu2 = mu2, rho = rho))
  }
  
  
  
  offset <- function(y, w) {
    if (!is.null(sigma1)) {
      RET <- log(sigma1)
    }
    else {
      RET <- log(sd(y[w!=0,1], na.rm = TRUE))
    }
    return(RET)
  }
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), 
                 offset = offset, name = "Bivariate Normal distribution: sigma1 (log link)")
}

GaussianBVRho <- function (mu1 = NULL, sigma1 = NULL ,mu2 = NULL, sigma2 = NULL, rho = NULL ,  scale.grad = TRUE) 
{

  loss <- function(mu2, sigma2,  mu1, sigma1,  y,  f = f){ 
    y1 <- y[,1]
    y2 <- y[,2]
    -(-log(2*pi) - log(sigma1) - log(sigma2) - .5*log(1 - (f / sqrt( 1 + f^2))^2) - 1/(2*(1 - (f / sqrt( 1 + f^2))^2)) *((y1 - mu1)^2/sigma1^2  +  (y2 - mu2)^2/sigma2^2  
                                                                                                                         - 2*(f / sqrt( 1 + f^2))*((y1 - mu1)*(y2 - mu2)/(sigma1*sigma2))))
  } 
  
  ngradient <- function(y, f, w = 1) {
    y1 <- y[,1]
    y2 <- y[,2]
    ngr <-   1/(1 - (f / sqrt( 1 + f^2))^2) * ( (f / sqrt( 1 + f^2)) / (1 + f^2)^(3/2)) - f*( (y1 - mu1)^2/sigma1^2  +  (y2 - mu2)^2/sigma2^2 ) +
      (1 + 2*f^2) / sqrt(1 + f^2)  * (y1 - mu1)/sigma1 * (y2 - mu2)/sigma2
    
    return(ngr)
  }
  
  
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu1 = mu1, sigma2 = sigma2, mu2 = mu2, sigma1 = sigma1))
  }
  
  
  
  offset <- function(y, w = 1) {
    if (!is.null(rho)) {
      RET <- (rho / sqrt(1 - rho^2))
    }
    else {
      rho = cor(y[w !=0,1], y[w!=0,2])
      RET <-  (rho / sqrt(1 - rho^2))
    }
    return(RET)
  }
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) (f / sqrt( 1 + f^2)), 
                 offset = offset, name = "Bivariate Normal distribution: rho (eta = (rho / sqrt(1 - rho^2)) )")
}


#------------------ complete gamboostLSS families

GaussianBV <- function (mu1 = NULL, sigma1 = NULL, mu2 = NULL, sigma2 = NULL, rho = NULL) 
{
  if ((!is.null(sigma1) && sigma1 <= 0)) 
    stop(sQuote("sigma1"), " must be greater than zero")
  if ((!is.null(sigma2) && sigma2 <= 0)) 
    stop(sQuote("sigma2"), " must be greater than zero")
  if ((!is.null(rho) && rho <= 0 && rho >= 1)) 
    stop(sQuote("rho"), " must be greater than 0 and smaller than 1")
  
  
  Families(   mu1 = GaussianBVMu1(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, rho = rho), 
              mu2 = GaussianBVMu2(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, rho = rho), 
              sigma1 = GaussianBVSigma1(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, rho = rho), 
              sigma2 = GaussianBVSigma2(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, rho = rho), 
              rho = GaussianBVRho(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, rho = rho),        
              name = "GaussianBV")
  
}
