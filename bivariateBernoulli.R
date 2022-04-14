

biP1 <- function(mu1 = NULL, mu2 = NULL, or = NULL){
  
  y_check <- function(y) {
    if ((is.matrix(y) && NCOL(y)!=2))
      stop("response should be a two-column matrix (y1 and y2) for this bivariate family")
  }
  
  # neg. log-likelihood
  loss <- function(mu2, or, y, f=f){
    y1 <- y[,1]
    y2 <- y[,2]
    
    el = exp(f)
    p1 = el/(1 + el)
    
    if(length(or) ==1) or = rep(or,length(y1))
    if(length(mu2) ==1) mu2 = rep(mu2,length(y1))
    
    
      psiminone = or - 1
      hilfs1 = 1 + (p1 + mu2)*psiminone
      hilfs2 = -4*or*psiminone*p1*mu2
      
      if(any(or == 1)){
        p11 <- p1 * mu2
        p11[ which(or != 1)]  = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      } else{
        p11  = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      }
      
      sapply(1:length(y1),function(i)(y1[i]==0&&y2[i]==0) * -log(1 + p11[i] -  mu2[i] - p1[i]) + (y1[i]==1 && y2[i]==0) * -log(p1[i] - p11[i]) + (y1[i]==0 && y2[i]==1) *  -log(mu2[i] - p11[i]) + (y1[i]==1 && y2[i]==1) * -log(p11[i]))
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(or = or, mu2 = mu2, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    y1 <- y[,1]
    y2 <- y[,2]
    
    el = exp(f)
    p1 = el/(1 + el)
    
    if(length(or) ==1) or = rep(or,length(y1))
    if(length(mu2) ==1) mu2 = rep(mu2,length(y1))
    
      psiminone = or - 1
      hilfs1 = 1 + (p1 + mu2)*psiminone
      hilfs2 = -4*or*psiminone*p1*mu2
      
      if(any(or == 1)){
        p11 <- p1 * mu2
        dp11  <-  mu2
        p11[ which(or != 1)]  <-  0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
        dp11[ which(or != 1)] <- 1/2 - (2*hilfs1*psiminone - 4*mu2*or*psiminone)/(4*sqrt(hilfs1^2 + hilfs2)*psiminone)
      } else{
        p11 <-  0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
        dp11 <- 1/2 - (2*hilfs1*psiminone - 4*mu2*or*psiminone)/(4*sqrt(hilfs1^2 + hilfs2)*psiminone)
      }
      
      sapply(1:length(y1), function(i)(y1[i]==0 && y2[i]==0) * (dp11[i] - 1)/(1 + p11[i] -  mu2[i] - p1[i]) + (y1[i]==1 && y2[i]==0) * (1 - dp11[i])/(p1[i] - p11[i]) + (y1[i]==0 && y2[i]==1) *  ( - dp11[i])/(mu2[i] - p11[i]) + (y1[i]==1 && y2[i]==1) * dp11[i]/p11[i] )
             
  }
  
  
  offset <- function(y, w) {
    if (!is.null(mu1)) {
      temp1 <- mu1
      temp1 <- log(temp1/(1-temp1))
      
      RET = temp1
    }else {

      a <- sum(rep(1, length(which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]==0)])) * w[which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]==0)]])
      b <- sum(rep(1, length(which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]!=0)])) * w[which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]!=0)]])
      c <- sum(rep(1, length(which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]==0)])) * w[which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]==0)]])
      d <- sum(rep(1, length(which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]!=0)])) * w[which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]!=0)]])
      
      temp1 <- (c+d)/sum(a+b+c+d)
      RET <- log(temp1/(1-temp1))
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f)/(1 + exp(f)), offset = offset,
                  name = "bivariate Bernoulli distribution: mu1 (logit link)")
}









biP2 <- function(mu1 = NULL, mu2 = NULL, or = NULL){

  
  # neg. log-likelihood
  loss <- function(mu1, or, y, f=f){
    y1 <- y[,1]
    y2 <- y[,2]
    
    el = exp(f)
    p2 = el/(1 + el)
    
    if(length(or) ==1) or = rep(or,length(y1))
    if(length(mu1) ==1) mu1 = rep(mu1,length(y1))
    
    psiminone = or - 1
    hilfs1 = 1 + (mu1 + p2)*psiminone
    hilfs2 = -4*or*psiminone*mu1*p2
    
    if(any(or == 1)){
      p11 <- mu1 * p2
      p11[ which(or != 1)]  = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
    } else{
      p11 = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
    }
    
    sapply(1:length(y1), function(i)(y1[i]==0&&y2[i]==0) * -log(1 + p11[i] -  p2[i] - mu1[i]) + (y1[i]==1 && y2[i]==0) * -log(mu1[i] - p11[i]) + (y1[i]==0 && y2[i]==1) *  -log(p2[i] - p11[i]) + (y1[i]==1 && y2[i]==1) * -log(p11[i]))
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(or = or, mu1 = mu1, y = y, f = f))
  }
  
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    y1 <- y[,1]
    y2 <- y[,2]
    
    el = exp(f)
    p2 = el/(1 + el)
    
    if(length(or) ==1) or = rep(or,length(y1))
    if(length(mu1) ==1) mu1 = rep(mu1,length(y1))
    
    psiminone = or - 1
    hilfs1 = 1 + (mu1 + p2)*psiminone
    hilfs2 = -4*or*psiminone*mu1*p2
    
    if(any(or == 1)){
      p11 <- mu1 * p2
      dp11 <-  mu1
      p11[ which(or != 1)]  <-  0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      dp11[ which(or != 1)] <- 1/2 - (2*hilfs1*psiminone - 4*mu1*or*psiminone)/(4*sqrt(hilfs1^2 + hilfs2)*psiminone)
    } else{
      p11 <-  0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      dp11 <- 1/2 - (2*hilfs1*psiminone - 4*mu1*or*psiminone)/(4*sqrt(hilfs1^2 + hilfs2)*psiminone)
    }
    
    sapply(1:length(y1),function(i) (y1[i]==0 && y2[i]==0) * (dp11[i] - 1)/(1 + p11[i] -  p2[i] - mu1[i]) + (y1[i]==1 && y2[i]==0) * (-dp11[i])/(mu1[i] - p11[i]) + (y1[i]==0 && y2[i]==1) *  (1 - dp11[i])/(p2[i] - p11[i]) + (y1[i]==1 && y2[i]==1) * dp11[i]/p11[i] )
           
  }
  
  offset <- function(y, w) {
    if (!is.null(mu2)) {
      temp2 <- mu2
      temp2 <- log(temp2/(1-temp2))
      RET <- temp2
    }  else {

      a <- sum(rep(1, length(which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]==0)])) * w[which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]==0)]])
      b <- sum(rep(1, length(which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]!=0)])) * w[which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]!=0)]])
      c <- sum(rep(1, length(which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]==0)])) * w[which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]==0)]])
      d <- sum(rep(1, length(which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]!=0)])) * w[which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]!=0)]])
      
      temp2 <- (b+d)/sum(a+b+c+d)
      RET <- log(temp2/(1-temp2))
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f)/(1+exp(f)), offset = offset,
                  name = "bivariate Bernoulli distribution: mu2 (logit link)")
}





biOR <- function(mu1 = NULL, mu2 = NULL, or = NULL){
  # f = exp(OR)

  # neg. log-likelihood
  loss <- function(mu1, mu2, y, f=f){
    y1 <- y[,1]
    y2 <- y[,2]
    
    OR = exp(f)
  
    if(length(mu2) ==1) mu2 = rep(mu2,length(y1))
    if(length(mu1) ==1) mu1 = rep(mu1,length(y1))
    
    
      psiminone = OR - 1
      hilfs1 = 1 + (mu1 + mu2)*psiminone
      hilfs2 = -4*OR*psiminone*mu1*mu2
      
      if(any(OR == 1)){  
        p11 <- mu1 * mu2
        p11[ which(OR != 1)] = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      } else{
        p11 = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      }
      

    sapply(1:length(y1),function(i)(y1[i]==0&&y2[i]==0) * -log(1 + p11[i] -  mu2[i] - mu1[i]) + (y1[i]==1 && y2[i]==0) * -log(mu1[i] - p11[i]) + (y1[i]==0 && y2[i]==1) *  -log(mu2[i] - p11[i]) + (y1[i]==1 && y2[i]==1) * -log(p11[i]) )
           
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu1 = mu1, mu2 = mu2, y = y, f = f))
  }

  
  # ngradient is the negative derivate 
  ngradient <- function(y, f, w = 1){
    y1 <- y[,1]
    y2 <- y[,2]
    
    OR = exp(f)
      
    
    if(length(mu2) ==1) mu2 = rep(mu2,length(y1))
    if(length(mu1) ==1) mu1 = rep(mu1,length(y1))
    
    
      psiminone = OR - 1
      hilfs1 = 1 + (mu1 + mu2)*psiminone
      hilfs2 = -4*OR*psiminone*mu1*mu2
      
      if(any(OR == 1)){ 
        p11 <- mu1 * mu2
        dp11 = 0
        
        p11[which(OR != 1)]  = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
        dp11[which(OR != 1)] <- (
          sqrt(hilfs1^2 + hilfs2) + 
            psiminone * (mu1 + mu2 -(2*hilfs1*(mu1 + mu2) - 4 * mu1 * mu2 * (2*OR -1))/(2*sqrt(hilfs1^2+hilfs2)) )
          - (mu1 + mu2) * (psiminone)- 1 
        )/ (2*psiminone^2)
      } else{
        p11  = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
        dp11 <- (
          sqrt(hilfs1^2 + hilfs2) + 
            psiminone * (mu1 + mu2 -(2*hilfs1*(mu1 + mu2) - 4 * mu1 * mu2 * (2*OR -1))/(2*sqrt(hilfs1^2+hilfs2)) )
          - (mu1 + mu2) * (psiminone)- 1 
        )/ (2*psiminone^2)
      }
      
      sapply(1:length(y1),function(i)(y1[i]==0&&y2[i]==0) * (dp11[i])/(1 + p11[i] -  mu2[i] - mu1[i]) + (y1[i]==1 && y2[i]==0) * (- dp11[i])/(mu1[i] - p11[i]) + (y1[i]==0 && y2[i]==1) * (- dp11[i])/(mu2[i] - p11[i]) + (y1[i]==1 && y2[i]==1) * dp11[i]/p11[i])
             
  }
  
  offset <- function(y, w) {
    if (!is.null(or)) {
      RET <- log(or)
      
    }  else {

            a <- sum(rep(1, length(which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]==0)])) * w[which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]==0)]])
      b <- sum(rep(1, length(which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]!=0)])) * w[which(y[,1] == 0)[which(y[,1]==0) %in% which(y[,2]!=0)]])
      c <- sum(rep(1, length(which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]==0)])) * w[which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]==0)]])
      d <- sum(rep(1, length(which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]!=0)])) * w[which(y[,1] != 0)[which(y[,1]!=0) %in% which(y[,2]!=0)]])
      
      RET <- log((a*d)/(b*c))
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), offset = offset,
                  name = "bivariate Bernoulli distribution: OR (log link)")
}



#------------------ complete gamboostLSS families
BernoulliBV <- function (mu1 = NULL,  mu2 = NULL, or = NULL) 
{

  Families(   mu1 = biP1(mu1 = mu1, mu2 = mu2, or = or), 
              mu2 = biP2(mu1 = mu1, mu2 = mu2, or = or), 
              or  = biOR(mu1 = mu1, mu2 = mu2, or = or), 
              name = "BernoulliBV")
}

