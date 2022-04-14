#########################################
require("extraDistr")  
library("BatchJobs")
library("gamboostLSS")
library('mvtnorm')
library('plyr')
library('pROC')
library('scoringRules')
#########################################

source('bivariatePoisson.R')

sim <- function(seed,n.train,p,n.test){
  
  loss <- function(lambda1, lambda2, lambda3, y){
    y2 <- y[,2]
    y1 <- y[,1]
    
    x <- sapply(1:length(y1), function(i) min(y1[i],y2[i]))
    B <- sapply(1:length(x),function(j){sum(sapply(0:x[j], function(a) {choose(y1[j],a) * choose(y2[j],a) * factorial(a) * (lambda3[j] / (lambda1[j]*lambda2[j]))^a}))}) 
    
    sapply(1:length(y1), function(i) lambda3[i] + lambda1[i] + lambda2[i] - y1[i] * log(lambda1[i]) + lfactorial(y1[i]) - y2[i] * log(lambda2[i]) + lfactorial(y2[i])- log(B[i])) 
    
  }
  
  
  TrueBeta <-  vector('list')
  TrueBeta$lambda1 <- paste(" sqrt(x.train[,1])* x.train[,1]")
  TrueBeta$lambda2 <-  paste("cos(X[,2])")
  TrueBeta$lambda3 <- paste("sin(X[,3])")
  
  set.seed(seed)
  
  n.mstop = 1500
  
  n = n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop))
  
  # - training data
  x.train <- matrix(runif(p * n, 0,1), n)
  x.train <- cbind(x.train, 1)

  eta.lambda1 <-  sqrt(x.train[,1])* x.train[,1]
  eta.lambda2 <-  cos(x.train[,2])
  eta.lambda3 <-  sin(x.train[,3])

  train.eta.lambda1<-  exp(eta.lambda1)
  train.eta.lambda2 <- exp(eta.lambda2)
  train.eta.lambda3 <- exp(eta.lambda3)

  y.train <- rbvpois(n, train.eta.lambda1, train.eta.lambda2, train.eta.lambda3)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  

  # - test data
  x.test <- matrix(runif(p * n.test, 0,1), n.test)
  x.test <- cbind(x.test, 1)
  
  eta.lambda1_test <-  sqrt(x.test[,1])*x.test[,1]
  eta.lambda2_test <-  cos(x.test[,2])
  eta.lambda3_test <-  sin(x.test[,3])
  
  test.eta.lambda1<-  exp(eta.lambda1_test)
  test.eta.lambda2 <- exp(eta.lambda2_test)
  test.eta.lambda3 <- exp(eta.lambda3_test)
  
  y.test <- rbvpois(n.test, test.eta.lambda1, test.eta.lambda2, test.eta.lambda2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  x <- paste(c(paste("bbs(X", 1:p, ")", sep = ""),paste("bols(X", (p+1) , ",intercept = F)",sep = "")), collapse = "+")
  form <- as.formula(( paste( "cbind(y1,y2) ~",  x)))
  
  # - Model
  bivPois <- gamboostLSS(form, data = dat.train, families = PoissonBV(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL), 
                         control = boost_control(mstop = 100, risk = 'oobag',nu  = 0.1), method = 'noncyclic', weights = weight.mstop)
  
  
  MSTOP <- which.min(risk(bivPois,merge = T))
  
  if(MSTOP >= 90){
    bivPois[500]
  }
  
  MSTOP <- which.min(risk(bivPois,merge = T))
  
  if(MSTOP >= 490){
    bivPois[100]
  }
  
  MSTOP <- which.min(risk(bivPois,merge = T))
  
  
  if(MSTOP >= 990){
    bivPois[1500]
  }
  
  
  MSTOP <- which.min(risk(bivPois,merge = T))
  oobag.risk <- risk(bivPois,merge = T)
  
  rm(bivPois) 
  dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  bivPois <- gamboostLSS(form, data = dat.train_biv, families = PoissonBV(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL), 
                         control = boost_control(mstop = MSTOP, nu  = 0.1), method = 'noncyclic')
  
  mstop.bivPois <-  vector('list')
  mstop.bivPois$mstop <- MSTOP
  mstop.bivPois$lambda1 <- bivPois$lambda1$mstop()
  mstop.bivPois$lambda2 <- bivPois$lambda2$mstop()
  mstop.bivPois$lambda3 <- bivPois$lambda3$mstop()
  nu  =  boost_control(bivPois)$nu
  
  coef.bivPois <- coef(bivPois, which = '')
  
  MSEPbivPois <-  vector('list')
  lik <- vector('list')
  plot<- vector('list')
  
  newX <- matrix(rep(seq(from = 0, to = 1,length.out =500),times = p),nrow =500 )
  plot$newX_lambda1 <- predict(bivPois, data.frame(newX), parameter = 'lambda1',which = 1, type = 'link' )
  plot$newX_lambda2 <- predict(bivPois, data.frame(newX), parameter = 'lambda2',which = 2, type = 'link' )
  plot$newX_lambda3 <- predict(bivPois, data.frame(newX), parameter = 'lambda3',which = 3, type = 'link' )
  
  pred.lambda1 <- as.numeric(predict(bivPois$lambda1, newdata = dat.test, type = 'response'))
  pred.lambda2 <- predict(bivPois$lambda2, newdata = dat.test, type = 'response')
  pred.lambda3 <- predict(bivPois$lambda3, newdata = dat.test, type = 'response')
  
  
  MSEPbivPois$lambda1 <- mean((pred.lambda1 - dat.test$y1)^2)
  MSEPbivPois$lambda2 <- mean((pred.lambda2 - dat.test$y2)^2)
  lik$biv <- sum(loss(lambda1 = pred.lambda1, lambda2 = pred.lambda2, lambda3 = pred.lambda3, y = y.test))
  
  #############################################
  ################ - univariat - ##############
  #############################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  MSEP.uni <- vector('list')
  
  ##########
  # --- mu1
  dat.train.lambda1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.lambda1  = dat.test[, !(names(dat.test) %in% "y2")]
  
  x <- paste(c(paste("bbs(X", 1:p, ")", sep = ""),paste("bols(X", (p+1) , ",intercept = F)",sep = "")), collapse = "+")
  form1 <- as.formula(( paste( "y1 ~",  x)))
  
 
  glm.uni.lambda1 <- gamboost(form1, data = dat.train.lambda1, family = Poisson(), control = boost_control(nu = 0.1, risk = 'oobag'), weights = weight.mstop)
 
  mstop.uni$lambda1 <- which.min(risk(glm.uni.lambda1))
  if(mstop.uni$lambda1 >= 90){
    glm.uni.lambda1[500]
  }
  mstop.uni$lambda1 <- which.min(risk(glm.uni.lambda1))
  
  if(mstop.uni$lambda1 >= 490){
    glm.uni.lambda1[1000]
  }
  mstop.uni$lambda1 <- which.min(risk(glm.uni.lambda1))
  
  dat.train_lambda1 <- dat.train.lambda1[weight.mstop == 1, ]
  
  glm.uni.lambda1 <- gamboost(form1, data = dat.train_lambda1, family = Poisson(), 
                              control = boost_control(mstop = mstop.uni$lambda1, nu  = 0.1))
  
  coef.uni$lambda1 <- coef(glm.uni.lambda1, which = '')
  pred.lambda1.uni <- as.numeric(predict(glm.uni.lambda1, newdata = dat.test.lambda1,type = "response"))
  MSEP.uni$lambda1 <- mean((pred.lambda1.uni-(as.numeric(dat.test.lambda1$y1)))^2)
  
  ##########
  # --- mu2
  dat.train.lambda2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.lambda2= dat.test[, !(names(dat.test) %in% "y1")]
  
  form2 <- as.formula(( paste( "y2 ~",  x)))
  
  glm.uni.lambda2 <- gamboost(form2, data = dat.train.lambda2, family = Poisson(), control = boost_control(nu = 0.1, risk = 'oobag'), weights = weight.mstop)
  
  mstop.uni$lambda2 <- which.min(risk(glm.uni.lambda2))
  if(mstop.uni$lambda2 >= 90){
    glm.uni.lambda2[500]
  }
  mstop.uni$lambda2 <- which.min(risk(glm.uni.lambda2))
  
  if(mstop.uni$lambda2 >= 490){
    glm.uni.lambda2[1000]
  }
  mstop.uni$lambda2 <- which.min(risk(glm.uni.lambda2))
  
  dat.train_lambda2 <- dat.train.lambda2[weight.mstop == 1, ]
  
  glm.uni.lambda2 <- gamboost(form2, data = dat.train_lambda2, family = Poisson(), 
                              control = boost_control(mstop = mstop.uni$lambda2, nu  = 0.1))
                              
                              
  coef.uni$lambda2 <- coef(glm.uni.lambda2, which = '')
  pred.lambda2.uni <- as.numeric(predict(glm.uni.lambda2, newdata = dat.test.lambda2,type = "response"))
  MSEP.uni$lambda2 <- mean((pred.lambda2.uni-(as.numeric(dat.test.lambda2$y2)))^2)

  
  
  # Likelihood
  lambda1.uni.loglik <- sum(-dpois(x = dat.test.lambda1$y1,lambda = pred.lambda1.uni, log = T))
  lambda2.uni.loglik <- sum(-dpois(x = dat.test.lambda2$y2, lambda = pred.lambda2.uni, log = T))
  
  lik$uni <- lambda1.uni.loglik + lambda2.uni.loglik
  
  lik$uni_lossPois <- sum(loss(lambda1 = predict(glm.uni.lambda1, newdata = dat.test.lambda1,type = "response"), lambda2 =  predict(glm.uni.lambda2, newdata = dat.test.lambda2,type = "response"), lambda3 = rep(0, dim(y.test)[1]), y = y.test))

  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- list(pred.lambda1, pred.lambda2, pred.lambda3)
  pred.uni <- list(pred.lambda1.uni, pred.lambda2.uni)
  
  plot_uni<- vector('list')
  
  newX <- matrix(rep(seq(from = 0, to = 1,length.out =500),times = p),nrow =500 )
  plot_uni$newX_lambda1 <- predict(glm.uni.lambda1, data.frame(newX),which=1)
  plot_uni$newX_lambda2 <- predict(glm.uni.lambda2, data.frame(newX),which=2)
  
  ################################
  ######## Energy Score ##########
  ################################
  
  
  es_biv <- vector()
  es_uni <- vector()
  
  for(i in 1:length(pred.lambda1.uni)){
    
    pred_sample_uni <- matrix(NA, nrow = 2, ncol = 10000)
    pred_sample_biv <- matrix(NA, nrow = 2, ncol = 10000)
    
    # univariate
    sample_uni <- rbvpois(10000, pred.lambda1.uni[i], pred.lambda2.uni[i], 0)
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    # bivariate
    sample_biv <- rbvpois(10000, pred.lambda1[i], pred.lambda2[i], pred.lambda3[i])
    pred_sample_biv[1, ] <- sample_biv[,1]
    pred_sample_biv[2, ] <- sample_biv[,2]
    
    es_biv[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_biv) 
    
  }
  
  energy_score <- list()
  energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)

  
  return(list(TrueBeta = TrueBeta, n = n.ges, p  = p, Coefficients = coef.bivPois, MSEPbivPois = MSEPbivPois, Likelihood = lik, oobag.risk = oobag.risk, nu = nu,
               plot = plot, plot_uni = plot_uni, predict = pred.ges, predict.uni = pred.uni, energy_score = energy_score, mstop = mstop.bivPois, Coefficients.uni = coef.uni, MSEP.uni = MSEP.uni, mstop.uni = mstop.uni))
  
}

n.train = 1000
p = 10
n.test = 1000

results = mclapply(1:100, sim, mc.cores = 10, n.train = n.train, p  = p, n.test = n.test, mc.preschedule = FALSE)


results$n.train= n.train
results$p = p
results$n.test = n.test

save(results, file="Sim_bivPoisson_nonlin.RData")
