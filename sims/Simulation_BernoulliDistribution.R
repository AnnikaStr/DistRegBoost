#########################################
library('VGAM')
library("gamboostLSS")
library('mvtnorm')
library('pROC')
library('scoringRules')
#########################################


source('bivariateBernoulli.R')

sim <- function(seed,n.train,p,n.test,corr){
  
  ## Likelihood 
  loss <- function(mu1, mu2, or, y){
    y1 <- y[,1]
    y2 <- y[,2]
    
    
    loss_vec <- c()
    for(i in 1:length(y1)){
      
      psiminone = or[i] - 1
      
      hilfs1 = 1 + (mu1[i] + mu2[i])*psiminone
      hilfs2 = -4*or[i]*psiminone*mu1[i]*mu2[i]
      
      if(or[i] == 1) {
        p11 <- mu1[i] * mu2[i]
      }else{
        p11 = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      }
      
      if(y1[i] == 0 && y2[i] == 0){
        loss_vec[i] <- -log(1 + p11 - mu2[i] - mu1[i])
      } else if(y1[i] == 0 && y2[i] == 1){
        loss_vec[i] <- -log(mu2[i] - p11)
      } else if(y1[i] == 1 && y2[i] == 0){
        loss_vec[i] <- -log(mu1[i] - p11)
      } else{
        loss_vec[i] <- -log(p11)
      }
      
      
    }
    loss_vec
  }
  
  

  mu1 = c(1, 1.5, -1, 1.5, 0, 0, rep(0,p-6))
  mu2 = c(2, -1, 1, 0, 0,rep(0,p-5)) 
  or  = c(0, 0, 0, 0, 1, 1.5, 0, 0,rep(0,p-8))
  

  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n.mstop = 1500
  n = n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  # - training data
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  train.eta.mu1_center<-  (x.train %*% mu1)
  train.eta.mu1<-  logitlink(train.eta.mu1_center, inv = TRUE)
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  logitlink(train.eta.mu2, inv = T)
  train.eta.or_center <- -1.5+ (x.train %*% or)
  train.eta.or <-  exp( train.eta.or_center)
  
  y.train <- rbinom2.or(n, mu1 = train.eta.mu1,  mu2 = train.eta.mu2, oratio = train.eta.or, exch = F)
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  # - test data
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  test.eta.mu1<-  logitlink(x.test %*% mu1 , inv = TRUE)
  test.eta.mu2 <-  logitlink(x.test %*% mu2, inv = T)
  test.eta.or <- exp(-1.5 + x.test %*% or)
  
  y.test <- rbinom2.or(n.test, mu1 = test.eta.mu1,  mu2 = test.eta.mu2, oratio = test.eta.or, exch = F)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  # - Model
  bivBern <- glmboostLSS(cbind(y1,y2)~., data = dat.train, families = BernoulliBV(mu1 = NULL, mu2 = NULL, or = NULL), 
                         control = boost_control(mstop = 100, risk = 'oobag', nu  = 0.1), method = 'noncyclic', weights = weight.mstop)
  
  MSTOP <- which.min(risk(bivBern,merge = T))

  if(MSTOP >= 90){
    mstop(bivBern) <- 300
  }

  MSTOP <- which.min(risk(bivBern,merge = T))

  if(MSTOP >= 290){
    bivBern[1000]
  }

  MSTOP <- which.min(risk(bivBern,merge = T))

  if(MSTOP >= 990){
    bivBern[2000]
  }
  MSTOP <- which.min(risk(bivBern,merge = T))
  oobag.risk <- risk(bivBern,merge = T)
  
  rm(bivBern)
  dat.train_biv <- dat.train[weight.mstop == 1, ]

  
  bivBern <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = BernoulliBV(mu1 = NULL, mu2 = NULL, or = NULL), 
                         control = boost_control(mstop = MSTOP, nu  = 0.1), method = 'noncyclic')
  
  
  mstop.bivBern <-  vector('list')
  mstop.bivBern$mstop <- MSTOP
  mstop.bivBern$mu1 <- bivBern$mu1$mstop()
  mstop.bivBern$mu2 <- bivBern$mu2$mstop()
  mstop.bivBern$or  <- bivBern$or$mstop()
  
  coef.bivBern <- coef(bivBern, which = "")
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')

  pred.mu1 <- predict(bivBern$mu1, newdata = dat.test, type = 'response')
  pred.mu2 <- predict(bivBern$mu2, newdata = dat.test, type = 'response')
  pred.or <- predict(bivBern$or, newdata = dat.test, type = 'response')
  
  AUCbivBern$mu1 <- roc(dat.test$y1, as.numeric(pred.mu1))$auc
  AUCbivBern$mu2 <- roc(dat.test$y2, as.numeric(pred.mu2))$auc
  
  BrierbivBern$mu1 <- mean((as.numeric(pred.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBern$mu2 <- mean((as.numeric(pred.mu2)-(as.numeric(dat.test$y2)))^2)
  
  lik$biv <- sum(loss(mu1 = pred.mu1, mu2 = pred.mu2, or = pred.or, y = y.test))

  #############################################
  ################ - univariat - ##############
  #############################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , data = dat.train.mu1, family = Binomial(), control = boost_control(risk = 'oobag'), weights = weight.mstop)

  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 90){
    glm.uni.mu1[500]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 490){
    glm.uni.mu1[1000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[1500]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , data = dat.train_mu1, family = Binomial(), 
                              control = boost_control(mstop = mstop.uni$mu1, nu  = 0.1))
    
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , data = dat.train.mu2, family = Binomial(), control = boost_control(risk = 'oobag'), weights = weight.mstop)
 
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 90){
    glm.uni.mu2[500]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 490){
    glm.uni.mu2[1000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 990){
    glm.uni.mu2[1500]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(), 
                          control = boost_control(mstop = mstop.uni$mu1, nu  = 0.1))
  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  lik$uni_lossBern <- sum(loss(mu1 = predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"), mu2 =  predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"), or = rep(1, dim(y.test)[1]), y = y.test))
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  
  ################################
  ######## Energy Score ##########
  ################################
  es_biv <- vector()
  es_uni<- vector()

  for(i in 1:length(pred.mu1.uni)){
    
    pred_sample_uni <- matrix(NA, nrow = 2, ncol = 10000)
    pred_sample_biv <- matrix(NA, nrow = 2, ncol = 10000)
    
    
    # univariate
    sample_uni <- rbinom2.or(10000, mu1 =  pred.mu1.uni[i], mu2 = pred.mu2.uni[i], oratio = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 

    # bivariate
    sample_biv <- rbinom2.or(10000, mu1 = pred.mu1[i],  mu2 = pred.mu2[i], oratio = pred.or[i], exch = F)

    pred_sample_biv[1, ] <- sample_biv[,1]
    pred_sample_biv[2, ] <- sample_biv[,2]
    
    es_biv[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_biv) 
    
  }
  
  energy_score <- list()
  energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)

  

  return(list(TrueBeta = TrueBeta, n = n.ges, p  = p, Coefficients = coef.bivBern, oobag.risk = oobag.risk, AUCbivBern = AUCbivBern, BrierbivBern = BrierbivBern, Likelihood = lik,
               predict = pred.ges, predict.uni = pred.uni, energy_score = energy_score,  mstop = mstop.bivBern, Coefficients.uni = coef.uni, AUC.uni = AUC.uni, Brier.uni = Brier.uni, mstop.uni = mstop.uni))
  
}


n.train = 1000
p = 10
n.test = 1000
corr = 0.5

results = mclapply(1:100, sim, mc.cores = 10, n.train = n.train, p  = p, corr = corr, n.test = n.test, mc.preschedule = FALSE)

ssresults$n.train = n.train
results$p = p
results$corr = corr
results$n.test = n.test

save(results, file="Sim_bivBernoulli.RData")
