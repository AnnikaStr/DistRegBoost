
#########################################
library("BatchJobs")
library("gamboostLSS")
library('mvtnorm')
library('plyr')
library("scoringRules")
library('R2BayesX')
library('BayesX')
#########################################



source('bivariateGaussian.R')

sim <- function(seed,n.train,p,n.test){
  
  ###### Maps
  # western germany: 327 regions
  map <- read.bnd("westerngermany.bnd")
  pmat <- read.gra("westerngermany.gra")
  
  centr <- as.data.frame((get.centroids(map)))
  helpdata <- data.frame(nr = 1:length(map))
  helpdata$region <- names(map)

  helpdata$xcoord <- centr$centroid_x
  helpdata$ycoord <- centr$centroid_y

  helpdata$xsd <- scale(centr$centroid_x)
  helpdata$ysd <- scale(centr$centroid_y)

  helpdata$falpha <-
    sin(helpdata$xsd) * cos(0.5 * helpdata$ysd)
  helpdata$falpha <-
    scale(helpdata$falpha, center = TRUE, scale = FALSE)
  helpdata <- helpdata[, -1]

  # get regions that are splitted
  regionid <-
    (row.names(centr))[duplicated(row.names(centr))]
  indremove <- c()
  for (i in 1:length(regionid)) {
    indexmap <- which(row.names(centr) == regionid[i])
    centr$centroid_x[indexmap] <-
      mean(centr$centroid_x[indexmap])
    centr$centroid_y[indexmap] <-
      mean(centr$centroid_y[indexmap])
    centr$area[indexmap] <- sum(centr$area[indexmap])
  }
  
  # regions
  findid <- paste0("X",unlist(dimnames(pmat)[2]))
  centroids <- data.frame(region = findid)
  
  for (i in 1:length(findid)) {
    index <- which(findid[i] == row.names(centr))
    centroids$area[i] <- centr$area[index[1]]
    centroids$centroid_x[i] <-
      centr$centroid_x[index[1]]
    centroids$centroid_y[i] <-
      centr$centroid_y[index[1]]
  }
  
  helpdata <- data.frame(nr = 1:length(findid))
  helpdata$region <- findid
  helpdata$xcoord <- centroids$centroid_x
  helpdata$ycoord <- centroids$centroid_y
  
  helpdata$xsd <-
    (centroids$centroid_x - mean(centroids$centroid_x)) / sd(centroids$centroid_x)
  helpdata$ysd <-
    (centroids$centroid_y - mean(centroids$centroid_y)) / sd(centroids$centroid_y)
  
  helpdata$falpha <-
    sin(helpdata$xsd) * cos(0.5 * helpdata$ysd) 
  helpdata$falpha <-
    helpdata$falpha - mean(helpdata$falpha)
  helpdata <- helpdata[, -1]
  
  region <- (unlist(dimnames(pmat)[2]))
  regdat <- (data.frame(region, helpdata$falpha))   
  
  # drawmap(map=map,regionvar="region",plotvar="helpdata.falpha",data=regdat)

  set.seed(seed)
  
  n.mstop = 1500
  
  n = n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  nregion <- floor(n/length(regdat[,1]))
  
  region <- c()
  if(nregion > 0) {
    for(i in 1:nregion) {
      region <- c(region,sample(as.character(regdat[,1]),size=length(regdat[,1]),replace=F))
    }
    region <- c(region,sample(as.character(regdat[,1]),size=(n-nregion*length(regdat[,1])),replace=T))
    
  } else {
    
    region <- c(region,sample(as.character(regdat[,1]),size=(n-nregion*length(regdat[,1])),replace=T)) # sample from regdata (here 1100)
  }
  region # subsample of regions
  
  
  # values for the specific region
  fspatalpha <- rep(0,n)
  for(i in 1:n){
    fspatalpha[i] <- helpdata$falpha[paste0("X",region[i])==helpdata$region] # values for the sample
  }

  
  ###### data generation
  # - training data
  x.train <- matrix(runif(p * n, 0,1), n)
  x.train <- cbind(x.train, 1)
  
  train.eta.mu1 <-  sin(2*x.train[,1])/0.5 + x.train[,6]
  train.eta.mu1 <-  train.eta.mu1 - mean(train.eta.mu1)
  train.eta.mu1 <- train.eta.mu1 + fspatalpha
  
  train.eta.mu2 <-  2 + 3*cos(2*x.train[,2]) + 0.5*x.train[,7]
  train.eta.mu2 <-  train.eta.mu2 - mean(train.eta.mu2)
  train.eta.mu2 <- train.eta.mu2 + fspatalpha
  
  train.eta.sigma1 <-  sqrt(x.train[,3]) * x.train[,3] - .5 *x.train[,8] 
  train.eta.sigma1 <-  train.eta.sigma1 - mean(train.eta.sigma1)
  train.eta.sigma1 <- train.eta.sigma1 + fspatalpha
  train.eta.sigma1 <- exp(train.eta.sigma1)
  
  
  train.eta.sigma2 <- (x.train[,4])*cos(x.train[,4]) + 0.25 *x.train[,9] 
  train.eta.sigma2 <- train.eta.sigma2 - mean(train.eta.sigma2)
  train.eta.sigma2 <- train.eta.sigma2 + fspatalpha
  train.eta.sigma2 <- exp(train.eta.sigma2)
  
  eta.rho <-   log(x.train[,5]^2) + x.train[,10]
  eta.rho <-  eta.rho - mean(eta.rho)
  eta.rho <- eta.rho + fspatalpha
  train.eta.rho <- eta.rho/sqrt(1+eta.rho^2)
  
  
  y.train <- NULL
  for(i in 1:n) {
    y.train <- rbind(y.train, rmvnorm(1,mean = c(train.eta.mu1[i], train.eta.mu2[i])  ,
                                      sigma = matrix(c(train.eta.sigma1[i]^2, train.eta.rho[i]*train.eta.sigma1[i]*train.eta.sigma2[i], train.eta.rho[i]*train.eta.sigma1[i]*train.eta.sigma2[i], train.eta.sigma2[i]^2), 2, 2)))
  }
  
  colnames(y.train) <- c('y1', 'y2')
  dat.train <- data.frame(y.train,x.train, region = as.factor(region))
  
  # - test data
  nregion.test <- floor(n.test/length(regdat[,1]))
  
  region.test <- c()

  if(nregion.test > 0) {
    for(i in 1:nregion.test) {
      region.test <- c(region.test,sample(as.character(regdat[,1]),size=length(regdat[,1]),replace=F))
    }
    region.test <- c(region.test,sample(as.character(regdat[,1]),size=(n.test-nregion.test*length(regdat[,1])),replace=T))
    
  } else {
    
    region.test <- c(region.test,sample(as.character(regdat[,1]),size=(n.test-nregion.test*length(regdat[,1])),replace=T))
  }
  region.test 
  
  fspatalpha.test <- rep(0,n.test)
  for(i in 1:n.test){
    fspatalpha.test[i] <- helpdata$falpha[paste0("X",region[i])==helpdata$region]
  }
  
  x.test <- matrix(runif(p * n.test, 0,1), n.test)
  x.test <- cbind(x.test, 1)
  
  test.eta.mu1 <-  sin(2 * x.test[,1])/0.5  + x.test[,6]
  test.eta.mu1 <-  test.eta.mu1 - mean(test.eta.mu1)
  test.eta.mu1 <- test.eta.mu1 + fspatalpha.test
  
  test.eta.mu2 <-  2 + 3*cos(2*x.test[,2]) + 0.5*x.test[,7] 
  test.eta.mu2 <-  test.eta.mu2 - mean(test.eta.mu2)
  test.eta.mu2 <- test.eta.mu2 + fspatalpha.test
  
  
  test.eta.sigma1 <-  sqrt(x.test[,3])*x.test[,3] - .5 *x.test[,8] 
  test.eta.sigma1 <-  test.eta.sigma1 - mean(test.eta.sigma1)
  test.eta.sigma1 <- test.eta.sigma1 + fspatalpha.test
  test.eta.sigma1 <- exp(test.eta.sigma1)
  
  test.eta.sigma2 <- cos(x.test[,4])*x.test[,4] + 0.25 *x.test[,9]
  test.eta.sigma2 <- test.eta.sigma2 - mean(test.eta.sigma2)
  test.eta.sigma2 <- test.eta.sigma2 + fspatalpha.test
  test.eta.sigma2 <- exp(test.eta.sigma2)
  
  
  eta.rho <-  log(x.test[,5]^2) + x.test[,10]
  eta.rho <-  eta.rho - mean(eta.rho)
  eta.rho <- eta.rho + fspatalpha.test
  test.eta.rho <- eta.rho/sqrt(1+eta.rho^2)
  
  y.test <- NULL
  for(i in 1:n.test) {
    y.test <- rbind(y.test, rmvnorm(1, mean = c(test.eta.mu1[i], test.eta.mu2[i])  ,
                                    sigma = matrix(c(test.eta.sigma1[i]^2, test.eta.rho[i]*test.eta.sigma1[i]*test.eta.sigma2[i], test.eta.rho[i]*test.eta.sigma1[i]*test.eta.sigma2[i], test.eta.sigma2[i]^2), 2, 2)))
    
  }
  
  colnames(y.test) <- c('y1', 'y2')
  dat.test <- data.frame(y.test,x.test, region = as.factor(region.test))
  #######################################
  
  # form high
  # x <- paste(c(paste("bbs(X", 1:5, ")", sep = ""),paste("bols(X", 6:10, ",intercept = F)", sep = ""),paste("bols(X", ((10+1):p) , ",intercept = F)",sep = ""), 'bmrf(region, bnd = map)'), collapse = "+")
  # form <- as.formula(( paste( "cbind(y1,y2) ~",  x)))
  
  # form low
  x <- paste(c(paste("bbs(X", 1:5, ")", sep = ""),paste("bols(X", 6:p, ",intercept = F)", sep = ""), paste("bols(X", (p+1) , ",intercept = F)",sep = ""), 'bmrf(region, bnd = map)'), collapse = "+")
  form <- as.formula(( paste( "cbind(y1,y2) ~",  x)))
  
  bivNorm <- gamboostLSS(form, data = dat.train, families = GaussianBV(mu1 = NULL, mu2 = NULL, sigma1 = NULL,  sigma2 = NULL, rho = NULL), 
                         control = boost_control(mstop = 1000, risk = 'oobag',nu  = 0.1), method = 'noncyclic', weights = weight.mstop)
  
  MSTOP <- which.min(risk(bivNorm,merge = T))
  
  oobag.risk <- risk(bivNorm,merge = T)
  
  rm(bivNorm) 
  dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  bivNorm <- gamboostLSS(form, data = dat.train_biv, families = GaussianBV(mu1 = NULL, mu2 = NULL, sigma1 = NULL,  sigma2 = NULL, rho = NULL), 
                         control = boost_control(mstop = MSTOP, nu  = 0.1), method = 'noncyclic')
  
  mstop.bivNorm <-  vector('list')
  mstop.bivNorm$mstop <- MSTOP
  mstop.bivNorm$mu1 <- bivNorm$mu1$mstop()
  mstop.bivNorm$mu2 <- bivNorm$mu2$mstop()
  mstop.bivNorm$sigma1 <- bivNorm$sigma1$mstop()
  mstop.bivNorm$sigma2 <- bivNorm$sigma1$mstop()
  mstop.bivNorm$rho <- bivNorm$rho$mstop()
  
  coef.bivNorm <- coef(bivNorm,which = '')
  coef.bivNorm_ow <- coef(bivNorm)
  
  
  ##### for plotting regions
  fitted_mu1 <- fitted(bivNorm, parameter = 'mu1', which = 'region', type = 'response')
  if(fitted_mu1 == 0){
    fitted_mu1 = 0
  }else{
    fitted_mu1 <- tapply(fitted_mu1, dat.train_biv$region, FUN = mean)
  }
  
  fitted_mu2 <- fitted(bivNorm, parameter = 'mu2', which = 'region', type = 'response')
  if(fitted_mu2 == 0){
    fitted_mu2 = 0
  }else{
    fitted_mu2 <- tapply(fitted_mu2, dat.train_biv$region, FUN = mean)
  }
  fitted_sigma1 <- fitted(bivNorm, parameter = 'sigma1', which = 'region', type = 'response')
  if(fitted_sigma1 == 0){
    fitted_sigma1 = 0 
  }else{
    fitted_sigma1 <- tapply(fitted_sigma1, dat.train_biv$region, FUN = mean)
  }
  
  fitted_sigma2 <- fitted(bivNorm, parameter = 'sigma2', which = 'region', type = 'response')
  if(fitted_sigma2 == 0){
    fitted_sigma2 = 0
  }else{
    fitted_sigma2 <- tapply(fitted_sigma2, dat.train_biv$region, FUN = mean)
  }
  
  fitted_rho <- fitted(bivNorm, parameter = 'rho', which = 'region', type = 'response')
  if(fitted_rho == 0){
    fitted_rho = 0
  }else{
    fitted_rho <- tapply(fitted_rho, dat.train_biv$region, FUN = mean)
  }
  
  plot.data <- list(region = names(fitted_mu1), 
                    mu1 = fitted_mu1, mu2 = fitted_mu2,
                    sigma1 = fitted_sigma1, sigma2 = fitted_sigma2,
                    rho  = fitted_rho) 
  
  MSEPbivNorm <-  vector('list')
  lik <- vector('list')
  plot<- list('list')
  
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = p), nrow = 500 )
  plot$newX_mu1 <- predict(bivNorm, data.frame(newX), parameter = 'mu1',which = 1, type = 'link' )
  plot$newX_mu2 <- predict(bivNorm, data.frame(newX), parameter = 'mu2',which = 2, type = 'link' )
  plot$newX_sigma1 <- predict(bivNorm, data.frame(newX), parameter = 'sigma1',which = 3, type = 'link' )
  plot$newX_sigma2 <- predict(bivNorm, data.frame(newX), parameter = 'sigma2',which = 4, type = 'link' )
  plot$newX_rho <- predict(bivNorm, data.frame(newX), parameter = 'rho',which = 5, type = 'link' )
  
  
  pred.mu1 <- as.numeric(predict(bivNorm$mu1, newdata = dat.test, type = 'response'))
  pred.mu2 <- predict(bivNorm$mu2, newdata = dat.test, type = 'response')
  pred.sigma1 <- predict(bivNorm$sigma1, newdata = dat.test, type = 'response')
  pred.sigma2 <- predict(bivNorm$sigma2, newdata = dat.test, type = 'response')
  pred.rho <- predict(bivNorm$rho, newdata = dat.test, type = 'response')
  
  loss <- function( mu1, sigma1, mu2, sigma2, rho,  y){ 
    y1 <- y[,1]
    y2 <- y[,2]
    
    -(-log(2*pi) - log(sigma1) - log(sigma2) - .5*log(1 - rho^2) - 1/(2*(1 - rho^2)) *((y1 - mu1)^2/sigma1^2  +  (y2 - mu2)^2/sigma2^2  
                                                                                       - 2*rho*((y1 - mu1)*(y2 - mu2)/(sigma1*sigma2))))
  } 
  
  MSEPbivNorm$mu1 <- mean((pred.mu1 - dat.test$y1)^2)
  MSEPbivNorm$mu2 <- mean((pred.mu2 - dat.test$y2)^2)
  lik$biv <- sum(loss(mu1 = pred.mu1, mu2 = pred.mu2, sigma1 = pred.sigma1, sigma2 = pred.sigma2, rho = pred.rho, y = y.test))
  
  
  #############################################
  ################ - univariat - ##############
  #############################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  MSEP.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[, !(names(dat.test) %in% "y2")]
  
  form1 <- as.formula(( paste( "y1 ~",  x)))
  dat.train.mu1$y1 <- dat.train.mu1$y1
  
  glm.uni.mu1 <- gamboostLSS(form1 , data = dat.train.mu1, families = GaussianLSS(), control = boost_control(nu = 0.1, risk = "oobag", mstop = 1000), method = 'noncyclic', weights = weight.mstop)
  
  
  mstop.uni$mstop_mu1 <-  which.min(risk(glm.uni.mu1,merge = T))
  dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  
  glm.uni.mu1 <- gamboostLSS(form1 , data = dat.train_mu1, families = GaussianLSS(), control = boost_control(nu = 0.1, mstop = mstop.uni$mstop_mu1), method = 'noncyclic')
  
  
  mstop.uni$mu1 <- glm.uni.mu1$mu$mstop()
  mstop.uni$mu1 <- glm.uni.mu1$sigma$mstop()
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')

  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2= dat.test[, !(names(dat.test) %in% "y1")]
  
  form2 <- as.formula(( paste( "y2 ~",  x)))
  dat.train.mu2$y2 <- dat.train.mu2$y2
  
  
  glm.uni.mu2 <- gamboostLSS(form2 , data = dat.train.mu2, families = GaussianLSS(), control = boost_control(nu = 0.1, risk = "oobag", mstop = 1000), method = 'noncyclic', weights = weight.mstop)
  
  mstop.uni$mstop_mu2 <-  which.min(risk(glm.uni.mu2,merge = T))
  dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  
  glm.uni.mu2 <- gamboostLSS(form2 , data = dat.train_mu2, families = GaussianLSS(), control = boost_control(nu = 0.1, mstop = mstop.uni$mstop_mu1), method = 'noncyclic')

  mstop.uni$mu2 <- glm.uni.mu2$mu$mstop()
  mstop.uni$mu2 <- glm.uni.mu2$sigma$mstop()
  
  coef.uni$mu2 <- coef(glm.uni.mu2, which = '')

  ####### Predictive Performance #######
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test, type = 'response'))
  pred.sigma1.uni <- predict(glm.uni.mu1$sigma, newdata = dat.test, type = 'response')
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test, type = 'response'))
  pred.sigma2.uni <- predict(glm.uni.mu2$sigma, newdata = dat.test, type = 'response')
  
  # MSEP
  MSEP.uni$mu1 <- mean((pred.mu1.uni-(as.numeric(dat.test.mu1$y1)))^2)
  MSEP.uni$mu2 <- mean((pred.mu2.uni-(as.numeric(dat.test.mu2$y2)))^2)
  

  mu1.uni.loglik <- -sum(dnorm(x = dat.test.mu1$y1,mean = pred.mu1.uni, sd = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- -sum(dnorm(x = dat.test.mu2$y2, mean = pred.mu2.uni, sd = pred.sigma2.uni, log = T))
  
  lik$uni <- mu1.uni.loglik + mu2.uni.loglik
  lik$uni_lossNorm <- sum(loss(mu1 = pred.mu1.uni, mu2 =  pred.mu2.uni, sigma1 = pred.sigma1.uni, sigma2  = pred.sigma2.uni, rho = rep(0,length(pred.mu1.uni)), y = y.test))
  
  plot_uni <- vector('list')
  newX <- matrix(rep(seq(from = 0, to = 1,length.out =500),times = p),nrow =500 )
  plot_uni$newX_mu1 <- predict(glm.uni.mu1$mu, data.frame(newX),which = 1)
  plot_uni$newX_sigma1 <- predict(glm.uni.mu1$sigma, data.frame(newX),which = 3)
  plot_uni$newX_mu2 <- predict(glm.uni.mu2, data.frame(newX),which=2)
  plot_uni$newX_sigma2 <- predict(glm.uni.mu2, data.frame(newX),which=4)
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- list(pred.mu1, pred.mu2, pred.sigma1, pred.sigma2, pred.rho)
  pred.uni <- list(pred.mu1.uni, pred.sigma1.uni, pred.mu2.uni, pred.sigma2.uni)
  
  
  fitted_mu1_uni <- fitted(glm.uni.mu1, parameter = 'mu', which = 'region', type = 'response')
  if(fitted_mu1_uni == 0){
    fitted_mu1_uni = 0
  }else{
    fitted_mu1_uni <- tapply(fitted_mu1_uni, dat.train_mu1$region, FUN = mean)
  }
  
  fitted_sigma1_uni <- fitted(glm.uni.mu1, parameter = 'sigma', which = 'region', type = 'response')
  if(fitted_sigma1_uni == 0){
    fitted_sigma1_uni = 0
  }else{
    fitted_sigma1_uni <- tapply(fitted_sigma1_uni, dat.train_mu1$region, FUN = mean)
  }
  
  fitted_mu2_uni <- fitted(glm.uni.mu2, parameter = 'mu', which = 'region', type = 'response')
  if(fitted_mu2_uni == 0){
    fitted_mu2_uni = 0 
  }else{
    fitted_mu2_uni <- tapply(fitted_mu2_uni, dat.train_mu2$region, FUN = mean)
  }
  
  fitted_sigma2_uni <- fitted(glm.uni.mu2, parameter = 'sigma', which = 'region', type = 'response')
  if(fitted_sigma2_uni == 0){
    fitted_sigma2_uni = 0
  }else{
    fitted_sigma2_uni <- tapply(fitted_sigma2_uni, dat.train_mu2$region, FUN = mean)
  }
  
  plot.data_uni <- list(region = names(fitted_mu1_uni),
                        mu1 = fitted_mu1_uni, sigma1 = fitted_sigma1_uni,
                        mu2 = fitted_mu1_uni, sigma2 = fitted_sigma1_uni)
  
  
  ################################
  ######## Energy Score ##########
  ################################
  es_biv <- vector()
  es_uni <- vector()
  
  for(i in 1:length(pred.mu1.uni)){
    
    pred_sample_uni <- matrix(NA, nrow = 2, ncol = 10000)
    pred_sample_biv <- matrix(NA, nrow = 2, ncol = 10000)
    
    # univariate
    sample_uni <- rmvnorm(10000, mean = c(pred.mu1.uni[i], pred.mu2.uni[i]) ,
                          sigma = matrix(c(pred.sigma1.uni[i]^2, 0, 
                                           0, pred.sigma2.uni[i]^2), 2, 2))
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    # bivariate
    sample_biv <- rmvnorm(10000, mean = c(pred.mu1[i], pred.mu2[i]) ,
                          sigma = matrix(c(pred.sigma1[i]^2, pred.rho[i]*pred.sigma1[i]*pred.sigma2[i], 
                                           pred.rho[i]*pred.sigma1[i]*pred.sigma2[i], pred.sigma2[i]^2), 2, 2))
    
    pred_sample_biv[1, ] <- sample_biv[,1]
    pred_sample_biv[2, ] <- sample_biv[,2]
    
    es_biv[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_biv) 
    
  }
  
  energy_score <- list()
  energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  
  return(list(n = n.ges, p  = p, Coefficients = coef.bivNorm, MSEPbivNorm = MSEPbivNorm, Likelihood = lik, oobag.risk = oobag.risk, plot = plot, plot_uni = plot_uni, energy_score = energy_score,
              predict = pred.ges, predict.uni = pred.uni, mstop = mstop.bivNorm, plot.data = plot.data, plot.data_uni = plot.data_uni, Coefficients.uni = coef.uni, MSEP.uni = MSEP.uni, mstop.uni = mstop.uni))
  
}


n.train = 1000
p = 10
n.test = 1000

results = mclapply(1:100, sim, mc.cores = 10, n.train = n.train, p  = p, n.test = n.test, mc.preschedule = FALSE)

results$n.train= n.train
results$p = p
results$n.test = n.test

save(results, file="Sim_bivGaussian.RData")


