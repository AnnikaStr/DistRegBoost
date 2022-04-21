#########################################
library('gamboostLSS')
library('BayesX')
source('bivariateGaussian.R')
#########################################

sim <- function(seed){

  load('Nigeria3.RData')
  neighborhood = read.bnd('nigeriamap.bnd')

  Nigeria3 <- Nigeria3[-which(Nigeria3$region == 0),]
  Nigeria3$cbirthorder <- as.factor(Nigeria3$cbirthorder)
  Nigeria3$region <- as.factor(Nigeria3$region)
  Nigeria3$subregion <- as.factor(Nigeria3$subregion)
  
  set.seed(seed)
  
  form = as.formula(cbind(stunting2, wasting2) ~ bbs(cage) + bbs(edupartner) + bbs(mage) + bbs(mbmi) 
                     + bols(bicycle) + bols(car) + bols(cbirthorder) + bols(csex) + bols(ctwin)
                     + bols(electricity) + bols(motorcycle) + bols(mresidence) + bols(munemployed)
                     + bols(radio) + bols(refrigerator) + bols(television) + bmrf(subregion, bnd = neighborhood))
  
  gam_nigeria <- gamboostLSS(formula =  form, data = Nigeria3, families = GaussianBV(),
                             control = boost_control(mstop = 1), method = 'noncyclic')
  
  cvr_nigeria3 <- cvrisk(gam_nigeria, folds = cv(model.weights(gam_nigeria)), grid = 1:10000, papply = mclapply)
  MSTOP <- mstop(cvr_nigeria3)
  
  mstop.bivNorm <-  vector('list')
  mstop.bivNorm$mstop <- MSTOP
  mstop.bivNorm$mu1 <- gam_nigeria$mu1$mstop()
  mstop.bivNorm$mu2 <- gam_nigeria$mu2$mstop()
  mstop.bivNorm$sigma1 <- gam_nigeria$sigma1$mstop()
  mstop.bivNorm$sigma2 <- gam_nigeria$sigma1$mstop()
  mstop.bivNorm$rho <- gam_nigeria$rho$mstop()
  
  gam_nigeria[MSTOP]
  
  coef.nigeria <- coef(gam_nigeria, which = '')
  
  ##### for plotting regions
  fitted_mu1 <- fitted(gam_nigeria, parameter = 'mu1', which = 'subregion', type = 'response')
  fitted_mu1 <- tapply(fitted_mu1, Nigeria3$subregion, FUN = mean)
  fitted_mu2 <- fitted(gam_nigeria, parameter = 'mu2', which = 'subregion', type = 'response')
  fitted_mu2 <- tapply(fitted_mu2, Nigeria3$subregion, FUN = mean)
  fitted_sigma1 <- fitted(gam_nigeria, parameter = 'sigma1', which = 'subregion', type = 'response')
  fitted_sigma1 <- tapply(fitted_sigma1, Nigeria3$subregion, FUN = mean)
  fitted_sigma2 <- fitted(gam_nigeria, parameter = 'sigma2', which = 'subregion', type = 'response')
  fitted_sigma2 <- tapply(fitted_sigma2, Nigeria3$subregion, FUN = mean)
  fitted_rho <- fitted(gam_nigeria, parameter = 'rho', which = 'subregion', type = 'response')
  fitted_rho <- tapply(fitted_rho, Nigeria3$subregion, FUN = mean)
  
  plot.data <- data.frame(region = names(fitted_mu1), 
                          mu1 = fitted_mu1, mu2 = fitted_mu2,
                          sigma1 = fitted_sigma1, sigma2 = fitted_sigma2,
                          rho  = fitted_rho) 

  predict_mu1 <- vector('list')
  predict_mu1$cage <- predict(gam_nigeria, Nigeria3, parameter = 'mu1', which = 'cage', type = 'link')
  predict_mu1$edupartner <- predict(gam_nigeria, Nigeria3, parameter = 'mu1', which = 'edupartner', type = 'link')
  predict_mu1$mage <- predict(gam_nigeria, Nigeria3, parameter = 'mu1', which = 'mage', type = 'link')
  predict_mu1$mbmi <- predict(gam_nigeria, Nigeria3, parameter = 'mu1', which = 'mbmi', type = 'link')
  
  predict_mu2 <- vector('list')
  predict_mu2$cage <- predict(gam_nigeria, Nigeria3, parameter = 'mu2', which = 'cage', type = 'link')
  predict_mu2$edupartner <- predict(gam_nigeria, Nigeria3, parameter = 'mu2', which = 'edupartner', type = 'link')
  predict_mu2$mage <- predict(gam_nigeria, Nigeria3, parameter = 'mu2', which = 'mage', type = 'link')
  predict_mu2$mbmi <- predict(gam_nigeria, Nigeria3, parameter = 'mu2', which = 'mbmi', type = 'link')
  
  predict_sigma1 <- vector('list')
  predict_sigma1$cage <- predict(gam_nigeria, Nigeria3, parameter = 'sigma1', which = 'cage', type = 'link')
  predict_sigma1$edupartner <- predict(gam_nigeria, Nigeria3, parameter = 'sigma1', which = 'edupartner', type = 'link')
  predict_sigma1$mage <- predict(gam_nigeria, Nigeria3, parameter = 'sigma1', which = 'mage', type = 'link')
  predict_sigma1$mbmi <- predict(gam_nigeria, Nigeria3, parameter = 'sigma1', which = 'mbmi', type = 'link')
  
  predict_sigma2 <- vector('list')
  predict_sigma2$cage <- predict(gam_nigeria, Nigeria3, parameter = 'sigma2', which = 'cage', type = 'link')
  predict_sigma2$edupartner <- predict(gam_nigeria, Nigeria3, parameter = 'sigma2', which = 'edupartner', type = 'link')
  predict_sigma2$mage <- predict(gam_nigeria, Nigeria3, parameter = 'sigma2', which = 'mage', type = 'link')
  predict_sigma2$mbmi <- predict(gam_nigeria, Nigeria3, parameter = 'sigma2', which = 'mbmi', type = 'link')
  
  predict_rho <- vector('list')
  predict_rho$cage <- predict(gam_nigeria, Nigeria3, parameter = 'rho', which = 'cage', type = 'link')
  predict_rho$edupartner <- predict(gam_nigeria, Nigeria3, parameter = 'rho', which = 'edupartner', type = 'link')
  predict_rho$mage <- predict(gam_nigeria, Nigeria3, parameter = 'rho', which = 'mage', type = 'link')
  predict_rho$mbmi <- predict(gam_nigeria, Nigeria3, parameter = 'rho', which = 'mbmi', type = 'link')
  
  return(list(Coefficients = coef.nigeria, MSTOP = mstop.bivNorm, plot.data = plot.data, 
              predict_mu1 = predict_mu1, predict_mu2 = predict_mu2, predict_sigma1 = predict_sigma1,
              predict_sigma2 = predict_sigma2, predict_rho = predict_rho))

}

results = mclapply(1, sim, mc.cores = 10, mc.preschedule = FALSE)

save(results, file="Nigeria_Gaussian.RData")
