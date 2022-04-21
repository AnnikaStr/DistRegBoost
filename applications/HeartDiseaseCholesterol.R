#########################################
library("gamboostLSS")
source('bivariateBernoulli.R')
#########################################

sim <- function(seed){
  
  load("genotype_data_chol_heart.RData")
  ukbb_i25_cholesterol <- read.delim("phenotypes_subset_chol_heart.txt")
  
  ukbb_i25_cholesterol$cholesterol_0 <- 0
  ukbb_i25_cholesterol$cholesterol_0[which(ukbb_i25_cholesterol$cholesterol > 6.19)]  <- 1
  
  ukbb_i25_cholesterol[,c('FID','IID', 'PC1','PC2', 'PC3', 'PC4','cholesterol','sex','age','bmi')]  <- list(NULL)
  genotype_data[,c('FID', 'IID')] <- list(NULL)
  genotype_data = data.frame(genotype_data)

  set.seed(seed)
  ALL <- data.frame(cbind(ukbb_i25_cholesterol, genotype_data))
  samp <- sample(1:30000, 10000)
  ALL_val <- ALL[samp,]
  ALL_train <- ALL[-samp,]
  weights_train <- rep(1,30000)  
  weights_train[samp] <- 0 
  
  perc <- table(ALL$I25_cases)

  glm_BV <- glmboostLSS(cbind(I25_cases,cholesterol_0) ~ ., data = ALL, families = BernoulliBV(), method = 'noncyclic', weights = weights_train,
                        control = boost_control(mstop = 20000, nu = list(mu1 = 0.1, mu2 = 0.1, or = 0.1), trace = T, risk = 'oobag'))
  
  MSTOP <- which.min(risk(glm_BV,merge = T))
  risk_b <- risk(glm_BV, merge = T)
  
  rm(bivBern)
  ALL <- ALL[weights_train == 1, ]
  
  glm_BV <- glmboostLSS(cbind(I25_cases,cholesterol_0) ~ ., data = ALL, families = BernoulliBV(), method = 'noncyclic', 
                        control = boost_control(mstop = MSTOP, nu = list(mu1 = 0.1, mu2 = 0.1, or = 0.1), trace= T))
  
  risk_a <- risk(glm_BV, merge = T)
  
  mstop.bivBern <-  vector('list')
  mstop.bivBern$mstop <- MSTOP
  mstop.bivBern$mu1 <- glm_BV$mu1$mstop()
  mstop.bivBern$mu2 <- glm_BV$mu2$mstop()
  mstop.bivBern$or <- glm_BV$or$mstop()
  
  coef.bivBern <- coef(glm_BV)

  return(list(Coefficients = coef.bivBern, fit_biv = fit_biv,risk_after = risk_a, risk_before = risk_b, mstop = mstop.bivBern, perc = perc, time.taken = time.taken))
  
}  

results = mclapply(1, sim, mc.cores = 10, mc.preschedule = FALSE)

save(results, file="HeartDiseaseCholesterol_Bernoulli.RData")
