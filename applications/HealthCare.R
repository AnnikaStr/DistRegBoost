#########################################
library("gamboostLSS")
source('bivariatePoisson.R')
#########################################


sim <- function(seed){
  
  load("ex3.health.rda")
  ex3.health$sex <- as.factor(ex3.health$sex)
  
  set.seed(seed)
  
  form.lambda1 <- as.formula(cbind(doctorco,prescrib) ~ bols(sex) + bbs(age) + bbs(income))
  form.lambda2 <- as.formula(cbind(doctorco,prescrib) ~ bols(sex) + bbs(age) + bbs(income))
  form.lambda3 <- as.formula(cbind(doctorco,prescrib) ~ bols(sex) + bbs(age) + bbs(income))

  form <- list(lambda1 = form.lambda1,
               lambda2 = form.lambda2,
               lambda3 = form.lambda3)
  
  glm.health <- gamboostLSS(form, data = ex3.health, families = PoissonBV(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL), method = 'noncyclic',
                            control = boost_control(mstop = 50))

  coef(glm.health)
  
  cvr.health <- cvrisk(glm.health, folds = cv(model.weights(glm.health)), grid = 1:2000, papply = mclapply)
  MSTOP <- mstop(cvr.health)
  
  glm.health[MSTOP]
  
  mstop.bivPois <-  vector('list')
  mstop.bivPois$mstop <- MSTOP
  mstop.bivPois$lambda1 <- glm.health$lambda1$mstop()
  mstop.bivPois$lambda2 <- glm.health$lambda2$mstop()
  mstop.bivPois$lambda3 <- glm.health$lambda3$mstop()

  coef.bivPois <- coef(glm.health)
  
  fit_lambda1 <- lapply(1:3, function(i)fitted(glm.health$lambda1, which = i))
  fit_lambda2 <- lapply(1:3, function(i)fitted(glm.health$lambda2, which = i))
  fit_lambda3 <- lapply(1:3, function(i)fitted(glm.health$lambda3, which = i))
  
  fit_biv <- list(fit_lambda1,fit_lambda2,fit_lambda3)
  
  predict_lambda2 <- vector('list')
  predict_lambda2$age <- predict(glm.health, ex3.health, parameter = 'lambda2', which = 'age', type = 'link')
  predict_lambda2$income <- predict(glm.health, ex3.health, parameter = 'lambda2', which = 'income', type = 'link')

  predict_lambda1 <- vector('list')
  predict_lambda1$age <- predict(glm.health, ex3.health, parameter = 'lambda1', which = 'age', type = 'link')
  predict_lambda1$income <- predict(glm.health, ex3.health, parameter = 'lambda1', which = 'income', type = 'link')
  
  predict_lambda3 <- vector('list')
  predict_lambda3$age <- predict(glm.health, ex3.health, parameter = 'lambda3', which = 'age', type = 'link')
  predict_lambda3$income <- predict(glm.health, ex3.health, parameter = 'lambda3', which = 'income', type = 'link')
  
  
  return(list(Coefficients = coef.bivPois, fit_biv = fit_biv, mstop = mstop.bivPois, predict_lambda1 =predict_lambda3, predict_lambda2 = predict_lambda2,
              predict_lambda3 = predict_lambda3))
  
}  

results = mclapply(1, sim, mc.cores = 1, mc.preschedule = FALSE)

save(results, file="Health_Pois.RData")

