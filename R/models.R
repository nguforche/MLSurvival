classifier.name <- function(method = "gbm"){
mod <- c("gbm", "ranger", "glm", "glmnet", "xgbTree", "svmRadial")
fun <- c("gbm", "ranger", "glm", "glmnet", "xgb.train", "ksvm")
names(fun) <- mod
fun[method]
}

Models <- function(x){
switch(x, 
gbm = {extra.para = GBM()}, 
ranger = {extra.para = RF()},
glm = {extra.para <- GLM()},
xgbTree = {xtra.para <- XGBoost()},
glmnet = {extra.para <- GLMNET()},
svmRadial = {extra.para <- SVM()},
stop(paste0(x," is not yet implemented"))
) 
}

