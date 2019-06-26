#' @title Function to implement proportional hazard model  
#' @description Train the  Cox model. Optionally, the regularize cox can also be trained based on the implementation in \code{glmnetUtils}.  
#' 
#' @param form survival formula 
#' @param dat  data frame 
#' @param trControl  list of control parameters:  
#'  \enumerate{
#'  \item maxit: number of iterations in \code{glmnet}  
#'  \item regularize: train regularize cox?
#'  \item nfolds :  number of folds in \code{glmnet} 
#'  \item lambda : numeric vector of lambda values in \code{glmnet} 
#'  \item alpha : numeric vector of alpha values in \code{cva.glmnet}    
#' }  
#' @return returns a list with items: 
#' \itemize{
#' \item{model: }{ object of class coxph}
#' \item{form: }{ survival formula} 
#' }


 COX <- function(form, dat, trControl = NULL){
  
  if(trControl$regularize){
    fit <- cva.glmnet(formula = form, data=dat, family="cox", alpha = trControl$alpha, 
                      lambda = trControl$lambda, maxit = trControl$maxit, nfolds = trControl$nfolds)
    mod <- fit$modlist
    alpha <- fit$alpha
    ## for each cv.glmnet model, select the lambda with the mimimum cross-validation loss
    err <- sapply(mod, function(xx){
      lambda <- xx$lambda
      lambda.min <- xx$lambda.min
      ix <- which(lambda == lambda.min)
      cvm <- xx$cvm
      cvm[ix]
    })
    afa <- alpha[which.min(err)]  ## select alpha corresponding to mimimum cross validation loss over all lambdas
    beta <- coef(fit, alpha = afa)[, 1]
    beta <- beta[beta != 0]
    rhs.vars <- names(beta)
    
    resp.vars <- all.vars(form[[2]]) ## 
    time = resp.vars[1]; status = resp.vars[2]
    
    form <- formula(paste("Surv(", time, ", ", status, ") ~", paste0(rhs.vars, collapse = " + ")))
    mod <- survival::coxph(form, data = dat, model = TRUE, init= beta,iter= 0, x=TRUE)
  } else 
    mod <- survival::coxph(formula=form, data = dat,model=TRUE)
  
  return(list(model = mod, form = form))
  }
  
  
  