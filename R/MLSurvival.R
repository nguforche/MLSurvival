#' @title MLSurvival   
#' @description Train and evelaute machine learning survival and classification models for time to event data
#
#' @name MLSurvival    
#' @param form survival formula 
#' @param dat  data frame 
#' @param method character verctor of machine learning algorithms. Implemented algorithms 
#'  \enumerate{
#'  \item glm  logistic regression 
#'  \item glmnet elastic net 
#'  \item gbm gradient boosting machine 
#'  \item ranger random forest 
#'  \item svmRadial  support vector machine with radial basis kernel 
#'  \item xgbTree extreem boosting machine  
#' }  
#' @param predict.times numeric vector containing the survival prediction times  
#' @param trControl  list of control parameters for caret and the ranger models
#' @param parallel run cross-validation in parallel? 
#' @param dummy.vars create  dummy variables/model.matrix  
#' @param mc.cores number of cores 
#' @param seed random seed
#' @param perf get performance metrics ? 
#' @param \dots further arguments passed to caret or other methods. 
#' @return returns a list with items: 
#' \itemize{
#' \item{model:}{ trained survival model}
#' \item{perf:}{ performance of models at each survival prediction time: PCC, AUC, sensitivity, specificity, g-mean etc. } 
#' \item{perf.ave:}{ average of perf with confidence intervals} 
#' }

# 
NULL 
#' @rdname MLSurvival
#' @export
MLSurvival <- function(x, ...) UseMethod("MLSurvival")
#' @rdname MLSurvival
#' @export
MLSurvival.default <-  function(x, y, method, predict.times, trControl, parallel = FALSE, dummy.vars = TRUE, mc.cores = 2, 
                        seed = 123, perf = TRUE, ...){				


}
#' @rdname MLSurvival
#' @export
MLSurvival.formula <-  function(form, dat, newdata = NULL, method, predict.times, trControl, parallel = FALSE, dummy.vars = TRUE, mc.cores = 2, 
                        seed = 123, perf = TRUE,  ...){				
  dummy.model <- NULL
  perf.ave <- NULL
  perf.res <- NULL 
  fac.vars <- NULL 
  rhs.vars <- NULL 
 
  if(dummy.vars){
  resp.vars <- all.vars(form[[2]]) ## 
  time = resp.vars[1]; status = resp.vars[2]
  rhs.vars <- all.vars(form[[3]])
  
  col.class <- sapply(dat[, rhs.vars], class)
  fac.vars <- names(col.class)[col.class%in%c("factor", "character")]
  other.vars <- setdiff(names(dat), fac.vars)
  if(length(fac.vars) > 0){
    dat[,  fac.vars]  <- lapply(dat[, fac.vars, drop = FALSE], factor)
    form <- as.formula(paste0("~", fac.vars))
    dummy.model <- dummyVars(formula = form, data = dat, fullRank = TRUE) ### model to convert factors to dummy variables 
    pred <- predict(dummy.model, newdata = dat)
    vars <- setdiff(rhs.vars, fac.vars)
    dat <- cbind(dat[, other.vars], pred)
    rhs.vars <- unique(c(vars, colnames(pred)))
    form <- formula(paste("Surv(", time, ", ", status, ") ~", paste0(rhs.vars, collapse = " + ")))
  }
  }

classifier <- c("gbm", "ranger", "glm", "glmnet", "xgbTree", "svmRadial")
method1 <- c("COX", "RSF", "Classifier")
names(method1) <- method1 

res <- lapply(method, function(mod){
  mod1 = mod;
  if(mod%in%classifier) {mod = "Classifier"}     
  
switch(method1[mod], 
COX = {model <- train_cox(form=form, dat=dat, predict.times=predict.times, trControl=trControl, 
                    parallel = parallel, mc.cores = mc.cores, seed = seed, ...)
#cat("done Cox \n")

}, 
RSF = {model <- train_RSF(form=form, dat=dat, predict.times=predict.times, trControl=trControl, 
                    parallel = parallel, mc.cores = mc.cores, seed = seed, ...)
                                        
#cat("done RSF \n")
}, 
Classifier = {model <- train_classifier(form=form, dat=dat, newdata = newdata, method = mod1, predict.times=predict.times, trControl=trControl, 
                    parallel = parallel, mc.cores = mc.cores, seed = seed, ...)	
#cat("done classifier \n")
}, 
stop(paste0(x," not yet implemented"))
)

if(perf){
perf.res <- Performance(object=model, method=mod)
rownames(perf.res) <- NULL 
perf.ave <- getResults.ci2(perf.res, alpha = 0.05, groups = c("method", "predict.times") )
rownames(perf.ave) <- NULL 
}
cat("done", mod1, "\n")

return(list(model = model, perf = perf.res, perf.ave = perf.ave))
})

names(res) <- method 
Res <- list(MLsurvival = res, dummy.model=dummy.model, factor.vars = fac.vars, formula = form, rhs.vars = rhs.vars, method = method)
class(Res) <- "MLSurvival"
Res 
}  

