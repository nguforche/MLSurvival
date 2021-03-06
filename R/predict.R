#' @title Get Predictions for new data  
#' @param oject oject 
#' @param newdata new data to get predictions
#' @export 
predict.MLSurvival <- function(object,  method = "gbm", newdata, predict.times){

if (!inherits(object, "MLSurvival")) stop("Object must be a \"MLSurvival \"'")
if(!method%in%object$method) stop("Method not known!")

if(!is.null(object$dummy.model)){
if(length(object$factor.vars) > 0) newdata[,  object$factor.vars]  <- lapply(newdata[, object$factor.vars, drop = FALSE], factor)

pred <- predict(object$dummy.model, newdata = newdata)


other.vars <- setdiff(object$rhs.vars, colnames(pred))
newdata <- cbind(newdata[, other.vars], pred)
}


classifier <- c("gbm", "ranger", "glm", "glmnet", "xgbTree", "svmRadial")
method1 <- c("COX", "RSF", "Classifier")
names(method1) <- method1 
mod <- method 
if(mod%in%classifier) {mod = "Classifier"}     

switch(method1[mod], 
COX = {
pred <- 1-predictSurvProb_Cox(object =object$MLsurvival$COX$model$finalModel, newdata = newdata, times = predict.times)
}, 
RSF = {pred <- 1-predictSurvProb_ranger(object = object$MLsurvival$RSF$model$finalModel, newdata = newdata, times= predict.times)
}, 
Classifier = {
if(!all(predict.times%in%object$predict.times)) stop("Predict times in model differ from provided times") 
pred <- do.call(cbind, lapply(object$MLsurvival[[method]]$model$finalModel, function(xx) predict(xx, newdata = newdata, type = "prob")[,"Yes"]))
colnames(pred) <- names(object)	
}, 
stop(paste0(x," not yet implemented"))
)
return(pred)
}



