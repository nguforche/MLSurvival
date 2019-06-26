#' @export
Performance <- function(object, method="gbm"){

threshold <- object$threshold
predict.times <- object$predict.times 
pred <- object$resamples 

classifier <- c("gbm", "ranger", "glm", "glmnet", "xgbTree", "svmRadial")
method1 <- c("COX", "RSF", "Classifier")
names(method1) <- method1 

if(method%in%classifier) {mod1 = method; method = "Classifier"}  


switch(method1[method], 
Classifier = {
perf <- ddply(pred, .variables = "Resample", .fun = function(xx){
ddply(xx, .variables = "time.id", .fun = function(yy){
ix <- unique(yy$time.id)
p <- yy$Yes 
obs <- ifelse(yy$obs == "Yes", 1, 0) 
cbind(method = object$method, predict.times = predict.times[ix], Performance.measures(pred= p, obs= obs, threshold= threshold[ix]))
})
})
perf$time.id <- NULL
perf$Resample <- NULL 
rownames(perf) <- NULL 
},
COX = {
perf <- do.call(rbind, lapply(pred, function(xx){
xx <- data.frame(xx)
obs <- xx[, "class"] 
ts <- xx[, "time"] 
nme  <- colnames(xx)[-c(1,2)]
names(predict.times) <- nme
names(threshold) <- nme
do.call(rbind, lapply(nme, function(ii){
p <- xx[, ii]
y = obs
y[(predict.times[ii] < ts)] <- 0
AUC <- Performance.measures(pred=p, obs= y, threshold= threshold[ii])

cbind(method = "COX", predict.times = predict.times[ii],  AUC)
}))
}))
rownames(perf) <- NULL 
}, 

RF = {
perf <- do.call(rbind, lapply(pred, function(xx){
xx <- data.frame(xx)
obs <- xx[, "class"] 
ts <- xx[, "time"] 
nme  <- colnames(xx)[-c(1,2)]
names(predict.times) <- nme 
names(threshold) <- nme
do.call(rbind, lapply(nme, function(ii){
p <- xx[, ii]
y = obs
y[(predict.times[ii] < ts) & (obs == 1)] <- 0
AUC <- Performance.measures(pred=1-p, obs= y, threshold= threshold[ii])
cbind(predict.times = predict.times[ii],  AUC)
}))
}))
rownames(perf) <- NULL 
}, 
RSF = {
perf <- do.call(rbind, lapply(pred, function(xx){
xx <- data.frame(xx)
obs <- xx[, "class"] 
ts <- xx[, "time"] 
nme  <- colnames(xx)[-c(1,2)]
names(predict.times) <- nme 
names(threshold) <- nme
do.call(rbind, lapply(nme, function(ii){
p <- xx[, ii]
y = obs
y[(predict.times[ii] < ts)] <- 0
AUC <- Performance.measures(pred=p, obs= y, threshold= threshold[ii])
cbind(method = "RSF", predict.times = predict.times[ii],  AUC)
}))
}))
rownames(perf) <- NULL 
}, 
stop(paste0(x," not yet implemented"))
)
return(perf)

}






