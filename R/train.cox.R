#' @title Cox proportional hazard model  
#' @description Train the Cox model through cross-validation and select the optimal survival classification threshold. 
#' A regularized Cox approach which performs feature selection is also implemeted. For regularize cox, the optimal set of variables is 
#' selected through cross-validation and used to train the final model on the complete data  
#' @name train.cox  
#' @param form survival formula 
#' @param dat  data frame 
#' @param predict.times survival prediction times 
#' @param trControl  list of control parameters:  
#'  \enumerate{
#'  \item number: number of cross-validations 
#'  \item regularize: train regularize cox?   
#' }  
#' @param parallel run cross-validation in parallel? Uses mclapply which works only on linux 
#' @param tuneLength  same as tuneLength in the caret package 
#' @param \dots further arguments passed to caret or other methods.  
#' @return returns a list with items: 
#' \itemize{
#' \item{finalModel: }{ final model trained on the complete data (dat) using optimal tuning paramters}
#' \item{fitted: }{ predictions on complete data (dat)} 
#' \item{threshold: }{ optimal classification threshold}
#' \item{resamples: }{ cross-validation results: predictions on resampled data } 
#' \item{predict.times: }{ survival prediction times}
#' \item{bestTune: }{ optimal tuning parameters} 
#' }
NULL 
#' @rdname cox  
#' @export
train_cox = function(form, dat, predict.times, trControl = NULL, parallel = FALSE, mc.cores = 2, 
                        seed = 123, ...){				
### 
predict.times <- sort(predict.times)
  if(parallel) {
    pfun <-  get("mclapply")
    set.seed(seed, "L'Ecuyer") 
  } else {
    set.seed(seed)
    pfun = get("lapply")
  }

if(is.null(trControl))
  trControl = list(number = 5, regularize = FALSE, lambda = seq(0.001,0.1,by = 0.001), maxit = 10000, 
                   alpha = 1, nfolds = 3)

resp.vars <- all.vars(form[[2]]) ## 
time = resp.vars[1]; status = resp.vars[2]

ix.cv <- createFolds(dat[, status], k= trControl$number, list=FALSE)
## cv 
resamples <- lapply(1:trControl$number, function(kk){      
    dat.trn  <- dat[kk!=ix.cv, ]
    dat.tst <-  dat[kk==ix.cv, ] 
    obs = dat.trn[, status]
    ts = dat.trn[, time]
    Obs = dat.tst[, status]
    Ts = dat.tst[, time]
#model <- survival::coxph(form, data = dat.trn, ...)
    
cox.mod <- COX(form, dat.trn, trControl = trControl)
model =cox.mod$model 
form = cox.mod$form 

pp <- predictSurvProb_Cox(object = model, newdata = dat.trn, times= predict.times)
PP <- predictSurvProb_Cox(object = model, newdata = dat.tst, times= predict.times)

perf <- lapply(1:ncol(pp), function(xx) {
y = dat.trn[, status]
y[(predict.times[xx] < ts)] <- 0

Y = dat.tst[, status]
Y[(predict.times[xx] < Ts)] <- 0

thresh <- opt.thresh(prob=1-pp[, xx], obs=y)
AUC = Performance.measures(pred=1-PP[, xx], obs= Y, threshold= thresh)$AUC

list(AUC = AUC, thresh = thresh)
} )

AUC <- sapply(perf, function(xx) xx$AUC)
threshold <- sapply(perf, function(xx) xx$thresh) 
colnames(PP) <- paste0("T", 1:length(predict.times))
prob <- cbind(class = dat.tst[, status],  time = dat.tst[, time], 1-PP)
rhs.vars <- all.vars(form[[3]])                                                     
res <- list(prob = prob, AUC=AUC, threshold = threshold, rhs.vars = rhs.vars)
return(res)
})

prob <- lapply(resamples, function(xx) xx$prob)
names(prob) <- paste0("cv", 1:trControl$number)

AUC <- do.call(rbind, lapply(resamples, function(xx) xx$AUC))
which.max1 <- function(x) {if(all(is.na(x))) return(NA) else return(which.max(x))}
ix <- as.numeric(apply(AUC, 2, which.max1))
threshold <- do.call(rbind, lapply(resamples, function(xx) xx$threshold))
threshold <- sapply(1:length(predict.times), function(ii) {ifelse(is.na(ix[ii]),  0.5, threshold[ix[ii], ii])})
names(threshold) <- paste0("T", 1:length(predict.times))

ix <- which.max(apply(AUC, 1, mean))
rhs.vars <- lapply(resamples, function(xx) xx$rhs.vars)[[ix]]

### train final model on all data 
form <- formula(paste("Surv(", time, ", ", status, ") ~", paste0(rhs.vars, collapse = " + ")))
trControl$regularize = FALSE
finalModel <- COX(form, dat, trControl = trControl) 

fitted <- 1-predictSurvProb_Cox(object =finalModel$model, newdata = dat, times= predict.times)
colnames(fitted) <- paste0("T", 1:length(predict.times))
fitted <- cbind(obs = dat[, status], fitted)

return(list(finalModel = finalModel, fitted = fitted, threshold = threshold, resamples = prob, 
predict.times =  predict.times, bestTune = NULL))
}







