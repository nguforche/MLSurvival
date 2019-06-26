suppressMessages(require(MLSurvival))
suppressMessages(require(plyr))
suppressMessages(require(caret))
suppressMessages(require(PresenceAbsence))
suppressMessages(require(randomForestSRC))
suppressMessages(require(ranger))
suppressMessages(require(rms))
suppressMessages(require(pec))
suppressMessages(require(missForest))
suppressMessages(library(glmnet))
suppressMessages(require(glmnetUtils))


data(pbc, package = "survival") 
dat <- pbc

rhs.vars <- c("age", "sex", "ascites", "hepato", "spiders", "edema", "bili", "chol", "albumin",
              "alk.phos", "ast", "platelet", "protime", "stage")
dat$status[dat$status==1] <- 0
dat$status[dat$status==2] <- 1


## convert to factors 
col.class <- sapply(dat[, rhs.vars], class)
fac.vars <- names(col.class)[col.class%in%c("factor", "character")]
dat[,  fac.vars]  <- lapply(dat[, fac.vars, drop = FALSE], factor)


### impute data 
#cores <-  length(rhs.vars)-1 
#registerDoParallel(cores=cores)
dat.imp = dat 
imp <- missForest(dat[, rhs.vars], maxiter = 10, ntree = 50, parallelize = 'no')$ximp 
dat.imp[, rhs.vars] <- imp  

status = "status"
time = "time" 
dat.imp$time <- dat.imp$time*0.00273973
predict.times <- c(2, 5, 8) 

form <- formula(paste("Surv(", "time", ", ", status, ") ~", paste0(rhs.vars, collapse = " + ")))

### input parameters 
trControl = list(
  method = "cv", 
  number = 5, # number of cross-validations 
  summaryFunction = twoClassSummary, 
  classProbs = TRUE, 
  savePredictions = "final", 
  allowParallel = FALSE,
  
  ntrees = 100, # ranger number of trees 
  importance = "none", # ranger feature importance 
  tuneLength = 5, ## tuning paramer grid size 
  regularize = TRUE, ## regularize cox model? 
  lambda = seq(0.001,0.1,by = 0.001), ## vector of lasso penalty parameters
  alpha = seq(0, 1, len = 11)^3, ## vector of elasticnet penalty parameters 
  maxit = 10000, 
  nfolds = 3
)



method = c("COX", "RSF", "gbm", "ranger", "glm", "glmnet", "svmRadial", "xgbTree")

mod <- MLSurvival(form = form, dat = dat.imp, method = method, predict.times = predict.times, trControl = trControl,  
                  parallel = FALSE, dummy.vars = TRUE)
res <- do.call(rbind, lapply(mod, function(xx)  xx$perf.ave))
res

















