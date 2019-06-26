#' @title Optimal Threshold 
#' @description Compute the optimal classifcation threshold based on the \code{optimal.thresholds} function in the 
#' \code{Presence.Absence} package 
#' @param prob predicted probabilities 
#' @param obs binary (0-1) ground truth 
#' @param opt.methods optimal threshold method. See \code{Presence.Absence} package 
#' @return optimal threshold 
#' @export 
opt.thresh <- function(prob, obs, opt.methods = 9){
  thresh = 0.5 
  if(length(unique(obs)) > 1){
    obs <- as.numeric(as.factor(obs))-1 
    SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, Predicted = prob)
    thresh <- optimal.thresholds(SIMDATA, threshold = 101, which.model = 1, opt.methods = opt.methods)
    thresh <- ifelse(length(thresh["Predicted"]) >= 1,as.numeric(thresh["Predicted"]), 0.5)
  }
  return(thresh)
}

#' @title Performance metrics  
#' @description Compute several performance metrics 
#' @param pred predicted probabilities 
#' @param obs binary (0-1) ground truth 
#' @param threshold optimal threshold.  
#' @return A data frame with performance metrics. 
#' @export 
Performance.measures <- function(pred, obs, threshold=NULL){
#  obs <- as.numeric(factor(obs))-1 
  nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
          "specificity", "specificity.sd")
  nme2 <- c("Pos Pred Value", "Neg Pred Value")#, "Balanced Accuracy")
  
  if(length(unique(obs)) <= 1){
    warning("No posivive cases in the test data, returning NA")
    res <- as.data.frame(matrix(NA, nrow = 1, ncol = 13))
    names(res) <- c(nme, "G.mean", "BER", nme2, "threshold")
    res$threshold <- 0.5 
  } else { 
  ## get best cut-off 
  if(is.null(threshold))
    threshold <- 0.5
  ### get these performance measures
  nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
          "specificity", "specificity.sd")
  xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)
  accuracy <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, nme]
  pred.prev <- predicted.prevalence(DATA=xx, threshold = threshold)[, c("Obs.Prevalence", "Predicted")]
  
  accuracy$G.mean <- sqrt(as.numeric(accuracy$sensitivity)*as.numeric(accuracy$specificity))
  accuracy$BER <- 1 - 0.5*(as.numeric(accuracy$sensitivity) + as.numeric(accuracy$specificity))  
  prevalence = as.numeric(pred.prev$Obs.Prevalence)
  obs <- factor(ifelse(obs == 1, "Yes", "No"), levels = c("Yes", "No"))
  pred <- factor(ifelse(pred >= threshold, "Yes", "No"), levels = c("Yes", "No"))
  cmx <- confusionMatrix(data=pred, reference=obs,  prevalence = prevalence)$byClass[nme2]
  res <- cbind.data.frame(accuracy, t(cmx),  threshold = threshold)
  }
  return(res)
}
#' @title Normalize data  
#' @description Normalize data to (0,1)
#' @param x data frame  
#' @return Normalized data with attributes min and max representing the min and max of each variable in \code{x}  
#' @export 
# normalize data to (0,1)
normalize <- function(x) { 
  x <- as.matrix(x)
  minAttr=apply(x, 2, min, na.rm=TRUE)
  maxAttr=apply(x, 2, max, na.rm=TRUE)
  x <- sweep(x, 2, minAttr, FUN="-") 
  x=sweep(x, 2,  maxAttr-minAttr, "/") 
  attr(x, 'min') = minAttr
  attr(x, 'max') = maxAttr
  return (x)
} 
#' @title Denormalized data   
#' @description Take the ouput of \code{normalize} and convert back to original scale.  
#' @param min minimum value of each variable in the original data. This value is stored as an attribute of \code{normalize} 
#' @param max maximum value of each variable in the original data. This value is stored as an attribute of \code{normalize}
#' @return Original un-normalized data  
#' @export 
# denormalize back to original scale
denormalize <- function (normalized, min, max) {
  if(dim(normalized)[2]!= length(min)) stop("length of min or max must equal number of columns of data ")
  nme <- colnames(normalized)
  if( !all(nme%in%names(min)) ) stop("colnames of data do not match names of min or max")
  sapply(nme, function(x)   normalized[, x] * (max[x]-min[x]) + min[x] )
}

#' @title Predict Cox   
#' @description Get predicted probabilities from the cox model for new data at different time points 
#' @param object trained cox model. Output of \code{train_Cox}  
#' @param newdata out of sample data 
#' @param times new time points 
#' @return predicted probabilities at eact time point in \code{times}
#' @export 
# prediction method for coxph model 
predictSurvProb_Cox <- function (object, newdata, times)
{
  survfit.object <- survival::survfit(object, newdata = newdata, se.fit = FALSE, conf.int = FALSE)
  if (is.null(attr(object$terms, "specials")$strata)) {
    inflated.pred <- summary(survfit.object, times = times)$surv
    if (is.null(inflated.pred)) {
      p = matrix(NA, ncol = length(times), nrow = NROW(newdata))
    }
    else {
      p <- t(inflated.pred)
      if ((beyond <- (length(times) - NCOL(p))) > 0)
        p <- cbind(p, matrix(NA, nrow = NROW(newdata),
                             ncol = beyond))
    }
  }
  else {
    inflated.pred <- summary(survfit.object, times = times)
    plist <- split(inflated.pred$surv, inflated.pred$strata)
    p <- do.call("rbind", lapply(plist, function(x) {
      beyond <- length(times) - length(x)
      c(x, rep(NA, beyond))
    }))
  }
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ",
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  p
}


#' @title Predict Ranger    
#' @description Get predicted probabilities from the ranger model for new data at different time points 
#' @param object trained ranger model. Output of \code{train_ranger}   
#' @param newdata out of sample data 
#' @param times new time points 
#' @param \dots further arguments passed to caret or other methods. 
#' @return predicted probabilities at eact time point in \code{times}
#' @export 
predictSurvProb_ranger <- function(object, newdata, times, ...){
  ptemp <- predict(object=object, data=newdata, ...)$survival   
  pos <- prodlim::sindex(jump.times=object$unique.death.times, eval.times=times)
  p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",
               NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ", NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  #colnames(p) = paste0("time_", times)
  p
}


sortImp <- function (object, top) 
{
    
    best <- "max"
    featureRank <- switch(best, max = rank(-apply(object, 1, max, na.rm = TRUE)), 
                                min = rank(apply(object, 1, min, na.rm = TRUE)), 
                                maxabs = rank(-apply(abs(object), 1, max, na.rm = TRUE)))

    tiedRanks <- as.numeric(names(table(featureRank)[table(featureRank) > 1]))
    if (length(tiedRanks) > 0) {
        for (i in seq(along = tiedRanks)) {
            tmp <- featureRank[featureRank == tiedRanks[i]]
            featureRank[featureRank == tiedRanks[i]] <- tmp + 
                runif(length(tmp), min = 0.001, max = 0.999)
        }
    }
    featureOrder <- order(featureRank)
    out <- object[featureOrder, , drop = FALSE]
    out <- out[1:top, , drop = FALSE]
    out
}
#' @title Plot variable importance     
#' @description Plot variable importance  
#' @param x data frame with variable importance    
#' @param top number of variables to plot 
#' @export 
VimPlot <- function (x, top = min(20, length(x$importance)), ...) 
{
    varSubset <- sortImp(x, top)
    plotObj <- stack(varSubset)
    if (dim(varSubset)[2] == 1) {
        plotObj <- varSubset
        names(plotObj) <- "values"
        plotObj$ind <- "Overall"
    }
    else plotObj <- stack(varSubset)
    plotObj$Var <- rep(rownames(varSubset), dim(varSubset)[2])
    plotObj$Var <- factor(plotObj$Var, levels = rev(rownames(varSubset)))
    if (dim(varSubset)[2] < 3) {
        if (dim(varSubset)[2] > 1) 
            plotObj <- plotObj[plotObj$ind == levels(plotObj$ind)[1], ]
#            out <- dotplot(Var ~ values, data = plotObj, as.table = TRUE,  xlab = "Importance", cex = 1.5, cex.labels = 1.5, pch=21, ...)
            par(mar=c(5,5.5,4,2) + 0.1)  
            out <- dotchart2(data=plotObj$values, labels=plotObj$Var, horizontal=TRUE, pch=19,col = "blue", 
			xlab="Importance", ylab="Variables", lty=1, lines=TRUE, dotsize = 1.2,
			cex = 1.2, cex.labels = 1.1, sort.=FALSE, add=FALSE, xaxis=TRUE, width.factor= 1,
			lcolor='gray', leavepar=FALSE, ...)
	    par()			    
#            out <- dotchart(plotObj$values, labels=plotObj$Var,  xlab = "Importance", ylab="Variables", lty=1, ...)
    }
    else {
        out <- dotplot(Var ~ values, data = plotObj, groups = plotObj$ind, 
            auto.key = list(columns = min(3, length(levels(plotObj$ind)))), 
            as.table = TRUE, xlab = "Importance", ...)
    }
    out
}
#' @title get Summary Results      
#' @description Takes a table of performance metrics, such as cross-validation results and compute summaries (mean and confidence interval) 
#' ready for publication.  
#' @param tab table with performance results
#' @param alpha confidence level
#' @return data frame with summaries (confidence interavals are represented in brackets)     
#' @export 
# get results with confidence intervals 
getResults.ci <- function(tab,  alpha = 0.05){
names(tab)[names(tab)%in%"n.dead"] <- "n.cases"
nme <- c("n", "n.cases", "PCC", "AUC",  "sensitivity",  "specificity", "G.mean", "BER", "Pos Pred Value")
nme2 <- c("PCC", "AUC",  "sensitivity",  "specificity", "G.mean", "BER", "Pos Pred Value")
tab$model <- factor(tab$model, levels = unique(tab$model))
tab$status <- factor(tab$status, levels = unique(tab$status))

mn = ddply(tab, .variables = c("model", "status"), numcolwise(mean), na.rm = TRUE)[, c("model", "status", nme)]
ci  <-   ddply(tab, .variables = c("model", "status"), numcolwise(quantile), 
probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)[,c("model", "status", nme)]
 
mn[, nme] <- format(round(mn[, nme], 2), nsmall = 2) 
ci[, nme] <- format(round(ci[, nme], 2), nsmall = 2) 

#### split by model, identify each status and merge 
df <- do.call(rbind.data.frame, lapply(unique(mn$model), function(mod){
xx <- mn[mn$model == mod, ]
yy <- ci[ci$model == mod, ]
do.call(rbind, lapply(unique(xx$status), function(ii){
x <- xx[xx$status == ii, ]
y <- yy[yy$status == ii, ]
n = x$n
cases = x$n.cases 
tb <- cbind.data.frame(mod, as.character(ii), n, cases, 
t(as.character(paste0(paste0(paste0(paste0(paste0(x[1, nme2], "(" ),  y[1, nme2]), ","), y[2, nme2]), ")"))))
colnames(tb) <- c("model", "status", nme)
tb 
}))
}))
df
}

#' @title get Summary Results 2      
#' @description Takes a table of performance metrics, such as cross-validation results and compute summaries (mean and confidence interval) 
#' ready for publication.  
#' @param tab table with performance results 
#' @param alpha confidence level
#' @param groups  variable in \code{tab} to group by
#' @return data frame with summaries (confidence interavals are represented in brackets)     
#' @export 
# get santized results for publication: include standard errors and confidence intervals 
# maximum grouping variables == 2 
## stats = the stats metrics to report 
getResults.ci2 <- function(tab,  alpha = 0.05, groups = c("model", "status"), 
        stats =c("PCC", "AUC",  "sensitivity",  "specificity", "G.mean", "BER", "Pos Pred Value")){
if(length(groups) > 2) stop("can only group by max of 2 variables")
sdError <- function(x) sd(x, na.rm = TRUE)/length(x)

tab[, groups] <- lapply(tab[, groups, drop = FALSE], function(xx) factor(xx, levels = unique(xx)))

mn = ddply(tab, .variables = groups, numcolwise(mean), na.rm = TRUE)
#sem <- ddply(tab, .variables = groups, numcolwise(sdError))
ci  <-   ddply(tab, .variables = groups, numcolwise(quantile), 
probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
mn[, stats] <- format(round(mn[, stats], 2), nsmall = 2) 
#sem[, stats] <- format(round(sem[, stats], 2), nsmall = 2) 
ci[, stats] <- format(round(ci[, stats], 2), nsmall = 2) 

#### split by grouping variables and sanitize
if(length(groups) == 2){  
df <- do.call(rbind.data.frame, lapply(unique(mn[, groups[1]]), function(mod){
xx <- mn[mn[, groups[1]] == mod, ]
#zz <- sem[sem[, groups[1]] == mod, ]
yy <- ci[ci[, groups[1]]  == mod, ]
do.call(rbind, lapply(unique(xx[, groups[2]]), function(ii){
x <- xx[xx[, groups[2]] == ii, ]
y <- yy[yy[, groups[2]] == ii, ]
#z = zz[zz[, groups[2]] == ii, ]
tb <- cbind.data.frame(mod, as.character(ii),
t(as.character(paste0(paste0(paste0(paste0(paste0(x[1, stats], "(" ),  y[1, stats]), ","), y[2, stats]), ")")))) 
#t(as.character(paste0(paste0(paste0(paste0(paste0(paste0(paste0(x[1,stats], " (" ), z[1,stats]), ","), y[1,stats]), "-"), 
#y[2,stats]), ")"))))
colnames(tb) <- c(groups, stats) 
tb 
}))
}))
} else {
df <- do.call(rbind.data.frame, lapply(unique(mn[, groups[1]]), function(mod){
x <- mn[mn[, groups[1]] == mod, ]
#z <- sem[sem[, groups[1]] == mod, ]
y <- ci[ci[, groups[1]]  == mod, ]
tb <- cbind.data.frame(mod,  
t(as.character(paste0(paste0(paste0(paste0(paste0(x[1, stats], "(" ),  y[1, stats]), ","), y[2, stats]), ")"))))
#t(as.character(paste0(paste0(paste0(paste0(paste0(paste0(paste0(x[1,stats], " (" ), z[1,stats]), ","), y[1,stats]), "-"), 
#y[2,stats]), ")"))))
colnames(tb) <- c(groups, stats)
tb 
}))
}
df
}






