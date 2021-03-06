% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train.RSF.R
\name{RSF}
\alias{RSF}
\alias{train_RSF}
\title{RSF: Random Survival Forest}
\usage{
train_RSF(form, dat, newdata = NULL, predict.times, trControl = NULL,
  seed = 123, parallel = FALSE, mc.cores = 2, ...)
}
\arguments{
\item{form}{survival formula}

\item{dat}{data frame}

\item{predict.times}{survival prediction times}

\item{trControl}{list of control parameters:
 \enumerate{
 \item ntrees: mumber of trees 
 \item number: number of cross-validations 
 \item tuneLength: tuning paramer grid size  
 \item importance: ranger variable importance  
}}

\item{parallel}{run cross-validation in parallel? Uses mclapply which works only on linux}

\item{\dots}{further arguments passed to caret or other methods.}

\item{tuneLength}{same as tuneLength in the caret package}
}
\value{
returns a list with items: 
\itemize{
\item{finalModel: }{ final model trained on the complete data (dat) using optimal tuning paramters}
\item{fitted: }{ predictions on complete data (dat)} 
\item{threshold: }{ optimal classification threshold}
\item{resamples: }{ cross-validation results: predictions on resampled data } 
\item{predict.times: }{ survival prediction times}
\item{bestTune: }{ optimal tuning parameters} 
}
}
\description{
Train the random survival forest through the ranger package. The optimal RSF tuning parameters: 
min.node.size,mtry, and splitrule can be selected through grid search.
}
