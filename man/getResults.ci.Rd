% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{getResults.ci}
\alias{getResults.ci}
\title{get Summary Results}
\usage{
getResults.ci(tab, alpha = 0.05)
}
\arguments{
\item{tab}{table with performance results}

\item{alpha}{confidence level}
}
\value{
data frame with summaries (confidence interavals are represented in brackets)
}
\description{
Takes a table of performance metrics, such as cross-validation results and compute summaries (mean and confidence interval) 
ready for publication.
}
