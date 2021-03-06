% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resultsSeeds.R
\name{hiddenInputs}
\alias{hiddenInputs}
\alias{hiddenInputs,list,numeric-method}
\alias{hiddenInputs,list,missing-method}
\alias{hiddenInputs,resultsSeeds,missing-method}
\title{Get the estimated hidden inputs}
\usage{
hiddenInputs(resultsSeeds, ind)

\S4method{hiddenInputs}{list,numeric}(resultsSeeds, ind)

\S4method{hiddenInputs}{list,missing}(resultsSeeds, ind)

\S4method{hiddenInputs}{resultsSeeds,missing}(resultsSeeds, ind)
}
\arguments{
\item{resultsSeeds}{A object of the class 'resultsSeeds', which is returned from the algorithms.}

\item{ind}{A numeric indicating the index of a 'resultsSeeds'-Object in a list. If not set the last listed object will be used.}
}
\value{
Dataframe containing the estimated hidden inputs
}
\description{
Get the estimated hidden inputs
}
\examples{

data(uvb_res)

hiddenInputs(res[[2]])

}
