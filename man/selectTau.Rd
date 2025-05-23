% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{selectTau}
\alias{selectTau}
\title{Randomly select the chromosome}
\usage{
selectTau(N, prange, minDist, Pb, mmax, lmax)
}
\arguments{
\item{N}{The length of time series.}

\item{prange}{A list object containing the possible range for other
pre-defined model parameters, i.e. AR/MA order of ARMA models.}

\item{minDist}{The minimum length between two adjacent changepoints.}

\item{Pb}{Same as \code{Pchangepoint}, the probability that a changepoint has occurred.}

\item{mmax}{The maximum possible number of changepoints in the data set.}

\item{lmax}{The maximum possible length of the chromosome representation.}
}
\value{
A single changepoint configuration format as above.
}
\description{
Randomly select the changepoint configuration for population initialization.
The selected changepoint configuration represents a changepoint chromosome.
The first element of the chromosome represent the number of changepoints
and the last non-zero element always equal to the length of time series + 1.
}
