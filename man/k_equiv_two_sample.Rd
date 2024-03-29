% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{k_equiv_two_sample}
\alias{k_equiv_two_sample}
\title{Calculate the factors for a two-sample acceptance test}
\usage{
k_equiv_two_sample(alpha, n, m)
}
\arguments{
\item{alpha}{the desired probability of Type 1 error}

\item{n}{the size of the qualification sample}

\item{m}{the size of the acceptance sample}
}
\value{
A vector of length 2 with the contents \code{c(k1, k2)}
}
\description{
Calculates the factors k1 and k2, which are used for setting acceptance
values for lot acceptance. These factors consider both the
size of the qualification sample (\code{n})
and the size of acceptance sample (\code{m}).
This test is detailed in a forthcoming paper.
}
\references{
Kloppenborg, S. (2023). Lot acceptance testing using sample mean and
extremum with finite qualification samples. Journal of Quality Technology,
https://doi.org/10.1080/00224065.2022.2147884
}
