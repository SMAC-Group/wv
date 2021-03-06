% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dwt.R
\name{summary.dwt}
\alias{summary.dwt}
\title{Summary Discrete Wavelet Transform}
\usage{
\method{summary}{dwt}(object, ...)
}
\arguments{
\item{object}{A \code{dwt} object}

\item{...}{additional arguments affecting the summary produced.}
}
\value{
Prints the dwt matrix decomposition
}
\description{
Prints DWT object in a concise format
}
\examples{
set.seed(999)
x = rnorm(2^8)
summary(dwt(x))
}
\author{
Nathanael Claussen and Justin Lee
}
\keyword{internal}
