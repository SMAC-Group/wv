% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covwv.R
\name{wccv_get_y}
\alias{wccv_get_y}
\title{Mapping to log10 scale}
\usage{
wccv_get_y(x, tick_y_min, tick_y_step)
}
\arguments{
\item{x}{A \code{vector} with dimensions J x 1.}

\item{tick_y_min}{A \code{negtive integer} the minimum power of 10, which corresponds to the smallest scale on y-axis.}

\item{tick_y_step}{An \code{integer} indicating the increment of the sequence.}
}
\value{
A \code{field<vec>} that contains values in log10 scale.
}
\description{
Map x to the value in log10 scale
}
\details{
\code{tick_y_min} is usually chosen as \eqn{floor(min(log10(abs(x))))}
}
\examples{
x = 2^(-1:-9)
y.min = floor(min(log10(abs(x))))
y.step = 2
wccv_get_y(x, y.min, y.step)
}
\author{
James Balamuta and Justin Lee
}
