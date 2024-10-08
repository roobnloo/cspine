% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cspine-methods.R
\name{coef.cspine}
\alias{coef.cspine}
\title{Extract estimated coefficients from a 'cspine' object}
\usage{
\method{coef}{cspine}(fit)
}
\arguments{
\item{fit}{An object of class [cspine].}
}
\value{
A list with components 'beta', an array of precision matrix
 components, and 'gamma', the estimate mean component.
}
\description{
Extract estimated coefficients from a 'cspine' object
}
