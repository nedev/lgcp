\name{neattable}
\alias{neattable}
\title{neattable function...}
\usage{neattable(mat, indent=0)}
\description{neattable function}
\details{Function to print right-aligned tables to the console.}
\value{prints to screen with specified indent}
\arguments{\item{mat}{a numeric or character matrix object}
\item{indent}{indent}}
\examples{mat <- rbind(c("one","two","three"),matrix(round(runif(9),3),3,3))
neattable(mat)}
