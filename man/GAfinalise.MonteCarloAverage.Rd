\name{GAfinalise.MonteCarloAverage}
\alias{GAfinalise.MonteCarloAverage}
\title{GAfinalise.MonteCarloAverage function}
\usage{
  \method{GAfinalise}{MonteCarloAverage} (F, ...)
}
\arguments{
  \item{F}{an object of class MonteCarloAverage}

  \item{...}{additional arguments}
}
\value{
  computes Monte Carlo averages
}
\description{
  Finalise a Monte Carlo averaging scheme. Divide the sum
  by the number of iterations.
}
\seealso{
  \link{MonteCarloAverage}, \link{setoutput},
  \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise},
  \link{GAreturnvalue}
}

