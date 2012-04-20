\name{mcmctrace.lgcpPredict}
\alias{mcmctrace.lgcpPredict}
\title{mcmctrace.lgcpPredict function}
\usage{
\method{mcmctrace}{lgcpPredict}(obj, ...)
}
\arguments{
  \item{obj}{an object of class lgcpPredict}

  \item{...}{additional arguments}
}
\value{
returns the saved MCMC chains in an object of class
\code{mcmcdiag}.
}
\description{
If \code{MCMCdiag} was positive when \code{lgcpPredict} was
called, then this retrieves information from the chains
stored.
}
\seealso{
\link{lgcpPredict}, \link{plot.mcmcdiag}
}

