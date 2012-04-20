\name{integerise.stppp}
\alias{integerise.stppp}
\title{integerise...}
\usage{\method{integerise}{stppp} (obj, ...)}
\description{integerise.stppp function}
\details{Function for converting the times and time limits of an stppp object into integer values. Do this before estimating mu(t), and hence
before creating the temporalAtRisk object. Not taking this step is possible in lgcp, but can cause minor complications connected with the scaling of mu(t).}
\value{The stppp object, but with integerised times.}
\arguments{\item{obj}{an stppp object}
\item{...}{additional arguments}}
