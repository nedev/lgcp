\name{discreteWindow.lgcpPredict}
\alias{discreteWindow.lgcpPredict}
\title{discreteWindow.lgcpPredict function}
\usage{
\method{discreteWindow}{lgcpPredict}(obj, inclusion = "touching", ...)
}
\arguments{
  \item{obj}{an lgcpPredict object}

  \item{inclusion}{criterion for cells being included into
  observation window. Either 'touching' or 'centroid'. The
  former includes all cells that touch the observation
  window, the latter includes all cells whose centroids are
  inside the observation window.}

  \item{...}{additional arguments}
}
\value{
...
}
\description{
A function for extracting the FFT discrete window from an
lgcpPredict object.
}

