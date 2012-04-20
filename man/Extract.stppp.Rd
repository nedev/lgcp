\name{Extract.stppp}
\alias{"[.stppp"}
\alias{Extract.stppp}
\title{Extract.stppp function}
\usage{
"x[subset]"
}
\arguments{
  \item{x}{an object of class stppp}

  \item{subset}{the subset to extract}
}
\value{
extracts subset of an stppp object
}
\description{
extracting subsets of an stppp object.
}
\examples{
\dontrun{xyt <- lgcpSim()}
\dontrun{xyt}
\dontrun{xyt[xyt$t>0.5]}
}

