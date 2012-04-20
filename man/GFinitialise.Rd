\name{GFinitialise}
\alias{GFinitialise}
\title{GFinitialise function}
\usage{
  GFinitialise(F, ...)
}
\arguments{
  \item{F}{an object}

  \item{...}{additional arguments}
}
\value{
  method GFinitialise
}
\description{
  Generic function defining the the initialisation step for
  the \code{gridFunction} class of objects. The function is
  called invisibly within \code{MALAlgcp} and facilitates
  the dumping of data to disk
}
\seealso{
  \link{setoutput}, \link{GFupdate}, \link{GFfinalise},
  \link{GFreturnvalue}
}

