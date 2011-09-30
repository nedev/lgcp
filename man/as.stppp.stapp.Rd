\name{as.stppp.stapp}
\alias{as.stppp.stapp}
\title{as...}
\usage{\method{as.stppp}{stapp} (obj, popden, n=100, dmin=0, check=TRUE, ...)}
\description{as.stppp.stapp function}
\details{A function to convert stapp objects to stppp objects for use in lgcpPredict. The regional counts in the stapp object are
assigned a random location within each areal region proportional to a population density (if that is available) else
the counts are distributed uniformly across the observation windows.}
\value{...}
\arguments{\item{obj}{an object of class stapp}
\item{popden}{a 'spatialAtRisk' of sub-class 'fromXYZ' object representing the population density, or for better results, lambda(s) can also
be used here. Cases are distributed across the spatial region according to popden. NULL by default, which has the effect of assigning counts uniformly.}
\item{n}{if popden is NULL, then this parameter controls the resolution of the uniform. Otherwise if popden is of class 'fromFunction', it controls the size of the imputation grid used for sampling. Default is 100.}
\item{dmin}{If any reginal counts are missing, then a set of polygonal 'holes' in the observation window will be computed for each. dmin is the parameter used to control the simplification of these holes (see ?simplify.owin). default is zero.}
\item{check}{logical. If any reginal counts are missing, then roughly speaking, check specifies whether to check the 'holes'.}
\item{...}{additional arguments}}
