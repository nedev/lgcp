\name{stapp.SpatialPolygonsDataFrame}
\alias{stapp.SpatialPolygonsDataFrame}
\title{stapp...}
\usage{\method{stapp}{SpatialPolygonsDataFrame} (obj, counts, t, tlim, window,
    ...)}
\description{stapp.SpatialPolygonsDataFrame function}
\details{Construct a space-time aggregated point-process (stapp) object from a SpatialPolygonsDataFrame (along with some other info)}
\value{an object of class stapp}
\arguments{\item{obj}{an SpatialPolygonsDataFrame object}
\item{counts}{a (length(t) by N) matrix containing aggregated case counts for each of the geographical regions defined by the SpatialPolygonsDataFrame, where N is the number of regions}
\item{t}{vector of times, for each element of t there should correspond a column in the matrix 'counts'}
\item{tlim}{vector giving the upper and lower bounds of the temporal observation window}
\item{window}{the observation window, of class owin, see ?owin}
\item{...}{additional arguments}}
