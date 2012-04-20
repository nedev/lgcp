\name{stapp.list}
\alias{stapp.list}
\title{stapp...}
\usage{\method{stapp}{list} (obj, ...)}
\description{stapp.list function}
\details{A wrapper function for stapp.SpatialPolygonsDataFrame

Construct a space-time aggregated point-process (stapp) object from a list object. The first element of the list 
should be a SpatialPolygonsDataFrame, the second element of the list a counts matrix, the third element of the list a
vector of times, the fourth element a vector giving the bounds of the temporal observation window and the fifth element
a spatstat owin object giving the spatial observation window.}
\value{an object of class stapp}
\arguments{\item{obj}{an list object as described above, see ?stapp.SpatialPolygonsDataFrame for further details on the requirements of the list}
\item{...}{additional arguments}}
