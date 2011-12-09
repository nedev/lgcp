\name{as.owinlist.stapp}
\alias{as.owinlist.stapp}
\title{as...}
\usage{\method{as.owinlist}{stapp} (obj, dmin=0, check=TRUE, ...)}
\description{as.owinlist.stapp function}
\details{A function to create a list of owin objects from a stapp}
\value{a list of owin objects corresponding to the constituent Polygons objects}
\arguments{\item{obj}{an stapp object}
\item{dmin}{purpose is to simplify the SpatialPolygons. A numeric value giving the smallest permissible length of an edge. See ? simplify.owin}
\item{check}{whether or not to use spatstat functions to check the validity of SpatialPolygons objects}
\item{...}{additional arguments}}
