\name{grid2spoly}
\alias{grid2spoly}
\title{grid2spoly function}
\usage{
  grid2spoly(xgrid, ygrid,
    proj4string = CRS(as.character(NA)))
}
\arguments{
  \item{xgrid}{vector of x centroids (equally spaced)}

  \item{ygrid}{vector of x centroids (equally spaced)}

  \item{proj4string}{proj 4 string: specify in the usual
  way}
}
\value{
  a SpatialPolygons object
}
\description{
  A function to convert a regular (x,y) grid of centroids
  into a SpatialPolygons object
}

