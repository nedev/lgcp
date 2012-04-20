\name{tempRaster}
\alias{tempRaster}
\title{tempRaster function}
\usage{
  tempRaster(mcens, ncens)
}
\arguments{
  \item{mcens}{vector of equally-spaced coordinates of cell
  centroids in x-direction}

  \item{ncens}{vector of equally-spaced coordinates of cell
  centroids in y-direction}
}
\value{
  an empty raster object
}
\description{
  A function to create a temporary raster object from an
  x-y regular grid of cell centroids. Useful for projection
  from one raster to another.
}

