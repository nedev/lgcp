\name{computeGradtruncSpatioTemporal}
\alias{computeGradtruncSpatioTemporal}
\title{computeGradtruncSpatioTemporal function...}
\usage{computeGradtruncSpatioTemporal(nsims=100, scale=1, nis, mu, rootQeigs,
    invrootQeigs, spatial, temporal, bt, cellarea)}
\description{computeGradtruncSpatioTemporal function}
\details{\bold{Advanced use only.} A function to compute a gradient truncation parameter for 'spatial only' MALA via simulation. The function
requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.}
\value{gradient truncation parameter}
\seealso{\link{fftgrid}}
\arguments{\item{nsims}{The number of simulations to use in computation of gradient truncation.}
\item{scale}{multiplicative scaling constant, returned value is scale*max(gradient over simulations). Default scale is 1.}
\item{nis}{cell counts on the extended grid}
\item{mu}{parameter of latent field, mu}
\item{rootQeigs}{root of eigenvalues of precision matrix of latent field}
\item{invrootQeigs}{reciprocal root of eigenvalues of precision matrix of latent field}
\item{spatial}{spatial at risk interpolated onto grid of requisite size}
\item{temporal}{fitted temporal values}
\item{bt}{vectoer of variances b(delta t) in Brix and Diggls 2001}
\item{cellarea}{cell area}}
