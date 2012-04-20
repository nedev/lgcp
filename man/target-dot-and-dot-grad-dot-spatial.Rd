\name{target.and.grad.spatial}
\alias{target.and.grad.spatial}
\title{target...}
\usage{target.and.grad.spatial(Gamma, nis, cellarea, rootQeigs, invrootQeigs,
    mu, spatial, logspat, scaleconst, gradtrunc)}
\description{target.and.grad.spatial function}
\details{A function to compute the target and gradient for 'spatial only' MALA}
\value{the back-transformed Y, its exponential, the log-target and gradient for use in MALAlgcpSpatial}
\arguments{\item{Gamma}{current state of the chain, Gamma}
\item{nis}{matrix of cell counts}
\item{cellarea}{area of cells, a positive number}
\item{rootQeigs}{square root of the eigenvectors of the precision matrix}
\item{invrootQeigs}{inverse square root of the eigenvectors of the precision matrix}
\item{mu}{parameter of the latent Gaussian field}
\item{spatial}{spatial at risk function, lambda, interpolated onto correct grid}
\item{logspat}{log of spatial at risk function, lambda*scaleconst, interpolated onto correct grid}
\item{scaleconst}{the expected number of cases}
\item{gradtrunc}{gradient truncation parameter}}
