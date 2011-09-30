\name{target.and.grad}
\alias{target.and.grad}
\title{target...}
\usage{target.and.grad(at.mat.list, theta, varvals, nicounts.mat.list,
    timediffs, car.mat, mu, spatial.mat, temporal, lam.tilde, cin,
    M.ext, N.ext, gradtrunc, estgradtrunc=FALSE)}
\description{target.and.grad function}
\value{The log target and gradient}
\arguments{\item{at.mat.list}{list of matrices}
\item{theta}{parameter beta, see Brix and Diggle (2001)}
\item{varvals}{variances, pre-computed}
\item{nicounts.mat.list}{list of matrices of ni counts}
\item{timediffs}{vector of delta ts}
\item{car.mat}{matrix}
\item{mu}{mean parameter, see Brix and Diggle}
\item{spatial.mat}{fixed spatial component, output from function fftgrid}
\item{temporal}{fixed temporal component}
\item{lam.tilde}{transpose of Fourier transform on extended grid, output from function fftgrid}
\item{cin}{cellInside, output from function fftgrid}
\item{M.ext}{number of x 'gridlines' in extended space}
\item{N.ext}{number of y 'gridlines' in extended space}
\item{gradtrunc}{truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Set to zero for no gradient truncation (default).}
\item{estgradtrunc}{set to true if estimating gradient truncation, see ?computeGradtrunc}}
