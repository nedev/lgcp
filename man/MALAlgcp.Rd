\name{MALAlgcp}
\alias{MALAlgcp}
\title{MALAlgcp function...}
\usage{MALAlgcp(mcmcloop, inits, adaptivescheme, THETA, vars, nis, tdiffs,
    cellArea.mat, MU, spatialvals, temporal.fitted, LAM.tildefft,
    cellInside, MCMCdiag, gradtrunc, gridfun, gridav)}
\description{MALAlgcp function}
\details{\bold{For Internal/Advanced Use Only.}}
\value{mean Gaussian field from MCMC, exceedances, output of MCMC chain  etc (if specified)}
\arguments{\item{mcmcloop}{an MCMC loop controller (see ?mcmcLoop)}
\item{inits}{optional initial point for MALA - useful if updating on a daily basis}
\item{adaptivescheme}{the type of adaptive mcmc to use, see ?constanth (constant h) or ?andrieuthomsh (adaptive MCMC of Andrieu and Thoms (2008))}
\item{THETA}{corresponds to parameter beta in Brix & Diggle (2001)}
\item{vars}{pre-computed coefficient (1-exp(-2*beta*delta_t)) in G(t)}
\item{nis}{these are the cell counts at each time lag}
\item{tdiffs}{time differences (note first entry is an NA)}
\item{cellArea.mat}{matrix of cell areas}
\item{MU}{estimate of parameter mu in Brix & Diggle (2001)}
\item{spatialvals}{fixed spatial component}
\item{temporal.fitted}{fixed temporal component}
\item{LAM.tildefft}{transpose of Fourier transform on extended grid, output from function fftgrid}
\item{cellInside}{(M-1)x(N-1) matrix of TRUE/FALSE indicating whether cell is inside observation window}
\item{MCMCdiag}{non-negative integer, if greater than zero saves information from the MCMC chain}
\item{gradtrunc}{truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Set to zero for no gradient truncation (default).}
\item{gridfun}{chosen grid function, see ?nullfunction (default) or ?dump2dir for example}
\item{gridav}{chosen grid average, see ?nullAverage (default) or ?MonteCarloAverage for example}}
