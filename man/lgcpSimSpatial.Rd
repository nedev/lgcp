\name{lgcpSimSpatial}
\alias{lgcpSimSpatial}
\title{lgcpSimSpatial function...}
\usage{lgcpSimSpatial(owin, spatial.intensity, expectednumcases=100,
    cellwidth=0.05, model.parameters=lgcppars(sigma = 2, phi = 0.2),
    spatial.covmodel="exponential", covpars=c(), ext=2, plot=FALSE)}
\description{lgcpSimSpatial function}
\details{A function to simulate from a log gaussian process}
\value{a ppp object containing the data}
\arguments{\item{owin}{observation window}
\item{spatial.intensity}{an object that can be coerced to one of class spatialAtRisk}
\item{expectednumcases}{the expected number of cases}
\item{cellwidth}{width of cells}
\item{model.parameters}{parameters of model, see ?lgcppars. Only set sigma and phi for spatial model.}
\item{spatial.covmodel}{spatial covariance function, default is exponential, see ?CovarianceFct}
\item{covpars}{vector of additional parameters for spatial covariance function, in order they appear in chosen model in ?CovarianceFct}
\item{ext}{how much to extend the parameter space by. Default is 2.}
\item{plot}{logical, whether to plot the latent field.}}
