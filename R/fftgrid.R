###
# FFT grid handling functions
###

##' fftgrid function
##'
##' \bold{Advanced use only.} Computes various quantities for use in \code{lgcpPredict},
##' \code{lgcpSim} .
##' 
##' @param xyt object of class stppp
##' @param M number of centroids in x-direction 
##' @param N number of centroids in y-direction
##' @param spatial an object of class spatialAtRisk
##' @param sigma scaling paramter for spatial covariance function, see Brix and Diggle (2001) 
##' @param phi scaling paramter for spatial covariance function, see Brix and Diggle (2001)
##' @param model correlation type see ?CovarianceFct
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @return fft objects for use in MALA
##' @examples
##' xyt <- stppp(ppp(),t=c(),tlim=c(0,1))
##' sar <- spatialAtRisk(function(x,y){return(1)},warn=FALSE)
##' grid <- fftgrid(xyt,64,64,sar,1,0.01,"exponential",covpars=c())
##' @export


fftgrid <- function(xyt,M,N,spatial,sigma,phi,model,covpars){

    verifyclass(xyt,"stppp")
    verifyclass(spatial,"spatialAtRisk")
    
    study.region <- xyt$window
	
	## DEFINE LATTICE & CENTROIDS ##
		
	del1 <- (study.region$xrange[2]-study.region$xrange[1])/M
	del2 <- (study.region$yrange[2]-study.region$yrange[1])/N 
	
	M.ext <- 2*M 	
	N.ext <- 2*N ##
	
	mcens <- study.region$xrange[1]+.5*del1+(0:(M.ext-1))*del1
	ncens <- study.region$yrange[1]+.5*del2+(0:(N.ext-1))*del2	
	
	## REQUIRED SIMULATION QUANTITIES ##
	
	cellArea.mat <- matrix(0,M.ext,N.ext)
	cellArea.mat[1:M,1:N] <- del1*del2
	
	cellInside <- inside.owin(x=sort(rep(mcens,N.ext)),y=rep(ncens,M.ext),w=study.region)
	cellInside <- matrix(as.logical(cellInside),M.ext,N.ext,byrow=T)[1:M,1:N]
	cellInsideBIG <- matrix(0,M.ext,N.ext)
	cellInsideBIG[1:M,1:N][cellInside] <- 1
	
	## OBTAIN SPATIAL VALS ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatialvals <- fftinterpolate(spatial,mcens,ncens)
	spatialvals <- spatialvals / (del1*del2*sum(spatialvals))
	
	#setting up axis-specific torus distance matrices
	
	d1il.mat <- matrix(NA,M.ext,M.ext)
	d2jk.mat <- matrix(NA,N.ext,N.ext)
	for(index in 1:M.ext){
	    abs.diff <- abs(mcens[index]-mcens)
	    tor.diff <- del1*M.ext - abs.diff
	    d1il.mat[index,] <- pmin(abs.diff,tor.diff)
	}
	for(index in 1:N.ext){
	    abs.diff <- abs(ncens[index]-ncens)
	    tor.diff <- del2*N.ext - abs.diff
	    d2jk.mat[index,] <- pmin(abs.diff,tor.diff)
	}
	
	C.tilde <- t(matrix(gu(u=d.func(mat1il=d1il.mat,mat2jk=d2jk.mat,i=1,j=1,l=sort(rep(1:M.ext,N.ext)),k=rep(1:N.ext,M.ext)),sigma=sigma,phi=phi,model=model,additionalparameters=covpars),N.ext,M.ext))
	LAM.tildefft <- Re(fft(C.tilde,inverse=T))
	eigs <-  Re(fft(C.tilde[1,]))
	return(list(LAM.tildefft=LAM.tildefft,cellArea.mat=cellArea.mat,cellInside=cellInside,gridvals=spatialvals,eigs=eigs,Q=fft(C.tilde),mcens=mcens,ncens=ncens))
}



##' fftinterpolate function
##'
##' Generic function used for computing interpolations used in the function \link{fftgrid}.
##'
##' @param spatial an object
##' @param ... additional arguments
##' @return method fftinterpolate
##' @seealso \link{fftgrid}
##' @export

fftinterpolate <- function(spatial,...){
    UseMethod("fftinterpolate")
}



##' interpolate.fromXYZ function
##'
##' This method performs interpolation within the function \code{fftgrid} for \code{fromXYZ} objects.
##'
##' @method fftinterpolate fromXYZ
##' @param spatial objects of class spatialAtRisk
##' @param mcens x-coordinates of interpolation grid in extended space
##' @param ncens y-coordinates of interpolation grid in extended space
##' @param ... additional arguments
##' @return matrix of interpolated values
##' @seealso \link{fftgrid}, \link{spatialAtRisk.fromXYZ}
##' @export

fftinterpolate.fromXYZ <- function(spatial,mcens,ncens,...){
    M.ext <- length(mcens)
    N.ext <- length(ncens)
    M <- M.ext/2
    N <- N.ext/2
	
	spatialvals <- matrix(0,M.ext,N.ext)
    spatialextend <- spatial$Zm
	spatialextend[is.na(spatialextend)] <- 0
	sv <- interp.im(im(t(spatialextend),xcol=spatial$X,yrow=spatial$Y),rep(mcens[1:M],N),rep(ncens[1:N],each=M))
	sv[is.na(sv)] <- 0
	spatialvals[1:M,1:N] <- matrix(sv,M,N)
	return(spatialvals)
}



##' fftinterpolate.fromFunction function
##'
##' This method performs interpolation within the function \code{fftgrid} for \code{fromFunction} objects.
##'
##' @method fftinterpolate fromFunction
##' @param spatial objects of class spatialAtRisk
##' @param mcens x-coordinates of interpolation grid in extended space
##' @param ncens y-coordinates of interpolation grid in extended space
##' @param ... additional arguments
##' @return matrix of interpolated values
##' @seealso \link{fftgrid}, \link{spatialAtRisk.function}
##' @export

fftinterpolate.fromFunction <- function(spatial,mcens,ncens,...){
    M.ext <- length(mcens)
    N.ext <- length(ncens)
    M <- M.ext/2
    N <- N.ext/2
    spatialvals <- matrix(0,M.ext,N.ext)
    xyvals <- matrix(cbind(rep(mcens[1:M],N),rep(ncens[1:N],each=M)),M*N,2)
    interp <- apply(xyvals,1,function(pt){return(spatial$f(pt[1],pt[2]))})
    spatialvals[1:M,1:N] <- matrix(interp,M,N)
    return(spatialvals)
}


##' fftinterpolate.fromSPDF function
##'
##' This method performs interpolation within the function \code{fftgrid} for \code{fromSPDF} objects.
##'
##' @method fftinterpolate fromSPDF
##' @param spatial objects of class spatialAtRisk
##' @param mcens x-coordinates of interpolation grid in extended space
##' @param ncens y-coordinates of interpolation grid in extended space
##' @param ... additional arguments
##' @return matrix of interpolated values
##' @seealso \link{fftgrid}, \link{spatialAtRisk.SpatialPolygonsDataFrame}
##' @export

fftinterpolate.fromSPDF<- function(spatial,mcens,ncens,...){
    M.ext <- length(mcens)
    N.ext <- length(ncens)
    M <- M.ext/2
    N <- N.ext/2
    spatialvals <- matrix(0,M.ext,N.ext)
    xyvals <- SpatialPoints(matrix(cbind(rep(mcens[1:M],N),rep(ncens[1:N],each=M)),M*N,2))
    interp <- overlay(spatial$spdf,xyvals)$atrisk    
    spatialvals[1:M,1:N] <- matrix(interp,M,N)
    spatialvals[is.na(spatialvals)] <- 0
    return(spatialvals)
}
