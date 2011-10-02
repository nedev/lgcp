###
# FFT grid handling functions
###

##' fftgrid function
##'
##' \bold{Advanced use only.} Computes various quantities for use in \code{lgcpPredict},
##' \code{lgcpSim} and \link{computeGradtrunc}.
##' 
##' @param xyt object of class stppp
##' @param M number of 'gridlines' in x-direction (= number of centroids + 1)
##' @param N number of 'gridlines' in y-direction
##' @param spatial an object of class spatialAtRisk
##' @param sigma scaling paramter for spatial covariance function, see Brix and Diggle (2001) 
##' @param phi scaling paramter for spatial covariance function, see Brix and Diggle (2001)
##' @param model correlation type see ?CovarianceFct
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @return fft objects for use in MALA
##' @seealso \link{computeGradtrunc}
##' @examples
##' xyt <- stppp(ppp(),t=c(),tlim=c(0,1))
##' sar <- spatialAtRisk(function(x,y){return(1)},warn=FALSE)
##' grid <- fftgrid(xyt,65,65,sar,1,0.01,"exponential",covpars=c())
##' @export


fftgrid <- function(xyt,M,N,spatial,sigma,phi,model,covpars){

    verifyclass(xyt,"stppp")
    verifyclass(spatial,"spatialAtRisk")
    
    study.region <- xyt$window
	
	## DEFINE LATTICE & CENTROIDS ##
		
	m <- 0:(M-1)
	n <- 0:(N-1)
	del1 <- (study.region$xrange[2]-study.region$xrange[1])/(M-1)
	del2 <- (study.region$yrange[2]-study.region$yrange[1])/(N-1) 
	
	M.ext <- 2*M-1 ## since we are working with centroids, it is M.ext-1 and N.ext-1 that must be powers of 2, NOT M.ext and N.ext alone as in m & w 2004	
	N.ext <- 2*N-1 ##
	
	mcens <- study.region$xrange[1]+.5*del1+(0:(M.ext-2))*del1
	ncens <- study.region$yrange[1]+.5*del2+(0:(N.ext-2))*del2	
	
	## REQUIRED SIMULATION QUANTITIES ##
	
	cellArea.mat <- matrix(0,M.ext-1,N.ext-1)
	cellArea.mat[1:(M-1),1:(N-1)] <- del1*del2
	
	cellInside <- inside.owin(x=sort(rep(mcens,N.ext-1)),y=rep(ncens,M.ext-1),w=study.region)
	cellInside <- matrix(as.logical(cellInside),M.ext-1,N.ext-1,byrow=T)[1:(M-1),1:(N-1)]
	cellInsideBIG <- matrix(0,M.ext-1,N.ext-1)
	cellInsideBIG[1:(M-1),1:(N-1)][cellInside] <- 1
	
	## OBTAIN SPATIAL VALS ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatialvals <- fftinterpolate(spatial,mcens,ncens)
	spatialvals <- spatialvals / (del1*del2*sum(spatialvals))
	
	#setting up axis-specific torus distance matrices
	
	d1il.mat <- matrix(NA,M.ext-1,M.ext-1)
	d2jk.mat <- matrix(NA,N.ext-1,N.ext-1)
	for(index in 1:(M.ext-1)){
	    abs.diff <- abs(mcens[index]-mcens)
	    tor.diff <- del1*(M.ext-1) - abs.diff
	    d1il.mat[index,] <- pmin(abs.diff,tor.diff)
	}
	for(index in 1:(N.ext-1)){
	    abs.diff <- abs(ncens[index]-ncens)
	    tor.diff <- del2*(N.ext-1) - abs.diff
	    d2jk.mat[index,] <- pmin(abs.diff,tor.diff)
	}
	
	C.tilde <- t(matrix(gu(u=d.func(mat1il=d1il.mat,mat2jk=d2jk.mat,i=1,j=1,l=sort(rep(1:(M.ext-1),(N.ext-1))),k=rep(1:(N.ext-1),(M.ext-1))),sigma=sigma,phi=phi,model=model,additionalparameters=covpars),N.ext-1,M.ext-1))
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
    M.ext <- length(mcens) + 1
    N.ext <- length(ncens) + 1
    M <- (M.ext+1)/2
    N <- (N.ext+1)/2
	
	spatialvals <- matrix(0,M.ext-1,N.ext-1)
    spatialextend <- spatial$Zm
	spatialextend[is.na(spatialextend)] <- 0
	sv <- interp.im(im(t(spatialextend),xcol=spatial$X,yrow=spatial$Y),rep(mcens[1:(M-1)],N-1),rep(ncens[1:(N-1)],each=M-1))
	sv[is.na(sv)] <- 0
	spatialvals[1:(M-1),1:(N-1)] <- matrix(sv,M-1,N-1)
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
    M.ext <- length(mcens) + 1
    N.ext <- length(ncens) + 1
    M <- (M.ext+1)/2
    N <- (N.ext+1)/2
    spatialvals <- matrix(0,M.ext-1,N.ext-1)
    xyvals <- matrix(cbind(rep(mcens[1:(M-1)],N-1),rep(ncens[1:(N-1)],each=M-1)),(M-1)*(N-1),2)
    interp <- apply(xyvals,1,function(pt){return(spatial$f(pt[1],pt[2]))})
    spatialvals[1:(M-1),1:(N-1)] <- matrix(interp,M-1,N-1)
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
    M.ext <- length(mcens) + 1
    N.ext <- length(ncens) + 1
    M <- (M.ext+1)/2
    N <- (N.ext+1)/2
    spatialvals <- matrix(0,M.ext-1,N.ext-1)
    xyvals <- SpatialPoints(matrix(cbind(rep(mcens[1:(M-1)],N-1),rep(ncens[1:(N-1)],each=M-1)),(M-1)*(N-1),2))
    interp <- overlay(spatial$spdf,xyvals)$atrisk    
    spatialvals[1:(M-1),1:(N-1)] <- matrix(interp,M-1,N-1)
    spatialvals[is.na(spatialvals)] <- 0
    return(spatialvals)
}
