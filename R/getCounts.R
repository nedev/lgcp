##' getCounts function
##'
##' This function is used to count the number of observations falling inside grid cells, the output
##' is used in the function \link{lgcpPredict}.
##'
##' @param xyt stppp or ppp data object
##' @param subset Logical vector. Subset of data of interest, by default this is all data.
##' @param M number of 'gridlines' in x-direction (= number of centroids + 1)
##' @param N number of 'gridlines' in y-direction
##' @return The number of observations in each grid cell returned on a grid suitable for use in the extended FFT space.
##' @seealso \link{lgcpPredict}
##' @examples
##' xyt <- stppp(ppp(runif(100),runif(100)),t=1:100,tlim=c(1,100))
##' cts <- getCounts(xyt,M=65,N=65) # gives an output grid of size 128 by 128
##' ctssub <- cts[1:64,1:64] # returns the cell counts in the observation
##'                          # window of interest
##' @export
getCounts <- function(xyt,subset=rep(TRUE,xyt$n),M,N){
	if(M<3 | N<3){stop("M and/or N too small")}
	
	test1 <- FALSE
	test2 <- FALSE
	for (i in 1:20){
		if (2^i==M-1){test1 <- TRUE}
		if (2^i==N-1){test2 <- TRUE}
	}
	if ((!test1)&(!test2)){
		stop("Both (M-1) and (N-1) must be a power of 2 and both M,N <= 2^20+1")
	}
	
	xran <- xyt$window$xrange
	yran <- xyt$window$yrange
	bbox <- rbind(xran,yran)
	xwd <- diff(bbox[1,])/(M-1)
	ywd	<- diff(bbox[2,])/(N-1)
	sg <- SpatialGrid(GridTopology(c(bbox[1,1]+xwd/2,bbox[2,1]+ywd/2),c(xwd,ywd),c(M-1,N-1)))
	sp <- SpatialPoints(matrix(cbind(xyt$x,xyt$y)[subset,],sum(subset),2))
	ol <- overlay(sg,sp)
	tol <- table(ol)
	idx <- as.numeric(rownames(tol))	
	smat <- matrix(0,M-1,N-1)
	smat[idx] <- tol
	smat <- smat[,(N-1):1] # change orientation 
	
	nis <- matrix(0,(2*M-2),(2*N-2))
	nis[1:(M-1),1:(N-1)] <- smat 
 
    return(nis)
}
