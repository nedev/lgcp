##' lgcpForecast function
##'
##' Function to produce forecasts for the mean field \eqn{Y}{Y} at times beyond the last time point in the 
##' analysis (given by the argument \code{T} in the function \code{lgcpPredict}).
##'
##' @param lg an object of class lgcpPredict
##' @param ptimes vector of time points for prediction. Must start strictly after last inferred time point.
##' @param spatial.intensity the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk
##' @param temporal.intensity the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk
##' @return forcasted relative risk, Poisson intensities and Y values over grid, together with approximate variance.
##' @references Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##' @seealso \link{lgcpPredict}
##' @export 

lgcpForecast <- function(lg,ptimes,spatial.intensity,temporal.intensity){

    verifyclass(lg,"lgcpPredict")
    
    maxt <- max(lg$aggtimes)
   
    if(any(ptimes <=maxt)){
        stop("required prediction times must start strictly after last inferred time point.")
    }    
    
    if(!inherits(temporal.intensity,"temporalAtRisk")){
        temporal.intensity <- temporalAtRisk(temporal.intensity,tlim=range(ptimes))
    }
    temporalfit <- sapply(ptimes,temporal.intensity)

    MN <- dim(lg$EY.mean$grid[[1]])
    cellarea <- diff(lg$mcens[1:2])*diff(lg$ncens[1:2]) 

    tdiffs <- ptimes - maxt    
    n <- length(tdiffs) 
    gt <- 1 - exp(-2*tdiffs*lg$theta)

    spatial <- spatial.intensity
    if(!inherits(spatial.intensity,"spatialAtRisk")){
        spatial <- spatialAtRisk(spatial.intensity)
    }    
    spatialvals <- fftinterpolate(spatial,c(lg$mcens,lg$mcens[length(lg$mcens)]+diff(lg$mcens[1:2])+lg$mcens-lg$mcens[1]),c(lg$ncens,lg$ncens[length(lg$ncens)]+diff(lg$ncens[1:2])+lg$ncens-lg$ncens[1]),ext=lg$ext)[1:lg$M,1:lg$N]   
    cellInside <- matrix(as.numeric(inside.owin(x=sort(rep(lg$mcens,lg$N)),y=rep(lg$ncens,lg$M),w=lg$xyt$window)),lg$M,lg$N,byrow=TRUE)
    spatialvals <- spatialvals*cellInside
    spatialvals <- spatialvals / (cellarea*sum(spatialvals)) 
    
    ymats <- list()
    yvars <- list()
    rr <- list()
    rr.var <- list()
    intens <- list()
    intens.var <- list()  		
    for(j in 1:n){ 
        ymats[[j]] <- at(tdiffs[j],lg$mu,lg$theta)+bt.scalar(tdiffs[j],lg$theta)*lg$y.mean$grid[[length(lg$y.mean$grid)]]
        yvars[[j]] <- (bt.scalar(tdiffs[j],lg$theta)^2)*lg$y.var$grid[[length(lg$y.var$grid)]] + gt[j]*lg$sigma^2 
        rr[[j]] <- exp(ymats[[j]])
        rr.var[[j]] <- exp(ymats[[j]]-yvars[[j]]/2) 
        intens[[j]] <- cellarea*temporalfit[j]*spatialvals*rr[[j]]
        intens.var[[j]] <- (cellarea*temporalfit[j]*spatialvals)^2*rr.var[[j]]
    }
    
    return(list(relrisk=lgcpgrid(rr),relrisk.var=lgcpgrid(rr.var),intensity=lgcpgrid(intens),intensity.var=lgcpgrid(intens.var),y.mean=lgcpgrid(ymats),y.var=lgcpgrid(yvars),mcens=lg$mcens,ncens=lg$ncens))
    
}
