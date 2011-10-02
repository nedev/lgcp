##' lgcpForecast function
##'
##' Function to produce forecasts for the mean field \eqn{Y}{Y} at times beyond the last time point in the 
##' analysis (given by the argument \code{T} in the function \code{lgcpPredict}).
##'
##' @param lg an object of class lgcpPredict
##' @param ptimes vector of time points for prediction. Must start strictly after last inferred time point.
##' @return forcasted Y values over grid
##' @references Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##' @seealso \link{lgcpPredict}
##' @export 

lgcpForecast <- function(lg,ptimes){

    verifyclass(lg,"lgcpPredict")
    
    maxt <- max(lg$aggtimes)
    
    if(any(ptimes <=maxt)){
        stop("required prediction times must start strictly after last inferred time point.")
    }    
    
    tdiffs <- ptimes - maxt    
    n <- length(tdiffs)      		
    
    ymats <- list()  		
    for(j in 1:n){ 
        ymats[[j]] <- at(tdiffs[j],lg$mu,lg$theta)+bt.scalar(tdiffs[j],lg$theta)*lg$y.mean[[length(lg$y.mean)]]
    }
    
    return(ymats)
    
}
