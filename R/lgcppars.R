##' lgcppars function
##'
##' A function for setting the parameters sigma, phi and theta for \code{lgcpPredict}. Note that the returned
##' set of parameters also features mu=-0.5*sigma^2, gives mean(exp(Y)) = 1.
##'
##' @param sigma sigma parameter
##' @param phi phi parameter
##' @param theta this is beta parameter in Brix and Diggle (2001)
##' @seealso \link{lgcpPredict}
##' @export

lgcppars <- function(sigma=NULL,phi=NULL,theta=NULL){
    
	return(list(sigma=sigma,phi=phi,mu=-0.5*sigma^2,theta=theta))	
	# NB choice of parameter mu=-0.5*sigma^2 gives mean(exp(Y)) = 1  					
}						
