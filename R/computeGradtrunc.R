##' computeGradtrunc function
##'
##' \bold{Advanced use only.} A function to compute a gradient truncation parameter for MALA via simulation. The function
##' requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.
##'
##' @param nsims The number of simulations to use in computation of gradient truncation.
##' @param scale multiplicative scaling constant, returned value is scale*max(gradient over simulations). Default scale is 1.
##' @param nis cell counts on the extended grid
##' @param theta temporal correlation parameter 
##' @param timediffs vector of time-differences between time intervals under consideration
##' @param mu mean parameter, see Brix and Diggle
##' @param temporal temporal fitted values
##' @param grid output from function fftgrid
##' @return gradient truncation parameter
##' @seealso \link{fftgrid}
##' @export 

computeGradtrunc <- function(	nsims=100,
                                scale=1,
                                nis,
                                theta,
    		                    timediffs,
    		                    mu,
    		                    temporal,
    		                    grid){
    
    M.ext <- dim(grid$LAM.tildefft)[1] + 1
    N.ext <- dim(grid$LAM.tildefft)[2] + 1
    
    cat("Computing gradient truncation ...\n")
    grds <- c()
    pb <- txtProgressBar(min=1,max=nsims,style=3)
    for (i in 1:nsims){
        at.mat.list <- list()
        for(j in 1:length(timediffs)){
            at.mat.list[[j]] <- matrix(rnorm((M.ext-1)*(N.ext-1)),M.ext-1,N.ext-1)
        }  		                    
        		                    
        grds[i] <- suppressWarnings(target.and.grad(at.mat.list=at.mat.list,
                                                    theta=theta,
                                                    varvals=1-exp(-2*theta*timediffs),
                                                    nicounts.mat.list=nis,
                                                    timediffs=timediffs,
                                                    car.mat=grid$cellArea.mat,
                                                    mu=mu,
                                                    spatial.mat=grid$gridvals,
                                                    temporal=temporal,
                                                    lam.tilde=grid$LAM.tildefft,
                                                    cin=grid$cellInside,
                                                    M.ext=M.ext,    
                                                    N.ext=N.ext,
                                                    gradtrunc=NA,
                                                    estgradtrunc=TRUE))
        setTxtProgressBar(pb,i)                                                        		                    
    }
    close(pb)
    gt <- floor(scale*max(grds,na.rm=TRUE))
    if (is.na(gt) | is.infinite(gt) | is.nan(gt)){
        stop("Could not compute gradient truncation. To set manually, see ?lgcpPredict") 
    }
    cat(paste("Using gradient truncation of ",gt,"\n",sep="")) 
    return(gt)
}