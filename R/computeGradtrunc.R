##' computeGradtrunc function
##'
##' NOTE THIS FUNCTiON IS NOW REDUNDANT AND HAS BEEN REPLACED BY computeGradtruncSpatioTemporal
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
        stop("Could not compute gradient truncation. To set manually, see gradtrunc argument of  ?lgcpPredict") 
    }
    cat(paste("Using gradient truncation of ",gt,"\n",sep="")) 
    return(gt)
}



##' computeGradtruncSpatial function
##'
##' \bold{Advanced use only.} A function to compute a gradient truncation parameter for 'spatial only' MALA via simulation. The function
##' requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.
##'
##' @param nsims The number of simulations to use in computation of gradient truncation.
##' @param scale multiplicative scaling constant, returned value is scale*max(gradient over simulations). Default scale is 1.
##' @param nis cell counts on the extended grid
##' @param mu parameter of latent field, mu 
##' @param rootQeigs root of eigenvalues of precision matrix of latent field 
##' @param invrootQeigs reciprocal root of eigenvalues of precision matrix of latent field
##' @param scaleconst expected number of cases, or ML estimate of this quantity
##' @param spatial spatial at risk interpolated onto grid of requisite size
##' @param cellarea cell area
##' @return gradient truncation parameter
##' @seealso \link{fftgrid}
##' @export 

computeGradtruncSpatial <- function(nsims=100,
                                    scale=1,
                                    nis,
                                    mu,
                                    rootQeigs,
                                    invrootQeigs,
                                    scaleconst,
                                    spatial,
                                    cellarea){
    
    cat("Computing gradient truncation ...\n")
    grds <- c()
    pb <- txtProgressBar(min=1,max=nsims,style=3)
    M <- nrow(nis)
    N <- ncol(nis)
    for (i in 1:nsims){
        Gamma <- matrix(rnorm(M*N),M,N)
        Y <- YfromGamma(Gamma=Gamma,invrootQeigs=invrootQeigs,mu=mu)
        expY <- exp(Y)                  		                    
        grds[i] <- max(expY) ####suppressWarnings(max((-1)*Gamma +(1/length(Y))*Re(fft(fft(nis-scaleconst*spatial*expY*cellarea,inverse=TRUE)*rootQeigs,inverse=TRUE))))
        setTxtProgressBar(pb,i)                                                        		                    
    }
    close(pb)
    gt <- floor(scale*max(grds,na.rm=TRUE))
    if (is.na(gt) | is.infinite(gt) | is.nan(gt)){
        stop("Could not compute gradient truncation. To set manually, see gradtrunc argument of ?lgcpPredictSpatial") 
    }
    cat(paste("Using gradient truncation of ",gt,"\n",sep="")) 
    return(gt)
}


##' computeGradtruncSpatioTemporal function
##'
##' \bold{Advanced use only.} A function to compute a gradient truncation parameter for 'spatial only' MALA via simulation. The function
##' requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.
##'
##' @param nsims The number of simulations to use in computation of gradient truncation.
##' @param scale multiplicative scaling constant, returned value is scale*max(gradient over simulations). Default scale is 1.
##' @param nis cell counts on the extended grid
##' @param mu parameter of latent field, mu 
##' @param rootQeigs root of eigenvalues of precision matrix of latent field 
##' @param invrootQeigs reciprocal root of eigenvalues of precision matrix of latent field
##' @param spatial spatial at risk interpolated onto grid of requisite size
##' @param temporal fitted temporal values
##' @param bt vectoer of variances b(delta t) in Brix and Diggls 2001
##' @param cellarea cell area
##' @return gradient truncation parameter
##' @seealso \link{fftgrid}
##' @export 

computeGradtruncSpatioTemporal <- function( nsims=100,
                                            scale=1,
                                            nis,
                                            mu,
                                            rootQeigs,
                                            invrootQeigs,
                                            spatial,
                                            temporal,
                                            bt,
                                            cellarea){
    
    cat("Computing gradient truncation ...\n")
    grds <- c()
    pb <- txtProgressBar(min=1,max=nsims,style=3)
    M <- nrow(nis[[1]])
    N <- ncol(nis[[1]])
    for (i in 1:nsims){
        Gamma <- list()                           
        
        lapply(1:length(nis),function(i){Gamma[[i]]<<-matrix(rnorm(M*N),M,N)})
        Y <- lapply(Gamma,YfromGamma,invrootQeigs=invrootQeigs,mu=mu)
        expY <- lapply(Y,exp)                 		                    
        
        grds[i] <- suppressWarnings(max(sapply(expY,max)))        
        
        setTxtProgressBar(pb,i)                                                        		                    
    }
    close(pb)
    gt <- floor(scale*max(grds,na.rm=TRUE))
    if (is.na(gt) | is.infinite(gt) | is.nan(gt)){
        stop("Could not compute gradient truncation. To set manually, see gradtrunc argument of ?lgcpPredictSpatial") 
    }
    cat(paste("Using gradient truncation of ",gt,"\n",sep="")) 
    return(gt)
}

