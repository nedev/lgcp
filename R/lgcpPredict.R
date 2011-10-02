##' lgcpPredict function
##'
##' The function \code{lgcpPredict} performs spatiotemporal prediction for log-Gaussian Cox Processes
##'
##' The following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s,t)}{\mathcal Y(s,t)} be a spatiotemporal Gaussian process, \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space and \eqn{T\subset R_{\geq 0}}{T\subset R_{\geq 0}} be an interval of time of interest. 
##' Cases occur at spatio-temporal positions \eqn{(x,t) \in W \times T}{(x,t) \in W \times T} 
##'  according to an inhomogeneous spatio-temporal Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x,t)}{R(x,t)},
##'   The number of cases, \eqn{X_{S,[t_1,t_2]}}{X_{S,[t_1,t_2]}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} during the interval \eqn{[t_1,t_2]\subseteq T}{[t_1,t_2]\subseteq T} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}}{%
##'    X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\mu(t)\exp\{\mathcal Y(s,t)\}.}{%
##'    R(s,t) = \lambda(s)\mu(t)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1,}{%
##'    \int_W\lambda(s)d s=1,}
##' whilst the fixed temporal component, 
##'  \eqn{\mu:R_{\geq 0}\mapsto R_{\geq 0}}{\mu:R_{\geq 0}\mapsto R_{\geq 0}}, is also a known function with
##' \deqn{\mu(t) \delta t = E[X_{W,\delta t}],}{%
##'    \mu(t) \delta t = E[X_{W,\delta t}],}
##' for \eqn{t}{t} in a small interval of time, \eqn{\delta t}{\delta t}, over which the rate of the process over \eqn{W}{W} can be considered constant.
##'
##' \bold{
##'     NOTE: the xyt stppp object can be recorded in continuous time, but for the purposes of prediciton,    
##'     discretisation must take place. For the time dimension, this is achieved invisibly by \code{as.integer(xyt$t)} and
##'     \code{as.integer(xyt$tlim)}. Therefore, before running an analysis please make sure that this is commensurate
##'     with the physical inerpretation and requirements of your output. The spatial discretisation is
##'     chosen with the argument cellwidth (or gridsize). If the chosen discretisation in time and space is too coarse for a
##'     given set of parameters (sigma, phi and theta) then the proper correlation structures implied by the model will not
##'     be captured in the output.
##' }
##'
##' Before calling this function, the user must decide on the time point of interest, the
##' number of intervals of data to use, the parameters, spatial covariance model, spatial discretisation,
##' fixed spatial (\eqn{\lambda(s)}{\lambda(s)}) and temporal (\eqn{\mu(t)}{\mu(t)}) components, mcmc parameters, and whether or not any output is
##' required.
##'
##' @param xyt a spatio-temporal point pattern object, see ?stppp
##' @param T time point of interest
##' @param laglength specifies lag window, so that data from and including  time (T-laglength) to time T is used in the MALA algorithm
##' @param model.parameters values for parameters, see ?lgcppars
##' @param spatial.covmodel correlation type see ?CovarianceFct 
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @param cellwidth width of grid cells on which to do MALA (grid cells are square). Note EITHER gridsize OR cellwidthe must be specified.
##' @param gridsize size of output grid required. Note EITHER gridsize OR cellwidthe must be specified.
##' @param spatial.intensity the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk
##' @param temporal.intensity the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput
##' @param autorotate logical: whether or not to automatically do MCMC on optimised, rotated grid.   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Set to NULL to estimate this automatically (default). Set to zero for no gradient truncation.
##' further notes on autorotate argument: If set to TRUE, and the argument spatial is not NULL, then the argument spatial must be computed in the original frame of reference (ie NOT in the rotated frame). 
##' Autorotate performs bilinear interpolation (via interp.im) on an inverse transformed grid; if there is no computational advantage in doing this, a warning message will be issued. Note that best accuracy 
##' is achieved by manually rotating xyt and then computing spatial on the transformed xyt and finally feeding these in as arguments to the function lgcpPredict. By default autorotate is set to FALSE.
##' @return the results of fitting the model in an object of class \code{lgcpPredict}
##' @references 
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{KinhomAverage}, \link{ginhomAverage}, \link{lambdaEst}, \link{muEst}, \link{spatialparsEst}, \link{thetaEst},  
##' \link{spatialAtRisk}, \link{temporalAtRisk}, \link{lgcppars}, \link{CovarianceFct}, \link{mcmcpars}, \link{setoutput} 
##' \link{print.lgcpPredict}, \link{xvals.lgcpPredict}, \link{yvals.lgcpPredict}, \link{plot.lgcpPredict}, \link{meanfield.lgcpPredict},
##' \link{rr.lgcpPredict}, \link{serr.lgcpPredict}, \link{intens.lgcpPredict},   
##' \link{varfield.lgcpPredict}, \link{gridfun.lgcpPredict}, \link{gridav.lgcpPredict}, \link{hvals.lgcpPredict}, \link{window.lgcpPredict},
##' \link{mcmctrace.lgcpPredict}, \link{plotExceed.lgcpPredict}, \link{quantile.lgcpPredict}, \link{identify.lgcpPredict}, \link{expectation.lgcpPredict},
##' \link{extract.lgcpPredict}, \link{showGrid.lgcpPredict}, \link{computeGradtrunc}
##' @export 
    
lgcpPredict <- function(xyt,
					    T,
					    laglength,
					    model.parameters=lgcppars(),
					    spatial.covmodel="exponential",
					    covpars=c(),
					    cellwidth=NULL,
					    gridsize=NULL,
					    spatial.intensity,
					    temporal.intensity,					
					    mcmc.control,
					    output.control=setoutput(),
					    autorotate=FALSE,
					    gradtrunc=NULL){
    
    starttime <- Sys.time()
    
    ###
    # Convert times into integer-valued vectors
    ###
    
    if (!inherits(T,"integer")){
	    warning("Converting T into integer value, see ?as.integer",immediate.=TRUE)
	    T <- as.integer(T) 
	}
	if (!inherits(laglength,"integer")){
	    warning("Converting laglength into integer values, see ?as.integer",immediate.=TRUE)
	    laglength <- as.integer(laglength) 
	}
	if (!inherits(xyt$tlim,"integer")){
	    warning("Converting xyt$tlim into integer values, see ?as.integer",immediate.=TRUE)
	    xyt$tlim <- as.integer(xyt$tlim) # convert times into integer values: they should already be in this form.
	}	
	if (!inherits(xyt$t,"integer")){
	    warning("Converting xyt$t into integer values, see ?as.integer",immediate.=TRUE)
	    xyt$t <- as.integer(xyt$t)
	}
	
	###
	# select cellwidth if gridsize specified
	###
	
	if(is.null(cellwidth) & is.null(gridsize)){
	    stop("Either cell width OR grid size must be specified")
	}
	if(!is.null(cellwidth) & !is.null(gridsize)){
	    stop("Either cell width OR grid size must be specified")
	}
	if (!all(sapply(gridsize,is.pow2))){
	    stop("All elements of gridsize must be a power of 2")
	}
	if(!is.null(gridsize) & autorotate==TRUE){
	    warning("In order to use autorotate, you must specify a cell width instead ... SETTING autorotate=FALSE.",immediate.=TRUE)
	    autorotate <- FALSE
	}
	if(!is.null(gridsize)){
	    approxcw <- diff(xyt$window$xrange)/gridsize[1] # approx cell width
	    cwseq <- seq(approxcw/2,2*approxcw,length.out=500)
	    cwfun <- function(cw){
	        ow <- selectObsWindow(xyt,cw)
	        return(c(ow$M-1,ow$N-1))
	    }
	    gsmat <- t(sapply(cwseq,cwfun))
	    tf <- apply(gsmat,1,function(x){return(all(x==gridsize))})
	    if(sum(tf)==0){
	        stop("No sensible observation window found: either change gridsize, or specify cellwidth instead")
	    }
	    else{
	        cellwidth <- cwseq[min(which(tf))]
	    }
	}
    
    ###
    # Perform basic checks 
    ###    					

    if (!is.null(gradtrunc)){
        if(gradtrunc<0){
            stop("gradtrunc must be non-negative")
        }
    }
		
	if (!inherits(temporal.intensity,"temporalAtRisk")){
	    temporal.intensity <- temporalAtRisk(temporal.intensity,tlim=xyt$tlim,xyt=xyt)
	}
	else{
	    if(!all(as.integer(xyt$tlim)==attr(temporal.intensity,"tlim"))){
	        stop("Incompatible temporal.intensity, integer time limits (xyt$tlim and temporal.intensity$tlim) do not match")
	    }
	}
	
	if(laglength==0){
	    stop("laglength must be >= 1")
	}
	
	if(mcmc.control$burnin>mcmc.control$mala.length){
		stop("Number of burnin iterations must be less than the total number of iterations")
	}
	
	aggtimes <- T - laglength:0
	nobser <- 0
	for (i in 1:(laglength+1)){
	    nobser <- nobser + sum(xyt$t==i)
	}
	if(nobser==0){
	    cat("NOTE: time data should be integer-valued.\n")
		stop("No data in chosen time interval")
	}
	
	temporalfit <- sapply(aggtimes,temporal.intensity)
	if (any(is.na(temporalfit))){
	    stop("Missing temporal fitted values")
	}
	
	###
	# Initialise
	###		

	sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu
	theta <- model.parameters$theta
	
	###
    # compute whether there is any efficiency gain in rotating window
    ###  
    
    if (!autorotate){
        test <- roteffgain(xyt,cellwidth)            
    }
    else{
        test <- roteffgain(xyt,cellwidth)
        if (!test){
            warning("There is no gain in efficiency by rotating, see ?roteffgain",immediate.=TRUE)
        }
        rotmat <- getRotation(xyt)$rotation
        xyt <- affine(xyt,mat=rotmat)
    }
    ow <- selectObsWindow(xyt,cellwidth)
	xyt <- ow$xyt
	M <- ow$M
	N <- ow$N
	
	if (M*N>=(256^2)){
	    Sys.sleep(1)
	    cat("\n")
	    warning("USING LARGE FFT GRID: COMPUTATION MAY BE SLOW ON SOME MACHINES ...",.immediate=TRUE)
	    cat("\n")
	}  
	
	cat(paste("FFT Grid size: [",2*(M-1)," , ",2*(N-1),"]\n",sep=""))
	Sys.sleep(1)
    rm(ow)
    
    ###
    # Deal with spatial component and rotate, if necessary
    ###
					
	if(!any(class(spatial.intensity)=="spatialAtRisk")){			
        spatial <- spatialAtRisk(spatial.intensity)
    }
    else{
        spatial <- spatial.intensity
    }
    
    if(autorotate){
        spatial <- affine(spatial,mat=rotmat)
    }

	grid <- fftgrid(xyt=xyt,
					M=M,
					N=N,
					spatial=spatial,
					sigma=sigma,
					phi=phi,
					model=spatial.covmodel,
					covpars=covpars)		
	
	nday <- length(aggtimes) # aggtimes defined above
	tdiffs <- c(NA,diff(aggtimes))	
	vars <- 1 - exp(-2*theta*tdiffs)
	vars[1] <- 1
	
    ###
    # Set up MCMC loop, required to compute nsamp, below
    ###
    mLoop = mcmcLoop(N=mcmc.control$mala.length,burnin=mcmc.control$burnin,thin=mcmc.control$retain,progressor=mcmcProgressTextBar)	
	
	# issue warning if dumping information to disc
	nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
	if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1]=="dump2dir"){
    	cat("WARNING: disk space required for saving is approximately ",round(nsamp*object.size(array(runif((M-1)*(N-1)*nday),dim=c((M-1),(N-1),nday)))/1024^2,2)," Mb, ",sep="")
        if (!output.control$gridfunction$forceSave){
            m <- menu(c("yes","no"),title="continue?")
            if(m==1){
                cat("Note: to bypass this menu, set forceSave=TRUE in dump2dir\n")
                Sys.sleep(2)
            }
            else{
                stop("Stopped")
            }
        }
    }
	
	nis <- list()
	for(i in 1:nday){
	    if (sum(xyt$t==aggtimes[i])>0){
		    nis[[i]] <- getCounts(xyt=xyt,subset=(xyt$t==aggtimes[i]),M=M,N=N)
		}
		else{
		    nis[[i]] <- matrix(0,2*M-2,2*N-2)
		}
		ct1 <- sum(nis[[i]])
		nis[[i]] <- nis[[i]] * (grid$gridvals>0)
		ct2 <- sum(nis[[i]])
		if(ct2<ct1){
		    warning(paste("Time ",aggtimes[i],": ",ct1-ct2," data points lost due to discretisation.",sep=""),immediate.=TRUE)
		}
	}
	
	###
	# Compute gradient truncation, if necessary
	###
	
	if(is.null(gradtrunc)){
        gradtrunc <- computeGradtrunc(  nis=nis,  
                                        theta=theta,
    	                                timediffs=tdiffs,
    		                            mu=mu,
    		                            temporal=temporalfit,
    		                            grid=grid)
    }
	
	###
	# Run MALA
	##
	
	gridfun <- output.control$gridfunction
	if (is.null(gridfun)){
	    gridfun <- nullFunction()
	}
    gridav <- output.control$gridmeans
	if (is.null(gridav)){
	    gridav <- nullAverage()
	}   
    
    lg <- MALAlgcp( mLoop,
                    inits=mcmc.control$inits,
                    adaptivescheme=mcmc.control$adaptivescheme,
                    THETA=theta,
                    vars=vars,
                    nis=nis,
                    tdiffs=tdiffs,
                    cellArea.mat=grid$cellArea.mat,
                    MU=mu,
                    spatialvals=grid$gridvals,
                    temporal.fitted=temporalfit,
                    LAM.tildefft=grid$LAM.tildefft,
                    cellInside=grid$cellInside,
                    MCMCdiag=mcmc.control$MCMCdiag,
                    gradtrunc=gradtrunc,
                    gridfun=gridfun,
                    gridav=gridav)
	
	del1 <- (xyt$window$xrange[2]-xyt$window$xrange[1])/(M-1)
	del2 <- (xyt$window$yrange[2]-xyt$window$yrange[1])/(N-1) 
	mcens <- xyt$window$xrange[1]+.5*del1+(0:(M-2))*del1
	ncens <- xyt$window$yrange[1]+.5*del2+(0:(N-2))*del2
	
	endtime <- Sys.time()
	timetaken <- endtime-starttime
	
	lg$xyt <- xyt
	if(autorotate){
	    lg$rotation <- rotmat # required for post processing
	}
	lg$M <- M
	lg$N <- N
	lg$aggtimes <- aggtimes
	lg$tdiffs <- tdiffs
	lg$vars <- vars
	lg$spatial <- spatial
	lg$temporal <- temporalfit
	lg$grid <- grid
	lg$nis <- lgcpgrid(nis)
	lg$mcens <- mcens
	lg$ncens <- ncens
	lg$sigma <- sigma
	lg$phi <- phi
	lg$mu <- mu
	lg$theta <- theta
	lg$mcmcpars <- mcmc.control
	lg$timetaken <- timetaken
	
	class(lg) <- c("lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}					
