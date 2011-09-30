##' lgcpsim function
##'
##' Approximate simulation from a spatiotemoporal log-Gaussian Cox Process. Returns an stppp object.
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
##' @param owin a spatstat owin object, possibly with polygonal boundary
##' @param tlim time interval on which to simulate data
##' @param spatial.intensity object that can be coerced into a spatialAtRisk object. if NULL then uniform spatial is chosen 
##' @param temporal.intensity the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk
##' @param cellwidth size of gris discretisation as in lgcpPredict
##' @param model.parameters parameters of hte Gaussian process, specify using function lgcppars,
##' @param spatial.covmodel spatial covariance model, default is "exponential"
##' @param covpars additional covariance parameters
##' @param progressbar logical, whether to print a progress bar. Default TRUE.
##' @param returnintensities logigal, whether to return the vector of simulation times and corresponding array of spatial intensities at each time. Default FALSE.
##' @param plot plot each time step? Only really useful for debugging purposes.
##' @return returns an object of class stppp containing the simulated data
##' @references 
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{lgcpPredict}, \link{showGrid.stppp}, \link{stppp}
##' @examples xyt <- lgcpSim()
##' @export



lgcpSim <- function(owin=NULL,
                    tlim=as.integer(c(0,10)),
                    spatial.intensity=NULL,
                    temporal.intensity=NULL,
                    cellwidth = 0.05,
                    model.parameters=lgcppars(sigma=2,phi=0.2,theta=3),
                    spatial.covmodel="exponential",
                    covpars=c(),
                    progressbar=TRUE,
                    returnintensities=FALSE,
                    plot=FALSE){
                    
    if (!inherits(tlim,"integer")){
	    warning("Converting tlim into integer values, see ?as.integer")
	    tlim <- as.integer(tlim) # convert times into integer values: they should already be in this form.
	}
	tlim <- sort(tlim)
	if (tlim[1]==tlim[2]){
	    stop("Length of time interval given by as.integer(tlim) must be >= 1")
	}
	toffset <- tlim[1]
	maxt <- tlim[2] - toffset
	                        
                              
    sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu
	theta <- model.parameters$theta

    if(is.null(owin)){
        owin <- owin()
    }    
    
    const0 <- 0.05 # level below which correlation function must drop to be considered not important
    const1 <- 10 # time discretisation
    const2 <- 100
    if (is.null(temporal.intensity)){
        temporal.intensity <- constantInTime(const2,tlim)
    }
    else{
        if (!inherits(temporal.intensity,"temporalAtRisk")){
            temporal.intensity <- temporalAtRisk(temporal.intensity,tlim)
        }
        if(!all(tlim==attr(temporal.intensity,"tlim"))){
	        stop("Incompatible temporal.intensity, integer time limits (tlim and temporal.intensity$tlim) do not match")
	    }
    }
    
    
    sleeplen <- 0.1
    
    # check / choose time discretisation
    c1 <- -(1/theta)*log(const0) # the number at which exp(-theta*x) = 0.05, use this to choose time discretisation 
    ndivs <- const1 * ceiling(maxt/c1)
    if (progressbar){
        pb <- txtProgressBar(min=1,max=ndivs,style=3)
    }
    tdiff = maxt/ndivs
    times <- tdiff/2 + tdiff*(0:(ndivs-1))    
    
    # check space discretisation
    c2 <- -phi*log(const0) # the number at which exp(-theta*x) = 0.05, use this to choose time discretisation 
    if (cellwidth>c2/2){
        warning(paste("cellwidth should be at least",c2/2,"to get accurate results."))
    }
     
    xyt <- ppp(window=owin)
    xyt <- stppp(xyt,t=c(),tlim=tlim) 
      
    ow <- selectObsWindow(xyt,cellwidth)
	xyt <- ow$xyt
	M <- ow$M
	N <- ow$N
	if(returnintensities){
	    intensities <- array(NA,c(M-1,N-1,ndivs))
	}
	cat(paste("Grid size: [",M-1," , ",N-1,"]\n",sep=""))
	if(is.null(spatial.intensity)){
        spatial <- spatialAtRisk(list(X=seq(1/(2*M),1,by=1/M),Y=seq(1/(2*N),1,by=1/N),Zm=matrix(1/(M*N),M,N)))
    }
    else{
        if(!any(class(spatial.intensity)=="spatialAtRisk")){		
            spatial <- spatialAtRisk(spatial.intensity)
        }
        else{
            spatial <- spatial.intensity
        }
    } 
    grid <- fftgrid(xyt=xyt,
					M=M,
					N=N,
					spatial=spatial,
					sigma=sigma,
					phi=phi,
					model=spatial.covmodel,
					covpars=covpars)
    
    spatialvals <- grid$gridvals[1:(M-1),1:(N-1)]
    LAM.tildefft <- grid$LAM.tildefft
    M.ext <- dim(LAM.tildefft)[1] + 1 
	N.ext <- dim(LAM.tildefft)[2] + 1
	
	xdiv <- diff(xyt$window$xrange) / (M-1)
	ydiv <- diff(xyt$window$yrange) / (N-1)
	gridx <- xyt$window$xrange[1] + xdiv/2 + xdiv*(0:(M-2))
	gridy <- xyt$window$yrange[1] + ydiv/2 + ydiv*(0:(N-2))
	xvals <- rep(gridx,N-1)
    yvals <- rep(gridy,each=M-1) # note use of 'rev' function here      
	
	###
	# Simulate first Y and observations  
	###
	obs <- matrix(rep(NA,3),1,3) # (x,y,t)
	GAM <- matrix(rnorm((M.ext-1)*(N.ext-1)),M.ext-1,N.ext-1)
	Y <- (Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(sqrt(LAM.tildefft)*((1/sqrt((M.ext-1)*(N.ext-1)))*fft(GAM)),inverse=T)) + mu)[1:(M-1),1:(N-1)]
    mut <- temporal.intensity(times[1]+toffset)
	rate <- tdiff * mut * spatialvals * cellwidth^2 * exp(Y)
	if(returnintensities){
	    intensities[,,1] <- rate
	}
    counts <- rpois((M-1)*(N-1),lambda=as.vector(rate))
    idx <- which(counts!=0)
    if(plot){
        rate[rate==0] <- NA
        image(gridx,gridy,rate^0.25)
    }
    for (j in idx){
        for (k in 1:counts[j]){
            obs <- rbind(obs,c(xvals[j]+runif(1,-xdiv/2,xdiv/2),yvals[j]+runif(1,-ydiv/2,ydiv/2),runif(1,0,times[1]+tdiff/2))) # note noise added here, hence 'approximate simulation'
            if(plot){
                points(xvals[j]+runif(1,-xdiv/2,xdiv/2),yvals[j]+runif(1,-ydiv/2,ydiv/2),pch="+",cex=0.5) # debugging only
            }
        }
    }
    if(plot){
        Sys.sleep(sleeplen)
    }
	
	###
	# Simulate remaining Ys and observations  
	###
	for(i in 2:ndivs){
	    mut <- temporal.intensity(times[i]+toffset)
	    GAM <- matrix(rnorm((M.ext-1)*(N.ext-1)),M.ext-1,N.ext-1)
	    Y <- at(tdiff,mu,theta)+bt.scalar(tdiff,theta)*Y + sqrt((1-exp(-2*theta*tdiff)))*Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(sqrt(LAM.tildefft)*((1/sqrt((M.ext-1)*(N.ext-1)))*fft(GAM)),inverse=T))[1:(M-1),1:(N-1)]
	    rate <- tdiff * mut * spatialvals * cellwidth^2 * exp(Y)
	    if(returnintensities){
    	    intensities[,,i] <- rate
    	}
        counts <- rpois((M-1)*(N-1),lambda=as.vector(rate))
        idx <- which(counts!=0)
        if(plot){
            rate[rate==0] <- NA
            image(gridx,gridy,rate^0.25)
        }
        for (j in idx){
            for (k in 1:counts[j]){
                obs <- rbind(obs,c(xvals[j]+runif(1,-xdiv/2,xdiv/2),yvals[j]+runif(1,-ydiv/2,ydiv/2),runif(1,times[i]-tdiff/2,times[i]+tdiff/2))) # note noise added here, hence 'approximate simulation'
                if(plot){
                    points(xvals[j]+runif(1,-xdiv/2,xdiv/2),yvals[j]+runif(1,-ydiv/2,ydiv/2),pch="+",cex=0.5) # debugging only
                }
            }
        }
        if(plot){
            Sys.sleep(sleeplen)
        }
        if (progressbar){
            setTxtProgressBar(pb,i)
        }
	}
	
	if (progressbar){
	    close(pb)
	}
	n <- dim(obs)[1]-1
	if(n==0){
	    stop("No data generated for chosen parameters")
	}
	obs <- matrix(obs[-1,],n,3)  
    
    idx <- inside.owin(x=obs[,1],y=obs[,2],w=xyt$window)
    if(sum(idx)==0){
	    stop("No data generated for chosen parameters")
	}
	obs <- obs[idx,]	
	
    xyt$x <- obs[,1]
    xyt$y <- obs[,2]
    xyt$t <- toffset + obs[,3]
    xyt$n <- length(xyt$x)
    attr(xyt,"xvals") <- gridx
    attr(xyt,"yvals") <- gridy
    attr(xyt,"grid") <- grid
    if(returnintensities){
        attr(xyt,"times") <- times + toffset
        attr(xyt,"intensities") <- intensities
    }
    return(xyt)
}
