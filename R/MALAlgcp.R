##' MALAlgcp function
##'
##' \bold{For Internal/Advanced Use Only.}
##'
##' @param mcmcloop an MCMC loop controller (see ?mcmcLoop)
##' @param inits optional initial point for MALA - useful if updating on a daily basis
##' @param adaptivescheme the type of adaptive mcmc to use, see ?constanth (constant h) or ?andrieuthomsh (adaptive MCMC of Andrieu and Thoms (2008))
##' @param THETA corresponds to parameter beta in Brix & Diggle (2001) 	
##' @param vars pre-computed coefficient (1-exp(-2*beta*delta_t)) in G(t)
##' @param nis these are the cell counts at each time lag
##' @param tdiffs time differences (note first entry is an NA)
##' @param cellArea.mat matrix of cell areas
##' @param MU estimate of parameter mu in Brix & Diggle (2001) 
##' @param spatialvals fixed spatial component
##' @param temporal.fitted fixed temporal component
##' @param LAM.tildefft transpose of Fourier transform on extended grid, output from function fftgrid
##' @param cellInside (M-1)x(N-1) matrix of TRUE/FALSE indicating whether cell is inside observation window
##' @param MCMCdiag non-negative integer, if greater than zero saves information from the MCMC chain
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Set to zero for no gradient truncation (default).
##' @param gridfun chosen grid function, see ?nullfunction (default) or ?dump2dir for example
##' @param gridav chosen grid average, see ?nullAverage (default) or ?MonteCarloAverage for example
##' @return mean Gaussian field from MCMC, exceedances, output of MCMC chain  etc (if specified)
##' @export 

MALAlgcp <-
function(   mcmcloop,
            inits,
            adaptivescheme,
            THETA,
            vars,
            nis,
            tdiffs,
            cellArea.mat,
            MU,
            spatialvals,
            temporal.fitted,
            LAM.tildefft,
            cellInside,
            MCMCdiag,
            gradtrunc,
            gridfun,
            gridav){
###

	time.block <- length(tdiffs) # note this assumes the first element of tdiffs is an NA
	M.ext <- dim(LAM.tildefft)[1] + 1 # this will need to be checked for general M and N 
	N.ext <- dim(LAM.tildefft)[2] + 1 # this will need to be checked for general M and N
	M <- (M.ext + 1)/2
	N <- (N.ext + 1)/2
	
	GFinitialise(gridfun) # note these two lines mst come after M and N have been computed
	GAinitialise(gridav) # note these two lines mst come after M and N have been computed
    
    h <- initialiseAMCMC(adaptivescheme)
    hrec <- h	
	
	GAM <- grad.means <- gamma.norms <- y.means <- ymats <- y.mean <- y.var <- EY.mean <- EY.var <- exy.means <- gam.sum <- y.sum <- jdmats <- Ysample <- GAMsample <- list()
	acceptances <- acceptance.probs <- denvec <- c()
	
	
	tarrec <- c()
	nsamp <- 0

    Ycount <- 1
	
	if(MCMCdiag>0){
		MCMCacc <- c()
		MCMCsamp <- c()
		ijcells <- cbind(sample(1:(M.ext-1),MCMCdiag,replace=TRUE),sample(1:(N.ext-1),MCMCdiag,replace=TRUE))
	}
	
	for(j in 1:time.block){
        EY.mean[[j]] <- EY.var[[j]] <- y.mean[[j]] <- y.var[[j]] <- gam.sum[[j]] <- y.sum[[j]] <- 0
        if(MCMCdiag>0){
        	jdmats[[j]] <- matrix(0,M.ext-1,N.ext-1)
        }
    }
    
    if (!is.null(inits)){
        cat("Initialising MALA using supplied initial values ...\n")
    	GAM <- inits
    }
    else{
        cat("Initialising MALA ...\n")            
        for(j in 1:time.block){
            GAM[[j]] <- matrix(0,M.ext-1,N.ext-1)
        }
    }
    cat("MALA successfully initialised ... please wait\n")   
    ## MALA now initialised
    
    GAMmean <- GAM
    
    oldtargrad <- target.and.grad(	at.mat.list=GAM,
                    				theta=THETA,
                    				varvals=vars,
                    				nicounts.mat.list=nis,
                    				timediffs=tdiffs,
                    				car.mat=cellArea.mat,
                    				mu=MU,
                    				spatial.mat=spatialvals,
                    				temporal=temporal.fitted,
                    				lam.tilde=LAM.tildefft,
                    				cin=cellInside,
                    				M.ext,
                    				N.ext,
                    				gradtrunc=gradtrunc)
    
                                       			
	while(nextStep(mcmcloop)){
	
        propmeans <- proposal.means(state.mat.list=GAM,
									h.var=h,
									del=oldtargrad$del)    									
		next.candidate <- list()
	    for(j in 1:length(GAM)){
	        next.candidate[[j]] <- matrix(rnorm((M.ext-1)*(N.ext-1),mean=as.vector(propmeans[[j]]),sd=sqrt(h)),M.ext-1,N.ext-1)
	    }       
                
        							
	
		newtargrad <- target.and.grad(	at.mat.list=next.candidate,
	                    				theta=THETA,
	                    				varvals=vars,
	                    				nicounts.mat.list=nis,
	                    				timediffs=tdiffs,
	                    				car.mat=cellArea.mat,
	                    				mu=MU,
	                    				spatial.mat=spatialvals,
	                    				temporal=temporal.fitted,
	                    				lam.tilde=LAM.tildefft,
	                    				cin=cellInside,
	                    				M.ext,
	                    				N.ext,
	                    				gradtrunc=gradtrunc)
	                    				
	                   				 
	    
	    # compute acceptance probability                				
	    
        revpropmeans <- proposal.means(	state.mat.list=next.candidate,
										h.var=h,
										del=newtargrad$del)
		
		# note that the object ac, defined on the next line is accessed by other functions via the parent.frame function																			
		ac <- exp(newtargrad$target - oldtargrad$target + (-1/(2*h))*diffnorm2(GAM,revpropmeans) - (-1/(2*h))*diffnorm2(next.candidate,propmeans))
		
		trigger <- FALSE
		if (newtargrad$target==-Inf | is.na(ac) | is.nan(ac)){ # gradient truncation insufficient, so reduce
	        gradtrunc <- gradtrunc/2
	        cat("Reducing gradient truncation to:",gradtrunc,"\n")
	        oldtargrad <- target.and.grad(	at.mat.list=GAM,
                    				theta=THETA,
                    				varvals=vars,
                    				nicounts.mat.list=nis,
                    				timediffs=tdiffs,
                    				car.mat=cellArea.mat,
                    				mu=MU,
                    				spatial.mat=spatialvals,
                    				temporal=temporal.fitted,
                    				lam.tilde=LAM.tildefft,
                    				cin=cellInside,
                    				M.ext,
                    				N.ext,
                    				gradtrunc=gradtrunc)
	        
	        if (!is.burnin(mcmcloop)){
	            cat("Gradient truncation currently",gradtrunc,"\n")
	            cat("Suggest reducing this further and setting gradtrunc manually in lgcpPredict.\n")
	            stop(paste("Problem with gradient truncation at iteration",iteration(mcmcloop),"acceptance probability =",ac))
	        }
	        ac <- 0 # don't accept the move if a suitable gradient truncation has not been found.
	        trigger <- TRUE # this is set to true so that the adaptive scheme is halted for one iteration if a suitable gradient trunctation has not been found
	    } 

	                       	
		
		if (ac>1){
			ac <- 1
		}							          
	    
	    
        if (MCMCdiag>0){
            if (is.retain(mcmcloop)){
                MCMCacc <- c(MCMCacc,ac)
                MCMCsamp <- rbind(MCMCsamp,as.vector(apply(ijcells,1,function(ij){GAM[[time.block]][ij[1],ij[2]]})))               
            }
        }
	    
	    if(ac>runif(1)){
	    	GAM <- next.candidate
	    	oldtargrad <- newtargrad
	    }
    	
        if (iteration(mcmcloop)>1){
	        hrec <- c(hrec,h)
	    }
	    
	    if (!trigger){ # ie if there was no problem with gradient truncation
	        h <- updateAMCMC(adaptivescheme)
        }	        
	    
	    ##print(h)
	    
	    ##print(c(ac,h,newtargrad$target))	
	    ##print(lapply(newtargrad$del,function(x){summary(as.vector(x))}))    
	    
	    if (is.retain(mcmcloop)){
    		nsamp <- nsamp + 1 # must be at start of loop so that exceedances are calculated correctly	    
	        ymats <- list()
	        ymats[[1]] <- (Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(sqrt(LAM.tildefft)*((1/sqrt((M.ext-1)*(N.ext-1)))*fft(GAM[[1]])),inverse=T)) + MU)[1:(M-1),1:(N-1)]
	        for(j in 2:time.block){
	            ymats[[j]] <- Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(sqrt(LAM.tildefft)*((1/sqrt((M.ext-1)*(N.ext-1)))*fft(GAM[[j]])),inverse=T))[1:(M-1),1:(N-1)] + at(tdiffs[j],MU,THETA)+bt.scalar(tdiffs[j],THETA)*ymats[[j-1]]
	        }
	        EY <- lapply(ymats,exp)
	        for(j in 1:time.block){
	        	y.mean[[j]] <- ((nsamp-1)/nsamp) * y.mean[[j]] + ymats[[j]]/nsamp
	        	EY.mean[[j]] <- ((nsamp-1)/nsamp) * EY.mean[[j]] + EY[[j]]/nsamp
	        	if (nsamp>1){
        			y.var[[j]] <- ((nsamp-2)/(nsamp-1))*y.var[[j]] + (nsamp/(nsamp-1)^2)*(y.mean[[j]]-ymats[[j]])^2
        			EY.var[[j]] <- ((nsamp-2)/(nsamp-1))*EY.var[[j]] + (nsamp/(nsamp-1)^2)*(EY.mean[[j]]-EY[[j]])^2
        		}                 
	        }
	        
			GFupdate(gridfun)
			GAupdate(gridav)		    
		}
				
	}
	
	retlist <- list(lasth=h,lastGAM=lgcpgrid(GAM))	
	
	GFfinalise(gridfun) # these two lines must appear after retlist has been initialised
	GAfinalise(gridav)  #  
	
	retlist$hrec <- hrec
	if(MCMCdiag>0){
		retlist$mcmcacc <- MCMCacc
		retlist$mcmcsamp <- MCMCsamp
		retlist$ijcells <- ijcells
	}
    retlist$y.mean <- lgcpgrid(y.mean)
    retlist$y.var <- lgcpgrid(y.var)
    retlist$EY.mean <- lgcpgrid(EY.mean)
    retlist$EY.var <- lgcpgrid(EY.var)
    retlist$gridfunction <- GFreturnvalue(gridfun)
    retlist$gridaverage <- GAreturnvalue(gridav)
    retlist$mcmcinfo <- mcmcloop
    retlist$gradtrunc <- gradtrunc	
	
	return(retlist)
}

