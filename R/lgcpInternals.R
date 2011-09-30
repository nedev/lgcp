##' at function
##' @param t change in time parameter, see Brix and Diggle (2001)
##' @param mu mean
##' @param theta parameter beta in Brix and Diggle
##' @return ...
at <- function(t,mu,theta){
    return(mu*(1-exp(-t*theta)))
}



##' bt.scalar function
##' @param t change in time, see Brix and Diggle (2001)
##' @param theta parameter beta in Brix and Diggle
##' @return ...
bt.scalar <- function(t,theta){
    return(exp(-t*theta))
}



##' d.func function
##' @param mat1il matrix 1
##' @param mat2jk matrix 2
##' @param i index matrix 1 number 1
##' @param j index matrix 2 number 1
##' @param l index matrix 1 number 2
##' @param k index matrix 2 number 2
##' @return ...
d.func <- function(mat1il,mat2jk,i,j,l,k){
# warning global variables
  return(sqrt(mat1il[i,l]^2+mat2jk[j,k]^2))
}



##' diffnorm2 function
##' @param list1 list 1
##' @param list2 list 2
##' @return ...
diffnorm2 <- function(list1,list2){
    result <- 0
    for(i in 1:length(list1)){
        result <- result + sum((list1[[i]]-list2[[i]])^2)
    }
    return(result)
}



##' evallogprobY function
##' @param mat.list list of Y matrices 
##' @param at (vector) component a_t from Brix and Diggle
##' @param bt (vector) component b_t from Brix and Diggle
##' @param gmult (vector) (1-exp(-2*beta*delta_t)) as in G(T); Brix and Diggle
##' @param logdetcovY log determinant of R matrix as in G(T); Brix and Diggle
##' @param Rinv inverse of R as in G(T); Brix and Diggle
##' @return ...
evallogprobY <- function(mat.list,at,bt,gmult,logdetcovY,Rinv){
    n <- dim(Rinv)[1]
    v <- as.vector(mat.list[[1]]) - rep(at[1],n)
    logpr <- -0.5*logdetcovY - 0.5*t(v)%*%Rinv%*%v
    if(length(mat.list)>=2){
        for(i in 2:length(mat.list)){
            v <- as.vector(mat.list[[i]]) - rep(at[i],n) - rep(bt[i],n)*as.vector(mat.list[[i-1]])
            logpr <- logpr -0.5*logdetcovY - 0.5*t(v)%*%((1/gmult[i])*Rinv)%*%v 
        }
    }
    return(logpr)
}



##' grad.comp function
##' @param m parameter
##' @param i index
##' @param theta beta parameter in Brix and Diggle (2001)
##' @param varvals pre computed variances
##' @param nicounts.mat.list matrix of nis
##' @param timediffs vector of delta ts
##' @param car.mat matrix
##' @param ys Gaussian field
##' @param spatial fixed spatial component
##' @param temporal fixed temporal component
##' @return gradient
grad.comp <- function(	m,
	                    i,
	                    theta,
	                    varvals,
	                    nicounts.mat.list,
	                    timediffs,
	                    car.mat,
	                    ys,
	                    spatial,
	                    temporal){
    	                             	
    p1<-bt.scalar(timediffs[m],theta)
    p2 <- nicounts.mat.list[[i]]-car.mat*exp(ys[[i]])*spatial*temporal[i]
    return(list(result=p1*p2,resids=p2))
}




##' gu function
##' @param u distance
##' @param sigma variance parameter, see Brix and Diggle (2001)
##' @param phi scale parameter, see Brix and Diggle (2001)
##' @param model correlation type, see ?CovarianceFct
##' @param additionalparameters vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @return this is just a wrapper for CovarianceFct
gu <- function(u,sigma,phi,model,additionalparameters){
        return(CovarianceFct(x=u,param=c(mean=0,variance=sigma^2,nugget=0,scale=phi,additionalparameters),model=model))

}



##' meanselect function
##' @param mymatrix matrix
##' @param selector selection
##' @return ...
meanselect <- function(mymatrix,selector){
    return(mean(mymatrix[selector]))
}



##' normsq function
##' @param mymatrix matrix
##' @param selector selection
##' @param M number of x 'gridlines'
##' @param N number of y 'gridlines'
##' @return ...
normsq <- function(mymatrix,selector,M,N){
    return(sum(mymatrix[1:(M-1),1:(N-1)][selector]^2))
}



##' proposal.means function
##' @param state.mat.list list of matrices
##' @param h.var variance of MALA proposal
##' @param del gradient
##' @return ...
proposal.means <- function(	state.mat.list,
			                h.var,
			                del){               	
    xilist <- list()
    for(i in 1:length(state.mat.list)){
        xilist[[i]] <- as.vector(state.mat.list[[i]]+0.5*h.var*del[[i]])
    }
    return(xilist)
}



##' target.and.grad function
##' @param at.mat.list list of matrices
##' @param theta parameter beta, see Brix and Diggle (2001)
##' @param varvals variances, pre-computed 
##' @param nicounts.mat.list list of matrices of ni counts
##' @param timediffs vector of delta ts
##' @param car.mat matrix
##' @param mu mean parameter, see Brix and Diggle
##' @param spatial.mat fixed spatial component, output from function fftgrid
##' @param temporal fixed temporal component
##' @param lam.tilde transpose of Fourier transform on extended grid, output from function fftgrid
##' @param cin cellInside, output from function fftgrid
##' @param M.ext number of x 'gridlines' in extended space
##' @param N.ext number of y 'gridlines' in extended space
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Set to zero for no gradient truncation (default).
##' @param estgradtrunc set to true if estimating gradient truncation, see ?computeGradtrunc
##' @return The log target and gradient
target.and.grad <- function(at.mat.list,
    		                theta,
    		                varvals,
    		                nicounts.mat.list,
    		                timediffs,
    		                car.mat,
    		                mu,
    		                spatial.mat,
    		                temporal,
    		                lam.tilde,
    		                cin,
    		                M.ext,
    		                N.ext,
    		                gradtrunc,
    		                estgradtrunc=FALSE){
    
    ys <- list()
    
    spatial.mat.subset <- spatial.mat[1:((M.ext-1)/2),1:((N.ext-1)/2)]
    car.mat.subset <- car.mat[1:((M.ext-1)/2),1:((N.ext-1)/2)]
    
    y.dummy <- sqrt(lam.tilde)*((1/sqrt((M.ext-1)*(N.ext-1)))*fft(at.mat.list[[1]]))
    ys[[1]] <- Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(y.dummy,inverse=T)) + mu

    for(i in 2:length(at.mat.list)){
        y.dummy <- sqrt(lam.tilde)*((1/sqrt((M.ext-1)*(N.ext-1)))*fft(at.mat.list[[i]]))
        ys[[i]] <- Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(y.dummy,inverse=T)) +at(timediffs[i],mu,theta)+bt.scalar(timediffs[i],theta)*ys[[i-1]]
    }

	# compute target
	    f1 <- 0
	    for(i in 1:length(at.mat.list)){
	        determ <- spatial.mat.subset*temporal[i] #'spatial' is now M*N interpolated density matrix
	        determ <- log(determ)
	        determ[is.infinite(determ)] <- 0
	        f1 <- f1 + sum((ys[[i]][1:((M.ext-1)/2),1:((N.ext-1)/2)][cin]+determ[cin])*nicounts.mat.list[[i]][1:((M.ext-1)/2),1:((N.ext-1)/2)][cin]-car.mat.subset[cin]*exp(ys[[i]][1:((M.ext-1)/2),1:((N.ext-1)/2)][cin])*exp(determ[cin]))  #fark... how to deal with deterministic grid...?
	    }
	    f2 <- 0
	    for(i in 2:length(at.mat.list)){
	        f2 <- f2 + (0.5/varvals[i])*sum(at.mat.list[[i]]^2)
	    }    
	    target <- f1-f2-0.5*sum(at.mat.list[[1]]^2)
	    
	# compute gradient
	    del <- gradchunkresidslist <- list()
	    for(j in 1:length(at.mat.list)){
	        del[[j]] <- (-1/varvals[j])*at.mat.list[[j]]
	        gradchunk <- gradchunkresids <- matrix(0,(M.ext-1),(N.ext-1))        
	        for(i in j:length(at.mat.list)){
	            gradbit <- gradresids <- matrix(1,(M.ext-1),(N.ext-1))
	            if(j+1>i){
	            	next
	            }
	            for(m in (j+1):i){
	                gc <- grad.comp(m=m,
	                                 i=i,
	                                 theta=theta,
	                                 varvals=varvals,
	                                 nicounts.mat.list=nicounts.mat.list,
	                                 timediffs=timediffs,
	                                 car.mat=car.mat,
	                                 ys=ys,
	                                 spatial=spatial.mat,
	                                 temporal=temporal)
	
	                componentQ <- sqrt(lam.tilde)*Re(((1/sqrt((M.ext-1)*(N.ext-1)))*fft(gc$result)))
	                componentQ <- Re((1/sqrt((M.ext-1)*(N.ext-1)))*fft(componentQ,inverse=T))              
	                gradbit <- gradbit*componentQ
	                ##print(c(j,i)) # debugging purposes
	                ##print(summary(as.vector(gc$result)))
	            	if(j==1){
	            		gradresids <- gradresids*gc$resids
	            	}
	            }
	            gradchunk <- gradchunk+gradbit
	            gradchunkresids <- gradchunkresids+gradresids
	        }
	        gradchunkresidslist[[j]] <- gradchunkresids
	        del[[j]] <- del[[j]]+gradchunk#t(t(gradchunk)%*%t(Q))
	    }
	    
	if (estgradtrunc){
	    return(max(sapply(del,max,na.rm=TRUE),na.rm=TRUE))
	}    
	
	for(i in 1:length(del)){ # now truncate gradient
	    del[[i]][del[[i]]>gradtrunc] <- gradtrunc
	}       
        
    return(list(del=del,resids=gradchunkresidslist,target=target)) # return of "resids=gradchunkresidslist" appears to be redundant in reset of code ??
}
