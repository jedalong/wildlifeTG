# ---- roxygen documentation ----
#' @title Numerically estimate c2 parameter
#'
#' @description
#'   Numerically estimate the \code{c2} parameter for \code{fbtgUD} using a data-driven approach.
#' @details
#' THe estimation of the \code{c2} parameter for \code{fbtgUD} model takes an identical approach to that used with Brownian bridges, as originally proposed by Horne et al. (2007). A leave-one-out estimation technique is used, which essentially removes fixes and then estimates the fbtgUD surface between the two adjacent fixes  -- termed a segment -- and computes the probability associated with the missing fix. The \code{c2} value is returned that numerically maximizes the log-likelihood (for a given \code{timefun} -- see details in \code{fbtgUD}) for \code{length.out} evenly spaced values of \code{c2} between (user-defined) \code{min} and \code{max}. With \code{plot = TRUE} the user can verify that the chosen range of potential \code{c2} values is appropriate, and if not, retry using a differnt range. The process is computationally demanding and is highly dependent on the number of 'segments' used. The parameter \code{rand} can be used to adjust the number of segments used to minimize computational time. If \code{rand = NA} it removes every second fix, and estimates \code{c2} based on these n/2 segments (the default). Otherwise \code{rand} can be passed in as an integer, and \code{rand} randomly selected segments will be chosen to estimate \code{c2}. This is beneficial for trajectories with many fixes, where it might be useful to choose \code{rand <<< n/2}. Parrallelization is possible to further decrease computational time. This can be implemented by choosing an appropriate integer value for the \code{parallel} parameter (implemented using the \code{foreach} package). 
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param tl a \code{TransitionLayer} object
#' @param timefun method for converting time into probability (see \code{fbtgUD}); one of:\cr 
#' -- \code{'inverse'} \eqn{= \frac{1}{c_2 * t}},\cr
#' -- \code{'inverse2'} \eqn{= \frac{1}{(c_2 * t)^2}},\cr
#' -- \code{'exp'} \eqn{= \exp{-c_2 * t}},\cr
#' -- \code{'norm'} \eqn{= \exp{-c_2 * t^2}},\cr
#' -- \code{'rootexp'} \eqn{= \exp{-c_2 * \sqrt{t}}},\cr
#' -- \code{'pareto'} \eqn{= \exp{-c_2 * \log{t}}},\cr
#' -- \code{'lognorm'}\eqn{= \exp{-c_2 * \log{t^2}}}.\cr
#' @param sigma location uncertainty parameter (see fbtgUD)
#' @param min lower bound for \code{c2} testing range (default = 0) 
#' @param min upper bound for \code{c2} testing range (default = 1)
#' @param length.out used to define precision of \code{c2} values tested, using \code{seq(min,max,length.out=length.out)}
#' @param rand if \code{NA} (the default) every second segment is evauluated (n/2), otherwise an integer indicating how many random segments to test.
#' @param parallel the number of parralel processors to use (see the \code{foreach} package), default = 1 (Not Used).
#' @param plot logical, whether or not to plot the log-likelihood curve.
#' 
#' @return
#'   This function returns a numerical value estimate for \code{c2} associated with the maximum of the log-likelihood.
#'
# @references
# @keywords 
#' @seealso fbtgUD
# @examples
#' @export
#
# ---- End of roxygen documentation ----
likec2 <- function(traj,tl,timefun,sigma=0,min=0,max=1,length.out=50,rand=NA,parallel=1,plot=TRUE){
  
  pckgs <- c('gdistance')
  
  #internal likelihood funciton 
  PiFun <- function(c2,x,zi,timefun){
    Pt <- internalPfun(x,timefun,c2)
    c1 <- 1/sum(Pt)
    Pi <- internalPfun(zi,timefun,c2)
    Pi <- c1*Pi
    return(Pi)
  }
  
  #use leave-one-out bootstrap to estimate theta in a similar fashion to proposed by Horne et al. 2007 as is commonly used with Brownian bridge, can speed up by chosing smaller number of random segments to test.
  x <- ld(traj)
  n <- dim(x)[1]
  if (is.na(rand)){
    ii <- seq(1,(n-2),by=2)
  } else {
    ii <- sample(1:(n-2),rand)
  }

  c2. <- seq(min,max,length.out=length.out)
  pi.mat <- matrix(nrow=length(ii),ncol=length.out)
  
  #Parallelize
  #cl<-makeCluster(parallel)
  #registerDoParallel(cl)
  
  tm <- transitionMatrix(tl)
  gr <- graph.adjacency(tm, mode="directed", weighted=TRUE)
  E(gr)$weight <- 1/E(gr)$weight
  
  #pi.mat <- foreach(i=1:length(ii),.combine=rbind,.packages=pckgs) %dopar% {
  for (i in 1:length(ii)){
    j <- ii[i]
    
    A <- as.numeric(x[j,c('x','y')])
    B <- as.numeric(x[j+1,c('x','y')])
    C <- as.numeric(x[j+2,c('x','y')])
    
    Ai <- cellFromXY(raster(tl),A)
    Bi <- cellFromXY(raster(tl),B)
    Ci <- cellFromXY(raster(tl),C)
    
    Tshort <- costDistance(tl,A,C)[1]
    
    #Compute Accumulated cost (i.e, time) for location A and B
    #This is the value input for various functions
    Tai <- distances(gr,v=Ai,to=V(gr),mode='out')
    Tib <- distances(gr,v=Ci,to=V(gr),mode='in')
    
    t1 <- x$dt[j]
    t2 <- x$dt[j+1]
    
    tt <- t1+t2
    
    #Compute the likelihood #highly inefficient nesting of for loops, how to fix?
    clipPPS=FALSE   #For speed reasons always use clipPPS=FALSE!
    for (k in 1:length(c2.)){
      pi.mat[i,k] <- internalTS(t1,tt,Tai,Tib,tl,Tshort,sigma=sigma,timefun,c2.[k],clipPPS)[Bi]
    }
    #compute the likelihood
    #pi.mat[i,] <- sapply(c2.,PiFun,Ti,Zi,timefun)
  }
  #stopCluster(cl)
  
  #Calculate the log-likelihood
  pi.mat <- log(pi.mat)
  LL <- colSums(pi.mat)
  ###Calculate max value
  imax <- which.max(LL)   #maximum value observed - good enough for this purpose
  
  if (plot){
    plot(c2.,LL,xlab='c2',ylab='log-likelihood',type='l')
    points(c2.[imax],LL[imax],pch=20,col='red')
    abline(v=c2.[imax],col='red')
  }

  return(c2.[imax])
}



