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
#' @param max upper bound for \code{c2} testing range (default = 1)
#' @param rand if \code{NA} (the default) every second segment is evauluated (n/2), otherwise an integer indicating how many random segments to test.
#' @param niter used to define maximum number of iterations of golden-search routine
#' @param tolerance used to define precision of golden search routine (i.e., routine stops when the absolute difference between two consecutive test points is below this value). 
#' @param plot logical, whether or not to plot the log-likelihood curve.
#' 
#' @return
#'   This function returns a numerical value estimate for \code{c2} associated with the maximum of the log-likelihood.
#'
# @references
# @keywords 
#' @seealso fbtgUD
# @examples
#' @importFrom graphics abline points
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom igraph E V distances graph.adjacency
#' @export
# ---- End of roxygen documentation ----
estc2 <- function(traj,tl,timefun='exp',sigma=0,min=0,max=1,rand=NA,niter=10,tolerance=0.01,plot=TRUE){
  
  #use leave-one-out bootstrap to estimate theta in a similar fashion to proposed by Horne et al. 2007 as is commonly used with Brownian bridge, can speed up by chosing smaller number of random segments to test.
  x <- ld(traj)
  n <- dim(x)[1]
  if (is.na(rand)){
    ii <- seq(1,(n-2),by=2)
  } else {
    ii <- sample(1:(n-2),rand)
  }

  tm <- transitionMatrix(tl)
  gr <- graph.adjacency(tm, mode="directed", weighted=TRUE)
  E(gr)$weight <- 1/E(gr)$weight
  
  
  ### Golden Search Routine
  
  #internal likelihood funciton
  c2func <- function(c2,ii,gr,tl,sigma,timefun){
    A <- SpatialPoints(x[ii,c('x','y')])
    B <- SpatialPoints(x[ii+1,c('x','y')])
    C <- SpatialPoints(x[ii+2,c('x','y')])
    Ai <- cellFromXY(raster(tl),A)
    Bi <- cellFromXY(raster(tl),A)
    Ci <- cellFromXY(raster(tl),A)
    
    Tshort <- diag(costDistance(tl,A,C))
    Tai <- distances(gr,v=Ai,to=V(gr),mode='out')
    Tib <- distances(gr,v=Ci,to=V(gr),mode='in')
    t1 <- x$dt[ii]
    t2 <- x$dt[ii+1]
    tt <- t1+t2
    pz <- 0*ii
    for (i in 1:length(ii)){
      #Compute the likelihood for each segment
      pz[i] <- internalTS(t1[i],tt[i],Tai[i,],Tib[i,],tl,Tshort[i],sigma=sigma,timefun,c2,clipPPS=FALSE)[Bi[i]]
    }
    #Calculate the negative of the log-likelihood - we are using a minimizing golden search function
    LLpz <- -log(pz) 
    #REMOVE -INFs
    LLpz[is.infinite(LLpz)] <- NA
    LL <- sum(LLpz,na.rm=T)
    return(LL)
  }
  
  #### GOLDEN SEARCH ROUTINE ###
  #Progress Bar
  cat('Iterations: \n')
  
  golden.ratio = 2/(sqrt(5) + 1)
  
  upper.bound <- max
  lower.bound <- min
  
  ### Evaluate the function at the extremes
  fmin = c2func(min,ii,gr,tl,sigma,timefun)
  cat('1 \n')
  fmax = c2func(max,ii,gr,tl,sigma,timefun)
  cat('2 \n')
  
  ### Use the golden ratio to set the initial test points
  x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  
  ### Evaluate the function at the first test points
  f1 <- c2func(x1,ii,gr,tl,sigma,timefun)
  cat('3 \n')
  f2 <- c2func(x2,ii,gr,tl,sigma,timefun)
  cat('4 \n')
  ### Output values storage
  c2.val <- c(min,max,x1,x2)
  LL.val <- c(fmin,fmax,f1,f2)
  
  #Search
  iteration = 0
  #Progress Bar
  
  while (iteration < (niter-4) & abs(upper.bound - lower.bound) > tolerance){
    iteration = iteration + 1
    if (f2 > f1){
      upper.bound = x2
      x2 = x1
      f2 = f1
      x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
      f1 = c2func(x1,ii,gr,tl,sigma,timefun)
      c2.val <- c(c2.val,x1)
      LL.val <- c(LL.val,f1)
    } else {
      lower.bound = x1
      x1 = x2
      f1 = f2
      x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
      f2 = c2func(x2,ii,gr,tl,sigma,timefun)
      c2.val <- c(c2.val,x2)
      LL.val <- c(LL.val,f2)
    }
    #update progress bar
    est.min = (lower.bound + upper.bound)/2
    cat(paste(iteration+4,' c2.est = ',est.min,'\n'))
  }
  #est.min = (lower.bound + upper.bound)/2

  if (plot){
    ord <- order(c2.val)
    plot(c2.val[ord],-LL.val[ord],xlab='theta',ylab='log-likelihood',pch=20)
    points(c2.val[ord],-LL.val[ord],type='l')
    abline(v=est.min,col='red')
  }

  return(est.min)
}



