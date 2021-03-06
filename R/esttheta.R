# ---- roxygen documentation ----
#' @title Numerically estimate theta parameter for rspUD
#'
#' @description
#'   Estimate theta parameter using leave-one-out estimation procedure for randomized shortest path utilization distributions.
#' @details
#' Still under development, currently it is WAY TO SLOW, looking at Panzachi and Saerens to see how to speed up (i.e., by not using passage function).
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param r a \code{RasterLayer} object describing the preference/permeability of the landscape. May be the result of a resource selection function, or other analyses. Note: Higher values should be associated higher preference or permeability.
#' @param rangetheta uppper and lower bound for \code{theta} testing range (default = [0,1]).
#' @param rand if \code{NA} (the default) every second segment is evauluated (n/2), otherwise an integer indicating how many random segments to test.
#' @param niter used to define maximum number of iterations of golden-search routine
#' @param tolerance used to define precision of golden search routine (i.e., routine stops when the absolute difference between two consecutive test points is below this value). 
#' @param dmin Use only segments where the movement distance is greater than dmin (default is NA or all segments).
#' @param dmax Use only segments where the movement distance is less than dmax (default is NA or all segments).
#' @param plot logical, whether or not to plot the likelihood curve.
#'
#' @return
#'   This function returns an estimate for the theta parameter which can be used with the rspUD function.
#'
# @references
# @keywords 
#' @seealso rspUD
# @examples
#' @importFrom graphics abline points
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#
# ---- End of roxygen documentation ----

esttheta <- function(traj,r,rangetheta=c(0,1),rand=NA,niter=10,tolerance=0.01,dmin=NA,dmax=NA,plot=TRUE){

  #Function to compute likelihood for a set of fixes and a given level of theta
  thetafunc <- function(theta,x,ii,tr){
    pz <- 0*ii
    w <- 0*ii
    mdt <- mean(x$dt[ii]+x$dt[ii+1])
    for (i in 1:length(ii)){
      j <- ii[i]
      sp1 <- SpatialPoints(x[j,c('x','y')])
      sp2 <- SpatialPoints(x[j+2,c('x','y')])
      spz <- SpatialPoints(x[j+1,c('x','y')])
      c1 <- raster::cellFromXY(tr,sp1)
      c2 <- raster::cellFromXY(tr,sp2)
      cz <- raster::cellFromXY(tr, spz)
      wi <- (mdt/(x$dt[j]+x$dt[j+1]))^2
      #It only makes sense to cross validate where Movement Occurs (at least from one cell to another)
      chck <- anyDuplicated(c(c1,c2,cz))
      if (chck == 0) {
        Pt <- passage(tr,sp1,sp2,theta=wi*theta,totalNet='net')
        Pt[c(c1,c2)] <- 0.5
        Pt <- Pt / cellStats(Pt,sum)    
        pz[i] <- Pt[cz]
      } else {
        pz[i] <- 1
      }

    }
    
    #Calculate the negative log likelihood - we are using a minimizing golden search function
    LL <- -sum(log(pz))
    return(LL)
  }
  
  
  #use leave-one-out bootstrap to estimate theta in a similar fashion to proposed by Horne et al. 2007 as is commonly used with Brownian bridge, can speed up by chosing smaller number of random segments to test.
  x <- ld(traj)
  n <- dim(x)[1]
  if (is.na(rand)){
    ii <- seq(1,(n-2),by=2)
  } else {
    ii <- sample(1:(n-2),rand)
  }
  #Can choose to only use movement fixes (based on dmin) - reduces number of segments in test. 
  if (!is.na(dmin)){
    ii <- ii[which(x$dist[ii] + x$dist[ii+1] >= dmin)]   
    print(paste('Using a dmin value of',dmin, ' ; ', length(ii), 'fixes will be used to estimate theta.'))
  }

  #Can choose to only use non-movement fixes (based on dmax) - reduces number of segments in test.
  if (!is.na(dmax)){
    ii <- ii[which(x$dist[ii] + x$dist[ii+1] <= dmax)]   
    print(paste('Using a dmax value of',dmax, ' ; ', length(ii), 'fixes will be used to estimate theta.'))
  }

  
  #Make *Symmetric* TransitionLayer from raster (e.g., avg cells)
  s1 <- function(x){x[1]}
  s2 <- function(x){x[2]}
  tr1 <- transition(r, s1, 8,symm=F)
  tr2 <- transition(r, s2, 8,symm=F)
  tr <- (tr1 + tr2)/2
  tr <- geoCorrection(tr,type='c',multpl=FALSE)

  
  #### GOLDEN SEARCH ROUTINE ###
  lower <- rangetheta[1]
  upper <- rangetheta[2]
  #Progress Bar
  cat('Iterations: \n')

  golden.ratio = 2/(sqrt(5) + 1)
  
  ### Evaluate the function at the extremes
  fmin = thetafunc(lower,x,ii,tr)
  cat('1 \n')
  fmax = thetafunc(upper,x,ii,tr)
  cat('2 \n')
                   
  ### Use the golden ratio to set the initial test points
  x1 = upper - golden.ratio*(upper - lower)
  x2 = lower + golden.ratio*(upper - lower)
  
  ### Evaluate the function at the first test points
  f1 = thetafunc(x1,x,ii,tr)
  cat('3 \n')
  f2 = thetafunc(x2,x,ii,tr)
  cat('4 \n')
  ### Output values storage
  theta.val <- c(lower,upper,x1,x2)
  LL.val <- c(fmin,fmax,f1,f2)
  
  #Search
  iteration = 0
  #Progress Bar

  while (iteration < (niter-4) & abs(upper - lower) > tolerance){
    iteration = iteration + 1
    if (f2 > f1){
      upper = x2
      x2 = x1
      f2 = f1
      x1 = upper - golden.ratio*(upper - lower)
      f1 = thetafunc(x1,x,ii,tr)
      theta.val <- c(theta.val,x1)
      LL.val <- c(LL.val,f1)
    } else {
      lower = x1
      x1 = x2
      f1 = f2
      x2 = lower + golden.ratio*(upper - lower)
      f2 = thetafunc(x2,x,ii,tr)
      theta.val <- c(theta.val,x2)
      LL.val <- c(LL.val,f2)
    }
    #update progress bar
    est.min = (lower + upper)/2
    cat(paste(iteration+4,' theta.est = ',est.min,'\n'))
  }
  #est.min = (lower + upper)/2
  
  if (plot){
    ord <- order(theta.val)
    plot(theta.val[ord],-LL.val[ord],xlab='theta',ylab='log-likelihood')
    points(theta.val[ord],-LL.val[ord],type='l')
    abline(v=est.min,col='red')
  }
  
  return(est.min)
}


