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
#' @param lower lower bound for \code{theta} testing range (default = 0) 
#' @param upper upper bound for \code{theta} testing range (default = 1)
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

esttheta <- function(traj,r,lower=0,upper=1,rand=NA,niter=10,tolerance=0.01,dmin=NA,dmax=NA,plot=TRUE){

  #Function to compute likelihood for a set of fixes and a given level of theta
  thetafunc <- function(theta,x,ii,tr,r,log){
    pz <- 0*ii
    for (i in 1:length(ii)){
      j <- ii[i]
      sp1 <- SpatialPoints(x[j,c('x','y')])
      sp2 <- SpatialPoints(x[j+2,c('x','y')])
      spz <- SpatialPoints(x[j+1,c('x','y')])
      c1 <- raster::cellFromXY(tr,sp1)
      c2 <- raster::cellFromXY(tr,sp2)
      cz <- raster::cellFromXY(tr, spz)
      chck <- anyDuplicated(c(c1,c2,cz))
      # if (c1 == c2){
      #   #Start and end pixel is the same which means no movement.
      #   # Arbitrarily set end location to the next pixel over (check if edge)
      #   ########################################
      #   ind <- adjacent(r,c1,pairs=FALSE,id=TRUE)
      #   c2 <- ind[which.max(r[ind])]
      #   sp2 <- raster::xyFromCell(r,c2,spatial=TRUE)
      #   sp2@proj4string <-sp1@proj4string   #this could cause an issue if traj and raster not in same projection
      # }
      if (chck == 0) {
        #Movement Occurs (at least from one cell to another)
        Pt <- passage(tr,sp1,sp2,theta=theta,totalNet='net')
        pz[i] <- Pt[cz]
      } else {
        pz[i] <- 0
      }

    }
    
    #Calculate the negative of the likelihood - we are using a minimizing golden search function
    LLpz <- -pz
    LL <- log(sum(LLpz,na.rm=T))
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
    ii <- ii[which(x$dist[ii] >= dmin & x$dist[ii+1] >= dmin)]   
    print(paste('Using a dmin value of',dmin, ' ; ', length(ii), 'fixes will be used to estimate theta.'))
  }

  #Can choose to only use non-movement fixes (based on dmax) - reduces number of segments in test.
  if (!is.na(dmax)){
    ii <- ii[which(x$dist[ii] <= dmax & x$dist[ii+1] <= dmax)]   
    print(paste('Using a dmax value of',dmax, ' ; ', length(ii), 'fixes will be used to estimate theta.'))
  }

  
  #Make *Symmetric* TransitionLayer from raster (e.g., avg cells)
  s1 <- function(x){x[1]}
  s2 <- function(x){x[2]}
  tr1 <- transition(r, s1, 8,symm=F)
  tr2 <- transition(r, s2, 8,symm=F)
  tr <- (tr1 + tr2)/2
  tr <- geoCorrection(tr,type='c',multpl=FALSE)
  r <- raster(tr)
  
  #### GOLDEN SEARCH ROUTINE ###
  #Progress Bar
  cat('Iterations: \n')

  golden.ratio = 2/(sqrt(5) + 1)
  
  ### Evaluate the function at the extremes
  fmin = thetafunc(lower,x,ii,tr,r)
  cat('1 \n')
  fmax = thetafunc(upper,x,ii,tr,r)
  cat('2 \n')
                   
  ### Use the golden ratio to set the initial test points
  x1 = upper - golden.ratio*(upper - lower)
  x2 = lower + golden.ratio*(upper - lower)
  
  ### Evaluate the function at the first test points
  f1 = thetafunc(x1,x,ii,tr,r)
  cat('3 \n')
  f2 = thetafunc(x2,x,ii,tr,r)
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
      f1 = thetafunc(x1,x,ii,tr,r)
      theta.val <- c(theta.val,x1)
      LL.val <- c(LL.val,f1)
    } else {
      lower = x1
      x1 = x2
      f1 = f2
      x2 = lower + golden.ratio*(upper - lower)
      f2 = thetafunc(x2,x,ii,tr,r)
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
    LL.val <- -log(LL.val)
    plot(theta.val[ord],LL.val[ord],xlab='theta',ylab='log-likelihood',type=n)
    ss <- smooth.spline(theta.val[ord],LL.val[ord],df=4)
    lines(ss)
    #points(theta.val[ord],-LL.val[ord],type='l')
    abline(v=est.min,col='red')
  }
  
  return(est.min)
}
