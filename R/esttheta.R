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
#' @param plot logical, whether or not to plot the log-likelihood curve.
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

esttheta <- function(traj,r,lower=0,upper=1,rand=NA,niter=10,tolerance = 0.01,plot=TRUE){
  

  #Function to alter tran matrix from gdistance
  tranSolid <- function(x){
    selection <- which(Matrix::rowMeans(gdistance::transitionMatrix(x,inflate=FALSE))>1e-300)
    x@transitionCells <- x@transitionCells[selection]
    x@transitionMatrix <- gdistance::transitionMatrix(x,inflate=FALSE)[selection,selection]
    return(x)
  }
  
  #Function to compute LL probability for a set of fixes and a given level of theta
  ### MODIFIED FROM passage function in gdistance
  thetafunc <- function(theta,x,ii,tm,tc,trR,P,Id,nr,nc){
    pz <- 0*ii
    W <- trR
    W@x <- exp(-theta * trR@x) #zero values are not relevant because of next step
    W <- W * P

    for (i in 1:length(ii)){
      j <- ii[i]
      a <- SpatialPoints(x[j,c('x','y')])
      b <- SpatialPoints(x[j+1,c('x','y')])
      c <- SpatialPoints(x[j+2,c('x','y')])
      cellnri <- raster::cellFromXY(tm, a)
      cellnrz <- raster::cellFromXY(tm, b)
      cellnrj <- raster::cellFromXY(tm, c)
      ci <- match(cellnri,tc)
      cz <- match(cellnrz,tc)
      cj <- match(cellnrj,tc)
      Ij <- Diagonal(nr)
      Ij[cbind(cj,cj)] <- 1 - 1 / length(cj)
      Wj <- Ij %*% W
      ei <- rep(0,times=nr)
      ei[ci] <- 1 / length(ci)
      ej <- rep(0,times=nr)
      ej[cj] <- 1 / length(cj)
      IdMinusWj <- as((Id - Wj), "dgCMatrix")
      zci <- Matrix::solve(t(IdMinusWj),ei) 
      zcj <- Matrix::solve(IdMinusWj, ej)
      zcij <- sum(ei*zcj)
      if(zcij < 1e-300){
        n <- rep(0,times=nc)
      } else {
        N <- (Diagonal(nr, as.vector(zci)) %*% Wj %*% Diagonal(nr, as.vector(zcj))) / zcij
        
        # #totalnet = 'total'
        # n <- pmax(rowSums(N),colSums(N)) #not efficient but effective
        
        #totalnet = 'net'
        nNet <- abs(skewpart(N))
        n <- pmax(rowSums(nNet),colSums(nNet))
        n[c(ci,cj)] <- 2 * n[c(ci,cj)]
      }
      pz[i] <- n[cz]  
    }
    #Calculate the negative of the log-likelihood - we are using a minimizing golden search function
    LLpz <- -log(pz) 
    #REMOVE -INFs
    LLpz[is.infinite(LLpz)] <- NA
    LL <- sum(LLpz,na.rm=T)
    return(LL)
  }
  
  
  #use leave-one-out bootstrap to estimate theta in a similar fashion to proposed by Horne et al. 2007 as is commonly used with Brownian bridge, can speed up by chosing smaller number of random segments to test.
  x <- ld(traj)
  n <- dim(x)
  if (is.na(rand)){
    ii <- seq(1,(n-2),by=2)
  } else {
    ii <- sample(1:(n-2),rand)
  }
  
  #dmin=res(r)[1]
  #ii <- ii[which(x$dist[ii] > dmin)]   #Only use movement fixes 
  
  #Make *Symmetric* TransitionLayer from raster (e.g., avg cells)
  s1 <- function(x){x[1]}
  s2 <- function(x){x[2]}
  tr1 <- transition(r, s1, 8,symm=F)
  tr2 <- transition(r, s2, 8,symm=F)
  tm <- (tr1 + tr2)/2
  tm <- geoCorrection(tm,type='c',multpl=FALSE)
  
  #Prepare transition matrix for PFUN
  ts <- tranSolid(tm)
  tc <- transitionCells(ts)
  tr <- transitionMatrix(ts,inflate=FALSE)
  trR <- tr
  trR@x <- 1 / trR@x 
  nr <- dim(tr)[1] 
  Id <- Diagonal(nr) 
  rs <- rowSums(tr)
  rs[rs>0] <- 1/rs[rs>0]
  P <- tr * rs
  nc <- ncell(ts)
  
  #### GOLDEN SEARCH ROUTINE ###
  #Progress Bar
  cat('Iterations: \n')

  golden.ratio = 2/(sqrt(5) + 1)
  
  ### Evaluate the function at the extremes
  fmin = thetafunc(lower,x,ii,tm,tc,trR,P,Id,nr,nc)
  cat('1 \n')
  fmax = thetafunc(upper,x,ii,tm,tc,trR,P,Id,nr,nc)
  cat('2 \n')
                   
  ### Use the golden ratio to set the initial test points
  x1 = upper - golden.ratio*(upper - lower)
  x2 = lower + golden.ratio*(upper - lower)
  
  ### Evaluate the function at the first test points
  f1 = thetafunc(x1,x,ii,tm,tc,trR,P,Id,nr,nc)
  cat('3 \n')
  f2 = thetafunc(x2,x,ii,tm,tc,trR,P,Id,nr,nc)
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
      f1 = thetafunc(x1,x,ii,tm,tc,trR,P,Id,nr,nc)
      theta.val <- c(theta.val,x1)
      LL.val <- c(LL.val,f1)
    } else {
      lower = x1
      x1 = x2
      f1 = f2
      x2 = lower + golden.ratio*(upper - lower)
      f2 = thetafunc(x2,x,ii,tm,tc,trR,P,Id,nr,nc)
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
    plot(theta.val[ord],-LL.val[ord],xlab='theta',ylab='log-likelihood',pch=20)
    points(theta.val[ord],-LL.val[ord],type='l')
    abline(v=est.min,col='red')
  }
  
  return(est.min)
}
