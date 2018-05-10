# ---- roxygen documentation ----
#' @title Randomised Shortest Path Utilization Distribution
#'
#' @description
#'   Calculate the utilization distribution of an animal based on randomised shortest paths between fixes.
#' @details
#' Still under development, but should work.
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param r a \code{RasterLayer} object describing the preference/affinity of the landscape. May be the result of a resource selection function, or other analyses. Note: Higher values should be associated higher affinity.
#' @param theta a global value for theta.
#' @param theta.col character string of the name of the column containing values of theta for each segment (only used if theta is not provided).
#'
#' @return
#'   This function returns a \code{RasterLayer} which can be used to estimate the UD of an animal.
#'
# @references
# @keywords 
#' @seealso ppa, fbtgUD
# @examples
#' @export
#
# ---- End of roxygen documentation ----

rspUD <- function(traj, r,theta,theta.col){
  #make a trajectory dataframe
  df <- ld(traj)
  n <- dim(df)[1]
  TT <- as.numeric(df$date[n]-df$date[1], units='secs')
  
  #Make *Symmetric* TransitionLayer from raster (e.g., avg cells)
  s1 <- function(x){x[1]}
  s2 <- function(x){x[2]}
  tr1 <- transition(r, s1, 8,symm=F)
  tr2 <- transition(r, s2, 8,symm=F)
  tr <- (tr1 + tr2)/2
  tr <- geoCorrection(tr,type='c',multpl=FALSE)
  
  #Define how theta is computed
  if (missing(theta)){
    if (missing(theta.col)){
      stop('theta or theta.col need to be provided.')
    } else {
      df$theta <- df[,theta.col]
    }
  } else {
    df$theta <- theta
  }
  
  #make an output raster
  Pi <- r*0
  #Progress Bar
  pb = txtProgressBar(min = 0, max = (n-1), initial = 0,style=3) 
  for (i in 1:(n-1)){   #Could parrallelize or apply?
    theta <- df$theta[i]
    dt <- df$dt[i]
    sp1 <- SpatialPoints(df[i,c('x','y')])
    sp2 <- SpatialPoints(df[i+1,c('x','y')])
    c1 <- cellFromXY(tr,sp1)
    c2 <- cellFromXY(tr,sp2)
    if (c1 == c2){
      #Start and end pixel is the same which means no movement. Need to adjust passage function.
      # Arbitrarily set end location to the next pixel over (check if edge)
      ########################################
      ind <- adjacent(r,c1,pairs=FALSE,id=TRUE)
      c2 <- ind[which.max(r[ind])]
      sp2 <- xyFromCell(r,c2,spatial=TRUE)
      sp2@proj4string <-sp1@proj4string
    } 
    ##########################################
    #Movement Occurs (at least from one cell to another)
    Pt <- passage(tr,sp1,sp2,theta=theta,totalNet='net')
    #Make sure Pi sum is proportional to the segment dt to account for unequally timed segments.
    Pt <- (Pt/cellStats(Pt,sum)) * dt/TT
    Pi <- Pi + Pt
    #update progress bar
    setTxtProgressBar(pb,i)
  }
  
  #Make 0's NA for easy plotting on output
  #Pi[Pi == 0] <- NA
  
  return(Pi)
}

