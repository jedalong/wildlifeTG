# ---- roxygen documentation ----
#' @title Randomised Shortest Path Utilization Distribution
#'
#' @description
#'   Calculate the utilization distribution of an animal based on randomised shortest paths between fixes.
#' @details
#' A randomised shortest path model is fit between every pair of fixes. The randomised shortest path model is derived from the \code{passage} function in the package \code{gdistance}. It uses the net number of packages (see \code{?passage}) to estimate the probability and then scales the values appropriately. An input \code{rasterLayer} object is required which defines the ability for movement through the landscape, which might typically be derived from a resource selection function, or be related to known barriers on the landscape. It requires only a single parameter \code{theta} which can be estimated from the data using the function \code{esttheta}. 
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param r a \code{RasterLayer} object describing the preference/affinity of the landscape. May be the result of a resource selection function, or other analyses. Note: Higher values should be associated higher affinity.
#' @param theta a global value for theta.
#' @param theta.col character string of the name of the column containing values of theta for each segment (only used if theta is not provided).
#' @param timescale (logical) whether or not to scale each segment so the sum of the output surface is equal to the duration of the segment (in seconds). Default = FALSE.
#'
#' @return
#'   This function returns a \code{RasterLayer} which can be used to estimate the UD of an animal.
#'
#' @references Long, J.A. Estimating wildlife utilization distributions using randomized shortest paths. (in Preparation)
# @keywords 
#' @seealso esttheta, volras
# @examples
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#
# ---- End of roxygen documentation ----

rspUD <- function(traj,r,theta,theta.col,timescale=FALSE){
  #make a trajectory dataframe
  df <- ld(traj)
  n <- dim(df)[1]
  TT <- sum(df$dt,na.rm=T)
  
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
    df$theta <- theta*(mean(df$dt,na.rm=T)/df$dt)^2
  }
  
  #make an output raster
  Pi <- r*0
  #Progress Bar
  pb = txtProgressBar(min = 0, max = (n-1), initial = 0,style=3) 
  for (i in 1:(n-1)){   #Could parrallelize or apply?
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
    Pt <- passage(tr,sp1,sp2,theta=df$theta[i],totalNet='net')
    Pt[c(c1,c2)] <- 0.5   #scale down prob 1 at end points to account for double counting fixes.
    Pt <- (Pt/cellStats(Pt,sum)) 
    if (timescale){
      #Make sure Pi sum is proportional to the segment dt to account for unequally timed segments.
      Pt <- Pt*(df$dt[i]/TT)
    }
    Pi <- Pi + Pt
    #update progress bar
    setTxtProgressBar(pb,i)
  }
  
  #Make 0's NA for easy plotting on output
  #Pi[Pi == 0] <- NA
  
  return(Pi)
}

