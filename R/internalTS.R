# ---- roxygen documentation ----
#' @title Time Slice Function for fbtg
#'
#' @description
#'   Internal function for computing time slices.
#' @details
#'   Used in the \code{fbtg} function. 
#'   
#' @param ta time point of the slice.
#' @param t total duration of the segment.
#' @param cai accumulated cost surface from anchor A to location i.
#' @param cib accumulated cost sufrace from location i to anchor B. 
#' @param sigma parameter for modelling location uncertainty (s.d. of Guassian location error).
#' @param timefun function for modelling probabilities from time, input into Pfun.
#' @param c2 coefficient for modelling probabilities from time, input into Pfun.
#' @param clipPPS whether or not to clip to the PPS (default = TRUE).
#' 
#' @return
#'   This function returns a \code{RasterLayer} of the probabilities for time slice ta.
#'
#' @keywords internal 
# ---- End of roxygen documentation ----

internalTS <- function(ta,t,cai,cib,Tshort,sigma,timefun,c2,clipPPS){
  tb <- t - ta
  #Compute ST cone from point a at tk
  m <- matrix(c(ta+1,Inf,Inf),ncol=3,byrow=T)
  cai <- reclassify(cai,m)
  #Compute ST cone from point b at time s
  m <- matrix(c(tb+1,Inf,Inf),ncol=3,byrow=T)
  cib <- reclassify(cib,m)
  
  #Compute the binary PPS
  PPS <- cai+cib
  PPS[!is.infinite(PPS)] <- 1
  #Check that there exists some locations in the PPS
  ## NEED TO ADD IN AN ERROR CHECKING MECHANISM HERE ##
  if (is.infinite(cellStats(PPS,'min'))){
    next
  }
  
  #Compute time fraction relative to Expected time based on LCP and clip to PPS
  Tfa <- abs(ta/t * Tshort - cai)*PPS
  Tfb <- abs(tb/t * Tshort - cib)*PPS
  
  #Compute the total time deviation from the Expected time based on LCP
  Tfab <- Tfa+Tfb
  #Find minimal time (should be ~0 and should be on LCP, may be more than one cell)
  min <- cellStats(Tfab,'min',na.rm=TRUE)
  #Set minimal time location(s) = 0
  Tfab <- Tfab - min
  #Compute probabilities
  dt <- Tshort/t
  Pt <- internalPfun(Tfab,timefun,c2,dt)
  
  #####===== Locational UNCERTAINTY ANALYSIS===
  #uses same formulation as in Brownian bridge
  sigma2 <- ((1-(ta/t))^2 + (ta/t)^2)*sigma^2
  #uses cosine function to adjust locational uncertainty (highest near telemetry fixes, lowest away from them)
  #sigma2 <-(sigma^2/2)*(cos(ta/t*(2*pi))+1)   #division by 2 due to upward shift of cosine function
  w <- focalWeight(Pt,sqrt(sigma2),type='Gauss')
  if (dim(w)[1]> 1){ Pt <- focal(Pt,w,sum,pad=T,padValue=0) }   #just check to see if it makes a viable Gaussian function
  ####==============================
  
  #If clipPPS = TRUE: set to zero outside of PPS
  if (clipPPS){ 
    PPS[is.infinite(PPS)] <- 0
    Pt <- Pt*PPS 
  }
  
  #Normalize so that probability at any given time sums to 1
  Pt <- Pt / cellStats(Pt,stat='sum')
  
  return(getValues(Pt))
}