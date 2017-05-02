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
#' @param Tai accumulated costs from anchor A to location i.
#' @param Tib accumulated costs from location i to anchor B. 
#' @param tl \code{TransitionLayer} object used to compute Tai, Tib 
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

internalTS <- function(ta,t,Tai,Tib,tl,Tshort,sigma,timefun,c2,clipPPS){
  
  tb <- t - ta
  
  #Compute Delta T_i,t
  Tfab <- abs(Tai - ta/t*Tshort) + abs(Tib - tb/t*Tshort)
  
  #Find minimal time (should be ~0 and should be on LCP, may be more than one cell)
  min <- min(Tfab,na.rm=TRUE)
  Tfab <- Tfab - min                            #Set minimal time location(s) == 0
  
  Pt <- internalPfun(Tfab,timefun,c2)           #Compute probabilities
  Pt <- setValues(raster(tl),Pt)                #Make probabilities raster

  #####===== Locational UNCERTAINTY ANALYSIS===
  #uses same formulation as in Brownian bridge
  sigma2 <- ((1-(ta/t))^2 + (ta/t)^2)*sigma^2
  #uses cosine function to adjust locational uncertainty (highest near telemetry fixes, lowest away from them)
  #sigma2 <-(sigma^2/2)*(cos(ta/t*(2*pi))+1)   #division by 2 due to upward shift of cosine function
  w <- focalWeight(r,sqrt(sigma2),type='Gauss')
  if (dim(w)[1]> 1){ Pt <- focal(r,w,sum,pad=T,padValue=0) }   #just check to see if it makes a viable Gaussian function
  ####==============================
  
  #If clipPPS = TRUE: set to zero outside of PPS
  if (clipPPS){
    #Compute forward ST cone from point a
    m <- matrix(c(0,ta,1,ta,Inf,0),ncol=3,byrow=T)
    Tai <- reclassify(Tai,m)
    #Compute backward ST cone from point b
    m <- matrix(c(0,tb,1,tb,Inf,0),ncol=3,byrow=T)
    Tib <- reclassify(Tib,m)
    #Compute the binary PPS
    PPS <- Tai*Tib
    #Check that there exists some locations in the PPS
    ## NEED TO ADD IN AN ERROR CHECKING MECHANISM HERE ##
    #if (cellStats(PPS,'max') == 0 ){
    #  next
    #}
    r <- r*PPS 
  }
  
  #C1 Parameter = Normalize so that probability at any given time sums to 1.
  r <- r / cellStats(r,stat='sum')
  
  return(getValues(r))
}