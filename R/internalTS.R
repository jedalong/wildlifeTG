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
  print(1)
  Pt <- setValues(raster(tl),Pt)                #Make probabilities raster

  #####===== Locational UNCERTAINTY ANALYSIS===
  #uses same formulation as in Brownian bridge
  sigma2 <- ((1-(ta/t))^2 + (ta/t)^2)*sigma^2
  w <- focalWeight(Pt,sqrt(sigma2),type='Gauss')
  if (dim(w)[1]> 1){ Pt <- focal(Pt,w,sum,pad=T,padValue=0) }   #just check to see if it makes a viable Gaussian function
  ####==============================
  
  
  #If clipPPS = TRUE: set to Inf outside of PPS, expand due to location uncertainty if sigma > 0.
  if (clipPPS){
    ppsa <- Tai*0
    ppsb <- Tib*0
    ppsa[which(Tai <=ta)] <- 1
    ppsb[which(Tib <=tb)] <- 1
    PPS <- Pt*0
    PPS <- setValues(PPS,ppsa*ppsb)
    #Include location uncertainty in clipping of PPS
    if (dim(w)[1]>1){
      kk <- trunc(dim(w)[1]/2)*mean(res(PPS))
      PPS[Which(PPS == 0)] <- NA
      PPS <- buffer(PPS,width=kk)
      PPS[is.na(PPS)] <- 0
    }
    Pt <- Pt*PPS
  }

  #C1 Parameter = Normalize so that probability at any given time sums to 1.
  Pt <- Pt / cellStats(Pt,stat='sum')
  
  return(getValues(Pt))
}