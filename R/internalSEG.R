# ---- roxygen documentation ----
#' @title Segment Function for fbtg
#'
#' @description
#'   Internal function for computing Pi for segments.
#' @details
#'   Used in the \code{fbtgUD/TS} functions. 
#'   
#' @param j index of the segment.
#' @param df trajectory dataframe (from \code{ld(traj)})
#' @param tl \code{TransitionLayer} modelling movement costs. 
#' @param k number of steps for use in \code{internalTS}.
#' @param sigma parameter for modelling location uncertainty (s.d. of Guassian location error).
#' @param timefun function for modelling probabilities from time, input into Pfun.
#' @param c2 coefficient for modelling probabilities from time, input into Pfun.
#' @param clipPPS whether or not to clip to the PPS (default = TRUE).
#' 
#' @return
#'   This function returns a \code{RasterLayer} of the probabilities for time slice ta.
#' @keywords internal 
# ---- End of roxygen documentation ----

#
internalSEG <- function(j, df, tl, k, sigma, timefun, c2, clipPPS){
  
  A <- c(df[j,1],df[j,2])         #origin point
  B <- c(df[j+1,1],df[j+1,2])     #destination point
  xy <- rbind(A,B)
  t <- df$dt[j]
  
  #Tai <- JaccCost(tl,A,'out')
  #Tib <- JaccCost(tl,B,'in')
  Ai <- cellFromXY(raster(tl),A)
  Bi <- cellFromXY(raster(tl),B)
  
  #Compute Accumulated cost (i.e, time) for location A and B
  #This is the value input for various functions
  tm <- transitionMatrix(tl)
  gr <- graph.adjacency(tm, mode="directed", weighted=TRUE)
  E(gr)$weight <- 1/E(gr)$weight
  Tai <- distances(gr,v=Ai,to=V(gr),mode='out')
  Tib <- distances(gr,v=Bi,to=V(gr),mode='in')
  
  #compute the cost of the shortest path
  Tshort <- costDistance(tl,A,B)[1]
  #time slices used to estimate prism
  tk <- seq(0,t,length.out=k+2)[2:(k+1)]
  
  #empty raster to populate
  P <- raster(tl)*0
  P.v <- getValues(P)
  #Compute the values for each Time Slice using internalTS function
  PP <- vapply(tk,internalTS,FUN.VALUE=P.v,t=t,Tai=Tai,Tib=Tib,tl=tl,Tshort=Tshort,sigma=sigma,timefun=timefun,c2=c2,clipPPS=clipPPS)
  P <- setValues(P, apply(PP,1,sum))
  
  
  #make an origin destination raster that is 1 for origin/destination pixels or 0 otherwise, and adjust for uncertainty
  xx <- rasterize(xy,raster(tl))*0+1 
  xx[is.na(xx)] <- 0
  if (sigma > 0){
    w <- focalWeight(xx,sigma,type='Gauss')
    if (dim(w)[1]>1){ xx <- focal(xx,w,sum,pad=T,padValue=0) }
  }
  #if (clipPPS){ xx <- xx*(rasterize(xy,raster(tl))*0+1) }
  
  #Numerical integration by trapezoidal rule approximation cumulative probability of prism for pixel.
  Pi <- t/(k+1) * (P + 0.5*xx)  #xx is anchors (multiply by 0.5 because each point is used twice in trajectory - except for boundary points, but ok)
  
  #Make sure Pi sums to the segment time budget to account for unequally timed segments  ### check necessary?
  Pi <- (Pi/cellStats(Pi,sum)) * t
  
  return(Pi)
}