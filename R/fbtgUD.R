# ---- roxygen documentation ----
#' @title ***UNDER DEVELOPMENT*** Field-based Time Geography Utilization Distribution
#'
#' @description
#'   Calculate the utilization distribution of an animal using the field-based time geographic model. ***** THIS FUNCTION IS STILL IN DEVELOPMENT *****
#' @details
#'   Calculates the field-based time geography utilization distribution (UD) for an animal. Field-based time geography is based on an underlying resistance surface, which constratins potential movement by the animal between two location fixes. This model is applied recursively over an entire trajectory in order to compute a UD for an animal. The UD inherently considers the movement limitations described by the underlying resistance surface, and thus is considered a landscape-based model for a UD. The landscape-based approach deviates from current models building upon random walks and diffusion processes. \cr
#'   The model requires that the resistance surface be directly related to an animals speed of passing through that environment, please see the vignette for more details as to how this might be constructed.
#'   The timefun parameter is used to choose the model for converting time into a probability based on the described function. The \code{c2} parameter can be used to tune these functions based on some fine-scale movement data if available. 
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param surf a \code{TransitionLayer} object
#' @param timefun method for convertingtime into probability; one of:\cr 
#' -- \code{'inverse'} \eqn{= \frac{1}{c_2 * t}},\cr
#' -- \code{'inverse2'} \eqn{= \frac{1}{c_2 * t^2}},\cr
#' -- \code{'exp'} \eqn{= \exp{-c_2 * t}},\cr
#' -- \code{'norm'} \eqn{= \exp{-c_2 * t^2}},\cr
#' -- \code{'rootexp'} \eqn{= \exp{-c_2 * \sqrt{t}}},\cr
#' -- \code{'pareto'} \eqn{= \exp{-c_2 * \log{t}}},\cr
#' -- \code{'lognorm'}\eqn{= \exp{-c_2 * \log{t^2}}}.\cr
#' @param c2 Parameter input into \code{timefun} default=1.
#' @param k number of time slices between fixes upon which to estimate the UD, default=10. 
#' @param spd.fac speed adjustment factor, can be used to globally adjust the resistance surface based on some factor such as season, time of day, beahvioural mode, etc. Takes as input as the name of the column storing the adjustment in the \code{InfoLocs} part of the \code{ltraj} object.
#' @param d.min minimum distance, below which the segment is removed from analysis. Can be used to focus the UD on only longer movement segments. Default is the pixel size of \code{surf}.
#' @param dt.max maximum time, above which the segment is removed from analysis. Can be used to remove segments with missing data from analysis, as segments with a very long time between fixes are problematic in time geographic analysis.


#' @return
#'   This function returns a \code{RasterLayer} which can be used to estimate the UD of an animal.
#'
# @references
# @keywords 
#' @seealso ppa
# @examples
#' @export
#
# ---- End of roxygen documentation ----

fbtgUD <- function(traj,surf,timefun='inverse',c2=1,k=10,spd.fac=FALSE,d.min='default',dt.max=FALSE){
  #get minimum movement distance that potentially starts and ends in the same cell.
  if (d.min == 'default'){
    res <- res(as(surf,'RasterLayer'))
    d.min <- sqrt(res[1]^2+res[2]^2)
  }
  
  df <- ld(traj)
  n <- dim(df)[1]
  r0 <- raster(surf)*0
  n.m <- 0                     #count of actual 'movement' segments
  
  #Get speed adjustment column, ideally something like output from dynvmax function.
  df$spd.fac <- 1
  if (spd.fac != FALSE){df$spd.fac <- df[,get('spd.fac')]}
  
  #create dt vector outside of loop
  dt <- df$dt
  #get time threshold to deal with long durations without data
  if (dt.max == FALSE){dt.max <- Inf}
  #if duation between fixes is above dt.max make dt = 0 and then it will be an 'id.fail'
  dt[which(dt > dt.max)] <- 0
  
  #IDs of non-movement segments  
  id.still <- NULL
  #IDs of fails because cost surface incorrectly modelled (i.e., not enough time) or due to long gaps.
  id.fail <- NULL
  
  #can we vectorize this?
  for (j in 1:(n-1)){
    
    #Right now just don't consider any segments where movement isn't above some baseline threshold (defaults to pixel size)
    #Will need to fix this eventually.
    if (df$dist[j] <= d.min){
      id.still <- c(id.still,j)
      next
    }
    #If it is a movement segment then do the analysis
    A <- c(df[j,1],df[j,2])         #origin point
    B <- c(df[j+1,1],df[j+1,2])     #destination point
    xy <- rbind(A,B)
    
    #make an origin destination raster that is 1 for origin/destination pixels or 0 otherwise
    #Can we adjust this for uncertainty?
    xx <- rasterize(xy,raster(surf))*0+1 
    xx[is.na(xx)] <- 0
    
    
    t <- dt[j]
    surf.j <- surf*df$spd.fac[j]  
    #Compute Accumulated cost (i.e, time) for location A and B
    cai <- JaccCost(surf.j,A,'out')
    cib <- JaccCost(surf.j,B,'in')
    ts <- seq(0,t,length.out=k+2)[2:(k+1)]
    P <- raster(surf)*0
    
    #can we vectorize this
    for (s in 1:k){
      chk <- 1
      ta <- ts[s]
      tb <- t - ts[s]
      
      #while loop to check if the time budget is feasible
      while (chk > 0){
        chk <- 0
        #Compute ST cone from point a at time s
        m <- matrix(c(ta+1,Inf,Inf),ncol=3,byrow=T)
        acs <- reclassify(cai,m)
        #Compute ST cone from point b at time s
        m <- matrix(c(tb+1,Inf,Inf),ncol=3,byrow=T)
        bcs <- reclassify(cib,m)
        #Compute PPS at time s
        PPS <- acs + bcs
        if (cellStats(PPS,'min')==Inf){
          ta <- ta*1.1
          tb <- tb*1.1
          chk <- 1
          id.fail <- c(id.fail,j)
        }
      }
      #Convert infinity to NA
      PPS[is.infinite(PPS)] <- NA
      
      #Inverse time to get probability and control for minimum amount of time to destination.
      Pt <- switch(timefun,
                   inverse = 1 / (c2*PPS),
                   inverse2 = 1 / (c2*(PPS^2)),
                   exp = exp(-c2*PPS),
                   norm = exp(-c2*(PPS^2)),
                   rootexp = exp(-c2*(PPS^0.5)),
                   pareto = exp(-c2*log(PPS)),
                   lognorm = exp(-c2*(log(PPS)^2)),
                   stop(paste('The time function',timefun,'does not exist.'))
      )
      
      #Convert NA to 0
      Pt[is.na(Pt)] <- 0

      
      #Normalize so that probability at any given time sums to 1
      Pt <- Pt / cellStats(Pt,stat='sum')
      P <- P + Pt
    }
    #Numerical integration by trapezoidal rule approximation cumulative probability of prism for pixel.
    Pi <- t/(k+1) * (P + 0.5*xx)  #xx is 1 for anchors, 0 otherwise. This may be biased by GPS error and pixel resolutions. Revisit.
    r0 <- r0 + Pi
    n.m <- n.m + 1
  }
  #sum raster, and compute normalized version
  #rr <- calc(rstack,sum)  
  #values(r0) <- values(r0)/n.m  
  #Make 0's NA for easy plotting on output
  r0[r0 == 0] <- NA
  if (!is.null(id.still)){print(paste('The no. of IDs where movement was below d.min:',length(id.still)))}
  if (!is.null(id.fail)){print(paste('The no. of IDs where the cost surface was too slow:',length(id.fail)))}
  return(r0)
}


