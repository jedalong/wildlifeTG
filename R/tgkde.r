# ---- roxygen documentation ----
#' @title Time geographic kernel density estimation
#'
#' @description
#'   Compute the time geographic kernel density estimate of an animals utilization distribution following the methods described in Downs et al. (2011). 
#' @details
#'   The function \code{tgkde} can be used to delineate an animals home range using the time geographic kernel density estimation procedure described by Downs et al. (2011). Specifically, it modifies the shape of the traditional kernel to consider the time geographic limits of movement opportunity - termed the geoellipse by Downs, which is analagous to the potential path area concept described by Long & Nelson (2012,2015). Several basic functions -- including inverse distance, inverse distance squared, exponential, and normal -- can be used to quantify movement probabilities within the geoellipse. The output is then the utilization distribution of the animal, confined to the accessibility space defined by the potential path area home range as in Long & Nelson (2012,2015).
#'   The function \code{volras} can be used to extract volume contours (e.g., 95% volume contour) which are commonly used to delineate home ranges based on the utilization distribution.
#'
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first object. Note this object must be a \code{type II ltraj} object. For more information on objects of this type see \code{ help(ltraj)}.
#' @param disfun one of { \code{'inv', 'inv2', 'exp', 'norm'} } representing the shape of the distance decay function used to model movement probabilities within the geoellipse (see Downs et al. (2011) for more information). 
#' @param grid spatial resolution (pixel size) of output utilization raster in appropriate units. Default is chosen based on the x and y range of the input telemetry data. Alternatively, a \code{RasterLayer} can be passed in upon which the UD is computed.
#' @param ... additional parameters to be passed to the function \code{dynvmax}. For example, should include options for \code{dynamic} and \code{method}; see the documentation for \code{dynvmax} for more detailed information on what to include here.
#'    
#' @return
#'   This function returns a \code{RasterLayer} representing the utilization distribution of the animal
#'
#' @references
#' Downs, J.A., Horner, M.W., Tucker, A.D. (2011) Time-geographic density estimation for home range analysis. Annals of GIS. 17(3): 163-171.
# @keywords 
#' @seealso dynvmax, dynppa, volras
# @examples
#' 
#' @export
#
# ---- End of roxygen documentation ----

tgkde <- function(traj,disfun='inv',grid=NA,...){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  #Append the local Vmax values to the trajectory using the dynvmax function
  traj <- dynvmax(traj, ...)
  
  #convert to dataframe
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #If grid is rasterLayer
  if (class(grid) == 'RasterLayer'){
    xy <- data.frame(coordinates(grid))
  } else {
    #-------------------------------------------
    # Create the output grid
    #-------------------------------------------
    x.range <- max(trDF$x) - min(trDF$x)
    y.range <- max(trDF$y) - min(trDF$y)
    
    if (is.na(grid)){
      grid <- 0.01*min(c(x.range,y.range)) 
    }
    xx <- seq(min(trDF$x)-0.1*x.range,max(trDF$x)+0.1*x.range,by=grid)
    yy <- seq(min(trDF$y)-0.1*y.range,max(trDF$y)+0.1*y.range,by=grid)
    xy <- expand.grid(x=xx,y=yy)
    #-------------------------------------------
  }
  xy$z <- 0

  #loop through Traj object and compute the tgkde estimate
  x <- xy$x
  y <- xy$y
  for (i in 1:(n-1)){
    sx <- trDF$x[i]
    sy <- trDF$y[i]
    ex <- trDF$x[i+1]
    ey <- trDF$y[i+1]
    dt <- trDF$dt[i]
    vmax <- trDF$dynVmax[i]
    if (is.na(vmax)){next}
    dd <- sqrt((sx-x)^2 + (sy-y)^2) + sqrt((x-ex)^2 + (y-ey)^2)
    dp <- sqrt((sx-ex)^2 + (sy-ey)^2)
    ind <- which(dd > dt*vmax)
    #insert function here perhaps using switch?
    g <- switch(disfun,
                inv = dp/dd,              #inverse distance
                inv2 = (dp/dd)^2,         #inverse distance ^2        
                exp = exp(-dd/dp),        #exponential function 
                norm = exp(-(dd/dp)^2),   #normal function
                stop(paste('The distance decay function',disfun,'does not exist.'))
    )
      
    g[ind] <- 0
    #In Theory, need to normalize here, because each ellipse should sum to dt
    g <- dt*g/sum(g)
    
    xy$z <- xy$z + g
  }
  
  #---- This all seems rather arbitrary? -----
  ##Normalize following the eqn. in Downs et al. (2011)
  ##get the 'average' vmax
  #vmax. <- mean(trDF$dynVmax,na.rm=TRUE)
  ##get the 'overall' time difference
  #DT <- as.numeric(difftime(trDF$date[n], trDF$date[1], units='secs'))
  ##Normalize the values
  #xy$z <- xy$z * (1 / ((n-1)*(DT*vmax.)^2))
  
  # Downs et al. 2011 method does not sum to 1... 
  # based on normalization of each ellipse, just dividing by n-1 should work
  # xy$z <- xy$z / (n-1)
  
  #----------------------------
  #  Format output to Raster
  #----------------------------
  coordinates(xy) = ~x+y
  gridded(xy) <- TRUE
  ras <- raster(xy)
  #----------------------------
  return(ras)
}
#=============DEMO==========================
# traj <- simm.crw(1:500)
# trDF <- ld(traj)
# outras <- tgkde(traj,grid=0.25)
# #Draw the plot
# library(fields)
# image(outras,xlab="X",ylab="Y",asp=1)
# lines(trDF$x,trDF$y)
#============================================