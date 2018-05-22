# ---- roxygen documentation ----
#' @title Dynamic PPA Measure of Animal Space Use
#'
#' @description
#'   The function \code{dynppa} computes the (dynamic) PPA measure of the accessibility space of an animal. 
#'   The PPA method can be thought as an alternative view on the home range; one that explicitly considers the
#'   spatial and temporal constraints on movement given known telemetry fixes, and a (dynamic) measure of maximum
#'   mobility - termed Vmax. The PPA method incorporates dynamic behaviour into the calculation of the vmax parameter 
#'   used to delineate the original version of the PPA method, but the original method is still an option here.
#' @details
#'   The function \code{dyn.ppa} represents an extension to an existing PPA method (Long and Nelson, 2012). 
#'   Dynamic calculation of the PPA method improves upon the original version by flexibly modelling the vmax 
#'   parameter according to wildlife behaviour. See the function \code{dyn.vmax} for more information on how 
#'   to incorporate dynamic behaviour into the vmax parameter estimation. 
#'   
#' @param traj an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param tol parameter used to filter out those segments where the time between fixes is overly 
#'    large (often due to irregular sampling or missing fixes); which leads to an overestimation of the 
#'    activity space via the PPA method. Default is the maximum sampling interval from \code{traj}.
#' @param dissolve (logical) whether or not to dissolve output elliplse polygons to create a single
#'    output polygon, or keep the individual segment PPA ellipses. Default = TRUE.
#' @param proj4string a string object containing the projection information to be passed included in the output 
#'    \code{SpatialPolygonsDataFrame} object. For more information see the \code{CRS-class} in the packages
#'    \code{sp} and \code{rgdal}. Default is \code{NA}.
#' @param ePoints number of vertices used to construct each PPA ellipse. More points will necessarily provide
#'    a more detailed ellipse shape, but will slow computation; default is 360.
#' @param ... additional parameters to be passed to the function \code{dynvmax}. For example, should include
#'    \code{method} and/or \code{dynamic} parameters, see the documentation for \code{dynvmax} for more detailed
#'    information on what to include here.
#'    
#' @return
#'   This function returns a \code{SpatialPolygonsDataFrame} representing the dynamic PPA measure of the accessibility
#'   space of an individual animal.
#'
#' @references
#'   Long, JA, Nelson, TA. (2012) Time geography and wildlife home range delineation. \emph{Journal of Wildlife
#'   Management}, 76(2):407-413.\cr \cr
#'   Long, JA, Nelson, TA. (2014) Home range and habitat analysis using dynamic time geography. \emph{Journal of 
#'   Wildlife Management}. Accepted: 2014-12-03.\cr
#'   
# @keywords 
#' @seealso dynvmax
#' @examples 
#' data(m3)
#' ppa1 <- dynppa(m3,method='vanderWatt')
#' ppa2 <- dynppa(m3,method='vanderWatt',dynamic='focal')
#' plot(ppa1)
#' plot(ppa2,add=T,border='red')
#' 
#' @export
#
# ---- End of roxygen documentation ----

dynppa <- function(traj, tol=max(ld(traj)$dt,na.rm=TRUE),dissolve=TRUE,proj4string=CRS(as.character(NA)),ePoints=360, ...){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}

  
  #Append the local Vmax values to the trajectory using the LocalVmax function
  traj <- dynvmax(traj, ...)
  
  #convert to dataframe
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #loop through trajectory dataset and calculate PPA
  polyList <- vector('list',length=(n-1))
  for (i in 1:(n-1)){ #replace with function and use lapply
    #check to see if the local Vmax value exists
    if (!is.na(trDF$dynVmax[i])){      
      #check to see if the time interval is below the tolerance level
      if (trDF$dt[i] <= tol){
        cpX <- (trDF$x[i] + trDF$x[i+1])/2
        cpY <- (trDF$y[i] + trDF$y[i+1])/2
        #check to see if the angle is NA, if it is make it zero (no movement)
        if (is.na(trDF$abs.angle[i]) == TRUE) {
          thetaRot <- 0
          } else {thetaRot <- trDF$abs.angle[i]}
        
        c <- trDF$dist[i]/2
        a <- (trDF$dt[i]*trDF$dynVmax[i])/2
        b <- sqrt((a^2)-(c^2))

        #Compute the ellipse
        polyList[[i]] <- ppaEllipse(cpX,cpY,a,b,thetaRot,ePoints)
        }
      }
    }
  
  #---------------------------------------------
  # Create spatial polygons from the output
  #---------------------------------------------
  ind.poly <- which(lapply(polyList, is.null) == FALSE)
  polyList <- polyList[ind.poly]
  polyList <- lapply(polyList,list)
  tempPoly <- mapply(Polygons, polyList, ind.poly)
  data <- data.frame(ind=ind.poly, date=trDF$date[ind.poly])
  spPoly <- SpatialPolygonsDataFrame(SpatialPolygons(tempPoly,proj4string=proj4string),data,match.ID=FALSE)
  
  #"union" multiple spatial polygons stored in single sp object
  #(dissolve in GIS) 
  if (dissolve == TRUE){
    spUnion <- gUnaryUnion(spPoly)
    spUnion <- SpatialPolygonsDataFrame(spUnion,data=data.frame(id=1:length(spUnion)))
    return(spUnion)
  } else {
    return(spPoly)
  }
}
  #End of Function
#===============================================================================



  



