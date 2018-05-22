# ---- roxygen documentation ----
#' @title Slow movement areas
#'
#' @description
#'   The function \code{sma} computes the areas representing slow movement areas as described in the paper Nelson et al. (2014). Slow movement areas represent areas of sustained or intense habitat use, related to slow movement behaviours. Slow movement areas are defined by counting consecutive fixes within time geographic ellipses, and represented spatially as spatial polygons that are the union of included telemetry fixes within the slow movement area.
#'   
#' @details
#'   The function \code{sma} can be used to map slow movement areas identifiable from wildlife telemetry data.of potential interaction between two animals. Slow movement areas can be ranked, according to their importance, which equates to consecutive time spent in an area. That is, the first slow movement area will be the area where the animal stayed the longest, and so on. Thus, slow movement areas can be useful for identifying where encamped behaviour or intensively exploited habitat on the landscape.
#'
#' @param traj an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the animal. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param sma.keep an integer value indicating the number of slow movement areas to delineate, default is 1.
#' @param sma.tol a value <= 1 indicating used when sma.keep > 1 to define how much overlap is allowed between SMA's, if sma.tol=1 no overlap is allowed, if sma.tol=0, any and all overlap is allowed. Typically something in between is most useful. Defaults to 1.
#' @param tol (optional) parameter used to filter out those segments where the time between fixes is overly 
#'    large (often due to missing fixes); which leads to an overestimation of the 
#'    activity space via the PPA method. Default is the maximum sampling interval from \code{traj1}.
#' @param proj4string a string object containing the projection information to be passed included in the output 
#'    \code{SpatialPolygonsDataFrame} object. For more information see the \code{CRS-class} in the packages
#'    \code{sp} and \code{rgdal}. Default is \code{NA}.
#' @param ePoints number of vertices used to construct each PPA ellipse. More points will necessarily provide
#'    a more detailed ellipse shape, but will slow computation; default is 360.
#' @param ... additional parameters to be passed to the function \code{dynvmax}. For example, should include options for
#'    \code{dynamic} and \code{method}; see the documentation for \code{dynvmax} for more detailed information 
#'    on what to include here.
#'    
#' @return
#'   This function returns a \code{SpatialPolygonsDataFrame} representing the joint accessibility space between
#'   the two animals.
#'
#' @references
#' Nelson, T.A., Long, J.A., Laberee, K., Stewart, B.P. (2015) A time geographic approach for delineating areas of sustained wildlife use. Annals of GIS. 21(1): 81-90.
# @keywords 
#' @seealso dynvmax, dynppa
#' @examples
#' data(m3)
#' sm1 <- sma(m3,method='vanderWatt')
#' sm2 <- sma(m3,sma.keep=2,method='vanderWatt')
#' 
#' @export
#
# ---- End of roxygen documentation ----


sma <- function(traj, 
                sma.keep=1, 
                sma.tol=1, 
                tol=max(ld(traj)$dt,na.rm=TRUE), 
                proj4string=CRS(as.character(NA)), 
                ePoints=360, 
                ...){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  
  #Append the local Vmax values to the trajectory using the LocalVmax function
  traj <- dynvmax(traj, ...)
  
  #convert to dataframe
  trDF <- ld(traj)
  n <- dim(trDF)[1]

  #loop through trajectory dataset and calculate PPA
  polyList <- vector('list',length=(n-1))
  smaList <- rep(0,(n-1))
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
        
        #If the ellipse is not NULL count successive intersecting points...
        if (!is.null(polyList[[i]])){
          #Create an SP object
          tempPoly <- SpatialPolygons(list(Polygons(polyList[i],ID=as.character(i))))
          #Count the number of consecutive telemetry points in the PPA ellipse
          pts <- SpatialPointsDataFrame(trDF[i:n,1:2],data=data.frame(ID=i:n))
          int <- gIntersects(tempPoly,pts,byid=T)
          bb <- pts$ID[int]
          cc <- which(diff(bb) != 1)[1]
          smaList[i] <- cc
        }  
      }
    }
  }

  #df is the data.frame associated with information of each sma. It sorts here by SMA.value
  # then it filters out overlap based on the parameter sma.tol {0,1}. Here we are trying to 
  # identify 'new' SMA's by eliminating those that overlap temporally with previous SMA's. Note
  # that spatial overlap is OK (and therefore not checked). The values of new.sma represent 'new'
  # SMA's identified with 0 indicating complete overlap with a previously defined SMA and 1 
  # indicating no overlap.
  
  df <- data.frame(sort(smaList,decreasing=TRUE,index.return=TRUE))
  names(df) <- c('SMA.value','i.start')
  df$i.end <- df$i.start + df$SMA.value-1
  df$t.start <- trDF$date[df$i.start]
  df$t.end <- trDF$date[df$i.end]
  df$new.sma <- 0
  df$new.sma[1] <- 1
  n <- dim(df)[1]
  ind <- NULL
  for (i in 2:n){
    ind <- unique(c(ind,df$i.start[i-1]:df$i.end[i-1]))
    ind. <- df$i.start[i]:df$i.end[i]
    n1 <- length(ind.) - length(which(ind. %in% ind))
    df$new.sma[i] <- n1/length(ind.)
  }
  df <- subset(df, df$new.sma >= sma.tol)
  
  #identify which indices to keep
  i.start <- df$i.start[1:sma.keep]
  i.end <- df$i.end[1:sma.keep]
  
  ###FAILING BELOW

  ###This section compiles the polygon ellipses associated with each SMA.
  for (i in 1:sma.keep){
    #extract the polys associated with each SMA
    # Note: the i.end - 1 is because ellipses are defined for two consecutive points.
    ind.poly <- i.start[i]:(i.end[i]-1)
    pList <- polyList[ind.poly]
    #Remove NULL polygons
    ind.null<- which(lapply(pList, is.null) == FALSE)
    pList <- pList[ind.null]
    ind.poly <- ind.poly[ind.null]
    pList <- lapply(pList,list)
    tempPoly <- mapply(Polygons, pList, ind.poly)
    poly <- SpatialPolygons(tempPoly,proj4string=proj4string)
    u.poly <- gUnaryUnion(poly,id=as.character(i.start[i]))
    if (i == 1){
      sma.p <- u.poly
    } else {
      sma.p <- rbind(sma.p,u.poly)
    }
  }

  #subset output attribute data
  df.out <- df[1:sma.keep,1:5]
  #fix projection issues
  #slot(sma.p,'proj4string') <- proj4string
  #output SpatialPolygonsDataFrame
  sma.poly <- SpatialPolygonsDataFrame(sma.p,data=df.out,match.ID=FALSE)
  
  return(sma.poly)
}

  