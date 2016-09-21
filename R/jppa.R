# ---- roxygen documentation ----
#' @title Joint Potential Path Area  of Two Animals
#'
#' @description
#'   The function \code{jppa} computes the joint accessibility space between two animals. It can be used
#'   to map (as a spatial polygon) the area that could have been jointly accessed by two individual animals 
#'   in space ant time. The jPPA represents a spatial measure of spatial-temporal interaction.
#' @details
#'   The function \code{jppa} can be used to map areas of potential interaction between two animals.
#'   Specifically, this represents a measure of spatial overlap that also considers the temporal sequencing 
#'   of telemetry points. In this respect it improves significantly over static measures of home range overlap, 
#'   often used to measure static interaction, and can be considered as a spatial measure of dynamic interaction.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param t.int (optional) time parameter (in seconds) used to determine the frequency of time slices 
#'    used to delineate the joint activity space. Default is 1/10th of the mode of the temporal sampling
#'    interval from \code{traj1}. Smaller values for \code{t.int} will result in smoother output polygons.
#' @param tol (optional) parameter used to filter out those segments where the time between fixes is overly 
#'    large (often due to irregular sampling or missing fixes); which leads to an overestimation of the 
#'    activity space via the PPA method. Default is the maximum sampling interval from \code{traj1}.
#' @param dissolve logical parameter indicating whether (\code{=TRUE}; the default) or not (\code{=FALSE})
#'    to return a spatially dissolved polygon of the joint activity space.
#' @param ePoints number of vertices used to construct each PPA ellipse. More points will necessarily provide
#'    a more detailed ellipse shape, but will slow computation; default is 360.
#' @param proj4string a string object containing the projection information to be passed included in the output 
#'    \code{SpatialPolygonsDataFrame} object. For more information see the \code{CRS-class} in the packages
#'    \code{sp} and \code{rgdal}. Default is \code{NA}.
#' @param ... additional parameters to be passed to the function \code{dynvmax}. For example, should include options for
#'    \code{dynamic} and \code{method}; see the documentation for \code{dynvmax} for more detailed information 
#'    on what to include here.
#'    
#' @return
#'   This function returns a \code{SpatialPolygonsDataFrame} representing the joint accessibility space between
#'   the two animals.
#'
#' @references
#' Long, J.A., Webb, S.L., Nelson, T.A., Gee, K. (2015) Mapping areas of spatial-temporal overlap from wildlife telemetry data. Movement Ecology. 3:38.
# @keywords 
#' @seealso dynvmax, dynppa
# @examples
#' 
#' @export
#
# ---- End of roxygen documentation ----
##Function for computing interaction range
jppa <- function(traj1,
                      traj2,
                      t.int=0.1*as.numeric(names(sort(-table(ld(traj1)$dt)))[1]),
                      tol=max(ld(traj1)$dt,na.rm=T), 
                      dissolve = TRUE, 
                      proj4string=CRS(as.character(NA)),
                      ePoints=360,           
                      ...){
  
  #EXTRA FUNCTIONS
  #====================================================================================  
  # Function get.anchors, gets two trajectory points surrounding a time point.
  get.anchors <- function(traj,indo){
    pF <- traj[indo,c('x','y','date','dynVmax','dt','dist','abs.angle')]
    pP <- traj[indo+1,c('x','y','date','dynVmax','dt','dist','abs.angle')]
    ret <- rbind(pF,pP)
    ret$ind <- c(indo,(indo+1))
    return(ret)
  }
  get.anchor <- function(traj,t.slice){
    indF <- max(which(difftime(traj$date,t.slice) <= 0))
    return(indF)
  }
  
  #function for identifying the prism slice instersection polygon
  prism.slice.int <- function(t,a1,a2,tol,ePoints){
    
    p1 <- prism.slice(a1,t,tol,ePoints)
    p2 <- prism.slice(a2,t,tol,ePoints)
    #prism slice intersection
    if (!is.null(p1) & !is.null(p2)) {
      pInt <- gIntersection(p1,p2,id=as.character(t))
      if (!is.null(pInt)) {
        jppa <- slot(pInt,'polygons')
        return(jppa)
      } else {return(NULL)}
    } else {return(NULL)}
  }
  
  # Function prism.slice computes the ST prism polygon slice for a given time point.
  prism.slice <- function(anchors,t.slice,tol,ePoints){
    if (anchors[1,'dt'] > tol || is.na(anchors[1,'dynVmax']) || is.na(anchors[1,'dt'])){return(NULL)}
    else{
      #vmax always associated with the forward point due to how it is defined for segment.
      vmax <- anchors[1,'dynVmax']
      #the number of points around the circle
      theta <- seq(0,2*pi,length.out=ePoints)
      #calculate the radius of the Future Cone
      tF <- difftime(t.slice,anchors[1,3],units="secs")
      rF <- vmax*tF  #vmax is anchors[1,4]
      #calculate the radius of the Past Cone
      tP <- difftime(anchors[2,3],t.slice,units="secs")
      rP <- vmax*tP
      #get x,y coords of Future Cone circle
      x1 <- anchors[1,1] + rF*cos(theta)
      y1 <- anchors[1,2] + rF*sin(theta)
      #get x,y coords of Past Cone circle
      x2 <- anchors[2,1] + rP*cos(theta)
      y2 <- anchors[2,2] + rP*sin(theta)
      #create spatial polygons and intersect them
      c1 <- Polygon(rbind(cbind(x1,y1),c(x1[1],y1[1])))
      c2 <- Polygon(rbind(cbind(x2,y2),c(x2[1],y2[1])))
      spc1 <- SpatialPolygons(list(Polygons(list(c1),ID="1")))
      spc2 <- SpatialPolygons(list(Polygons(list(c2),ID="2")))
      slicePoly <- gIntersection(spc1,spc2) 
      return(slicePoly)
    }
  }
  
  #### NEW ####
  #Function that computes single ppa to speed up performance
  single.ppa <- function(aa,ePoints){
    cpX <- (aa$x[1] + aa$x[2])/2
    cpY <- (aa$y[1] + aa$y[2])/2
    #check to see if the angle is NA, if it is make it zero (no movement)
    if (is.na(aa$abs.angle[1]) == TRUE) {
      thetaRot <- 0
    } else {thetaRot <- aa$abs.angle[1]}
    
    c <- aa$dist[1]/2
    a <- (aa$dt[1]*aa$dynVmax[1])/2
    b <- sqrt((a^2)-(c^2))
    
    #Compute the ellipse
    poly <- ppaEllipse(cpX,cpY,a,b,thetaRot,ePoints)
  }
  #===========================================================================
  
  traj1 <- dynvmax(traj1,...)
  traj2 <- dynvmax(traj2,...)
  tr1 <- ld(traj1)
  tr2 <- ld(traj2)
  
  if (t.int < 0){stop(paste('Time interval is too small: ',t.int))}
  t.start <- max(c(min(tr1$date),min(tr2$date)))
  t.end <- min(c(max(tr1$date),max(tr2$date)))
  if (t.start > t.end){stop(paste('Error: The two ltraj objects have no temporal overlap.'))}
  tt <- seq(t.start+t.int, t.end-t.int, by=t.int)
  
  #print(paste('Estimated time for processing is:',length(tt)*0.0125/60,'minutes.'))
  anchor.frame <- data.frame(tt=tt,A1=0,A2=0,chunk=0)
  anchor.frame$A1 <- unlist(lapply(tt,get.anchor,traj=tr1))
  anchor.frame$A2 <- unlist(lapply(tt,get.anchor,traj=tr2))
  anchor.frame$chunk <- as.integer(as.numeric(factor(with(anchor.frame, paste(A1,'.', A2)))))
  
  #How can this be sped up?
  chunks <- unique(anchor.frame$chunk)
  poly.list <- vector('list',length(chunks))
  data <- data.frame(a1.time=NULL,a2.time=NULL)
  for (i in 1:length(chunks)){
    ind.c <- which(anchor.frame$chunk == chunks[i])
    ind.a1 <- anchor.frame$A1[ind.c]
    ind.a2 <- anchor.frame$A2[ind.c]
    a1 <- get.anchors(traj=tr1,indo=ind.a1[1])
    a2 <- get.anchors(traj=tr2,indo=ind.a2[1])
    
    ### NEW ### Speed up by checking if PPA overlaps
    #---------------------------------------------------
    pp1 <- SpatialPolygons(list(Polygons(list(single.ppa(a1,ePoints)),ID='1')),proj4string=proj4string)
    pp2 <- SpatialPolygons(list(Polygons(list(single.ppa(a2,ePoints)),ID='2')),proj4string=proj4string)
    if (gIntersects(pp1,pp2)){
      pol.list <- sapply(tt[ind.c],prism.slice.int,a1=a1,a2=a2,tol=tol,ePoints=ePoints)
      ind <- which(sapply(pol.list,is.null,simplify=TRUE,USE.NAMES=FALSE)==FALSE)
      if (length(ind)>0){
        pol.list <- unlist(pol.list[ind])
        for (j in 1:length(pol.list)){
          slot(pol.list[[j]],'ID') <- as.character(j)
        }
        sp.raw <- SpatialPolygons(pol.list)
        poly.hull <- gConvexHull(sp.raw,id=as.character(i))@polygons
        poly.list[i] <- poly.hull
        data <- rbind(data,c(a1$date[1],a2$date[1]))
      }
    }
    #------------ end speed up addition ---------------
  poly.list <- poly.list[!sapply(poly.list, is.null)]
  }
  
  
  if (length(poly.list) > 0){
    sp.time <- SpatialPolygonsDataFrame(SpatialPolygons(poly.list,proj4string=proj4string),data=data,match.ID=FALSE)
    spdf.diss <- gUnaryUnion(sp.time)
    spdf.diss2 <- SpatialPolygonsDataFrame(spdf.diss,data=data.frame(id=1:length(spdf.diss)),match.ID=F)
    if (dissolve == TRUE){
      return(spdf.diss2)
    } else {
      return(sp.time)
    }
  }else {return(NULL)}
}

#=====================================================================