# ---- roxygen documentation ----
#' @title Volume contour from Raster
#'
#' @description
#'   Compute a percent volume contour from a raster UD.
#' @details
#'   The volras function simply modifies the getvolumeUD function from the package \code{adehabitatHR} developed by C. Calenge. Basically, I was getting errors when using the getvolumeUD function, so I modified that function, and added the simplify polygon functionality.
#'   
#' @param x a \code{RasterLayer}
#' @param percent a percent value to get the volume contour, e.g., 95.
#' @param maxpixels numeric value, how many pixels to be included in countour line. More pixels results in a more detailed output polygon, but slower calcultion. 
#' 
#' @return
#'   A \code{SpatialPolygonsDataFrame}.
#'
# @references
# @keywords internal 
# @seealso 
# @examples
#' @export
#
# ---- End of roxygen documentation ----

volras <- function(x,percent=95,maxpixels=10000){
  
  #Function taken from Stack Overflow modified slightly
  ## http://stackoverflow.com/questions/14379828/how-does-one-turn-contour-lines-into-filled-contours#
  raster2contourPolys <- function(r, levels = NULL,percent=NULL,maxpixels=maxpixels) {
    
    ## set-up levels
    levels <- sort(levels)
    percent <- sort(percent)
    plevels <- c(min(values(r), na.rm=TRUE), levels, max(values(r), na.rm=TRUE)) # pad with raster range
    llevels <- paste(plevels[-length(plevels)], plevels[-1], sep=" - ")  
    llevels[1] <- paste("<", min(levels))
    llevels[length(llevels)] <- paste(">", max(levels))
    
    ## convert raster object to matrix so it can be fed into contourLines
    xmin <- extent(r)@xmin
    xmax <- extent(r)@xmax
    ymin <- extent(r)@ymin
    ymax <- extent(r)@ymax
    rx <- seq(xmin, xmax, length.out=ncol(r))
    ry <- seq(ymin, ymax, length.out=nrow(r))
    rz <- t(as.matrix(r))
    rz <- rz[,ncol(rz):1] # reshape
    
    ## get contour lines and convert to SpatialLinesDataFrame
    cat("Converting to contour lines...\n")
    cl <- rasterToContour(r,levels=levels,maxpixels=maxpixels) 
    cl <- SpatialLinesDataFrame(cl,data=data.frame(percent=percent),match.ID = F)
    
    ## extract coordinates to generate overall boundary polygon
    xy <- coordinates(r)[which(!is.na(values(r))),]
    i <- chull(xy)
    b <- xy[c(i,i[1]),]
    b <- SpatialPolygons(list(Polygons(list(Polygon(b, hole = FALSE)), "1")))
    
    ## add buffer around lines and cut boundary polygon
    cat("Converting contour lines to polygons...\n")
    bcl <- gBuffer(cl, width = 0.0001) # add small buffer so it cuts bounding poly
    cp <- gDifference(b, bcl)
    
    ## restructure and make polygon number the ID
    polys <- list() 
    for(j in seq_along(cp@polygons[[1]]@Polygons)) {
      polys[[j]] <- Polygons(list(cp@polygons[[1]]@Polygons[[j]]),j)
    }
    cp <- SpatialPolygons(polys)
    cp <- SpatialPolygonsDataFrame(cp, data.frame(id=seq_along(cp)))
    
    ## cut the raster by levels
    rc <- cut(r, breaks=plevels)
    
    ## loop through each polygon, create internal buffer, select points and define overlap with raster
    cat("Adding attributes to polygons...\n")
    l <- character(length(cp))
    for(j in seq_along(cp)) {
      p <- cp[cp$id==j,] 
      bp <- gBuffer(p, width = -max(res(r))) # use a negative buffer to obtain internal points
      if(!is.null(bp)) {
        xy <- SpatialPoints(coordinates(bp@polygons[[1]]@Polygons[[1]]))[1]
        l[j] <- llevels[extract(rc,xy)]
      } 
      else { 
        xy <- coordinates(gCentroid(p)) # buffer will not be calculated for smaller polygons, so grab centroid
        l[j] <- llevels[extract(rc,xy)]
      } 
    }
    
    ## assign level to each polygon
    cp$level <- factor(l, levels=llevels)
    cp$min <- plevels[-length(plevels)][cp$level]
    cp$max <- plevels[-1][cp$level]  
    cp <- cp[!is.na(cp$level),] # discard small polygons that did not capture a raster point
    df <- unique(cp@data[,c("level","min","max")]) # to be used after holes are defined
    df <- df[order(df$min),]
    row.names(df) <- df$level
    llevels <- df$level
    df$percent <- percent
    
    ## define depressions in higher levels (ie holes)
    cat("Defining holes...\n")
    spolys <- list()
    p <- cp[cp$level==llevels[1],] # add deepest layer
    p <- gUnaryUnion(p)
    spolys[[1]] <- Polygons(p@polygons[[1]]@Polygons, ID=llevels[1])
    for(i in seq(length(llevels)-1)) {
      p1 <- cp[cp$level==llevels[i+1],] # upper layer
      p2 <- cp[cp$level==llevels[i],] # lower layer
      x <- numeric(length(p2)) # grab one point from each of the deeper polygons
      y <- numeric(length(p2))
      id <- numeric(length(p2))
      for(j in seq_along(p2)) {
        xy <- coordinates(p2@polygons[[j]]@Polygons[[1]])[1,]
        x[j] <- xy[1]; y[j] <- xy[2]
        id[j] <- as.numeric(p2@polygons[[j]]@ID)
      }
      xy <- SpatialPointsDataFrame(cbind(x,y), data.frame(id=id))
      holes <- over(xy, p1)$id
      holes <- xy$id[which(!is.na(holes))]
      if(length(holes)>0) {
        p2 <- p2[p2$id %in% holes,] # keep the polygons over the shallower polygon
        p1 <- gUnaryUnion(p1) # simplify each group of polygons
        p2 <- gUnaryUnion(p2)
        p <- gDifference(p1, p2) # cut holes in p1      
      } else { p <- gUnaryUnion(p1) }
      spolys[[i+1]] <- Polygons(p@polygons[[1]]@Polygons, ID=llevels[i+1]) # add level 
    }
    cp <- SpatialPolygons(spolys, pO=seq_along(llevels), proj4string=CRS(proj4string(r))) # compile into final object
    cp <- SpatialPolygonsDataFrame(cp, df)
    cat("Done!")
    cp
    
  }
  
  x[is.na(x)] <- 0
  pfs <- proj4string(x)
  
  ## standardize it so that the total volume is 1 over the area
  v <- as.vector(values(x))
  index<-1:length(v)
  vord<-v[order(v, decreasing=TRUE)]
  vsu<-cumsum(vord)
  
  confun <- function(percent,vsu){
    cont <- which(vsu > (percent/100)*max(vsu))[1]
    cont.val <- vord[cont]
    return(cont.val)
  }
  
  cont.lev <- sapply(percent,confun,vsu)

  sl <- raster2contourPolys(x,levels=cont.lev,percent=percent,maxpixels=maxpixels)

  #return the volume contour
  return(sl)
}

  
         
         
         