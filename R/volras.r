# ---- roxygen documentation ----
#' @title Volume contour from Raster
#'
#' @description
#'   Compute a percent volume contour polygon from a raster UD.
#' @details
#'   The volras function is a simpler version of the getvolumeUD function from the package \code{adehabitatHR} developed by C. Calenge. It allows the output to be a 'raster looking' polygon (i.e., the cells that are within the UD) or a simplified (smoothed) polygon.
#'   
#' @param x a \code{RasterLayer}
#' @param percent a percent value to get the volume contour, e.g., 95. Note: This is a simple function and only accepts one value at a time.
#' @param simplify (logical; default = TRUE) whether or not to simplify the output home range polygon using \code{gSimplify} with a tolerance value of 1.5 times the spatial resolution of the UD. 
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

volras <- function(x,percent=95,simplify=TRUE){
  
  x[is.na(x)] <- 0
  pfs <- proj4string(x)
  
  ## standardize it so that the total volume is 1 over the area
  v <- as.vector(values(x))
  index<-1:length(v)
  vord<-v[order(v, decreasing=TRUE)]
  vsu<-cumsum(vord)
  
  cont <- which(vsu > (percent/100)*max(vsu))[1]
  cont.lev <- vord[cont]
  
  #Get all the cells above the cont.lev and make 1
  m <- c(0, cont.lev, 0, cont.lev,cellStats(x,'max'),1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  x2 <- reclassify(x, rclmat)
  
  #Convert to polygon and simplify if desired
  hr <- rasterToPolygons(x2,fun=function(x){x==1},n=4,na.rm=TRUE,dissolve=TRUE)
  if(simplify){
    hr <- gSimplify(hr,tol=res(x2)[1]*1.5)
  }
  
  #return the volume contour polygon
  return(hr)
}

  
         
         
         