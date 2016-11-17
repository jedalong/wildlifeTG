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
#' @param simplify (logical) whether or not to simplify the output contour (uses \code{gSimplify}). 
#' 
#' @return
#'   A \code{SpatialPolygonsDataFrame}.
#'
# @references
#' @keywords internal 
# @seealso ppa, fbtgTS, volras
# @examples
#' @export
#
# ---- End of roxygen documentation ----

volras <- function(x,percent,simplify=FALSE){
  x[is.na(x)] <- 0
  pfs <- proj4string(x)
  sg <- as(x, "SpatialPixelsDataFrame")
  sg <- as(sg, "SpatialGridDataFrame")
  gr <- gridparameters(sg)
  uu <- names(sg)
  gri <- as.image.SpatialGridDataFrame(sg)
  xyok <- expand.grid(gri$y,gri$x)[,2:1]
  asc <- gri$z
  
  cs <- gr[1, 2]
  asc <- asc/(sum(asc)*cs*cs)
  
  ## computes the volume for each pixel
  ## thanks to a call to the C function calcvolume
  v<-.C("calcvolume", as.double(t(asc)), as.integer(ncol(asc)),
        as.integer(nrow(asc)), as.double(cs), PACKAGE="adehabitatHR")[[1]]
  
  ## standardize it so that the total volume is 1 over the area
  index<-1:length(v)
  vord<-v[order(v, decreasing=TRUE)]
  indord<-index[order(v, decreasing=TRUE)]
  vsu<-cumsum(vord)
  
  ## output
  vreord<-data.frame(n=vsu[order(indord)]*100)
  coordinates(vreord) <- xyok
  gridded(vreord) <- TRUE
  
  
  x<- as(vreord,'RasterLayer')
  #m <- matrix(c(0,percent,1,percent,100,NA),ncol=3,byrow=T)
  #x2 <- reclassify(x,m)
  if (simplify==TRUE){
    options(max.contour.segments=10000) 
    cc <- rasterToContour(x,levels=percent)
    c2 <- gSimplify(cc,tol=2*res(x)[1])
  } else {
    c2 <- rasterToContour(x,levels=percent)
  }
  return(c2)
}

  
         
         
         