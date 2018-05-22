#' Land cover data for northern BC, Canada
#'
#' Land cover data to be used in examples and vignettes alongside the caribou telemetry data. 
#' 
#' Landcover data see Wulder et al. (2008) that has been re-sampled using a modal filter to a 100 m resolution and clipped to the area where the telemetry data is present.
#'
#'
#' @docType data
#' @keywords datasets
#' @format A \code{RasterLayer} with 275 rows and 347 columns.
#' @references 
#' Wulder, M., White, J., Cranny, M., Hall, R., Luther, J., Beaudoin, A. (2008) Monitoring Canada's forests-Part 1: Completion of the EOSD land cover project. Canadian Journal of Remote Sensing. 34(6): 549-562.
#' @name eosd
#' @examples
#' data(eosd)
#' raster::plot(eosd)
#-----------------------------------
NULL