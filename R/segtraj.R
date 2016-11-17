# ---- roxygen documentation ----
#' @title Simple trajectory segmentation
#'
#' @description
#'   Using the \code{classIntervals} function from the \code{classInt} package perform trajectory segmentation based on a variety of simple algorithms. 
#' @details
#'   The segementation is based on a single trajectory attribute (default is the movement distance). The default method for performing the segmentation is the Fisher-Jenks algorithm, but many are available (see the full list in the documentation for the \code{classIntervals} function).  
#'   
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first object. Note this object must be a \code{type II ltraj} object. For more information on objects of this type see \code{help(ltraj)},
#' @param k (integer) the number of output classes desired from the segmentation (default = 2),
#' @param style (character)the algorithm for identifying class breaks (default = 'fisher'), 
#' @param col the name (as a character string) of the column upon which to perform the segmentation (default = 'dist')
#' 
#' @return
#'   An \code{ltraj} object with a new infolocs column 'class'.
#'
# @references
# @keywords
# @seealso 
# @examples
#' @export
#
# ---- End of roxygen documentation ----
segtraj <- function(traj,k=2,style='fisher',col='dist'){
  trj <- ld(traj)
  x <- trj[,col]
  a <- classIntervals(x,n=k,style=style,na.rm=T)
  fac <- cut(x,breaks=a$brks,include.lowest=TRUE)
  levels(fac) <- 1:k #lowest to highest values
  trj$class <- fac
  traj <- dl(trj)
  return(traj)
}
