# ---- roxygen documentation ----
#' @title PPA Ellipse
#'
#' @description
#'   Internal ellipse calculation function.
#' @details
#'   Internal function for calculating ellipses in time geographic analysis.
#' @param x first coordinate
#' @param y second coordinate
#' @param a semi-major axis
#' @param b semi-minor axis
#' @param theta rotation angle of the ellipse (in radians)
#' @param steps number of segments, from ePoints parameter in \code{dyn.ppa.hr}

#' @return
#'   This function returns a polygon ellipse.
#'
# @references
#' @keywords internal 
#' @seealso dynppa
# @examples
# @export
#
# ---- End of roxygen documentation ----
#===============================================================================
# function: ppaEllipse
# purpose: Calculate Time geography ellipses
#---------------------------------------------------------------
#x,y are center points
#a,b are semi-major and semi-minor axis
#theta is rotation angle (in radians)
#steps is number of segments to draw (default: 360 from ePoints parameter)
# This formulation uses the general parametric form of an ellipse
ppaEllipse <- function(x,y,a,b,theta,steps){
  X=rep(0,steps);Y=rep(0,steps)
  sinTheta <- sin(theta)
  cosTheta <- cos(theta)
  
  for (i in 1:steps){
    alpha <- i*(360/steps)*(pi/180)
    sinAlpha <- sin(alpha)
    cosAlpha <- cos(alpha)
    X[i] = x + (a*cosAlpha*cosTheta - b*sinAlpha*sinTheta)
    Y[i] = y + (a*cosAlpha*sinTheta + b*sinAlpha*cosTheta)
  }
  #save X,Y as Points Data frame
  ptDF <- data.frame(x=X,y=Y)
  ptDF[steps+1,] <- ptDF[1,]     #add first point to end to close ellipse
  #turn ellipses from points to polygons
  ellipsePoly <- Polygon(ptDF,hole=F)
  return(ellipsePoly)
}
#End of Function  
#===============================================================================