# ---- roxygen documentation ----
#' @title Functions for converting time to probabilities
#'
#' @description
#'   Internal function for computing movement probabilities.
#' @details
#'   Used in the fbtg function. 
#'   
#' @param x a numeric object (i.e., a \code{RasterLayer}) upon which to compute the probabilities.
#' @param timefun method for convertingtime into probability; one of:\cr 
#' -- \code{'inverse'} \eqn{= \frac{1}{c_2 * t}},\cr
#' -- \code{'inverse2'} \eqn{= \frac{1}{(c_2 * t)^2}},\cr
#' -- \code{'exp'} \eqn{= \exp{-c_2 * t}},\cr
#' -- \code{'norm'} \eqn{= \exp{-c_2 * t^2}},\cr
#' -- \code{'rootexp'} \eqn{= \exp{-c_2 * \sqrt{t}}},\cr
#' -- \code{'pareto'} \eqn{= \exp{-c_2 * \log{t}}},\cr
#' -- \code{'lognorm'}\eqn{= \exp{-c_2 * \log{t^2}}}.\cr
#' @param c2 Parameter input into \code{timefun} default=1.
#' 
#' @return
#'   This function returns an object identical to \code{x} with \code{timefun} applied to its values.
#'
#' @keywords internal 
# ---- End of roxygen documentation ----

internalPfun <- function(x,timefun,c2,dt){
  
  #Check if c2 is
  if (is.na(c2)){
    c2 <- switch(timefun,
                 inverse = 1,
                 inverse2 = 1,
                 exp = dt^(-1),
                 norm = dt^(-2),
                 rootexp = dt^(-0.5),
                 pareto = log(dt)^(-1),
                 lognorm = log(dt)^(-2),
                 stop(paste('The time function',timefun,'does not exist.'))
    )
  }
  
  #Functions to get probability relative to an expected continuous movement along the LCP
  Pt <- switch(timefun,
               inverse = 1 / (x+c2),
               inverse2 = 1 / ((x+c2)^2),
               exp = exp(-c2*(x)),
               norm = exp(-c2*((x)^2)),
               rootexp = exp(-c2*((x)^0.5)),
               pareto = exp(-c2*log(x)),
               lognorm = exp(-c2*(log(x)^2)),
               stop(paste('The time function',timefun,'does not exist.'))
  )
  return(Pt)
}