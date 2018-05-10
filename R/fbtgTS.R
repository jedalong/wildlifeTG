# ---- roxygen documentation ----
#' @title Field-based Time Geography Utilization Distribution
#'
#' @description
#'   Calculate the utilization distribution of an animal using the field-based time geographic model.
#' @details
#'   Calculates the field-based time geography probability surface for an animal for a chosen time - often called a time slice. Field-based time geography is based on an underlying resistance surface, which constratins potential movement by the animal between two location fixes. The output probability surface represents the probability of finding an animal in a given location at the input time, based on the field-based time geography movement model.  
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param tk time to compute the probability surface, \code{ta < tk < tb}.
#' @param tl a \code{TransitionLayer} object
#' @param timefun method for convertingtime into probability; one of:\cr 
#' -- \code{'inverse'} \eqn{= \frac{1}{c_2 * t}},\cr
#' -- \code{'inverse2'} \eqn{= \frac{1}{(c_2 * t)^2}},\cr
#' -- \code{'exp'} \eqn{= \exp{-c_2 * t}},\cr
#' -- \code{'norm'} \eqn{= \exp{-c_2 * t^2}},\cr
#' -- \code{'rootexp'} \eqn{= \exp{-c_2 * \sqrt{t}}},\cr
#' -- \code{'pareto'} \eqn{= \exp{-c_2 * \log{t}}},\cr
#' -- \code{'lognorm'}\eqn{= \exp{-c_2 * \log{t^2}}}.\cr
#' @param sigma locational uncertainty parameter for Guassian error kernel, default=0.
#' @param c2 Parameter input into \code{timefun} default=1, see details.
#' @param clipPPS (logical) whether or not the output probabilities should be clipped to the potential path space, default=TRUE.
#'
#' @return
#'   This function returns a \code{RasterLayer} which can be used to estimate the UD of an animal.
#'
# @references
# @keywords 
#' @seealso ppa, fbtgTS, volras
# @examples
#' @export
#
# ---- End of roxygen documentation ----

fbtgTS <- function(traj,tk,tl,timefun='inverse',sigma=0,c2=1,clipPPS=TRUE){
  
  df <- ld(traj)
  
  #get index associated with tk segment
  ind <- tail(which(df$date < tk),n=1)
  
  if (length(ind) == 0){ stop('tk value was not appropriately defined.') }
  if (ind == dim(df)[1]){ stop('tk value was not appropriately defined.') }
  
  A <- as.numeric(df[ind,1:2])
  B <- as.numeric(df[ind+1,1:2])
  t <- as.numeric(df$dt[ind],units='secs')
  ta <- as.numeric(tk-df$date[1],units='secs')
  tb <- as.numeric(df$date[2]-tk,units='secs')

  #Tai <- JaccCost(tl,A,'out')
  #Tib <- JaccCost(tl,B,'in')
  Ai <- cellFromXY(raster(tl),A)
  Bi <- cellFromXY(raster(tl),B)
  
  #Compute Accumulated cost (i.e, time) for location A and B
  #This is the value input for various functions
  tm <- transitionMatrix(tl)
  gr <- graph.adjacency(tm, mode="directed", weighted=TRUE)
  E(gr)$weight <- 1/E(gr)$weight
  Tai <- distances(gr,v=Ai,to=V(gr),mode='out')
  Tib <- distances(gr,v=Bi,to=V(gr),mode='in')
  
  #compute the cost of the shortest path
  Tshort <- costDistance(tl,A,B)[1]
  
  #Compute the time slice using internalTS function
  Pt <- internalTS(ta,t,Tai,Tib,raster(tl),Tshort,sigma,timefun,c2,clipPPS)
  #Make raster
  Pt <- setValues(raster(tl),Pt)          
  
  #Pt[Pt == 0] <- NA
  return(Pt)
}


