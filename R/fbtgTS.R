# ---- roxygen documentation ----
#' @title Field-based Time Geography Utilization Distribution
#'
#' @description
#'   Calculate the utilization distribution of an animal using the field-based time geographic model.
#' @details
#'   Calculates the field-based time geography probability surface for an animal for a chosen time - often called a time slice. Field-based time geography is based on an underlying resistance surface, which constratins potential movement by the animal between two location fixes. The output probability surface represents the probability of finding an animal in a given location at the input time, based on the field-based time geography movement model.  
#' @param A coordinates of origin point of a movement segment in the form \code{c(x,y,t)}.
#' @param B coordinates of end point of a movement segment in the form \code{c(x,y,t)}.
#' @param tk time to compute the probability surface, \code{ta < tk < tb}.
#' @param timefun method for converting movement time into probability; one of: \code{'inverse' (default),'inverse2','exp','norm','rootexp','pareto','lognorm'}, see details.
#' @param c2 Parameter input into \code{timefun} default=1, see details.
#' @param k number of time slices between fixes upon which to estimate the UD, default=10.
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

fbtgTS <- function(traj,tk,surf,sigma,timefun='inverse',c2=1,clipPPS=TRUE){
  
  df <- ld(traj)
  
  #get index associated with tk segment
  ind <- tail(which(df$date < tk),n=1)
  
  if (length(ind) == 0){ stop('tk value was not appropriately defined.') }
  if (ind == dim(df)[1]){ stop('tk value was not appropriately defined.') }
  
  A <- as.numeric(df[ind,1:2])
  B <- as.numeric(df[ind+1,1:2])
  t <- as.numeric(df$dt[ind],units='secs')

  cai <- JaccCost(surf,A,'out')
  cib <- JaccCost(surf,B,'in')
  ta <- as.numeric(tk-df$date[1],units='secs')
  tb <- as.numeric(df$date[2]-tk,units='secs')
   
  #compute the cost of the shortest path
  Tshort <- costDistance(surf,A,B)[1]
  
  #Compute the time slice using internalTS function
  Pt <- internalTS(ta,t,cai,cib,Tshort,sigma,timefun,c2,clipPPS)
  
  Pt <- setValues(raster(surf),Pt)
  
  #Pt[Pt == 0] <- NA
  return(Pt)
}


