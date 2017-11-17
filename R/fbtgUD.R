# ---- roxygen documentation ----
#' @title Field-based Time Geography Utilization Distribution
#'
#' @description
#'   Calculate the utilization distribution of an animal using the field-based time geographic model.
#' @details
#'   Calculates the field-based time geography utilization distribution (UD) for an animal. Field-based time geography is based on an underlying resistance surface, which constratins potential movement by the animal between two location fixes. This model is applied recursively over an entire trajectory in order to compute a UD for an animal. The UD inherently considers the movement limitations described by the underlying resistance surface, and thus is considered a landscape-based model for a UD. The landscape-based approach deviates from current models building upon random walks and diffusion processes. \cr
#'   The model requires that the resistance surface be directly related to an animals speed of passing through that environment.
#'   The timefun parameter is used to choose the model for converting time into a probability based on the described function. The \code{c2} parameter can be used to tune these functions based on some fine-scale movement data if available, but defaults to a value of 1. The sigma parameter represents the locational uncertainty, which can be interpreted as the standard deviation of the location error in a similar fashion to what is done in Brownian bridge models. The parameter k, which defines how many 'time slices' are to be computed, is the most significant influencer of computational time. Lower values will result speed up computations, but result in less-smooth output UD surfaces. 
#'   
#' @param traj animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}
#' @param tl a \code{TransitionLayer} object
#' @param timefun method for convertingtime into probability; one of:\cr 
#' -- \code{'inverse'} \eqn{= \frac{1}{c_2 * t}},\cr
#' -- \code{'inverse2'} \eqn{= \frac{1}{(c_2 * t)^2}},\cr
#' -- \code{'exp'} \eqn{= \exp{-c_2 * t}},\cr
#' -- \code{'norm'} \eqn{= \exp{-c_2 * t^2}},\cr
#' -- \code{'rootexp'} \eqn{= \exp{-c_2 * \sqrt{t}}},\cr
#' -- \code{'pareto'} \eqn{= \exp{-c_2 * \log{t}}},\cr
#' -- \code{'lognorm'}\eqn{= \exp{-c_2 * \log{t^2}}}.\cr
#' @param c2 Parameter input into \code{timefun} (see \code{likec2}) default=1.
#' @param sigma locational uncertainty parameter for Guassian error kernel, default=0.
#' @param k number of time slices between fixes upon which to estimate the UD, default=100. 
#' @param clipPPS (logical) whether or not the output probabilities should be clipped to the potential path space, default=TRUE.
#' @param d.min minimum distance, below which the segment is removed from analysis. Can be used to focus the UD on only longer movement segments. Default is the pixel size of \code{surf}.
#' @param dt.max maximum time, above which the segment is removed from analysis. Can be used to remove segments with missing data from analysis, as segments with a very long time between fixes can be problematic in time geographic analysis.


#' @return
#'   This function returns a \code{RasterLayer} which can be used to estimate the UD of an animal.
#'
# @references
# @keywords 
#' @seealso ppa, likec2
# @examples
#' @export
#
# ---- End of roxygen documentation ----

fbtgUD <- function(traj,tl,timefun='inverse',c2=1,sigma=0,k=100,clipPPS=TRUE,d.min=0,dt.max=Inf){
  
  #make a trajectory dataframe
  df <- ld(traj)
  
  #make an output raster
  Pi <- raster(tl)*0
  
  #Get only segments that meet the following criteria:
  # 1. dist > d.min  - minimum distance of movement (i.e., ignore non-movement segments)
  # 2. dt < dt.max   - maximum time difference (i.e., ignore segments with long durations)
  #ind <- which(df$dist > d.min & df$dt < dt.max)

  
  #Compute the values for each Time Slice using internalTS function
  n <- dim(df)[1]-1
  pb <- txtProgressBar(min=0,max=n,style=3)
  for (i in 1:n){
    Pi <- Pi + internalSEG(i, df, tl, k, sigma, timefun, c2, clipPPS)
    #update progress
    setTxtProgressBar(pb,i)
  }
  #P.V <- getValues(Pi)
  #P.I <- vapply(ind,internalSEG,FUN.VALUE=P.V, df, gr, k, sigma, timefun, c2, clipPPS)
  #Pi <- setValues(Pi, apply(P.I,1,sum))
    
  
  #Make 0's NA for easy plotting on output
  Pi[Pi == 0] <- NA
  
  return(Pi)
}


