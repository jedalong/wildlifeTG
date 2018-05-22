# ---- roxygen documentation ----
#' @title Dynamic Calculation of the Vmax Parameter
#'
#' @description
#'   The function \code{dynvmax} computes a dynamic version of the Vmax parameter for the PPA method. It can be used to incorporate changes in animal movement behaviour into the PPA method caluculation to better model that area accessible to an individual animal given the set of known telemetry locations in space and time.
#'   
#' @details
#'   The function \code{dynvmax} represents an intermediary function used to extend and improve upon an existing PPA home range method (Long and Nelson, 2012) as described in the paper (Long and Nelson, 2014). Four options are available for computing the vmax parameter dynamically and are passed into the \code{dynvmax} function using \code{dynamic} option.\cr 
#'  \cr 1) \code{NA} -- if \code{dynamic} = \code{'NA'} (the default) the function estimates the original, 
#'          non-dynamic estimate of Vmax which is a global estimate, as per Long & Nelson (2012). 
#'  \cr 2) \code{focal} -- a moving window approach whereby a window of size \code{w} is moved along the 
#'        trajectory and vmax computed dynamically within each window and assigned to the central segment.
#'  \cr 3) \code{cumulative} -- A moving window of size \code{w} is again used, only in this case the value 
#'        is assigned to the end segment. This represents the vmax calculation of the previous \code{w} segments.
#'  \cr 4) \code{class} -- A priori analysis (e.g., obtained via state-space models, or from expert knowledge) 
#'        is used to identify discrete behavioural states in the telemetry data and these stored in a column 
#'        which is then passed into the function.\cr\cr 
#'The \code{class} method is the preferred choice, as it allows the use of more sophisticated models for identifying behavioural shifts in telemetry data where we would expect to see clear differences in the Vmax parameter based on changing movement behaviour.\cr\cr 
#'The use of the \code{'focal'} or \code{'cumulative'} dynamic methods uses a moving window approach, which is sensitive to edge effects at the initial and ending times of the trajectory. Thus, the dynamic Vmax parameter is only computed for those segments that have a valid window and the dataset is shrunk by \code{w-1} segments.
#'   
#' @param traj an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param dynamic one of \code{'NA'}, \code{'focal'}, \code{'cumulative'}, or \code{'class'}; which signifies 
#'    whether or how to dynamically compute the Vmax parameter. See \bold{Details} for more information on 
#'    each of the choices.
#' @param method method for computing the Vmax parameter dynamically; can be one of several options:
#'    -- \code{"Robson"} for the Robson & Whitlock (1964) method,\cr
#'    -- \code{"RobsonLL"} for the R & W (1964) lower \eqn{(1-\alpha)*100\%} C.I. limit, \cr
#'    -- \code{"RobsonUL"} for the R & W (1964) upper \eqn{(1-\alpha)*100\%} C.I. limit, \cr
#'    -- \code{"vanderWatt"} for the van der Watt (1980) method, \cr
#'    -- \code{"vanderWattLL"} for the van der Watt (1980) lower \eqn{(1-\alpha)*100\%} C.I. limit, \cr
#'    -- \code{"vanderWattUL"} for the van der Watt (1980) upper \eqn{(1-\alpha)*100\%} C.I. limit. \cr
#' @param w (optional) window size (only used with \code{dynamic = 'focal'} or \code{'cumulative'}). 
#' @param class.col (optional) character indicating the name of the column in the \code{infolocs} dataframe
#'    of \code{traj} containing the categorized behavioural states of the animal (which can be stored as
#'    a character or numeric column).
#' @param k (optional) value for the \emph{k} parameter in the van der watt (1980) method; default is 5.
#' @param alpha (optional) value for the \eqn{\alpha} parameter if using upper or lower C.I. methods; default is 0.05.
#' @param manualVmax (optional) Character name of column in \code{traj} storing user input column of vmax values (typically call the column dynVmax).
#' @param vmaxtrunc (optional) due to irregular sampling intervals, or errors in GPS location, or other
#'    effects, the calculation of the vmax parameter through the statistical methods outlined above can be
#'    heavily influenced by high outliers. Thus, it may be useful to exclude those segments from calculation 
#'    of the dynamic Vmax parameter. Default is \code{NA}.
#' 
#' @return
#'   This function returns the original \code{traj} object with a new column -- \code{dynVmax} in the \code{infolocs} dataframe
#'   containing the dynamic vmax parameter for each trajectory segment.
#'
#' @references
#'   Long, JA, Nelson, TA. (2012) Time geography and wildlife home range delineation. \emph{Journal 
#'   of Wildlife Management}. 76(2):407-413.\cr\cr
#'   Long, JA, Nelson, TA. (2015) Home range and habitat analysis using dynamic time geography. \emph{Journal of -
#'   Wildlife Management}. 79(3):481-490.\cr\cr
#'   Robson, DS, Whitlock, JH. (1964) Estimation of a truncation point. \emph{Biometrika}
#'   51:33-39.\cr\cr
#'   van der Watt, P. (1980) A note on estimation bounds of random variables. \emph{Biometrika}
#'   67(3):712-714.\cr
# @keywords 
#' @seealso dynppa
#' @examples
#' m3R <- dynvmax(m3,dynamic='focal',method='Robson')
#' m3V <- dynvmax(m3,dynamic='focal',method='vanderWatt')
#' m3c <- dynvmax(m3,dynamic='cumulative')
#' 
#' @export
#
# ---- End of roxygen documentation ----

dynvmax <- function(traj,dynamic="NA",w=9,class.col="dt",method="Robson",k=5,alpha=0.05,manualVmax=NA,vmaxtrunc=NA){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  #check to see if k < w iff the vanderWatt methods are used
  #if (method == "vanderWatt" || method == "vanderWattLL" || method == "vanderWattUL"){
  #  if(k > w) {stop("The k parameter is greater than the w parameter.")}}

  #store trajectory data as a dataframe for indexing
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #compute segment velocities
  trDF$vi <- c(trDF$dist[1:(n-1)] / trDF$dt[1:(n-1)],NA)
  
  if (is.na(vmaxtrunc)){
    vmaxtrunc <- max(trDF$vi,na.rm=TRUE)
  }  
  
  trDF$dynVmax <- NA
  
  #original vmax calculation
  if (dynamic == "NA"){
    v.temp <- trDF$vi
    v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
    v <- sort(v.temp,decreasing=T)
    trDF$dynVmax <- switch(method,
      #Robson & Whitlock (1964) method --- main estimate and UL come
      # directly from the paper, however the LL was inferred
      Robson = v[1] + (v[1]-v[2]),
      RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
      RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
      #van der Watt (1980) method --- all come directly from the paper.
      vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
      vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
      vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
      stop(paste("The vMax method is spelt incorrectly: ",method)))
  }
  
  #Focal local Vmax
  if (dynamic == "focal"){
    w1 <- trunc(w/2)
    for (i in (w1+1):(n-1-w1)){
      v.temp <- trDF$vi[(i-w1):(i+w1)]
      v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
      v <- sort(v.temp,decreasing=T)
      
      trDF$dynVmax[i] = switch(method,
          #Robson & Whitlock (1964) method --- main estimate and UL come
          # directly from the paper, however the LL was inferred
          Robson = v[1] + (v[1]-v[2]),
          RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
          RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
          #van der Watt (1980) method --- all come directly from the paper.
          vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
          vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
          vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
          stop(paste("The vMax method is spelt incorrectly: ",method)))
      }
    }
    
  #cumulative local Vmax
  if (dynamic == "cumulative"){
    for (i in w:n){
      v.temp <- trDF$vi[(i-w+1):i]
      v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
      v <- sort(v.temp,decreasing=T)
      
      trDF$dynVmax[i] = switch(method,
          #Robson & Whitlock (1964) method --- main estimate and UL come
          # directly from the paper, however the LL was inferred
          Robson = v[1] + (v[1]-v[2]),
          RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
          RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
          #van der Watt (1980) method --- all come directly from the paper.
          vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
          vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
          vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
          stop(paste("The vMax method is spelt incorrectly: ",method)))
      }
    }

  #based on behaviour classes
  if (dynamic == "class"){
    #Check to see the input field matches one from the column
    colID <- which(names(trDF) == class.col)[1]
    if (is.na(colID)==TRUE) {stop(paste("The column name used for the class.col argument does not exist."))}
    #get the unique values of from class column
    classes <- unique(trDF[,colID])
    for (class in classes){
      ind <- which(trDF[,colID] == class)
      v.temp <- trDF$vi[ind]
      v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
      v <- sort(v.temp,decreasing=T)
      
      trDF$dynVmax[ind] = switch(method,
          #Robson & Whitlock (1964) method --- main estimate and UL come
          # directly from the paper, however the LL was inferred
          Robson = v[1] + (v[1]-v[2]),
          RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
          RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
          #van der Watt (1980) method --- all come directly from the paper.
          vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
          vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
          vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
          stop(paste("The vMax method is spelt incorrectly: ",method)))
      }
    }
  #User-defined vector of Vmax values as column in ltraj object.
  if (!is.na(manualVmax)){
    col <- which(names(trDF)==manualVmax)
    trDF$dynVmax <- trDF[,col]
  }
  
  #-----------------------------------------------------
  ## Process vmaxtrunc
  vmax.inflate <- 1.05   #could be an optional variable to pass in
  ind <- which(trDF$vi > vmaxtrunc)
  trDF$dynVmax[ind] <- trDF$vi[ind]*vmax.inflate
  
  #-----------------------------------------------------
  traj <- dl(trDF)
  return(traj)
  }

#============= End =============================================================