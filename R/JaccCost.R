# ---- roxygen documentation ----
#' @title Accumulated cost function - Modified by Jed
#'
#' @description
#'   Modified Accumulated cost function.
#' @details
#'   Modified Accumulated Cost function from the gdistance accCost function. 
#'   
#' @param x a \code{TransitionLayer} object
#' @param fromCoords location coordinates for accumulated cost calculation
#' @param mode one of \code{'in'} or \code{'out'} signifying if the direction of accumulation is towards or away from the location. 
#' 
#' @return
#'   This function returns a \code{RasterLayer} which has the accumulated costs.
#'
# @references
#' @keywords internal 
# @seealso ppa, fbtgTS, volras
# @examples
# @export
#
# ---- End of roxygen documentation ----

setGeneric("JaccCost", function(x, fromCoords,mode) standardGeneric("JaccCost"))

setMethod("JaccCost", signature(x = "TransitionLayer", fromCoords = "Coords",mode="character"), def = function(x, fromCoords, mode)
	{
  coordsToMatrix <- function(Coords){
    if(class(Coords) == "numeric")
    {
      if(length(Coords) == 2) {Coords <- t(as.matrix(Coords))} 
      else{stop("coordinates given as a vector, but the vector does not have a length of two")}
    }
    
    if(class(Coords) == "matrix")
    {
      if(!(ncol(Coords) == 2)){stop("coordinates given as a matrix, but the matrix does not have two columns")}
    }	
    
    if(inherits(Coords, "SpatialPoints"))  
    {
      Coords <- coordinates(Coords)
    }
    return(Coords)
  }
    fromCoords <- coordsToMatrix(fromCoords) 
		fromCells <- cellFromXY(x, fromCoords)
		if(!all(!is.na(fromCells))){
			warning("some coordinates not found and omitted")
			fromCells <- fromCells[!is.na(fromCells)]
		}
		tr <- transitionMatrix(x)

		adjacencyGraph <- graph.adjacency(tr, mode="directed", weighted=TRUE)
		E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight		
		
		dd <- distances(adjacencyGraph,v=fromCells,to=V(adjacencyGraph),mode=mode)

		result <- as(x, "RasterLayer")
		result <- setValues(result, dd[1,])

		return(result)
	}
)




