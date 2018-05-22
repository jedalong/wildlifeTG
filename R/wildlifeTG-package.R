# This is package documentation for wildlifeTG.
# roxygen will use this file to create a NAMESPACE file.
# Of importance is the @import command, as it lists package dependencies.

#' wildlifeTG - Time Geographic Analysis of Wildlife Telemetry Data
#'
#' The package \code{wildlifeTG} provides tools for performing time geographic analysis of animal tracking data. The functions provide useful tools for examining wildflife movement from the context of accessibility. The time geographic framework is an alternative view on typical home range estimation procedures. Currently, the package provides functions that faciliate the calculation of the PPA and dynPPA measures of an individuals accessibility space along with some newer methods for mapping utilization distributions that simultaneously consider barriers to movement. Please see the wildlifeTG on GitHub \url{https://github.com/jedalong/wildlifeTG} for the latest information, example code/vignettes, and the most up-to-date version of the package.
#'
#' \code{wildlifeTG}'s functions utilize the \code{ltraj} objects from the package \code{adehabitat}. 
#'
#' @author Jed Long
# @references
#' @import methods sp adehabitatLT adehabitatHR rgeos classInt Matrix devtools
#' @rawNamespace import(gdistance, except = normalize)
#' @rawNamespace import(raster, except = union)
#' @rawNamespace import(igraph, except = union)
#' @docType package
#' @name wildlifeTG-package
NULL
