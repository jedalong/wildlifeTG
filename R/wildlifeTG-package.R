# This is package documentation for wildlifeTG.
# roxygen will use this file to create a NAMESPACE file.
# Of importance is the @import command, as it lists package dependencies.

#' wildlifeTG - Time Geographic Analysis of Wildlife Telemetry Data
#'
#' The package \code{wildlifeTG} provides tools for performing time geographic analysis of animal tracking data. The functions provide useful tools for examining wildflife movement from the context of accessibility. The time geographic framework is an alternative view on typical home range estimation procedures. Currently, the package provides functions that faciliate the calculation of the PPA and dynPPA measures of an individuals accessibility space. Please note that the package is still under development. Please see the wildlifeTG website at http://jedalong.github.io/wildlifeTG for the latest information and most up-to-date version of the package.
#'
#' \code{wildlifeTG}'s functions utilize the \code{ltraj} objects from the package \code{adehabitat}. 
#'
#' @author Jed Long
# @references
#'
#' @import methods sp gdistance adehabitatLT adehabitatHR rgeos raster classInt foreach doParallel
#' @docType package
#' @name wildlifeTG-package
NULL
