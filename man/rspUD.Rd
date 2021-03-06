% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rspUD.R
\name{rspUD}
\alias{rspUD}
\title{Randomised Shortest Path Utilization Distribution}
\usage{
rspUD(traj, r, theta, theta.col, timescale = FALSE)
}
\arguments{
\item{traj}{animal movement trajectory in the form of an \code{ltraj} object, see package \code{adehabitatLT}}

\item{r}{a \code{RasterLayer} object describing the preference/affinity of the landscape. May be the result of a resource selection function, or other analyses. Note: Higher values should be associated higher affinity.}

\item{theta}{a global value for theta.}

\item{theta.col}{character string of the name of the column containing values of theta for each segment (only used if theta is not provided).}

\item{timescale}{(logical) whether or not to scale each segment so the sum of the output surface is equal to the duration of the segment (in seconds). Default = FALSE.}
}
\value{
This function returns a \code{RasterLayer} which can be used to estimate the UD of an animal.
}
\description{
Calculate the utilization distribution of an animal based on randomised shortest paths between fixes.
}
\details{
A randomised shortest path model is fit between every pair of fixes. The randomised shortest path model is derived from the \code{passage} function in the package \code{gdistance}. It uses the net number of packages (see \code{?passage}) to estimate the probability and then scales the values appropriately. An input \code{rasterLayer} object is required which defines the ability for movement through the landscape, which might typically be derived from a resource selection function, or be related to known barriers on the landscape. It requires only a single parameter \code{theta} which can be estimated from the data using the function \code{esttheta}.
}
\references{
Long, J.A. Estimating wildlife utilization distributions using randomized shortest paths. (in Preparation)
}
\seealso{
esttheta, volras
}
