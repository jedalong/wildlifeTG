# ---- roxygen documentation ----
#' @title Potential path volume: 3D home range
#'
#' @description
#'   The function \code{ppv} computes the potential path volume (PPV) 3D home range. It is also a measure of the accessibility space of an animal in 3D. The PPV method considers the spatial and temporal constraints on movement given known telemetry fixes, and a (dynamic) measure of maximum mobility - termed Vmax.
#'   
#' @details
#'   The function \code{ppv} represents an extension to an existing 2D home range method (Long and Nelson, 2012) for 3D. ADD MORE HERE.
#'   
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first object. Note this object must be a \code{type II ltraj} object. It also requires a column in the \code{infolocs} which provides the 3rd Dimension (height) infromation. For more information on objects of this type see \code{help(ltraj)}.
#' @param zcol (character string) giving the name of the column in the \code{infolocs} where the height attribute data is stored.
#' @param tol parameter used to filter out those segments where the time between fixes is overly large (often due to irregular sampling or missing fixes); which leads to an overestimation of the activity space via the PPV method. Default is the maximum sampling interval from \code{traj}.
#' @param vmax either a single vmax value giving the vmax value, or a vector, length n-1, giving the vmax value for each segment.
#' @param vox size of voxels in map units (e.g., if {x,y,z} are in meters then vox should be in meters).
#'    
#' @return
#'   This function returns a XX representing the ppv 3D home range measure, which can be interpreted as the 3D accessibility space of an individual animal.
#'
#' @references
#'   Demsar, U, Long, JA (in review) Potential Path Volume (PPV): a geometric estimator for space use in 3D. In review at Movement Ecology. 
#'   
#'   Long, JA, Nelson, TA. (2012) Time geography and wildlife home range delineation. \emph{Journal of Wildlife
#'   Management}, 76(2):407-413.\cr \cr
#'   
# @keywords 
#' @seealso dynvmax dynppa
# @examples 
#' 
#' @export
# ---- End of roxygen documentation ----

ppv <- function(traj,zcol='z',vmax,vox){
  df <- ld(traj)
  
  n <- dim(df)[1]  #number of fixes
  
  # Set extent of the volume 
  minXcoord <- min(df$x) # westernmost point 
  maxXcoord <- max(df$x) # easternmost point
  minYcoord <- min(df$y) # southernmost point
  maxYcoord <- max(df$y) # northernmost point
  minZcoord <- min(df[,zcol]) # lowest point
  maxZcoord <- max(df[,zcol]) # highest point
  
  #distance and velocity in 3D
  df$dist3D <- NA
  for (i in 1:(n-1)){
    df$dist3D[i] <- sqrt( (df$x[i]-df$x[i+1])^2 + (df$y[i]-df$y[i+1])^2 + (df[i,zcol]-df[i+1,zcol])^2 )
  }
  df$v3D <- df$dist3D / df$dt # segment velocity in 3D
  
  
  # Ellipsoid parameters
  df$a <- (vmax * df$dt)/2 
  df$b <- sqrt(df$a^2-((df$dist3D^2)/4))   #throws warning
  
  # Coordinates of central point on segment P_C: xc_i, yc_i, zc_i and rotations
  df$xc <- 0           # central point on segment P_C, X coordinate 
  df$yc <- 0           # dentral point on segmentvP_C, Y coordinate 
  df$zc <- 0           # central point on segment P_C, Z coordinate 
  df$alpha <- 0        # first rotation angle, alpha_i 
  df$beta <- 0         # second rotation angle, beta_i 
  df$direction <- 0    # +/- for direction angle
  df$dZ <- 0
  df$dX <- 0
  
  # ----------------------------------------
  # Calculate the largest extent, depending on the largest b
  
  # Add a buffer in distance to extent to fit ellipsiods around end points
  # Buffer is the size of the largest b across all trajectories and segments
  
  distBuffer <- max(df$b,na.rm=T)
  
  distX <- maxXcoord - minXcoord
  distY <- maxYcoord - minYcoord
  distZ <- maxZcoord - minZcoord
  minXcoord <- minXcoord - distBuffer
  maxXcoord <- maxXcoord + distBuffer
  minYcoord <- minYcoord - distBuffer
  maxYcoord <- maxYcoord + distBuffer
  minZcoord <- minZcoord - distBuffer
  maxZcoord <- maxZcoord + distBuffer
  
  ##Automatically choose a default vox here
  
  # Build the global volume
  x <- seq(minXcoord,maxXcoord,by=vox)
  y <- seq(minYcoord,maxYcoord,by=vox)
  z <- seq(minZcoord,maxZcoord,by=vox)
  m <- expand.grid(x=x,y=y,z=z,p=0)
  
  
  #Loop through and process each segment - could we do better here!?!?
  for (i in 1:(n-1)){
    df$xc[i] <- (df$x[i] + df$x[i+1])/2 
    df$yc[i] <- (df$y[i] + df$y[i+1])/2
    df$zc[i] <- (df[i,zcol] + df[i+1,zcol])/2
    # Rotation angles alpha and beta 
    df$alpha[i] <- atan((df$y[i+1]-df$y[i]) / (df$x[i+1]-df$x[i]))   
    df$beta[i] <-  asin((df[i+1,zcol]-df[i,zcol]) / df$dist3D[i])
    
    # Rotation direction for beta depends on the direction of movement (up or down) along the segment and dx
    df$dZ[i] <- df[i+1,zcol] - df[i,zcol]
    df$dX[i] <- df$x[i+1] - df$x[i]
    df$direction[i] <- df$dZ[i] * df$dX[i]
    
    # Transform the coordinates by 1) translation to PC, 2) rotation for alpha, 3)rotation for beta
    
    y_new <- - sin(df$alpha[i])*(m$x - df$xc[i]) + cos(df$alpha[i])*(m$y - df$yc[i])
      
    if(df$direction[i]>=0) { # dz*dx is more than zero -> we turn up or down based on dz
      
      if (df$dZ[i] >0) { # Moving up: z(i)>z(i-1) -> Rotation 1
        x_new <- cos(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) + sin(df$beta[i])*(m$z - df$zc[i])
        z_new <- - sin(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) + cos(df$beta[i])*(m$z - df$zc[i])
      } else if (df$dZ[i]<0) { # Moving down: z(i)<z(i-1) -> Rotation 2
        x_new <- cos(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) - sin(df$beta[i])*(m$z - df$zc[i])
        z_new <- sin(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) + cos(df$beta[i])*(m$z - df$zc[i])
      } else { # Same level: z(i)=z(i-1) 
        x_new <- cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i])
        z_new <- (m$z - df$zc[i])
      }
      
      
    } else if (df$direction[i]<0) { # dz*dx is less than zero -> we turn down or up based on dz
      
      if (df$dZ[i]>0) { # Moving up: z(i)>z(i-1) -> Rotation 2
        x_new <- cos(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) - sin(df$beta[i])*(m$z - df$zc[i])
        z_new <- sin(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) + cos(df$beta[i])*(m$z - df$zc[i])
      } else if (df$dZ[i]<0) { # Moving down: z(i)<z(i-1) -> Rotation 1
        x_new <- cos(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) + sin(df$beta[i])*(m$z - df$zc[i])
        z_new <- - sin(df$beta[i])*( cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i]) ) + cos(df$beta[i])*(m$z - df$zc[i])
      } else { # Same level: z(i)=z(i-1) 
        x_new <- cos(df$alpha[i])*(m$x - df$xc[i]) + sin(df$alpha[i])*(m$y - df$yc[i])
        z_new <- (m$z - df$zc[i])
      }
    }
    #------------------------------------------------------------------------
    # Is the new point (x_new, y_new, z_new) inside the ellipsoid?
    inside <- (x_new/df$a[i])^2+(y_new/df$b[i])^2+(z_new/df$b[i])^2
    
    # It is, if inside <= 1, then we assign value 1 to the voxel in totalPPV: this will
    # produce a union of all the new voxels that were inside the PPV for this trajectory:
    # add them to the totalPPV volume:
    m$p[which(inside <= 1)] <- 1
  }

  #return the voxel dataframe - m
  return(m)
}

