# TO DO LIST
# 3. Add more flexibility for definining output grid.
# 4. Add in tolerance? move to dynvmax?


#=========================================
# library(wildlifeTG)
# library(raster)
# library(adehabitatLT)
# library(sp)
# #=========================================

tgkde <- function(traj,disfun='inv',grid=NA,...){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  #Append the local Vmax values to the trajectory using the dynvmax function
  traj <- dynvmax(traj, ...)
  
  #convert to dataframe
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #-------------------------------------------
  # Create the output grid
  #-------------------------------------------
  x.range <- max(trDF$x) - min(trDF$x)
  y.range <- max(trDF$y) - min(trDF$y)
  
  if (is.na(grid)){
    grid <- 0.01*min(c(x.range,y.range)) 
  }
  xx <- seq(min(trDF$x)-0.1*x.range,max(trDF$x)+0.1*x.range,by=grid)
  yy <- seq(min(trDF$y)-0.1*y.range,max(trDF$y)+0.1*y.range,by=grid)
  xy <- expand.grid(x=xx,y=yy)
  xy$z <- 0
  #-------------------------------------------

  #loop through Traj object and compute the tgkde estimate
  x <- xy$x
  y <- xy$y
  for (i in 1:(n-1)){
    sx <- trDF$x[i]
    sy <- trDF$y[i]
    ex <- trDF$x[i+1]
    ey <- trDF$y[i+1]
    dt <- trDF$dt[i]
    vmax <- trDF$dynVmax[i]
    if (is.na(vmax)){next}
    dd <- sqrt((sx-x)^2 + (sy-y)^2) + sqrt((x-ex)^2 + (y-ey)^2)
    dp <- sqrt((sx-ex)^2 + (sy-ey)^2)
    ind <- which(dd > dt*vmax)
    #insert function here perhaps using switch?
    g <- switch(disfun,
                inv = dp/dd,              #inverse distance
                inv2 = (dp/dd)^2,         #inverse distance ^2        
                exp = exp(-dd/dp),        #exponential function 
                norm = exp(-(dd/dp)^2),   #normal function
                stop(paste('The distance decay function',disfun,'does not exist.'))
    )
      
    g[ind] <- 0
    #In Theory, need to normalize here, because each ellipse should sum to dt
    #g <- dt*g/sum(g)
    
    xy$z <- xy$z + g
  }
  
  #---- This all seems rather arbitrary? -----
  #Normalize following the eqn. in Downs et al. (2011)
  #get the 'average' vmax
  vmax. <- mean(trDF$dynVmax,na.rm=TRUE)
  #get the 'overall' time difference
  DT <- as.numeric(difftime(trDF$date[n], trDF$date[1], units='secs'))
  #Normalize the values
  xy$z <- xy$z * (1 / ((n-1)*(DT*vmax.)^2))
  
  # Downs et al. 2011 method does not sum to 1... 
  # based on normalization of each ellipse, just dividing by n-1 should work
  # xy$z <- xy$z / (n-1)
  
  #----------------------------
  #  Format output to Raster
  #----------------------------
  coordinates(xy) = ~x+y
  gridded(xy) <- TRUE
  ras <- raster(xy)
  #----------------------------
  return(ras)
}
#=============DEMO==========================
# traj <- simm.crw(1:500)
# trDF <- ld(traj)
# outras <- tgkde(traj,grid=0.25)
# #Draw the plot
# library(fields)
# image(outras,xlab="X",ylab="Y",asp=1)
# lines(trDF$x,trDF$y)
#============================================