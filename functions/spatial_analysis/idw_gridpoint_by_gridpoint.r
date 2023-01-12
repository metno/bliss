# Inverse Distance Weighting interpolation
idw_gridpoint_by_gridpoint<-function( i,
                                      pow_par=2) {
#------------------------------------------------------------------------------
  # define vectors
  dist <- sqrt((envtmp$x[i]-envtmp$obs_x)**2 + (envtmp$y[i]-envtmp$obs_y)**2)

  # some of the distances are equal to 0
  if ( (p <- length( ixa <- which(dist==0))) > 0) {
    xa <- mean(envtmp$obs_value[ixa])
  } else {
    w  <- 1 / (dist)**pow_par
    xa <- sum(w*envtmp$obs_value) / sum(w)
  }

  xa
}
