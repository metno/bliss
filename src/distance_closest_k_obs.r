#+ mean radial distance between an observation and its k-th closest obs
distance_closest_k_obs<-function(i,
                                 k=10,
                                 sort_decreasing=FALSE) {
#------------------------------------------------------------------------------
# mode == "analysis"
# global variables
#  obs_x/y. observation coordinates
#  grid_x/y. grid point coordinates
#------------------------------------------------------------------------------
  ifelse(length(obs_x)<k,NA,
         sqrt(sort( (obs_x-grid_x[i])**2 + (obs_y-grid_y[i])**2, 
              decreasing = sort_decreasing)[k]) )
}

