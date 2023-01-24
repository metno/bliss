#+ Define correlations 
corr1d <- function( values, par, label) {
# correlations when input vectors are referred to a specific point (e.g. distances from a point)
#
# par is the vector of parameters (1 for each vector of values to be processed)
# label is the vector of labels identifying the specific elaboration
# values is the vector/matrix of values to be processed (one vector is the column of an array)
#-----------------------------------------------------------------------------

  res <- rep( 1, length(value[,1]))

  for (i in 1:length(par)) {

    if (label[i] == "gaussian") {
      res <- res * exp( -0.5* (value[,i]*value[,i]) / (par[i]*par[i]) )

    } else if (label[i] == "soar")  {
      res <- res * (1+value[,i]/par[i])*exp(-value[,i]/par[i])

    } else if (label[i] == "powerlaw")  {
      res <- res * 1 / (1 + 0.5*(value[,i]*value[,i])/(par[i]*par[i]))

    } else if (label[i] == "toar")  {
      res <- res * (1 + value[,i]/par[i] + (value[,i]*value[,i])/(3*par[i]*par[i])) * exp(-value[,i]/par[i])

    } else if (label[i] == "linear")  {
      res <- res * (1 - (1-par[i]) * abs(value[,i]))

    } else {
      res <- rep( NA, length(value[,i]))

    } # end if

  } # end for

  res
}

#  if ( length(par) == 1) {
#    dist <- values
#    dh   <- par
#    if (label == "gaussian") {
#      res <- exp( -0.5* (dist*dist) / (dh*dh) )
#    } else if (label == "soar")  {
#      res <- (1+dist/dh)*exp(-dist/dh)
#    } else if (label == "powerlaw")  {
#      res <- 1 / (1 + 0.5*(dist*dist)/(dh*dh))
#    } else if (label == "toar")  {
#      res <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
#    } else {
#      res <- rep( NA, length(dist))
#    }
#  } else if ( length(par) == 2 ) {
#    dist <- values[,1]
#    dh   <- par[1]
#    elev_diff <- values[,2]
#    dz   <- par[2]
#    laf_diff <- values[,3]
#    lafmn   <- par[3]
#
#  } else if ( length(par) == 3 ) {
#
#  }
#  res
#}
