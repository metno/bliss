#+ Define correlations
corr2d <- function( values, par, label, values_as_globalVar=F, 
                    values2=NULL, values3=NULL) {
# correlations between pair of points
#------------------------------------------------------------------------------



  for (i in 1:length(par)) {
    
    if (values_as_globalVar) {
      if (i==1) { 
        mat2 <- envtmp$dist2 
      } else if (i==2) { 
        mat2 <- envtmp$dist2_z 
      } else if (i==3) { 
        mat2 <- envtmp$dist2_laf 
      }
    } else if (i==1) {
      mat2 <- values 
    } else if (i==2) {
      mat2 <- values2 
    } else if (i==3) {
      mat2 <- values3
    } # end if 

    if (i==1) { res <- mat2; res[] <- 1 }

    if (label[i] == "gaussian") {
      res <- res * exp( -0.5* mat2 / (par[i]*par[i]) )

    } else if (label[i] == "soar")  {
      mat_norm <- sqrt(mat2) / par[i]
      res <- res * (1+mat_norm) * exp(-mat_norm)

    } else if (label[i] == "powerlaw")  {
      res <- res * 1 / (1 + 0.5 * mat2 / (par[i]*par[i]))

    } else if (label[i] == "toar")  {
      mat <- sqrt(mat2)
      res <- res * (1 + mat/par[i] + mat2/(3*par[i]*par[i])) * exp(-mat/par[i])

    } else if (label[i] == "linear")  {
      mat <- sqrt(mat2)
      res <- res * (1 - (1-par[i]) * sqrt(mat2))

    } else {
      res <- rep( NA, length(value[,i]))

    } # end if

  } # end for

  res
}

#
#  if ( length( par == 1)) {
#    dh <- par
#    if (dist2_is_global) {
#      dist2 <- envtmp$dist2
#    } else {
#      x  <- values[,1]
#      y  <- values[,2]
#      dist2 <- outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2
#    }
#    if (label == "gaussian") {
#      res <- exp(-0.5 * dist2 / (dh*dh))
#    } else if (label == "soar")  {
#      distnorm <- sqrt(dist2) / dh 
#      res <- (1+distnorm)*exp(-distnorm)
#    } else if (label == "powerlaw")  {
#      res <- 1 / (1 + 0.5*dist2 / (dh*dh))
#    } else if (label == "toar")  {
#      dist <- sqrt(dist2)
#      res <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
#    }
#  }
#  res
#}
