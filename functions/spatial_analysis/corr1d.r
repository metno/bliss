#+ Define correlations 
corr1d <- function( values, par, label, values2=NULL, values3=NULL) {
# correlations when input vectors are referred to a specific point (e.g. distances from a point)
#
# par is the vector of parameters (1 for each vector of values to be processed)
# label is the vector of labels identifying the specific elaboration
# values is the vector/matrix of values to be processed (one vector is the column of an array)
#-----------------------------------------------------------------------------
  res <- rep( 1, length(values))

  for (i in 1:length(par)) {
    if (i==2) values <- values2
    if (i==3) values <- values3
    if (label[i] == "gaussian") {
      res <- res * exp( -0.5* (values*values) / (par[i]*par[i]) )

    } else if (label[i] == "soar")  {
      res <- res * (1+values/par[i])*exp(-values/par[i])

    } else if (label[i] == "powerlaw")  {
      res <- res * 1 / (1 + 0.5*(values*values)/(par[i]*par[i]))

    } else if (label[i] == "toar")  {
      res <- res * (1 + values/par[i] + (values*values)/(3*par[i]*par[i])) * exp(-values/par[i])

    } else if (label[i] == "linear")  {
      res <- res * (1 - (1-par[i]) * abs(values))

    } else {
      res <- rep( NA, length(values))

    } # end if

  } # end for

  res
}
