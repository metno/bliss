#+ nearest neighbour interpolation 
ngb_gridpoint_by_gridpoint<-function( x_from, y_from, val_from, x_to, y_to,
                                      buffer=100000000) {
#------------------------------------------------------------------------------

  # initialization
  n_to <- length(x_to)
  val_to <- vector( mode="numeric", length= n_to)
  val_to[] <- NA

  # 
  nn2 <- nn2( cbind( x_from, y_from), 
              query = cbind( x_to, y_to), 
              k = 1, searchtype = "radius", 
              radius = buffer)
  ix <- which(nn2$nn.idx > 0)
  val_to[ix] <- val_from[nn2$nn.idx[ix]]

  # 
  val_to
}

