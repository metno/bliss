# Adapted From Yue (Michael) Ying repository on github
# https://github.com/myying/QG_Multiscale_DA/blob/master/util.py
# Main modifications are for making this work in R and for using raster
# GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

warp <- function( r, u, v, method="simple") {
  xy <- xyFromCell( r, 1:ncell(r))
  x <- xy[,1] + getValues(u)
  y <- xy[,2] + getValues(v)
  s <- r 
  s[] <- extract( r, cbind(x,y), method=method)
  s
}

coarsen <- function( r, level) {
  if (level < 1) return(r)
  for (k in 1:level) {
    s <- aggregate( r, fact=2, fun=mean, expand=T, na.rm=T)
    r <- s
  }
  r
}

sharpen <- function( r_from, r_to, method="ngb") {
  resample( r_from, r_to, method=method )
}

deriv_x <- function(mat) {
#  t( 0.5 * apply( cbind( mat[,1], mat, mat[,ncol(mat)]), 1, diff, lag=2))
  0.5 * (shift_left(mat) - shift_right(mat))
}

deriv_y <- function(mat) {
#  0.5 * apply( rbind( mat[1,], mat, mat[nrow(mat),]), 2, diff, lag=2)
  0.5 * (shift_dw(mat) - shift_up(mat))
}

laplacian <- function(mat) {
  ( shift_left(mat) + shift_right(mat) + shift_dw(mat) + shift_up(mat)) / 6 +
  ( shift_left( shift_dw(mat)) +  shift_right( shift_dw(mat)) + shift_left( shift_up(mat)) + shift_right( shift_up(mat))) / 12 - mat
}

deriv_xy <- function(mat) {
  0.25 *  ( shift_dw( shift_left(mat)) + shift_up( shift_right(mat)) - shift_up( shift_left(mat)) - shift_dw( shift_right(mat)))
}

deriv_xx <- function(mat) {
  ( shift_left(mat) + shift_right(mat)) / 2 - mat
}

deriv_yy <- function(mat) {
  ( shift_dw(mat) + shift_up(mat)) / 2 - mat 
}

shift_left <- function(mat) {
  cbind( mat[,2:ncol(mat)], mat[,ncol(mat)])
}

shift_right <- function(mat) {
  cbind( mat[,1],mat[,1:(ncol(mat)-1)])
}

shift_up <- function(mat) {
  rbind( mat[2:nrow(mat),], mat[nrow(mat),])
}

shift_dw <- function(mat) {
  rbind( mat[1,], mat[1:(nrow(mat)-1),])
}
