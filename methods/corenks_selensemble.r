#+ Select the background fields using an observed field as the reference
corenks_selensemble <- function( argv, fg_env, env, plot=F, dir_plot=NA) {
#
# output
#  fg_env$score - score for each potential background field (higher the better)
#  fg_env$ixf - data source for each potential background field
#  fg_env$ixe - ensemble member number for each potential background field
#  fg_env$ixs - selection of the k background fields (indexes wrt fg_env$score)
#  y_env$yo(v)$value_fgsel - matrix( y_env$yo(v)$n, env$k_dim)
#------------------------------------------------------------------------------
  options( warn = 2)

  cat( "-- corenks selection of ensemble members --\n")

  if ( fg_env$nfg == 0) return( FALSE)

  # compute score for each potential background field
  cat (" compute score for each potential background field \n")
  # reference observed field
#  uo <- getValues( u_env$uo[[1]]$r_main)
#  uo[uo<u_env$rain] <- 0
#  uo[uo>u_env$rain] <- 1
  uo <- env$mergeobs$value
  uo[uo<y_env$rain] <- 0
  uo[uo>=y_env$rain] <- 1
  score <- numeric(0)
  ixf <- integer(0)
  ixe <- integer(0)
  j <- 0
  for (f in 1:fg_env$nfg) {
    cat( paste("  data source",f,"ensembles"))
    for (e in 1:fg_env$fg[[f]]$k_dim) {
      r <- getValues( subset( fg_env$fg[[f]]$r_main, subset=e))
#      r[r<u_env$rain] <- 0
#      r[r>u_env$rain] <- 1
      r[r<y_env$rain] <- 0
      r[r>=y_env$rain] <- 1
      a <- as.numeric( length( which(r==1 & uo==1)))
      b <- as.numeric( length( which(r==1 & uo==0)))
      c <- as.numeric( length( which(r==0 & uo==1)))
      d <- as.numeric( length( which(r==0 & uo==0)))
      if (argv$corenks_selens_mode == "ets") {
        a_random <- as.numeric( (a+c)*(a+b) / (a+b+c+d))
        score <- c( score, (a-a_random) / (a+c+b-a_random))
      } else if (argv$corenks_selens_mode == "maxoverlap") {
        score <- c( score, a/length(r))
      }
      ixf <- c( ixf, f)
      ixe <- c( ixe, e)
      j <- j+1
      cat( ".")
    }
    cat( "\n")
  } 

  cat( paste(" number of background ensemble members (tot) =", env$k_dim,"(",length(ixf),")\n"))
#  # select only the k fields with the highest scores
  if (argv$corenks_selens_mode == "identity") {
    ixs <- 1:env$k_dim
  } else {
    ixs <- order( score, decreasing=T)[1:env$k_dim]
  }

  # output
  fg_env$score <- score
  fg_env$ixf <- ixf
  fg_env$ixe <- ixe
  fg_env$ixs <- ixs

  # add info to y-structures
  if (env$cv_mode | env$cv_mode_random) {
    y_env$yov$value_fgsel  <- array( data=NA, dim=c( y_env$yov$n, env$k_dim))
    for (e in 1:env$k_dim) { 
      i <- fg_env$ixs[e]
      y_env$yov$value_fgsel[,e] <- extract( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), cbind( y_env$yov$x, y_env$yov$y)) 
    }
  }

  y_env$yo$value_fgsel  <- array( data=NA, dim=c( y_env$yo$n, env$k_dim))
  for (e in 1:env$k_dim) { 
    i <- fg_env$ixs[e]
    y_env$yo$value_fgsel[,e] <- extract( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), cbind( y_env$yo$x, y_env$yo$y)) 
  }

  return( TRUE)

}
