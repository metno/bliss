#+ Preprocessing step: Select background fields using reference gridded observations
preproc_selensemble <- function( argv, fg_env, env) {
#
# output
#  fg_env$score - score for each potential background field (higher the better)
#  fg_env$ixf - data source for each potential background field
#  fg_env$ixe - ensemble member number for each potential background field
#  fg_env$ixs - selection of the k background fields (indexes wrt fg_env$score)
#  y_env$yo(v)$value_fgsel - matrix( y_env$yo(v)$n, env$k_dim)
#------------------------------------------------------------------------------
  options( warn = 2)

  cat( "-- preproc selection of ensemble members --\n")

  if ( fg_env$nfg == 0) return( FALSE)
  
  ixf <- integer(0)
  ixe <- integer(0)

  # no need to compute score, take the first k_dim
  if (argv$selens_mode == "identity") {
    for (f in 1:fg_env$nfg) {
      for (e in 1:fg_env$fg[[f]]$k_dim) {
        ixf <- c( ixf, f)
        ixe <- c( ixe, e)
      }
    }
    ixs <- 1:env$k_dim
    score <- NA

  # compute score for each potential background field
  } else if (argv$selens_mode %in% c( "ets", "maxoverlap")) {
    cat (" compute score for each potential background field \n")
    # reference observed field
    if (argv$preproc_mergeobs) {
      uo <- env$mergeobs$value
      uo[uo<y_env$rain]  <- 0
      uo[uo>=y_env$rain] <- 1
    } else {
      yo <- y_env$yo$value
      yo[yo<y_env$rain]  <- 0
      yo[yo>=y_env$rain] <- 1
    }
    score <- numeric(0)
    for (f in 1:fg_env$nfg) {
      cat( paste("  data source",f,"ensembles"))
      for (e in 1:fg_env$fg[[f]]$k_dim) {
        xb <- getValues( subset( fg_env$fg[[f]]$r_main, subset=e))
        xb[xb<y_env$rain] <- 0
        xb[xb>=y_env$rain] <- 1
        if (argv$preproc_mergeobs) {
          ix <- which( !is.na(xb) & !is.na(uo) & is.finite(xb) & is.finite(uo))
          a <- as.numeric( length( which(xb[ix]==1 & uo[ix]==1)))
          b <- as.numeric( length( which(xb[ix]==1 & uo[ix]==0)))
          c <- as.numeric( length( which(xb[ix]==0 & uo[ix]==1)))
          d <- as.numeric( length( which(xb[ix]==0 & uo[ix]==0)))
        } else {
          yb <- extract( subset( fg_env$fg[[f]]$r_main, subset=e), 
                         cbind(y_env$yo$x, y_env$yo$y))
          ix <- which( !is.na(yb) & !is.na(yo) & is.finite(yb) & is.finite(yo))
          a <- as.numeric( length( which(yb[ix]==1 & yo[ix]==1)))
          b <- as.numeric( length( which(yb[ix]==1 & yo[ix]==0)))
          c <- as.numeric( length( which(yb[ix]==0 & yo[ix]==1)))
          d <- as.numeric( length( which(yb[ix]==0 & yo[ix]==0)))
        }
        if (argv$selens_mode == "ets") {
          a_random <- as.numeric( (a+c)*(a+b) / (a+b+c+d))
          score <- c( score, (a-a_random) / (a+c+b-a_random))
        } else if (argv$selens_mode == "maxoverlap") {
          score <- c( score, mean( a, na.rm=T))
        }
        ixf <- c( ixf, f)
        ixe <- c( ixe, e)
        cat( ".")
      } # loop over ensembles in one file
      cat( "\n")
    } # loop over files

    # select only the k fields with the highest scores
    ixs <- order( score, decreasing=T)[1:env$k_dim]
    cat( paste( " number of background ensemble members (tot) =", 
                env$k_dim,"(",length(ixf),")\n"))
  } # end selection of ensembles 

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
      y_env$yov$value_fgsel[,e] <- extract( 
       subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), 
       cbind( y_env$yov$x, y_env$yov$y) ) 
    }
  }

  y_env$yo$value_fgsel  <- array( data=NA, dim=c( y_env$yo$n, env$k_dim))
  for (e in 1:env$k_dim) { 
    i <- fg_env$ixs[e]
    y_env$yo$value_fgsel[,e] <- extract( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), cbind( y_env$yo$x, y_env$yo$y)) 
  }

  return( TRUE)

}
