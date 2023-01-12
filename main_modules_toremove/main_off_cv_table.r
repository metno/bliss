#+
write_cv_table <- function( argv, y_env) {
#------------------------------------------------------------------------------
    
    y_env$yov$y      <- data$y[ixcv] 
    y_env$yov$x_orig <- data$x_orig[ixcv]
    y_env$yov$y_orig <- data$y_orig[ixcv]
    y_env$yov$lat    <- data$lat[ixcv]
    y_env$yov$lon    <- data$lon[ixcv]
    y_env$yov$z      <- data$z[ixcv]
    y_env$yov$souid  <- data$sourceId[ixcv]
    y_env$yov$value  <- data$value[ixcv]
    if (file.exists(argv$iff_laf)) y_env$yov$laf <-extract( env$rlaf, cbind( y_env$yov$x, y_env$yov$y), na.rm=T)
    y_env$yov$prid   <- data$prId[ixcv]
    if ( !is.na( argv$rrinf)) {
      y_env$yov$ixwet <- which( y_env$yov$value >= argv$rrinf)
      y_env$yov$ixdry <- which( y_env$yov$value <  argv$rrinf)
      y_env$yov$nwet  <- length( y_env$yov$ixwet)
      y_env$yov$ndry  <- length( y_env$yov$ixdry)
    }
    y_env$yov$nuo  <- u_env$nuo 
    if (u_env$nuo > 0) {
      y_env$yov$value_uo  <- array( data=NA, dim=c( y_env$yov$n, u_env$nuo))
      for (i in 1:u_env$nuo) y_env$yov$value_uo[,i]  <- data$value_uo[ixcv,i]
    }
    if (fg_env$ktot_dim > 0) {
      y_env$yov$value_fg  <- array( data=NA, dim=c( y_env$yov$n, fg_env$ktot_dim))
      for (i in 1:fg_env$ktot_dim) y_env$yov$value_fg[,i]  <- data$value_fg[ixcv,i]
    }

  if (!file.exists(argv$iff_laf)) y_env$yov$laf <- rep(NA,y_env$yov$n)
  if (u_env$nuo == 0) y_env$yov$value_uo <- rep(NA,y_env$yov$value_uo)

  names_out <- c("date","x","y","x_orig","y_orig","lat","lon","z","souid","laf","prid","value","value_uo")
  out <- data.frame( rep(argv$date_out, y_env$yov$n),
                     y_env$yov$x,
                     y_env$yov$y,
                     y_env$yov$x_orig,
                     y_env$yov$y_orig,
                     y_env$yov$lat,
                     y_env$yov$lon,
                     y_env$yov$z,
                     y_env$yov$souid,
                     y_env$yov$laf,
                     y_env$yov$prid,
                     y_env$yov$value,
                     y_env$yov$value_uo)
  for (i in 1:fg_env$ktot_dim) {
    
    out <- cbind( out, y_env$yov$value_fg[,i])
  }

  cat("date;sourceId;x;y;z;yo_cv;yb_cv;ya_cv;yidi_cv;dqc;\n",
      file=argv$off_cv_table,append=F)
  cat(paste(argv$date_out,
            formatC(VecS_cv,format="f",digits=0),
            formatC(VecX_cv,format="f",digits=0),
            formatC(VecY_cv,format="f",digits=0),
            formatC(VecZ_cv,format="f",digits=0),
            formatC(yo_cv,format="f",digits=2),
            formatC(yb_cv,format="f",digits=2),
            formatC(ya_cv,format="f",digits=2),
            formatC(yidi_cv,format="f",digits=4),
            rep(0,length(VecS_cv)),
            "\n",sep=";"),
      file=argv$off_cv_table,append=T)
  print(paste("output saved on file",argv$off_cv_table))
}
