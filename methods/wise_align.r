#+ Select the background fields using an observed field as the reference
#wise_align <- function( argv, fg_env, u_env, env, plot=F, dir_plot=NA) {
wise_align <- function( argv, fg_env, env, plot=F, dir_plot=NA) {
#
# output
#  fg_env$score score for each potential background field (higher the better)
#  fg_env$ixf data source for each potential background field
#  fg_env$ixe ensemble member number for each potential background field
#  fg_env$ixs selection of the k background fields (indexes wrt fg_env$score)
#
#------------------------------------------------------------------------------
  options( warn = 2)

  cat( "-- main_wise_align --\n")

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
      if (argv$wise_align_mode == "ets") {
        a_random <- as.numeric( (a+c)*(a+b) / (a+b+c+d))
        score <- c( score, (a-a_random) / (a+c+b-a_random))
      } else if (argv$wise_align_mode == "maxoverlap") {
        score <- c( score, a/length(r))
      }
      ixf <- c( ixf, f)
      ixe <- c( ixe, e)
      j <- j+1
      cat( ".")

      if (plot) {
        fout<-file.path(dir_plot,paste0("wise_align_f",f,"_e",formatC(e,width=2,flag="0"),".png"))
        f1 <- file.path(dir_plot,"f1.png")
        f2 <- file.path(dir_plot,"f2.png")

        proj4.lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
        b_utm33<-readOGR("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_UTM33/TM_WORLD_BORDERS_UTM33-0.2.shp","TM_WORLD_BORDERS_UTM33-0.2",verbose=F)
        b<-spTransform(b_utm33,CRS(proj4.lcc))
        rm(b_utm33,proj4.lcc)

#        ra <- u_env$uo[[1]]$r_main
#        rb <- u_env$uo[[1]]$r_main
        ra <- env$rmaster
        rb <- env$rmaster
        png(file=f1,width=800,height=800)
        ra[]<-uo
        image(ra,breaks=c(-1,0.5,1.5),col=c("gray","cornflowerblue"))
        plot(b,add=T)
        dev.off()
        png(file=f2,width=800,height=800)
        rb[]<-r
        rb<-mask(rb,ra)
        image(rb,breaks=c(-1,0.5,1.5),col=c("gray","cornflowerblue"),main=paste("score ",round(score[j],4)))
        plot(b,add=T)
        dev.off()
        system( paste0("convert +append ",f1," ",f2," ",fout))
        system( paste0("rm ",f1," ",f2))
        rm(ra,rb)
      }

    }
    cat( "\n")
  } 

  cat( paste(" number of background ensemble members (tot) =", env$k_dim,"(",length(score),")\n"))
#  # select only the k fields with the highest scores
  ixs <- order( score, decreasing=T)[1:env$k_dim]

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

  # diagnostics
  if (plot) {
    cat( "diagnostics ")
    j <- 0
    for (i in ixs) {
      j <- j + 1
      cat(".")
      r <- subset( fg_env$fg[[ixf[i]]]$r_main, subset=ixe[i])
      fout<-file.path(dir_plot,paste0("wise_align_sel_e",formatC(j,width=2,flag="0"),".png"))
#      ra <- u_env$uo[[1]]$r_main
      ra <- env$rmaster
      ra[] <- env$mergeobs$value
      png(file=f1,width=800,height=800)
      image(ra,breaks=c(0,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(7))))
      plot(b,add=T)
      dev.off()
      png(file=f2,width=800,height=800)
      rb<-mask(r,ra)
      image(rb,breaks=c(0,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(7))),
            main=paste0("data source=",ixf[i]," ensemble=",ixe[i]))
      plot(b,add=T)
      dev.off()
      system( paste0("convert +append ",f1," ",f2," ",fout))
      system( paste0("rm ",f1," ",f2))
      rm(ra,rb)
    }
    cat( "\n")
  }

  return( TRUE)

}
