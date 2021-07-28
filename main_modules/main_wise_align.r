#+ Select the background fields using an observed field as the reference
main_wise_align <- function( argv, fg_env, u_env, env,
                             plot=F, dir_plot=NA) {
#
# output
#  fg_env$ets ETS for each potential background field
#  fg_env$ixf data source for each potential background field
#  fg_env$ixe ensemble member number for each potential background field
#  fg_env$ixs selection of the k background fields (indexes wrt fg_env$ets)
#
#------------------------------------------------------------------------------
  options( warn = 2)

  cat( "-- main_wise_align --\n")

  if ( ( nfg <- length( fg_env$fg)) == 0) return( FALSE)

  # compute ETS for each potential background field
  cat (" compute ETS for each potential background field \n")
  # reference observed field
  uo <- getValues( u_env$uo[[1]]$r_main)
  uo[uo<u_env$rain] <- 0
  uo[uo>u_env$rain] <- 1
  ets <- numeric(0)
  ixf <- integer(0)
  ixe <- integer(0)
  j <- 0
  for (f in 1:nfg) {
    cat( paste("  data source",f,"ensembles"))
    ef <- nlayers( fg_env$fg[[f]]$r_main)
    for (e in 1:ef) {
      r <- getValues( subset( fg_env$fg[[f]]$r_main, subset=e))
      r[r<u_env$rain] <- 0
      r[r>u_env$rain] <- 1
      a <- as.numeric( length( which(r==1 & uo==1)))
      b <- as.numeric( length( which(r==1 & uo==0)))
      c <- as.numeric( length( which(r==0 & uo==1)))
      d <- as.numeric( length( which(r==0 & uo==0)))
      a_random <- as.numeric( (a+c)*(a+b) / (a+b+c+d))
      ets <- c( ets, (a-a_random) / (a+c+b-a_random))
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

        ra <- u_env$uo[[1]]$r_main
        rb <- u_env$uo[[1]]$r_main
        png(file=f1,width=800,height=800)
        ra[]<-uo
        image(ra,breaks=c(-1,0.5,1.5),col=c("gray","cornflowerblue"))
        plot(b,add=T)
        dev.off()
        png(file=f2,width=800,height=800)
        rb[]<-r
        rb<-mask(rb,ra)
        image(rb,breaks=c(-1,0.5,1.5),col=c("gray","cornflowerblue"),main=paste("ets ",round(ets[j],4)))
        plot(b,add=T)
        dev.off()
        system( paste0("convert +append ",f1," ",f2," ",fout))
        system( paste0("rm ",f1," ",f2))
        rm(ra,rb)
      }

    }
    cat( "\n")
  } 
  
  # select only the k fileds with the highest ETSs
  ixs <- order( ets, decreasing=T)[1:env$k_dim]

  # output
  fg_env$ets <- ets
  fg_env$ixf <- ixf
  fg_env$ixe <- ixe
  fg_env$ixs <- ixs

  # diagnostics
  if (plot) {
    j <- 0
    for (i in ixs) {
      j <- j + 1
      print(j)
      r <- subset( fg_env$fg[[ixf[i]]]$r_main, subset=ixe[i])
      fout<-file.path(dir_plot,paste0("wise_align_sel_e",formatC(j,width=2,flag="0"),".png"))
      ra <- u_env$uo[[1]]$r_main
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
  }

  return( TRUE)

}
