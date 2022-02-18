#+ Read input observations
read_obs <- function( argv, env, y_env) {
#
#==============================================================================

  #
  # Read data from text file
  dat<-read.table( file             = argv$iff_obs,
                   header           = T,
                   sep              = argv$iff_obs.sep,
                   stringsAsFactors = F,
                   strip.white      = T)

  #
  # read x_orig,y_orig,value & set x,y 
  varidxtmp <- match( c( argv$iff_obs.x,
                         argv$iff_obs.y,
                         argv$iff_obs.value),
                         names(dat))
  if (any(is.na(varidxtmp))) {
    print("ERROR in the specification of the variable names")
    print(paste("        x=",argv$iff_obs.x))
    print(paste("        y=",argv$iff_obs.y))
    print(paste("    value=",argv$iff_obs.value))
    print("header of input file:")
    print(argv$argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }

  #
  # x_orig, y_orig and value are the info we always need
  data <- data.frame( dat[,varidxtmp])
  names(data) <- c( "x_orig", "y_orig", "value")
  data$x_orig <- suppressWarnings( as.numeric( data$x_orig))
  data$y_orig <- suppressWarnings( as.numeric( data$y_orig))
  data$value  <- suppressWarnings( as.numeric( data$value))
  ndata <- length( data$x_orig)

  #
  # re-project (x,y) if needed  
  if (argv$iff_obs.proj4!=argv$grid_master.proj4) {
    xymaster <- spTransform( SpatialPoints( cbind( data$x_orig, data$y_orig),
                                            proj4string = CRS(argv$iff_obs.proj4)),
                          CRS(argv$grid_master.proj4))
    data$x <- attr( xymaster, "coords")[,1]
    data$y <- attr( xymaster, "coords")[,2]
    rm(xymaster)
  } else {
    data$x <- data$x_orig
    data$y <- data$y_orig
  }

  #
  # lat-lon coords are required by verif
  if (argv$iff_obs.proj4!=env$proj4.llwgs84) {
    xyll <- spTransform( SpatialPoints( cbind( data$x_orig, data$y_orig),
                                    proj4string=CRS(argv$iff_obs.proj4)) ,
                      CRS(env$proj4.llwgs84))
    data$lon<-attr(xyll,"coords")[,1]
    data$lat<-attr(xyll,"coords")[,2]
    rm(xyll)
  } else {
    data$lon<-data$x_orig
    data$lat<-data$y_orig
  }
  
  #
  # read z (vertical coordinate / altitude / elevation)
  if (argv$iff_obs.z!="none") {
    varidxtmp <- which( argv$iff_obs.z == names(dat))
    if ( length(varidxtmp) == 0) {
      print("ERROR in the specification of the variable names")
      print(paste("     z=",argv$iff_obs.z))
      print("header of input file:")
      print(argv$iff_obs)
      print(names(dat))
      quit(status=1)
    }
    data$z <- dat[,varidxtmp]
  } else {
    data$z <- rep(0,length(data$x))
  }

  #
  # read sourceId 
  if (argv$iff_obs.sourceId!="none") {
    varidxtmp<-which(argv$iff_obs.sourceId==names(dat))
    if (length(varidxtmp)==0) {
      print("ERROR in the specification of the variable names")
      print(paste("     sourceId=",argv$iff_obs.sourceId))
      print("header of input file:")
      print(argv$iff_obs)
      print(names(dat))
      quit(status=1)
    }
    data$sourceId<-dat[,varidxtmp]
  } else {
    data$sourceId<-1:length(data$x)
  }

  #
  # read prId 
  if (argv$iff_obs.prId!="none") {
    varidxtmp<-which(argv$iff_obs.prId==names(dat))
    if (length(varidxtmp)==0) {
      print("ERROR in the specification of the variable names")
      print(paste("     prId=",argv$iff_obs.prId))
      print("header of input file:")
      print(argv$iff_obs)
      print(names(dat))
      quit(status=1)
    }
    data$prId<-dat[,varidxtmp]
  } else {
    data$prId<-rep(0,length(data$x))
  }
  
  #
  # read dqc flag 
  if (argv$iff_obs.dqc!="none") {
    varidxtmp<-which(argv$iff_obs.dqc==names(dat))
    if (length(varidxtmp)==0) {
      print("ERROR in the specification of the variable names")
      print(paste("     dqc=",argv$iff_obs.dqc))
      print("header of input file:")
      print(argv$iff_obs)
      print(names(dat))
      quit(status=1)
    }
    data$dqc<-dat[,varidxtmp]
  } else {
    data$dqc<-rep(0,length(data$x))
  }

  # so far:
  #   data, dataframe: x,y,z,x_orig,y_orig,value,sourceId,prId,dqc

  #
  # blacklist
  if (file.exists(argv$iff_black)) {
    bstid <- read.csv( argv$iff_black, header=T, stringsAsFactors=F,strip.white=T)
  } else {
    bstid <- integer(0)
  }

  #
  # extract aux info from u_env
  flag_in_uo <- rep( T, ndata)
  if ( u_env$nuo > 0) {
    k <- 0 
    for (i in 1:u_env$nuo) {
      for (j in 1:nlayers(u_env$uo[[i]]$r_main)) {
        if ( class(u_env$uo[[i]]$r_main)[1] == "RasterStack") {
          r <- raster(u_env$uo[[i]]$r_main, layer=j)
        } else if ( class(u_env$uo[[i]]$r_main)[1] == "RasterLayer") {
          r <- u_env$uo[[i]]$r_main
        } else {
          return(FALSE)
        }
        if ( k == 0) {
          data$value_uo <- array( data=NA, dim=c(ndata, 1))
          data$value_uo[,1] <- extract( r, cbind( data$x, data$y)) }
        else {
          data$value_uo <- cbind(  data$value_uo, extract( r, cbind( data$x, data$y)))
        }
        k <- k + 1
        flag_in_uo <- flag_in_uo & !is.na(data$value_uo[,k])
      }
    }
  }

  #
  # extract aux info from fg_env
  flag_in_fg <- rep( T, ndata)
  if ( fg_env$nfg > 0) {
    data$value_fg <- array( data=NA, dim=c(ndata, fg_env$ktot_dim))
    # loop over the background files
    for (i in 1:fg_env$nfg) {
      # loop over the backg fields of the j-th backg file 
      for (j in 1:fg_env$fg[[i]]$k_dim) {
        if ( class(fg_env$fg[[i]]$r_main)[1] == "RasterStack") {
          r <- raster(fg_env$fg[[i]]$r_main, layer=j)
        } else if ( class(fg_env$fg[[i]]$r_main)[1] == "RasterLayer") {
          r <- fg_env$fg[[i]]$r_main
        } else {
          return(FALSE)
        }
        data$value_fg[,j] <- extract( r, cbind( data$x, data$y))
        flag_in_fg <- flag_in_fg & !is.na(data$value_fg[,j])
      } # end loop over the backg fields of the j-th backg file
    } # end loop over the background files
  }

  #
  # select only observation within the master grid
  flag_in_master <- !is.na( extract( env$rmaster, cbind( data$x, data$y)))

  # on-the-fly dqc, used for testing
  #  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) &
  #              data$value > (extract(rfg,cbind(data$x,data$y))-0.2*extract(rfg,cbind(data$x,data$y)))

  #
  # CVmode based on data provider
  if (env$cv_mode | env$cv_mode_random) {

    auxflag <- rep( T, ndata)
    if ( any( !is.na( argv$prId.exclude)))
      auxflag <- !(data$prId %in% argv$prId.exclude)  

    if (env$cv_mode) {
      # prId=1 MET-WMO stations
      ixcv <- which( auxflag &
                     data$dqc == 0               & 
                     data$prId %in% argv$prId.cv & 
                     !is.na( data$value)         &
                     flag_in_master              &
                     flag_in_fg )
    # CVmode based random selection of stations (predefined distance between them)
    } else if (env$cv_mode_random) {
      # define the set of indices for the cv stations
      stmp <- vector()
      xtmp <- vector()
      ytmp <- vector()
      j  <- 0
      if (!is.na(argv$cv_mode_random_setseed)) set.seed(argv$cv_mode_random_setseed)
      randomize <- sample( ndata, replace=FALSE)
      for (i in randomize) {
        if ( !(auxflag[i]) | is.na(data$value[i]) | data$dqc[i]!=0 | !flag_in_master[i] | !flag_in_fg[i] ) next
        if ( j == 0) {
          j <- j + 1
          stmp[j] <- data$sourceId[i]
          xtmp[j] <- data$x[i]
          ytmp[j] <- data$y[i]
        } else {
          dtmp <- sqrt( (data$x[i]-xtmp[1:j])**2 + (data$y[i]-ytmp[1:j])**2 )
          if ( !any( dtmp < argv$d.cv)) {
            j <- j + 1
            stmp[j] <- data$sourceId[i]
            xtmp[j] <- data$x[i]
            ytmp[j] <- data$y[i]
          }
        }
      }
      ixcv <- which( data$sourceId %in% stmp[1:j])
    }

    if ( ( y_env$yov$n <- length(ixcv)) == 0) boom("ERROR \"cv_mode\" running without CV-observations ")
    y_env$yov$x      <- data$x[ixcv] 
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
    data$value[ixcv] <- NA
    data$dqc[ixcv]   <- 999
  } # end cv_mode and cv_mode_random

  #
  # Select stations to be used as input (exclude some stations if required)
  if ( any( !is.na( argv$prId.exclude))) {
    ix0 <- which( data$dqc == 0                       &  
                  !(data$prId %in% argv$prId.exclude) &
                  flag_in_master                      &
                  flag_in_fg                          &
                  !is.na(data$value)                  &
                  !is.nan(data$value) )
  } else {
    ix0 <- which( data$dqc==0        &
                  flag_in_master     &
                  flag_in_fg         &
                  !is.na(data$value) &
                  !is.nan(data$value) )
  }

  #
  # definitive station list
  if ( ( y_env$yo$n <- length(ix0)) == 0) boom("ERROR No valid input observations.")
  y_env$yo$x      <- data$x[ix0] 
  y_env$yo$y      <- data$y[ix0] 
  y_env$yo$x_orig <- data$x_orig[ix0]
  y_env$yo$y_orig <- data$y_orig[ix0]
  y_env$yo$lat    <- data$lat[ix0]
  y_env$yo$lon    <- data$lon[ix0]
  y_env$yo$z      <- data$z[ix0]
  y_env$yo$souid  <- data$sourceId[ix0]
  y_env$yo$value  <- data$value[ix0]
  if (file.exists(argv$iff_laf)) y_env$yo$laf <- extract( env$rlaf, cbind( y_env$yo$x, y_env$yo$y), na.rm=T)
  y_env$yo$prid   <- data$prId[ix0]
  y_env$yo$dqc_flag <- rep( 0, length=y_env$yo$n)
  if ( !is.na( argv$rrinf)) {
    y_env$yo$ixwet <- which( y_env$yo$value >= argv$rrinf)
    y_env$yo$ixdry <- which( y_env$yo$value <  argv$rrinf)
    y_env$yo$nwet  <- length( y_env$yo$ixwet)
    y_env$yo$ndry  <- length( y_env$yo$ixdry)
  }
  y_env$yo$nuo  <- u_env$nuo 
  if (u_env$nuo > 0) {
    y_env$yo$value_uo  <- array( data=NA, dim=c( y_env$yo$n, u_env$nuo))
    for (i in 1:u_env$nuo) y_env$yo$value_uo[,i]  <- data$value_uo[ix0,i]
  }
  y_env$yo$nfg  <- fg_env$nfg 
  if (fg_env$ktot_dim > 0) {
    y_env$yo$value_fg  <- array( data=NA, dim=c( y_env$yo$n, fg_env$ktot_dim))
    for (i in 1:fg_env$ktot_dim) y_env$yo$value_fg[,i]  <- data$value_fg[ix0,i]
  }
  rm(data)

  if (argv$cv_mode_calcidiv) {
    y_env$yov$idi <- rep( 0, y_env$yov$n)
    suppressPackageStartupMessages( library( "RANN"))
    nn2 <- nn2( cbind( y_env$yo$x, y_env$yo$y), 
                query = cbind( y_env$yov$x, y_env$yov$y),
                k = argv$calcidiv_nobs, 
                searchtype = "radius", radius = argv$calcidiv_radius)
    mat <- nn2[[1]]
    c_xy <- which( ( aux <- rowSums( mat)) > 0 )
    if ( length(c_xy) > 0) {
      dh2 <- argv$calcidiv_dh * argv$calcidiv_dh
      mapply_idi <- function(i) {
        j <- c_xy[i]
        n <- length(which(mat[j,]!=0))
        vx <- y_env$yo$x[mat[j,1:n]]
        vy <- y_env$yo$y[mat[j,1:n]]
        deltax <- abs(y_env$yov$x[j]-vx)
        deltay <- abs(y_env$yov$y[j]-vy)
        rloc<-exp( -0.5 * (deltax*deltax+deltay*deltay) / dh2 )
        S<-exp(-0.5*(outer(vx,vx,FUN="-")**2. + outer(vy,vy,FUN="-")**2)/dh2)
        SRinv<-chol2inv(chol( (S+diag(x=0.1,n)) ))
        return( sum(rloc*as.vector(rowSums(SRinv))))
      }
      y_env$yov$idi[c_xy] <- t( mapply( mapply_idi, 1:length(c_xy), SIMPLIFY = T))
    }
    t1 <- Sys.time()
    cat( paste("cv_mode_calcidiv time=", round(t1-t0,1), attr(t1-t0,"unit")))
  }

#  #
#  # Aggregate observations
#  # NOTE: we take into account X and Y dimensions, but not the Z
#  if (argv$obspp.agg) {
#
#    # temporary data structures
#    v1    <- numeric(0)
#    x1    <- numeric(0)
#    y1    <- numeric(0)
#    prid1 <- integer(0)
#
#    # loop over the provider ids
#    for (i in 1:length(argv$obspp.agg_prId)) {
#      
#      # super-observations (superobs) for the i-th provider are obtained as:
#      #  a) generate a predefined regular grid
#      #  b) for each grid box, take the average of all obs within that box
#
#      # select observations belonging to the i-th provider
#      if ( is.na( argv$obspp.agg_prId[i])) {
#        ix <- 1:y_env$yo$n
#      } else {
#        ix <- which( y_env$yo$prid == argv$obspp.agg_prId[i])
#      }
#      if ( length(ix) == 0) next
#      # generate the predef grid
#      if ( argv$obspp.agg_fact[i] > 1) {
#        r <- aggregate(rmaster, fact=argv$obspp.agg_fact[i], expand=T, na.rm=T)
#      } else {
#        r <- rmaster
#      }
#      # obtain the superobs
#      r1  <- rasterize( x     = cbind( y_env$yo$x[ix], y_env$yo$y[ix]),
#                        y     = r,
#                        field = y_env$yo$value[ix],
#                        fun   = mean,
#                        na.rm = T)
#      if ( ( n1 <- length( ix1 <- which( !is.na( getValues(r1))))) == 0) next
#      # superobs values and coordinates
#      v1 <- c( v1, getValues(r1)[ix1])
#      x1 <- c( x1, xyFromCell(r1,ix1)[,1])
#      y1 <- c( y1, xyFromCell(r1,ix1)[,2])
#      prid1 <- c( prid1, rep( argv$obspp.agg_prId[i], n1))
#    } # end loop over provider ids
#
#    # we need coordinates in the original CRS (could be different from the one of the grid)
#    if ( argv$iff_obs.proj4 != argv$grid_master.proj4) {
#      xymaster <- spTransform( SpatialPoints( cbind( x1, y1), proj4string=CRS(argv$grid_master.proj4)), CRS(argv$iff_obs.proj4))
#      x1_orig <- attr( xymaster, "coords")[,1]
#      y1_orig <- attr( xymaster, "coords")[,2]
#      rm(xymaster)
#    } else {
#      x1_orig <- x1
#      y1_orig <- y1
#    }
#    # we need coordinates in the geographical CRS
#    if ( env$proj4.llwgs84 != argv$grid_master.proj4) {
#      xyll <- spTransform( SpatialPoints( cbind( x1, y1), proj4string=CRS(argv$grid_master.proj4)), CRS(proj4.llwgs84))
#      x1_ll <- attr( xyll, "coords")[,1]
#      y1_ll <- attr( xyll, "coords")[,2]
#      rm(xyll)
#    } else {
#      x1_ll <- x1
#      y1_ll <- y1
#    }
#
#    # backup of the original obs
#    y_env$yo_backup <-  y_env$yo
#
#    #
#    # replace the definitive station list
#    if ( ( y_env$yo$n <- length(v1)) == 0) boom("ERROR No valid input observations.")
#    y_env$yo$x      <- x1 
#    y_env$yo$y      <- y1 
#    y_env$yo$x_orig <- x1_orig
#    y_env$yo$y_orig <- y1_orig 
#    y_env$yo$lat    <- x1_ll 
#    y_env$yo$lon    <- y1_ll
#    y_env$yo$z      <- rep(  0, y_env$yo$n)
#    y_env$yo$souid  <- rep( NA, y_env$yo$n)
#    y_env$yo$value  <- v1
#    y_env$yo$prid   <- prid1
#    if (file.exists(argv$iff_laf)) y_env$yo$laf <- extract( env$rlaf, cbind( y_env$yo$x, y_env$yo$y), na.rm=T)
#    y_env$yo$prid   <- data$prId[ix0]
#    y_env$yo$dqc_flag <- rep( 0, length=y_env$yo$n)
#    if ( !is.na( argv$rrinf)) {
#      y_env$yo$ixwet <- which( y_env$yo$value >= argv$rrinf)
#      y_env$yo$ixdry <- which( y_env$yo$value <  argv$rrinf)
#      y_env$yo$nwet  <- length( y_env$yo$ixwet)
#      y_env$yo$ndry  <- length( y_env$yo$ixdry)
#    }
#
#    if ( y_env$yo$nuo > 0) 
#      for (i in 1:u_env$nuo)
#        y_env$yo$value_uo[,i] <- extract( u_env$uo[[i]]$r_main, cbind( y_env$yo$x, y_env$yo$y))
#
#    if ( y_env$yo$nfg > 0) 
#      for (i in 1:fg_env$nfg)
#        y_env$yo$value_fg[,i] <- extract( fg_env$fg[[i]]$r_main, cbind( y_env$yo$x, y_env$yo$y))
#
#  } # end observation aggregation

  #
  # observation error variance correction factor
  ovarc <- rep( 1, y_env$yo$n)
  if ( any( !is.na( argv$ovarc.prId))) {
    for (i in 1:length(argv$ovarc.prId)) {
      if ( any( prId == argv$ovarc.prId[i])) 
        ovarc[which(prId==argv$ovarc.prId[i])] <- argv$ovarc[i]
    }
  }


  if (!is.na(argv$off_obspp)) {
    dataout<-array(data=NA, dim=c(length(VecLat),8))
    names(dataout)<-c("lat","lon","elev","value","prid","dqc","sct","rep")
    dataout[,1]<-round(VecLat,6)
    dataout[,2]<-round(VecLon,6)
    dataout[,3]<-round(VecZ,0)
    dataout[,4]<-round(yo,3)
    dataout[,5]<-prId
    dataout[,6]<-round(ydqc.flag,0)
    dataout[,7]<-0
    dataout[,8]<-round(eps2[-ixsus],6)
    write.table(file=argv$off_obspp,
                dataout,
                quote=F,
                col.names=c("lat","lon","elev","value","prid","dqc","sct","rep"),
                row.names=F,
                dec=".",
                sep=";")
    print(paste("output saved on file",argv$off_obspp))
    quit(status=0)
  }

  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    if (!is.na(argv$rrinf)) {
      print(paste("#observations (wet/dry) =",y_env$yo$n,"(",y_env$yo$nwet,"/",y_env$yo$ndry,")"))
      if (env$cv_mode | env$cv_mode_random) {
        print(paste("#cv-observations (wet/dry) =",y_env$yov$n,"(",y_env$yov$nwet,"/",y_env$yov$ndry,")"))
      }
    } else {
      print(paste("#observations =",y_env$yo$n))
    }
    print("+...............................................................+")
  }

  return (TRUE)

}
