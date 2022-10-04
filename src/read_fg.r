#+ Read first-guess fields 
read_fg <- function( argv, fg_env, u_env, env) {
#==============================================================================

  t0a<-Sys.time()

  if ( ( fg_env$nfg <- length( fg_env$fg)) == 0) return( FALSE)

  if (argv$verbose) cat("+----------------------------------------------------+\n")

  #
  #----------------------------------------------------------------------------
  # u_env$uo configuration file to read observations (not in-situ)
  if ( ( u_env$nuo <- length( u_env$uo)) == 1) {

    if (argv$verbose) cat("Read observations for alignment ...")
 
    # initializations

    u_env$uo[[1]]$r_main <- NULL

    if ( is.null( u_env$uo[[1]]$main.proj4)) 
      u_env$uo[[1]]$main.proj4 <- ""
    if ( is.null( u_env$uo[[1]]$main.proj4_var)) 
      u_env$uo[[1]]$main.proj4_var <- ""
    if ( is.null( u_env$uo[[1]]$main.proj4_att)) 
      u_env$uo[[1]]$main.proj4_att <- ""
    if ( is.null( u_env$uo[[1]]$main.t)) 
      u_env$uo[[1]]$main.t <- NA
    if ( is.null( u_env$uo[[1]]$main.acc)) 
      u_env$uo[[1]]$main.acc <- F

    # read the main file
    if ( !is.null( u_env$uo[[1]]$main.file)) {
      if ( file.exists( u_env$uo[[1]]$main.file)) {
        first <- T
        if ( is.null( u_env$uo[[1]]$main.epos)) {
          ei <- 0
        } else {
          if ( is.null( u_env$uo[[1]]$main.e)) {
            ei <- nc4.getDim( u_env$uo[[1]]$main.file, 
                              varid = u_env$uo[[1]]$main.dimnames[u_env$uo[[1]]$main.epos])
          } else {
            ei <- u_env$uo[[1]]$main.e
          }
        }
        if ( is.null( u_env$uo[[1]]$main.t) | is.na(u_env$uo[[1]]$main.t)) {
            ti <- nc4.getTime( u_env$uo[[1]]$main.file)
        } else {
            ti <- u_env$uo[[1]]$main.t
        }
        for (ens in 1:length(ei)) {
          if( ei[ens] == 0) { nc_e <- NA} else { nc_e <- ei[ens]}
          res <- read_and_regrid_nc( nc.file    = u_env$uo[[1]]$main.file,
                                     nc.varname = u_env$uo[[1]]$main.varname,
                                     topdown    = u_env$uo[[1]]$main.topdown,
                                     out.dim    = list( ndim  = u_env$uo[[1]]$main.ndim,
                                                        tpos  = u_env$uo[[1]]$main.tpos,
                                                        epos  = u_env$uo[[1]]$main.epos,
                                                        names = u_env$uo[[1]]$main.dimnames),
                                     proj4      = u_env$uo[[1]]$main.proj4,
                                     nc.proj4   = list( var = u_env$uo[[1]]$main.proj4_var, 
                                                        att = u_env$uo[[1]]$main.proj4_att),
                                     selection  = list( t      = ti,
                                                        format = "%Y%m%d%H%M",
                                                        e      = nc_e),
                                     adjfact=u_env$uo[[1]]$main.cfact,
                                     adjval=u_env$uo[[1]]$main.offset,
                                     rmaster = env$rmaster,
                                     grid_master.proj4 = as.character( env$rmaster),
                                     nc.varname_lat ="none",
                                     nc.varname_lon ="none",
                                     out.dim_ll = list( ndim  = NA,
                                                        tpos  = NA,
                                                        epos  = NA,
                                                        names = NA),
                                     upscale = F,
                                     upscale_fun = "mean",
                                     projectraster_method = "bilinear") 
#b<-readOGR("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_LATLON/TM_WORLD_BORDERS-0.2.shp","TM_WORLD_BORDERS-0.2")
#proj4ll<-"+proj=longlat +datum=WGS84"
#proj4lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
#aux<-which(b@data$LON>=-10 & b@data$LON<=90 & b@data$LAT>=0 & b@data$LAT<=85)
##145<-Antarctica
#subset<-b[aux,]
##ETRS89 / ETRS-LAEA
##EPSG:32633 +proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs 
## Russia not included
#blcc<-spTransform(subset, CRS(proj4lcc))
#png(file="test.png", height=1200, width=1200)
#plot(res$raster)
#plot(blcc,add=T)
#dev.off()
#q()
          if ( !is.null( res)) {
            if ( first) {
              u_env$uo[[1]]$r_main <- u_env$uo[[1]]$main.offset + res$raster * u_env$uo[[1]]$main.cfact 
            } else {
              u_env$uo[[1]]$r_main <- stack( u_env$uo[[1]]$r_main,
                                             u_env$uo[[1]]$main.offset + res$raster * u_env$uo[[1]]$main.cfact)
            }
            first <- F
          }
          rm(res)
        } # end loop over ensemble 
      } # end if main file exists 
    } # end read the main file

    if (argv$verbose) cat("done!\n")

  }

  #
  #----------------------------------------------------------------------------
  if (argv$verbose) cat("Read background\n")

  for (f in 1:fg_env$nfg) {
    if (argv$verbose) cat(paste("file",f,"ensembles "))

    # initializations

    fg_env$fg[[f]]$r_main <- NULL

    if ( is.null( fg_env$fg[[f]]$main.proj4)) 
      fg_env$fg[[f]]$main.proj4 <- ""
    if ( is.null( fg_env$fg[[f]]$main.proj4_var)) 
      fg_env$fg[[f]]$main.proj4_var <- ""
    if ( is.null( fg_env$fg[[f]]$main.proj4_att)) 
      fg_env$fg[[f]]$main.proj4_att <- ""
    if ( is.null( fg_env$fg[[f]]$main.t)) 
      fg_env$fg[[f]]$main.t <- NA
    if ( is.null( fg_env$fg[[f]]$main.acc)) 
      fg_env$fg[[f]]$main.acc <- F

    # read the main file
    if ( !is.null( fg_env$fg[[f]]$main.file)) {
      if ( file.exists( fg_env$fg[[f]]$main.file)) {
        first <- T
        if ( is.null( fg_env$fg[[f]]$main.epos)) {
          ei <- 0
        } else {
          if ( is.null( fg_env$fg[[f]]$main.e)) {
            ei <- nc4.getDim( fg_env$fg[[f]]$main.file, 
                              varid = fg_env$fg[[f]]$main.dimnames[fg_env$fg[[f]]$main.epos])
          } else {
            ei <- fg_env$fg[[f]]$main.e
          }
        }
        if ( is.null( fg_env$fg[[f]]$main.t) | is.na(fg_env$fg[[f]]$main.t)) {
            ti <- nc4.getTime( fg_env$fg[[f]]$main.file)
        } else {
            ti <- fg_env$fg[[f]]$main.t
        }
        for (ens in 1:length(ei)) {
          if( ei[ens] == 0) { nc_e <- NA} else { nc_e <- ei[ens]}
          if (argv$verbose) cat(".")
          res <- read_and_regrid_nc( nc.file    = fg_env$fg[[f]]$main.file,
                                     nc.varname = fg_env$fg[[f]]$main.varname,
                                     topdown    = fg_env$fg[[f]]$main.topdown,
                                     out.dim    = list( ndim  = fg_env$fg[[f]]$main.ndim,
                                                        tpos  = fg_env$fg[[f]]$main.tpos,
                                                        epos  = fg_env$fg[[f]]$main.epos,
                                                        names = fg_env$fg[[f]]$main.dimnames),
                                     proj4      = fg_env$fg[[f]]$main.proj4,
                                     nc.proj4   = list( var = fg_env$fg[[f]]$main.proj4_var, 
                                                        att = fg_env$fg[[f]]$main.proj4_att),
                                     selection  = list( t      = ti,
                                                        format = "%Y%m%d%H%M",
                                                        e      = nc_e),
                                     adjfact=fg_env$fg[[f]]$main.cfact,
                                     adjval=fg_env$fg[[f]]$main.offset,
                                     rmaster = env$rmaster,
                                     grid_master.proj4 = as.character( env$rmaster),
                                     nc.varname_lat ="none",
                                     nc.varname_lon ="none",
                                     out.dim_ll = list( ndim  = NA,
                                                        tpos  = NA,
                                                        epos  = NA,
                                                        names = NA),
                                     upscale = F,
                                     upscale_fun = "mean",
                                     projectraster_method = "bilinear") 
#b<-readOGR("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_LATLON/TM_WORLD_BORDERS-0.2.shp","TM_WORLD_BORDERS-0.2")
#proj4ll<-"+proj=longlat +datum=WGS84"
#proj4lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
#aux<-which(b@data$LON>=-10 & b@data$LON<=90 & b@data$LAT>=0 & b@data$LAT<=85)
##145<-Antarctica
#subset<-b[aux,]
##ETRS89 / ETRS-LAEA
##EPSG:32633 +proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs 
## Russia not included
#blcc<-spTransform(subset, CRS(proj4lcc))
#png(file="test.png", height=1200, width=1200)
#plot(res$raster)
#plot(blcc,add=T)
#dev.off()
#q()
          if ( !is.null( res)) {
            if ( first) {
              fg_env$fg[[f]]$r_main <- fg_env$fg[[f]]$main.offset + res$raster * fg_env$fg[[f]]$main.cfact 
            } else {
              fg_env$fg[[f]]$r_main <- stack( fg_env$fg[[f]]$r_main,
                                              fg_env$fg[[f]]$main.offset + res$raster * fg_env$fg[[f]]$main.cfact)
            }
            first <- F
          }
          rm(res)
        } # end loop over ensemble
        if (argv$verbose) cat("done!\n")
        fg_env$fg[[f]]$k_dim <- nlayers(fg_env$fg[[f]]$r_main)
      } # end if main file exists 
    } # end read the main file
  } # end loop over fg files
  
  # total number of background fields
  fg_env$ktot_dim <- 0
  for (f in 1:fg_env$nfg)
    fg_env$ktot_dim <- fg_env$ktot_dim + fg_env$fg[[f]]$k_dim 
 
  #----------------------------------------------------------------------------
  # print info
  if (argv$verbose) {
    t1a<-Sys.time()

    for (f in 1:fg_env$nfg) {
      if (!is.null( fg_env$fg[[f]]$r_main)) {
        for (ens in 1:nlayers(fg_env$fg[[f]]$r_main)) {
          if ( class( fg_env$fg[[f]]$r_main) == "RasterLayer") {
            r <- fg_env$fg[[f]]$r_main
          } else {
            r <- raster( fg_env$fg[[f]]$r_main, ens)
          }
          if (ens == 1) {
            cat( paste( "FG", formatC( f, width=2, flag="0"),
                        "ensemble member", formatC( ens, width=2, flag="0"), ".",
                        "dim =", ncell( r), ".", 
                        "range =", round( min( getValues( r), na.rm=T), 2), 
                                   round( max( getValues( r), na.rm=T), 2), "\n"))
          } else {
            cat( paste( "      ensemble member", formatC( ens, width=2, flag="0"), ".",
                        "dim =", ncell( r), ".", 
                        "range =", round( min( getValues( r), na.rm=T), 2), 
                                   round( max( getValues( r), na.rm=T), 2), "\n"))
          }
        }
      }
      if (!is.null( fg_env$fg[[f]]$r_aux1)) {
        r <- fg_env$fg[[f]]$r_aux1
        cat( paste( "                   aux 1 .",
                    "dim =", ncell( r), ".", 
                    "range =", round( min( getValues( r), na.rm=T), 2), 
                               round( max( getValues( r), na.rm=T), 2), "\n"))
      }
    }
    cat( paste( "total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))
    cat( "+---------------------------------+\n")
  }

  return( TRUE)
}
