  if (argv$verbose) {
    cat("+---------------------------------------------------------------+\n")
    cat(paste("+ first guess ",argv$iff_fg,"\n"))
  }
  iff_is_ens<-F
  # First-guess is not an ensemble
  if (is.null(argv$iff_fg.epos)) {
    res <- read_and_regrid_nc( 
      nc.file    = argv$iff_fg,
      nc.varname = argv$iff_fg.varname,
      topdown    = argv$iff_fg.topdown,
      out.dim    = list( ndim=argv$iff_fg.ndim,
                         tpos=argv$iff_fg.tpos,
                         epos=argv$iff_fg.epos,
                         names=argv$iff_fg.names),
      proj4      = argv$iff_fg.proj4,
      nc.proj4   = list( var=NULL,
                         att=NULL),
      selection  = list( t=argv$iff_fg.t,
                         format=argv$iff_fg.tfmt,
                         e=NULL),
      adjfact    = argv$iff_fg.adjfact,
      adjval     = argv$iff_fg.adjval,
      rmaster    = rmaster,
      grid_master.proj4 = argv$grid_master.proj4 ,
      nc.varname_lat = "none",
      nc.varname_lon = "none",
      out.dim_ll     = list(ndim=argv$iff_fg.ndim_ll,
                            tpos=argv$iff_fg.tpos_ll,
                            epos=NULL,
                            names=argv$iff_fg.names_ll),
      upscale = argv$fg_upscale) 
    rfg <- mask( res$raster, rmaster)
    rm( res)
    xb_tmp <- getValues(rfg)
    if ( !is.na( argv$min_plaus_val)) 
      xb_tmp[ which( xb_tmp < argv$min_plaus_val)] <- argv$min_plaus_val
    rfg[] <- xb_tmp
    aix <- which(!is.na(xb_tmp))
    xb  <- xb_tmp[aix]
    rm( xb_tmp)
  # First-guess is an ensemble
  } else {
    iff_is_ens<-T
    nens<-0 # number of ensemble members having at least one value different from NA
    for (e in 1:length(argv$iff_fg.e)) {
      res <- read_and_regrid_nc( 
        nc.file    = argv$iff_fg,
        nc.varname = argv$iff_fg.varname,
        topdown    = argv$iff_fg.topdown,
        out.dim    = list( ndim=argv$iff_fg.ndim,
                           tpos=argv$iff_fg.tpos,
                           epos=argv$iff_fg.epos,
                           names=argv$iff_fg.names),
        proj4      = argv$iff_fg.proj4,
        nc.proj4   = list( var=NULL,
                           att=NULL),
        selection  = list( t=argv$iff_fg.t,
                           format=argv$iff_fg.tfmt,
                           e=argv$iff_fg.e[e]),
        adjfact    = argv$iff_fg.adjfact,
        adjval     = argv$iff_fg.adjval,
        rmaster    = rmaster,
        grid_master.proj4 = argv$grid_master.proj4 ,
        nc.varname_lat = "none",
        nc.varname_lon = "none",
        out.dim_ll     = list(ndim=argv$iff_fg.ndim_ll,
                              tpos=argv$iff_fg.tpos_ll,
                              epos=NULL,
                              names=argv$iff_fg.names_ll),
        upscale = argv$fg_upscale) 
      rfg <- mask( res$raster, rmaster)
      rm( res)
      xb_tmp <- getValues(rfg)
      aix <- which(!is.na(xb_tmp))
      if ( length( aix)>0) {
        if ( !exists("xb_ens_tmp")) xb_ens_tmp <- array(data=NA,
                                          dim=c( length(aix),length(argv$iff_fg.e)))
        if ( length(aix) != dim(xb_ens_tmp)[1]) boom("ERROR while reading the background file")
        xb_ens_tmp[,e] <- xb_tmp[aix]
        nens <- nens+1
        if (!exists("valens")) valens <- integer(0)
        valens <- c( valens, argv$iff_fg.e[e])
      }
      rm(xb_tmp)
    }
    if ( nens==0) boom("ERROR while reading the background file")
    xb <- array( data=NA, dim=c(dim(xb_ens_tmp)[1], nens))
    e  <- 0
    for (i in 1:length(argv$iff_fg.e)) {
      if ( any( !is.na( xb_ens_tmp[,i]))) {
        e <- e+1
        xb[,e] <- xb_ens_tmp[,i]
      }
    }
    rm( xb_ens_tmp)
    if ( !is.na( argv$min_plaus_val)) 
      xb[ which( xb < argv$min_plaus_val)] <- argv$min_plaus_val
    if ( nens == 1) { xb <- drop(xb); iff_is_ens<-F }
  }
  if (argv$verbose) {
    cat(paste("# grid points (not NAs)=",length(aix),"\n"))
    cat("+...............................................................+\n")
  }

