#+ Read Land-Area Fraction
read_dem <- function( argv, env) {
#==============================================================================

  if ( file.exists( argv$iff_dem)) {

    if (argv$verbose) {
      cat("+---------------------------------------------------------------+\n")
      cat(paste("+ digital elevation model",argv$iff_dem,"\n"))
    }
    res <- read_and_regrid_nc( 
      nc.file    = argv$iff_dem,
      nc.varname = argv$iff_dem.varname,
      topdown    = argv$iff_dem.topdown,
      out.dim    = list( ndim=argv$iff_dem.ndim,
                         tpos=argv$iff_dem.tpos,
                         epos=argv$iff_dem.epos,
                         names=argv$iff_dem.names),
      proj4      = argv$iff_dem.proj4,
      nc.proj4   = list( var=NULL,
                         att=NULL),
      selection  = list( t=argv$iff_dem.t,
                         format=argv$iff_dem.tfmt,
                         e=argv$iff_dem.e),
      adjfact    = argv$iff_dem.adjfact,
      adjval     = argv$iff_dem.adjval,
      rmaster    = rmaster,
      grid_master.proj4 = argv$grid_master.proj4 ,
      nc.varname_lat = "none",
      nc.varname_lon = "none",
      out.dim_ll     = list(ndim=argv$iff_dem.ndim_ll,
                            tpos=argv$iff_dem.tpos_ll,
                            epos=NULL,
                            names=argv$iff_dem.names_ll)) 
    rdem <- res$raster
    dem  <- res$values
    rm( res)
    env$rdem_exists <- TRUE
    env$rdem <- rdem

  } else {

    env$rdem_exists <- FALSE
  }

}
