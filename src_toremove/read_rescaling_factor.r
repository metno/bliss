#+ Read rescaling factor
read_rescaling_factor <- function( argv, env) {
#=============================================================================

  if (file.exists(argv$iff_rf)) {

    if (argv$verbose) {
      cat("+---------------------------------------------------------------+\n")
      cat(paste("+ rescaling factor",argv$iff_rf,"\n"))
    }
    res <- read_and_regrid_nc( 
      nc.file    = argv$iff_rf,
      nc.varname = argv$iff_rf.varname,
      topdown    = argv$iff_rf.topdown,
      out.dim    = list( ndim=argv$iff_rf.ndim,
                         tpos=argv$iff_rf.tpos,
                         epos=argv$iff_rf.epos,
                         names=argv$iff_rf.names),
      proj4      = argv$iff_rf.proj4,
      nc.proj4   = list( var=NULL,
                         att=NULL),
      selection  = list( t=argv$iff_rf.t,
                         format=argv$iff_rf.tfmt,
                         e=argv$iff_rf.e),
      adjfact    = argv$iff_rf.adjfact,
      adjval     = argv$iff_rf.adjval,
      rmaster    = rmaster,
      grid_master.proj4 = argv$grid_master.proj4 ,
      nc.varname_lat = argv$iff_rf.varname_lat,
      nc.varname_lon = argv$iff_rf.varname_lon,
      out.dim_ll     = list(ndim=argv$iff_rf.ndim_ll,
                            tpos=argv$iff_rf.tpos_ll,
                            epos=NULL,
                            names=argv$iff_rf.names_ll)) 
    rrf <- res$raster
    rf  <- res$values
    rm( res)

    env$rrf_exists <- TRUE
    env$rrf <- rrf

  } else {

    env$rrf_exists <- FALSE

  }
}
