#+ Read Land-Area Fraction
main_laf <- function( argv, env) {
#==============================================================================

  if ( file.exists( argv$iff_laf)) {

    if (argv$verbose) {
      cat("+---------------------------------------------------------------+\n")
      cat(paste("+ land area fraction",argv$iff_laf,"\n"))
    }
    res <- read_and_regrid_nc( 
      nc.file    = argv$iff_laf,
      nc.varname = argv$iff_laf.varname,
      topdown    = argv$iff_laf.topdown,
      out.dim    = list( ndim=argv$iff_laf.ndim,
                         tpos=argv$iff_laf.tpos,
                         epos=argv$iff_laf.epos,
                         names=argv$iff_laf.names),
      proj4      = argv$iff_laf.proj4,
      nc.proj4   = list( var=NULL,
                         att=NULL),
      selection  = list( t=argv$iff_laf.t,
                         format=argv$iff_laf.tfmt,
                         e=argv$iff_laf.e),
      adjfact    = argv$iff_laf.adjfact,
      adjval     = argv$iff_laf.adjval,
      rmaster    = rmaster,
      grid_master.proj4 = argv$grid_master.proj4 ,
      nc.varname_lat = "none",
      nc.varname_lon = "none",
      out.dim_ll     = list(ndim=argv$iff_laf.ndim_ll,
                            tpos=argv$iff_laf.tpos_ll,
                            epos=NULL,
                            names=argv$iff_laf.names_ll)) 
    rlaf <- res$raster
    laf  <- res$values
    env$rlaf_exists <- TRUE
    env$rlaf <- rlaf

  } else {

    env$rlaf_exists <- FALSE
  }
}

