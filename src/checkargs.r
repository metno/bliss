#+ Check input arguments and set some flags for elaborations
checkargs <- function( argv, env) {

#------------------------------------------------------------------------------

#
#------------------------------------------------------------------------------
# check input arguments
#

if (!(file.exists(argv$iff_obs))) boom(paste0("file not found ",argv$iff_obs))

if (argv$mode=="OI_multiscale") {
  if (argv$verbose) {
    if (!(file.exists(argv$iff_rf))) 
      print("warning: file not found",argv$iff_rf)
  }
} else if (argv$mode=="OI_firstguess") {
  if (!(file.exists(argv$iff_fg))) 
    boom(paste0("file not found ",argv$iff_fg))
} else if (argv$mode=="OI_twosteptemperature") {
  if (argv$verbose) {
    if (!(file.exists(argv$iff_dem))) 
      print("warning: file not found",argv$iff_dem)
  }
} else if ( argv$mode=="hyletkf") {
  print("chosen mode is Hybrid Local Ensemble Transform Kalman Filter")
} else if ( argv$mode=="letkf") {
  print("chosen mode is Local Ensemble Transform Kalman Filter")
} else if ( argv$mode=="OI_Bratseth") {
  print("chosen mode is OI_Bratseth")
} else if ( argv$mode=="ensigap") {
  print("chosen mode is ensigap")
} else if ( argv$mode=="wise") {
  print("chosen mode is wise")
} else if ( argv$mode=="rasterize") {
  print("chosen mode is rasterize")
} else {
  boom("error statistical interpolation scheme undefined")
}

# define/check paths and load external functions
#if ( !(file.exists(argv$path2src)) ) 
#  ext<-boom("path not found")
#
argv$iff_rf.adjfact  <- as.numeric( gsub( "_", "-", argv$iff_rf.adjfact))
argv$iff_rf.adjval   <- as.numeric( gsub( "_", "-", argv$iff_rf.adjval))
argv$iff_laf.adjfact <- as.numeric( gsub( "_", "-", argv$iff_laf.adjfact))
argv$iff_laf.adjval  <- as.numeric( gsub( "_", "-", argv$iff_laf.adjval))
argv$iff_dem.adjfact <- as.numeric( gsub( "_", "-", argv$iff_dem.adjfact))
argv$iff_dem.adjval  <- as.numeric( gsub( "_", "-", argv$iff_dem.adjval))
argv$iff_fg.adjfact  <- as.numeric( gsub( "_", "-", argv$iff_fg.adjfact))
argv$iff_fg.adjval   <- as.numeric( gsub( "_", "-", argv$iff_fg.adjval))

# parameter used in output session, select between deterministic/ensemble
argv$iff_fg.epos     <- set_NAs_to_NULL( argv$iff_fg.epos)

#
# load external C functions
#dyn.load(file.path(argv$path2src,"oi_rr_first.so"))
#dyn.load(file.path(argv$path2src,"oi_rr_fast.so"))
#dyn.load(file.path(argv$path2src,"oi_rr_var.so"))
#dyn.load(file.path(argv$path2src,"oi_t_xb_upd.so"))
#dyn.load(file.path(argv$path2src,"obsop_LapseRateConst.so"))
#if (!is.na(argv$cores)) {
#  suppressPackageStartupMessages(library("parallel"))
#  if (argv$cores==0) argv$cores <- detectCores()
#  print(paste("--> multi-core run, cores=",argv$cores))
#}
#if ( !file.exists( file.path( argv$path2src, "read_and_regrid_nc.r")))
#  boom( paste( "file not found", file.path(argv$path2src, "read_and_regrid_nc.r")))
#source( file.path( argv$path2src, "read_and_regrid_nc.r"))

#
# define elaboration modes
if (!is.na(argv$off_x)) {
  assign( "x_elab",  TRUE, envir = .GlobalEnv)
  env$x_elab <- T
} else {
  assign( "x_elab", FALSE, envir = .GlobalEnv)
  env$x_elab <- F
}

if (!is.na(argv$off_y_table)   | !is.na(argv$off_yt_table)   |
    !is.na(argv$off_y_verif_a) | !is.na(argv$off_yt_verif_a) |
    !is.na(argv$off_y_verif_b) | !is.na(argv$off_yt_verif_b) ) {
  assign( "y_elab",  TRUE, envir = .GlobalEnv)
  env$y_elab <- T
} else {
  assign( "y_elab", FALSE, envir = .GlobalEnv)
  env$y_elab <- F
}

assign( "cv_mode_random", FALSE, envir = .GlobalEnv)
env$cv_mode_random <- F

assign( "cv_mode", FALSE, envir = .GlobalEnv)
env$cv_mode <- F

if ( argv$cv_mode | 
     !is.na( argv$off_cv_table)   | !is.na( argv$off_cvt_table)   |
     !is.na( argv$off_cv_verif_a) | !is.na( argv$off_cvt_verif_a) |
     !is.na( argv$off_cv_verif_b) | !is.na( argv$off_cvt_verif_b)) {
  assign( "cv_elab", TRUE, envir = .GlobalEnv)
  env$cv_elab <- T
  if ( any( !is.na( argv$prId.cv))) {
    assign( "cv_mode", TRUE, envir = .GlobalEnv)
    env$cv_mode <- T
  } else {
    assign( "cv_mode_random", TRUE, envir = .GlobalEnv)
    env$cv_mode_random <- T
  }
} else {
  assign( "cv_elab", FALSE, envir = .GlobalEnv)
  env$cv_elab <- F
}

if (!is.na(argv$off_lcv_table)   | !is.na(argv$off_lcvt_table)   |
    !is.na(argv$off_lcv_verif_a) | !is.na(argv$off_lcvt_verif_a) |
    !is.na(argv$off_lcv_verif_b) | !is.na(argv$off_lcvt_verif_b) ) {
  assign( "loocv_elab", TRUE, envir = .GlobalEnv)
  env$loocv_elab <- T
} else {
  assign( "loocv_elab", FALSE, envir = .GlobalEnv)
  env$loocv_elab <- F
}

#

return(argv)

}
