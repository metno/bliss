#+ Check input arguments and set some flags for elaborations
checkargs <- function( argv, env) {

#------------------------------------------------------------------------------

cat("Checking command line arguments: ")

#
#------------------------------------------------------------------------------
# check input arguments
#

if (!(file.exists(argv$iff_obs))) boom(paste0("file not found ",argv$iff_obs))

if ( argv$mode == "oi_multiscale_senorge_prec") {
  cat("chosen mode is OI multiscale developed for the spatial analysis of seNorge precipitation\n")
} else if ( argv$mode == "oi_twostep_senorge_temperature") {
  cat("chosen mode is OI two-step developed for the spatial analysis of seNorge temperature\n")
} else if ( argv$mode == "hyletkf") {
  cat("chosen mode is Hybrid Local Ensemble Transform Kalman Filter\n")
} else if ( argv$mode == "letkf") {
  cat("chosen mode is Local Ensemble Transform Kalman Filter\n")
} else if ( argv$mode == "successivecorrections") {
  cat("chosen mode is Successive Corrections\n")
} else if ( argv$mode == "oi") {
  cat("chosen mode is oi\n")
} else if ( argv$mode == "ensi") {
  cat("chosen mode is ensi\n")
} else if ( argv$mode == "ensigap") {
  cat("chosen mode is ensigap\n")
} else if ( argv$mode == "corensi") {
  cat("chosen mode is corensi\n")
} else if ( argv$mode == "msa") {
  cat("chosen mode is msa\n")
} else if ( argv$mode == "msaensi") {
  cat("chosen mode is msaensi\n")
} else if ( argv$mode == "msaensi_dev") {
  cat("chosen mode is msaensi_dev\n")
} else if ( argv$mode == "rasterize") {
  cat("chosen mode is rasterize\n")
} else {
  boom("error statistical interpolation scheme undefined\n")
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
#if (argv$mode %in% c("oi_multiscale_senorge_prec")) {
#  if ( !( file.exists(argv$path2src))) ext <- boom("path2src not found")
#  # remember: to get .so you have to "$>R CMD SHLIB oi_rr_fast.c"
#  dyn.load( file.path( argv$path2src, "oi_rr_first.so"))
#  dyn.load( file.path( argv$path2src, "oi_rr_fast.so"))
#}
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

# cross-validation
assign( "cv_mode_random", FALSE, envir = .GlobalEnv)
env$cv_mode_random <- F

assign( "cv_mode", FALSE, envir = .GlobalEnv)
env$cv_mode <- F

if ( argv$cv_mode | argv$cv_mode_random | 
     !is.na( argv$off_cv_table)   | !is.na( argv$off_cvt_table)   |
     !is.na( argv$off_cv_verif_a) | !is.na( argv$off_cvt_verif_a) |
     !is.na( argv$off_cv_verif_b) | !is.na( argv$off_cvt_verif_b)) {
  assign( "cv_elab", TRUE, envir = .GlobalEnv)
  env$cv_elab <- T
  if ( argv$cv_mode & any( !is.na( argv$prId.cv))) {
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

# leave-one-out cross-validation
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
