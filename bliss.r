#!/usr/bin/env Rscript
# --~- bliss.r -~--
# Bayesian statisticaL Interpolation for Spatial analySis
# See the software repository here: 
#..............................................................................
#Copyright and license
# Copyright (C) 2018 MET Norway. The software is licensed under GPL version 3 
# or (at your option) any later version.
# https://www.gnu.org/licenses/gpl-3.0.en.html
# 
# History:
# 26.10.2018 - Cristian Lussana. Original code.
# -----------------------------------------------------------------------------
#
rm(list=ls())
#
# -----------------------------------------------------------------------------
# Libraries
suppressPackageStartupMessages( library( "argparser"))
suppressPackageStartupMessages( library( "sp"))
suppressPackageStartupMessages( library( "raster"))
suppressPackageStartupMessages( library( "igraph"))
suppressPackageStartupMessages( library( "rgdal"))
suppressPackageStartupMessages( library( "ncdf4"))
suppressPackageStartupMessages( library( "dotnc"))
options( warn = 2, scipen = 999)
#options(scipen = 999)
# 
#..............................................................................
# Functions

#+ manage fatal error
boom <- function( str=NA, code=NA) {
  cat("Fatal Error ")
  if ( !is.na(code)) {
    if ( code == 1) cat("file not found ")
  }
  if ( !is.na(str)) cat( str)
  cat("\n")
  quit( status= 1 )
}

#+ the end 
rip <- function( str=NA, code=NA, t0=NA) {
  cat( "the End : ")
  if ( !is.na(code) ) {
    if ( code == 0 ) cat( "normal exit : ")
  }
  if ( !is.na(t0)) {
    t1 <- Sys.time()
    cat( paste( "total time=", round(t1-t0,1), attr(t1-t0,"unit")))
  }
  if ( !is.na(str)) cat( str)
  cat("\n")
  quit( status= 0 )
}

#
#==============================================================================
# MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN -
#==============================================================================
t0 <- Sys.time() # ladies and gentlemen, start your engines
#
#-----------------------------------------------------------------------------
# path to the titan functions is stored in the enviroment var BLISS_PATH 
bliss_fun_path <- file.path( Sys.getenv( "BLISS_PATH"), "src")
bliss_mod_path <- file.path( Sys.getenv( "BLISS_PATH"), "main_modules")
#
#-----------------------------------------------------------------------------
# load functions
fun_list <- c( "oi_var_gridpoint_by_gridpoint.r",
               "oi_bratseth_gridpoint_by_gridpoint.r",
               "oivar.r",
               "gamma_anamorphosis.r",
               "gamma_get_shape_rate_from_dataset.r",
               "oi_rr_fast.r",
               "henoi.r",
               "henoi_alt.r",
               "obsop_precip.r",
               "boxcox.r",
               "oi_temperature_util.r",
               "superobs.r",
               "misc_util.r",
               "read_and_regrid_nc.r",
               "debug_util.r")
for (fun in fun_list) {
  if ( !file.exists(file.path( bliss_fun_path, fun)))
    boom( file.path( bliss_fun_path, fun), code=1)
  source( file.path( bliss_fun_path, fun))
}
rm( fun_list, fun)               
#
#-----------------------------------------------------------------------------
# check all main modules exists
mod_list <- c( "main_constants.r", "main_argparser.r", "main_checkargs.r",
               "main_mastergrid.r", "main_emptygrid.r", "main_iff_rf.r",
               "main_laf.r", "main_dem.r", "main_iff_fg.r", "main_iff_obs.r",
               "main_oi_multiscale.r", "main_oi_multicale_setpar.r",
               "main_rasterize.r",
               "main_oi_fg.r",
               "main_oi_twosteps.r",
               "main_oi_bratseth.r",
               "main_sc_barnes.r",
#               "main_hyletkf.r",
               "main_ensi.r",
               "main_ensigap.r",
#               "main_wise.r",
               "main_off_y_table.r", "main_off_yt_table.r",
               "main_off_cv_table.r", "main_off_cvt_table.r",
               "main_off_lcv_table.r", "main_off_lcvt_table.r",
               "main_off_verif.r",
               "main_off_rdata.r",
               "main_off_x.r")
for (mod in mod_list) {
  if ( !file.exists(file.path( bliss_mod_path, mod)))
    boom( file.path( bliss_mod_path, mod), code=1)
}
rm( mod_list, mod)               
source( file.path( bliss_mod_path, "main_iff_fg_wise.r"))
source( file.path( bliss_mod_path, "main_wise_align.r"))
source( file.path( bliss_mod_path, "main_wise_analysis.r"))
source( file.path( bliss_mod_path, "main_wise_analysis_loop.r"))
source( file.path( bliss_mod_path, "main_wise_analysis_loop_alignment.r"))
source( file.path( bliss_mod_path, "main_wise_sampling_postpdf.r"))
source( file.path( bliss_mod_path, "main_wise_plot.r"))
source( file.path( bliss_mod_path, "main_argparser.r"))

#
#-----------------------------------------------------------------------------
# define constants
source( file.path( bliss_mod_path, "main_constants.r"))

#
#-----------------------------------------------------------------------------
# define environments

# main environment
env <- new.env( parent = emptyenv())

# first guess (background) environment
fg_env    <- new.env( parent=emptyenv())
fg_env$fg <- list()

# observation enviroment
y_env    <- new.env( parent=emptyenv())
y_env$yo <- list()

# observation environment (alignment)
u_env    <- new.env( parent=emptyenv())
u_env$uo <- list()

#
#-----------------------------------------------------------------------------
# read command line arguments and/or configuration file
argv <- argparser( env, fg_env, y_env, u_env)

#
#-----------------------------------------------------------------------------
# Multi-cores run
if ( !is.na( argv$cores)) {
  suppressPackageStartupMessages( library( "parallel"))
  if ( argv$cores==0) argv$cores <- detectCores()
  cat( paste( "--> multi-core run, cores=", argv$cores, "\n"))
}

#-----------------------------------------------------------------------------
# checks on input arguments
source( file.path( bliss_mod_path, "main_checkargs.r"))

#
#------------------------------------------------------------------------------
# Create master grid
source( file.path( bliss_mod_path, "main_mastergrid.r"))
env$rmaster <- rmaster

#
#------------------------------------------------------------------------------
# Empty grid on a gridded output
if ( argv$empty_grid & !is.na( argv$off_x) ) {
  source( file.path( bliss_mod_path, "main_emptygrid.r"))
  rip( code=0, t0=t0)
}
#
#------------------------------------------------------------------------------
# read rescaling factor
if (file.exists(argv$iff_rf)) source( file.path( bliss_mod_path, "main_iff_rf.r"))
#
#------------------------------------------------------------------------------
# read land area fraction
if (file.exists(argv$iff_laf)) source( file.path( bliss_mod_path, "main_laf.r"))
#
#------------------------------------------------------------------------------
# read digital elevation model 
if (file.exists(argv$iff_dem)) source( file.path( bliss_mod_path, "main_dem.r"))
#
#------------------------------------------------------------------------------
# -~- read first guess -~-
# output data structures:
# - xb, deterministic (vector) or ensemble (array = dim( ngrid,nens) field(s)
#    xb values are stored only for not NAs points, masked using rmaster
#    xb values less than min_plaus_val are set to min_plaus_val
# - aix, vector. indices of background-valid gridpoints wrt the original grid 
#    note: ngrid = length(aix)
# - rfg, raster file with the masked background field 
#        (if ensemble, then is the last one)
# - valens, index to the valid ensemble members 
# - iff_is_ens. is iff_fg an ensemble? used when writing output
#
if (file.exists(argv$iff_fg)) source( file.path( bliss_mod_path, "main_iff_fg.r")) 
if (argv$mode=="wise") {
  dir_plot<-"/home/cristianl/data/wise"
  dir_plot <- file.path(dir_plot,paste0("case_",argv$date_out))
  load_if_present <- T
  ffff<- file.path(dir_plot,paste0("tmp_wise_input_",argv$date_out,".rdata"))
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    res <- main_iff_fg_wise( argv, fg_env, u_env, env)
    save(file=ffff,argv, fg_env, u_env, env)
  }
}
#
#------------------------------------------------------------------------------
# Read observations
source( file.path( bliss_mod_path, "main_iff_obs.r"))
if (argv$mode=="wise") {
  y_env$yo$x <- VecX
  y_env$yo$y <- VecY
  y_env$yo$value <- yo
  ffff<- file.path(dir_plot,paste0("tmp_wise_input_",argv$date_out,".rdata"))
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }
}
#
#------------------------------------------------------------------------------
# Set the OI multi-scale parameters
if (argv$mode=="OI_multiscale") 
  source( file.path( bliss_mod_path, "main_oi_multicale_setpar.r"))
#
#------------------------------------------------------------------------------
# compute Disth (symmetric) matrix: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
if ( argv$mode != "hyletkf" & n0 < argv$maxobs_for_matrixInv ) {
  Disth <- matrix( ncol=n0, nrow=n0, data=0.)
  Disth <- ( outer(VecY,VecY,FUN="-")**2.+
             outer(VecX,VecX,FUN="-")**2. )**0.5/1000.
}
#
#------------------------------------------------------------------------------
# ANALYSIS
if (argv$verbose) 
  cat("+-------------------------------------------------------+\nAnalysis\n")
#..............................................................................
# ===> Rasterize  <===
if (argv$mode=="rasterize") {
  source( file.path( bliss_mod_path, "main_rasterize.r"))
#..............................................................................
# ===> OI with deterministic background  <===
} else if (argv$mode=="OI_firstguess") {
  source( file.path( bliss_mod_path, "main_oi_fg.r"))
#..............................................................................
# ===>  OI multiscale (without background)   <===
} else if (argv$mode=="OI_multiscale") {
  source( file.path( bliss_mod_path, "main_oi_multiscale.r"))
#..............................................................................
# ===>  OI two-step spatial interpolation (without background)   <===
} else if (argv$mode=="OI_twosteptemperature") {
  source( file.path( bliss_mod_path, "main_oi_twosteps.r"))
#..............................................................................
# ===>  OI Bratseth, Brathset's iterative method for OI (with background)  <===
} else if (argv$mode=="OI_Bratseth") {
  source( file.path( bliss_mod_path, "main_oi_bratseth.r"))
#..............................................................................
# ===>  SCB, Successive Correction Barnes' scheme (with background)  <===
} else if (argv$mode=="SC_Barnes") {
  source( file.path( bliss_mod_path, "main_sc_barnes.r"))
#..............................................................................
# ===>  Hybrid Local Ensemble Transform Kalman Filter  <===
#} else if (argv$mode=="hyletkf") {
#  source( file.path( bliss_mod_path, "main_hyletkf.r"))
#..............................................................................
# ===>  Ensemble-based Statistical Interpolation  <===
} else if (argv$mode=="ensi") {
  source( file.path( bliss_mod_path, "main_ensi.r"))
#..............................................................................
# ===>  EnSI with Gaussian Anamorphosis for Precipitation  <===
} else if (argv$mode=="ensigap") {
  source( file.path( bliss_mod_path, "main_ensigap.r"))
#..............................................................................
# ===>  Wavelet statistical interpolation  <===
} else if (argv$mode=="wise") {
  ffff<- file.path(dir_plot,paste0("tmp_wise_align_",argv$date_out,".rdata"))
  load_if_present<-T
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    res <- main_wise_align( argv, fg_env, u_env, env, plot=F, dir_plot=dir_plot)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }
  ffff<- file.path(dir_plot,paste0("tmp_wise_analysis_",argv$date_out,".rdata"))
  load_if_present<-F
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    res <- main_wise_analysis_loop( argv, y_env, fg_env, env, seed=1, obs_k_dim=30, plot=T, dir_plot=dir_plot)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }
  ffff<- file.path(dir_plot,paste0("tmp_wise_postpdf_",argv$date_out,".rdata"))
q()
  load_if_present<-T
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    res <- main_wise_sampling_postpdf( argv, y_env, fg_env, u_env, env, seed=100, plot=T, dir_plot=dir_plot)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }
  res <- main_wise_plot( argv, y_env, fg_env, env, dir_plot=dir_plot)
q()
} # end if selection among spatial analysis methods
#
#------------------------------------------------------------------------------
if (argv$verbose) 
  cat("+-------------------------------------------------------+\nOutput\n")
#
# -- text files --
# table
if (!is.na(argv$off_y_table)) 
  source( file.path( bliss_mod_path, "main_off_y_table.r"))
# table (transformed data)
if (!is.na(argv$off_yt_table)) 
  source( file.path( bliss_mod_path, "main_off_yt_table.r"))
# table cross-validation
if (!is.na(argv$off_cv_table)) 
  source( file.path( bliss_mod_path, "main_off_cv_table.r"))
# table cross-validation (transformed data)
if (!is.na(argv$off_cvt_table)) 
  source( file.path( bliss_mod_path, "main_off_cvt_table.r"))
# table leave-one-out cross-validation
if (!is.na(argv$off_lcv_table)) 
  source( file.path( bliss_mod_path, "main_off_lcv_table.r"))
# table leave-one-out cross-validation (transformed data)
if (!is.na(argv$off_lcvt_table)) 
  source( file.path( bliss_mod_path, "main_off_lcvt_table.r"))
#
# -- netcdf - verif --
source( file.path( bliss_mod_path, "main_off_verif.r"))
#
# -- rdata --
if ( !is.na( argv$off_rdata)) 
  source( file.path( bliss_mod_path, "main_off_rdata.r"))
#
# -- netcdf - gridded output --
if ( !is.na( argv$off_x)) 
  source( file.path( bliss_mod_path, "main_off_x.r"))
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# normal exit
rip( code=0, t0=t0)
