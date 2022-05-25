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
# get path from system environment
bliss_path <- Sys.getenv( "BLISS_PATH")

# load the function we need to read the command line arguments
source( file.path( bliss_path, "src", "argparser.r"))

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
y_env$yo  <- list()
y_env$yov <- list() # list for cv

# observation environment (alignment)
u_env    <- new.env( parent=emptyenv())
u_env$uo <- list()

#
#-----------------------------------------------------------------------------
# read command line arguments and/or configuration file
argv <- argparser( env, fg_env, y_env, u_env)

#
#------------------------------------------------------------------------------
# Load functions
for (file in list.files(path = file.path( bliss_path, "src"), pattern = ".r$", full.names=T)) 
  source(file)
for (file in list.files(path = file.path( bliss_path, "methods"), pattern = ".r$", full.names=T)) 
  source(file)
rm(file)

#
#-----------------------------------------------------------------------------
# define constants
define_constants( env)

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
argv <- checkargs( argv, env)

#
#------------------------------------------------------------------------------
# Create master grid (output env$rmaster)
define_mastergrid( argv, env)

#
#------------------------------------------------------------------------------
# Just write an empty grid and exit
if ( argv$empty_grid & !is.na( argv$off_x) ) {
  main_emptygrid( argv, env)
  rip( code=0, t0=t0)
}

#
#------------------------------------------------------------------------------
# read rescaling factor
read_rescaling_factor( argv, env)

#
#------------------------------------------------------------------------------
# read land area fraction
# output: env$rlaf_exists, env$rlaf
read_laf( argv, env)

#
#------------------------------------------------------------------------------
# read digital elevation model 
# output: env$rlaf_exists, env$rlaf
read_dem( argv, env)

#
#------------------------------------------------------------------------------
# Read background fields
# fg_env$nfg = number of files where the background fields are stored
# fg_env$fg[[f]] = f-th element of the list, corresponding to the f-th file,
#                  where the background data and metadata are stored

dir_plot <- "/home/cristianl/data/wise/pngs"
ffff<- file.path(dir_plot,paste0("tmp_fg_",argv$date_out,".rdata"))
load_if_present<-F
if (file.exists(ffff) & load_if_present) {
  load(ffff)
} else {
  res <- read_fg( argv, fg_env, u_env, env)
  if (!res) boom( code=1, str="ERROR problems while reading the background")
  save(file=ffff, argv, fg_env, u_env, env, y_env)
}

#
#------------------------------------------------------------------------------
# Read observations

ffff<- file.path(dir_plot,paste0("tmp_obs_",argv$date_out,".rdata"))
load_if_present<-F
if (file.exists(ffff) & load_if_present) {
  load(ffff)
} else {
  res <- read_obs( argv, env, y_env)
  if (!res) boom( code=1, str="ERROR problems while reading the observations")
#  png(file="tmp_obs.png",width=1200,height=1200)
#  plot(y_env$yo$x,y_env$yo$y,pch=21,bg="gray",col="darkgray")
#  points(y_env$yov$x,y_env$yov$y,pch=21,bg="red",col="darkred")
#  dev.off()
  save(file=ffff, argv, fg_env, u_env, env, y_env)
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
if ( !(argv$mode %in% c( "hyletkf", "wise", "oi")) & y_env$yo$n < argv$maxobs_for_matrixInv ) {
  Disth <- matrix( ncol=y_env$yo$n, nrow=y_env$yo$n, data=0.)
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
} else if (argv$mode=="oi") {
  suppressPackageStartupMessages( library( "RANN"))
  envtmp <- new.env( parent = emptyenv())
  res <- fg_u_align( argv, fg_env, u_env, env)
  res <- oi_driver( argv, y_env, fg_env, env)
#..............................................................................
# ===>  EnSI with Gaussian Anamorphosis for Precipitation  <===
} else if (argv$mode=="ensigap") {
  source( file.path( bliss_mod_path, "main_ensigap.r"))
#..............................................................................
# ===>  Wavelet statistical interpolation  <===
} else if (argv$mode=="wise") {

  suppressPackageStartupMessages( library( "waveslim"))
  suppressPackageStartupMessages( library( "RANN"))

  ffff<- file.path(dir_plot,paste0("tmp_wise_align_",argv$date_out,".rdata"))
  load_if_present<-F
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    t00<-Sys.time()
    res <- wise_align( argv, fg_env, u_env, env, plot=F, dir_plot=dir_plot)
    print(Sys.time()-t00)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }

  ffff<- file.path(dir_plot,paste0("tmp_wise_analysis_",argv$date_out,".rdata"))
  load_if_present<-F
  plot<-F
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    t00<-Sys.time()
    res <- wise_analysis_loop( argv, y_env, fg_env, env, 
                               supob_nobs=argv$wise_supob_nobs,
                               supob_radius=argv$wise_supob_radius,
                               supob_q=argv$wise_supob_q,
                               max_it=argv$wise_opt_maxit,
                               opttol=argv$wise_opt_opttol,
                               En2_adj_fun=argv$wise_En2_adj_fun,
                               En2_adj_min=argv$wise_En2_adj_min,
                               rescale_min_obs=argv$wise_rescale_min_obs,
                               rescale_min_cells=rgv$wise_rescale_min_cells,
                               plot=plot, dir_plot=dir_plot)
    print(Sys.time()-t00)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }

  ffff<- file.path(dir_plot,paste0("tmp_wise_postpdf_",argv$date_out,".rdata"))
  load_if_present<-F
  plot<-F
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    t00<-Sys.time()
    res <- wise_sampling_postpdf( argv, y_env, fg_env, u_env, env,
                                  resample=argv$wise_a_resample, 
                                  seed=argv$wise_a_resample_setseed,
                                  plot=plot, dir_plot=dir_plot)
    print(Sys.time()-t00)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }

  ffff<- file.path(dir_plot,paste0("tmp_wise_agg_",argv$date_out,".rdata"))
  load_if_present<-F
  plot<-F
  if (file.exists(ffff) & load_if_present) {
    load(ffff)
  } else {
    t00<-Sys.time()
    res <- wise_aggregation( argv, y_env, fg_env, u_env, env,
                             plot=plot, dir_plot=dir_plot)
    print(Sys.time()-t00)
    save(file=ffff, argv, fg_env, u_env, env, y_env)
  }

#  res <- main_wise_plot( argv, y_env, fg_env, env, dir_plot=dir_plot)

} # end if selection among spatial analysis methods
#
#------------------------------------------------------------------------------
if (argv$verbose) 
  cat("+-------------------------------------------------------+\nOutput\n")
#
# -- text files --
# table
#if (!is.na(argv$off_y_table)) 
#  source( file.path( bliss_mod_path, "main_off_y_table.r"))
## table (transformed data)
#if (!is.na(argv$off_yt_table)) 
#  source( file.path( bliss_mod_path, "main_off_yt_table.r"))
## table cross-validation
#if (!is.na(argv$off_cv_table)) 
#  source( file.path( bliss_mod_path, "main_off_cv_table.r"))
## table cross-validation (transformed data)
#if (!is.na(argv$off_cvt_table)) 
#  source( file.path( bliss_mod_path, "main_off_cvt_table.r"))
## table leave-one-out cross-validation
#if (!is.na(argv$off_lcv_table)) 
#  source( file.path( bliss_mod_path, "main_off_lcv_table.r"))
## table leave-one-out cross-validation (transformed data)
#if (!is.na(argv$off_lcvt_table)) 
#  source( file.path( bliss_mod_path, "main_off_lcvt_table.r"))
#
# -- netcdf - verif --
#source( file.path( bliss_mod_path, "main_off_verif.r"))

#
# -- rdata --
if ( !is.na( argv$off_yenv_rdata))  write_rdata( file=argv$off_yenv_rdata,  arguments=argv, environment=y_env) 
if ( !is.na( argv$off_env_rdata))   write_rdata( file=argv$off_env_rdata,   arguments=argv, environment=env) 
if ( !is.na( argv$off_fgenv_rdata)) write_rdata( file=argv$off_fgenv_rdata, arguments=argv, environment=fg_env) 
if ( !is.na( argv$off_uenv_rdata))  write_rdata( file=argv$off_uenv_rdata,  arguments=argv, environment=u_env) 

#
# -- netcdf - gridded output --
if ( !is.na( argv$off_x))  write_off_x_nc(  argv,  y_env, fg_env, u_env, env) 
if ( !is.na( argv$off_xb)) write_off_xb_nc( argv,  y_env, fg_env, u_env, env) 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# normal exit
rip( code=0, t0=t0)
