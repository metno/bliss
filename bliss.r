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
suppressPackageStartupMessages( library( "RANN"))
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

#--------debug or test
dirdeb <- "/home/cristianl/data/msaensi/debug"
ffdeb  <- file.path( dirdeb, paste0("debtest_beforemsaensi_", argv$date_out, ".rdata"))
if ( file.exists( ffdeb)) {
  load(ffdeb)
} else {

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
res <- read_fg( argv, fg_env, u_env, env)
if (!res) boom( code=1, str="ERROR problems while reading the background")
#
#------------------------------------------------------------------------------
# Read observations
res <- read_obs( argv, env, y_env)
if (!res) boom( code=1, str="ERROR problems while reading the observations")
#
#------------------------------------------------------------------------------
# Set the OI multi-scale parameters
if (argv$mode=="OI_multiscale") 
  source( file.path( bliss_mod_path, "main_oi_multicale_setpar.r"))
#
#------------------------------------------------------------------------------
# compute Disth (symmetric) matrix: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
if ( !(argv$mode %in% c( "hyletkf", "oi", "corens", "msaensi")) & y_env$yo$n < argv$maxobs_for_matrixInv ) {
  Disth <- matrix( ncol=y_env$yo$n, nrow=y_env$yo$n, data=0.)
  Disth <- ( outer(VecY,VecY,FUN="-")**2.+
             outer(VecX,VecX,FUN="-")**2. )**0.5/1000.
}
#------debug and test
save( file=ffdeb, argv, env, y_env, fg_env, u_env)
}
#
#------------------------------------------------------------------------------
# ANALYSIS
if (argv$verbose) 
  cat("+-------------------------------------------------------+\nAnalysis\n")
#..............................................................................
# ===> Rasterize  <===
if (argv$mode=="rasterize") { # still to test
  res <- rasterize_with_bliss( argv, y_env, env)
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
# ===>  Ensemble-based Statistical Interpolation  <===
} else if (argv$mode=="oi") {
  envtmp <- new.env( parent = emptyenv())
  res <- fg_u_align( argv, fg_env, u_env, env)
  res <- oi_driver( argv, y_env, fg_env, env)
#..............................................................................
# ===>  EnSI with Gaussian Anamorphosis for Precipitation  <===
} else if (argv$mode=="ensigap") {
  source( file.path( bliss_mod_path, "main_ensigap.r"))
#..............................................................................
# ===>  Change-of-Resolution Ensemble Rauch-Tung-Striebel smoother  <===
} else if (argv$mode=="corens") {
#next 4 lines are debug/test
ffdeb <- file.path( dirdeb, paste0("debtest_corens_", argv$date_out, ".rdata"))
if ( file.exists( ffdeb)) {
  load(ffdeb)
} else {
  envtmp <- new.env( parent = emptyenv())
  res <- corens_mergeobs( argv, y_env, u_env, env)
  rm(envtmp)
  res <- corens_selensemble( argv, fg_env, env)
  envtmp <- new.env( parent = emptyenv())
  res <- corens( argv, y_env, fg_env, env)
  rm( envtmp)
# next 2 lines are debug/test
save( file=ffdeb, argv, env, y_env, fg_env, u_env)
}
#..............................................................................
# ===>  Multiscale Alignment Ensemble Statistical Interpolation  <===
} else if (argv$mode=="msaensi") {
#next 4 lines are debug/test
ffdeb <- file.path( dirdeb, paste0("debtest_msaensi_", argv$date_out, ".rdata"))
if ( file.exists( ffdeb)) {
  load(ffdeb)
  envtmp <- new.env( parent = emptyenv())
  res <- msaensi( argv, y_env, fg_env, env)
} else {
  envtmp <- new.env( parent = emptyenv())
  res <- preproc_mergeobs( argv, y_env, u_env, env)
  res <- preproc_selensemble( argv, fg_env, env)
  save( file=ffdeb, argv, env, y_env, fg_env, u_env)
q()
#  res <- msaensi( argv, y_env, fg_env, env)
}
rm( envtmp)
q()
} # end if selection among spatial analysis methods
#
#------------------------------------------------------------------------------
if (argv$verbose) 
  cat("+-------------------------------------------------------+\nOutput\n")

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
