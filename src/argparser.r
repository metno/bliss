#+
argparser<-function( env, fg_env, y_env, u_env) {
# create parser object
p <- arg_parser("bliss")
#------------------------------------------------------------------------------
# miscellaneous 
p <- add_argument(p, "--verbose",
                  help="verbose mode",
                  flag=T,
                  short="-v")
p <- add_argument(p, "--debug",
                  help="debug mode",
                  flag=T)
p <- add_argument(p, "--rrinf",
                  help="precipitation yes/no threshold",
                  type="numeric",
                  default=0.1)
#.............................................................................. 
#
# first-guess configuration files
p <- add_argument(p, "--fg.files",
                  help="information used to read first-guess nc-file(s) (use conf <- list( var1=..., var2=..., ... )",
                  type="character",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--fg.filenames",
                  help="file names of the first-guess nc-files. It is an optional argument, if specified: i) must have the same length of --fg.files ii) override the main.file arguments in the fg.files",
                  type="character",
                  default=NULL,
                  nargs=Inf)
#------------------------------------------------------------------------------
p <- add_argument(p, "--uo.file",
                  help="information used to read observation file used for alignment (use conf <- list( var1=..., var2=..., ... )",
                  type="character",
                  default=NA)
p <- add_argument(p, "--uo.filename",
                  help="file names of the observation nc-file used for alignment. It is an optional argument, if specified it overrides the main.file argument in the uo.file",
                  type="character",
                  default=NA)
#------------------------------------------------------------------------------
# cross-validation mode
p <- add_argument(p, "--cv_mode",
                  help="standard cross-validation mode (exclude one observation provider from the analysis and use it for cv, see also ''prId.cv'')",
                  flag=T)
p <- add_argument(p, "--prId.cv",
                  help="observation provider identifiers to reserve for cross-validation",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--cv_mode_random",
                  help="standard cross-validation mode (random selection of stations, pre-set minimum distance between them is ''d_cv'')",
                  flag=T)
p <- add_argument(p, "--d.cv",
                  help="distance used for ''cv_mode_random'' (km)",
                  type="numeric",
                  default=50)
p <- add_argument(p, "--cv_mode_random_setseed",
                  help="set the seed used by the random number generator, this way one always get the same set of stations. Useful for testing/debugging.",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--loocv_mode",
                  help="leave-one-out cross-validation mode",
                  type="logical",
                  default=F)
p <- add_argument(p, "--cv_mode_calcidiv",
                  help="compute IDI-cv",
                  type="logical",
                  default=F)
p <- add_argument(p, "--calcidiv_nobs",
                  help="IDI-cv maximum number of observations to use locally",
                  type="integer",
                  default=50)
p <- add_argument(p, "--calcidiv_radius",
                  help="IDI-cv radius of the neighbourhood",
                  type="numeric",
                  default=30000)
p <- add_argument(p, "--calcidiv_dh",
                  help="IDI-cv decorrelation lenght",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--idiv_instead_of_elev",
                  help="compute IDI-cv",
                  type="logical",
                  default=F)
p <- add_argument(p, "--twostep_superobbing",
                  help="superobbing (used only if \"OI_twosteptemperature\")",
                  flag=T)
p <- add_argument(p, "--twostep_nogrid",
                  help="calculation only for station points",
                  flag=T)
#------------------------------------------------------------------------------
# statistical interpolation mode
p <- add_argument(p, "--mode",
#                  help="statistical interpolation scheme (\"rasterize\",\"OI_multiscale\",\"OI_firstguess\",\"OI_twosteptemperature\",\"SC_Barnes\",\"OI_Bratseth\",\"hyletkf\",\"letkf\",\"ensip\")",
                  help="statistical interpolation scheme (\"rasterize\",\"OI_multiscale\",\"OI_firstguess\",\"OI_twosteptemperature\",\"SC_Barnes\",\"OI_Bratseth\",\"ensi\",\"ensigap\", \"wise\", \"oi\")",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# time-related variables
p <- add_argument(p, "--date_out",
                  help="date to write in output files (nc, %Y%m%d%H%M)",
                  type="character",
                  default=NULL)
p <- add_argument(p, "--date_out_fmt",
                  help="date out format",
                  type="character",
                  default="%Y%m%d%H%M")
p <- add_argument(p, "--time_bnds_string",
                  help="time bounds with respect to date_out (e.g., \"-1 day\" \"-1 min\")",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# rasterize
# output variables in the netcdf are "mean_raster", "sd_raster",
#  "n_raster", "q_raster_XX" (e.g. XX=01 stands for 1st percentile).
# q_raster_XX must have the same XX value in rasterize_q to actually generate
# the output field (e.g. XX=01 corresponds to rasterize_q=0.01)
p <- add_argument(p, "--rasterize_nmin",
                  help="set to NA if less than this number of observations are within a box",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--rasterize_q",
                  help="quantiles to be computed",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
#------------------------------------------------------------------------------
# OI shared
p <- add_argument(p, "--corrfun",
                  help="correlation functions (\"gaussian\",\"soar\")",
                  type="character",
                  default="gaussian")
p <- add_argument(p, "--pmax",
                  help="maximum number of observations in the neighbourhood of a gridpoint for OI",
                  type="numeric",
                  default=200)
# OI_multiscale / OI_firstguess parameters
p <- add_argument(p, "--eps2",
                  help="ratio of observation to background error covariance",
                  type="numeric",
                  default=1)
p <- add_argument(p, "--Dh",
                  help="horizontal de-corellation length scale (km)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--oimult.eps2_idi",
                  help="OI multiscale, optional, eps2 for the IDI",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--oimult.Dh_idi",
                  help="OI multiscale, optional, horizontal de-corellation length scale for IDI (km)",
                  type="numeric",
                  default=NULL)
# NOTE: ovarc is somewhat equivalent to oifg.eps2_r...
p <- add_argument(p, "--ovarc",
                  help="observation error variance correction factor",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--ovarc.prId",
                  help="observation error variance correction factor (prId)",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--prId.exclude",
                  help="observation provider identifiers to exclude",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
# OI_firstguess
p <- add_argument(p, "--oifg.eps2",
                  help="eps2 values for \"OI_firstguess\" mode, prId- and observation- dependent",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--oifg.eps2_prId",
                  help="provided Id associated with eps2, NA stands for default",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--oifg.eps2_r",
                  help="threshold for observed value (e.g., two different eps2 for observations [0,1) and [1,1000) is oifg.eps2_r=c(0,1,1000))",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--oifg.Dh",
                  help="horizontal decorrelation lenght for \"OI_firstguess\" mode (km)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--oifg.xta_errvar_smooth",
                  help="smoothing length for the analysis error variance (m), used only in case of data transformation",
                  type="numeric",
                  default=50000)

# additional OI parameters two-step temperature
p <- add_argument(p, "--nmaxo",
                  help="number of closest observations to be used in the OI to adjust background values at a grid point",
                  type="numeric",
                  default=20)
p <- add_argument(p, "--dz",
                  help="background error covariance matrix, vertical decorrelation length scale (m)",
                  type="numeric",
                  default=600)
p <- add_argument(p, "--lafmin",
                  help="background error covariance matrix, land area fraction parameter",
                  type="numeric",
                  default=0.5)
# Background parameters
p <- add_argument(p, "--grid.bg",
                  help="nrow ncol (i.e. number_of_rows number_of_columns) used to define grid of sub-regions.",
                  type="integer",
                  nargs=2,
                  default=c(5,5))
p <- add_argument(p, "--vmin",
                  help="minimum allowed value [units of the variable specified]",
                  type="numeric",
                  default=-50)
p <- add_argument(p, "--vmax",
                  help="maximum allowed value [units of the variable specified]",
                  type="numeric",
                  default=40)
p <- add_argument(p, "--gamma.standard",
                  help="standard value for the moist adiabatic temperature lapse rate dT/dz (degC/m)",
                  type="numeric",
                  default=-0.0065)
p <- add_argument(p, "--n.bg",
                  help="minimum number of observations allowed in a sub-region",
                  type="numeric",
                  default=50)
p <- add_argument(p, "--dz.bg",
                  help="elevation difference. If the elevation range (=elev 95th-perc minus elev 5th-perc) is less than this value, then fit a linear profile of temperature. Otherwise, fit a non-linear profile.",
                  type="numeric",
                  default=30)
p <- add_argument(p, "--nclose.bg",
                  help="n-th closest observation to consider in the calculation of the OI de-correlation length",
                  type="numeric",
                  default=4)
p <- add_argument(p, "--maxboxl",
                  help="maximum length (m) of the box used to define a sub-region",
                  type="numeric",
                  default=250000)
p <- add_argument(p, "--obs.outbuffer",
                  help="distance (m) defining the \"buffer\" region outside the masked region where to consider the observation, so to reduce border effects",
                  type="numeric",
                  default=50000)
#------------------------------------------------------------------------------
# OI_Bratseth
p <- add_argument(p, "--oibr.nSCloops",
                  help="number of Successive Corrections loops (\"OI_Bratset\")",
                  type="numeric",
                  default=10)
p <- add_argument(p, "--oibr.dh",
                  help="horizontal decorellation length scale (m) (\"OI_Bratset\")",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--oibr.delta_hor",
                  help="horizontal displacement (m) used for adaptive estimation of dh (\"OI_Bratset\")",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--oibr.box_o_nearest_halfwidth",
                  help="search neighbours within a box (m) (\"OI_Bratset\")",
                  type="numeric",
                  default=200000)
p <- add_argument(p, "--oibr.dz",
                  help="vertical decorellation length scale (m) (\"OI_Bratset\")",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--oibr.lafmin",
                  help="land-area fraction minimum (\"OI_Bratset\")",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--oibr.dh_adaptive",
                  help="adaptive estimation of dh (\"OI_Bratset\")",
                  flag=T)
p <- add_argument(p, "--oibr.dh_adaptive_min",
                  help="minimum allowed distance (m) (\"OI_Bratset\")",
                  type="numeric",
                  default=0)
p <- add_argument(p, "--oibr.dh_adaptive_max",
                  help="maximum allowed distance (m) (\"OI_Bratset\")",
                  type="numeric",
                  default=0)
#------------------------------------------------------------------------------
# hyletkf
p <- add_argument(p, "--hyletkf.eps2_prec_default",
                  help="LETKF default value for the ratio between obs error variance and backg error variance in case of precipitation",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--hyletkf.eps2_noprec_default",
                  help="LETKF default value for the ratio between obs error variance and backg error variance in case of no-precipitation",
                  type="numeric",
                  default=0.05)
p <- add_argument(p, "--hyletkf.Dh",
                  help="horizontal de-correlation length for the LETKF localization",
                  type="numeric",
                  default=10)
p <- add_argument(p, "--hyletkf.pmax",
                  help="maximum number of observations in the neighbourhood of a gridpoint for LETKF",
                  type="numeric",
                  default=200)
p <- add_argument(p, "--hyletkf.sigma2_min",
                  help="minimum allowed background variance (in the transformed space) for LETKF",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--hyletkf.eps2_prId",
                  help="provider identifier corresponding to the specified ratio between obs error variance and backg error variance (LETKF)",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--hyletkf.eps2",
                  help="observation-provider dependent ratio between obs error variance and backg error variance (LETKF)",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--hyletkf.Dh_oi",
                  help="horizontal de-correlation length for the OI step of the HyLETKF (km)",
                  type="numeric",
                  default=10)
p <- add_argument(p, "--hyletkf.eps2_oi",
                  help="ratio between obs error variance and backg error variance for the OI step of the HyLETKF",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--hyletkf.rloc_min",
                  help="use an observation in LETKF only if the localization weight is greater than the specified minimum value",
                  type="numeric",
                  default=0.0013)
#------------------------------------------------------------------------------
# ensip - ensemble-based statistical interpolation of precipitation
p <- add_argument(p, "--ensip.pmax",
                  help="max number of neighbouring observations to consider",
                  type="numeric",
                  default=200)
p <- add_argument(p, "--ensip.rloc_min",
                  help="do not use observations if they are too far from a location. pmax is still they strongest constrain. Set this parameter to 0 to use always pmax observations.",
                  type="numeric",
                  default=0)
p <- add_argument(p, "--ensip.henoi_par_notadaptive",
                  help="if true use the default values for henoi parameters",
                  flag=T)
p <- add_argument(p, "--ensip.henoi_Dh_loc",
                  help="Reference lenght scale for localization of ensemble-based background error covariance matrix (m, or same unit as grid/observation coords).",
                  type="numeric",
                  default=20000)
p <- add_argument(p, "--ensip.henoi_Dh",
                  help="Reference lenght scale for static background error covariance matrix (m, or same unit as grid/observation coords).",
                  type="numeric",
                  default=20000)
p <- add_argument(p, "--ensip.henoi_eps2",
                  help="eps2",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--ensip.henoi_alpha",
                  help="alpha",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--ensip.var_o_coeff",
                  help="coefficient(s) for the observation error variance (default)",
                  type="numeric",
                  default=1)
p <- add_argument(p, "--ensip.var_o_coeff_special",
                  help="coefficient(s) for the observation error variance (special)",
                  type="numeric",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--ensip.var_o_coeff_special_prId",
                  help="coefficient(s) for the observation error variance (special, Id)",
                  type="numeric",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--ensip.henoi_reflen_min",
                  help="HEnOI adaptive parameter estimation (Dh and Dh_loc are set to the same value). Minimum value for the reference lenght scale for localization of ensemble-based background error covariance matrix (m, or same unit as grid/observation coords).",
                  type="numeric",
                  default=3000)
p <- add_argument(p, "--ensip.henoi_reflen_max",
                  help="HEnOI adaptive parameter estimation (Dh and Dh_loc are set to the same value). Maximum value for the reference lenght scale for localization of ensemble-based background error covariance matrix (m, or same unit as grid/observation coords).",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--ensip.henoi_reflen_aggfact",
                  help="HEnOI adaptive parameter estimation. The estimation needs a coarser grid and this is the aggregation factor applied to the original grid (number of grid cells).",
                  type="numeric",
                  default=5)
p <- add_argument(p, "--ensip.henoi_reflen_k",
                  help="HEnOI adaptive parameter estimation. The reflen is estimated for each grid point on a coarser grid (\"henoi_reflen_aggfact\") such that it coincides with the distance to the closest \"k\" observations, provided that is within a predetermined range (\"ensip.henoi_reflen_min\",\"ensip.henoi_reflen_max\") (k units = number of observations).",
                  type="numeric",
                  default=10)
p <- add_argument(p, "--ensip.no_data_transf",
                  help="do not transform data",
                  flag=T)
p <- add_argument(p, "--ensip.shape",
                  help="Gamma distribution shape parameter",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--ensip.rate",
                  help="Gamma distribution rate parameter",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--ensip.henoi_statcov_backg",
                  help="HEnOI static error covariance matrix function (exp or gauss)",
                  type="character",
                  default="gauss")
p <- add_argument(p, "--ensip.backg_best",
                  help="strategy to get the best background",
                  type="character",
                  default="mean")
#------------------------------------------------------------------------------
# paths
p <- add_argument(p, "--path2src",
                  help="path to the shared objects (.so files)",
                  type="character",
                  default="none")
p <- add_argument(p, "--path2bliss",
                  help="path to the bliss main directory",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# input files
p <- add_argument(p, "--config.file",
                  help="configuration file",
                  type="character",
                  default=NULL,
                  short="cf")
p <- add_argument(p, "--iff_obs",
                  help="full file name for the observations (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg",
                  help="full file name for the first-guess field (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf",
                  help="full file name for the rescaling factor (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem",
                  help="full file name for the digital elevation model (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf",
                  help="full file name for the land area fraction (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask",
                  help="full file name for the mask (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_black",
                  help="blacklist file",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# output files
p <- add_argument(p, "--empty_grid",
                  help="create an output file on the grid with all NAs",
                  flag=T)
p <- add_argument(p, "--off_x",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_xb",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_rdata",
                  help="full file name for output (rdata)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_yenv_rdata",
                  help="full file name for output at station locations (rdata, y_env environment)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_env_rdata",
                  help="full file name for output (rdata, env environment)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_fgenv_rdata",
                  help="full file name for output (rdata, fg_env environment)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_uenv_rdata",
                  help="full file name for output (rdata, u_env environment)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_y_table",
                  help="full file name for output at station locations (txt)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_yt_table",
                  help="full file name for output at station locations, transformed values (txt)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_cv_table",
                  help="full file name for output at station locations, cross-validation (txt)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_cvt_table",
                  help="full file name for output at station locations, cross-validation, transformed values (txt)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_lcv_table",
                  help="full file name for output at station locations, leave-one-out cross-validation (txt)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_lcvt_table",
                  help="full file name for output at station locations, leave-one-out cross-validation, transformed values (txt)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_y_verif_a",
                  help="full file name for output at station locations analysis vs obs, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_y_verif_av",
                  help="full file name for output at station locations cv-analysis vs obs, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_y_verif_b",
                  help="full file name for output at station locations background vs obs, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_yt_verif_a",
                  help="full file name for output at station locations analysis vs obs, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_yt_verif_av",
                  help="full file name for output at station locations cv-analysis vs obs, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_yt_verif_b",
                  help="full file name for output at station locations background vs obs, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_cv_verif_a",
                  help="full file name for output at station locations analysis vs obs, cross-validation, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_cv_verif_b",
                  help="full file name for output at station locations background vs obs, cross-validation, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_cvt_verif_a",
                  help="full file name for output at station locations analysis vs obs, cross-validation, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_cvt_verif_b",
                  help="full file name for output at station locations background vs obs, cross-validation, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_lcv_verif_a",
                  help="full file name for output at station locations analysis vs obs, leave-one-out cross-validation, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_lcv_verif_b",
                  help="full file name for output at station locations background vs obs, leave-one-out cross-validation, verif format (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_lcvt_verif_a",
                  help="full file name for output at station locations analysis vs obs, leave-one-out cross-validation, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--off_lcvt_verif_b",
                  help="full file name for output at station locations background vs obs, leave-one-out cross-validation, verif format, transformed values (nc)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--verif_var",
                  help="name of the verif variable",
                  type="character",
                  default="Precipitation")
p <- add_argument(p, "--verif_units",
                  help="units of the verif variable",
                  type="character",
                  default="mm")
p <- add_argument(p, "--off_vernc.varname",
                  help="name of the verif variable for off_vernc files",
                  type="character",
                  default="Precipitation")
p <- add_argument(p, "--off_vernc.stdvarname",
                  help="standard name of the verif variable for off_vernc files",
                  type="character",
                  default="precipitation_amount")
p <- add_argument(p, "--off_vernc.varunits",
                  help="units of the verif variable for off_vernc files",
                  type="character",
                  default="mm")
p <- add_argument(p, "--off_vernc.thresholds",
                  help="thresholds for the off_vernc files",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--off_vernc.quantile",
                  help="quantiles for the off_vernc files",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
#------------------------------------------------------------------------------
# customization of the observation file
p <- add_argument(p, "--iff_obs.sep",
                  help="separator character",
                  type="character",
                  default=";")
p <- add_argument(p, "--iff_obs.x",
                  help="easting coordinate name",
                  type="character",
                  default="x")
p <- add_argument(p, "--iff_obs.y",
                  help="northing coordinate name",
                  type="character",
                  default="y")
p <- add_argument(p, "--iff_obs.z",
                  help="elevation name",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.value",
                  help="variable name",
                  type="character",
                  default="value")
p <- add_argument(p, "--iff_obs.proj4",
                  help="proj4 string for the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.sourceId",
                  help="station identifier",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.prId",
                  help="provider identifier",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.dqc",
                  help="data quality control",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# Master grid definition
p <- add_argument(p, "--grid_master.x1",
                  help="easting coordinate of the first gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.xn",
                  help="easting coordinate of the last gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.y1",
                  help="northing coordinate of the first gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.yn",
                  help="northing coordinate of the last gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.resx",
                  help="grid spacing along the easting coordinate",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.resy",
                  help="grid spacing along the northing coordinate",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.proj4",
                  help="proj4 string for the master grid",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# Mask file netcdf parameters
p <- add_argument(p, "--iff_mask.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_mask.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_mask.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# Rescaling factor file netcdf parameters
p <- add_argument(p, "--iff_rf.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_rf.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.ndim_ll",
                  help="number of dimensions for the lat-lon variables",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.tpos_ll",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_rf.names_ll",
                  help="dimension names for the lat-lon variables",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_rf.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.adjfact",
                  help="correction factor",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_rf.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
#------------------------------------------------------------------------------
# Land area fraction file netcdf parameters
p <- add_argument(p, "--iff_laf.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_laf.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_laf.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.adjfact",
                  help="correction factor (laf should be 0-1)",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_laf.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
#------------------------------------------------------------------------------
# Land area fraction file netcdf parameters
p <- add_argument(p, "--iff_dem.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_dem.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_dem.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.adjfact",
                  help="correction factor",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_dem.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
#------------------------------------------------------------------------------
# first-guess file netcdf parameters
p <- add_argument(p, "--fg_upscale",
                  help="upscale first-guess on master grid",
                  flag=T)
p <- add_argument(p, "--iff_fg.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_fg.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_fg.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.adjfact",
                  help="correction factor",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_fg.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
#------------------------------------------------------------------------------
# output file netcdf parameters
p <- add_argument(p, "--off_x.grid",
                  help="grid type",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.variables",
                  help="type of variable (analysis,background,idi,...)",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_x.varname",
                  help="variable name",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_x.varlongname",
                  help="variable long name",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.standardname",
                  help="variable standard name",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_x.varversion",
                  help="variable version",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.varunit",
                  help="variable unit",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_x.timesunit",
                  help="time unit",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.reference",
                  help="references",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.write_lonlat",
                  help="add latitude and longitude variables",
                  type="logical",
                  default=F)
p <- add_argument(p, "--off_x.diground",
                  help="rounding digits",
                  type="numeric",
                  default=3)
p <- add_argument(p, "--off_x.summary",
                  help="summary",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.sourcestring",
                  help="source string",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.title",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.comment",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_x.cell_methods",
                  help="title",
                  type="character",
                  nargs=Inf,
                  default=NA)
#------------------------------------------------------------------------------
# output file netcdf parameters
p <- add_argument(p, "--off_xb.grid",
                  help="grid type",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.variables",
                  help="type of variable (analysis,background,idi,...)",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_xb.varname",
                  help="variable name",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_xb.varlongname",
                  help="variable long name",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.standardname",
                  help="variable standard name",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_xb.varversion",
                  help="variable version",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.varunit",
                  help="variable unit",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_xb.timesunit",
                  help="time unit",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.reference",
                  help="references",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.write_lonlat",
                  help="add latitude and longitude variables",
                  type="logical",
                  default=F)
p <- add_argument(p, "--off_xb.diground",
                  help="rounding digits",
                  type="numeric",
                  default=3)
p <- add_argument(p, "--off_xb.summary",
                  help="summary",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.sourcestring",
                  help="source string",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.title",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.comment",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_xb.cell_methods",
                  help="title",
                  type="character",
                  nargs=Inf,
                  default=NA)
#------------------------------------------------------------------------------
# gaussian anamorphosis
p <- add_argument(p, "--transf",
                  help="transformation used in the gaussian anamorphosis (\"none\",\"Box-Cox\")",
                  type="character",
                  default="none")
p <- add_argument(p, "--transf.boxcox_lambda",
                  help="Box-Cox power parameter",
                  type="numeric",
                  default=NA)
#------------------------------------------------------------------------------
# Miscellaneous
p <- add_argument(p, "--maxobs_for_matrixInv",
                  help="Box-Cox power parameter",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--cores",
                  help="set the number of cores for parallel runs. Rpackage \"parallel\" required. 0 stands for \"use detectCores\". Default do not use it.",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--min_plaus_val",
                  help="minimum plausible value (e.g., 0 for prec)",
                  type="numeric",
                  default=NA)
# observation post-processing
p <- add_argument(p, "--obspp.agg",
                  help="observation post-processing, aggregation",
                  flag=T)
p <- add_argument(p, "--obspp.agg_prId",
                  help="observation post-processing, list of provider Ids (NA stands for all)",
                  type="numeric",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--obspp.agg_fact",
                  help="observation post-processing, list of aggregation factors (with respect to the master grid)",
                  type="numeric",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--obspp.dqcpuddle",
                  help="observation post-processing, data quality control puddle check",
                  flag=T)
p <- add_argument(p, "--obspp.dqcpuddle_fact",
                  help="observation post-processing, dqc puddle aggregation factor (with respect to the master grid)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--off_obspp",
                  help="output full file name for the post-processed observations. write file in TITAN format, then exit",
                  type="character",
                  default=NA)

p <- add_argument(p, "--oi_k_dim",
                  help="number of background ensemble members",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--oi_a_dim",
                  help="number of background ensemble members",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--oi_rain_uo",
                  help="rain yes/no threshold for alignment (mm)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--oi_align_mode",
                  help="strategy used for alignment ('ets','maxoverlap')",
                  type="character",
                  default="maxoverlap")
p <- add_argument(p, "--oi_rain_yo",
                  help="rain yes/no threshold for interpolation (mm)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--oi_pmax",
                  help="maximum number of observations to use in the surrounding of a grid point",
                  type="integer",
                  default=50)
p <- add_argument(p, "--oi_range",
                  help="range of allowed values",
                  type="numeric",
                  nargs=2,
                  default=c(NA,NA))
p <- add_argument(p, "--oi_dh",
                  help="horizontal decorrelation lenght scale",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--oi_eps2",
                  help="ration observation and background error variances",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--oi_corrfun",
                  help="correlation functions (\"gaussian\",\"soar\",\"toar\",\"powerlaw\")",
                  type="character",
                  default="gaussian")

p <- add_argument(p, "--wise_oi_pmax",
                  help="maximum number of observations to use in the surrounding of a grid point",
                  type="integer",
                  default=50)
p <- add_argument(p, "--wise_oi_range",
                  help="range of allowed values",
                  type="numeric",
                  nargs=2,
                  default=c(NA,NA))
p <- add_argument(p, "--wise_oi_dh",
                  help="horizontal decorrelation lenght scale",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--wise_oi_eps2",
                  help="ration observation and background error variances",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--wise_oi_corrfun",
                  help="correlation functions (\"gaussian\",\"soar\",\"toar\",\"powerlaw\")",
                  type="character",
                  default="gaussian")


#------------------------------------------------------------------------------
# gaussian anamorphosis
p <- add_argument(p, "--wise_rain_uo",
                  help="rain yes/no threshold for alignment (mm)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--wise_align_mode",
                  help="strategy used for alignment ('ets','maxoverlap')",
                  type="character",
                  default="maxoverlap")
p <- add_argument(p, "--wise_rain_yo",
                  help="rain yes/no threshold for interpolation (mm)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--wise_k_dim",
                  help="number of background ensemble members",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--wise_a_resample",
                  help="should we resample the analysis members from the posterior PDF? If yes, then use wise_a_dim to define the number of anlysis members, otherwise use wise_k_dim",
                  flag=T)
p <- add_argument(p, "--wise_a_resample_setseed",
                  help="set seed when resampling (useful when testing)",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--wise_a_dim",
                  help="number of analysis ensemble members",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--wise_wf",
                  help="wavelet type",
                  type="character",
                  default="bl20")
p <- add_argument(p, "--wise_boundary",
                  help="strategy adopted to avoid boundary effects (\"periodic\", \"reflection\"",
                  type="character",
                  default="periodic")
p <- add_argument(p, "--wise_n_levs_mx",
                  help="number of coarsest level to consider",
                  type="integer",
                  default=8)
p <- add_argument(p, "--wise_n_levs_mn",
                  help="number of finest level to consider",
                  type="integer",
                  default=1)
p <- add_argument(p, "--wise_supob_nobs",
                  help="Wise Superobbing, max number of observations considered in the neighbourhhod of a grid point",
                  type="integer",
                  default=50)
p <- add_argument(p, "--wise_supob_radius",
                  help="Wise Superobbing, radius defining the neighbourhhod of a grid point",
                  type="numeric",
                  default=1500)
p <- add_argument(p, "--wise_supob_mode",
                  help="Wise Superobbing strategy ('quantile', 'mean')",
                  type="character",
                  default="mean")
p <- add_argument(p, "--wise_supob_q",
                  help="Wise Superobbing, quantile used for superobbing",
                  type="numeric",
                  default=0.99)
p <- add_argument(p, "--wise_opt_maxit",
                  help="Wise Optimization, maximum number of iterations in the main loop",
                  type="integer",
                  default=100)
p <- add_argument(p, "--wise_opt_opttol",
                  help="Wise Optimization, tolerance threshold used to break out from the main loop",
                  type="numeric",
                  default=0.02)
p <- add_argument(p, "--wise_En2_adj_fun",
                  help="Wise adjustment function to avoid abrupt variations is the energy specturm between two consecutive iterations of the analysis loop (\"Gaussian\",\"SOAR\",\"powerlaw\",\"TOAR\")",
                  type="character",
                  default="Gaussian")
p <- add_argument(p, "--wise_En2_adj_min",
                  help="lower boundary for the adjustment function",
                  type="numeric",
                  default=0.)
p <- add_argument(p, "--wise_rescale_min_obs",
                  help="minimum number of observations in an event to apply wise",
                  type="integer",
                  default=100)
p <- add_argument(p, "--wise_rescale_min_cells",
                  help="minimum number of cells for the extent of a precipitation event to be considered significative enough to apply wise",
                  type="integer",
                  default=100)


#------------------------------------------------------------------------------
#
argv <- parse_args(p)
#
#------------------------------------------------------------------------------
# read configuration file
if (!is.na(argv$config.file)) {
  if (file.exists(argv$config.file)) {
    source(argv$config.file)
    argv_tmp<-append(argv,conf)
    names_argv_tmp<-names(argv_tmp)
    argv_def<-list()
    names_argv_def<-integer(0)
    k<-0
    for (i in 1:length(argv_tmp)) {
      if (names_argv_tmp[i] %in% names_argv_def) next
      k<-k+1
      j<-which(names_argv_tmp==names_argv_tmp[i])
      argv_def[[k]]<-argv_tmp[[j[length(j)]]]
      names_argv_def<-c(names_argv_def,names_argv_tmp[i])
    }
    names(argv_def)<-names_argv_def
    rm(argv_tmp,names_argv_tmp,names_argv_def)
    rm(argv)
    argv<-argv_def
    rm(argv_def)
  } else {
    print("WARNING: config file not found")
    print(argv$config.file)
  }
}


#
#-----------------------------------------------------------------------------
# read fg config files
if ( any( !is.na( argv$fg.files))) {
  for (f in 1:length(argv$fg.files)) {
    if ( file.exists( argv$fg.files[f])) {
      source( argv$fg.files[f], local=T)
      fg_env$fg[[f]] <- conf
      rm( conf)
    } else {
      print( "WARNING: config file not found")
      print( argv$fg.files[f])
    }
  }
}

# one can specify file names outside the fg config files
if ( any( !is.na( argv$fg.filenames))) {
  if ( length( argv$fg.files) != length( argv$fg.filenames))
    boom( "since --fg.files and --fg.filenames are both defined, they must have the same length")
  for (f in 1:length(argv$fg.files)) {
    if ( file.exists( argv$fg.filenames[f])) {
      fg_env$fg[[f]]$main.file <- argv$fg.filenames[f]
    } else {
      boom( paste( "ERROR: file not foud", argv$fg.filenames[f]))
    }
  }
}

#
#-----------------------------------------------------------------------------
# read uo config file
if ( !is.na( argv$uo.file)) {
  if ( file.exists( argv$uo.file)) {
    source( argv$uo.file, local=T)
    u_env$uo[[1]] <- conf
    rm( conf)
  } else {
    print( "WARNING: config file not found")
    print( argv$uo.file)
  }
}

# one can specify file names outside the uo config file
if ( !is.na( argv$uo.filename)) {
  if ( file.exists( argv$uo.filename)) {
    u_env$uo[[1]]$main.file <- argv$uo.filename
  } else {
    boom( paste( "ERROR: file not foud", argv$uo.filename))
  }
}

#
#-----------------------------------------------------------------------------
# set variables of the env environment
if (argv$mode=="wise") {
  env$k_dim <- argv$wise_k_dim
  env$a_dim <- argv$wise_a_dim
  if (is.na(env$a_dim) || !argv$wise_a_resample) env$a_dim <- env$k_dim
  env$wf <- argv$wise_wf
  env$boundary <- argv$wise_boundary
  env$n_levs_mx <- argv$wise_n_levs_mx
  env$n_levs_mn <- argv$wise_n_levs_mn
  u_env$rain <- argv$wise_rain_uo
  y_env$rain <- argv$wise_rain_yo
} else if (argv$mode=="oi") {
  env$k_dim <- argv$oi_k_dim
  env$a_dim <- argv$oi_a_dim
  if ( is.na(env$a_dim)) env$a_dim <- env$k_dim
  u_env$rain <- argv$oi_rain_uo
  y_env$rain <- argv$oi_rain_yo
}

#
return(argv)
}
