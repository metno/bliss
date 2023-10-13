#+
argparser<-function( env, fg_env, y_env, u_env) {

cat("Read command line arguments...")

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
                  help="file names of the first-guess nc-files. It is an optional argument, if specified: i) must have the same length of --fg.files ii) it overrides the main.file arguments of the fg.files",
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
p <- add_argument(p, "--oi2step_superobbing",
                  help="superobbing (not yet implemented, originally was used for \"oi_twostep_senorge_temperature\")",
                  flag=T)
#------------------------------------------------------------------------------
# statistical interpolation mode
p <- add_argument(p, "--mode",
                  help="statistical interpolation scheme (\"rasterize\",\"oi_multiscale_senorge_prec\",\"oi_twostep_senorge_temperature\",\"successivecorrections\",\"ensi\",\"ensigap\", \"corensi\", \"oi\")",
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
# preproc - mergeobs
p <- add_argument(p, "--preproc_mergeobs",
                  help="preprocess the observations yes/no",
                  type="logical",
                  default=F)
p <- add_argument(p, "--mergeobs_eps2",
                  help="eps2 for mergeobs",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--mergeobs_pmax",
                  help="max number of observations for mergeobs",
                  type="integer",
                  default=30)
p <- add_argument(p, "--mergeobs_dh",
                  help="horizontal decorrelation length scale for mergeobs",
                  type="numeric",
                  default=3)
p <- add_argument(p, "--mergeobs_corrfun",
                  help="correlation function for mergeobs (gaussian, soar, powerlaw, toar)",
                  type="character",
                  default="toar")
p <- add_argument(p, "--mergeobs_range",
                  help="range check for mergeobs",
                  type="numeric",
                  nargs=2,
                  default=c(NA,NA))
#------------------------------------------------------------------------------
# preproc - selensemblemble
p <- add_argument(p, "--selens_mode",
                  help="selection of ensemble members (ets, maxoverlap, identity)",
                  type="character",
                  default="identity")
#------------------------------------------------------------------------------
p <- add_argument(p, "--obs_perturb",
                  help="should we perturb the observations?",
                  flag=T)
p <- add_argument(p, "--obs_perturb_rmin",
                  help="Perturb the observations with random factors from a uniform distribution",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--obs_perturb_rmax",
                  help="Perturb the observations with random factors from a uniform distribution",
                  type="numeric",
                  default=1.5)

#------------------------------------------------------------------------------
# parameters common to several interpolation methods
p <- add_argument(p, "--k_dim",
                  help="number of background ensemble members",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--k_dim_corr",
                  help="number of background ensemble members",
                  type="integer",
                  default=NA)
p <- add_argument(p, "--rain_uo",
                  help="rain yes/no threshold for gridded observations (mm)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--rain_yo",
                  help="rain yes/no threshold for in-situ observations (mm)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--range",
                  help="range check",
                  type="numeric",
                  nargs=2,
                  default=c(NA,NA))
p <- add_argument(p, "--corrfun",
                  help="correlation functions (\"gaussian\",\"soar\",\"toar\",\"powerlaw\")",
                  type="character",
                  default="gaussian")
p <- add_argument(p, "--pmax",
                  help="maximum number of observations in the neighbourhood of a gridpoint for OI",
                  type="numeric",
                  default=200)
p <- add_argument(p, "--analysis_range",
                  help="range of allowed values for the analysis",
                  type="numeric",
                  nargs=2,
                  default=c(NA,NA))
p <- add_argument(p, "--area_small_clumps",
                  help="max area (units m**2) of those small clumps of connected cells that we want to remove from the analysis (e.g. small precipitation events, most likely interpolation artifacts)",
                  type="numeric",
                  default=6250000)
p <- add_argument(p, "--alpha",
                  help="parameter used to weight the static and dynamical contributions in the definition of the hybrid correlations",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--dh",
                  help="horizontal decorrelation lenght scale",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--dhloc",
                  help="horizontal decorrelation lenght scale (used for localization)",
                  type="numeric",
                  default=100000)
p <- add_argument(p, "--dh_idi",
                  help="horizontal decorrelation lenght scale used only for IDI calculations",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--eps2",
                  help="ratio between observation and background error variances",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--use_relativeAnomalies",
                  help="use relative anomalies (observation/background)",
                  flag=T)
p <- add_argument(p, "--nSCloops",
                  help="number of loops used in successive corrections",
                  type="integer",
                  default=100)
#------------------------------------------------------------------------------
# MSA
p <- add_argument(p, "--msa_ididense",
                  help="MSA IDI threshold for defining data dense regions (IDI is defined with respect to preproc-mergeobs OI",
                  type="numeric",
                  default=0.8)
p <- add_argument(p, "--msa_eps2",
                  help="MSA ratio between observation and background error variances (0.1-1)",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--dh_vec_min",
                  help="define a sequence of spatial scales, this is the smallest scale",
                  type="numeric",
                  default=2500)
p <- add_argument(p, "--dh_vec_max",
                  help="define a sequence of spatial scales, this is the largest scale",
                  type="numeric",
                  default=200000)
p <- add_argument(p, "--dh_vec_by",
                  help="define a sequence of spatial scales, this is the largest scale",
                  type="numeric",
                  default=5000)
p <- add_argument(p, "--dh_r_max",
                  help="maximum value of the horizontal decorrelation lenght scale used in the Multiscal alignement procedure, analysis step",
                  type="numeric",
                  default=10000)
p <- add_argument(p, "--dh_obs_fact",
                  help="factor used to optimize calculation times. Given the value of Dh, the grid is defined with resolution of Dh/dh_obs_fact and it is used in the Multiscal alignement procedure, analysis step",
                  type="numeric",
                  default=3)
p <- add_argument(p, "--delta_sbe_mean_r",
                  help="SBE threshold used to stop the optimization. The score used in the optimization is the SBE computed over the sequence of spatial scales defined by dh_vec_min, dh_vec_max, dh_vec_by. Specifically, the Multiscal alignement optimization stops when the difference between two consecutive optimization iteration is less than delta_sbe_mean_r.",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--sbe_filter",
                  help="number of time steps used to smooth the SBE vector (e.g. 3 means that the SBE used for optimization at the i-th spatial scale is mean(sbe(i-1),sbe(i),sbe(i+1)).",
                  type="integer",
                  default=3)
p <- add_argument(p, "--m_step",
                  help="The analysis is processed in batches of m_steps points. This parameter can be adjusted to keep memory usage within acceptable limits.",
                  type="integer",
                  default=10000)
#------------------------------------------------------------------------------
# MSA-EnSI
p <- add_argument(p, "--msaensi_ididense",
                  help="MSA-EnSI IDI threshold for defining data dense regions (IDI is defined with respect to preproc-mergeobs OI",
                  type="numeric",
                  default=0.8)
p <- add_argument(p, "--msaensi_eps2",
                  help="MSA-EnSI ratio between observation and background error variances (0.1-1)",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--msaensi_wf",
                  help="MSA-EnSI name of the wavelet filter to use in the decomposition",
                  type="character",
                  default="haar")
p <- add_argument(p, "--msaensi_boundary",
                  help="MSA-EnSI boundary",
                  type="character",
                  default="periodic")
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

# additional OI parameters two-step temperature
# Background parameters
p <- add_argument(p, "--oi2step.bg_nrnc",
                  help="nrow ncol (i.e. \"nr\" number_of_rows, \"nc\" number_of_columns) used to define grid of sub-regional centroids (each node of this grid is a candidate centroid)",
                  type="integer",
                  nargs=2,
                  default=c(5,5))
p <- add_argument(p, "--oi2step.bg_centroids_buffer",
                  help="this is the radius of a circular region (m) around each candidate centroid where we check if (a) there is at least one grid point of the master grid that is not masked (b) there are at least \"oi2step.bg_centroids_nobsmin\" observations. If the two conditions hold true, then we have found a suitable centroid",
                  type="numeric",
                  default=100000)
p <- add_argument(p, "--oi2step.bg_centroids_nobsmin",
                  help="see the help for \"oi2step.bg_centroids_buffer\"",
                  type="integer",
                  default=1)
p <- add_argument(p, "--oi2step.bg_vertprof_gamma",
                  help="standard value for the moist adiabatic temperature lapse rate dT/dz (degC/m)",
                  type="numeric",
                  default=-0.0065)
p <- add_argument(p, "--oi2step.bg_vertprof_vmin",
                  help="minimum allowed value [units of the variable specified]",
                  type="numeric",
                  default=-50)
p <- add_argument(p, "--oi2step.bg_vertprof_vmax",
                  help="maximum allowed value [units of the variable specified]",
                  type="numeric",
                  default=40)
p <- add_argument(p, "--oi2step.bg_vertprof_dzmin",
                  help="elevation difference. If the elevation range (=elev 95th-perc minus elev 5th-perc) is less than this value, then fit a linear profile of temperature. Otherwise, fit a non-linear profile.",
                  type="numeric",
                  default=30)
p <- add_argument(p, "--oi2step.bg_blending_corr",
                  help="correlation function used when blending the sub-regional profiles",
                  type="character",
                  default="soar")
p <- add_argument(p, "--oi2step.bg_blending_dh",
                  help="horizontal decorellation length scale (m) used when blending the sub-regional profiles",
                  type="numeric",
                  default=100000)
p <- add_argument(p, "--oi2step.loop_deltam",
                  help="number of gridpoint simultaneously used when computing distances (optimization of memory usage and computing time)",
                  type="integer",
                  default=10000)
p <- add_argument(p, "--oi2step.analysis_dz",
                  help="vertical decorrelation length scale(s) used in the analysis loop",
                  type="numeric",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--oi2step.analysis_eps2",
                  help="ratio(s) between obs error variance and backg error variance used in the analysis loop",
                  type="numeric",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--oi2step.analysis_corr",
                  help="correlation function(s) used in the analysis loop",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--oi2step.analysis_lafmin",
                  help="minimum correlation value for the LAF correlation function used in the analysis loop",
                  type="numeric",
                  default=0)
p <- add_argument(p, "--oi2step.analysis_nclose",
                  help="number of closest observations (to each analysis point) used in the analysis loop",
                  type="integer",
                  default=50)
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
p <- add_argument(p, "--iff_obs.dqc_flags_to_keep",
                  help="data quality control flags to keep",
                  type="integer",
                  nargs=Inf,
                  default=NA)
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

#------------------------------------------------------------------------------
# Change-of-Resolution Ensemble Kalman Smoother

p <- add_argument(p, "--corensi_ididense",
                  help="threshold used to define dense station regions",
                  type="numeric",
                  default=0.8)
p <- add_argument(p, "--corensi_eps2",
                  help="ratio between observation and background error variances (0.1-1)",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--corensi_alpha",
                  help="parameter used to weight the static and dynamical contributions in the definition of the hybrid correlations within corensi_up sweep",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--corensi_k_dim_corr",
                  help="number of ensemble members considered for the computation of the correlations within corensi_up sweep (default is k_dim)",
                  type="integer",
                  default=NA)

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
      boom( paste( "ERROR: file not found", argv$fg.filenames[f]))
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
y_env$rain <- NA
if ( argv$mode %in% c( "corensi", "msa", "msaensi_dev", "msaensi", "oi", "ensi")) {
  env$k_dim <- argv$k_dim
  u_env$rain <- argv$rain_uo
  y_env$rain <- argv$rain_yo
  if ( is.na(argv$k_dim_corr)) argv$k_dim_corr <- argv$k_dim
  if (argv$mode == "corensi") {
    if ( is.na(argv$corensi_k_dim_corr)) argv$corensi_k_dim_corr <- argv$k_dim
  }
} else if (argv$mode %in% c("oi_multiscale_senorge_prec")) {
  y_env$rain <- argv$rain_yo
}

#
cat("ok\n")
return(argv)
}
