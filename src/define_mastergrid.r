#+ Define master grid
define_mastergrid <- function( argv, env) {

#------------------------------------------------------------------------------

xmn<-as.numeric(argv$grid_master.x1)-as.numeric(argv$grid_master.resx)/2
xmx<-as.numeric(argv$grid_master.xn)+as.numeric(argv$grid_master.resx)/2
ymn<-as.numeric(argv$grid_master.y1)-as.numeric(argv$grid_master.resy)/2
ymx<-as.numeric(argv$grid_master.yn)+as.numeric(argv$grid_master.resy)/2
rmaster<-raster(extent(xmn,xmx,ymn,ymx),
                res=c(as.numeric(argv$grid_master.resx),
                      as.numeric(argv$grid_master.resy)),
                crs=argv$grid_master.proj4)
rmaster[]<-1
# use mask if provided
if (file.exists(argv$iff_mask)) {
  argv$iff_mask.epos<-set_NAs_to_NULL(argv$iff_mask.epos)
  argv$iff_mask.tpos<-set_NAs_to_NULL(argv$iff_mask.tpos)
  argv$iff_mask.e<-set_NAs_to_NULL(argv$iff_mask.e)
  if (!is.null(argv$iff_mask.tpos) & argv$iff_mask.t=="none") 
    argv$iff_mask.t<-nc4.getTime(argv$iff_mask)[1]
  raux<-try(read_dotnc(nc.file=argv$iff_mask,
                       nc.varname=argv$iff_mask.varname,
                       topdown=argv$iff_mask.topdown,
                       out.dim=list(ndim=argv$iff_mask.ndim,
                                    tpos=argv$iff_mask.tpos,
                                    epos=argv$iff_mask.epos,
                                    names=argv$iff_mask.names),
                       proj4=argv$iff_mask.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_mask.t,
                                      format=argv$iff_mask.tfmt,
                                      e=argv$iff_mask.e)))
  if (is.null(raux)) 
    boom("error reading the mask file")
  rmask<-raux$stack; rm(raux)
  if (!rasters_match(rmask,rmaster)) rmask<-projectRaster(rmask,rmaster)
  rmaster<-mask(rmaster,rmask); rm(rmask)
}
#
nx<-ncol(rmaster)
ny<-nrow(rmaster)
xmn<-xmin(rmaster)
xmx<-xmax(rmaster)
ymn<-ymin(rmaster)
ymx<-ymax(rmaster)

#
# extract all the cell values: zvalues[1] contains the rmaster[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues<-getValues(rmaster)
storage.mode(zvalues)<-"numeric"
xy<-xyFromCell(rmaster,1:ncell(rmaster))
x<-sort(unique(xy[,1]))
y<-sort(unique(xy[,2]),decreasing=T)
xgrid<-xy[,1]
ygrid<-xy[,2]
mask<-which(!is.na(zvalues))
ng<-length(xgrid)
ngrid<-length(mask)

# clean memory
rm(zvalues,xy)

# debug info
if (argv$verbose) {
  cat( "+---------------------------------------------------------------+\n")
  cat( "+ master grid parameters\n")
  cat( paste("nx ny dx dy", as.integer(nx), as.integer(ny),
       round(xres(rmaster),2), round(yres(rmaster),2),"\n"))
  cat( paste("xmn xmx ymn ymx", round(xmn,2), round(xmx,2),
       round(ymn,2), round(ymx,2),"\n"))
  cat( paste("# grid points=",as.integer(ng),"\n"))
  cat( paste("# unmasked grid points=",as.integer(ngrid),"\n"))
}

env$rmaster <- rmaster
env$x <- x
env$y <- y
env$xgrid <- xgrid
env$ygrid <- ygrid
env$mask <- mask
env$ng <- ng
env$ngrid <- ngrid
env$nx <- nx
env$ny <- ny
env$xmn <- xmn
env$xmx <- xmx
env$ymn <- ymn
env$ymx <- ymx

}
