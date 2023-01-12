#+
read_and_regrid_nc <- function( nc.file,
                                nc.varname,
                                topdown,
                                out.dim,
                                proj4,
                                nc.proj4,
                                selection,
                                adjfact,
                                adjval,
                                rmaster,
                                nc.varname_lat="none",
                                nc.varname_lon="none",
                                out.dim_ll,
                                upscale=F,
                                upscale_fun="mean",
                                projectraster_method="bilinear") {
#------------------------------------------------------------------------------
  grid_master.proj4 <- projection(rmaster, asText=TRUE)
  out.dim$epos    <- set_NAs_to_NULL(out.dim$epos)
  out.dim$tpos    <- set_NAs_to_NULL(out.dim$tpos)
  out.dim_ll$tpos <- set_NAs_to_NULL(out.dim_ll$tpos)
  selection$e     <- set_NAs_to_NULL(selection$e)
  if ( any( out.dim$names == "time")) {
    if ( selection$t == "none") selection$t <- nc4.getTime( nc.file)[1]
  } else {
    selection$t <- NULL
  }
  raux <- try(read_dotnc(nc.file=nc.file,
                         nc.varname=nc.varname,
                         topdown=topdown,
                         out.dim=list(ndim=out.dim$ndim,
                                      tpos=out.dim$tpos,
                                      epos=out.dim$epos,
                                      names=out.dim$names),
                         proj4=proj4,
                         nc.proj4=nc.proj4,
                         selection=list(t=selection$t,
                                        format=selection$format,
                                        e=selection$e)))
  if (is.null(raux)) boom(paste("error reading",nc.file))
  if (adjfact!=1 | adjval!=0) 
    raux$stack[]<-getValues(raux$stack)*adjfact+adjval
  r<-raux$stack; rm(raux)
  val<-getValues(r)
  if (!rasters_match(r,rmaster)) {
    if (nc.varname_lat=="none") {
      # Upscale to coarser grid
      if (upscale) {
        ix<-which(!is.na(val))
        if (proj4==grid_master.proj4) {
          coord.new<-xyFromCell(r,ix)
        } else {
          coord.new<-spTransform( 
                      SpatialPoints(xyFromCell(r,ix),
                                     proj4string=CRS(proj4)) 
                                                ,CRS(grid_master.proj4))
        }
        if (upscale_fun=="mean") {
          r1<-rasterize(x=coord.new, y=rmaster, field=val[ix], fun=mean)
        } else if (upscale_fun=="max") {
          r1<-rasterize(x=coord.new, y=rmaster, field=val[ix], fun=max)
        } else if (upscale_fun=="min") {
          r1<-rasterize(x=coord.new, y=rmaster, field=val[ix], fun=min)
        }
        r<-r1
        rm( r1, coord.new)
        if (!any(!is.na(val<-getValues(r)))) 
          cat( paste( "warning: all NAs after upscale\n"))
      } else {
        if ( projection(r) != projection(rmaster)) {
          r <- projectRaster(r,rmaster,method=projectraster_method)
        } else {
          r <- resample(r,rmaster,method=projectraster_method)
        }
        val<-getValues(r)
      }
    } else {
      raux <- try(read_dotnc(nc.file=nc.file,
                             nc.varname=nc.varname_lat,
                             topdown=topdown,
                             out.dim=out.dim_ll,
                             proj4=proj4.llwgs84,
                             nc.proj4=list(var=NULL, att=NULL),
                             selection=list(t=selection$t,
                                            format=selection$format,
                                            e=NULL)))
      if (is.null(raux)) boom(paste("error reading (lat) ",nc.file))
      lat<-raux$data; rm(raux)
      raux <- try(read_dotnc(nc.file=nc.file,
                             nc.varname=nc.varname_lon,
                             topdown=topdown,
                             out.dim=out.dim_ll,
                             proj4=proj4.llwgs84,
                             nc.proj4=list(var=NULL, att=NULL),
                             selection=list(t=selection$t,
                                            format=selection$format,
                                            e=NULL)))
      if (is.null(raux)) boom(paste("error reading (lon) ",nc.file))
      lon<-raux$data; rm(raux)
      coord.new<-spTransform(SpatialPoints(cbind(lon,lat),
                                           proj4string=CRS(proj4.llwgs84)),
                             CRS(grid_master.proj4))
      #
      val<-getValues(r)
      ragg<-rasterize(coord.new,
                        aggregate(rmaster,fact=4),
                        val)
      r<-mask( crop( disaggregate(ragg,fact=4,method="bilinear"),
                       rmaster),
                 rmaster)
      val<-getValues(r)
    } # end if read lat lon
  } # end if raster does not match rmaster
  return( list(raster = r, values = val))
}

