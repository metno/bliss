  if (!exists("r.list")) r.list<-list()
  r<-rmaster; r[]<-NA
  for (i in 1:length(argv$off_x.variables)) {
    r.list[[i]]<-matrix(data=getValues(r),
                        ncol=length(y),
                        nrow=length(x))
  }
  # define time for output
  tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                    format="%Y%m%d%H%M",tz="GMT")
  time_bnds<-array(format(rev(seq(strptime(argv$date_out,argv$date_out_fmt),
                                           length=2,by=argv$time_bnds_string)),
                   format="%Y%m%d%H%M",tz="GMT"),dim=c(1,2))
  out<-write_dotnc(grid.list=r.list,
                   file.name=argv$off_x,
                   grid.type=argv$off_x.grid,
                   x=x,
                   y=y,
                   var.name=argv$off_x.varname,
                   var.longname=argv$off_x.varlongname,
                   var.standardname=argv$off_x.varstandardname,
                   var.version=argv$off_x.varversion,
                   var.unit=argv$off_x.varunit,
                   times=tstamp_nc,
                   times.unit=argv$off_x.timesunit,
                   reference=argv$off_x.reference,
                   proj4.string=argv$grid_master.proj4,
                   lonlat.out=argv$off_x.write_lonlat,
                   round.dig=argv$off_x.diground,
                   summary=argv$off_x.summary,
                   source.string=argv$off_x.sourcestring,
                   title=argv$off_x.title,
                   comment=argv$off_x.comment,
                   atts.var.add=NULL,
                   var.cell_methods=argv$off_x.cell_methods,
                   time_bnds=time_bnds,
                   cf_1.7=T)
  print(paste("output saved on file",argv$off_x))

