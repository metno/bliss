  if (nwet>0) {
    # define sequence of nearest station values
    if (n0<100) {
      kseq<-rev(unique(c(2,3,4,seq(5,n0,by=5),n0)))
    } else {
      kseq<-rev(unique(c(2,3,4,seq(5,100,by=5),seq(100,n0,by=200),n0)))
    }
    # average distance between a station and its nearest kseq[1], kseq[2],...
    # NOTE: dobs_fun with k=1 returns 0 km, so k=2 is the distance to closest station
    vecd_tmp<-vector()
    vecf<-vector()
    for (i in 1:length(kseq)) 
      vecd_tmp[i]<-round(dobs_fun(obs=data.frame(x=VecX,y=VecY),k=kseq[i])/1000,0)
    vecd_tmp<-unique(sort(c(vecd_tmp,2:100),decreasing=T))
    vecf_tmp<-pmin(min(c(nx,ny))/3,pmax(1,round(vecd_tmp/2,0)))
    vecd<-vecd_tmp[which(!duplicated(vecf_tmp,fromLast=T))]
    vecd<-vecd_tmp[which(!duplicated(vecf_tmp,fromLast=F))]
    vecd<-unique(sort(c(vecd_tmp[which(!duplicated(vecf_tmp,fromLast=T) & vecd_tmp>=2)],
                        vecd_tmp[which(!duplicated(vecf_tmp,fromLast=F) & vecd_tmp>=2)]),
                      decreasing=T))
    # aggregation factor is half the horizontal decorrelation length
    vecf<-round(vecd/2,0)
    # same weight to the obs and the background from previous iteration
    vece<-rep(argv$eps2,length(vecf))
    nl<-length(vecd)
    vece[]<-1
    rm(vecf_tmp,vecd_tmp)
    if (argv$verbose) {
      print("+---------------------------------------------------------------+")
      print("vecd vecf vece")
      print(" km   -    -")
      print(cbind(vecd,vecf,vece))
      print("+---------------------------------------------------------------+")
    }
  }

