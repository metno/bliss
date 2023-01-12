#
debug_plots<-function() {
  png(file=paste0("rb_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  image(rb,breaks=c(0,seq(0.001,0.8,length=25)),
        col=c("gray",rev(rainbow(24))))
  points(VecX,VecY,pch=19,col="black",cex=0.5)
  dev.off()
  png(file=paste0("innoh_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  d<-yo_relan[ixwet]-yb[ixwet]
  hist(d,breaks=seq(-1000.025,1000,by=.05),xlim=c(-1,1))
  dev.off()
  png(file=paste0("plot_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  plot(yo_relan,yb,pch=19,col="blue",xlim=c(0,1),ylim=c(0,1))
  ix9<-which(prId==9)
  points(yo_relan[ix9],yb[ix9],pch=19,col="red")
  lines(-1000:1000,-1000:1000,col="red",lwd=2)
  dev.off()
}

#
debug_01_plots<-function() {
  png(file=paste0("ra.png"),width=800,height=800)
  image(ra,breaks=c(0,seq(0.1,50,length=45)),col=c("gray",rev(rainbow(44))))
  dev.off()
  png(file=paste0("ra1.png"),width=800,height=800)
  image(ra,breaks=c(0,0.09999,1000),col=c("gray","blue"))
  points(VecX[ixwet],VecY[ixwet],pch=19,col="cyan",cex=0.5)
  points(VecX[ixdry],VecY[ixdry],pch=19,col="red",cex=0.5)
  dev.off()
}

