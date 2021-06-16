library(waveslim)

xmn<-as.numeric(argv$grid_master.x1)-as.numeric(argv$grid_master.resx)/2
xmx<-as.numeric(argv$grid_master.xn)+as.numeric(argv$grid_master.resx)/2
ymn<-as.numeric(argv$grid_master.y1)-as.numeric(argv$grid_master.resy)/2
ymx<-as.numeric(argv$grid_master.yn)+as.numeric(argv$grid_master.resy)/2

n<-ceiling(log(max(nx,ny),2))

r<-rmaster

rdyad <- raster( extent( xmn, xmx, ymn, ymx), ncol=2**n, nrow=2**n, crs=argv$grid_master.proj4)

dwt<-list()

xb[xb<0] <- 0
for (i in 1:10) {
  print(i)
  r[]<-xb[,i]
  s<-resample(r,rdyad,method="bilinear")
  s[s<0] <- 0
  dwt[[i]]<-dwt.2d(as.matrix(s),wf="haar",J=12)
}

save.image("tmp.rdata")

for (i in 1:length(dwt[[1]])) {
  vec <- array( data=NA, dim=c(length(dwt[[1]][[i]]),10) )
  for (j in 1:10) vec[,j] <- as.vector(dwt[[j]][[i]])
  var <- mean( apply( vec, MAR=1, FUN= var))
  mu <- abs( mean( rowMeans(vec)))
  print( paste( i, round(var,3), round(mu,3)))
}

save.image("tmp.rdata")

q()
