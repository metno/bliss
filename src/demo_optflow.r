source("~/projects/bliss/src/optflow_util.r")
source("~/projects/bliss/src/optflow_msa.r")
x <- y <- matrix(0, 100, 100)
x[11:20,31:40] <- 1
x[61:70,61:70] <- 1
y[15:25,41:51] <- 1
#y[11:21,41:51] <- 1
y[51:61,61:71] <- 1
r1<-r2<-raster(ncol=100,nrow=100)
r1[]<-x
r2[]<-y
#res<-optical_flow_HS(r1,r2,niter=1000,nlevel=4,w1=1,w2=1)
res<-optical_flow_HS(r1,r2,niter=1000,nlevel=4,w1=0.25,w2=0.03)

rfig <- r1+r2
rmod <- warp( r1, -res$u, -res$v)
image(rfig+rmod)
plot_arrows(res$u, res$v, fact=0.2, length=0.03)
