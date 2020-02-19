#+
`OI_RR_var`<-function(yo,
                      yb,
                      xb,
                      gx,
                      gy,
                      ox,
                      oy,
                      Dh,
                      eps2
                      ) {
#------------------------------------------------------------------------------
  no<-length(yo)
  ng<-length(gx)
  xa<-vector(mode="numeric",length=ng)
  ya<-vector(mode="numeric",length=no)
  xa_errvar<-vector(mode="numeric",length=ng)
  ya_errvar<-vector(mode="numeric",length=no)
  o_errvar<-0
  xa[]<-0
  ya[]<-0
  xa_errvar[]<-0
  ya_errvar[]<-0
  out<-.C("oi_rr_var",ng=as.integer(ng),
                      no=as.integer(no),
                      SRinv=as.numeric(InvD),
                      eps2=as.double(eps2),
                      Dh=as.double(Dh),
                      gx=as.double(gx),
                      gy=as.double(gy),
                      ox=as.double(ox),
                      oy=as.double(oy),
                      yo=as.double(yo),
                      yb=as.double(yb),
                      xb=as.double(xb),
                      xa=as.double(xa),
                      ya=as.double(ya),
                      xa_errvar=as.double(xa_errvar),
                      ya_errvar=as.double(ya_errvar),
                      o_errvar=as.double(o_errvar))
  return( list(xa= out$xa[1:ng],
               ya= out$ya[1:no],
               xa_errvar= out$xa_errvar[1:ng],
               ya_errvar= out$ya_errvar[1:no],
               o_errvar= out$o_errvar))
}

