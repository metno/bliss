#include <math.h>
#include<stdio.h>

#define PI 3.14159265

/* Update the regional temperature background field estimate (bg) over a set of ng points G with coordinates (xg,yg,zg)  given that this update is based on the observations measured at a set of no points O with coordinates (xo,yo). The set O defines a sub-region inside G. The local profile valid for the sub-region is obtained using the tpar parameters. The regional background is obtained through the incremental calculation of weighted mean over all the local backgrounds, with weigths wg that are computed by means of OI-IDI concepts (Dh and vec1). */
void oi_t_xb_upd(int *ng,
                 int *no, 
                 double *xg,
                 double *yg,
                 double *zg,
                 double *bg,
                 double *wg,
                 double *dg,
                 double *dg_new,
                 double *xo,
                 double *yo,
                 double *Dh,
                 double *vec1,
                 double *tpar,
                 double *na) {

/*  int ngv = ng[0];
  int nov = no[0]; 
  double nav = na[0];
  double Dhv = Dh[0];*/
  int i,j;
  double g,Dh2,hd2,vd2;
  double bg_new,wg_new;

  Dh2=Dh[0]*Dh[0];

/* tpar basic 0=method(1) 1=t0 2=gamma 3=na 4=na 5=na  */
/* tpar Frei  0=method(2) 1=t0 2=gamma 3= a 4=h0 5=h1i */
  if (tpar[5]<0) {
    tpar[5]=-tpar[5];
  }
  /* loop over G */
  for (i=0;i<ng[0];i++) {
    /* temperature at zg[i] */
    if (tpar[0]==1) {
      bg_new=tpar[1]+tpar[2]*zg[i];
    } else {
      if (zg[i]<=tpar[4]) {
        bg_new=tpar[1]+tpar[2]*zg[i]-tpar[3];
      } else if (zg[i]>=(tpar[4]+tpar[5])) {
        bg_new=tpar[1]+tpar[2]*zg[i];
      } else {
        bg_new=tpar[1]+tpar[2]*zg[i]-tpar[3]/2*(1+cos(PI*(zg[i]-tpar[4])/tpar[5]));
      }
    }
    /* loop over O, get the current weight */
    wg_new=0;
    for (j=0;j<no[0];j++) { 
      hd2=((xg[i]-xo[j])*(xg[i]-xo[j])+(yg[i]-yo[j])*(yg[i]-yo[j]));
      if (hd2<(49*Dh2)) {
        g=exp(-0.5*(hd2/Dh2));
        wg_new=wg_new+g*vec1[j];
      }
    }
    /* incremental calculation of the weighted mean */
    if (bg[i]==na[0]) {
      wg[i]=wg_new;
      bg[i]=bg_new;
      dg[i]=dg_new[0];
    } else if (wg_new>0.) {
      wg[i]=wg[i]+wg_new;
      bg[i]=bg[i]+wg_new/wg[i]*(bg_new-bg[i]);
      dg[i]=dg[i]+wg_new/wg[i]*(dg_new[0]-dg[i]);
    }
  }
  return;
}
