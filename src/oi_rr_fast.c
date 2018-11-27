#include <math.h>
#include<stdio.h>


/* optimal interpolation for precipitation.
 ng = number of grid points
 no = number of observations
 (xo,yo,zo) = observation coordinates and elevations
 (xg,yg,zg) = grid point coordinates and elevations
 xb = background at grid point
 vec = (S+R)^-1 * (yo-yb)
 Dh,Dz = OI parameters
 xa = analysis over the grid points
*/
void oi_rr_fast(int *ng, 
                int *no, 
                double *xg,
                double *yg,
                double *zg,
                double *xo,
                double *yo,
                double *zo,
                double *Dh, 
                double *Dz,
                double *xb, 
                double *vec,
                double *xa,
                double *zero) {
  int i,j;
  double g,Dh2,Dz2,hd2,vd2;

  Dh2=Dh[0]*Dh[0];
  Dz2=Dz[0]*Dz[0];
  for (i=0;i<ng[0];i++) {
    xa[i]=xb[i];
    for (j=0;j<no[0];j++) { 
      hd2=((xg[i]-xo[j])*(xg[i]-xo[j])+(yg[i]-yo[j])*(yg[i]-yo[j])) / (1000.*1000.);
      if (hd2<(49*Dh2)) { 
        vd2=(zg[i]-zo[j])*(zg[i]-zo[j]);
        g=exp(-0.5*(hd2/Dh2+vd2/Dz2));
        xa[i]=xa[i]+g*vec[j];
      }
    }
    if (xa[i]<zero[0]) xa[i]=zero[0]; 
  }

  return;
}
