#include <math.h>
#include<stdio.h>


/* optimal interpolation for precipitation.
 Input:
 ng = number of grid points
 no = number of observations
 (ox,oy) = observation coordinates
 (gx,gy) = grid point coordinates
 yo = observations
 xb = background at grid points
 yb = background at observation points
 Dh,eps2 = OI parameters ( eps2=o_errvar/b_errvar )
 SRinv = (S+R)^-1, where R = eps2 * I & S = f(Dh)

 Output:
 xa = analysis at grid points
 ya = analysis at observation points
 xa_errvar = analysis error variance at grid points
 ya_errvar = analysis error variance at observation points
 o_errvar = best estimate of observation error variance

Method:
xa[i] = xb[i] + sum_j{ K_ij %*% (yo[j]-yb[j]) }
K_ij = sum_l{ G_il ((S+R)^-1)_lj }
G_il background error correlation between i-th grid point and l-th observation point

xa_errvar[i] = b_errvar * [ 1 - (K%*%Gt)_ii]
(K%*%Gt)_ii = sum_j{ K_ij * Gt_ji }
            = sum_j{ K_ij * G_ij }

ya[i] = yb[i] + sum_j{ W_ij %*% (yo[j]-yb[j]) }
W_ij = sum_l{ S_il ((S+R)^-1)_lj }

ya_errvar[i] = b_errvar * [ 1 - (W%*%St)_ii]
(W%*%St)_ii = sum_j{ W_ij * St_ji }
            = sum_j{ W_ij * S_ij }

o_errvar=mean{ (yo-yb) * (yo-ya) }
b_errvar=o_errvar / eps2

*/
void oi_rr_var(int *ng, 
               int *no,
               double *SRinv,
               double *eps2,
               double *Dh, 
               double *gx,
               double *gy,
               double *ox,
               double *oy,
               double *yo, 
               double *yb, 
               double *xb, 
               double *xa,
               double *ya,
               double *xa_errvar,
               double *ya_errvar,
               double *o_errvar) {
  int i, j, k, l, m, k1, thr, no_i;
  double Dh2, hd2_ij, b_errvar;
  double KGt_ii, WSt_ii, K_ij, W_ij;
  double vec[no[0]],G_i[no[0]],S_i[no[0]];
  int sel_i[no[0]];

  Dh2= Dh[0]*Dh[0];

  /* initialize error variances (b_errvar set to 1 temporarly) */
  o_errvar[0]= 0;
  b_errvar= 1;

  /* loop over observation points */
  for (i=0;i<no[0];i++) {
    /* initialization */
    ya_errvar[i]= 0;
    ya[i]= yb[i];
    WSt_ii= 0;
    /* analysis loop */
    k= 0; /* k in the counter for the SRinv matrix */
    k1= 0; /* k in the counter for the SRinv matrix */
    for (j=0;j<no[0];j++) {
      W_ij= 0;
      for (l=0;l<no[0];l++) {
        /* compute the K_ij (W_ij) */
        if (j==0) {
          /* (a) compute the ith-row G and S elements, just once for every "i"
             (b) compute the vector "vec" of (S+R)^-1 %*% (yo-yb), just once */
          hd2_ij=((ox[i]-ox[l])*(ox[i]-ox[l])+(oy[i]-oy[l])*(oy[i]-oy[l])) / (1000.*1000.);
          S_i[l]= exp(-0.5*(hd2_ij/Dh2));
          if (i==0) {
            vec[l]=0;
            for (m=0;m<no[0];m++) {
              vec[l]= vec[l]+SRinv[k1]*(yo[m]-yb[m]);
              k1++;
            }
          }
        }
        W_ij= W_ij+SRinv[k]*S_i[l];
        k++;
      }
      /* analysis and (K%*%Gt)_ii ( (W%*%St)_ii ) */
      ya[i]= ya[i]+S_i[j]*vec[j];
      WSt_ii= WSt_ii+W_ij*S_i[j];
    } /* end of the analysis loop*/
    o_errvar[0]= o_errvar[0]+(yo[i]-yb[i])*(yo[i]-ya[i]);
    ya_errvar[i]= b_errvar*(1-WSt_ii);
  }
  /* error variances */
  /* NOTE: before i==(n[0]-1) ya_errvar are wrong by a factor of b_errvar */
  /*       they are adjusted now */
  o_errvar[0]= o_errvar[0] / no[0];
  b_errvar= o_errvar[0]/eps2[0];
  for (j=0;j<no[0];j++) {
    ya_errvar[j]= b_errvar*ya_errvar[j];
  }
  thr=7*7;
  /* loop over grid points & observation points at the same time */
  for (i=0;i<ng[0];i++) {
/*    if ( i % 10000 == 0) {
      printf("%d / %d \n", i,ng[0]);
    }*/
    no_i=0;
    for (j=0;j<no[0];j++) {
      hd2_ij=((gx[i]-ox[j])*(gx[i]-ox[j])+(gy[i]-oy[j])*(gy[i]-oy[j])) / (1000.*1000.);
      if ((hd2_ij/Dh2)<thr) {
        G_i[no_i]= exp(-0.5*(hd2_ij/Dh2));
        sel_i[no_i]=j;
        no_i++;
      }
    }
    /* initialization */
    xa[i]= xb[i];
    KGt_ii= 0;
    if (no_i>0) {
      /* analysis loop */
      k= 0; /* k in the counter for the SRinv matrix */
      for (j=0;j<no_i;j++) {
        K_ij= 0;
        for (l=0;l<no_i;l++) {
          K_ij= K_ij+SRinv[(sel_i[j]*no[0]+sel_i[l])]*G_i[l];
        }
        /* analysis and (K%*%Gt)_ii ( (W%*%St)_ii ) */
        xa[i]= xa[i]+G_i[j]*vec[sel_i[j]];
        KGt_ii= KGt_ii+K_ij*G_i[j];
      } /* end of the analysis loop*/
    }
    xa_errvar[i]= b_errvar*(1-KGt_ii);
  } /* end loop over grid points & observation points at the same time */
  return;
}
