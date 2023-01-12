#include <math.h>
#include <stdio.h>


void obsop_LapseRateConst(int *s, int *m,
                  double *xm, double *ym, double *zm,
                  double *xs, double *ys, double *zs,
                  double *ifield,
                  double *oval,
                  double *mMinElevDiff,
                  double *mMinGradient,
                  double *mMaxGradient,
                  double *mSearchRadius) {
/*=============================================================================
LETKF matrix multiplications
Input
 s: integer. number of (all the available) observations;
 xs;ys;zs: double vectors [s]. station x-y coord (m) and elevation (m)
 m: integer. numbver of gridpoints.
 xm;ym;zm: double vectors [m]. x-y coord (m) and elevation (m)
 ifield: double vector [m]. field on a regular grid.
 mMinElevDiff: double. 
 mMinGradient: double. 
 mMaxGradient: double. 
 mSearchRadius: double.
Output
 oval: double vector [s]. output.
-------------------------------------------------------------------------------
mat:
 1,1 1,2 1,3 ... 1,k
 2,1 2,2 2,3 ... 2,k
 ...
 s,1 s,2 s,3 ... s,k
vec:
l      0   1   2      s-1  s  s+1 s+2   2*s-1 ...3*s-1 ...k*s+-1
(i,j) 1,1 2,1 3,1 ... s,1 1,2 2,2 3,2 ... s,2 ... s,3  ... s,k
mat[i,j] -> vec[l], l=(j-1)*s+(i-1)
=============================================================================*/
int s0 = s[0];
int m0 = m[0];
double Hdist,Zdist;
int s1,m1;
int counter = 0;
double MV = -9999.;
double meanXY  = 0.; // elev*T
double meanX   = 0.; // elev
double meanY   = 0.; // T
double meanXX  = 0.; // elev*elev
double min = MV;
double max = MV;
double mMinElevDiff0=mMinElevDiff[0];
double mMinGradient0=mMinGradient[0];
double mMaxGradient0=mMaxGradient[0];
double mConstantGradient;
//double HdistMin;
//double nearestElev;
double mSearchRadius0=mSearchRadius[0];
double gradient;
double elevDiff;
double dElev;
double x,dx;
double y,dy;

//  mMinElevDiff = 30; // m
  mConstantGradient  = -0.0065;    // K/m
  for (s1=0;s1<s0;s1++) {
    /* Compute the model's gradient:
       The gradient is computed by using linear regression on forecast ~ elevation
       using all forecasts within a neighbourhood. To produce stable results, there
       is a requirement that the elevation within the neighbourhood has a large
       range (see mMinElevDiff).

       For bounded variables (e.g. wind speed), the gradient approach could cause
       forecasts to go outside its domain (e.g. negative winds). If this occurs,
       the nearest neighbour is used.
    */
    meanXY  = 0; // elev*T
    meanX   = 0; // elev
    meanY   = 0; // T
    meanXX  = 0; // elev*elev
    counter = 0;
    min = MV;
    max = MV;
//    HdistMin = MV;
    for (m1=0;m1<m0;m1++) {
      dx=fabs(xm[m1]-xs[s1])/1000.;
      dy=fabs(ym[m1]-ys[s1])/1000.;
      if (dx>mSearchRadius0 || dy>mSearchRadius0) continue;
      Hdist=pow( dx*dx + dy*dy,0.5); /* Km */
      Zdist=fabs(zm[m1]-zs[s1]); /* m */
      if (Hdist<=mSearchRadius0) {
//        if (HdistMin==MV || Hdist<HdistMin)
//          nearestElev = zm[m1];
        // Y=field; X=elev
        x = zm[m1];
        y = ifield[m1];
//        if (fabs(y)>100) printf("@@ %d %7.2f %7.2f \n",m1,x,y);
        meanXY += x*y;
        meanX  += x;
        meanY  += y;
        meanXX += x*x;
        counter++;
        // Found a new min
        if(min==MV || x < min)
           min = x;
        // Found a new max
        if(max==MV || x > max)
           max = x;
      }
    } // ==> end for grid points
    if (counter==0) {
      oval[s1] = MV;
//      printf("%7.2f\n",oval[s1]);
    } else {
      gradient = MV;
      // Compute elevation difference within neighbourhood
      elevDiff = MV;
      if (min!=MV && max!=MV) 
        elevDiff = max - min;
//      printf("%7.2f \n",elevDiff);
      // Use model gradient if:
      // 1) sufficient elevation difference in neighbourhood
      // 2) regression parameters are stable enough
      meanX  /= counter;
      meanY  /= counter;
      if (counter > 0 && elevDiff!=MV && elevDiff >= mMinElevDiff0 && meanXX != meanX*meanX) {
        // Estimate lapse rate
        meanXY /= counter;
        meanXX /= counter;
        gradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
      }
//      printf("%7.2f %7.2f %10.7f  \n",meanX,meanY,gradient);
      // Check against minimum and maximum gradients
      if(gradient==MV || gradient < mMinGradient0)
        gradient = mConstantGradient;
      if(gradient==MV || gradient > mMaxGradient0)
        gradient = mConstantGradient;
      dElev = zs[s1] - meanX;
      oval[s1] = meanY + dElev * gradient;
//      if (fabs(oval[s1])>50) printf("@@@@ %7.2f %7.2f %7.2f %10.7f \n",oval[s1],meanY,dElev,gradient);
    }
//    printf("---> %d %d %10.7f %7.2f \n",s1,counter,gradient,oval[s1]);
  } // end for stations
  return;
}
