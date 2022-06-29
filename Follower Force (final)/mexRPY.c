#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define PI 3.14159265358979323846264338327950288419716939937510

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  // Input declarations

  int N;
  double *mu;
  double *X, *Forces, *Torques, *R;

  // Output declarations

  double *V, *Omega;

  // Internal variable declarations

  int m, n;
  double am, an, amax, amin, r, a0, a1, a2, dot, A, B, heaviside_val, fx, fy,
  fz, tx, ty, tz, rhatx, rhaty, rhatz, xm, ym, zm, xn, yn, zn, vx, vy, vz, wx, wy, wz;

  // Assign input and output

  X = mxGetPr(prhs[0]);
  Forces = mxGetPr(prhs[1]);
  Torques = mxGetPr(prhs[2]);
  R = mxGetPr(prhs[3]);
  mu = mxGetPr(prhs[4]);
  N = mxGetN(prhs[0])/3;

  plhs[0] = mxCreateDoubleMatrix(1,3*N,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,3*N,mxREAL);

  V = mxGetPr(plhs[0]);
  Omega = mxGetPr(plhs[1]);

  // Begin the RPY loop

  for (m=0; m<N; m++){

    int mid = 3*m;

    xm = X[mid];
    ym = X[mid+1];
    zm = X[mid+2];

    am = R[m];

    vx = 0;
    vy = 0;
    vz = 0;

    wx = 0;
    wy = 0;
    wz = 0;

    for (n=0; n<N; n++){

      int nid = 3*n;

      xn = X[nid];
      yn = X[nid+1];
      zn = X[nid+2];

      an = R[n];

      fx = Forces[nid];
      fy = Forces[nid+1];
      fz = Forces[nid+2];

      tx = Torques[nid];
      ty = Torques[nid+1];
      tz = Torques[nid+2];

      r = pow((xm-xn)*(xm-xn) + (ym-yn)*(ym-yn) + (zm-zn)*(zm-zn),0.5) + pow(10,-14);

      rhatx = (xm-xn)/r;
      rhaty = (ym-yn)/r;
      rhatz = (zm-zn)/r;

      if (am>an){amax = am; amin = an;}else{amax = an; amin = am;}

      if (r > am+an){

        a0 = 1/(8*PI*mu[0]*r);
        a1 = 1 + (am*am + an*an)/(3*r*r);
        a2 = 1 - (am*am + an*an)/(r*r);

        dot = fx*rhatx + fy*rhaty + fz*rhatz;

        vx += a0*(a1*fx + a2*dot*rhatx);
        vy += a0*(a1*fy + a2*dot*rhaty);
        vz += a0*(a1*fz + a2*dot*rhatz);

        a0 = 1/(16*PI*mu[0]*r*r*r);

        dot = tx*rhatx + ty*rhaty + tz*rhatz;

        wx += a0*(3*dot*rhatx - tx);
        wy += a0*(3*dot*rhaty - ty);
        wz += a0*(3*dot*rhatz - tz);

        a0 = 1/(8*PI*mu[0]*r*r);

        vx += a0*(ty*rhatz - tz*rhaty);
        vy += a0*(tz*rhatx - tx*rhatz);
        vz += a0*(tx*rhaty - ty*rhatx);

        wx += a0*(fy*rhatz - fz*rhaty);
        wy += a0*(fz*rhatx - fx*rhatz);
        wz += a0*(fx*rhaty - fy*rhatx);

      }else if (r > (amax - amin)){

        a0 = 1/(6*PI*mu[0]*am*an);
        a1 = (16*(am+an)*r*r*r - pow((am-an)*(am-an) + 3*r*r,2))/(32*r*r*r);
        a2 = (3*pow((am-an)*(am-an) - r*r,2))/(32*r*r*r);

        dot = fx*rhatx + fy*rhaty + fz*rhatz;

        vx += a0*(a1*fx + a2*dot*rhatx);
        vy += a0*(a1*fy + a2*dot*rhaty);
        vz += a0*(a1*fz + a2*dot*rhatz);

        a0 = 1/(8*PI*mu[0]*pow(am*an,3));
        A = (5*pow(r,6) - 27*pow(r,4)*(am*am + an*an) + 32*pow(r,3)*(am*am*am + an*an*an) - 9*r*r*pow(am*am - an*an,2) - pow(am-an,4)*(am*am + 4*am*an + an*an))/(64*r*r*r);
        B = 3*pow((am-an)*(am-an) - r*r,2)*(am*am + 4*am*an + an*an - r*r)/(64*r*r*r);

        dot = tx*rhatx + ty*rhaty + tz*rhatz;

        wx += a0*(A*tx + B*dot*rhatx);
        wy += a0*(A*ty + B*dot*rhaty);
        wz += a0*(A*tz + B*dot*rhatz);

        a0 = (an*an + 2*an*(am+r) - 3*pow(am-r,2))*pow(am-an+r,2)/(8*16*r*r*PI*mu[0]*an*am*am*am);

        vx += a0*(ty*rhatz - tz*rhaty);
        vy += a0*(tz*rhatx - tx*rhatz);
        vz += a0*(tx*rhaty - ty*rhatx);

        wx += a0*(fy*rhatz - fz*rhaty);
        wy += a0*(fz*rhatx - fx*rhatz);
        wz += a0*(fx*rhaty - fy*rhatx);

      }else{

        vx += fx/(6*PI*mu[0]*amax);
        vy += fy/(6*PI*mu[0]*amax);
        vz += fz/(6*PI*mu[0]*amax);

        wx += tx/(8*PI*mu[0]*pow(amax,3));
        wy += ty/(8*PI*mu[0]*pow(amax,3));
        wz += tz/(8*PI*mu[0]*pow(amax,3));

        if (am-an<0){heaviside_val = 0;}else if (am-an == 0){heaviside_val = 0.5;}else{heaviside_val = 1;}

        a0 = heaviside_val*r/(8*PI*mu[0]*am*am*am);

        vx += a0*(ty*rhatz - tz*rhaty);
        vy += a0*(tz*rhatx - tx*rhatz);
        vz += a0*(tx*rhaty - ty*rhatx);

        wx += a0*(fy*rhatz - fz*rhaty);
        wy += a0*(fz*rhatx - fx*rhatz);
        wz += a0*(fx*rhaty - fy*rhatx);

      }

    }

    V[mid] = vx;
    V[mid+1] = vy;
    V[mid+2] = vz;

    Omega[mid] = wx;
    Omega[mid+1] = wy;
    Omega[mid+2] = wz;

  }

}
