#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define PI 3.14159265358979323846264338327950288419716939937510

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    // Input declarations
    
    int N;
    double *mu;
    double *X, *Forces, *Torques, *a;
    
    // Output declarations
    
    double *V, *Omega;
    
    // Internal variable declarations
    
    int m, n;
    double r, a0, a1, a2, a3, a4, h, dot, fx, fy, rhatx, rhaty, rhatz,
            fz, tx, ty, tz, ex, ey, ez, xm, ym, zm, xn, yn, zn, vx, vy, vz, wx, wy, wz;
    
    // Assign input and output
    
    X = mxGetPr(prhs[0]);
    Forces = mxGetPr(prhs[1]);
    Torques = mxGetPr(prhs[2]);
    a = mxGetPr(prhs[3]);
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
            
            fx = Forces[nid];
            fy = Forces[nid+1];
            fz = Forces[nid+2];
            
            tx = Torques[nid];
            ty = Torques[nid+1];
            tz = Torques[nid+2];
            
            if (m==n){
                
                // We include the effect of the wall in the self mobility
                
                h = zm/a[0];
                
                a0 = (1 - (9/h - 2*pow(h,-3) + pow(h,-5))/16)/(6*PI*mu[0]*a[0]);
                a1 = (1 - (9/h - 4*pow(h,-3) + pow(h,-5))/8)/(6*PI*mu[0]*a[0]);
                a2 = 3*pow(h,-4)/(32*6*PI*mu[0]*a[0]);
                
                vx += a0*fx - a2*ty;
                vy += a0*fy + a2*tx;
                vz += a1*fz;
                
                a0 = 1/(8*PI*mu[0]*pow(a[0],3)) - 15*pow(h,-3)/(64*6*PI*mu[0]*a[0]);
                a1 = 1/(8*PI*mu[0]*pow(a[0],3)) - 3*pow(h,-3)/(32*6*PI*mu[0]*a[0]);
                a2 = 3*pow(h,-4)/(32*6*PI*mu[0]*a[0]);
                
                wx += a0*tx + a2*fy;
                wy += a0*ty - a2*fx;
                wz += a1*tz;
                
            }else{
                
                // Start with the usual terms
                
                r = pow((xm-xn)*(xm-xn) + (ym-yn)*(ym-yn) + (zm-zn)*(zm-zn),0.5);// + pow(10,-14);
            
                rhatx = (xm-xn)/r;
                rhaty = (ym-yn)/r;
                rhatz = (zm-zn)/r;
        
                a0 = 1/(8*PI*mu[0]*r);
                a1 = 1 + (2*a[0]*a[0])/(3*r*r);
                a2 = 1 - (2*a[0]*a[0])/(r*r);
                
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
                
                // Now we add the wall reflection terms
                
                ex = (xm - xn)/a[0];
                ey = (ym - yn)/a[0];
                ez = (zm + zn)/a[0];
                
                r = pow(ex*ex + ey*ey + ez*ez,0.5);// + pow(10,-14);
                
                ex /= r;
                ey /= r;
                ez /= r;
                
    //            h = zn/(r*a[0]);
                
                h = zn/(zm+zn);
                
                a0 = -0.25*(3*(1 + 2*h*(1-h)*ez*ez)/r + 2*(1 - 3*ez*ez)*pow(r,-3) - 2*(1 - 5*ez*ez)*pow(r,-5));
                a1 = -0.25*(3*(1 - 6*h*(1-h)*ez*ez)/r - 6*(1 - 5*ez*ez)*pow(r,-3) + 10*(1 - 7*ez*ez)*pow(r,-5))*(ex*fx + ey*fy + ez*fz);
                a2 = 0.5*ez*(3*h*(1 - 6*(1-h)*ez*ez)/r - 6*(1 - 5*ez*ez)*pow(r,-3) + 10*(2 - 7*ez*ez)*pow(r,-5))*fz;
                a3 = 0.5*ez*(3*h/r - 10*pow(r,-5))*(ex*fx + ey*fy + ez*fz);
                a4 = -((3*h*h/r + 3*pow(r,-3))*ez*ez + (2 - 15*ez*ez)*pow(r,-5))*fz;
                
                vx += (a0*fx + (a1+a2)*ex)/(6*PI*mu[0]*a[0]);
                vy += (a0*fy + (a1+a2)*ey)/(6*PI*mu[0]*a[0]);
                vz += (a0*fz + (a1+a2)*ez + a3*ez + a4)/(6*PI*mu[0]*a[0]);
                
                a0 = 3*(1 - 6*ez*ez)*pow(r,-3)/8;
                a1 = -9*pow(r,-3)*(ex*tx + ey*ty + ez*tz)/8;
                a2 = 9*ez*pow(r,-3)*(ex*tx + ey*ty + ez*tz)/4;
                a3 = 9*pow(r,-3)*(ex*ty - ey*tx)/4;
                
                wx += (a0*tx + a1*ex - a3*ey)/(6*PI*mu[0]*a[0]);
                wy += (a0*ty + a1*ey + a3*ex)/(6*PI*mu[0]*a[0]);
                wz += (a0*tz + a1*ez + a2)/(6*PI*mu[0]*a[0]);
                
                a0 = 0.75*pow(r,-2);
                a1 = 1.5*fz*(6*h*ez*ez*pow(r,-2) + (1 - 10*ez*ez)*pow(r,-4));
                a2 = -1.5*ez*(ex*fx + ey*fy + ez*fz)*(3*h*pow(r,-2) - 5*pow(r,-4));
                a3 = -1.5*ez*(h*pow(r,-2) - pow(r,-4));
                
                wx += (a0*(fy*ez - fz*ey) - a1*ey - a2*ey + a3*fy)/(6*PI*mu[0]*a[0]);
                wy += (a0*(fz*ex - fx*ez) + a1*ex + a2*ex - a3*fx)/(6*PI*mu[0]*a[0]);
                wz += (a0*(fx*ey - fy*ex) )/(6*PI*mu[0]*a[0]);
                
                h = 1-h;
                ex *= -1.0;
                ey *= -1.0;
                
                a0 = -0.75*pow(r,-2);
                a1 = 1.5*(6*h*ez*ez*pow(r,-2) + (1 - 10*ez*ez)*pow(r,-4))*(ex*ty - ey*tx);
                a2 = -1.5*ez*(3*h*pow(r,-2) - 5*pow(r,-4))*(ex*ty - ey*tx);
                a3 = 1.5*ez*(h*pow(r,-2) - pow(r,-4));
                
                vx += (a0*(ty*ez - tz*ey) + a2*ex + a3*ty)/(6*PI*mu[0]*a[0]);
                vy += (a0*(tz*ex - tx*ez) + a2*ey - a3*tx)/(6*PI*mu[0]*a[0]);
                vz += (a0*(tx*ey - ty*ex) + a1 + a2*ez)/(6*PI*mu[0]*a[0]);

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
