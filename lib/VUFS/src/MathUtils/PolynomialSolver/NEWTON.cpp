/******************************************************************/
/*                                                                */
/* file:          newton.c                                        */
/*                                                                */
/* main function: newton()                                        */
/*                                                                */
/* version:       1.2                                             */
/*                                                                */
/* author:        B. Frenzel                                      */
/*                                                                */
/* date:          Jan 7 1993                                      */ 
/*                                                                */
/* input:         p[]     coefficient vector of the original      */
/*                        polynomial                              */
/*                n       the highest exponent of the original    */
/*                        polynomial                              */  
/*                ns      determined root with Muller method      */
/*                                                                */
/* output:        *dxabs error of determined root                 */
/*                                                                */
/* return:        xmin    determined root                         */ 
/* subroutines:   fdvalue()                                       */
/*                                                                */
/* description:                                                   */
/* newton() determines the root of a polynomial with complex      */
/* coefficients; the initial estimation is the root determined    */
/* by Muller's method                                             */
/*                                                                */
/* Copyright:                                                     */
/* Lehrstuhl fuer Nachrichtentechnik Erlangen                     */
/* Cauerstr. 7, 8520 Erlangen, FRG, 1993                          */
/* e-mail: rg_INT@nt.e-technik.uni-erlangen.de                       */
/*                                                                */
/******************************************************************/

#define  NEWTON
#include "PolynomialSolver.h"









/***** main routine of the Newton method *****/
dcomplex newton(dcomplex *p,rg_INT n,dcomplex ns,rg_REAL *dxabs,
                unsigned char flag)
/*dcomplex *p,         coefficient vector of the original polynomial   */
/*         ns;         determined root with Muller method              */
/*rg_INT      n;          highest exponent of the original polynomial     */
/*rg_REAL   *dxabs;     dxabs  = |P(x0)/P'(x0)|                         */
/*unsigned char flag;  flag = rg_TRUE  => complex coefficients            */
/*                     flag = rg_FALSE => real    coefficients            */
{
     dcomplex x0,                /* iteration variable for x-value     */
              xmin,              /* best x determined in newton()      */
              f,                 /* f       = P(x0)                    */
              df,                /* df      = P'(x0)                   */
              dx,                /* dx      = P(x0)/P'(x0)             */
              dxh;               /* help variable dxh = P(x0)/P'(x0)   */
     rg_REAL   fabsmin=FVALUE,    /* fabsmin = |P(xmin)|                */
              eps=DBL_EPSILON;   /* routine ends, when estimated dist. */
                                 /* between x0 and root is less or     */
                                 /* equal eps                          */
     rg_INT      iter   =0,         /* counter                            */
              noise  =0;         /* noisecounter                       */

     x0   = ns;             /* initial estimation = root determined    */
                            /* with Muller method                      */
     xmin = x0;             /* initial estimation for the best x-value */
     dx   = Complex(1.,0.); /* initial value: P(x0)/P'(x0)=1+j*0       */
     *dxabs = Cabs(dx);     /* initial value: |P(x0)/P'(x0)|=1         */
     /* printf("%8.4e %8.4e\n",xmin.r,xmin.i); */
     for (iter=0;iter<ITERMAX;iter++) { /* main loop                   */
          fdvalue(p,n,&f,&df,x0,rg_TRUE);  /* f=P(x0), df=P'(x0)          */

          if (Cabs(f)<fabsmin) {  /* the new x0 is a better            */
               xmin = x0;         /* approximation than the old xmin   */
               fabsmin = Cabs(f); /* store new xmin and fabsmin        */
               noise = 0;         /* reset noise counter               */
          }

          if (Cabs(df)!=0.) {     /* calculate new dx                  */
               dxh=Cdiv(f,df);
               if (Cabs(dxh)<*dxabs*FACTOR) { /* new dx small enough?  */
                    dx = dxh;          /* store new dx for next        */
                    *dxabs = Cabs(dx); /* iteration                    */
               }
          }

          if (Cabs(xmin)!=0.) {
               if (*dxabs/Cabs(xmin)<eps || noise == NOISEMAX) { 
                                                  /* routine ends      */ 
                    if (fabs(xmin.i)<BOUND && flag==0) {
                                     /* define determined root as real,*/ 
                         xmin.i=0.;  /* if imag. part<BOUND            */
		       }
                    *dxabs=*dxabs/Cabs(xmin); /* return relative error */ 
                    return xmin;     /* return best approximation      */
               }
          } 
          
          x0 = Csub(x0,dx); /* main iteration: x0 = x0 - P(x0)/P'(x0)  */

          noise++;  /* increase noise counter                          */ 
     }


     if (fabs(xmin.i)<BOUND && flag==0) /* define determined root      */
          xmin.i=0.;         /* as real, if imag. part<BOUND           */
     if (Cabs(xmin)!=0.)
          *dxabs=*dxabs/Cabs(xmin); /* return relative error           */ 
                             /* maximum number of iterations exceeded: */
     return xmin;            /* return best xmin until now             */
}


