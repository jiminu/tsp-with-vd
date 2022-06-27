/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : EigenDecomposition.h
//	  
//    DESCRIPTION : collections of function related with Eigen Decomposition
//                          
//	  CLASS NAME  : 
//
//    BASE CLASS  : 
//      
//    AUTHOR      : Ryu, JungHyun and William H. Press et al.
//
//    START DATE  :     
//    
//    HISTORY :
//
//					 
//
//           Copyright ¨Ï 1998 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _EIGENDECOMPOSITION_IS_INCLUDED
#define _EIGENDECOMPOSITION_IS_INCLUDED
#include "rg_Const.h"
#include <math.h>
#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

void elmhes(rg_REAL **a, rg_INT n);
//Reduction to Hessenberg form by the elimination method. The real, nonsymmetric matrix
//a[1..n][1..n] is replaced by an upper Hessenberg matrix with identical eigenvalues. Rec-
//ommended, but not required, is that this routine be preceded by balanc. On output, the
//Hessenberg matrix is in elements a[i][j] with i . j+1. Elements with i > j+1 are to be
//thought of as zero, but are returned with random values.


void hqr(rg_REAL **a, rg_INT n, rg_REAL wr[], rg_REAL wi[]);
//Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On input a can be
//exactly as output from elmhes x 11.5; on output it is destroyed. The real and imaginary parts
//of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.

#endif

