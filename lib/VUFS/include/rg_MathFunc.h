//********************************************************************
//
//	  FILENAME    : rg_MathFunc.h
//	  
//    DESCRIPTION : 
//           This is the interface of various function for mathmatics
//                          
//
//    AUTHOR      : Deok-Soo Kim, Sun-woong Lee
//    START DATE  : 21 Jun 1996    
//
//    HISTORY     :
//				1. 1998.2.7  Tae-Bum Jang revise the following function 
//                   rg_REAL         combination(const rg_REAL &n, const rg_REAL &i);
//				     rg_REAL         factorial(const rg_INT &n);
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************


#ifndef _MATHFUNC_EXT
#define _MATHFUNC_EXT

#include <math.h>
#include "rg_Const.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"
#include <list>
using namespace std;


class rg_MathFunc
{
public:
    static rg_REAL           cubicRoot(const rg_REAL &r);
    static rg_ComplexNumber *cubicRoot(const rg_ComplexNumber &c);
    static rg_ComplexNumber *root(const rg_ComplexNumber &c, const rg_INT &n);
    //rg_Point2D       *solveQuadraticEq(const rg_Point2D &a, const rg_Point2D &b, const rg_Point2D &c);
    static rg_REAL          *solveQuadraticEq(const rg_REAL &a, const rg_REAL &b, const rg_REAL &c);
    static rg_ComplexNumber *solveQuadraticEq(rg_ComplexNumber &a, rg_ComplexNumber &b, rg_ComplexNumber &c,  rg_INT& num);
    //rg_REAL         crossProduct(const rg_Point2D &p1, const rg_Point2D &p2);
    static rg_REAL         minValue(const rg_REAL &a, const rg_REAL &b, const rg_REAL &c);
    static rg_REAL         maxValue(const rg_REAL &a, const rg_REAL &b, const rg_REAL &c);
    static rg_REAL         norm(const rg_REAL *vec, const rg_INT &vectorSize, const rg_INT &k);
    static rg_INT         combination(const rg_INT &n, const rg_INT &i);
    static rg_REAL         factorial(const rg_INT &n);

    static void            compute_mean_and_standard_deviation(const list<rg_REAL>& observations, rg_REAL& mean, rg_REAL& standardDeviation);

    static inline rg_REAL   round(const rg_REAL& x)
    {
        return (rg_REAL)(x > 0 ? rg_INT(x + 0.5) : rg_INT(x - 0.5));
    }

    // Generating all the possible sequences that consists of 0's and 1's
    static rg_INT** enumerateSubset(rg_INT* set, rg_INT sizeOfSet, rg_INT sizeOfSubset);
    static rg_INT** enumerateZeroOneSequence(rg_INT noOfZero, rg_INT noOfOne);
    static rg_INT** enumerateZeroOneSequenceRevised(rg_INT noOfZero, rg_INT noOfOne);

    static rg_REAL randomNumber(const rg_REAL& min, const rg_REAL& max);
    static void    shuffleSequence(rg_INT*  inputSequence,
                                   rg_INT*& shuffledSequence, 
                                   const rg_INT& seqSize);

// Numerical computations
static rg_REAL compute_root_of_polynomial_via_Newton_method_with_initial_solution(rg_Polynomial& polynomial, const rg_REAL& initialSolution, const rg_REAL& terminatingTolerance = rg_MATH_RES, const rg_INT& numIterations = 40);
};
#endif

