#ifndef _SPLINECOEFFFORFITTINGLJPOTENERGY_H_
#define _SPLINECOEFFFORFITTINGLJPOTENERGY_H_

#include "ForceField.h"
#include "rg_RelativeOp.h"

// test
#include <iostream>

const rg_INT  DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY          = 2;
const rg_INT  NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY   = 100;
const rg_REAL MIN_DISTANCE_BETWEEN_ATOM_CENTERS               = 0.0;  // unit: Angstrom
const rg_REAL MAX_DISTANCE_BETWEEN_ATOM_CENTERS               = 100.0; // unit: Angstrom
const rg_REAL DIFF_MAX_MIN_DISTANCE_BETWEEN_ATOM_CENTERS      = MAX_DISTANCE_BETWEEN_ATOM_CENTERS - MIN_DISTANCE_BETWEEN_ATOM_CENTERS;
const rg_REAL NUM_INTERVALS_SPLINE_OVER_DIFF_MAX_MIN_DISTANCE = NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY / DIFF_MAX_MIN_DISTANCE_BETWEEN_ATOM_CENTERS;

//const rg_INT  DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY        = 2;
//const rg_INT  NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY = 90;
//const rg_REAL MIN_DISTANCE_BETWEEN_ATOM_CENTERS             = 1.0;  // unit: Angstrom
//const rg_REAL MAX_DISTANCE_BETWEEN_ATOM_CENTERS             = 10.0; // unit: Angstrom

struct SplineCoeffForFittingLJPotEnergy
{
	rg_REAL m_coeff[NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY + 1];

	SplineCoeffForFittingLJPotEnergy()
	{
		rg_INDEX i;
		for (i = 0;i < NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY;i++)
		{
			m_coeff[ i ][ 0 ] = m_coeff[ i ][ 1 ] = m_coeff[ i ][ 2 ] = 0.0;
		}
	}

	inline rg_REAL evaluatePtOnSpline(const rg_REAL& distBTWCenters)
	{
		//if(rg_GE(distBTWCenters, MAX_DISTANCE_BETWEEN_ATOM_CENTERS))
		if(distBTWCenters >= MAX_DISTANCE_BETWEEN_ATOM_CENTERS - rg_MATH_RES)
		{
			// test
			//cout << NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY - 1 << endl;
			
			rg_INDEX indexOfInterval = NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY - 1;
			// evaluate by Horner's method			
			//rg_REAL ptOnSpline = m_coeff[indexOfInterval][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY];
			//rg_INDEX i;
			//for (i = DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY - 1;i >=0;i--)
			//{
			//	ptOnSpline = ptOnSpline * MAX_DISTANCE_BETWEEN_ATOM_CENTERS + m_coeff[indexOfInterval][ i ];
			//}
			
			//return ptOnSpline;

			return (m_coeff[indexOfInterval][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY] * MAX_DISTANCE_BETWEEN_ATOM_CENTERS + m_coeff[indexOfInterval][ 1 ]) * MAX_DISTANCE_BETWEEN_ATOM_CENTERS + m_coeff[indexOfInterval][ 0 ];
			// return zero ??? DKIM
		}
		//else if(rg_LT(distBTWCenters, MIN_DISTANCE_BETWEEN_ATOM_CENTERS))

		// remove this block....DKIM
		//else if(distBTWCenters < MIN_DISTANCE_BETWEEN_ATOM_CENTERS - rg_MATH_RES)
		//{
			// test
			//cout << 0 << endl;

			//rg_INDEX indexOfInterval = 0;
			// evaluate by Horner's method
			//rg_REAL ptOnSpline = m_coeff[indexOfInterval][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY];
			//rg_INDEX i;
			//for (i = DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY - 1;i >=0;i--)
			//{
			//	ptOnSpline = ptOnSpline * MIN_DISTANCE_BETWEEN_ATOM_CENTERS + m_coeff[indexOfInterval][ i ];
			//}
			
			//return ptOnSpline;

			//return (m_coeff[ 0 ][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY] * MIN_DISTANCE_BETWEEN_ATOM_CENTERS + m_coeff[ 0 ][ 1 ]) * MIN_DISTANCE_BETWEEN_ATOM_CENTERS + m_coeff[ 0 ][ 0 ];
		//}
		// MIN_DISTANCE_BETWEEN_ATOM_CENTERS <= distBTWCenters < MAX_DISTANCE_BETWEEN_ATOM_CENTERS
		else
		{
			//rg_INDEX indexOfInterval = getIndexOfIntervalForSplineCorrTo(distBTWCenters);
			// test
			//cout << indexOfInterval << endl;

			// evaluate by Horner's method
			//rg_REAL ptOnSpline = m_coeff[indexOfInterval][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY];
			//rg_INDEX i;
			//for (i = DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY - 1;i >=0;i--)
			//{
			//	ptOnSpline = ptOnSpline * distBTWCenters + m_coeff[indexOfInterval][ i ];
			//}

			//return ptOnSpline;

			rg_INDEX indexOfInterval = (rg_INDEX) distBTWCenters;
			return (m_coeff[indexOfInterval][DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY] * distBTWCenters + m_coeff[indexOfInterval][ 1 ]) * distBTWCenters + m_coeff[indexOfInterval][ 0 ];			
		}		
	}

	inline rg_INDEX getIndexOfIntervalForSplineCorrTo(const rg_REAL& distBTWCenters)
	{
		//rg_INDEX indexOfInterval = (rg_INDEX) ( ( NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY ) * (distBTWCenters - MIN_DISTANCE_BETWEEN_ATOM_CENTERS) / (MAX_DISTANCE_BETWEEN_ATOM_CENTERS - MIN_DISTANCE_BETWEEN_ATOM_CENTERS) );
		//return indexOfInterval;

		rg_INDEX indexOfInterval = (rg_INDEX) ( distBTWCenters );
		// In principle, we have to use the statement "(rg_INDEX) ( NUM_INTERVALS_SPLINE_OVER_DIFF_MAX_MIN_DISTANCE * ( distBTWCenters - MIN_DISTANCE_BETWEEN_ATOM_CENTERS ) )"
		// However, we need not do that because of the following two reason:		
		// (1) MIN_DISTANCE_BETWEEN_ATOM_CENTERS = 0.0, MAX_DISTANCE_BETWEEN_ATOM_CENTERS = 100.0
		//     , and NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY = 100
		// (2) Hence, NUM_INTERVALS_SPLINE_OVER_DIFF_MAX_MIN_DISTANCE = 100 / 100.0 = 1
		// NOTE: Currently, we use squared distance instead of distance.
		return indexOfInterval;
	}

	//inline rg_REAL evaluateDerivativeOfLJPotEnergy(const rg_REAL& distBTWCenters)
	//{
	//	rg_REAL derivative = 0.0;

	//	return derivative;
	//}
};

#endif