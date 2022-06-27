#ifndef _XVOLNPOTENERGYVALUEPAIRS_H_
#define _XVOLNPOTENERGYVALUEPAIRS_H_

//#include "rg_Const.h"
#include "ForceField.h"
#include "rg_RelativeOp.h"

// test
#include <iostream>
#include <cmath>

const rg_INT NUM_XVOLLJPOTENERGYVAL_PAIRS = 100;

#include "ConstForXVolLJPotEnergy.h"

struct XVolNLJPotEnergyValuePairs
{
	rg_REAL m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS][ 2 ];
	// For i = 1, 2, ..., NUM_XVOLLJPOTENERGYVAL_PAIRS
	// m_XVolEnergyPairs[ i-1 ][ 0 ]: XVolume   w.r.t. distance = m_sumOfTwoRadii * i / NUM_XVOLLJPOTENERGYVAL_PAIRS
	// m_XVolEnergyPairs[ i-1 ][ 1 ]: LJ_Energy w.r.t. distance = m_sumOfTwoRadii * i / NUM_XVOLLJPOTENERGYVAL_PAIRS

	// For i = 0, 1, ..., NUM_XVOLLJPOTENERGYVAL_PAIRS - 1
	// m_XVolEnergyPairs[ i ][ 0 ]: XVolume   w.r.t. distance = m_sumOfTwoRadii * (1 - i / NUM_XVOLLJPOTENERGYVAL_PAIRS)
	// m_XVolEnergyPairs[ i ][ 1 ]: LJ_Energy w.r.t. distance = m_sumOfTwoRadii * (1 - i / NUM_XVOLLJPOTENERGYVAL_PAIRS)
	
	rg_REAL m_sumOfTwoRadii;                                    // R = r1 + r2, summation of radii for two atoms (r1 + r2)
	rg_REAL m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS;  // R / NUM_XVOLLJPOTENERGYVAL_PAIRS
	rg_REAL m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_diffMinMaxDist; // NUM_XVOLLJPOTENERGYVAL_PAIRS / (R - R / NUM_XVOLLJPOTENERGYVAL_PAIRS)	
	rg_REAL m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_sumOfTwoRadii;  // NUM_XVOLLJPOTENERGYVAL_PAIRS /  R

	XVolNLJPotEnergyValuePairs()
	{
		rg_INT i;
		for (i = 0;i < NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
		{
			m_XVolEnergyPairs[ i ][ 0 ] = m_XVolEnergyPairs[ i ][ 1 ] = 0.0;
		}
		m_sumOfTwoRadii = 0.0;
		m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS = 0.0;
		m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_diffMinMaxDist = 0.0;
		m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_sumOfTwoRadii = 0.0;
	}

	inline rg_REAL estimateLJPotEnergyValueOfXVolume(const rg_REAL& distBTWCenters, const rg_REAL& XVolume)
	{

#ifdef _DEBUG
        GaussianScaleFunc(XVolume, 1.0);
#endif  
		//if(rg_GE(distBTWCenters, m_sumOfTwoRadii))
		if(distBTWCenters >= m_sumOfTwoRadii - rg_MATH_RES)
		{
			// test
			//cout << NUM_XVOLLJPOTENERGYVAL_PAIRS - 1 << endl;

			return m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 1 ];
            //return linearScaleFunc(m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 0 ], m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 1 ]);
		}
		//else if(rg_LT(distBTWCenters, m_sumOfTwoRadii / NUM_XVOLLJPOTENERGYVAL_PAIRS))
		else if(distBTWCenters < m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS - rg_MATH_RES)
		{
			// test
			//cout << 0 << endl;

			return m_XVolEnergyPairs[ 0 ][ 1 ];
            //return linearScaleFunc(m_XVolEnergyPairs[ 0 ][ 0 ], m_XVolEnergyPairs[ 0 ][ 1 ]);
		}
		// m_sumOfTwoRadii / NUM_XVOLLJPOTENERGYVAL_PAIRS <= distBTWCenters < m_sumOfTwoRadii
		else
		{
			// construct LJEnergy by linear interpolation
			rg_INDEX indexOfInterval = getIndexOfIntervalContainingXVolEnergyPairCorrTo(distBTWCenters);
			rg_REAL  currXVolume     = m_XVolEnergyPairs[indexOfInterval    ][ 0 ];
			rg_REAL  nextXVolume     = m_XVolEnergyPairs[indexOfInterval + 1][ 0 ];
            rg_REAL  currLJEnergy    = m_XVolEnergyPairs[indexOfInterval    ][ 1 ];
            rg_REAL  nextLJEnergy    = m_XVolEnergyPairs[indexOfInterval + 1][ 1 ];
			//rg_REAL  currLJEnergy    = linearScaleFunc(m_XVolEnergyPairs[indexOfInterval    ][ 0 ], m_XVolEnergyPairs[indexOfInterval    ][ 1 ]);
			//rg_REAL  nextLJEnergy    = linearScaleFunc(m_XVolEnergyPairs[indexOfInterval + 1][ 0 ], m_XVolEnergyPairs[indexOfInterval + 1][ 1 ]);

			rg_REAL  paramOfLinearInterpolation = 0.0;
			if(fabs(currXVolume - nextXVolume) > rg_MATH_RES)
			// We use the intersection volume in order to construct the parameter of linear interpolation
				paramOfLinearInterpolation = (currXVolume - XVolume) / (currXVolume - nextXVolume);			
			else
			// When distBTWCenters < |r1 - r2|, this case occurs
			// In this case, we use the distance between atom centers in order to get parameter of linear interpolation.
			
			// Let d= distBTWCenters, R = sumOfTwoRadii, N = NUM_XVOLLJPOTENERGYVAL_PAIRS
			// Because "indexOfInterval" is the index corresponding to "distBTWCenters",
			// "(indexOfInterval + 1) * R / N" is the distance corresponding to "currXVolume" and
			// "(indexOfInterval + 2) * R / N" is the distance corresponding to "nextXVolume" 
				paramOfLinearInterpolation = (distBTWCenters - m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS * (indexOfInterval + 1) ) / m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS;
			// Hence,
			// paramOfLinearInterpolation = (d - R / N * (indexOfInterval + 1)) / ((indexOfInterval + 2) * R / N - (indexOfInterval + 1) * R / N)
			//                            = (d - R / N * (indexOfInterval + 1)) / (R / N)


			rg_REAL  LJEnergy = 0.0;
			LJEnergy = (1 - paramOfLinearInterpolation) * currLJEnergy + paramOfLinearInterpolation * nextLJEnergy;

			// test			
			//if(rg_NEG(paramOfLinearInterpolation) || rg_GT(paramOfLinearInterpolation, 1.0))
			//	cout << index << paramOfLinearInterpolation << " wired ";
			//else
			//	cout << index << " ";

			// test
			//cout << indexOfInterval << endl;

			return LJEnergy;
		}
	}

	inline rg_INDEX getIndexOfIntervalContainingXVolEnergyPairCorrTo(const rg_REAL& distBTWCenters)
	{
		// (1) The following formula is derived from: 
		// index = (NUM_XVOLLJPOTENERGYVAL_PAIRS - 1) * (m_sumOfTwoRadii - distBTWCenters) / (m_sumOfTwoRadii - m_sumOfTwoRadii / NUM_XVOLLJPOTENERGYVAL_PAIRS)
		// or
		// index = (NUM_XVOLLJPOTENERGYVAL_PAIRS - 1) - (NUM_XVOLLJPOTENERGYVAL_PAIRS - 1) * (distBTWCenters - m_sumOfTwoRadii / NUM_XVOLLJPOTENERGYVAL_PAIRS) / (m_sumOfTwoRadii - m_sumOfTwoRadii / NUM_XVOLLJPOTENERGYVAL_PAIRS)
		//rg_INDEX index = (rg_INDEX) ( NUM_XVOLLJPOTENERGYVAL_PAIRS * (m_sumOfTwoRadii - distBTWCenters) / m_sumOfTwoRadii );		

		//rg_INDEX indexOfInterval = (rg_INDEX) ( NUM_XVOLLJPOTENERGYVAL_PAIRS * (NUM_XVOLLJPOTENERGYVAL_PAIRS * distBTWCenters - m_sumOfTwoRadii) / (m_sumOfTwoRadii * (NUM_XVOLLJPOTENERGYVAL_PAIRS - 1)) );

		//rg_INDEX indexOfInterval = (rg_INDEX) (m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_diffMinMaxDist * (distBTWCenters - m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS));

		rg_INDEX indexOfInterval = (rg_INDEX) (m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_sumOfTwoRadii * distBTWCenters - 1);
		// The above formula is derived as follows:
		// Let d= distBTWCenters, R = sumOfTwoRadii, N = NUM_XVOLLJPOTENERGYVAL_PAIRS, and
		// (N - 1): # of intervals.
		// Then, 
		// indexOfInterval = (d - R/N) / (R - R/N) * (N - 1)
		//                 = (N*d - R) / (N*R - R) * (N - 1)
		//                 = (N/R * d - 1).
		return indexOfInterval;
	}

    // test scaling functions
    inline rg_REAL constantScaleFunc(const rg_REAL& energy)
    {
        return 10.0 * energy;
    }

    inline rg_REAL linearScaleFunc(const rg_REAL& xvol, const rg_REAL& energy)
    {
        return 10000000000.0 * xvol * energy;
    }

    inline rg_REAL GaussianScaleFunc(const rg_REAL& xvol, const rg_REAL& energy)
    {
        rg_REAL mu    = (4*m_XVolEnergyPairs[ 0 ][ 0 ] + m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 0 ]) / 5.0;
        rg_REAL sigma = (m_XVolEnergyPairs[ 0 ][ 0 ] - m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 0 ]) / 10.0; // 10sigma
        rg_REAL delta = pow(10.0, 8.0);
        //rg_REAL maxFuncValue = 10000000.0 * m_XVolEnergyPairs[ 0 ][ 0 ];
        //rg_REAL scale = maxFuncValue * sigma * sqrt( 2.0 * rg_PI );
        rg_REAL scale = sigma * sqrt( 2.0 * rg_PI ) * (delta*10.0*10.0 - delta);

        //return exp( -(xvol - mu)*(xvol - mu)/(2.0 * sigma*sigma) ) / ( sigma * sqrt( 2.0 * rg_PI ) ) * energy;
        //return exp( -(xvol - mu)*(xvol - mu)/(2.0 * sigma*sigma) ) * maxFuncValue * energy;
        return ( exp( -(xvol - mu)*(xvol - mu)/(2.0 * sigma*sigma) ) * (delta*10.0*10.0 - delta) + delta ) * energy;        
    }	
};
//struct XVolNLJPotEnergyValuePairs
//{
//	rg_REAL m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS][ 2 ];
//	rg_REAL m_lengthOfXVolumeRange;
//
//	XVolNLJPotEnergyValuePairs()
//	{
//		rg_INT i;
//		for (i = 0;i < NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
//		{
//			m_XVolEnergyPairs[ i ][ 0 ] = m_XVolEnergyPairs[ i ][ 1 ] = 0.0;
//		}
//		m_lengthOfXVolumeRange = 0.0;
//	}
//
//	inline rg_REAL estimateLJPotEnergyValueOfXVolume(const rg_REAL& XVolume)
//	{
//		if(rg_GE(XVolume, m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 0 ]))
//			return m_XVolEnergyPairs[NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 1 ];
//
//		rg_INDEX index = getIndexOfIntervalContainingXVol(XVolume);
//		rg_REAL  currXVolume    = m_XVolEnergyPairs[index    ][ 0 ];
//		rg_REAL  nextXVolume    = m_XVolEnergyPairs[index + 1][ 0 ];
//		rg_REAL  currLJEnergy   = m_XVolEnergyPairs[index    ][ 1 ];
//		rg_REAL  nextLJEnergy   = m_XVolEnergyPairs[index + 1][ 1 ];
//
//		rg_REAL  paramOfLinearInterpolation = (XVolume - currXVolume) / (nextXVolume - currXVolume);
//		rg_REAL  LJEnergy = 0.0;
//		LJEnergy = (1 - paramOfLinearInterpolation) * currLJEnergy + paramOfLinearInterpolation * nextLJEnergy;
//
//		return LJEnergy;
//	}
//
//	inline rg_INDEX getIndexOfIntervalContainingXVol(const rg_REAL& XVolume)
//	{
//		rg_INDEX index = (rg_INDEX)((NUM_XVOLLJPOTENERGYVAL_PAIRS - 1) * (XVolume - m_XVolEnergyPairs[ 0 ][ 0 ]) / m_lengthOfXVolumeRange);
//		return index;
//	}
//};
#endif
