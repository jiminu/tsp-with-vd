#ifndef _FUNCTIONSFORSPLINEFITTINGLJPOTENERGY_H_
#define _FUNCTIONSFORSPLINEFITTINGLJPOTENERGY_H_

#include "SplineCoeffForFittingLJPotEnergy.h"
#include "ForceField.h"
#include "ConstForXVolLJPotEnergy.h"

//#include <time.h>

class FunctionsForSplineFittingLJPotEnergy
{
public:
	static SplineCoeffForFittingLJPotEnergy m_splineCoefTablefForFittingLJPotEnergy[NUM_ATOM_TYPE_IN_AMBER][NUM_ATOM_TYPE_IN_AMBER];
	static SplineCoeffForFittingLJPotEnergy m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[NUM_ATOM_TYPE_IN_AMBER][NUM_ATOM_TYPE_IN_AMBER];

	// test
	//static rg_REAL m_time4EvalSplineOfLJEnergy;
	//static rg_REAL m_time4EvalSplineOfLJDerivativeOfEnergy;

private:
    static rg_BOOL                    m_bSplineFittingCoeffLoaded;
	static void readSplineCoeffForAtomPair(const string& fileNameWithPathForAtomPair, SplineCoeffForFittingLJPotEnergy& splineCoeff);
	//static void readSplineCoeffForFittingDerivativeOfLJPotEnergyForAtomPair(const string& fileNameWithPathForAtomPair, SplineCoeffForFittingLJPotEnergy& splineCoeff);
	//static string getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(const AmberAtomTypes& atomType1, const AmberAtomTypes& atomType2);
    
public:
	static inline rg_REAL evaluateLJPotEnergyBySplineFitting(const AmberAtomTypes& atomType1, 
		                                              const AmberAtomTypes& atomType2, 
		                                              const rg_REAL& distBTWCenters   )
	{
		//rg_INDEX index1 = (rg_INDEX) atomType1;
		//rg_INDEX index2 = (rg_INDEX) atomType2;

		//rg_REAL energy = 0.0;
		//energy = m_splineCoefTablefForFittingLJPotEnergy[index1][index2].evaluatePtOnSpline(distBTWCenters);

		//return energy;

		return m_splineCoefTablefForFittingLJPotEnergy[atomType1][atomType2].evaluatePtOnSpline(distBTWCenters);
	}

	static inline rg_REAL evaluateDerivativeOfLJPotEnergyBySplineFitting(const AmberAtomTypes& atomType1, 
		                                                          const AmberAtomTypes& atomType2, 
		                                                          const rg_REAL& distBTWCenters   )
	{
		//rg_INDEX index1 = (rg_INDEX) atomType1;
		//rg_INDEX index2 = (rg_INDEX) atomType2;

		//rg_REAL derivative = 0.0;
		//derivative = m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[index1][index2].evaluatePtOnSpline(distBTWCenters);

		//return derivative;

		return m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[atomType1][atomType2].evaluatePtOnSpline(distBTWCenters);
	}

	static void writeSplineCoeffTables(const string& fileNameWithPathOfSplineCoeffTableForFittingLJPotEnergy, const string& fileNameWithPathOfSplineCoeffTableForFittingDerivativeOfLJPotEnergy);
		static void writeSplineCoeffTableForFittingLJPotEnergyIntoBinaryFile(const string& fileNameWithPath);
		static void writeSplineCoeffTableForFittingDerivativeOfLJPotEnergyIntoBinaryFile(const string& fileNameWithPath);
	static void loadSplineCoeffTables(const string& fileNameWithPathOfSplineCoeffTableForFittingLJPotEnergy, const string& fileNameWithPathOfSplineCoeffTableForFittingDerivativeOfLJPotEnergy);
		static void loadSplineCoeffTableForFittingLJPotEnergy(const string& fileNameWithPath);
		static void loadSplineCoeffTableForFittingDerivativeOfLJPotEnergy(const string& fileNameWithPath);

            static string getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(const AmberAtomTypes& atomType1, const AmberAtomTypes& atomType2);

};

#endif