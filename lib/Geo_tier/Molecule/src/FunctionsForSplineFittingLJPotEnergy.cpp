#include "FunctionsForSplineFittingLJPotEnergy.h"
#include <iostream>
#include <fstream>

SplineCoeffForFittingLJPotEnergy FunctionsForSplineFittingLJPotEnergy::m_splineCoefTablefForFittingLJPotEnergy[NUM_ATOM_TYPE_IN_AMBER][NUM_ATOM_TYPE_IN_AMBER];
SplineCoeffForFittingLJPotEnergy FunctionsForSplineFittingLJPotEnergy::m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[NUM_ATOM_TYPE_IN_AMBER][NUM_ATOM_TYPE_IN_AMBER];
rg_BOOL                          FunctionsForSplineFittingLJPotEnergy::m_bSplineFittingCoeffLoaded = rg_FALSE;

// test
//rg_REAL                          FunctionsForSplineFittingLJPotEnergy::m_time4EvalSplineOfLJEnergy = 0.0;
//rg_REAL                          FunctionsForSplineFittingLJPotEnergy::m_time4EvalSplineOfLJDerivativeOfEnergy = 0.0;

void FunctionsForSplineFittingLJPotEnergy::readSplineCoeffForAtomPair(const string& fileNameWithPathForAtomPair, SplineCoeffForFittingLJPotEnergy& splineCoeff)
{
	ifstream fin( fileNameWithPathForAtomPair.c_str(), ios::in );

	if(! fin)
	{
		cout << "cannot open Spline Coefficient file of atom pair for fitting LJ Potential Energy. \n";
	}

	rg_INDEX i;
	for (i = 0;i < NUM_INTERVALS_SPLINE_FOR_FITTING_LJ_POTENERGY;i++)
	{
		rg_INDEX j;
		for (j = DEGREE_SPLINE_FOR_FITTING_LJ_POTENERGY;j >= 0;j--)
		{
			fin >> splineCoeff.m_coeff[ i ][ j ];
		}		
	}

	fin.close();

	if(! fin.good() )
	{
		cout << "A file error occurred.";
	}
}

//void FunctionsForSplineFittingLJPotEnergy::readSplineCoeffForFittingDerivativeOfLJPotEnergyForAtomPair(const string& fileNameWithPathForAtomPair, SplineCoeffForFittingLJPotEnergy& splineCoeff)
//{
//
//}


// moved to header file
//rg_REAL FunctionsForSplineFittingLJPotEnergy::evaluateLJPotEnergyBySplineFitting(const AmberAtomTypes& atomType1, 
//	                                                                             const AmberAtomTypes& atomType2, 
//	                                                                             const rg_REAL& distBTWCenters   )
//{
//	rg_INDEX index1 = (rg_INDEX) atomType1;
//	rg_INDEX index2 = (rg_INDEX) atomType2;
//
//	rg_REAL energy = 0.0;
//	energy = m_splineCoefTablefForFittingLJPotEnergy[index1][index2].evaluatePtOnSpline(distBTWCenters);
//
//	return energy;
//}

//rg_REAL FunctionsForSplineFittingLJPotEnergy::evaluateDerivativeOfLJPotEnergyBySplineFitting(const AmberAtomTypes& atomType1, 
//	                                                                                         const AmberAtomTypes& atomType2, 
//	                                                                                         const rg_REAL& distBTWCenters   )
//{
//	rg_INDEX index1 = (rg_INDEX) atomType1;
//	rg_INDEX index2 = (rg_INDEX) atomType2;
//
//	rg_REAL derivative = 0.0;
//	//derivative = m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[index1][index2].evaluateDerivativeOfLJPotEnergy(distBTWCenters);
//	derivative = m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[index1][index2].evaluatePtOnSpline(distBTWCenters);
//
//	return derivative;
//}

void FunctionsForSplineFittingLJPotEnergy::writeSplineCoeffTables(const string& fileNameWithPathOfSplineCoeffTableForFittingLJPotEnergy, const string& fileNameWithPathOfSplineCoeffTableForFittingDerivativeOfLJPotEnergy)
{
	writeSplineCoeffTableForFittingLJPotEnergyIntoBinaryFile(fileNameWithPathOfSplineCoeffTableForFittingLJPotEnergy);
	writeSplineCoeffTableForFittingDerivativeOfLJPotEnergyIntoBinaryFile(fileNameWithPathOfSplineCoeffTableForFittingDerivativeOfLJPotEnergy);
}

void FunctionsForSplineFittingLJPotEnergy::writeSplineCoeffTableForFittingLJPotEnergyIntoBinaryFile(const string& fileNameWithPath)
{
	ofstream fout( fileNameWithPath.c_str(), ios::out | ios::binary );

	if(! fout)
	{
		cout << "cannot open file. \n";
	}

	rg_INDEX i;
	for (i = 0;i < NUM_ATOM_TYPE_IN_AMBER;i++)
	{
		rg_INDEX j;
		for (j = 0;j < NUM_ATOM_TYPE_IN_AMBER;j++)
		{
			AmberAtomTypes atomType1 = (AmberAtomTypes) i;
			AmberAtomTypes atomType2 = (AmberAtomTypes) j;

			// test
			//if(atomType1 == CA_ATM && atomType2 == CT_ATM)
			//	int here = 1;

			AtomCode atomCode1 = ATOMCODE_FOR_AMBER_ATOM_TYPE[atomType1];
			AtomCode atomCode2 = ATOMCODE_FOR_AMBER_ATOM_TYPE[atomType2];
						
			SplineCoeffForFittingLJPotEnergy splineCoeff;
			if(atomCode1 != UNK_ATOM && atomCode2 != UNK_ATOM)
			{
				// get the spline coefficients from the corresponding file
				//string fileNameWithPathForAtomPair("LJ_coeff_");
				string fileNameWithPathForAtomPair("LJ_sq_dist_coeff_");
				string characterPairsOfAtomCode = getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(atomType1, atomType2);
				fileNameWithPathForAtomPair = fileNameWithPathForAtomPair + characterPairsOfAtomCode + string(".txt");			
				readSplineCoeffForAtomPair(fileNameWithPathForAtomPair, splineCoeff);
			}

			fout.write((const char*) &splineCoeff, sizeof(SplineCoeffForFittingLJPotEnergy) );
		}
	}

	fout.close();

	if(! fout.good() )
	{
		cout << "A file error occurred.";
	}
}

void FunctionsForSplineFittingLJPotEnergy::writeSplineCoeffTableForFittingDerivativeOfLJPotEnergyIntoBinaryFile(const string& fileNameWithPath)
{
	ofstream fout( fileNameWithPath.c_str(), ios::out | ios::binary );

	if(! fout)
	{
		cout << "cannot open file. \n";
	}

	rg_INDEX i;
	for (i = 0;i < NUM_ATOM_TYPE_IN_AMBER;i++)
	{
		rg_INDEX j;
		for (j = 0;j < NUM_ATOM_TYPE_IN_AMBER;j++)
		{
			AmberAtomTypes atomType1 = (AmberAtomTypes) i;
			AmberAtomTypes atomType2 = (AmberAtomTypes) j;

			AtomCode atomCode1 = ATOMCODE_FOR_AMBER_ATOM_TYPE[atomType1];
			AtomCode atomCode2 = ATOMCODE_FOR_AMBER_ATOM_TYPE[atomType2];

			SplineCoeffForFittingLJPotEnergy splineCoeff;
			if(atomCode1 != UNK_ATOM && atomCode2 != UNK_ATOM)
			{
				// get the spline coefficients from the corresponding file
				//string fileNameWithPathForAtomPair("LJ_der_coeff_");
				string fileNameWithPathForAtomPair("LJ_der_sq_dist_coeff_");
				string characterPairsOfAtomCode = getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(atomType1, atomType2);
				fileNameWithPathForAtomPair = fileNameWithPathForAtomPair + characterPairsOfAtomCode + string(".txt");			
				readSplineCoeffForAtomPair(fileNameWithPathForAtomPair, splineCoeff);
			}

			fout.write((const char*) &splineCoeff, sizeof(SplineCoeffForFittingLJPotEnergy) );
		}
	}

	fout.close();

	if(! fout.good() )
	{
		cout << "A file error occurred.";
	}
}

void FunctionsForSplineFittingLJPotEnergy::loadSplineCoeffTables(const string& fileNameWithPathOfSplineCoeffTableForFittingLJPotEnergy, const string& fileNameWithPathOfSplineCoeffTableForFittingDerivativeOfLJPotEnergy)
{
    if(m_bSplineFittingCoeffLoaded)
        return;
    else
        m_bSplineFittingCoeffLoaded = rg_TRUE;

	loadSplineCoeffTableForFittingLJPotEnergy(fileNameWithPathOfSplineCoeffTableForFittingLJPotEnergy);
	loadSplineCoeffTableForFittingDerivativeOfLJPotEnergy(fileNameWithPathOfSplineCoeffTableForFittingDerivativeOfLJPotEnergy);
}

void FunctionsForSplineFittingLJPotEnergy::loadSplineCoeffTableForFittingLJPotEnergy(const string& fileNameWithPath)
{
	ifstream fin( fileNameWithPath.c_str(), ios::in | ios::binary );

	if(! fin)
	{
		cout << "cannot open Spline Coefficient Table file for fitting LJ Potential Energy. \n";
	}

	rg_INDEX i;
	for (i = 0;i < NUM_ATOM_TYPE_IN_AMBER;i++)
	{
		rg_INDEX j;
		for (j = 0;j < NUM_ATOM_TYPE_IN_AMBER;j++)
		{
			// test
			//AmberAtomTypes atomType1 = (AmberAtomTypes) i;
			//AmberAtomTypes atomType2 = (AmberAtomTypes) j;
			//if(atomType1 == CA_ATM && atomType2 == CT_ATM)
			//	int here = 1;

			fin.read((char*) &m_splineCoefTablefForFittingLJPotEnergy[ i ][ j ], sizeof(SplineCoeffForFittingLJPotEnergy) );
		}		
	}

	fin.close();

	if(! fin.good() )
	{
		cout << "A file error occurred.";
	}
}

void FunctionsForSplineFittingLJPotEnergy::loadSplineCoeffTableForFittingDerivativeOfLJPotEnergy(const string& fileNameWithPath)
{
	ifstream fin( fileNameWithPath.c_str(), ios::in | ios::binary );

	if(! fin)
	{
		cout << "cannot open Spline Coefficient Table file for fitting Derivative of LJ Potential Energy. \n";
	}

	rg_INDEX i;
	for (i = 0;i < NUM_ATOM_TYPE_IN_AMBER;i++)
	{
		rg_INDEX j;
		for (j = 0;j < NUM_ATOM_TYPE_IN_AMBER;j++)
		{
			fin.read((char*) &m_splineCoefTablefForFittingDerivativeOfLJPotEnergy[ i ][ j ], sizeof(SplineCoeffForFittingLJPotEnergy) );
		}		
	}

	fin.close();

	if(! fin.good() )
	{
		cout << "A file error occurred.";
	}
}

string FunctionsForSplineFittingLJPotEnergy::getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(const AmberAtomTypes& atomType1, const AmberAtomTypes& atomType2)
{
    char oneCharacterForAtomCode1 = ONE_CHARACTER_OF_ATOMCODE[atomType1];
    char oneCharacterForAtomCode2 = ONE_CHARACTER_OF_ATOMCODE[atomType2];
    if(oneCharacterForAtomCode1 > oneCharacterForAtomCode2)
    {
        char temp = oneCharacterForAtomCode1;
        oneCharacterForAtomCode1 = oneCharacterForAtomCode2;
        oneCharacterForAtomCode2= temp;
    }

    string characterPairOfAtomCodePair;
    characterPairOfAtomCodePair += oneCharacterForAtomCode1;
    characterPairOfAtomCodePair += oneCharacterForAtomCode2;

    return characterPairOfAtomCodePair;
}