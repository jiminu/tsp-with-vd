#include "FunctionsForMappingFromXVolToLJPotEnergy.h"
#include "rg_GeoFunc.h"
#include <iostream>
#include <fstream>

#include "FunctionsForSplineFittingLJPotEnergy.h"

XVolNLJPotEnergyValuePairs FunctionsForMappingFromXVolToLJPotEnergy::m_XVolNLJPotEnergyValuePairsTable[NUM_ATOM_TYPE_IN_AMBER][NUM_ATOM_TYPE_IN_AMBER];
rg_BOOL                    FunctionsForMappingFromXVolToLJPotEnergy::m_bXVOL2LJEnergyTableLoaded = rg_FALSE;

void FunctionsForMappingFromXVolToLJPotEnergy::computeXVolNLJPotEnergyPairs(const AtomCode& atomCode1, 
	                                                                        const AtomCode& atomCode2, 
																		    const rg_REAL L_J_POTENTIAL_COEFFS[2], 
																		    XVolNLJPotEnergyValuePairs& xVolNLJPotEnergyValuePairs)
{
	rg_REAL radius1 = ATOM_FEATURES[atomCode1].radius;
	rg_REAL radius2 = ATOM_FEATURES[atomCode2].radius;
	rg_REAL sumOfTwoRadii = radius1 + radius2;
	rg_REAL distBTWCenters = 0.0;

	rg_INDEX i;
	//for (i = 0;i < NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
	for (i = 1;i <= NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
	{
		//distBTWCenters = sumOfTwoRadii * (1.0 - (rg_REAL) i / NUM_XVOLLJPOTENERGYVAL_PAIRS);
		distBTWCenters = sumOfTwoRadii * (rg_REAL) i / NUM_XVOLLJPOTENERGYVAL_PAIRS;

		//xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ i ][ 0 ] = computeXVolume(radius1, radius2, distBTWCenters);
		//xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ i ][ 1 ] = computeLJPotEnergy(L_J_POTENTIAL_COEFFS, distBTWCenters);
		xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ i-1 ][ 0 ] = computeXVolume(radius1, radius2, distBTWCenters);
		xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ i-1 ][ 1 ] = computeLJPotEnergy(L_J_POTENTIAL_COEFFS, distBTWCenters);

		// i,  interval, R = sumOfTwoRadii, N = NUM_XVOLLJPOTENERGYVAL_PAIRS
		// 0,  [     R/N, 2R/N) --> m_XVolEnergyPairs[ 0 ][ 0 ] = xvol for distance =  R/N
		// 1,  [    2R/N, 3R/N) --> m_XVolEnergyPairs[ 1 ][ 0 ] = xvol for distance = 2R/N
		// 2,  [    3R/N, 4R/N) --> m_XVolEnergyPairs[ 2 ][ 0 ] = xvol for distance = 3R/N
		// .
		// .
		// .
		// N-2 [(N-1)R/N,  R  ) --> m_XVolEnergyPairs[N-2][ 0 ] = xvol for distance = (N-1)R/N
		//                      --> m_XVolEnergyPairs[N-1][ 0 ] = xvol for distance =      R
	}
	xVolNLJPotEnergyValuePairs.m_sumOfTwoRadii = sumOfTwoRadii;
	xVolNLJPotEnergyValuePairs.m_sumOfTwoRadii_over_NUM_XVOLLJPOTENERGYVAL_PAIRS = sumOfTwoRadii / NUM_XVOLLJPOTENERGYVAL_PAIRS;
	xVolNLJPotEnergyValuePairs.m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_diffMinMaxDist = (NUM_XVOLLJPOTENERGYVAL_PAIRS * NUM_XVOLLJPOTENERGYVAL_PAIRS) / (sumOfTwoRadii * (NUM_XVOLLJPOTENERGYVAL_PAIRS - 1));
	xVolNLJPotEnergyValuePairs.m_NUM_XVOLLJPOTENERGYVAL_PAIRS_over_sumOfTwoRadii = NUM_XVOLLJPOTENERGYVAL_PAIRS / sumOfTwoRadii;
	//rg_REAL lowerBndOfXVol = xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ 0 ][ 0 ];
	//rg_REAL upperBndOfXVol = xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ NUM_XVOLLJPOTENERGYVAL_PAIRS - 1][ 0 ];
	//xVolNLJPotEnergyValuePairs.m_lengthOfXVolumeRange = (upperBndOfXVol - lowerBndOfXVol);
}

rg_REAL FunctionsForMappingFromXVolToLJPotEnergy::computeXVolume(const rg_REAL& radius1,
	                                                             const rg_REAL& radius2,
	                                                             const rg_REAL& distBTWCenters)
{
	return rg_GeoFunc::computeXVolumeOfTwoSpheres(radius1, radius2, distBTWCenters);
}

rg_REAL FunctionsForMappingFromXVolToLJPotEnergy::computeLJPotEnergy(const rg_REAL L_J_POTENTIAL_COEFFS[], const rg_REAL& distBTWCenters)
{
	rg_REAL squaredDistance = 0.0;
	squaredDistance = distBTWCenters * distBTWCenters;
	rg_REAL sixthPoweredDistance = 0.0;
	sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

	rg_REAL LJPotEnergy = 0.0;
	LJPotEnergy = ( L_J_POTENTIAL_COEFFS[0] / (sixthPoweredDistance*sixthPoweredDistance) ) - ( L_J_POTENTIAL_COEFFS[1] / sixthPoweredDistance ) ;

	return LJPotEnergy;
}


//rg_REAL FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(const AmberAtomTypes& atomType1, const AmberAtomTypes& atomType2, const rg_REAL& XVolume)
//{
//	rg_INDEX index1 = (rg_INDEX) atomType1;
//	rg_INDEX index2 = (rg_INDEX) atomType2;
//
//	rg_REAL LJPotEnergyValue = 0.0;
//	LJPotEnergyValue = m_XVolNLJPotEnergyValuePairsTable[index1][index2].estimateLJPotEnergyValueOfXVolume( XVolume );
//
//	return LJPotEnergyValue;
//}



// moved to header file
//rg_REAL FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(const AmberAtomTypes& atomType1, 
//	                                                                            const AmberAtomTypes& atomType2, 
//																				const rg_REAL& distBTWCenters  ,
//																				const rg_REAL& XVolume          )
//{
//		rg_INDEX index1 = (rg_INDEX) atomType1;
//		rg_INDEX index2 = (rg_INDEX) atomType2;
//	
//		rg_REAL LJPotEnergyValue = 0.0;
//		LJPotEnergyValue = m_XVolNLJPotEnergyValuePairsTable[index1][index2].estimateLJPotEnergyValueOfXVolume( distBTWCenters, XVolume );
//	
//		return LJPotEnergyValue;
//}

void    FunctionsForMappingFromXVolToLJPotEnergy::writeXVolNLJPotEnergyValuePairTableIntoBinaryFile(const string& fileNameWithPath)
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

#ifdef _DEBUG
            if(atomCode1 == C_ATOM && atomCode2 == O_ATOM)
                int here = 1;
#endif

			XVolNLJPotEnergyValuePairs xVolNLJPotEnergyValuePairs;
			if(atomCode1 != UNK_ATOM && atomCode2 != UNK_ATOM)
			{
				computeXVolNLJPotEnergyPairs(atomCode1,
					                         atomCode2,
											 L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[ i ][ j ], 
											 xVolNLJPotEnergyValuePairs);
			}

			fout.write((const char*) &xVolNLJPotEnergyValuePairs, sizeof(XVolNLJPotEnergyValuePairs) );
		}
	}	

	fout.close();

	if(! fout.good() )
	{
		cout << "A file error occurred.";
	}
}

void    FunctionsForMappingFromXVolToLJPotEnergy::loadXVolNLJPotEnergyValuePairTables(const string& fileNameWithPath)
{
    if(m_bXVOL2LJEnergyTableLoaded)
        return;
    else 
        m_bXVOL2LJEnergyTableLoaded = rg_TRUE;

	ifstream fin( fileNameWithPath.c_str(), ios::in | ios::binary );

	if(! fin)
	{
		cout << "cannot open XVolumeNLJPotENergyTagble file. \n";
		// May 26, 2016 by Joonghyun
		exit( 1 );
	}

	rg_INDEX i;
	for (i = 0;i < NUM_ATOM_TYPE_IN_AMBER;i++)
	{
		rg_INDEX j;
		for (j = 0;j < NUM_ATOM_TYPE_IN_AMBER;j++)
		{
			fin.read((char*) &m_XVolNLJPotEnergyValuePairsTable[ i ][ j ], sizeof(XVolNLJPotEnergyValuePairs) );
		}		
	}

	fin.close();

	if(! fin.good() )
	{
		//cout << "A file error occurred.";
		// May 26, 2016 by Joonghyun
		cout << "An error occurred during loading intersection volume and LJ-energy level pair table.";
		exit( 1 );
	}
}

void    FunctionsForMappingFromXVolToLJPotEnergy::writeXVolNLJPotEnergyValuePairsIntoSeparateFiles(const string& pathForFiles, const rg_BOOL& bIncreasingOrderOfDistance)
{
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

            if(atomCode1 != UNK_ATOM && atomCode2 != UNK_ATOM)
            {
                string fileNameWithPathForAtomPair("XVol2EnergyTable_");
                fileNameWithPathForAtomPair = pathForFiles + fileNameWithPathForAtomPair;

                string characterPairsOfAtomCode = FunctionsForSplineFittingLJPotEnergy::getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(atomType1, atomType2);
                fileNameWithPathForAtomPair = fileNameWithPathForAtomPair + characterPairsOfAtomCode + string(".txt");
                writeXVolNLJPotEnergyValuePairs(atomType1, atomType2, fileNameWithPathForAtomPair, bIncreasingOrderOfDistance);
            }
        }
    }
}

void    FunctionsForMappingFromXVolToLJPotEnergy::writeXVolNLJPotEnergyValuePairs(const AmberAtomTypes& atomType1, 
                                                                                  const AmberAtomTypes& atomType2, 
                                                                                  const string& fileNameWithPath,
                                                                                  const rg_BOOL& bIncreasingOrderOfDistance)
{
    rg_INDEX index1 = (rg_INDEX) atomType1;
    rg_INDEX index2 = (rg_INDEX) atomType2;

    ofstream fout( fileNameWithPath.c_str());

    if(! bIncreasingOrderOfDistance)
    {
        rg_INDEX i;
        for(i = NUM_XVOLLJPOTENERGYVAL_PAIRS - 1;i >= 0;i--)
            fout << m_XVolNLJPotEnergyValuePairsTable[index1][index2].m_XVolEnergyPairs[ i ][ 0 ] << " "
            << m_XVolNLJPotEnergyValuePairsTable[index1][index2].m_XVolEnergyPairs[ i ][ 1 ] << endl;
    }
    else
    {
        rg_INDEX i;
        for(i = 0;i < NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
            fout << m_XVolNLJPotEnergyValuePairsTable[index1][index2].m_XVolEnergyPairs[ i ][ 0 ] << " "
            << m_XVolNLJPotEnergyValuePairsTable[index1][index2].m_XVolEnergyPairs[ i ][ 1 ] << endl;
    }


    fout.close();
}

//void    FunctionsForMappingFromXVolToLJPotEnergy::writeXVolNLJPotEnergyValuePairs(const AmberAtomTypes& atomType1, 
//	                                                                              const AmberAtomTypes& atomType2, 
//																				  const string& fileNameWithPath)
//{
//	rg_INDEX index1 = (rg_INDEX) atomType1;
//	rg_INDEX index2 = (rg_INDEX) atomType2;
//
//	ofstream fout( fileNameWithPath.c_str());
//
//	rg_INDEX i;
//	for(i = 0;i < NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
//		fout << m_XVolNLJPotEnergyValuePairsTable[index1][index2].m_XVolEnergyPairs[ i ][ 0 ] << " "
//		     << m_XVolNLJPotEnergyValuePairsTable[index1][index2].m_XVolEnergyPairs[ i ][ 1 ] << endl;
//
//	fout.close();
//}

void    FunctionsForMappingFromXVolToLJPotEnergy::writeXVolNLJPotEnergyValuePairs(const XVolNLJPotEnergyValuePairs& xVolNLJPotEnergyValuePairs,
	                                                                              const string& fileNameWithPath)
{
	ofstream fout( fileNameWithPath.c_str());

	rg_INDEX i;
	for(i = 0;i < NUM_XVOLLJPOTENERGYVAL_PAIRS;i++)
		fout << xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ i ][ 0 ] << " "
		     << xVolNLJPotEnergyValuePairs.m_XVolEnergyPairs[ i ][ 1 ] << endl;

	fout.close();
}

void   FunctionsForMappingFromXVolToLJPotEnergy::writeLJCoeffForAtomPairs(const string& fileNameWithPath)
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

            if(atomCode1 != UNK_ATOM && atomCode2 != UNK_ATOM)
            {
                string characterPairsOfAtomCode = FunctionsForSplineFittingLJPotEnergy::getCharacterPiarOfAtomCodePairForReadingSplineCoeffFileCorrTo(atomType1, atomType2);
                fout << atomType1 << " " << atomType2 << " " << characterPairsOfAtomCode << " " << L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[ i ][ j ][ 0 ] << " " << L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[ i ][ j ][ 1 ] << endl;
            }
        }
    }	

    fout.close();

    if(! fout.good() )
    {
        cout << "A file error occurred.";
    }
}
	