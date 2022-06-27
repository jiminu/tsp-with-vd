#ifndef _FUNCTIONSFORMAPPINGFROMXVOLTOLJPOTENERGY_H_
#define _FUNCTIONSFORMAPPINGFROMXVOLTOLJPOTENERGY_H_

#include "XVolNLJPotEnergyValuePairs.h"
#include "ForceField.h"

class FunctionsForMappingFromXVolToLJPotEnergy
{
public:
	static XVolNLJPotEnergyValuePairs m_XVolNLJPotEnergyValuePairsTable[NUM_ATOM_TYPE_IN_AMBER][NUM_ATOM_TYPE_IN_AMBER];    

private:
    static rg_BOOL                    m_bXVOL2LJEnergyTableLoaded;
	static void computeXVolNLJPotEnergyPairs(const AtomCode& atomCode1, 
		                                     const AtomCode& atomCode2, 
											 const rg_REAL L_J_POTENTIAL_COEFFS[2], 
											 XVolNLJPotEnergyValuePairs& xVolNLJPotEnergyValuePairs);
		static rg_REAL computeXVolume(const rg_REAL& radius1,
			                          const rg_REAL& radius2,
			                          const rg_REAL& distBTWCenters);
		static rg_REAL computeLJPotEnergy(const rg_REAL L_J_POTENTIAL_COEFFS[], const rg_REAL& distBTWCenters);

public:
	//static rg_REAL getLJPotEnergyValueCorrToXVol(const AmberAtomTypes& atomType1, const AmberAtomTypes& atomType2, const rg_REAL& XVolume);
	static inline rg_REAL getLJPotEnergyValueCorrToXVol(const AmberAtomTypes& atomType1, 
		                                                const AmberAtomTypes& atomType2, 
												        const rg_REAL& distBTWCenters  ,
												        const rg_REAL& XVolume)
	{
		rg_INDEX index1 = (rg_INDEX) atomType1;
		rg_INDEX index2 = (rg_INDEX) atomType2;

		rg_REAL LJPotEnergyValue = 0.0;
		LJPotEnergyValue = m_XVolNLJPotEnergyValuePairsTable[index1][index2].estimateLJPotEnergyValueOfXVolume( distBTWCenters, XVolume );

		return LJPotEnergyValue;
	}

	static void    writeXVolNLJPotEnergyValuePairTableIntoBinaryFile(const string& fileNameWithPath);
	static void    loadXVolNLJPotEnergyValuePairTables(const string& fileNameWithPath);
	
    // For testing
    static void    writeXVolNLJPotEnergyValuePairsIntoSeparateFiles(const string& pathForFiles, const rg_BOOL& bIncreasingOrderOfDistance);

    static void    writeXVolNLJPotEnergyValuePairs(const AmberAtomTypes& atomType1, 
                                                   const AmberAtomTypes& atomType2, 
                                                   const string& fileNameWithPath,
                                                   const rg_BOOL& bIncreasingOrderOfDistance);

    //static void    writeXVolNLJPotEnergyValuePairs(const AmberAtomTypes& atomType1, 
    //	                                           const AmberAtomTypes& atomType2, 
    //											   const string& fileNameWithPath);


	static void    writeXVolNLJPotEnergyValuePairs(const XVolNLJPotEnergyValuePairs& xVolNLJPotEnergyValuePairs,
		                                           const string& fileNameWithPath);

    static void    writeLJCoeffForAtomPairs(const string& fileNameWithPath);

};

#endif