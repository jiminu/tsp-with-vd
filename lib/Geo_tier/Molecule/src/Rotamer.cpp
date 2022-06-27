#include "Rotamer.h"
#include "Residue.h"
#include "FunctionsForMolecule.h"
#include "FunctionsForRotamerLibrary.h"
#include "rg_TMatrix3D.h"

void Rotamer::setDihedralAngles(V::GeometryTier::Residue* residue)
{
	if(residue == rg_NULL)
		return;

	destroyDihedralAngles();

	m_numDihedralAngles
		= FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(m_residue->getResidueCode());

	if (m_numDihedralAngles > 0)
	{
		m_dihedralAngles = new rg_REAL[m_numDihedralAngles];
		residue->getSideChainDihedralAngles(m_dihedralAngles);
	}
}

Rotamer::Rotamer()
{
	setInitialValues();
}

Rotamer::Rotamer(V::GeometryTier::Residue* residue)
{
	setInitialValues();
	setResidue( residue );
}

Rotamer::Rotamer(const Rotamer& rotamer)
{
	setInitialValues();
	m_residue = rotamer.m_residue;
	setDihedralAngles(rotamer.m_dihedralAngles);
}

Rotamer::~Rotamer()
{
	m_residue    = rg_NULL;
	destroyDihedralAngles();
}

rg_REAL Rotamer::computeEnergyWithOtherRotamer( Rotamer& theOtherRotamer )
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);

	rg_INT numTheOtherRotamerAtoms = 0;
    V::GeometryTier::Atom** theOtherRotamerAtoms = rg_NULL;
	numTheOtherRotamerAtoms = theOtherRotamer.m_residue->getAtomsOnSidechain(theOtherRotamerAtoms);
	//numTheOtherRotamerAtoms = theOtherRotamer.getAtomsOnSidechain(theOtherRotamerAtoms);

	
    rg_INT numXAtoms = 0;
    // test
    energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, numRotamerAtoms, theOtherRotamerAtoms, numTheOtherRotamerAtoms);
        
    //energy = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(rotamerAtoms, numRotamerAtoms, theOtherRotamerAtoms, numTheOtherRotamerAtoms, numXAtoms);

    //energy =  computeNewScoreFunction(rotamerAtoms, numRotamerAtoms, theOtherRotamerAtoms, numTheOtherRotamerAtoms, numXAtoms);
    
	//energy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting(rotamerAtoms, numRotamerAtoms, theOtherRotamerAtoms, numTheOtherRotamerAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	if(theOtherRotamerAtoms != rg_NULL)
		delete [] theOtherRotamerAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithBackbone    ( Backbone     & backbone     )
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numBackboneAtoms = 0;
    V::GeometryTier::Atom** backboneAtoms = rg_NULL;
	numBackboneAtoms = backbone.getBackboneAtoms( backboneAtoms );

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);

    // test
    rg_INT numXAtoms = 0;
    //energy = computeNewScoreFunction(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithBackbone    ( Backbone     * backbone     )
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numBackboneAtoms = 0;
    V::GeometryTier::Atom** backboneAtoms = rg_NULL;
	numBackboneAtoms = backbone->getBackboneAtoms( backboneAtoms );

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);

	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	return energy;
}


rg_REAL Rotamer::computeEnergyWithOtherBackbone(V::GeometryTier::Residue* residue)
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);
	numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);

	rg_INT numTheOtherResidueBackboneAtoms = 0;
    V::GeometryTier::Atom** theOtherResidueBackboneAtoms = rg_NULL;
	numTheOtherResidueBackboneAtoms = residue->getAtomsOnBackbone(theOtherResidueBackboneAtoms);

    // test
    energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, 
		                                                numRotamerAtoms, 
														theOtherResidueBackboneAtoms, 
														numTheOtherResidueBackboneAtoms);

    rg_INT numXAtoms = 0;

    //energy = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(rotamerAtoms, 
    //                                                                             numRotamerAtoms, 
    //                                                                             theOtherResidueBackboneAtoms, 
    //                                                                             numTheOtherResidueBackboneAtoms,
    //                                                                             numXAtoms);

    //energy =  computeNewScoreFunction(rotamerAtoms, 
    //                                  numRotamerAtoms, 
    //                                  theOtherResidueBackboneAtoms, 
    //                                  numTheOtherResidueBackboneAtoms,
    //                                  numXAtoms);

	//energy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting(rotamerAtoms,
	//	                                                               numRotamerAtoms,
	//																   theOtherResidueBackboneAtoms,
	//																   numTheOtherResidueBackboneAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	if(theOtherResidueBackboneAtoms != rg_NULL)
		delete [] theOtherResidueBackboneAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithOtherResidue     (V::GeometryTier::Residue* residue)
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);
	numRotamerAtoms = m_residue->getAtomsOnSidechain( rotamerAtoms );

	rg_INT numTheOtherResidueAtoms = 0;
    V::GeometryTier::Atom** theOtherResidueAtoms = rg_NULL;
	numTheOtherResidueAtoms = residue->getAtoms(theOtherResidueAtoms);
	//numTheOtherResidueAtoms = FunctionsForRotamerLibrary::getAtomsOfResidue(residue, theOtherResidueAtoms);


    // test
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, numRotamerAtoms, theOtherResidueAtoms, numTheOtherResidueAtoms);
    
    rg_INT numXAtoms = 0;
    //energy = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(rotamerAtoms, numRotamerAtoms, theOtherResidueAtoms, numTheOtherResidueAtoms, numXAtoms);

    //energy = computeNewScoreFunction(rotamerAtoms, numRotamerAtoms, theOtherResidueAtoms, numTheOtherResidueAtoms, numXAtoms);

	//energy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting(rotamerAtoms, numRotamerAtoms, theOtherResidueAtoms, numTheOtherResidueAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	if(theOtherResidueAtoms != rg_NULL)
		delete [] theOtherResidueAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithItsBackbone      ()
{
    rg_REAL energy = 0.0;

    ResidueCode code = m_residue->getResidueCode();
    if(code == GLY_AMINO_RESIDUE)
        return energy;

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);
	numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);

	rg_INT numBackboneAtoms = 0;
    V::GeometryTier::Atom** backboneAtoms = rg_NULL;
	numBackboneAtoms = m_residue->getAtomsOnBackbone(backboneAtoms);


    // test
	
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms);

    rg_INT numXAtoms = 0;

    //energy = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);

    //energy = computeNewScoreFunction(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);

	//energy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	if(backboneAtoms != rg_NULL)
		delete [] backboneAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithItsBackboneAtomsNOtherResidue     (V::GeometryTier::Residue* residue)
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	//numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);
	numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);
	
	rg_dList<V::GeometryTier::Atom*> backboneAtoms;
	m_residue->getAtomsOnBackbone(&backboneAtoms);
	rg_dList<V::GeometryTier::Atom*> theOtherResidueAtoms;
	//FunctionsForRotamerLibrary::getAtomsOfResidue(residue, theOtherResidueAtoms);
	residue->getAtoms(theOtherResidueAtoms);
	theOtherResidueAtoms.append( backboneAtoms );

	rg_INT numBackboneNTheOtherResidueAtoms = 0;
	numBackboneNTheOtherResidueAtoms = theOtherResidueAtoms.getSize();
    V::GeometryTier::Atom** backboneNOtherResidueAtoms = new V::GeometryTier::Atom* [numBackboneNTheOtherResidueAtoms];

	rg_INDEX index = 0;
	theOtherResidueAtoms.reset4Loop();
	while (theOtherResidueAtoms.setNext4Loop())
	{
        V::GeometryTier::Atom* currAtom = theOtherResidueAtoms.getEntity();
		backboneNOtherResidueAtoms[ index ] = currAtom;
		index++;
	}

	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(rotamerAtoms, numRotamerAtoms, backboneNOtherResidueAtoms, numBackboneNTheOtherResidueAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	if(backboneNOtherResidueAtoms != rg_NULL)
		delete [] backboneNOtherResidueAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithOtherRotamerBySCWRL3( Rotamer& theOtherRotamer )
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	//numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);
	numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);

	rg_INT numTheOtherRotamerAtoms = 0;
    V::GeometryTier::Atom** theOtherRotamerAtoms = rg_NULL;
	//numTheOtherRotamerAtoms = theOtherRotamer.m_residue->getAtomsOnSidechain(theOtherRotamerAtoms);
	numTheOtherRotamerAtoms = theOtherRotamer.getAtomsOnSidechain(theOtherRotamerAtoms);

    rg_INT numXAtoms = 0;
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3(rotamerAtoms, numRotamerAtoms, theOtherRotamerAtoms, numTheOtherRotamerAtoms, numXAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	return energy;
}

rg_REAL Rotamer::computeEnergyWithBackboneBySCWRL3    ( Backbone     & backbone     )
{
	rg_REAL energy = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return energy;

	rg_INT numBackboneAtoms = 0;
    V::GeometryTier::Atom** backboneAtoms = rg_NULL;
	numBackboneAtoms = backbone.getBackboneAtoms( backboneAtoms );

	rg_INT numRotamerAtoms = 0;
    V::GeometryTier::Atom** rotamerAtoms = rg_NULL;
	//numRotamerAtoms = m_residue->getAtomsOnSidechain(rotamerAtoms);
	numRotamerAtoms = getAtomsOnSidechain(rotamerAtoms);

    rg_INT numXAtoms = 0;
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3(rotamerAtoms, numRotamerAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);

	if(rotamerAtoms != rg_NULL)
		delete [] rotamerAtoms;

	return energy;
}


rg_REAL Rotamer::computeRotamerPreferenceEnergyBySCWRL3(const rg_INDEX& rotLibID, const rg_INDEX& highestProbabilityRotLibID)
{
	ResidueCode code = getResidue()->getResidueCode();
	if (code == ALA_AMINO_RESIDUE || code == GLY_AMINO_RESIDUE)
	{
		return 0.0;
	}

	rg_REAL scalingFactor = 1.0;
	if(code == HIS_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE)
	{
		scalingFactor = 2.0;
	}

	rg_INT numDihedralAngles = getNumDihedralAngles();	
	rg_REAL* dihedralAngles  = getDihedralAngles();

	rg_REAL devFromRotLib = 0.0;
	rg_INDEX i;
	for (i = 0;i < numDihedralAngles;i++)
	{		
		devFromRotLib += ( ((dihedralAngles[ i ] - FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].chi[ i ]) / FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].stdDeviation[ i ]) * 
			               ((dihedralAngles[ i ] - FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].chi[ i ]) / FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].stdDeviation[ i ])   );
	}

	rg_REAL probability        = FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].probability;
	rg_REAL highestProbability = FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[highestProbabilityRotLibID].probability;
		
	rg_REAL logNormalizedProbability = 0.0;
	if(probability > 0)
		logNormalizedProbability = -log(probability / highestProbability);

	rg_REAL energy = 0.0;
	const rg_REAL k_tor = 0.4;
	energy = scalingFactor * (logNormalizedProbability + k_tor * devFromRotLib);

	return energy;
}

void Rotamer::evaluateGradientVectorOfRotamerPreferenceEnergyForDihedralAnglesBySCWRL3( Vector& givenSolVector, const rg_INDEX& rotLibID, Vector& gradientVector )
{
	rg_INT numDihedralAngles = getNumDihedralAngles();
	gradientVector.setSize( numDihedralAngles );

	ResidueCode code = getResidue()->getResidueCode();
	if (code == ALA_AMINO_RESIDUE || code == GLY_AMINO_RESIDUE)
	{
		return;
	}

	rg_REAL scalingFactor = 1.0;
	if(code == HIS_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE)
	{
		scalingFactor = 2.0;
	}

	const rg_REAL two_times_k_tor_scalingFactor = 2.0 * 0.4 * scalingFactor;
	rg_REAL* dihedralAngles  = getDihedralAngles();

	rg_INDEX i;
	for (i = 0;i < numDihedralAngles;i++)
	{
		rg_REAL stdDeviation = FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].stdDeviation[ i ];
		gradientVector[ i ] = two_times_k_tor_scalingFactor * (dihedralAngles[ i ] - FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[rotLibID].chi[ i ]) / (stdDeviation * stdDeviation);
	}
}

void Rotamer::computeMinimumEnclosingSphere(EnclosingSphereOfSpheres& MES_R)
{
	rg_dList<V::GeometryTier::Atom*> atomsOfRotamer;
	ResidueCode code = m_residue->getResidueCode();
	if (code != ALA_AMINO_RESIDUE && code != GLY_AMINO_RESIDUE)
	{		
		//getAtomsOnSidechain(atomsOfRotamer);		
		m_residue->getAtomsOnSidechain(atomsOfRotamer);
	}
	else if(code == ALA_AMINO_RESIDUE)
	{
        V::GeometryTier::Atom* betaCarbon = m_residue->getBetaCarbonInSideChainOfAminoResidue();
		atomsOfRotamer.add(betaCarbon);
	}
	else if(code == GLY_AMINO_RESIDUE)
	{		
        V::GeometryTier::Atom* alphaCarbon = m_residue->getAlphaCarbonOfAminoResidue();
		atomsOfRotamer.add(alphaCarbon);		
	}
	
	MES_R.setSpheres(atomsOfRotamer);
	MES_R.computeEnclosingSphere();
}

void Rotamer::computeMinimumEnclosingSphereOfResidue(EnclosingSphereOfSpheres& MES_R)
{
	rg_dList<V::GeometryTier::Atom*> atomsOfRotamer;
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	//FunctionsForRotamerLibrary::getAtomsOfResidue(m_residue, atomsOfRotamer);
	m_residue->getAtoms(atomsOfRotamer);
	MES_R.setSpheres(atomsOfRotamer);
	MES_R.computeEnclosingSphere();
}

rg_REAL Rotamer::computeSumOfPairwiseXVolumeWithOtherRotamer(Rotamer* rotamer)
{
	rg_REAL sumOfPairwiseXVolume = 0.0;

	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	//getAtomsOnSidechain(targetAtomsOfSidechain);
	m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	rg_dList<V::GeometryTier::Atom*> atomsOfSidechain;
	rotamer->getAtomsOnSidechain(atomsOfSidechain);

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{
		rg_REAL xVolume = 0.0;
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();
		xVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfSidechain);
		sumOfPairwiseXVolume += xVolume;
	}
	return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseXVolumeWithOtherRotamer(Rotamer* rotamer, rg_INT& numXAtoms)
{
    V::GeometryTier::Atom** targetAtomsOfSidechain = rg_NULL;
    rg_INT numTargetSidechainAtoms = m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);
    V::GeometryTier::Atom** atomsOfSidechain = rg_NULL;
    rg_INT numSidechainAtoms = rotamer->m_residue->getAtomsOnSidechain(atomsOfSidechain);

    rg_REAL sumOfPairwiseXVolume = 0.0;
    numXAtoms = 0;

    //test
    //sumOfPairwiseXVolume = computeVanDerWaalsEnergyBetweenTwoAtomSets(targetAtomsOfSidechain, numTargetSidechainAtoms, atomsOfSidechain, numSidechainAtoms);
    sumOfPairwiseXVolume = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(targetAtomsOfSidechain, numTargetSidechainAtoms, atomsOfSidechain, numSidechainAtoms, numXAtoms);

    if(targetAtomsOfSidechain != rg_NULL)
        delete [] targetAtomsOfSidechain;
    if(atomsOfSidechain != rg_NULL)
        delete [] atomsOfSidechain;
    
    return sumOfPairwiseXVolume;

	//rg_dList<Atom*> targetAtomsOfSidechain;
	////getAtomsOnSidechain(targetAtomsOfSidechain);
	//m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	//rg_dList<Atom*> atomsOfSidechain;
	////rotamer->getAtomsOnSidechain(atomsOfSidechain);
	//rotamer->m_residue->getAtomsOnSidechain(atomsOfSidechain);

	//rg_REAL sumOfPairwiseXVolume = 0.0;
	//numXAtoms = 0;

	//targetAtomsOfSidechain.reset4Loop();
	//while (targetAtomsOfSidechain.setNext4Loop())
	//{
	//	rg_REAL currXVolume = 0.0;
	//	rg_INT  currNumXAtoms = 0;

	//	Atom* currAtom = targetAtomsOfSidechain.getEntity();
	//	currXVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfSidechain, currNumXAtoms);

	//	sumOfPairwiseXVolume += currXVolume;
	//	numXAtoms += currNumXAtoms;
	//}
	//return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseXVolumeWithOtherBackbone(V::GeometryTier::Residue* residue, rg_INT& numXAtoms)
{
    V::GeometryTier::Atom** targetAtomsOfSidechain = rg_NULL;
    rg_INT numTargetAtoms = m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);
    V::GeometryTier::Atom** atomsOfBackbone = rg_NULL;
    rg_INT numBackboneAtoms = residue->getAtomsOnBackbone(atomsOfBackbone);

    rg_REAL sumOfPairwiseXVolume = 0.0;
    numXAtoms = 0;

    //test
    //sumOfPairwiseXVolume = computeVanDerWaalsEnergyBetweenTwoAtomSets(targetAtomsOfSidechain, numTargetAtoms, atomsOfBackbone, numBackboneAtoms);
    sumOfPairwiseXVolume = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(targetAtomsOfSidechain, numTargetAtoms, atomsOfBackbone, numBackboneAtoms, numXAtoms);

    if(targetAtomsOfSidechain != rg_NULL)
        delete [] targetAtomsOfSidechain;
    if(atomsOfBackbone != rg_NULL)
        delete [] atomsOfBackbone;

    return sumOfPairwiseXVolume;

	//rg_dList<Atom*> targetAtomsOfSidechain;
	////getAtomsOnSidechain(targetAtomsOfSidechain);
	//m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	//rg_dList<Atom*> atomsOfBackbone;
	//residue->getAtomsOnBackbone(&atomsOfBackbone);

	//rg_REAL sumOfPairwiseXVolume = 0.0;
	//numXAtoms = 0;

	//targetAtomsOfSidechain.reset4Loop();
	//while (targetAtomsOfSidechain.setNext4Loop())
	//{		
	//	Atom* currAtom = targetAtomsOfSidechain.getEntity();
	//	rg_REAL currXVolume = 0.0;
	//	rg_INT  currNumXAtoms = 0;
	//	currXVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfBackbone, currNumXAtoms);
	//	sumOfPairwiseXVolume += currXVolume;
	//	numXAtoms += currNumXAtoms;
	//}
	//return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseAtomXVolumesWithOtherResidue(V::GeometryTier::Residue* residue)
{
	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	//getAtomsOnSidechain(targetAtomsOfSidechain);
	m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	rg_dList<V::GeometryTier::Atom*> atomsOfResidue;
	//FunctionsForRotamerLibrary::getAtomsOfResidue(residue, atomsOfResidue);
	residue->getAtoms(atomsOfResidue);

	rg_REAL sumOfPairwiseXVolume = 0.0;

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{		
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();
		rg_REAL currXVolume = 0.0;
		currXVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfResidue);
		sumOfPairwiseXVolume += currXVolume;
	}
	return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseAtomXVolumesWithOtherResidue(V::GeometryTier::Residue* residue, rg_INT& numXAtoms)
{
	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	//getAtomsOnSidechain(targetAtomsOfSidechain);
	m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	rg_dList<V::GeometryTier::Atom*> atomsOfResidue;
	//FunctionsForRotamerLibrary::getAtomsOfResidue(residue, atomsOfResidue);
	residue->getAtoms(atomsOfResidue);

	rg_REAL sumOfPairwiseXVolume = 0.0;
	numXAtoms = 0;

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{		
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();
		rg_REAL currXVolume = 0.0;
		rg_INT  currNumXAtoms = 0;
		currXVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfResidue, currNumXAtoms);
		sumOfPairwiseXVolume += currXVolume;
		numXAtoms += currNumXAtoms;
	}
	return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseAtomXVolumesWithItsBackbone()
{
	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	//getAtomsOnSidechain(targetAtomsOfSidechain);
	m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	rg_dList<V::GeometryTier::Atom*> backboneAtoms;
	m_residue->getAtomsOnBackbone(&backboneAtoms);

	rg_REAL sumOfPairwiseXVolume = 0.0;

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();

		rg_REAL currXVolume = 0.0;
		currXVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, backboneAtoms);
		sumOfPairwiseXVolume += currXVolume;
	}
	return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseAtomXVolumesWithItsBackbone(rg_INT& numXAtoms)
{
    V::GeometryTier::Atom** targetAtomsOfSidechain = rg_NULL;
    rg_INT numTargetSidechainAtoms = m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

    V::GeometryTier::Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms = m_residue->getAtomsOnBackbone(backboneAtoms);

    rg_REAL sumOfPairwiseXVolume = 0.0;
    numXAtoms = 0;

    //test
    //sumOfPairwiseXVolume = computeVanDerWaalsEnergyBetweenTwoAtomSets(targetAtomsOfSidechain, numTargetSidechainAtoms, backboneAtoms, numBackboneAtoms);
    sumOfPairwiseXVolume = computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy(targetAtomsOfSidechain, numTargetSidechainAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);

    if(targetAtomsOfSidechain != rg_NULL)
        delete [] targetAtomsOfSidechain;
    if(backboneAtoms != rg_NULL)
        delete [] backboneAtoms;

    return sumOfPairwiseXVolume;
    
	//rg_dList<Atom*> targetAtomsOfSidechain;
	////getAtomsOnSidechain(targetAtomsOfSidechain);
	//m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	//rg_dList<Atom*> backboneAtoms;
	//m_residue->getAtomsOnBackbone(&backboneAtoms);

	//rg_REAL sumOfPairwiseXVolume = 0.0;
	//numXAtoms = 0;

	//targetAtomsOfSidechain.reset4Loop();
	//while (targetAtomsOfSidechain.setNext4Loop())
	//{
	//	Atom* currAtom = targetAtomsOfSidechain.getEntity();

	//	rg_REAL currXVolume = 0.0;
	//	rg_INT currNumXAtoms = 0;
	//	currXVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, backboneAtoms, currNumXAtoms);
	//	sumOfPairwiseXVolume += currXVolume;
	//	numXAtoms += currNumXAtoms;
	//}
	//return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseXVolumeWithItsBackboneAtomsNOtherResidue(V::GeometryTier::Residue* residue)
{
	rg_REAL sumOfPairwiseXVolume = 0.0;

	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	//getAtomsOnSidechain(targetAtomsOfSidechain);
	m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);

	// intersection volume with backbone
	rg_dList<V::GeometryTier::Atom*> atomsOfBackbone;
	m_residue->getAtomsOnBackbone(&atomsOfBackbone);

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{
		rg_REAL xVolume = 0.0;
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();
		xVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfBackbone);
		sumOfPairwiseXVolume += xVolume;
	}

	// intersection volume with other residue
	rg_dList<V::GeometryTier::Atom*> atomsOfResidue;
	//FunctionsForRotamerLibrary::getAtomsOfResidue(residue, atomsOfResidue);
	residue->getAtoms(atomsOfResidue);
	//atomsOfResidue.append( atomsOfBackbone );

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{
		rg_REAL xVolume = 0.0;
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();
		xVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfResidue);
		sumOfPairwiseXVolume += xVolume;
	}

	return sumOfPairwiseXVolume;
}

rg_REAL Rotamer::computeSumOfPairwiseXVolumeWithProteinBackbone(Backbone& backbone)
{
	rg_REAL sumOfPairwiseXVolume = 0.0;

	ResidueCode code = m_residue->getResidueCode();
	if(code == GLY_AMINO_RESIDUE)
		return sumOfPairwiseXVolume;

	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	//getAtomsOnSidechain(targetAtomsOfSidechain);
	m_residue->getAtomsOnSidechain(targetAtomsOfSidechain);
		
	rg_dList<V::GeometryTier::Atom*> atomsOfBackbone;
	backbone.getBackboneAtoms(atomsOfBackbone);

	targetAtomsOfSidechain.reset4Loop();
	while (targetAtomsOfSidechain.setNext4Loop())
	{
		rg_REAL xVolume = 0.0;
        V::GeometryTier::Atom* currAtom = targetAtomsOfSidechain.getEntity();
		xVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfBackbone);
		sumOfPairwiseXVolume += xVolume;
	}

	return sumOfPairwiseXVolume;
}

// NOTE: moved to header file
// because this makes LNK1120 
//inline Residue* Rotamer::getResidue()
//{
//	return m_residue;
//}

//inline rg_REAL* Rotamer::getDihedralAngles()
//{
//	return m_dihedralAngles;
//}

rg_INT  Rotamer::getDihedralAngles(rg_REAL*& dihedralAngles)
{
	dihedralAngles = new rg_REAL[ m_numDihedralAngles ];
	rg_INDEX i;
	for (i = 0;i < m_numDihedralAngles;i++)
	{
		dihedralAngles[ i ] = m_dihedralAngles[ i ];
	}

	return m_numDihedralAngles;
}

rg_INT Rotamer::getAtomsOnSidechain(rg_dList<V::GeometryTier::Atom*>& atomsOnSidechain )
{
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	//return FunctionsForRotamerLibrary::getAtomsOnSidechain(m_residue, atomsOnSidechain);
	return m_residue->getAtomsOnSidechain(atomsOnSidechain);
}

rg_INT Rotamer::getAtomsOnSidechain(V::GeometryTier::Atom**& atomsOnSidechain )
{
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	//return FunctionsForRotamerLibrary::getAtomsOnSidechain(m_residue, atomsOnSidechain);
	return m_residue->getAtomsOnSidechain( atomsOnSidechain );
}

rg_INT Rotamer::getAtomsOfResidue( rg_dList<V::GeometryTier::Atom*>& atomsOfResidue )
{
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	//return FunctionsForRotamerLibrary::getAtomsOfResidue(m_residue, atomsOfResidue);
	m_residue->getAtoms(atomsOfResidue);
	return atomsOfResidue.getSize();
}

void Rotamer::setDihedralAngles(const rg_REAL* dihedralAngles)
{
	if(dihedralAngles == rg_NULL)
		return;

	destroyDihedralAngles();

	m_numDihedralAngles 
		= FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(m_residue->getResidueCode());

	if (m_numDihedralAngles > 0)
	{
		m_dihedralAngles = new rg_REAL[m_numDihedralAngles];

		rg_INDEX i;
		for (i = 0;i < m_numDihedralAngles;i++)
		{
			m_dihedralAngles[ i ] = dihedralAngles[ i ];
		}
	}
}


void  Rotamer::setResidue(V::GeometryTier::Residue* residue)
{
	m_residue    = residue;

	destroyDihedralAngles();

	m_numDihedralAngles
		= FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(m_residue->getResidueCode());
}


void Rotamer::setResidueWithDihedralAngles(V::GeometryTier::Residue* residue)
{
	m_residue    = residue;
	setDihedralAngles( m_residue );
}


void    Rotamer::updateSidechainAtomCoordinates()
{
	//if(m_residue->hasSideChainDihedralAngle())
    if(m_numDihedralAngles > 0)
	{
		FunctionsForRotamerLibrary::transformAtomCentersOfSidechainUsingDihedralAngle(m_residue, m_dihedralAngles, m_numDihedralAngles);

		// if the side chain comes from one of the following 5 amino acids
		// then glue the remaining rigid part into the corresponding side chain
		// HIS_AMINO_RESIDUE, PHE_AMINO_RESIDUE, TRP_AMINO_RESIDUE, TYR_AMINO_RESIDUE, ARG_AMINO_RESIDUE	
		
        //FunctionsForRotamerLibrary::glueFixedStructureIntoResidue(m_residue);
        FunctionsForRotamerLibrary::glueFixedStructureIntoResidue2(m_residue);
        // fix bug for 5 types of residue (ARG, HIS, TRP, TYR, PHE)

		return;
	}
	else
	{
		// do nothing
		// because this side chain does not have any dihedral angle
		return;
	}
}


void Rotamer::updateSidechainAtomCoordinatesByApplyingRotLibID(const rg_INDEX& rotLibID)
{
	if (m_numDihedralAngles > 0)
	{
#ifdef _DEBUG
if(rotLibID == 477005)
	int here = 1;
#endif
        if(rotLibID != FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX && 
		   rotLibID != USE_CURRENT_CONFORMATION &&
		   rotLibID != UNKNOWN_ROT_INDEX)
        {
		    rg_REAL* dihedralAngles = rg_NULL;
		    dihedralAngles = new rg_REAL[m_numDihedralAngles];
		    FunctionsForRotamerLibrary::getDihedralAngles(rotLibID, dihedralAngles);
		    updateSidechainAtomCoordinatesByApplyingDihedralAngles( dihedralAngles );
		    if(dihedralAngles != rg_NULL)
			    delete [] dihedralAngles;
        }
	}
}

void Rotamer::updateSidechainAtomCoordinatesByApplyingDihedralAngles(const rg_REAL* dihedralAngles)
{
	setDihedralAngles( dihedralAngles );
	updateSidechainAtomCoordinates();
}

//void Rotamer::perturbSidechainAtomsByRotation(const Vector3D& rotAxis,
//                                              const rg_REAL& angle)
//{
//    rg_TMatrix3D rotMat;
//    rotMat.rotateArbitraryAxis(rotAxis, rg_PI * angle / 180.);
//
//    rg_dList<Atom*> sideChainAtoms;
//    //getAtomsOnSidechain(targetAtomsOfSidechain);
//    m_residue->getAtomsOnSidechain(sideChainAtoms);
//
//    sideChainAtoms.reset4Loop();
//    while (sideChainAtoms.setNext4Loop())
//    {
//        Atom* currAtom = sideChainAtoms.getEntity();
//        rg_Point3D currCenter = currAtom->getpAtomBall()->getCenter();
//        currAtom->getpAtomBall()->setCenter( rotMat * currCenter );
//    }
//}

Rotamer& Rotamer::operator=(const Rotamer& rotamer)
{
	if(this == &rotamer)
		return *this;

	m_residue = rotamer.m_residue;
	setDihedralAngles(rotamer.m_dihedralAngles);

	return *this;
}

void Rotamer::destroyDihedralAngles()
{
	if(m_numDihedralAngles > 0 && m_dihedralAngles != rg_NULL)
	{
		delete [] m_dihedralAngles;
		m_dihedralAngles = rg_NULL;
		m_numDihedralAngles = 0;
	}
}