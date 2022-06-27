#ifndef _POTENTIALENERGY_H
#define _POTENTIALENERGY_H

#include "rg_RelativeOp.h"
#include "ConstForPotentialEnergy.h"
#include "rg_Molecule.h"
#include "HydrogenBond.h"
#include "ForceField.h"



namespace V {
namespace GeometryTier {



class PotentialEnergy
{
private:
    rg_INT                  m_flagsOfAppliedEnergyTerms;

public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    PotentialEnergy();
    PotentialEnergy( const PotentialEnergy& potentialEnergy );
    ~PotentialEnergy();

    //  GET FUNCTION
    rg_INT  getFlagsOfAppliedEnergyTerms() const;

    //  SET FUNCTION
    void    setFlagsOfAppliedEnergyTerms( const rg_INT& flagsOfEnergyTerms );

    //  COMPUTATION FUNCTION
    
    rg_REAL	computeAndGetPotentialEnergyForTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule );	
    void    computeAndGetPotentialEnergyForTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule, rg_REAL* energy );
    void    computeAndGetPotentialEnergyForTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol, rg_REAL* energy );

	// jhryu
	rg_REAL computeAndGetPotentialEnergyForTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol);
    
    rg_REAL	computeVanDerWaalsBetweenTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule );
    rg_REAL	computeVanDerWaalsBetweenTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol );
    rg_REAL	computeElectrostaticBetweenTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule );
    rg_REAL	computeElectrostaticBetweenTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol );
    rg_REAL	computeHydrogenBondBetweenTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule );
    rg_REAL	computeHydrogenBondBetweenTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol );


    // ycho
    rg_REAL	computeAndGetPotentialEnergyForTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule);
    void    computeAndGetPotentialEnergyForTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule, rg_REAL* energy);

    rg_REAL	computeVanDerWaalsBetweenTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule);
    rg_REAL	computeElectrostaticBetweenTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule);


    void    getListOfHydrogenBondCandidates( Molecule* moleculeForHydrogens, Molecule* moleculeForAcceptors, rg_dList<HydrogenBond>& listOfHBondCandidates );
    void    getListOfHydrogenBondCandidates( Atom** atomsInMolForHydrogens,  const rg_INT& numOfAtomsInMolForHydrogens, 
                                             Atom*  atomsInMolForAcceptors,  const rg_INT& numOfAtomsInMolForAcceptors, rg_dList<HydrogenBond>& listOfHBondCandidates );
    void    getListOfHydrogenBondCandidates( Atom*  atomsInMolForHydrogens,  const rg_INT& numOfAtomsInMolForHydrogens, 
                                             Atom** atomsInMolForAcceptors,  const rg_INT& numOfAtomsInMolForAcceptors, rg_dList<HydrogenBond>& listOfHBondCandidates );


    void    getListOfHydrogenBond( const rg_INT& numOfAtomsInMolForHydrogens, const rg_INT& numOfAtomsInMolForAcceptors, rg_dList<HydrogenBond>* listOfHBondCandidates, rg_dList<HydrogenBond>& listOfHydrogenBond );


    void    addHBondCandidateToListWithDistanceOrder( rg_dList<HydrogenBond>& listHBondCandidates, HydrogenBond& hBondCandidate );

    //  OPERATOR OVERLOADING
    PotentialEnergy& operator =( const PotentialEnergy& potentialEnergy );
};



} // namespace GeometryTier
} // namespace V


#endif

