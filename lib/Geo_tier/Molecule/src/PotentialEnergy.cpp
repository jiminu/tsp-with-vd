#include "PotentialEnergy.h"
using namespace V::GeometryTier;



PotentialEnergy::PotentialEnergy()
: m_flagsOfAppliedEnergyTerms(0)
{   
}



PotentialEnergy::PotentialEnergy( const PotentialEnergy& potentialEnergy )
: m_flagsOfAppliedEnergyTerms(potentialEnergy.m_flagsOfAppliedEnergyTerms)
{
}



PotentialEnergy::~PotentialEnergy()
{
}



rg_INT PotentialEnergy::getFlagsOfAppliedEnergyTerms() const
{
    return m_flagsOfAppliedEnergyTerms;    
}



void PotentialEnergy::setFlagsOfAppliedEnergyTerms( const rg_INT& flagsOfEnergyTerms )
{
    m_flagsOfAppliedEnergyTerms = flagsOfEnergyTerms;    
}



rg_REAL PotentialEnergy::computeAndGetPotentialEnergyForTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule )
{
	rg_REAL energy = .0;

	if( m_flagsOfAppliedEnergyTerms&VDW_TERM_ON ) {
        energy  += computeVanDerWaalsBetweenTwoMolecules( fixedMolecule, movingMolecule );
    }
	if( m_flagsOfAppliedEnergyTerms&ELEC_TERM_ON ) {
        energy  += computeElectrostaticBetweenTwoMolecules( fixedMolecule, movingMolecule );
    }
	if( m_flagsOfAppliedEnergyTerms&HBOND_TERM_ON ) {
        energy  += computeHydrogenBondBetweenTwoMolecules( fixedMolecule, movingMolecule ); 
    }

//     rg_dList<Atom>* listOfAtomsInFixedMol  = fixedMolecule->getAtoms();
//     rg_dList<Atom>* listOfAtomsInMovingMol = movingMolecule->getAtoms();
//     
//     listOfAtomsInFixedMol->reset4Loop();
//     while( listOfAtomsInFixedMol->setNext4Loop() ) {
//         Atom*      currAtomInFixedMol = listOfAtomsInFixedMol->getpEntity();
//         rg_Point3D centerOfFixedAtom = currAtomInFixedMol->getAtomBall().getCenter();
//         rg_INT     atomTypeInAmberForFixedAtom = currAtomInFixedMol->getChemicalProperties().getAtomTypeInAmber();
//         rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getCharge();
// 
//         listOfAtomsInMovingMol->reset4Loop();
//         while( listOfAtomsInMovingMol->setNext4Loop() ) {
//             Atom*      currAtomInMovingMol = listOfAtomsInMovingMol->getpEntity();
//             rg_Point3D centerOfMovingAtom = currAtomInMovingMol ->getAtomBall().getCenter();
//             rg_INT     atomTypeInAmberForMovingAtom = currAtomInMovingMol->getChemicalProperties().getAtomTypeInAmber();
//             rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getCharge();
// 
//             rg_REAL    denom = ( centerOfMovingAtom - centerOfFixedAtom ).squaredMagnitude();
// 
//             if( m_flagsOfAppliedEnergyTerms&VDW_TERM_ON ) {
//                 energy += (   L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][0] / pow(denom, 6.0)
//                             - L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][1] / pow(denom, 3.0) );
//             }
//             if( m_flagsOfAppliedEnergyTerms&ELEC_TERM_ON ) {
//                 energy += 0.5 * chargeForFixedAtom * chargeForMovingAtom / denom;
//             }
//             if( m_flagsOfAppliedEnergyTerms&HBOND_TERM_ON ) {
//                 //
//             }
//         }
//     }

    return energy;
}



void PotentialEnergy::computeAndGetPotentialEnergyForTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule, rg_REAL* energy )
{
	rg_REAL VDWEnergy   = 0.0;
    rg_REAL ElecEnergy  = 0.0;
    rg_REAL HBondEnergy = 0.0;

	if( m_flagsOfAppliedEnergyTerms & VDW_TERM_ON ) {
        VDWEnergy  += computeVanDerWaalsBetweenTwoMolecules( fixedMolecule, movingMolecule );
    }
	if( m_flagsOfAppliedEnergyTerms & ELEC_TERM_ON ) {
        ElecEnergy  += computeElectrostaticBetweenTwoMolecules( fixedMolecule, movingMolecule );
    }
	if( m_flagsOfAppliedEnergyTerms & HBOND_TERM_ON ) {
        HBondEnergy  += computeHydrogenBondBetweenTwoMolecules( fixedMolecule, movingMolecule ); 
    }

    energy[ID_VDW_TERM]   = VDWEnergy;
    energy[ID_ELEC_TERM]  = ElecEnergy;
    energy[ID_HBOND_TERM] = HBondEnergy;
}


// FUNCTION IN USE...
void PotentialEnergy::computeAndGetPotentialEnergyForTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol, rg_REAL* energy )
{
	if( m_flagsOfAppliedEnergyTerms & VDW_TERM_ON ) {
        energy[ID_VDW_TERM] += computeVanDerWaalsBetweenTwoMolecules( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol );
    }
	if( m_flagsOfAppliedEnergyTerms & ELEC_TERM_ON ) {
        energy[ID_ELEC_TERM] += computeElectrostaticBetweenTwoMolecules( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol );
    }
	if( m_flagsOfAppliedEnergyTerms & HBOND_TERM_ON ) {
        energy[ID_HBOND_TERM] += computeHydrogenBondBetweenTwoMolecules( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol );
    }
}


rg_REAL PotentialEnergy::computeAndGetPotentialEnergyForTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol)
{
	rg_REAL energy = .0;

	if( m_flagsOfAppliedEnergyTerms & VDW_TERM_ON ) {
		energy += computeVanDerWaalsBetweenTwoMolecules( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol );
	}
	if( m_flagsOfAppliedEnergyTerms & ELEC_TERM_ON ) {
		energy += computeElectrostaticBetweenTwoMolecules( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol );
	}
	if( m_flagsOfAppliedEnergyTerms & HBOND_TERM_ON ) {
		energy += computeHydrogenBondBetweenTwoMolecules( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol );
	}

	return energy;
}

rg_REAL PotentialEnergy::computeVanDerWaalsBetweenTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule )
{
    rg_REAL energy = 0.0;

    rg_dList<Atom>* listOfAtomsInFixedMol  = fixedMolecule->getAtoms();
    rg_dList<Atom>* listOfAtomsInMovingMol = movingMolecule->getAtoms();
    

    Atom*      currAtomInFixedMol          = rg_NULL;
    rg_Point3D centerOfFixedAtom;
    rg_INT     atomTypeInAmberForFixedAtom = -1;

    Atom*      currAtomInMovingMol          = rg_NULL;
    rg_Point3D centerOfMovingAtom;
    rg_INT     atomTypeInAmberForMovingAtom = -1;

    listOfAtomsInFixedMol->reset4Loop();
    while( listOfAtomsInFixedMol->setNext4Loop() ) {
        currAtomInFixedMol = listOfAtomsInFixedMol->getpEntity();
        centerOfFixedAtom  = currAtomInFixedMol->getAtomBall().getCenter();
        atomTypeInAmberForFixedAtom = currAtomInFixedMol->getChemicalProperties().getAtomTypeInAmber();

        listOfAtomsInMovingMol->reset4Loop();
        while( listOfAtomsInMovingMol->setNext4Loop() ) {
            currAtomInMovingMol = listOfAtomsInMovingMol->getpEntity();
            centerOfMovingAtom  = currAtomInMovingMol ->getAtomBall().getCenter();
            atomTypeInAmberForMovingAtom = currAtomInMovingMol->getChemicalProperties().getAtomTypeInAmber();

            rg_REAL squaredDistance      = ( centerOfMovingAtom - centerOfFixedAtom ).squaredMagnitude();
            rg_REAL sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

            energy += (   ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][0] 
                            / (sixthPoweredDistance*sixthPoweredDistance) )
                        - ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][1] 
                            / sixthPoweredDistance ) );
        }
    }
    
    return energy;
}


// FUNCTION IN USE...
rg_REAL PotentialEnergy::computeVanDerWaalsBetweenTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol )
{
    rg_REAL energy = 0.0;

    rg_Point3D centerOfFixedAtom;
    rg_INT     atomTypeInAmberForFixedAtom  = -1;

    rg_Point3D centerOfMovingAtom;
    rg_INT     atomTypeInAmberForMovingAtom = -1;

    rg_REAL squaredDistance      = 0.0;
    rg_REAL sixthPoweredDistance = 0.0;

    for( rg_INT i_atomsInFixedMol=0; i_atomsInFixedMol<numOfAtomsInFixedMol; i_atomsInFixedMol++ ) {
        centerOfFixedAtom = atomsInFixedMolecule[i_atomsInFixedMol]->getpAtomBall()->getCenter();
        atomTypeInAmberForFixedAtom = (rg_INT)atomsInFixedMolecule[i_atomsInFixedMol]->getpChemicalProperties()->getAtomTypeInAmber();

//         if( atomsInFixedMolecule[i_atomsInFixedMol]->getID() == 58 )
//             int aaa =0;

        for( rg_INT j_atomsInMovingMol=0; j_atomsInMovingMol<numOfAtomsInMovingMol; j_atomsInMovingMol++ ) {
            centerOfMovingAtom = atomsInMovingMolecule[j_atomsInMovingMol].getpAtomBall()->getCenter();
            atomTypeInAmberForMovingAtom = (rg_INT)atomsInMovingMolecule[j_atomsInMovingMol].getpChemicalProperties()->getAtomTypeInAmber();

            squaredDistance = ( centerOfMovingAtom - centerOfFixedAtom ).squaredMagnitude();
            sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

            energy += (   ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][0] 
                            / (sixthPoweredDistance*sixthPoweredDistance) )
                        - ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][1] 
                            / sixthPoweredDistance ) );
        }
    }

    return energy;    
}



rg_REAL PotentialEnergy::computeElectrostaticBetweenTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule )
{
    rg_REAL energy = .0;

    rg_dList<Atom>* listOfAtomsInFixedMol  = fixedMolecule->getAtoms();
    rg_dList<Atom>* listOfAtomsInMovingMol = movingMolecule->getAtoms();
    
    listOfAtomsInFixedMol->reset4Loop();
    while( listOfAtomsInFixedMol->setNext4Loop() ) {
        Atom*      currAtomInFixedMol = listOfAtomsInFixedMol->getpEntity();
        rg_Point3D centerOfFixedAtom  = currAtomInFixedMol->getAtomBall().getCenter();
        ////rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getCharge();
        //// jhryu
        //rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getChargeInAmber();
        rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getChargeInAmber();

        listOfAtomsInMovingMol->reset4Loop();
        while( listOfAtomsInMovingMol->setNext4Loop() ) {
            Atom*      currAtomInMovingMol = listOfAtomsInMovingMol->getpEntity();
            rg_Point3D centerOfMovingAtom  = currAtomInMovingMol ->getAtomBall().getCenter();
            ////rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getCharge();
            ////jhryu
            //rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getChargeInAmber();
            rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getCharge();

            energy += (chargeForFixedAtom * chargeForMovingAtom) / ( (centerOfMovingAtom - centerOfFixedAtom).squaredMagnitude() );
        }
    }
    
    return energy;
}



rg_REAL PotentialEnergy::computeElectrostaticBetweenTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol )
{
    rg_REAL energy = .0;

    rg_Point3D centerOfFixedAtom;
    rg_REAL    chargeForFixedAtom = 0.0;

    rg_Point3D centerOfMovingAtom;
    rg_REAL    chargeForMovingAtom = 0.0;

    for( rg_INT i_atomsInFixedMol=0; i_atomsInFixedMol<numOfAtomsInFixedMol; i_atomsInFixedMol++ ) {
        centerOfFixedAtom = atomsInFixedMolecule[i_atomsInFixedMol]->getpAtomBall()->getCenter();
        ////chargeForFixedAtom = atomsInFixedMolecule[i_atomsInFixedMol]->getpChemicalProperties()->getCharge();
        //// jhryu
        //chargeForFixedAtom = atomsInFixedMolecule[i_atomsInFixedMol]->getpChemicalProperties()->getChargeInAmber();

        //  ycho 2017.12.29
        chargeForFixedAtom = atomsInFixedMolecule[i_atomsInFixedMol]->getpChemicalProperties()->getChargeInAmber();
        for( rg_INT j_atomsInMovingMol=0; j_atomsInMovingMol<numOfAtomsInMovingMol; j_atomsInMovingMol++ ) {
            centerOfMovingAtom = atomsInMovingMolecule[j_atomsInMovingMol].getpAtomBall()->getCenter();
            ////jhryu
            ////chargeForMovingAtom = atomsInMovingMolecule[j_atomsInMovingMol].getpChemicalProperties()->getCharge();
            //chargeForMovingAtom = atomsInMovingMolecule[j_atomsInMovingMol].getpChemicalProperties()->getChargeInAmber();

            //  ycho 2017.12.29
            chargeForMovingAtom = atomsInMovingMolecule[j_atomsInMovingMol].getpChemicalProperties()->getCharge();
            energy += chargeForFixedAtom * chargeForMovingAtom / (centerOfMovingAtom-centerOfFixedAtom).squaredMagnitude();
        }
    }
    
    return energy;    
}



rg_REAL PotentialEnergy::computeHydrogenBondBetweenTwoMolecules( Molecule* fixedMolecule, Molecule* movingMolecule )
{
    rg_REAL energy = .0;

    rg_INT numOfAtomsInFixedMolecule  = fixedMolecule->getAtoms()->getSize();
    rg_INT numOfAtomsInMovingMolecule = movingMolecule->getAtoms()->getSize();

    rg_dList<HydrogenBond> listOfHBondCandidates;
    getListOfHydrogenBondCandidates( fixedMolecule, movingMolecule, listOfHBondCandidates );

    rg_dList<HydrogenBond> listOfHydrogenBond;
    getListOfHydrogenBond( numOfAtomsInFixedMolecule, numOfAtomsInMovingMolecule, &listOfHBondCandidates, listOfHydrogenBond );



    listOfHBondCandidates.removeAll();
    getListOfHydrogenBondCandidates( movingMolecule, fixedMolecule, listOfHBondCandidates );

    getListOfHydrogenBond( numOfAtomsInMovingMolecule, numOfAtomsInFixedMolecule, &listOfHBondCandidates, listOfHydrogenBond );


    listOfHydrogenBond.reset4Loop();
    while( listOfHydrogenBond.setNext4Loop() ) {
        HydrogenBond* currHBond = listOfHydrogenBond.getpEntity();
        energy += currHBond->getHydrogenBondEnergy();
    }

    return energy;
}



rg_REAL PotentialEnergy::computeHydrogenBondBetweenTwoMolecules( Atom** atomsInFixedMolecule, const rg_INT& numOfAtomsInFixedMol, Atom* atomsInMovingMolecule, const rg_INT& numOfAtomsInMovingMol )
{
    rg_REAL energy = .0;

    rg_dList<HydrogenBond> listOfHBondCandidates;
    getListOfHydrogenBondCandidates( atomsInFixedMolecule, numOfAtomsInFixedMol, atomsInMovingMolecule, numOfAtomsInMovingMol, listOfHBondCandidates );

    rg_dList<HydrogenBond> listOfHydrogenBond;
    getListOfHydrogenBond( numOfAtomsInFixedMol, numOfAtomsInMovingMol, &listOfHBondCandidates, listOfHydrogenBond );

    listOfHBondCandidates.removeAll();
    getListOfHydrogenBondCandidates( atomsInMovingMolecule, numOfAtomsInMovingMol, atomsInFixedMolecule, numOfAtomsInFixedMol, listOfHBondCandidates );

    getListOfHydrogenBond( numOfAtomsInMovingMol, numOfAtomsInFixedMol, &listOfHBondCandidates, listOfHydrogenBond );

    listOfHydrogenBond.reset4Loop();
    while( listOfHydrogenBond.setNext4Loop() ) {
        HydrogenBond* currHBond = listOfHydrogenBond.getpEntity();
        energy += currHBond->getHydrogenBondEnergy();
    }

    return energy;
}


rg_REAL PotentialEnergy::computeAndGetPotentialEnergyForTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule)
{
    rg_REAL energy = .0;

    if (m_flagsOfAppliedEnergyTerms&VDW_TERM_ON) {
        energy += computeVanDerWaalsBetweenTwoMolecules_IgnoredHydrogen(fixedMolecule, movingMolecule);
    }
    if (m_flagsOfAppliedEnergyTerms&ELEC_TERM_ON) {
        energy += computeElectrostaticBetweenTwoMolecules_IgnoredHydrogen(fixedMolecule, movingMolecule);
    }
    if (m_flagsOfAppliedEnergyTerms&HBOND_TERM_ON) {
        energy += 0.0;
    }

    return energy;
}



void PotentialEnergy::computeAndGetPotentialEnergyForTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule, rg_REAL* energy)
{
    rg_REAL VDWEnergy = 0.0;
    rg_REAL ElecEnergy = 0.0;
    rg_REAL HBondEnergy = 0.0;

    if (m_flagsOfAppliedEnergyTerms & VDW_TERM_ON) {
        VDWEnergy += computeVanDerWaalsBetweenTwoMolecules_IgnoredHydrogen(fixedMolecule, movingMolecule);
    }
    if (m_flagsOfAppliedEnergyTerms & ELEC_TERM_ON) {
        ElecEnergy += computeElectrostaticBetweenTwoMolecules_IgnoredHydrogen(fixedMolecule, movingMolecule);
    }
    if (m_flagsOfAppliedEnergyTerms & HBOND_TERM_ON) {
        HBondEnergy += 0.0;
    }

    energy[ID_VDW_TERM] = VDWEnergy;
    energy[ID_ELEC_TERM] = ElecEnergy;
    energy[ID_HBOND_TERM] = HBondEnergy;
}



rg_REAL PotentialEnergy::computeVanDerWaalsBetweenTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule)
{
    rg_REAL energy = 0.0;

    rg_dList<Atom>* listOfAtomsInFixedMol = fixedMolecule->getAtoms();
    rg_dList<Atom>* listOfAtomsInMovingMol = movingMolecule->getAtoms();


    Atom*      currAtomInFixedMol = rg_NULL;
    rg_Point3D centerOfFixedAtom;
    rg_INT     atomTypeInAmberForFixedAtom = -1;

    Atom*      currAtomInMovingMol = rg_NULL;
    rg_Point3D centerOfMovingAtom;
    rg_INT     atomTypeInAmberForMovingAtom = -1;

    listOfAtomsInFixedMol->reset4Loop();
    while (listOfAtomsInFixedMol->setNext4Loop()) {
        currAtomInFixedMol = listOfAtomsInFixedMol->getpEntity();
        if (currAtomInFixedMol->getAtomCode() == H_ATOM) {
            continue;
        }

        centerOfFixedAtom = currAtomInFixedMol->getAtomBall().getCenter();
        atomTypeInAmberForFixedAtom = currAtomInFixedMol->getChemicalProperties().getAtomTypeInAmber();

        listOfAtomsInMovingMol->reset4Loop();
        while (listOfAtomsInMovingMol->setNext4Loop()) {
            currAtomInMovingMol = listOfAtomsInMovingMol->getpEntity();
            if (currAtomInMovingMol->getAtomCode() == H_ATOM) {
                continue;
            }

            centerOfMovingAtom = currAtomInMovingMol->getAtomBall().getCenter();
            atomTypeInAmberForMovingAtom = currAtomInMovingMol->getChemicalProperties().getAtomTypeInAmber();

            rg_REAL squaredDistance = (centerOfMovingAtom - centerOfFixedAtom).squaredMagnitude();
            rg_REAL sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

            energy += ((L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][0]
                / (sixthPoweredDistance*sixthPoweredDistance))
                - (L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForFixedAtom][atomTypeInAmberForMovingAtom][1]
                    / sixthPoweredDistance));
        }
    }

    return energy;
}



rg_REAL PotentialEnergy::computeElectrostaticBetweenTwoMolecules_IgnoredHydrogen(Molecule* fixedMolecule, Molecule* movingMolecule)
{
    rg_REAL energy = .0;

    rg_dList<Atom>* listOfAtomsInFixedMol = fixedMolecule->getAtoms();
    rg_dList<Atom>* listOfAtomsInMovingMol = movingMolecule->getAtoms();

    listOfAtomsInFixedMol->reset4Loop();
    while (listOfAtomsInFixedMol->setNext4Loop()) {
        Atom*      currAtomInFixedMol = listOfAtomsInFixedMol->getpEntity();
        if (currAtomInFixedMol->getAtomCode() == H_ATOM) {
            continue;
        }

        rg_Point3D centerOfFixedAtom = currAtomInFixedMol->getAtomBall().getCenter();
        ////rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getCharge();
        //// jhryu
        //rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getChargeInAmber();
        rg_REAL    chargeForFixedAtom = currAtomInFixedMol->getChemicalProperties().getChargeInAmber();

        listOfAtomsInMovingMol->reset4Loop();
        while (listOfAtomsInMovingMol->setNext4Loop()) {
            Atom*      currAtomInMovingMol = listOfAtomsInMovingMol->getpEntity();
            if (currAtomInMovingMol->getAtomCode() == H_ATOM) {
                continue;
            }

            rg_Point3D centerOfMovingAtom = currAtomInMovingMol->getAtomBall().getCenter();
            ////rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getCharge();
            ////jhryu
            //rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getChargeInAmber();
            rg_REAL    chargeForMovingAtom = currAtomInMovingMol->getChemicalProperties().getCharge();

            energy += (chargeForFixedAtom * chargeForMovingAtom) / ((centerOfMovingAtom - centerOfFixedAtom).squaredMagnitude());
        }
    }

    return energy;
}



void PotentialEnergy::getListOfHydrogenBondCandidates( Molecule* moleculeForHydrogens, Molecule* moleculeForAcceptors, rg_dList<HydrogenBond>& listOfHBondCandidates )
{
    rg_dList<Atom*> listOfHydrogenAtoms;
    rg_dList<Atom*> listOfHAcceptorAtoms;

    moleculeForHydrogens->getListOfHAtomsFromHDonor( &listOfHydrogenAtoms );
    moleculeForAcceptors->getListOfHAcceptorAtoms( &listOfHAcceptorAtoms );

    listOfHydrogenAtoms.reset4Loop();
    while( listOfHydrogenAtoms.setNext4Loop() ) {
        Atom* currHydrogen = listOfHydrogenAtoms.getEntity();
        listOfHAcceptorAtoms.reset4Loop();
        while( listOfHAcceptorAtoms.setNext4Loop() ) {
            Atom* currAcceptor = listOfHAcceptorAtoms.getEntity();
            HydrogenBond hBondCandidate( currHydrogen, currAcceptor );
            addHBondCandidateToListWithDistanceOrder( listOfHBondCandidates, hBondCandidate );
        }
    }
}



void PotentialEnergy::getListOfHydrogenBondCandidates( Atom** atomsInMolForHydrogens, const rg_INT& numOfAtomsInMolForHydrogens, Atom* atomsInMolForAcceptors, const rg_INT& numOfAtomsInMolForAcceptors, rg_dList<HydrogenBond>& listOfHBondCandidates )
{
    /////////////////////////////////// getListOfHAtomsFromHDonor( &listOfHydrogenAtoms ) ///////////////////////////////////
    //                                                                                                                     //
    rg_dList<Atom*> listOfHydrogenFromHDonor;

    Atom* firstBondedAtom  = rg_NULL;
    Atom* secondBondedAtom = rg_NULL;

    rg_INT i_atoms=0;
    for( i_atoms=0; i_atoms<numOfAtomsInMolForHydrogens; i_atoms++ ) {
        if( atomsInMolForHydrogens[i_atoms]->getAtomCode() == H_ATOM && atomsInMolForHydrogens[i_atoms]->getListChemicalBond()->getSize() != 0 ) {
            atomsInMolForHydrogens[i_atoms]->getListChemicalBond()->getFirstEntity()->getAtoms( firstBondedAtom, secondBondedAtom );

            if( ( firstBondedAtom->getAtomCode()  == O_ATOM || firstBondedAtom->getAtomCode()  == N_ATOM ) ||
                ( secondBondedAtom->getAtomCode() == O_ATOM || secondBondedAtom->getAtomCode() == N_ATOM )  ) {
                listOfHydrogenFromHDonor.addTail( atomsInMolForHydrogens[i_atoms] );
            }
        }
    }
    //                                                                                                                     //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////// getListOfHAcceptorAtoms( &listOfHydrogenAtoms ) ////////////////////////////////////
    //                                                                                                                     //
    rg_dList<Atom*> listOfHAcceptorAtoms;

    for( i_atoms=0; i_atoms<numOfAtomsInMolForAcceptors; i_atoms++ ) {
        if( atomsInMolForAcceptors[i_atoms].getpChemicalProperties()->getNumOfAcceptableHydrogen() >= 1 )
            listOfHAcceptorAtoms.addTail( &atomsInMolForAcceptors[i_atoms] );
    }
    //                                                                                                                     //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    listOfHydrogenFromHDonor.reset4Loop();
    while( listOfHydrogenFromHDonor.setNext4Loop() ) {
        Atom* currHydrogen = listOfHydrogenFromHDonor.getEntity();
        listOfHAcceptorAtoms.reset4Loop();
        while( listOfHAcceptorAtoms.setNext4Loop() ) {
            Atom* currAcceptor = listOfHAcceptorAtoms.getEntity();
            HydrogenBond hBondCandidate( currHydrogen, currAcceptor );
            addHBondCandidateToListWithDistanceOrder( listOfHBondCandidates, hBondCandidate );
        }
    }

}



void PotentialEnergy::getListOfHydrogenBondCandidates( Atom* atomsInMolForHydrogens, const rg_INT& numOfAtomsInMolForHydrogens, Atom** atomsInMolForAcceptors, const rg_INT& numOfAtomsInMolForAcceptors, rg_dList<HydrogenBond>& listOfHBondCandidates )
{
    /////////////////////////////////// getListOfHAtomsFromHDonor( &listOfHydrogenAtoms ) ///////////////////////////////////
    //                                                                                                                     //
    rg_dList<Atom*> listOfHydrogenFromHDonor;

    Atom* firstBondedAtom  = rg_NULL;
    Atom* secondBondedAtom = rg_NULL;

    rg_INT i_atoms=0;
    for( i_atoms=0; i_atoms<numOfAtomsInMolForHydrogens; i_atoms++ ) {
        if( atomsInMolForHydrogens[i_atoms].getAtomCode() == H_ATOM && atomsInMolForHydrogens[i_atoms].getListChemicalBond()->getSize() != 0 ) {
            atomsInMolForHydrogens[i_atoms].getListChemicalBond()->getFirstEntity()->getAtoms( firstBondedAtom, secondBondedAtom );

            if( ( firstBondedAtom->getAtomCode()  == O_ATOM || firstBondedAtom->getAtomCode()  == N_ATOM ) ||
                ( secondBondedAtom->getAtomCode() == O_ATOM || secondBondedAtom->getAtomCode() == N_ATOM )  )
                listOfHydrogenFromHDonor.addTail( &atomsInMolForHydrogens[i_atoms] );
        }
    }
    //                                                                                                                     //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////// getListOfHAcceptorAtoms( &listOfHydrogenAtoms ) ////////////////////////////////////
    //                                                                                                                     //
    rg_dList<Atom*> listOfHAcceptorAtoms;

    for( i_atoms=0; i_atoms<numOfAtomsInMolForAcceptors; i_atoms++ ) {
        if( atomsInMolForAcceptors[i_atoms]->getpChemicalProperties()->getNumOfAcceptableHydrogen() >= 1 )
            listOfHAcceptorAtoms.addTail( atomsInMolForAcceptors[i_atoms] );
    }
    //                                                                                                                     //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    listOfHydrogenFromHDonor.reset4Loop();
    while( listOfHydrogenFromHDonor.setNext4Loop() ) {
        Atom* currHydrogen = listOfHydrogenFromHDonor.getEntity();
        listOfHAcceptorAtoms.reset4Loop();
        while( listOfHAcceptorAtoms.setNext4Loop() ) {
            Atom* currAcceptor = listOfHAcceptorAtoms.getEntity();
            HydrogenBond hBondCandidate( currHydrogen, currAcceptor );
            addHBondCandidateToListWithDistanceOrder( listOfHBondCandidates, hBondCandidate );
        }
    }

}



void PotentialEnergy::getListOfHydrogenBond( const rg_INT& numOfAtomsInMolForHydrogens, const rg_INT& numOfAtomsInMolForAcceptors, rg_dList<HydrogenBond>* listOfHBondCandidates, rg_dList<HydrogenBond>& listOfHydrogenBond )
{
    rg_INT  i=0;

    rg_FLAG* isHydrogenBonded          = new rg_FLAG[ numOfAtomsInMolForHydrogens ];
    for( i=0; i<numOfAtomsInMolForHydrogens; i++ )
        isHydrogenBonded[i] = rg_FALSE;
    
    rg_INT* numOfAcceptedHydrogensForAcceptors = new rg_INT[ numOfAtomsInMolForAcceptors ];
    for( i=0; i<numOfAtomsInMolForAcceptors; i++ )
        numOfAcceptedHydrogensForAcceptors[i] = 0;


    listOfHBondCandidates->reset4Loop();
    while ( listOfHBondCandidates->getSize() != 0 ) {
        listOfHBondCandidates->setNext4Loop();

        HydrogenBond* currHBondCandidate = listOfHBondCandidates->getpEntity();
        Atom*  currAcceptor              = currHBondCandidate->getAcceptor();
        Atom*  currHydrogen              = currHBondCandidate->getHydrogen();
        rg_INT numOfAcceptableHydrogen   = currAcceptor->getpChemicalProperties()->getNumOfAcceptableHydrogen();

        if( isHydrogenBonded[ currHydrogen->getID() ] == rg_FALSE &&
            numOfAcceptedHydrogensForAcceptors[ currAcceptor->getID() ] < numOfAcceptableHydrogen ) {
            listOfHydrogenBond.addTail( *currHBondCandidate );
            numOfAcceptedHydrogensForAcceptors[ currAcceptor->getID() ] = numOfAcceptedHydrogensForAcceptors[ currAcceptor->getID() ] + 1;
            isHydrogenBonded[ currHydrogen->getID() ] = rg_TRUE;
        }
        listOfHBondCandidates->killCurrent();
    }

    delete [] isHydrogenBonded;
    delete [] numOfAcceptedHydrogensForAcceptors;
}



void PotentialEnergy::addHBondCandidateToListWithDistanceOrder( rg_dList<HydrogenBond>& listHBondCandidates, HydrogenBond& hBondCandidate )
{
    const rg_REAL UPPER_BOUND_DIST_HYDROGEN_ACCEPTOR        = 2.5;
    const rg_REAL UPPER_BOUND_DIST_DONOR_ACCEPTOR           = 3.9;
    const rg_REAL LOWER_BOUND_ANGLE_DONOR_HYDROGEN_ACCEPTOR = 90;

    if( ( hBondCandidate.getDistanceBetweenHydrogenAndAcceptor() <= UPPER_BOUND_DIST_HYDROGEN_ACCEPTOR        ) &&
        ( hBondCandidate.getDistanceBetweenDonorAndAcceptor()    <= UPPER_BOUND_DIST_DONOR_ACCEPTOR           ) &&
        ( hBondCandidate.getAngleBetweenDonorHydrogenAcceptor()  >= LOWER_BOUND_ANGLE_DONOR_HYDROGEN_ACCEPTOR )  ) {

        if( listHBondCandidates.getSize() == 0 ) {
            listHBondCandidates.addTail( hBondCandidate );
            return;
        }
        else if( listHBondCandidates.getSize() == 1 ) {
            rg_REAL distForNewCandidate   = hBondCandidate.getDistanceBetweenHydrogenAndAcceptor();
            rg_REAL distForFirstCandidate = listHBondCandidates.getFirstEntity().getDistanceBetweenHydrogenAndAcceptor();

            if( distForNewCandidate < distForFirstCandidate )
                listHBondCandidates.addHead( hBondCandidate );
            else
                listHBondCandidates.addTail( hBondCandidate );
            return;
        }
        else {
            rg_REAL distForNewCandidate   = hBondCandidate.getDistanceBetweenHydrogenAndAcceptor();
            rg_REAL distForFirstCandidate = listHBondCandidates.getFirstEntity().getDistanceBetweenHydrogenAndAcceptor();

            if( distForNewCandidate < distForFirstCandidate ) {
                listHBondCandidates.addHead( hBondCandidate );
                return;
            }

            rg_FLAG isCandidateAdded = rg_FALSE;
            listHBondCandidates.reset4Loop();
            while( listHBondCandidates.setNext4Loop() ) {
                rg_REAL distForPrevCandidate = listHBondCandidates.getPrevEntity().getDistanceBetweenHydrogenAndAcceptor();
                rg_REAL distForCurrCandidate = listHBondCandidates.getEntity().getDistanceBetweenHydrogenAndAcceptor();            
                if( ( distForNewCandidate >= distForPrevCandidate ) && ( distForNewCandidate < distForCurrCandidate ) ) {
                    listHBondCandidates.setCurrentPrev();
                    listHBondCandidates.add( hBondCandidate );
                    isCandidateAdded = rg_TRUE;
                    break;
                }
            }
            if( isCandidateAdded == rg_FALSE )
                listHBondCandidates.addTail( hBondCandidate );
        }
    }
    else
    {}   
}



PotentialEnergy& PotentialEnergy::operator=( const PotentialEnergy& potentialEnergy )
{
    if( this == &potentialEnergy )
        return *this;

    m_flagsOfAppliedEnergyTerms       = potentialEnergy.m_flagsOfAppliedEnergyTerms;

    return *this;    
}
