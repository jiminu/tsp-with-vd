#include "Residue.h"
#include "rg_Atom.h"
#include "Chain.h"
#include "ChemicalBond.h"
#include "rg_GeoFunc.h"
#include "ConstForRotamerLibrary.h"
#include "FunctionsForRotamerLibrary.h"
#include "FunctionsForMolecule.h"
//using namespace V::GeometryTier;

#include "BackboneDihedralAnglePair.h"
#include "SidechainDihedralAngleSet.h"

#include "ConstForMoleculeIOFunctions.h"


V::GeometryTier::Residue::Residue()
: m_ID(0)
, m_residueCode(UNK_RESIDUE)
, m_chain(rg_NULL)
, m_sequenceNumber(0)
{   
}



V::GeometryTier::Residue::Residue( const rg_INT& ID )
: m_ID(ID)
, m_residueCode(UNK_RESIDUE)
, m_chain(rg_NULL)
, m_sequenceNumber(0)
{
}



V::GeometryTier::Residue::Residue( const rg_INT& ID, const V::GeometryTier::ResidueCode& aResidueCode )
: m_ID(ID)
, m_residueCode(aResidueCode)
, m_chain(rg_NULL)
, m_sequenceNumber(0)
{
}



V::GeometryTier::Residue::Residue( const rg_INT& ID, const V::GeometryTier::ResidueCode& aResidueCode, const rg_INT& sequenceNumber )
: m_ID(ID)
, m_residueCode(aResidueCode)
, m_chain(rg_NULL)
, m_sequenceNumber(sequenceNumber)
{
}



V::GeometryTier::Residue::Residue( const V::GeometryTier::Residue& aResidue )
{
    m_ID                  = aResidue.m_ID;
    m_residueCode         = aResidue.m_residueCode;
    m_atoms               = aResidue.m_atoms;
    m_chain               = aResidue.m_chain;
    m_sequenceNumber      = aResidue.m_sequenceNumber;
    m_residueName         = aResidue.m_residueName;
}



V::GeometryTier::Residue::~Residue()
{
}


void    V::GeometryTier::Residue::getAtoms(rg_dList<V::GeometryTier::Atom*>& atomsOfResidue)
{
	m_atoms.reset4Loop();
	while (m_atoms.setNext4Loop())
	{
		Atom* currAtom = m_atoms.getEntity();
		atomsOfResidue.add(currAtom);
	}
}

rg_INT  V::GeometryTier::Residue::getAtoms(V::GeometryTier::Atom**& atomsOfResidue )
{
	rg_dList<V::GeometryTier::Atom*> residueAtomsInList;
	getAtoms(residueAtomsInList);
	rg_INT numAtoms = residueAtomsInList.getSize();

	atomsOfResidue = new V::GeometryTier::Atom* [numAtoms];

	rg_INDEX index = 0;
	residueAtomsInList.reset4Loop();
	while (residueAtomsInList.setNext4Loop())
	{
		atomsOfResidue[ index ] = residueAtomsInList.getEntity();
		index++;
	}

	return numAtoms;
}

V::GeometryTier::Atom*   V::GeometryTier::Residue::getAtom( const string& atomNameFromInputFile )
{
    V::GeometryTier::Atom* targetAtom = rg_NULL;

    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        if ( currAtom->getAtomNameFromInputFile().compare( atomNameFromInputFile ) == 0 ) {
            targetAtom = currAtom;
            break;
        }
    }
    return targetAtom;
}



string  V::GeometryTier::Residue::getResidueNameInPDBFormat() const
{
    string residueNameInPDBFormat;

    if ( isStandardResidue() ) {
        residueNameInPDBFormat = (string)V::GeometryTier::RESIDUE_FEATURES[m_residueCode].threeCodeName;
    }
    else {
        if (m_residueName.length() >= PDB_ATOM_RECORD_RESNAME_LENGTH) {
            residueNameInPDBFormat = m_residueName.substr(0, 3);
        }
        else {
            residueNameInPDBFormat = m_residueName;
            int numSpaces = PDB_ATOM_RECORD_RESNAME_LENGTH - m_residueName.length();
            for (int i = 0; i < numSpaces; ++i) {
                residueNameInPDBFormat += " ";
            }
        }
    }

    return residueNameInPDBFormat;
}



V::GeometryTier::Atom*   V::GeometryTier::Residue::getOxygenInCarboxylGroupOfAminoResidue() const
{
    V::GeometryTier::Atom* targetAtom = rg_NULL;

    targetAtom = getAtom( O_ATOM, UNK_REMOTE, UNK_BRANCH );

    if( targetAtom == rg_NULL )
        targetAtom = getAtom( O_ATOM, TERMINATE_REMOTE, FOURTH_BRANCH );


    return targetAtom;
}



V::GeometryTier::Atom*   V::GeometryTier::Residue::getAtom( const AtomCode& typeOfAtom, const RemoteIndicator& remote, const BranchDesignator& branch ) const
{
    V::GeometryTier::Atom* targetAtom = rg_NULL;
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        
        if( currAtom->getAtomCode() == typeOfAtom &&
            currAtom->getpChemicalProperties()->getRemoteIndicator()  == remote &&
            currAtom->getpChemicalProperties()->getBrangeDesignator() == branch   ) {
            targetAtom = currAtom;
            break;
        }
    }

    return targetAtom;
}



V::GeometryTier::Atom*   V::GeometryTier::Residue::getAtom( const RemoteIndicator& remote, const BranchDesignator& branch ) const
{
    V::GeometryTier::Atom* targetAtom = rg_NULL;
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        
        if( currAtom->getpChemicalProperties()->getRemoteIndicator()  == remote &&
            currAtom->getpChemicalProperties()->getBrangeDesignator() == branch   ) 
		{
            targetAtom = currAtom;
            break;
        }
    }

        
    return targetAtom;
}

void    V::GeometryTier::Residue::getAtoms( const V::GeometryTier::RemoteIndicator& remote, rg_dList<V::GeometryTier::Atom*>& atomsWithThisRemote ) const
{	
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) 
	{
        V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        V::GeometryTier::RemoteIndicator currRemote = currAtom->getpChemicalProperties()->getRemoteIndicator();
        
        if( currRemote == remote ) 
		{
			atomsWithThisRemote.add(currAtom);
        }
		
		if( currRemote > remote )
		{
			return;
		}
    }
}



void    V::GeometryTier::Residue::getAtomsOnBackbone( rg_dList<V::GeometryTier::Atom*>& atomsOnBackbone )
{
	m_atoms.reset4Loop();
	while ( m_atoms.setNext4Loop() ) 
	{
        V::GeometryTier::Atom* currAtom = m_atoms.getEntity();

		if( currAtom->getpChemicalProperties()->isOnBackBone() == rg_TRUE ) {
			atomsOnBackbone.addTail( currAtom );
		}
	}
}



void    V::GeometryTier::Residue::getAtomsOnBackbone( rg_dList<V::GeometryTier::Atom*>* listOfAtomsOnBackbone )
{
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        
        if( currAtom->getpChemicalProperties()->isOnBackBone() == rg_TRUE ) {
            listOfAtomsOnBackbone->addTail( currAtom );
        }
    }
}



rg_INT  V::GeometryTier::Residue::getAtomsOnBackbone(V::GeometryTier::Atom**& atomsOnBackbone )
{	
	rg_dList<V::GeometryTier::Atom*> listOfAtomsOnBackbone;
	getAtomsOnBackbone( &listOfAtomsOnBackbone );
	rg_INT numBackboneAtoms = listOfAtomsOnBackbone.getSize();

	atomsOnBackbone = new V::GeometryTier::Atom*[numBackboneAtoms];

	rg_INDEX i = 0;
	listOfAtomsOnBackbone.reset4Loop();
	while (listOfAtomsOnBackbone.setNext4Loop())
	{
		atomsOnBackbone[ i ] = listOfAtomsOnBackbone.getEntity();
		i++;
	}

	return numBackboneAtoms;
}



void    V::GeometryTier::Residue::getAtomsOnSideChain( rg_dList<V::GeometryTier::Atom*>* listOfAtomsOnSideChain )
{
    if ( isStandardResidue() == rg_FALSE ) {
        return;
    }
    else {
        m_atoms.reset4Loop();
        while ( m_atoms.setNext4Loop() ) {
            V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        
            if( currAtom->getpChemicalProperties()->isOnBackBone() == rg_FALSE &&
                currAtom->getAtomNameFromInputFile() != " O  " ) {
                listOfAtomsOnSideChain->addTail( currAtom );
            }
        }    
    }    
}



rg_INT  V::GeometryTier::Residue::getAtomsOnSidechain( rg_dList<V::GeometryTier::Atom*>& atomsOnSidechain )
{
	if ( isStandardResidue() == rg_FALSE ) 
	{
		return 0;
	}
	else 
	{
		m_atoms.reset4Loop();
		while ( m_atoms.setNext4Loop() ) {
            V::GeometryTier::Atom* currAtom = m_atoms.getEntity();

			if( currAtom->getpChemicalProperties()->isOnBackBone() == rg_FALSE &&
				currAtom->getAtomNameFromInputFile() != " O  " ) 
			{
					atomsOnSidechain.addTail( currAtom );
			}
		}
		return atomsOnSidechain.getSize();
	}
}



rg_INT  V::GeometryTier::Residue::getAtomsOnSidechain(V::GeometryTier::Atom**& atomsOnSidechain )
{
	rg_dList<V::GeometryTier::Atom*> sidechainAtomsInList;
	rg_INT numAtoms = getAtomsOnSidechain(sidechainAtomsInList);

	atomsOnSidechain = new V::GeometryTier::Atom* [numAtoms];

	rg_INDEX index = 0;
	sidechainAtomsInList.reset4Loop();
	while (sidechainAtomsInList.setNext4Loop())
	{
		atomsOnSidechain[ index ] = sidechainAtomsInList.getEntity();
		index++;
	}

	return numAtoms;
}


void    V::GeometryTier::Residue::getAtomsNotOnSideChain( rg_dList<V::GeometryTier::Atom*>* listOfAtomsNotOnSideChain )
{
    if ( isStandardResidue() == rg_FALSE ) {
        return;
    }
    else {
        m_atoms.reset4Loop();
        while ( m_atoms.setNext4Loop() ) {
            V::GeometryTier::Atom* currAtom = m_atoms.getEntity();
        
            if( currAtom->getpChemicalProperties()->isOnBackBone() == rg_TRUE ||
                currAtom->getAtomNameFromInputFile() == " O  " ) {
                listOfAtomsNotOnSideChain->addTail( currAtom );
            }
        }    
    }    
}



V::GeometryTier::Residue*    V::GeometryTier::Residue::getResidueConnectedByNitrogenInAminoGroupOfAminoResidue() const
{
    V::GeometryTier::Atom* nitrogenInAminoGroup = getNitrogenInAminoGroupOfAminoResidue();
	rg_dList<V::GeometryTier::ChemicalBond*>* bondList = nitrogenInAminoGroup->getListChemicalBond();
	bondList->reset4Loop();
	while (bondList->setNext4Loop())
	{
        V::GeometryTier::ChemicalBond* currBond = bondList->getEntity();
        V::GeometryTier::Atom* currBondedAtom = currBond->getBondedAtom(nitrogenInAminoGroup);
		if(currBondedAtom->getResidue()->getID() != m_ID)
			return currBondedAtom->getResidue();
	}
	return rg_NULL;
}



V::GeometryTier::Residue*    V::GeometryTier::Residue::getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue() const
{
    V::GeometryTier::Atom* carbonInCarboxylGroup = getCarbonInCarboxylGroupOfAminoResidue();
	rg_dList<V::GeometryTier::ChemicalBond*>* bondList = carbonInCarboxylGroup->getListChemicalBond();
	bondList->reset4Loop();
	while (bondList->setNext4Loop())
	{
        V::GeometryTier::ChemicalBond* currBond = bondList->getEntity();
        V::GeometryTier::Atom* currBondedAtom = currBond->getBondedAtom(carbonInCarboxylGroup);
		if(currBondedAtom->getResidue()->getID() != m_ID)
			return currBondedAtom->getResidue();
	}
	return rg_NULL;
}



rg_BOOL     V::GeometryTier::Residue::isResidueType( const V::GeometryTier::ResidueType& typeOfResidue ) const
{
    rg_BOOL isCorrectResidueType = rg_FALSE;

    switch ( typeOfResidue ) {
    case AMINO_RESIDUE_TYPE:
        isCorrectResidueType = isAminoResidue();
        break;
    case DNA_RESIDUE_TYPE:
        isCorrectResidueType = isDNAResidue();
        break;
    case RNA_RESIDUE_TYPE:
        isCorrectResidueType = isRNAResidue();
        break;
    case NEUCLEIC_RESIDUE_TYPE:
        isCorrectResidueType = isNeucleicAcidResidue();
        break;
    case STANDARD_RESIDUE_TYPE:
        isCorrectResidueType = isStandardResidue();
        break;
    case NONE_STD_RESIDUE_TYPE:
        if ( !isStandardResidue() )
            isCorrectResidueType = rg_TRUE;
        break;
    default:
        break;
    }

    return isCorrectResidueType;
}



rg_BOOL     V::GeometryTier::Residue::isAminoResidue() const
{
    if( (rg_INT)m_residueCode > 0 && (rg_INT)m_residueCode <= NUM_OF_AMINOACIDS )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_BOOL     V::GeometryTier::Residue::isNeucleicAcidResidue() const
{
    if( (rg_INT)m_residueCode > NUM_OF_AMINOACIDS && (rg_INT)m_residueCode <= NUM_OF_AMINOACIDS+NUM_OF_NUCLEICACIDS )
        return rg_TRUE;
    else
        return rg_FALSE;

}


rg_BOOL     V::GeometryTier::Residue::isDNAResidue() const
{
    if( (rg_INT)m_residueCode > NUM_OF_AMINOACIDS && (rg_INT)m_residueCode <= NUM_OF_AMINOACIDS+NUM_OF_DNAS )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_BOOL     V::GeometryTier::Residue::isRNAResidue() const
{
    if( (rg_INT)m_residueCode > NUM_OF_AMINOACIDS+NUM_OF_DNAS && (rg_INT)m_residueCode <= NUM_OF_AMINOACIDS+NUM_OF_NUCLEICACIDS )
        return rg_TRUE;
    else
        return rg_FALSE;    
}



rg_BOOL     V::GeometryTier::Residue::isStandardResidue() const
{
    if( (rg_INT)m_residueCode > 0 && (rg_INT)m_residueCode <= NUM_OF_AMINOACIDS+NUM_OF_NUCLEICACIDS )
        return rg_TRUE;
    else
        return rg_FALSE;
}


rg_BOOL     V::GeometryTier::Residue::getBackBoneDihedralAngles(rg_REAL& phi, rg_REAL& psi) const
{
	rg_BOOL bNotTerminal = rg_TRUE;

	if(!isN_Terminal() && !isC_Terminal())
	{
		V::GeometryTier::Atom* C_Atom_prev = getResidueConnectedByNitrogenInAminoGroupOfAminoResidue()->getCarbonInCarboxylGroupOfAminoResidue();
		V::GeometryTier::Atom* N_Atom      = getNitrogenInAminoGroupOfAminoResidue();
		V::GeometryTier::Atom* CA_Atom     = getAlphaCarbonOfAminoResidue();
		V::GeometryTier::Atom* C_Atom      = getCarbonInCarboxylGroupOfAminoResidue();
		//V::GeometryTier::Atom* O_Atom      = getOxygenInCarboxylGroupOfAminoResidue();
		V::GeometryTier::Atom* N_Atom_next = getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue()->getNitrogenInAminoGroupOfAminoResidue();

		rg_Point3D C_center_prev = C_Atom_prev->getAtomBall().getCenter();
		rg_Point3D N_center      = N_Atom->getAtomBall().getCenter();
		rg_Point3D CA_center     = CA_Atom->getAtomBall().getCenter();
		rg_Point3D C_center      = C_Atom->getAtomBall().getCenter();
		//rg_Point3D O_center      = O_Atom->getAtomBall().getCenter();
		rg_Point3D N_center_next = N_Atom_next->getAtomBall().getCenter();

		phi = rg_GeoFunc::computeSignedDihedralAngle(C_center_prev, N_center, CA_center, C_center);
		phi = 180.0 / rg_PI * phi;
		//psi = rg_GeoFunc::computeSignedDihedralAngle(N_center, CA_center, C_center, O_center);
		psi = rg_GeoFunc::computeSignedDihedralAngle(N_center, CA_center, C_center, N_center_next);
		psi = 180.0 / rg_PI * psi;
	}
	else if(isN_Terminal())
	{
		V::GeometryTier::Atom* N_Atom      = getNitrogenInAminoGroupOfAminoResidue();
		V::GeometryTier::Atom* CA_Atom     = getAlphaCarbonOfAminoResidue();
		V::GeometryTier::Atom* C_Atom      = getCarbonInCarboxylGroupOfAminoResidue();
		//Atom* O_Atom      = getOxygenInCarboxylGroupOfAminoResidue();

        V::GeometryTier::Atom* N_Atom_next = getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue()->getNitrogenInAminoGroupOfAminoResidue();

		rg_Point3D N_center      = N_Atom->getAtomBall().getCenter();
		rg_Point3D CA_center     = CA_Atom->getAtomBall().getCenter();
		rg_Point3D C_center      = C_Atom->getAtomBall().getCenter();
		//rg_Point3D O_center      = O_Atom->getAtomBall().getCenter();
		rg_Point3D N_center_next = N_Atom_next->getAtomBall().getCenter();

		phi = DEFAULT_PHI_N_TERMINAL;
		//psi = rg_GeoFunc::computeSignedDihedralAngle(N_center, CA_center, C_center, O_center);
		psi = rg_GeoFunc::computeSignedDihedralAngle(N_center, CA_center, C_center, N_center_next);
		psi = 180.0 / rg_PI * psi;

		bNotTerminal = rg_FALSE;
	}
	else
	// isC_Terminal()
	{
		V::GeometryTier::Atom* C_Atom_prev = getResidueConnectedByNitrogenInAminoGroupOfAminoResidue()->getCarbonInCarboxylGroupOfAminoResidue();
		V::GeometryTier::Atom* N_Atom      = getNitrogenInAminoGroupOfAminoResidue();
		V::GeometryTier::Atom* CA_Atom     = getAlphaCarbonOfAminoResidue();
		V::GeometryTier::Atom* C_Atom      = getCarbonInCarboxylGroupOfAminoResidue();

		rg_Point3D C_center_prev = C_Atom_prev->getAtomBall().getCenter();
		rg_Point3D N_center      = N_Atom->getAtomBall().getCenter();
		rg_Point3D CA_center     = CA_Atom->getAtomBall().getCenter();
		rg_Point3D C_center      = C_Atom->getAtomBall().getCenter();

		phi = rg_GeoFunc::computeSignedDihedralAngle(C_center_prev, N_center, CA_center, C_center);
		phi = 180.0 / rg_PI * phi;
		psi = DEFAULT_PSI_C_TERMINAL;

		bNotTerminal = rg_FALSE;
	}

	return bNotTerminal;
}


rg_INT  V::GeometryTier::Residue::getSideChainDihedralAngles(rg_REAL*& dihedralAngles)
{
	rg_INT numDihedralAngles = 0;
	numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(m_residueCode);

	if( numDihedralAngles > 0 )
	{
		// initialize the array of atom pointers which contains the consecutive four atom pointers on side chain
        V::GeometryTier::Atom* consecutiveFourAtoms[ 4 ]
		= {rg_NULL,
			getNitrogenInAminoGroupOfAminoResidue(),
			getAlphaCarbonOfAminoResidue(),
			getBetaCarbonInSideChainOfAminoResidue()};

		rg_Point3D consecutiveFourAtomCenters[ 4 ] 
		= {rg_Point3D(0.0, 0.0, 0.0),
		   consecutiveFourAtoms[ 1 ]->getAtomBall().getCenter(),
		   consecutiveFourAtoms[ 2 ]->getAtomBall().getCenter(),
		   consecutiveFourAtoms[ 3 ]->getAtomBall().getCenter()};

		rg_INT dihedralAngleIndex = -1;
				
		// get atoms on side chain with remote indicator greater than equal to start remote in breath first order
		rg_dList<V::GeometryTier::Atom*> atomsOnSideChainInBreathFirstOrder;
        V::GeometryTier::RemoteIndicator  startRemote = GAMMA_REMOTE;
		getAtomsOnSideChainInBreathFirstOrder(this, atomsOnSideChainInBreathFirstOrder, startRemote);
		//rg_dList<Atom*> atomsOnSideChainInBreathFirstOrder;
		//RemoteIndicator  endRemote = UNK_REMOTE;
		//BranchDesignator endBranch = UNK_BRANCH;
		//FunctionsForRotamerLibrary::getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(getResidueCode(), endRemote, endBranch);
		//RemoteIndicator  startRemote = GAMMA_REMOTE;
		//getAtomsOnSideChainInBreathFirstOrder(this, atomsOnSideChainInBreathFirstOrder, startRemote, endRemote, endBranch);

        V::GeometryTier::RemoteIndicator prevRemote  = BETA_REMOTE;
		atomsOnSideChainInBreathFirstOrder.reset4Loop();

		while(atomsOnSideChainInBreathFirstOrder.setNext4Loop())
		{
            V::GeometryTier::Atom* currAtom = atomsOnSideChainInBreathFirstOrder.getEntity();
            V::GeometryTier::RemoteIndicator currRemote = currAtom->getpChemicalProperties()->getRemoteIndicator();
			// if current atom has a different (deeper) remote
			if(prevRemote < currRemote)
			{
				// get next dihedral angle
				dihedralAngleIndex++;
				// no more dihedral angle left
				if(dihedralAngleIndex >= numDihedralAngles)
					break;

				// shift elements in atom pointer array
				rg_INDEX i;
				for(i = 0;i < 3;i++)
				{
					consecutiveFourAtoms[ i ] = consecutiveFourAtoms[i + 1];
					consecutiveFourAtomCenters[ i ] = consecutiveFourAtomCenters[i + 1];
				}
				// replace fourth atom with new one
				consecutiveFourAtoms[ 3 ] = currAtom;
				consecutiveFourAtomCenters[ 3 ] = consecutiveFourAtoms[ 3 ]->getAtomBall().getCenter();
				// construct dihedral angle one by one
				dihedralAngles[dihedralAngleIndex] =
					180 * rg_GeoFunc::computeSignedDihedralAngle(consecutiveFourAtomCenters[ 0 ],
														         consecutiveFourAtomCenters[ 1 ],
														         consecutiveFourAtomCenters[ 2 ],
														         consecutiveFourAtomCenters[ 3 ]) / rg_PI;
				prevRemote = currRemote;
			}
			// same remote yet different branch
			// The case (prevRemote == currRemote) is only possible
			else
			{
				// second branch atom		
			}
		}
	}

	return numDihedralAngles;
}


rg_BOOL     V::GeometryTier::Residue::hasSideChainDihedralAngle() const
{
	// AMINO_RESIDUE_TYPE excluding ALA and GLY	
	if(isAminoResidue()) 
	{
		if(m_residueCode != ALA_AMINO_RESIDUE && m_residueCode != GLY_AMINO_RESIDUE)
		{
			return rg_TRUE;
		}
		else
		{
			return rg_FALSE;
		}
	}
	else
	{
		return rg_FALSE;
	}
}



rg_BOOL     V::GeometryTier::Residue::isN_Terminal() const
{
    V::GeometryTier::Residue* residueConnectedByNitrogenInAminogroup = getResidueConnectedByNitrogenInAminoGroupOfAminoResidue();
	if(residueConnectedByNitrogenInAminogroup != rg_NULL)
		return rg_FALSE;
	else
		return rg_TRUE;
}



rg_BOOL     V::GeometryTier::Residue::isC_Terminal() const
{
    V::GeometryTier::Residue* residueConnectedByCarbonInCarboxylgroup  = getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue();
	if(residueConnectedByCarbonInCarboxylgroup != rg_NULL)
		return rg_FALSE;
	else
		return rg_TRUE;
}



void    V::GeometryTier::Residue::changeConformationWith(const BackboneDihedralAnglePair& backboneDihedralAnglePair, const SidechainDihedralAngleSet& sideChainDihedralAngleSet)
{

}



rg_REAL     V::GeometryTier::Residue::computeSumOfPairwiseAtomXVolumesWithOtherResidue(V::GeometryTier::Residue* residue)
{
	rg_REAL sumOfPairwiseXVolume = 0.0;

	rg_dList<V::GeometryTier::Atom*> atomsOfTargetResidue;
	FunctionsForRotamerLibrary::getAtomsOfResidue(this, atomsOfTargetResidue);

	rg_dList<V::GeometryTier::Atom*> atomsOfOtherResidue;
	FunctionsForRotamerLibrary::getAtomsOfResidue(residue, atomsOfOtherResidue);

	atomsOfTargetResidue.reset4Loop();
	while (atomsOfTargetResidue.setNext4Loop())
	{
		rg_REAL xVolume = 0.0;
        V::GeometryTier::Atom* currAtom = atomsOfTargetResidue.getEntity();
		xVolume = computeSumPairwiseXVolOfAtomWithAtomSet(currAtom, atomsOfOtherResidue);
		sumOfPairwiseXVolume += xVolume;
	}
	return sumOfPairwiseXVolume;
}



rg_REAL V::GeometryTier::Residue::computeSumOfPairwiseAtomXVolumesOfSidechainWithOtherResidue(V::GeometryTier::Residue* residue)
{
	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	getAtomsOnSidechain(targetAtomsOfSidechain);

	rg_dList<V::GeometryTier::Atom*> atomsOfResidue;
	FunctionsForRotamerLibrary::getAtomsOfResidue(residue, atomsOfResidue);

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



rg_REAL V::GeometryTier::Residue::computeSumOfPairwiseAtomXVolumesOfSidechainWithItsBackbone()
{
	rg_dList<V::GeometryTier::Atom*> targetAtomsOfSidechain;
	getAtomsOnSidechain(targetAtomsOfSidechain);

	rg_dList<V::GeometryTier::Atom*> backboneAtoms;
	getAtomsOnBackbone(&backboneAtoms);

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



V::GeometryTier::Residue& V::GeometryTier::Residue::operator=( const V::GeometryTier::Residue& aResidue )
{
    if (this != &aResidue) {
        m_ID = aResidue.m_ID;
        m_residueCode = aResidue.m_residueCode;
        m_atoms = aResidue.m_atoms;
        m_chain = aResidue.m_chain;
        m_sequenceNumber = aResidue.m_sequenceNumber;
        m_residueName = aResidue.m_residueName;
    }

    return *this;
}


