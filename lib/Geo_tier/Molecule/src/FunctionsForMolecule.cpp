#include "FunctionsForMolecule.h"
#include "rg_Atom.h"
#include "Residue.h"
#include "rg_Molecule.h"
#include "rg_TMatrix3D.h"
#include "RotamerSetOfResidue.h"
#include "FunctionsForRotamerLibrary.h"
#include "BackboneDependentRotamerLibrary.h"
#include "PotentialEnergy.h"
#include "StringFunctions.h"
using namespace V::GeometryTier;


#ifdef BETASCP_STAND_ALONE_VERSION
// test
//const rg_REAL ALPHA_COEFF_COMBINATORIAL_SEARCH = 100.0;
//const rg_REAL RADIUS_INFLATION_RATE = 0.1;
//#include "ParamsForExperiments.h"
// May 11, 2016 by Joonghyun
//extern rg_REAL ALPHA_COEFF_COMBINATORIAL_SEARCH;
extern rg_REAL LAMBDA_COEFF_COMBINATORIAL_SEARCH;
extern rg_REAL RADIUS_INFLATION_RATE;
#endif


#include "FunctionsForMappingFromXVolToLJPotEnergy.h"
#include "FunctionsForSplineFittingLJPotEnergy.h"
#include "rg_GeoFunc.h"
#include <time.h>

#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;


rg_INT  V::GeometryTier::compareResidueBySequenceNumber( const void* residueA, const void* residueB )
{
    Residue* residue_A = *(Residue**) residueA;
    Residue* residue_B = *(Residue**) residueB;

    rg_INT sequneceNumberOfResidueA = residue_A->getSequenceNumber();
    rg_INT sequneceNumberOfResidueB = residue_B->getSequenceNumber();

    if ( sequneceNumberOfResidueA > sequneceNumberOfResidueB )
        return 1;
    else if ( sequneceNumberOfResidueA < sequneceNumberOfResidueB )
        return -1;
    else 
        return 0;
}



rg_INT V::GeometryTier::compareResidueByChainIDSequenceNumber( const void* residueA, const void* residueB )
{
	Residue* residue_A = *(Residue**) residueA;
	Residue* residue_B = *(Residue**) residueB;

	rg_INT chainIDSeqNumberA = getChainIDSeqNumOfResidue(residue_A);
	rg_INT chainIDSeqNumberB = getChainIDSeqNumOfResidue(residue_B);

	if(chainIDSeqNumberA > chainIDSeqNumberB)
		return 1;
	else if(chainIDSeqNumberA < chainIDSeqNumberB)
		return -1;
	else
		return 0;
}



rg_BOOL     V::GeometryTier::isThereBondBetween(V::GeometryTier::Atom& atom1, V::GeometryTier::Atom& atom2)
{
	rg_dList<ChemicalBond*>* chemicalBondList1 = atom1.getListChemicalBond();

	rg_BOOL bShareBond = rg_FALSE;
	chemicalBondList1->reset4Loop();
	while (chemicalBondList1->setNext4Loop())
	{
		ChemicalBond* currBond = chemicalBondList1->getEntity();
		if(currBond->getBondedAtom(&atom1) == &atom2)
		{
			bShareBond = rg_TRUE;
			break;
		}
	}
	return bShareBond;
}



void V::GeometryTier::getAtomsOnBackboneCorrToResidues(rg_dList<V::GeometryTier::Residue*>& targetListOfResidues, rg_dList<V::GeometryTier::Atom*>& listOfAtomsOnBackbone )
{
    targetListOfResidues.reset4Loop();
    while ( targetListOfResidues.setNext4Loop() ) 
	{        
        rg_dList<Atom*> atomsOnBackboneInResidue;
        targetListOfResidues.getEntity()->getAtomsOnBackbone( &atomsOnBackboneInResidue );
        
        atomsOnBackboneInResidue.reset4Loop();
        while ( atomsOnBackboneInResidue.setNext4Loop() ) {
            listOfAtomsOnBackbone.addTail( atomsOnBackboneInResidue.getEntity() );
        }
    }
}



void V::GeometryTier::getAtomsOnBackboneCorrTo18ResiduesExcludingGLYNALASortedBySequenceNumber(Molecule& molecule, rg_dList<V::GeometryTier::Atom*>& listOfAtomsOnBackbone )
{
    listOfAtomsOnBackbone.removeAll();

    rg_dList<Residue*> targetListOfResidues;
    get18ResiduesExcludingGLYNALASortedBySequenceNumber(molecule, targetListOfResidues );

    targetListOfResidues.reset4Loop();
    while ( targetListOfResidues.setNext4Loop() ) 
	{        
        rg_dList<Atom*> atomsOnBackboneInResidue;
        targetListOfResidues.getEntity()->getAtomsOnBackbone( &atomsOnBackboneInResidue );
        
        atomsOnBackboneInResidue.reset4Loop();
        while ( atomsOnBackboneInResidue.setNext4Loop() ) {
            listOfAtomsOnBackbone.addTail( atomsOnBackboneInResidue.getEntity() );
        }
    }
}



void V::GeometryTier::get18ResiduesExcludingGLYNALASortedBySequenceNumber(Molecule& molecule, rg_dList<Residue*>& targetListOfResidues )
{
	rg_dList<Residue>* residues = molecule.getResidues();
    rg_INT numOfResidues = residues->getSize();

    Residue** sortedResidues = new Residue*[numOfResidues];

    rg_INT i_residue = 0;
    residues->reset4Loop();
    while ( residues->setNext4Loop() ) {
        sortedResidues[i_residue] = residues->getpEntity();
        i_residue++;
    }

    qsort( (void *)sortedResidues, numOfResidues, sizeof(Residue*), compareResidueBySequenceNumber );
	
    targetListOfResidues.removeAll();
    
    for ( i_residue=0; i_residue<numOfResidues; i_residue++ ) 
	{
		ResidueCode code = sortedResidues[i_residue]->getResidueCode();
		if(code == VAL_AMINO_RESIDUE || code == LEU_AMINO_RESIDUE || code == ILE_AMINO_RESIDUE || code == MET_AMINO_RESIDUE || code == PRO_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE ||
		   code == TRP_AMINO_RESIDUE || code == SER_AMINO_RESIDUE || code == CYS_AMINO_RESIDUE || code == THR_AMINO_RESIDUE || code == GLN_AMINO_RESIDUE || code == ASN_AMINO_RESIDUE ||
		   code == HIS_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == GLU_AMINO_RESIDUE || code == ASP_AMINO_RESIDUE || code == LYS_AMINO_RESIDUE || code == ARG_AMINO_RESIDUE)
		{
			targetListOfResidues.addTail( sortedResidues[i_residue] );
		}
    }
    
    delete [] sortedResidues;
}



void V::GeometryTier::get19ResiduesExcludingGLYSortedBySequenceNumber(Molecule& molecule, rg_dList<Residue*>& targetListOfResidues )
{
	molecule.getResiduesSortedBySequenceNumber(targetListOfResidues);
	// remove NON_AMINO_RESIDUE_TYPE or GLY_RESIDUE from the list
	targetListOfResidues.reset4Loop();
	while(targetListOfResidues.setNext4Loop())
	{
		Residue* currResidue = targetListOfResidues.getEntity();
		// AMINO_RESIDUE_TYPE excluding GLY	
		if(!currResidue->isAminoResidue())
		{			
			targetListOfResidues.killCurrent();
		}
		else
		{
			if(currResidue->getResidueCode() 
				== GLY_AMINO_RESIDUE)
				targetListOfResidues.killCurrent();
		}
	}
}



rg_INT V::GeometryTier::deleteAtomsInSideChainsCorrTo18ResiduesExcludingGLYNALA(Molecule& molecule, rg_BOOL replaceResidueCodeToUnknown, const RemoteIndicator& startRemoteIndicator)
{
	rg_INT numDeletedResidues = 0;
	rg_dList<Residue>* residues = molecule.getResidues();
	residues->reset4Loop();
    while ( residues->setNext4Loop() ) 
	{
        Residue* currResidue = residues->getpEntity();		
		if( currResidue->isAminoResidue() )
		{
			ResidueCode code = currResidue->getResidueCode();
			if( code != ALA_AMINO_RESIDUE && code != GLY_AMINO_RESIDUE )
		//if(code == VAL_AMINO_RESIDUE || code == LEU_AMINO_RESIDUE || code == ILE_AMINO_RESIDUE || code == MET_AMINO_RESIDUE || code == PRO_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE ||
		//   code == TRP_AMINO_RESIDUE || code == SER_AMINO_RESIDUE || code == CYS_AMINO_RESIDUE || code == THR_AMINO_RESIDUE || code == GLN_AMINO_RESIDUE || code == ASN_AMINO_RESIDUE ||
		//   code == HIS_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == GLU_AMINO_RESIDUE || code == ASP_AMINO_RESIDUE || code == LYS_AMINO_RESIDUE || code == ARG_AMINO_RESIDUE)
			{
			   molecule.deleteAtomsInSideChain( replaceResidueCodeToUnknown, startRemoteIndicator, currResidue );
			   numDeletedResidues++;
			}
		}
	}
	return numDeletedResidues;
}



rg_INT V::GeometryTier::recoverMissingAtomsIn18ResiduesExcludingGLYNALAWithoutCoordinates(Molecule& molecule, rg_dList<Residue*>& recoveredResidues)
{
	rg_INT numRecoveredResidues = 0;
	rg_dList<Residue>* residues = molecule.getResidues();

    residues->reset4Loop();
    while ( residues->setNext4Loop() ) 
	{
        Residue* currResidue = residues->getpEntity();
		if( currResidue->isAminoResidue() )
		{
			ResidueCode code = currResidue->getResidueCode();
			if( code != ALA_AMINO_RESIDUE && code != GLY_AMINO_RESIDUE )
			//if(code == VAL_AMINO_RESIDUE || code == LEU_AMINO_RESIDUE || code == ILE_AMINO_RESIDUE || code == MET_AMINO_RESIDUE || code == PRO_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE ||
			//   code == TRP_AMINO_RESIDUE || code == SER_AMINO_RESIDUE || code == CYS_AMINO_RESIDUE || code == THR_AMINO_RESIDUE || code == GLN_AMINO_RESIDUE || code == ASN_AMINO_RESIDUE ||
			//   code == HIS_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == GLU_AMINO_RESIDUE || code == ASP_AMINO_RESIDUE || code == LYS_AMINO_RESIDUE || code == ARG_AMINO_RESIDUE)
			{
				molecule.recoverMissingAtomsInResidueWithoutCoordinates(currResidue);
				recoveredResidues.add(currResidue);
				numRecoveredResidues++;
			}
		}
	}
	return numRecoveredResidues;
}



rg_INT V::GeometryTier::recoverMissingAtomsIn18ResiduesExcludingGLYNALAWithoutCoordinates(Molecule& molecule)
{
	rg_INT numRecoveredResidues = 0;
	rg_dList<Residue>* residues = molecule.getResidues();

    residues->reset4Loop();
    while ( residues->setNext4Loop() ) 
	{
		Residue* currResidue = residues->getpEntity();
		if( currResidue->isAminoResidue() )
		{
			ResidueCode code = currResidue->getResidueCode();
			if( code != ALA_AMINO_RESIDUE && code != GLY_AMINO_RESIDUE )
				//if(code == VAL_AMINO_RESIDUE || code == LEU_AMINO_RESIDUE || code == ILE_AMINO_RESIDUE || code == MET_AMINO_RESIDUE || code == PRO_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE ||
				//   code == TRP_AMINO_RESIDUE || code == SER_AMINO_RESIDUE || code == CYS_AMINO_RESIDUE || code == THR_AMINO_RESIDUE || code == GLN_AMINO_RESIDUE || code == ASN_AMINO_RESIDUE ||
				//   code == HIS_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == GLU_AMINO_RESIDUE || code == ASP_AMINO_RESIDUE || code == LYS_AMINO_RESIDUE || code == ARG_AMINO_RESIDUE)
			{
				molecule.recoverMissingAtomsInResidueWithoutCoordinates(currResidue);
				numRecoveredResidues++;
			}
		}
	}
	return numRecoveredResidues;
}



rg_INT V::GeometryTier::getAtomsOfResidueWithGivenRemoteIndicatorNBranchDesignator(Residue* targetResidue, rg_dList<Atom*>& atomsOfResidueUnnecessaryAtomsRemoved, const RemoteIndicator& endRemote, const BranchDesignator& maxBranchDesignator)
{
	// get atoms of residue
	rg_dList<Atom*> atomsOfResidue;
	targetResidue->getAtoms(atomsOfResidue);

	atomsOfResidue.reset4Loop();
	while(atomsOfResidue.setNext4Loop())
	{
		Atom* currAtom = atomsOfResidue.getEntity();
		RemoteIndicator currRemoteIndicator = UNK_REMOTE;
		currRemoteIndicator = currAtom->getpChemicalProperties()->getRemoteIndicator();
		BranchDesignator currBranchDesignator = UNK_BRANCH;
		currBranchDesignator = currAtom->getpChemicalProperties()->getBrangeDesignator();


		if(currRemoteIndicator > endRemote)
		{
			atomsOfResidue.killCurrent();
		}
		else if(currBranchDesignator > maxBranchDesignator)
		{
			atomsOfResidue.killCurrent();
		}
		else
		{
			atomsOfResidueUnnecessaryAtomsRemoved.add( currAtom );
		}

		// Because some pdb file do not follow the order of remote indicator and branch designator,
		// we should go through the entire residue atoms
		//if(currRemoteIndicator == endRemote && currBranchDesignator == endBranchDesignator)
		//{
		//	break;
		//}		
	}

	return atomsOfResidueUnnecessaryAtomsRemoved.getSize();
}



rg_INT V::GeometryTier::getAtomsOnSidechainWithGivenRemoteIndicatorNBranchDesignator(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChainUnnecessaryAtomsRemoved, const RemoteIndicator& startRemote, const RemoteIndicator& endRemote, const BranchDesignator& maxBranchDesignator)
{
	// get atoms on side chain
	rg_dList<Atom*> atomsOnSideChain;
	targetResidue->getAtomsOnSideChain(&atomsOnSideChain);

	atomsOnSideChain.reset4Loop();
	while(atomsOnSideChain.setNext4Loop())
	{
		Atom* currAtom = atomsOnSideChain.getEntity();
		RemoteIndicator currRemoteIndicator = UNK_REMOTE;
		currRemoteIndicator = currAtom->getpChemicalProperties()->getRemoteIndicator();
		BranchDesignator currBranchDesignator = UNK_BRANCH;
		currBranchDesignator = currAtom->getpChemicalProperties()->getBrangeDesignator();

		if(currRemoteIndicator < startRemote)
		{
			atomsOnSideChain.killCurrent();
		}
		else if(currRemoteIndicator > endRemote)
		{
			atomsOnSideChain.killCurrent();
		}
		else if(currBranchDesignator > maxBranchDesignator)
		{
			atomsOnSideChain.killCurrent();
		}
		else
		{
			atomsOnSideChainUnnecessaryAtomsRemoved.add( currAtom );
		}

		// Because some pdb file do not follow the order of remote indicator and branch designator,
		// we should go through the entire sidechain atoms
//		if(currRemoteIndicator == endRemote && currBranchDesignator == endBranchDesignator)
//		{
//			break;
//		}
	}

	return atomsOnSideChainUnnecessaryAtomsRemoved.getSize();
}



rg_INT V::GeometryTier::getAtomsOnSideChainInBreathFirstOrder(
    Residue* targetResidue, 
    rg_dList<Atom*>& atomsOnSideChainInBreathFirstOrder, 
    const RemoteIndicator& startRemote /*= BETA_REMOTE*/)
{
	//RemoteIndicator  endRemote = UNK_REMOTE;
	//BranchDesignator endBranch = UNK_BRANCH;
	//FunctionsForRotamerLibrary::getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(targetResidue->getResidueCode(), endRemote, endBranch);
	RemoteIndicator  endRemote = UNK_REMOTE;
	BranchDesignator maxBranch = UNK_BRANCH;
	FunctionsForRotamerLibrary::getEndRemoteIndicatorNMaximumBranchDesignatorOfSidechainForResidue(targetResidue->getResidueCode(), endRemote, maxBranch);
	getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsOnSideChainInBreathFirstOrder, startRemote, endRemote, maxBranch);

	return atomsOnSideChainInBreathFirstOrder.getSize();

	/*
	// get atoms on side chain
	rg_dList<Atom*> atomsOnSideChain;
	targetResidue->getAtomsOnSideChain(&atomsOnSideChain);	

	// test
	if (targetResidue->getResidueCode() == GLN_AMINO_RESIDUE && targetResidue->getSequenceNumber() == 38)
	{
		int here = 1;
	}

	atomsOnSideChain.reset4Loop();
	while(atomsOnSideChain.setNext4Loop())
	{
		Atom* currAtom = atomsOnSideChain.getEntity();
		if(currAtom->getpChemicalProperties()->getRemoteIndicator() < startRemote)
		{
			atomsOnSideChain.killCurrent();
		}
		//else
		//	break;
		// 해당 residue의 side-chain에 해당하는 원자와 함께 수소원자도 갖고 있는 PDB file이 있기 때문에
		// atom list 전체를 확인해 보아야 한다.
	}

	return getAtomsOnSideChainInBreathFirstOrder(atomsOnSideChain, atomsOnSideChainInBreathFirstOrder);
	*/
}



/*
rg_INT V::GeometryTier::getAtomsOnSideChainInBreathFirstOrder(
                                V::GeometryTier::Residue* targetResidue, 
                                rg_dList<V::GeometryTier::Atom*>& atomsOnSideChainInBreathFirstOrder, 
                                const RemoteIndicator& startRemote)
{
	// get atoms on side chain
	rg_dList<Atom*> atomsOnSideChain;
	targetResidue->getAtomsOnSideChain(&atomsOnSideChain);
	atomsOnSideChain.reset4Loop();

	while(atomsOnSideChain.setNext4Loop())
	{
		Atom* currAtom = atomsOnSideChain.getEntity();
		if(currAtom->getpChemicalProperties()->getRemoteIndicator() < startRemote)
		{
			atomsOnSideChain.killCurrent();
		}
		else
			break;
	}

	// sort atoms by considering their Branch Designators
	
	// copy atom pointers
	rg_INT size = atomsOnSideChain.getSize();
	Atom** arrayOfAtoms = new Atom*[size];
	rg_INDEX i = 0;
	atomsOnSideChain.reset4Loop();
	while(atomsOnSideChain.setNext4Loop())
	{
		Atom* currAtom = atomsOnSideChain.getEntity();
		arrayOfAtoms[ i++ ] = currAtom;
	}

	// selection sort
	for(i = 0;i < size;i++)
	{
		rg_INDEX j;
		rg_INDEX sIndex = i + 1;
		rg_INDEX minIndex = i;
		for(j = sIndex;j < size;j++)
		{
			if(compareAtomsByRemoteIndicatorBranchDesignator(arrayOfAtoms[minIndex], 
				                                             arrayOfAtoms[ j ]      ) == 1)
			{
				minIndex = j;
			}
		}
		if(minIndex != i)
		{
			Atom* temp = arrayOfAtoms[ i ];
			arrayOfAtoms[ i ] = arrayOfAtoms[minIndex];
			arrayOfAtoms[minIndex] = temp;			
		}
	}

	// copy atom pointers into linked list
	for(i = 0;i < size;i++)
	{
		atomsOnSideChainInBreathFirstOrder.add(arrayOfAtoms[ i ]);
	}
	
	if(arrayOfAtoms != rg_NULL)
		delete [] arrayOfAtoms;

	return atomsOnSideChainInBreathFirstOrder.getSize();
}
*/



rg_INT V::GeometryTier::getAtomsOnSideChainInBreathFirstOrder(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChainInBreathFirstOrder, const RemoteIndicator& startRemote, const RemoteIndicator& endRemote, const BranchDesignator& maxBranchDesignator)
{
	rg_dList<Atom*> atomsOnSideChainUnnecessaryAtomsRemoved;
	getAtomsOnSidechainWithGivenRemoteIndicatorNBranchDesignator(targetResidue, atomsOnSideChainUnnecessaryAtomsRemoved, startRemote, endRemote, maxBranchDesignator);

	return getAtomsOnSideChainInBreathFirstOrder(atomsOnSideChainUnnecessaryAtomsRemoved, atomsOnSideChainInBreathFirstOrder);
}



rg_INT V::GeometryTier::getAtomsOnSideChainInBreathFirstOrder(rg_dList<Atom*>& atomsOnSideChain, rg_dList<Atom*>& atomsOnSideChainInBreathFirstOrder)
{
	// sort atoms by considering their Branch Designators

	// copy atom pointers
	rg_INT size = atomsOnSideChain.getSize();
	Atom** arrayOfAtoms = new Atom*[size];
	rg_INDEX i = 0;
	atomsOnSideChain.reset4Loop();
	while(atomsOnSideChain.setNext4Loop())
	{
		Atom* currAtom = atomsOnSideChain.getEntity();
		arrayOfAtoms[ i++ ] = currAtom;
	}

	// selection sort
	for(i = 0;i < size;i++)
	{
		rg_INDEX j;
		rg_INDEX sIndex = i + 1;
		rg_INDEX minIndex = i;
		for(j = sIndex;j < size;j++)
		{
			if(compareAtomsByRemoteIndicatorBranchDesignator(arrayOfAtoms[minIndex], 
				arrayOfAtoms[ j ]      ) == 1)
			{
				minIndex = j;
			}
		}
		if(minIndex != i)
		{
			Atom* temp = arrayOfAtoms[ i ];
			arrayOfAtoms[ i ] = arrayOfAtoms[minIndex];
			arrayOfAtoms[minIndex] = temp;			
		}
	}

	// copy atom pointers into linked list
	for(i = 0;i < size;i++)
	{
		atomsOnSideChainInBreathFirstOrder.add(arrayOfAtoms[ i ]);
	}

	if(arrayOfAtoms != rg_NULL)
		delete [] arrayOfAtoms;

	return atomsOnSideChainInBreathFirstOrder.getSize();
}



rg_INT V::GeometryTier::compareAtomsByRemoteIndicatorBranchDesignator(Atom* atom1, Atom* atom2)
{
	ChemicalPropertiesOfAtom* chemProp1 = atom1->getpChemicalProperties();
	ChemicalPropertiesOfAtom* chemProp2 = atom2->getpChemicalProperties();

	// compare remote indicators
	RemoteIndicator remote1 = chemProp1->getRemoteIndicator();
	RemoteIndicator remote2 = chemProp2->getRemoteIndicator();
	if(remote1 > remote2)
		return 1;
	if(remote1 < remote2)
		return -1;

	// if remote indicators are equal to each other
	// we need to check out branch designators
	BranchDesignator branch1 = chemProp1->getBrangeDesignator();
	BranchDesignator branch2 = chemProp2->getBrangeDesignator();
	if(branch1 > branch2)
		return 1;
	// branch1 < branch2
	// branch1 can never be equal to branch2
	else
		return -1;
}



// currently not called from anywhere

rg_INT V::GeometryTier::getAtomsOnSideChainInBreathFirstOrder(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChain)
{
	rg_BOOL bFound = rg_FALSE;
	RemoteIndicator currRemote = BETA_REMOTE;
	do
	{
		rg_dList<Atom*> atomsWithThisRemote;
		bFound = getAtomsSortedByBranchDesignator(currRemote, targetResidue, atomsWithThisRemote);
		currRemote = getNextRemoteIndicator(currRemote);
		atomsOnSideChain.append(atomsWithThisRemote);

	}while(bFound && currRemote != UNK_REMOTE);

	return atomsOnSideChain.getSize();
}



rg_BOOL  V::GeometryTier::getAtomsSortedByBranchDesignator(const RemoteIndicator& remote, Residue* targetResidue, rg_dList<Atom*>& atomsWithThisRemote)
{
	rg_dList<Atom*> atomsWithSameRemote;
	targetResidue->getAtoms(remote, atomsWithSameRemote);
	rg_INT size = atomsWithSameRemote.getSize();

	if(size == 0)
		return rg_FALSE;

	if(size == 1)
	{
		atomsWithThisRemote.add(atomsWithSameRemote.getFirstEntity());
		return rg_TRUE;
	}

	// sort atoms by considering their Branch Designators
	
	// copy atom pointers
	Atom** arrayOfAtoms = new Atom*[size];
	rg_INDEX i = 0;
	atomsWithSameRemote.reset4Loop();
	while(atomsWithSameRemote.setNext4Loop())
	{
		Atom* currAtom = atomsWithSameRemote.getEntity();
		arrayOfAtoms[ i++ ] = currAtom;
	}

	// selection sort
	for(i = 0;i < size;i++)
	{
		rg_INDEX j;
		rg_INDEX sIndex = i + 1;
		rg_INDEX minIndex = i;
		for(j = sIndex;j < size;j++)
		{
			if(arrayOfAtoms[minIndex]->getpChemicalProperties()->getBrangeDesignator() >
			   arrayOfAtoms[ j ]->getpChemicalProperties()->getBrangeDesignator()  )
			{
				minIndex = j;
			}
		}
		if(minIndex != i)
		{
			Atom* temp = arrayOfAtoms[ i ];
			arrayOfAtoms[ i ] = arrayOfAtoms[minIndex];
			arrayOfAtoms[minIndex] = temp;			
		}
	}

	// copy atom pointers into linked list
	for(i = 0;i < size;i++)
	{
		atomsWithThisRemote.add(arrayOfAtoms[ i ]);
	}
	
	if(arrayOfAtoms != rg_NULL)
		delete [] arrayOfAtoms;

	return rg_TRUE;
	
}



V::GeometryTier::Atom*  V::GeometryTier::getNextAtomOnSideChainInBreathFirstOrder(V::GeometryTier::Atom* targetAtom)
{
	Residue*                  targetResidue  = targetAtom->getResidue();
	ChemicalPropertiesOfAtom* targetChemProp = targetAtom->getpChemicalProperties();

	if(targetChemProp->isOnBackBone())
		return rg_NULL;

	RemoteIndicator  targetRemote = targetChemProp->getRemoteIndicator();
	BranchDesignator targetBranch = targetChemProp->getBrangeDesignator();

	if(targetRemote == TERMINATE_REMOTE)
		return rg_NULL;

	// get atoms which share remote indicator with target atom
	rg_dList<Atom*> atomsWithSameRemote;
	targetResidue->getAtoms(targetRemote, atomsWithSameRemote);
	rg_INT size = atomsWithSameRemote.getSize();

	// if only target atom has this remote indicator
	// then try to find the next remote atom with FIRST_BRANCH or UNK_BRANCH
	// NOTE: In this case, there could be no next atom
	if(size == 1)
	{
		return getFirstAtomOnSideChainWithThisRemoteIndicatorInBreathFirstOrder(getNextRemoteIndicator(targetRemote), targetResidue);
	}
	// size >= 2
	else
	{
		// check out the max. branch in this remote
		BranchDesignator maxBranch = UNK_BRANCH;
		atomsWithSameRemote.reset4Loop();
		while(atomsWithSameRemote.setNext4Loop())
		{
			Atom* currAtom = atomsWithSameRemote.getEntity();
			BranchDesignator currBranchDesg = currAtom->getpChemicalProperties()->getBrangeDesignator();
			if(maxBranch < currBranchDesg)
				maxBranch = currBranchDesg;
		}
		// if there exists another branch
		// then get the next branch atom
		if(targetBranch < maxBranch)
		{
			return targetResidue->getAtom(targetRemote, getNextBranchDesignator(targetBranch));
		}
		// targetBranch == maxBranch
		// target atom is the last branch of this remote
		// hence, try to find the next remote atom with FIRST_BRANCH or UNK_BRANCH
		// NOTE: In this case, there could be no next atom
		else
		{
			return getFirstAtomOnSideChainWithThisRemoteIndicatorInBreathFirstOrder(getNextRemoteIndicator(targetRemote), targetResidue);
		}
	}
}



V::GeometryTier::Atom* V::GeometryTier::getFirstAtomOnSideChainWithThisRemoteIndicatorInBreathFirstOrder(const RemoteIndicator& remote, V::GeometryTier::Residue* targetResidue)
{
	rg_dList<Atom*> atomsWithThisRemote;
	targetResidue->getAtoms(remote, atomsWithThisRemote);
	rg_INT size = atomsWithThisRemote.getSize();

	if(size >= 1)
	{
		// find the atom with smallest Branch designator
		// this code is required due to Tryptophan
		Atom* atomWithSmallestBranchDesignator = atomsWithThisRemote.popFront();
		BranchDesignator minBranch = atomWithSmallestBranchDesignator->getpChemicalProperties()->getBrangeDesignator();
		
		atomsWithThisRemote.reset4Loop();
		while(atomsWithThisRemote.setNext4Loop())
		{
			Atom* currAtom = atomsWithThisRemote.getEntity();
			BranchDesignator currBranch = currAtom->getpChemicalProperties()->getBrangeDesignator();
			if(minBranch > currBranch)
			{
				minBranch = currBranch;
				atomWithSmallestBranchDesignator = currAtom;
			}
		}
		return atomWithSmallestBranchDesignator;
	}
	return rg_NULL;
}



BranchDesignator V::GeometryTier::getNextBranchDesignator(const BranchDesignator& target)
{
	const rg_INT numDesig = 5;
	BranchDesignator branchDesg[ numDesig ] = { UNK_BRANCH, 
		                                        FIRST_BRANCH, 
												SECOND_BRANCH, 
												THIRD_BRANCH, 
												FOURTH_BRANCH };
	rg_INT nextBranchIndex = 0;
	rg_INT index;
	for(index = 0;index < numDesig;index++)
	{
		if(branchDesg[ index ] == target)
		{
			nextBranchIndex = index + 1;
			break;
		}
	}

	if(nextBranchIndex < numDesig)
		return branchDesg[nextBranchIndex];
	else
		return branchDesg[ 0 ];
}



RemoteIndicator  V::GeometryTier::getNextRemoteIndicator(const RemoteIndicator& target)
{
	const rg_INT numInd = 9;
	RemoteIndicator remoteInd[ numInd ] = { UNK_REMOTE,   
		                                    ALPHA_REMOTE, 
											BETA_REMOTE, 
											GAMMA_REMOTE, 
											DELTA_REMOTE, 
											EPSILON_REMOTE, 
											ZETA_REMOTE,   
											ETA_REMOTE, 
											TERMINATE_REMOTE };

	rg_INT nextRemoteIndex = 0;
	rg_INT index;
	for(index = 0;index < numInd;index++)
	{
		if(remoteInd[ index ] == target)
		{
			nextRemoteIndex = index + 1;
			break;
		}
	}

	if(nextRemoteIndex < numInd)
		return remoteInd[nextRemoteIndex];
	else
		return remoteInd[ 0 ];
}



void V::GeometryTier::getAtomNamesOfResidue(const ResidueCode& code, rg_dList<string>& atomNames)
{	
	// get start and end indices
	rg_INDEX start_index = -1;
	rg_INDEX end_index   = -1;

	if(code  == PRO_AMINO_RESIDUE || 
		code == PR0_AMINO_RESIDUE ||
		code == PRZ_AMINO_RESIDUE   )
	{
		//start_index = 2;
		start_index = 0;
		end_index = RESIDUE_FEATURES[code].numOfBonds - 3;
	}
	else
	{
		//start_index = 4;
		start_index = 0;
		end_index = RESIDUE_FEATURES[code].numOfBonds - 2;
	}

	// collect the types of atoms in target residue	
	rg_INDEX i;
	for(i = start_index;i <= end_index;i++)
	{
		string firstAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]);
		string secondAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]);

		//if( !doesStringContainHydrogenAtomName(firstAtomName) )
			atomNames.addWithoutSame( firstAtomName );

		//if( !doesStringContainHydrogenAtomName(secondAtomName) )
			atomNames.addWithoutSame( secondAtomName );
	}
}



void V::GeometryTier::getAtomNamesOfResidueExceptHydrogenNOXT(const ResidueCode& code, rg_dList<string>& atomNames)
{
    // get start and end indices
    rg_INDEX start_index = -1;
    rg_INDEX end_index   = -1;

    if(code  == PRO_AMINO_RESIDUE || 
        code == PR0_AMINO_RESIDUE ||
        code == PRZ_AMINO_RESIDUE   )
    {
        //start_index = 2;
        start_index = 0;
        end_index = RESIDUE_FEATURES[code].numOfBonds - 3;
    }
    else
    {
        //start_index = 4;
        start_index = 0;
        end_index = RESIDUE_FEATURES[code].numOfBonds - 2;
    }

    // collect the types of atoms in target residue	
    string OXT(" OXT");
    rg_INDEX i;
    for(i = start_index;i <= end_index;i++)
    {
        string firstAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]);
        string secondAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]);

        string firstCharacterOfFirst = firstAtomName.substr(0, 1);
        string firstCharacterOfSecond = secondAtomName.substr(0, 1);
        string secondCharacterOfFirst = firstAtomName.substr(1, 1);
        string secondCharacterOfSecond = secondAtomName.substr(1, 1);

        // "RESIDUE_FEATURES"에서 각 residue별로 구성하는 atom의 종류를 파악한 후
        // hydrogen인 경우와 OXT인 경우만 제외시킴
        rg_BOOL bFirstAtomHydrogen = rg_FALSE;
        if(firstCharacterOfFirst == string("H") || secondCharacterOfFirst == string("H"))
            bFirstAtomHydrogen = rg_TRUE;

        rg_BOOL bSecondAtomHydrogen = rg_FALSE;
        if(firstCharacterOfSecond == string("H") || secondCharacterOfSecond == string("H"))
            bSecondAtomHydrogen = rg_TRUE;

        //if( !doesStringContainHydrogenAtomName(firstAtomName)  &&  firstAtomName != OXT)
        //    atomNames.addWithoutSame( firstAtomName );

        //if( !doesStringContainHydrogenAtomName(secondAtomName) &&  secondAtomName != OXT)
        //    atomNames.addWithoutSame( secondAtomName );
        if( !bFirstAtomHydrogen  &&  firstAtomName != OXT)
            atomNames.addWithoutSame( firstAtomName );

        if( !bSecondAtomHydrogen &&  secondAtomName != OXT)
            atomNames.addWithoutSame( secondAtomName );
    }
}



void V::GeometryTier::getAtomNamesOfSidechain(const ResidueCode& code, rg_dList<string>& atomNames)
{
	// get start and end indices
	rg_INDEX start_index = -1;
	rg_INDEX end_index   = -1; 

	if( code == PRO_AMINO_RESIDUE || 
	    code == PR0_AMINO_RESIDUE ||
	    code == PRZ_AMINO_RESIDUE   )
	{
		start_index = 6;
		end_index = RESIDUE_FEATURES[code].numOfBonds - 3;
	}
	else
	{
		start_index = 8;
		end_index = RESIDUE_FEATURES[code].numOfBonds - 2;
	}

	// collect the types of atoms in target residue	
	rg_INDEX i;
	for(i = start_index;i <= end_index;i++)
	{
		string firstAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]);
		string secondAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]);

		//if( !doesStringContainHydrogenAtomName(firstAtomName) )
			atomNames.addWithoutSame( firstAtomName );

		//if( !doesStringContainHydrogenAtomName(secondAtomName) )
			atomNames.addWithoutSame( secondAtomName );
	}

	if(code == ALA_AMINO_RESIDUE)
	{
		atomNames.addWithoutSame( string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ 7 ][ 1 ]) );
	}
}



void V::GeometryTier::getAtomNamesOfSidechainExceptHydrogenNOXT(const ResidueCode& code, rg_dList<string>& atomNames)
{
    // get start and end indices
    rg_INDEX start_index = -1;
    rg_INDEX end_index   = -1; 

    if( code == PRO_AMINO_RESIDUE || 
        code == PR0_AMINO_RESIDUE ||
        code == PRZ_AMINO_RESIDUE   )
    {
        start_index = 6;
        end_index = RESIDUE_FEATURES[code].numOfBonds - 3; // OXT를 제외한 마지막 atom index
    }
    else
    {
        start_index = 8;
        end_index = RESIDUE_FEATURES[code].numOfBonds - 2; // OXT를 제외한 마지막 atom index
    }

    // collect the types of atoms in target residue	
    string OXT(" OXT");
    rg_INDEX i;
    for(i = start_index;i <= end_index;i++)
    {
        string firstAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]);
        string secondAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]);

        //if( !doesStringContainHydrogenAtomName(firstAtomName) )
        //atomNames.addWithoutSame( firstAtomName );

        //if( !doesStringContainHydrogenAtomName(secondAtomName) )
        //atomNames.addWithoutSame( secondAtomName );
        
        string firstCharacterOfFirst = firstAtomName.substr(0, 1);
        string firstCharacterOfSecond = secondAtomName.substr(0, 1);
        string secondCharacterOfFirst = firstAtomName.substr(1, 1);
        string secondCharacterOfSecond = secondAtomName.substr(1, 1);

        // "RESIDUE_FEATURES"에서 각 residue별로 구성하는 atom의 종류를 파악한 후
        // hydrogen인 경우와 OXT인 경우만 제외시킴
        rg_BOOL bFirstAtomHydrogen = rg_FALSE;
        if(firstCharacterOfFirst == string("H") || secondCharacterOfFirst == string("H"))
            bFirstAtomHydrogen = rg_TRUE;

        rg_BOOL bSecondAtomHydrogen = rg_FALSE;
        if(firstCharacterOfSecond == string("H") || secondCharacterOfSecond == string("H"))
            bSecondAtomHydrogen = rg_TRUE;

        //if( !doesStringContainHydrogenAtomName(firstAtomName)  &&  firstAtomName != OXT)
        //    atomNames.addWithoutSame( firstAtomName );

        //if( !doesStringContainHydrogenAtomName(secondAtomName) &&  secondAtomName != OXT)
        //    atomNames.addWithoutSame( secondAtomName );
        if( !bFirstAtomHydrogen  &&  firstAtomName != OXT)
            atomNames.addWithoutSame( firstAtomName );

        if( !bSecondAtomHydrogen &&  secondAtomName != OXT)
            atomNames.addWithoutSame( secondAtomName );

    }

    if(code == ALA_AMINO_RESIDUE)
    {
        atomNames.addWithoutSame( string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ 7 ][ 1 ]) );
    }
}



void V::GeometryTier::getAtomNamesOfResidue(Residue* residue, rg_dList<string>& atomNames)
{
	rg_dList<Atom*> atomsOfResidue;
	residue->getAtoms(atomsOfResidue);

	atomsOfResidue.reset4Loop();
	while (atomsOfResidue.setNext4Loop())
	{
		Atom* currAtom = atomsOfResidue.getEntity();
		string currAtomName = currAtom->getAtomNameFromInputFile();
		atomNames.add(currAtomName);
	}
}



void V::GeometryTier::getAtomNames(rg_dList<Atom*>& atoms, rg_dList<string>& atomNames)
{
    atoms.reset4Loop();
    while (atoms.setNext4Loop())
    {
        Atom* currAtom = atoms.getEntity();
        string currAtomName = currAtom->getAtomNameFromInputFile();
        atomNames.add(currAtomName);
    }
}



rg_BOOL V::GeometryTier::doesStringContainHydrogenAtomName( const string& atomName )
{
    if (atomName.empty()) {
        return rg_FALSE;
    }


    if (atomName[0] == 'H') {
        bool hasNotAnyAlpha = true;
        int  length = atomName.length();
        for (int i = 1; i < length; ++i) {
            if (isalpha(atomName[i])) {
                hasNotAnyAlpha = false;
            }
        }

        if (hasNotAnyAlpha) {
            return rg_TRUE;
        }
    }

    string secondCharacter = atomName.substr(1, 1);
    if(secondCharacter == string("H"))
        return rg_TRUE;
    else
        return rg_FALSE;

    /*
    string firstCharacter = atomName.substr(0, 1);
    string secondCharacter = atomName.substr(1, 1);
    if( (atomName.find(" H") != string::npos) ||
    //(atomName.find("HG") == string::npos  && firstCharacter.find("H") != string::npos) )
    (atomName.find("HG") == string::npos  && secondCharacter.find("H") != string::npos) )
    return rg_TRUE;
    else
    return rg_FALSE;
    */
}



rg_BOOL V::GeometryTier::doesStringContainHydrogenAtomName(const string& atomName, Atom* targetAtom)
{
    if (atomName.empty()) {
        return rg_FALSE;
    }

    string chemicalSymbol = atomName.substr(0, 2);
    if (chemicalSymbol == " H") {
        return rg_TRUE;
    }

    if (atomName[0] == 'H') {
        if (targetAtom->getResidue()->isAminoResidue()) {
            return rg_TRUE;
        }

        bool hasNotAnyAlpha = true;
        int  length = atomName.length();
        for (int i = 1; i < length; ++i) {
            if (isalpha(atomName[i])) {
                hasNotAnyAlpha = false;
            }
        }

        if (hasNotAnyAlpha) {
            return rg_TRUE;
        }
    }

    string secondCharacter = atomName.substr(1, 1);
    if (secondCharacter == string("H"))
        return rg_TRUE;
    else
        return rg_FALSE;

    /*
    string firstCharacter = atomName.substr(0, 1);
    string secondCharacter = atomName.substr(1, 1);
    if( (atomName.find(" H") != string::npos) ||
    //(atomName.find("HG") == string::npos  && firstCharacter.find("H") != string::npos) )
    (atomName.find("HG") == string::npos  && secondCharacter.find("H") != string::npos) )
    return rg_TRUE;
    else
    return rg_FALSE;
        */
}



rg_INT  V::GeometryTier::getChainIDSeqNumOfResidue(Residue* residue)
{
	rg_INT chainID = residue->getChain()->getID();
	rg_INT seqNum  = residue->getSequenceNumber();
	rg_INT chainSeqID = (chainID + 1) * 100000 + seqNum;
	return chainSeqID;
}



rg_INT  V::GeometryTier::getChainIDSeqNumOfResidue(const rg_INT& chainID, const rg_INT& seqNum)
{
    rg_INT chainSeqID = (chainID + 1) * 100000 + seqNum;
    return chainSeqID;
}

void    V::GeometryTier::getChainIDNSeqNumFromChainSeqID(const rg_INT& chainSeqID, rg_INT& chainID, rg_INT& seqNum)
{
	chainID = chainSeqID / 100000 - 1;
	seqNum  = chainSeqID % 100000    ;
}



void    V::GeometryTier::getResidueCorrToChainIDSeqNum(Molecule* protein, const rg_INT& chainIDSeqNum, Residue*& targetResidue)
{
    Residue** aminoResidues = rg_NULL;
    rg_INT numResidus = protein->getAminoResidues( aminoResidues );

    targetResidue = rg_NULL;
    rg_INT i;
    for(i = 0;i < numResidus;i++)
    {
        rg_INT currChainIDSeqNum = getChainIDSeqNumOfResidue(aminoResidues[ i ]);
        if(chainIDSeqNum == currChainIDSeqNum)
        {
            targetResidue = aminoResidues[ i ];
            break;
        }
    }

    if(aminoResidues != rg_NULL)
        delete [] aminoResidues;
}



void    V::GeometryTier::sortBackboneDihedralAnglePairs(rg_dList<BackboneDihedralAnglePair>& backboneDihedralAnglePairs)
{
	vector<BackboneDihedralAnglePair> sortedBackboneDihedralAnglePairs;
	backboneDihedralAnglePairs.reset4Loop();
	while (backboneDihedralAnglePairs.setNext4Loop())
	{
		BackboneDihedralAnglePair currAnglePair = backboneDihedralAnglePairs.popFront();
		sortedBackboneDihedralAnglePairs.push_back( currAnglePair );
	}
	
	// sort
	sort(sortedBackboneDihedralAnglePairs.begin(), sortedBackboneDihedralAnglePairs.end(), backboneDihedralAnglePairSortCriterion);
	// get result
	vector<BackboneDihedralAnglePair>::iterator itr;
	for(itr = sortedBackboneDihedralAnglePairs.begin(); itr != sortedBackboneDihedralAnglePairs.end();itr++)
		backboneDihedralAnglePairs.add( *itr );
}



bool    V::GeometryTier::backboneDihedralAnglePairSortCriterion(const BackboneDihedralAnglePair& bdp1, const BackboneDihedralAnglePair& bdp2)
{
	return bdp1.getChainIDSeqNum() < bdp2.getChainIDSeqNum();
}



void    V::GeometryTier::sortSidechainDihedralAngleSets(rg_dList<SidechainDihedralAngleSet>& sidechainDihedralAngleSets)
{
	vector<SidechainDihedralAngleSet> sortedSidechainDihedralAngleSets;
	sidechainDihedralAngleSets.reset4Loop();
	while (sidechainDihedralAngleSets.setNext4Loop())
	{
		SidechainDihedralAngleSet currSet = sidechainDihedralAngleSets.popFront();
		sortedSidechainDihedralAngleSets.push_back( currSet );
	}

	// sort
	sort(sortedSidechainDihedralAngleSets.begin(), sortedSidechainDihedralAngleSets.end(), sidechainDihedralAngleSetSortCriterion);
	// get result
	vector<SidechainDihedralAngleSet>::iterator itr;
	for(itr = sortedSidechainDihedralAngleSets.begin(); itr != sortedSidechainDihedralAngleSets.end(); itr++)
		sidechainDihedralAngleSets.add( *itr );
}



bool    V::GeometryTier::sidechainDihedralAngleSetSortCriterion(const SidechainDihedralAngleSet& sds1, const SidechainDihedralAngleSet& sds2)
{
	return sds1.getChainIDSeqNum() < sds2.getChainIDSeqNum();
}



void    V::GeometryTier::sortChainIDSeqNumRotLibIDPairSet(rg_dList<ChainIDSeqNumRotLibIDPair>& chainIDSeqNumRotLibIDPairSet)
{
	vector<ChainIDSeqNumRotLibIDPair> sortedChainIDSeqNumRotLibIDPairs;
	chainIDSeqNumRotLibIDPairSet.reset4Loop();
	while (chainIDSeqNumRotLibIDPairSet.setNext4Loop())
	{
		ChainIDSeqNumRotLibIDPair currPair = chainIDSeqNumRotLibIDPairSet.popFront();
		sortedChainIDSeqNumRotLibIDPairs.push_back( currPair );
	}

	// sort
	sort(sortedChainIDSeqNumRotLibIDPairs.begin(), sortedChainIDSeqNumRotLibIDPairs.end(), chainIDSeqNumRotLibIDPairSortCriterion);
	// get result
	vector<ChainIDSeqNumRotLibIDPair>::iterator itr;
	for(itr = sortedChainIDSeqNumRotLibIDPairs.begin(); itr != sortedChainIDSeqNumRotLibIDPairs.end(); itr++)
		chainIDSeqNumRotLibIDPairSet.add( *itr );
}



bool    V::GeometryTier::chainIDSeqNumRotLibIDPairSortCriterion(const ChainIDSeqNumRotLibIDPair& csr1, const ChainIDSeqNumRotLibIDPair& csr2)
{
	return csr1.getChainIDSeqNum() < csr2.getChainIDSeqNum();
}



rg_REAL     V::GeometryTier::computeVDWEnergyBetweenTwoAtoms(V::GeometryTier::Atom* firstAtom, V::GeometryTier::Atom* secondAtom)
{
    Sphere firstBall  = firstAtom->getAtomBall();
    Sphere secondBall = secondAtom->getAtomBall();

    rg_Point3D center1 = firstBall.getCenter();
    rg_Point3D center2 = secondBall.getCenter();

    rg_INT atomTypeInAmberForAtom1  = (rg_INT)firstAtom->getpChemicalProperties()->getAtomTypeInAmber();
    rg_INT atomTypeInAmberForAtom2 = (rg_INT)secondAtom->getpChemicalProperties()->getAtomTypeInAmber();

    rg_REAL squaredDistance = center2.squaredDistance(center1);
    rg_REAL sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

    rg_REAL energy = 0.0;
    energy = (   ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][0] / (sixthPoweredDistance*sixthPoweredDistance) )
        - ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][1] /  sixthPoweredDistance )                       );

    return energy;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyOfProtein(Molecule& protein)
{
	// get backbone atoms
	Atom** backboneAtoms = rg_NULL;
	rg_INT numBackboneAtoms
		= protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

	// get target residues
	Residue** sortedResidues = rg_NULL;
	rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

	// energy between sidechains and backbone
	rg_REAL nonBondedEnergyBetweenSidechainsNBackbone = 0.0;
	nonBondedEnergyBetweenSidechainsNBackbone = computeVDWESHBEnergyBetweenSidechainsNBackbone(sortedResidues, numResidues, backboneAtoms, numBackboneAtoms);

	// energy between sidechain pairs
	rg_REAL nonBondedEnergyBetweenSidechainPairs = 0.0;
	nonBondedEnergyBetweenSidechainPairs = computeVDWESHBEnergyBetweenSidechainPairs(sortedResidues, numResidues);

	rg_REAL nonBondedPotentialEnergy = 0.0;
	nonBondedPotentialEnergy = nonBondedEnergyBetweenSidechainsNBackbone + nonBondedEnergyBetweenSidechainPairs;

	if(backboneAtoms != rg_NULL)
		delete[] backboneAtoms;

	if(sortedResidues != rg_NULL)
		delete [] sortedResidues;

	return nonBondedPotentialEnergy;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainsNBackbone(V::GeometryTier::Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms)
{
	Atom* atomsOnBackbone = rg_NULL;
	atomsOnBackbone = new Atom[numBackboneAtoms];
	rg_INDEX i;
	for (i = 0;i < numBackboneAtoms;i++)
	{
		atomsOnBackbone[ i ] = *backboneAtoms[ i ];
	}

	const rg_INT flagsOfEnergyTerms = VDW_TERM_ON | ELEC_TERM_ON | HBOND_TERM_ON;
	PotentialEnergy PE;
	PE.setFlagsOfAppliedEnergyTerms(flagsOfEnergyTerms);

	rg_REAL energyBetweenSidechainsNBackbone = 0.0;
		
	for(i = 0;i < numResidues;i++)
	{
		rg_INT numSidechainAtoms = 0;
		Atom** sidechainAtoms = rg_NULL;		
		numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);
	
		rg_REAL energyBetweenCurrentSidechainNBackbone = 0.0;
		energyBetweenCurrentSidechainNBackbone 
			= PE.computeAndGetPotentialEnergyForTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
		energyBetweenSidechainsNBackbone += energyBetweenCurrentSidechainNBackbone;

		if(sidechainAtoms != rg_NULL)
			delete [] sidechainAtoms;
	}

	if(atomsOnBackbone != rg_NULL)
		delete [] atomsOnBackbone;

	return energyBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainPairs(V::GeometryTier::Residue** residues, const rg_INT& numResidues)
{
	const rg_INT flagsOfEnergyTerms = VDW_TERM_ON | ELEC_TERM_ON | HBOND_TERM_ON;
	PotentialEnergy PE;
	PE.setFlagsOfAppliedEnergyTerms(flagsOfEnergyTerms);

	rg_REAL energyBetweenSidechainPairs = 0.0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_INT numSidechainAtoms_i = 0;
		Atom** sidechainAtoms_i = rg_NULL;
		numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

		Atom* atomsOnSidechain_i = rg_NULL;
		atomsOnSidechain_i = new Atom[numSidechainAtoms_i];
		rg_INDEX k;
		for (k = 0;k < numSidechainAtoms_i;k++)
		{
			atomsOnSidechain_i[ k ] = *sidechainAtoms_i[ k ];
		}

		rg_INDEX j;
		for(j = i + 1;j < numResidues;j++)
		{
			rg_INT numSidechainAtoms_j = 0;
			Atom** sidechainAtoms_j = rg_NULL;
			numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

			rg_REAL energyBetweenCurrentSidechainPairs = 0;
			energyBetweenCurrentSidechainPairs 
				= PE.computeAndGetPotentialEnergyForTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
			energyBetweenSidechainPairs += energyBetweenCurrentSidechainPairs;

			if(sidechainAtoms_j != rg_NULL)
				delete [] sidechainAtoms_j;
		}

		if(sidechainAtoms_i != rg_NULL)
			delete [] sidechainAtoms_i;

		if(atomsOnSidechain_i != rg_NULL)
			delete [] atomsOnSidechain_i;
	}

	return energyBetweenSidechainPairs;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyOfComputedProteinStructureForSCPWithMissingSideChainAtoms(Molecule& targetProtein, Molecule& referenceProtein, rg_REAL& VDWEnergy, rg_REAL& electrostaticEnergy, rg_REAL& hydrogenBondingEnergy)
{
    // get target residues
    Residue** sortedTargetResidues = rg_NULL;
    rg_INT numTargetResidues = targetProtein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedTargetResidues);

	// get reference residues
    Residue** sortedReferenceResidues = rg_NULL;
    rg_INT numReferenceResidues = referenceProtein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedReferenceResidues);

	// get backbone atoms corresponding to reference protein (reference protein might have some missing backbone atoms)
	rg_dList<Atom*> commonBackboneAtoms;
	rg_INDEX i;
	for(i = 0;i < numTargetResidues;i++)
	{
		rg_dList<Atom*> atomsOnBackboneInTargetResidue;
		sortedTargetResidues[ i ]->getAtomsOnBackbone(&atomsOnBackboneInTargetResidue);

		rg_dList<Atom*> atomsOnBackboneInReferenceResidue;
		sortedReferenceResidues[ i ]->getAtomsOnBackbone(&atomsOnBackboneInReferenceResidue);
				
		atomsOnBackboneInTargetResidue.reset4Loop();
		atomsOnBackboneInReferenceResidue.reset4Loop();
		while(atomsOnBackboneInTargetResidue.setNext4Loop() && atomsOnBackboneInReferenceResidue.setNext4Loop())
		{
			commonBackboneAtoms.add(atomsOnBackboneInTargetResidue.getEntity());
		}
	}
	rg_INT numBackboneAtoms = 0;
	numBackboneAtoms = commonBackboneAtoms.getSize();
	Atom** backboneAtoms = rg_NULL;
	backboneAtoms = new Atom*[numBackboneAtoms];
	i = 0;
	commonBackboneAtoms.reset4Loop();
	while(commonBackboneAtoms.setNext4Loop())
	{
		backboneAtoms[ i++ ] = commonBackboneAtoms.getEntity();
	}

    //Atom** backboneAtoms = rg_NULL;
    //rg_INT numBackboneAtoms
    //    = targetProtein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // energy between sidechains and backbone
    rg_REAL nonBondedEnergyBetweenSidechainsNBackbone = 0.0;
    rg_REAL VDWEnergyBetweenSidechainsNBackbone = 0.0;
    rg_REAL electrostaticEnergyBetweenSidechainsNBackbone = 0.0;
    rg_REAL hydrogenBondingEnergyBetweenSidechainsNBackbone = 0.0;
    nonBondedEnergyBetweenSidechainsNBackbone = computeVDWESHBEnergyBetweenSidechainsNBackbone(sortedTargetResidues, 
                                                                                               numTargetResidues, 
																							   sortedReferenceResidues,
																							   numReferenceResidues,
                                                                                               backboneAtoms, 
                                                                                               numBackboneAtoms,
                                                                                               VDWEnergyBetweenSidechainsNBackbone,
                                                                                               electrostaticEnergyBetweenSidechainsNBackbone,
                                                                                               hydrogenBondingEnergyBetweenSidechainsNBackbone);

    // energy between sidechain pairs
    rg_REAL nonBondedEnergyBetweenSidechainPairs = 0.0;
    rg_REAL VDWEnergyBetweenSidechainPairs = 0.0;
    rg_REAL electrostaticEnergyBetweenSidechainPairs = 0.0;
    rg_REAL hydrogenBondingEnergyBetweenSidechainPairs = 0.0;
    nonBondedEnergyBetweenSidechainPairs = computeVDWESHBEnergyBetweenSidechainPairs(sortedTargetResidues, 
                                                                                     numTargetResidues,
																					 sortedReferenceResidues,
																					 numReferenceResidues,
                                                                                     VDWEnergyBetweenSidechainPairs,
                                                                                     electrostaticEnergyBetweenSidechainPairs,
                                                                                     hydrogenBondingEnergyBetweenSidechainPairs);

    VDWEnergy = electrostaticEnergy = hydrogenBondingEnergy = 0.0;
    VDWEnergy = VDWEnergyBetweenSidechainsNBackbone + VDWEnergyBetweenSidechainPairs;
    electrostaticEnergy = electrostaticEnergyBetweenSidechainsNBackbone + electrostaticEnergyBetweenSidechainPairs;
    hydrogenBondingEnergy = hydrogenBondingEnergyBetweenSidechainsNBackbone + hydrogenBondingEnergyBetweenSidechainPairs;

    rg_REAL nonBondedPotentialEnergy = 0.0;
    nonBondedPotentialEnergy = nonBondedEnergyBetweenSidechainsNBackbone + nonBondedEnergyBetweenSidechainPairs;

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    if(sortedTargetResidues != rg_NULL)
        delete [] sortedTargetResidues;

	if(sortedReferenceResidues != rg_NULL)
		delete [] sortedReferenceResidues;

    return nonBondedPotentialEnergy;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainsNBackbone(
                                                        Residue** sortedTargetResidues,
                                                        const rg_INT& numTargetResidues, 
													    Residue** sortedReferenceResidues,
													    const rg_INT& numReferenceResidues,
                                                        Atom** backboneAtoms, 
                                                        const rg_INT& numBackboneAtoms,
                                                        rg_REAL& VDWEnergyBetweenSidechainsNBackbone,
                                                        rg_REAL& electrostaticEnergyBetweenSidechainsNBackbone,
                                                        rg_REAL& hydrogenBondingEnergyBetweenSidechainsNBackbone)
{
    Atom* atomsOnBackbone = rg_NULL;
    atomsOnBackbone = new Atom[numBackboneAtoms];
    rg_INDEX i;
    for (i = 0;i < numBackboneAtoms;i++)
    {
        atomsOnBackbone[ i ] = *backboneAtoms[ i ];
    }

    PotentialEnergy PE;
    
    VDWEnergyBetweenSidechainsNBackbone = electrostaticEnergyBetweenSidechainsNBackbone = hydrogenBondingEnergyBetweenSidechainsNBackbone = 0.;
    for(i = 0;i < numTargetResidues;i++)
    {
		// get only atoms that are shared by target and reference protein structures ///////////////////////////
		rg_dList<Atom*> sidechainAtomsOfTargetResidue;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedTargetResidues[ i ], sidechainAtomsOfTargetResidue );
		rg_dList<Atom*> sidechainAtomsOfReferenceResidue;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedReferenceResidues[ i ], sidechainAtomsOfReferenceResidue );

		rg_dList<Atom*> commonAtoms;
		sidechainAtomsOfTargetResidue.reset4Loop();
		sidechainAtomsOfReferenceResidue.reset4Loop();
		while (sidechainAtomsOfTargetResidue.setNext4Loop() && sidechainAtomsOfReferenceResidue.setNext4Loop())
		{
			Atom* currAtom = rg_NULL;
			currAtom = sidechainAtomsOfTargetResidue.getEntity();
			commonAtoms.add( currAtom );
		}
		rg_INT numSidechainAtoms = commonAtoms.getSize();
        Atom** sidechainAtoms = rg_NULL;
		sidechainAtoms = new Atom*[numSidechainAtoms];

		rg_INDEX j = 0;
		commonAtoms.reset4Loop();
		while(commonAtoms.setNext4Loop())
		{
			sidechainAtoms[ j++ ] = commonAtoms.getEntity();
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		
        VDWEnergyBetweenSidechainsNBackbone 
            //+= PE.computeVanDerWaalsBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
            += computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms);
        electrostaticEnergyBetweenSidechainsNBackbone
            += PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);

        //hydrogenBondingEnergyBetweenSidechainsNBackbone
        //    += PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }

    if(atomsOnBackbone != rg_NULL)
        delete [] atomsOnBackbone;

    rg_REAL energyBetweenSidechainsNBackbone = 0.0;
    energyBetweenSidechainsNBackbone = VDWEnergyBetweenSidechainsNBackbone + electrostaticEnergyBetweenSidechainsNBackbone + hydrogenBondingEnergyBetweenSidechainsNBackbone;

    return energyBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainPairs(
                                                  Residue** sortedTargetResidues,
                                                  const rg_INT& numTargetResidues,
												  Residue** sortedReferenceResidues,
												  const rg_INT& numReferenceResidues,
                                                  rg_REAL& VDWEnergyBetweenSidechainPairs,
                                                  rg_REAL& electrostaticEnergyBetweenSidechainPairs,
                                                  rg_REAL& hydrogenBondingEnergyBetweenSidechainPairs)
{
    PotentialEnergy PE;

    VDWEnergyBetweenSidechainPairs = electrostaticEnergyBetweenSidechainPairs = hydrogenBondingEnergyBetweenSidechainPairs = 0.0;

    rg_INDEX i;
    for(i = 0;i < numTargetResidues;i++)
    {
		// get only atoms that are shared by target and reference protein structures ///////////////////////////
		rg_dList<Atom*> sidechainAtomsOfTargetResidue_i;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedTargetResidues[ i ], sidechainAtomsOfTargetResidue_i );
		rg_dList<Atom*> sidechainAtomsOfReferenceResidue_i;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedReferenceResidues[ i ], sidechainAtomsOfReferenceResidue_i );

		rg_dList<Atom*> commonAtoms_i;
		sidechainAtomsOfTargetResidue_i.reset4Loop();
		sidechainAtomsOfReferenceResidue_i.reset4Loop();
		while (sidechainAtomsOfTargetResidue_i.setNext4Loop() && sidechainAtomsOfReferenceResidue_i.setNext4Loop())
		{
			Atom* currAtom = rg_NULL;
			currAtom = sidechainAtomsOfTargetResidue_i.getEntity();
			commonAtoms_i.add( currAtom );
		}
		rg_INT numSidechainAtoms_i = commonAtoms_i.getSize();
        Atom** sidechainAtoms_i = rg_NULL;
		sidechainAtoms_i = new Atom*[numSidechainAtoms_i];

		rg_INDEX k = 0;
		commonAtoms_i.reset4Loop();
		while(commonAtoms_i.setNext4Loop())
		{
			sidechainAtoms_i[ k++ ] = commonAtoms_i.getEntity();
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////

        Atom* atomsOnSidechain_i = rg_NULL;
        atomsOnSidechain_i = new Atom[numSidechainAtoms_i];
        rg_INDEX l;
        for (l = 0;l < numSidechainAtoms_i;l++)
        {
            atomsOnSidechain_i[ l ] = *sidechainAtoms_i[ l ];
        }

        rg_INDEX j;
        for(j = i + 1;j < numTargetResidues;j++)
        {
			// get only atoms that are shared by target and reference protein structures ///////////////////////////
			rg_dList<Atom*> sidechainAtomsOfTargetResidue_j;
			FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedTargetResidues[ j ], sidechainAtomsOfTargetResidue_j );
			rg_dList<Atom*> sidechainAtomsOfReferenceResidue_j;
			FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedReferenceResidues[ j ], sidechainAtomsOfReferenceResidue_j );

			rg_dList<Atom*> commonAtoms_j;
			sidechainAtomsOfTargetResidue_j.reset4Loop();
			sidechainAtomsOfReferenceResidue_j.reset4Loop();
			while (sidechainAtomsOfTargetResidue_j.setNext4Loop() && sidechainAtomsOfReferenceResidue_j.setNext4Loop())
			{
				Atom* currAtom = rg_NULL;
				currAtom = sidechainAtomsOfTargetResidue_j.getEntity();
				commonAtoms_j.add( currAtom );
			}
			rg_INT numSidechainAtoms_j = commonAtoms_j.getSize();
			Atom** sidechainAtoms_j = rg_NULL;
			sidechainAtoms_j = new Atom*[numSidechainAtoms_j];

			rg_INDEX k = 0;
			commonAtoms_j.reset4Loop();
			while(commonAtoms_j.setNext4Loop())
			{
				sidechainAtoms_j[ k++ ] = commonAtoms_j.getEntity();
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
            VDWEnergyBetweenSidechainPairs 
                //+= PE.computeVanDerWaalsBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
                += computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms_j, numSidechainAtoms_j, sidechainAtoms_i, numSidechainAtoms_i);
            electrostaticEnergyBetweenSidechainPairs
                += PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
            //hydrogenBondingEnergyBetweenSidechainPairs
            //    += PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;

        if(atomsOnSidechain_i != rg_NULL)
            delete [] atomsOnSidechain_i;
    }

    rg_REAL energyBetweenSidechainPairs = 0.0;
    energyBetweenSidechainPairs = VDWEnergyBetweenSidechainPairs + electrostaticEnergyBetweenSidechainPairs + hydrogenBondingEnergyBetweenSidechainPairs;

    return energyBetweenSidechainPairs;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyOfProtein(Molecule& protein, rg_REAL& VDWEnergy, rg_REAL& electrostaticEnergy, rg_REAL& hydrogenBondingEnergy)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** sortedResidues = rg_NULL;
    rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

    // energy between sidechains and backbone
    rg_REAL nonBondedEnergyBetweenSidechainsNBackbone = 0.0;
    rg_REAL VDWEnergyBetweenSidechainsNBackbone = 0.0;
    rg_REAL electrostaticEnergyBetweenSidechainsNBackbone = 0.0;
    rg_REAL hydrogenBondingEnergyBetweenSidechainsNBackbone = 0.0;
    nonBondedEnergyBetweenSidechainsNBackbone = computeVDWESHBEnergyBetweenSidechainsNBackbone(sortedResidues, 
                                                                                               numResidues, 
                                                                                               backboneAtoms, 
                                                                                               numBackboneAtoms,
                                                                                               VDWEnergyBetweenSidechainsNBackbone,
                                                                                               electrostaticEnergyBetweenSidechainsNBackbone,
                                                                                               hydrogenBondingEnergyBetweenSidechainsNBackbone);

    // energy between sidechain pairs
    rg_REAL nonBondedEnergyBetweenSidechainPairs = 0.0;
    rg_REAL VDWEnergyBetweenSidechainPairs = 0.0;
    rg_REAL electrostaticEnergyBetweenSidechainPairs = 0.0;
    rg_REAL hydrogenBondingEnergyBetweenSidechainPairs = 0.0;
    nonBondedEnergyBetweenSidechainPairs = computeVDWESHBEnergyBetweenSidechainPairs(sortedResidues, 
                                                                                     numResidues,
                                                                                     VDWEnergyBetweenSidechainPairs,
                                                                                     electrostaticEnergyBetweenSidechainPairs,
                                                                                     hydrogenBondingEnergyBetweenSidechainPairs);

    VDWEnergy = electrostaticEnergy = hydrogenBondingEnergy = 0.0;
    VDWEnergy = VDWEnergyBetweenSidechainsNBackbone + VDWEnergyBetweenSidechainPairs;
    electrostaticEnergy = electrostaticEnergyBetweenSidechainsNBackbone + electrostaticEnergyBetweenSidechainPairs;
    hydrogenBondingEnergy = hydrogenBondingEnergyBetweenSidechainsNBackbone + hydrogenBondingEnergyBetweenSidechainPairs;

    rg_REAL nonBondedPotentialEnergy = 0.0;
    nonBondedPotentialEnergy = nonBondedEnergyBetweenSidechainsNBackbone + nonBondedEnergyBetweenSidechainPairs;

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    if(sortedResidues != rg_NULL)
        delete [] sortedResidues;

    return nonBondedPotentialEnergy;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainsNBackbone(
                                                       Residue** residues,
                                                       const rg_INT& numResidues, 
                                                       Atom** backboneAtoms, 
                                                       const rg_INT& numBackboneAtoms,
                                                       rg_REAL& VDWEnergyBetweenSidechainsNBackbone,
                                                       rg_REAL& electrostaticEnergyBetweenSidechainsNBackbone,
                                                       rg_REAL& hydrogenBondingEnergyBetweenSidechainsNBackbone)
{
    Atom* atomsOnBackbone = rg_NULL;
    atomsOnBackbone = new Atom[numBackboneAtoms];
    rg_INDEX i;
    for (i = 0;i < numBackboneAtoms;i++)
    {
        atomsOnBackbone[ i ] = *backboneAtoms[ i ];
    }

    PotentialEnergy PE;
    
    VDWEnergyBetweenSidechainsNBackbone = electrostaticEnergyBetweenSidechainsNBackbone = hydrogenBondingEnergyBetweenSidechainsNBackbone = 0.;

    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        VDWEnergyBetweenSidechainsNBackbone 
            //+= PE.computeVanDerWaalsBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
            += computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms);
        electrostaticEnergyBetweenSidechainsNBackbone
            += PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);

        //hydrogenBondingEnergyBetweenSidechainsNBackbone
        //    += PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }

    if(atomsOnBackbone != rg_NULL)
        delete [] atomsOnBackbone;

    rg_REAL energyBetweenSidechainsNBackbone = 0.0;
    energyBetweenSidechainsNBackbone = VDWEnergyBetweenSidechainsNBackbone + electrostaticEnergyBetweenSidechainsNBackbone + hydrogenBondingEnergyBetweenSidechainsNBackbone;

    return energyBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainPairs(
                                                  Residue** residues,
                                                  const rg_INT& numResidues,
                                                  rg_REAL& VDWEnergyBetweenSidechainPairs,
                                                  rg_REAL& electrostaticEnergyBetweenSidechainPairs,
                                                  rg_REAL& hydrogenBondingEnergyBetweenSidechainPairs)
{
    PotentialEnergy PE;

    VDWEnergyBetweenSidechainPairs = electrostaticEnergyBetweenSidechainPairs = hydrogenBondingEnergyBetweenSidechainPairs = 0.0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        Atom* atomsOnSidechain_i = rg_NULL;
        atomsOnSidechain_i = new Atom[numSidechainAtoms_i];
        rg_INDEX k;
        for (k = 0;k < numSidechainAtoms_i;k++)
        {
            atomsOnSidechain_i[ k ] = *sidechainAtoms_i[ k ];
        }

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);
                        
            VDWEnergyBetweenSidechainPairs 
                //+= PE.computeVanDerWaalsBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
                += computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms_j, numSidechainAtoms_j, sidechainAtoms_i, numSidechainAtoms_i);
            electrostaticEnergyBetweenSidechainPairs
                += PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
            //hydrogenBondingEnergyBetweenSidechainPairs
            //    += PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;

        if(atomsOnSidechain_i != rg_NULL)
            delete [] atomsOnSidechain_i;
    }

    rg_REAL energyBetweenSidechainPairs = 0.0;
    energyBetweenSidechainPairs = VDWEnergyBetweenSidechainPairs + electrostaticEnergyBetweenSidechainPairs + hydrogenBondingEnergyBetweenSidechainPairs;

    return energyBetweenSidechainPairs;
}



rg_INT V::GeometryTier::computeVDWESHBEnergyOfResiduesWithBackboneNOtherSidechains(
                                                                  Molecule& protein,
                                                                  rg_REAL*& engergyWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                  rg_REAL*& engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                  rg_REAL&  energyOfProtein)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** aminoResidues_orderedByChainIDSeqNum = rg_NULL;
    rg_INT numAminoResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResidues_orderedByChainIDSeqNum);

    engergyWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    rg_INDEX i;
    for (i = 0;i < numAminoResidues;i++)
    {
        engergyWithBackboneOfResidue_orderedByChainIDSeqNum[ i ]        = 0.0;
        engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0.0;
    }

    // energy between sidechains and backbone
    rg_REAL nonBondedEnergyBetweenSidechainsNBackbone = 0.0;
    nonBondedEnergyBetweenSidechainsNBackbone = computeVDWESHBEnergyBetweenSidechainsNBackbone(aminoResidues_orderedByChainIDSeqNum, 
                                                                                               numAminoResidues, 
                                                                                               backboneAtoms, 
                                                                                               numBackboneAtoms,
                                                                                               engergyWithBackboneOfResidue_orderedByChainIDSeqNum);

    // energy between sidechain pairs
    rg_REAL nonBondedEnergyBetweenSidechainPairs = 0.0;
    nonBondedEnergyBetweenSidechainPairs = computeVDWESHBEnergyBetweenSidechainPairs(aminoResidues_orderedByChainIDSeqNum, 
                                                                                     numAminoResidues,
                                                                                     engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum);

    rg_REAL nonBondedPotentialEnergy = 0.0;
    nonBondedPotentialEnergy = nonBondedEnergyBetweenSidechainsNBackbone + nonBondedEnergyBetweenSidechainPairs;
    energyOfProtein = nonBondedPotentialEnergy;

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    if(aminoResidues_orderedByChainIDSeqNum != rg_NULL)
        delete [] aminoResidues_orderedByChainIDSeqNum;

    return numAminoResidues;
}



rg_REAL V::GeometryTier::computeVDWESHBEnergyBetweenSidechainsNBackbone(
                                                       Residue** residues,
                                                       const rg_INT& numResidues, 
                                                       Atom** backboneAtoms, 
                                                       const rg_INT& numBackboneAtoms,
                                                       rg_REAL*& engergyWithBackboneOfResidue_orderedByChainIDSeqNum)
{
    Atom* atomsOnBackbone = rg_NULL;
    if(numBackboneAtoms > 0)
        atomsOnBackbone = new Atom[numBackboneAtoms];
    rg_INDEX i;
    for (i = 0;i < numBackboneAtoms;i++)
    {
        atomsOnBackbone[ i ] = *backboneAtoms[ i ];
    }

    //const rg_INT flagsOfEnergyTerms = VDW_TERM_ON | ELEC_TERM_ON | HBOND_TERM_ON;
    PotentialEnergy PE;
    //PE.setFlagsOfAppliedEnergyTerms(flagsOfEnergyTerms);

    rg_REAL energyBetweenSidechainsNBackbone = 0.0;

    for(i = 0;i < numResidues;i++)
    {

        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        rg_REAL VDWEnergyBetweenCurrentSidechainsNBackbone 
            //= PE.computeAndGetPotentialEnergyForTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
            = computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms);
        rg_REAL electrostaticEnergyBetweenCurrentSidechainsNBackbone
            = PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
        //rg_REAL hydrogenBondingEnergyBetweenSidechainPairs
        //    = PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);

        rg_REAL energyBetweenCurrentSidechainNBackbone = VDWEnergyBetweenCurrentSidechainsNBackbone + electrostaticEnergyBetweenCurrentSidechainsNBackbone;

        energyBetweenSidechainsNBackbone += energyBetweenCurrentSidechainNBackbone;

        engergyWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = energyBetweenCurrentSidechainNBackbone;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }

    if(atomsOnBackbone != rg_NULL)
        delete [] atomsOnBackbone;

    return energyBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeVDWESHBEnergyBetweenSidechainPairs(
                                                  Residue** residues,
                                                  const rg_INT& numResidues,
                                                  rg_REAL*& engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum)
{
    //const rg_INT flagsOfEnergyTerms = VDW_TERM_ON | ELEC_TERM_ON | HBOND_TERM_ON;
    PotentialEnergy PE;
    //PE.setFlagsOfAppliedEnergyTerms(flagsOfEnergyTerms);

    rg_REAL energyBetweenSidechainPairs = 0.0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        Atom* atomsOnSidechain_i = rg_NULL;
        if(numSidechainAtoms_i > 0)
            atomsOnSidechain_i = new Atom[numSidechainAtoms_i];
        rg_INDEX k;
        for (k = 0;k < numSidechainAtoms_i;k++)
        {
            atomsOnSidechain_i[ k ] = *sidechainAtoms_i[ k ];
        }

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

            rg_REAL VDWEnergyBetweenCurrentSidechainPairs
                //= PE.computeAndGetPotentialEnergyForTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
                = computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms_j, numSidechainAtoms_j, sidechainAtoms_i, numSidechainAtoms_i);
            rg_REAL electrostaticEnergyBetweenCurrentSidechainPairs
                = PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
            //rg_REAL hydrogenBondingEnergyBetweenCurrentSidechainPairs
            //    = PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
            
            rg_REAL energyBetweenCurrentSidechainPairs = VDWEnergyBetweenCurrentSidechainPairs + electrostaticEnergyBetweenCurrentSidechainPairs;
            energyBetweenSidechainPairs += energyBetweenCurrentSidechainPairs;

            engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += energyBetweenCurrentSidechainPairs / 2.0;
            engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += energyBetweenCurrentSidechainPairs / 2.0;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;

        if(atomsOnSidechain_i != rg_NULL)
            delete [] atomsOnSidechain_i;
    }

    return energyBetweenSidechainPairs;
}



rg_INT  V::GeometryTier::computeXVolumeOfResidueWithBackboneNOtherSidechains(
                                                           Molecule& protein,
                                                           rg_REAL*& xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                           rg_REAL*& xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                           rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                           rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                           rg_REAL&  xvolumeOfProtein,
                                                           rg_INT&   numXOfProtein)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** aminoResidues_orderedByChainIDSeqNum = rg_NULL;
    rg_INT numAminoResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResidues_orderedByChainIDSeqNum);

    xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    numXWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_INT[numAminoResidues];
    numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_INT[numAminoResidues];
    rg_INDEX i;
    for (i = 0;i < numAminoResidues;i++)
    {
        xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum[ i ]        = 0.0;
        xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0.0;
        numXWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = 0;
        numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0;
    }

    // xvolume between sidechains and backbone
    rg_REAL xvolumeBetweenSidechainsNBackbone = 0.0;
    rg_INT numXBetweenSidechainsNBackbone = 0;
    xvolumeBetweenSidechainsNBackbone = computeXVolumeBetweenSidechainsNBackbone(aminoResidues_orderedByChainIDSeqNum, 
                                                                                 numAminoResidues, 
                                                                                 backboneAtoms, 
                                                                                 numBackboneAtoms,
                                                                                 xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                                 numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                                 numXBetweenSidechainsNBackbone);

    // xvolume between sidechain pairs
    rg_REAL xvolumeBetweenSidechainPairs = 0.0;
    rg_INT  numXBetweenSidechainPairs = 0;
    xvolumeBetweenSidechainPairs = computeXVolumeBetweenSidechainPairs(aminoResidues_orderedByChainIDSeqNum, 
                                                                       numAminoResidues,
                                                                       xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                       numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                       numXBetweenSidechainPairs);
        
    xvolumeOfProtein = xvolumeBetweenSidechainsNBackbone + xvolumeBetweenSidechainPairs;
    numXOfProtein = numXBetweenSidechainsNBackbone + numXBetweenSidechainPairs;

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    if(aminoResidues_orderedByChainIDSeqNum != rg_NULL)
        delete [] aminoResidues_orderedByChainIDSeqNum;

    return numAminoResidues;
}



rg_REAL V::GeometryTier::computeXVolumeBetweenSidechainsNBackbone(
                                                 V::GeometryTier::Residue** residues,
                                                 const rg_INT& numResidues, 
                                                 V::GeometryTier::Atom** backboneAtoms,
                                                 const rg_INT& numBackboneAtoms,
                                                 rg_REAL*& xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                 rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                 rg_INT& numXBetweenSidechainsNBackbone)
{
    rg_REAL xvolumeBetweenSidechainsNBackbone = 0.0;
    numXBetweenSidechainsNBackbone = 0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        rg_REAL xvolumeBetweenCurrentSidechainNBackbone = 0.0;
        rg_INT  numXAtoms = 0;
        xvolumeBetweenCurrentSidechainNBackbone 
            = computeSumPairwiseXVolBetweenTwoAtomSets(sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);
        xvolumeBetweenSidechainsNBackbone += xvolumeBetweenCurrentSidechainNBackbone;
        numXBetweenSidechainsNBackbone += numXAtoms;

        xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = xvolumeBetweenCurrentSidechainNBackbone;
        numXWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = numXAtoms;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }    
    
    return xvolumeBetweenSidechainsNBackbone;
}



rg_REAL V::GeometryTier::computeXVolumeBetweenSidechainPairs(
                                            V::GeometryTier::Residue** residues,
                                            const rg_INT& numResidues,
                                            rg_REAL*& xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                            rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                            rg_INT& numXBetweenSidechainPairs)
{
    rg_REAL xvolumeBetweenSidechainPairs = 0.0;
    numXBetweenSidechainPairs = 0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

            rg_REAL xvolumeBetweenCurrentSidechainPairs = 0;
            rg_INT  numXAtoms = 0;
            xvolumeBetweenCurrentSidechainPairs 
                = computeSumPairwiseXVolBetweenTwoAtomSets(sidechainAtoms_j, numSidechainAtoms_j, sidechainAtoms_i, numSidechainAtoms_i, numXAtoms);
            xvolumeBetweenSidechainPairs += xvolumeBetweenCurrentSidechainPairs;
            numXBetweenSidechainPairs += numXAtoms;

            xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += xvolumeBetweenCurrentSidechainPairs / 2.0;
            xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += xvolumeBetweenCurrentSidechainPairs / 2.0;
            numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += numXAtoms / 2.0;
            numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += numXAtoms / 2.0;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;
    }


    return xvolumeBetweenSidechainPairs;
}



rg_INT  V::GeometryTier::computeXAreaOfResidueWithBackboneNOtherSidechains(
                                                         Molecule& protein,
                                                         rg_REAL*& xAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_REAL*& xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                         rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                         rg_REAL& xAreaOfProtein,
                                                         rg_INT&  numXOfProtein)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** aminoResidues_orderedByChainIDSeqNum = rg_NULL;
    rg_INT numAminoResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResidues_orderedByChainIDSeqNum);

    xAreaWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    numXWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_INT[numAminoResidues];
    numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_INT[numAminoResidues];
    rg_INDEX i;
    for (i = 0;i < numAminoResidues;i++)
    {
        xAreaWithBackboneOfResidue_orderedByChainIDSeqNum[ i ]        = 0.0;
        xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0.0;
        numXWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = 0;
        numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0;
    }

    // xvolume between sidechains and backbone
    rg_REAL xAreaBetweenSidechainsNBackbone = 0.0;
    rg_INT numXBetweenSidechainsNBackbone = 0;
    xAreaBetweenSidechainsNBackbone = computeXAreaBetweenSidechainsNBackbone(aminoResidues_orderedByChainIDSeqNum, 
                                                                             numAminoResidues, 
                                                                             backboneAtoms, 
                                                                             numBackboneAtoms,
                                                                             xAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                             numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                             numXBetweenSidechainsNBackbone);

    // xvolume between sidechain pairs
    rg_REAL xAreaBetweenSidechainPairs = 0.0;
    rg_INT  numXBetweenSidechainPairs = 0;
    xAreaBetweenSidechainPairs = computeXAreaBetweenSidechainPairs(aminoResidues_orderedByChainIDSeqNum, 
                                                                     numAminoResidues,
                                                                     xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                     numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                     numXBetweenSidechainPairs);

    xAreaOfProtein = xAreaBetweenSidechainsNBackbone + xAreaBetweenSidechainPairs;
    numXOfProtein = numXBetweenSidechainsNBackbone + numXBetweenSidechainPairs;

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    if(aminoResidues_orderedByChainIDSeqNum != rg_NULL)
        delete [] aminoResidues_orderedByChainIDSeqNum;

    return numAminoResidues;
}



rg_REAL     V::GeometryTier::computeXAreaBetweenSidechainsNBackbone(
                                                V::GeometryTier::Residue** residues,
                                                const rg_INT& numResidues, 
                                                V::GeometryTier::Atom** backboneAtoms,
                                                const rg_INT& numBackboneAtoms,
                                                rg_REAL*& xAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                rg_INT& numXBetweenSidechainsNBackbone)
{
    rg_REAL xAreaBetweenSidechainsNBackbone = 0.0;
    numXBetweenSidechainsNBackbone = 0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        rg_REAL xAreaBetweenCurrentSidechainNBackbone = 0.0;
        rg_INT  numXAtoms = 0;
        xAreaBetweenCurrentSidechainNBackbone 
            = computeSumPairwiseXAreaBetweenTwoAtomSets(sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);
        xAreaBetweenSidechainsNBackbone += xAreaBetweenCurrentSidechainNBackbone;
        numXBetweenSidechainsNBackbone += numXAtoms;

        xAreaWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = xAreaBetweenCurrentSidechainNBackbone;
        numXWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = numXAtoms;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }    

    return xAreaBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeXAreaBetweenSidechainPairs(
                                          V::GeometryTier::Residue** residues,
                                          const rg_INT& numResidues,
                                          rg_REAL*& xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                          rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                          rg_INT& numXBetweenSidechainPairs)
{
    rg_REAL xAreaBetweenSidechainPairs = 0.0;
    numXBetweenSidechainPairs = 0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

            rg_REAL xAreaBetweenCurrentSidechainPairs = 0;
            rg_INT  numXAtoms = 0;
            xAreaBetweenCurrentSidechainPairs 
                = computeSumPairwiseXAreaBetweenTwoAtomSets(sidechainAtoms_j, numSidechainAtoms_j, sidechainAtoms_i, numSidechainAtoms_i, numXAtoms);
            xAreaBetweenSidechainPairs += xAreaBetweenCurrentSidechainPairs;
            numXBetweenSidechainPairs += numXAtoms;

            xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += xAreaBetweenCurrentSidechainPairs / 2.0;
            xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += xAreaBetweenCurrentSidechainPairs / 2.0;
            numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += numXAtoms / 2.0;
            numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += numXAtoms / 2.0;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;
    }


    return xAreaBetweenSidechainPairs;
}




rg_INT  V::GeometryTier::computeXCircleAreaOfResidueWithBackboneNOtherSidechains(
                                                                Molecule& protein,
                                                                rg_REAL*& xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL*& xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL*& xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL*& xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL& sumXCircleAreaOfProtein,
                                                                rg_REAL& sumXCircleRadiusOfProtein,
                                                                rg_INT&  numXOfProtein)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** aminoResidues_orderedByChainIDSeqNum = rg_NULL;
    rg_INT numAminoResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResidues_orderedByChainIDSeqNum);

    xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_REAL[numAminoResidues];
    numXWithBackboneOfResidue_orderedByChainIDSeqNum = new rg_INT[numAminoResidues];
    numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum = new rg_INT[numAminoResidues];
    rg_INDEX i;
    for (i = 0;i < numAminoResidues;i++)
    {
        xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum[ i ]        = 0.0;
        xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0.0;
        xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = 0.0;
        xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0.0;
        numXWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = 0;
        numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] = 0;
    }

    // xCircle between sidechains and backbone
    rg_REAL sumxCircleAreaBetweenSidechainsNBackbone = 0.0;
    rg_REAL sumXCircleRadiusBetweenSidechainsNBackbone = 0.0;
    rg_INT numXBetweenSidechainsNBackbone = 0;
    sumxCircleAreaBetweenSidechainsNBackbone = computeXCircleAreaBetweenSidechainsNBackbone(aminoResidues_orderedByChainIDSeqNum, 
                                                                                   numAminoResidues, 
                                                                                   backboneAtoms, 
                                                                                   numBackboneAtoms,
                                                                                   xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                                   xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                                   numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                                   sumXCircleRadiusBetweenSidechainsNBackbone,
                                                                                   numXBetweenSidechainsNBackbone);

    // xCircle between sidechain pairs
    rg_REAL sumXCircleAreaBetweenSidechainPairs = 0.0;
    rg_REAL sumXCircleRadiusBetweenSidechainPairs = 0.0;
    rg_INT  numXBetweenSidechainPairs = 0;
    sumXCircleAreaBetweenSidechainPairs = computeXCircleAreaBetweenSidechainPairs(aminoResidues_orderedByChainIDSeqNum, 
                                                                               numAminoResidues,
                                                                               xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                               xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                               numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                               sumXCircleRadiusBetweenSidechainPairs,
                                                                               numXBetweenSidechainPairs);

    sumXCircleAreaOfProtein = sumxCircleAreaBetweenSidechainsNBackbone + sumXCircleAreaBetweenSidechainPairs;
    sumXCircleRadiusOfProtein = sumXCircleRadiusBetweenSidechainsNBackbone + sumXCircleRadiusBetweenSidechainPairs;

    numXOfProtein = numXBetweenSidechainsNBackbone + numXBetweenSidechainPairs;

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    if(aminoResidues_orderedByChainIDSeqNum != rg_NULL)
        delete [] aminoResidues_orderedByChainIDSeqNum;

    return numAminoResidues;
}



rg_REAL     V::GeometryTier::computeXCircleAreaBetweenSidechainsNBackbone(
                                                        V::GeometryTier::Residue** residues,
                                                        const rg_INT& numResidues, 
                                                        V::GeometryTier::Atom** backboneAtoms,
                                                        const rg_INT& numBackboneAtoms,
                                                        rg_REAL*& xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                        rg_REAL*& xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                        rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                        rg_REAL& sumXCircleRadiusBetweenSidechainsNBackbone,
                                                        rg_INT& numXBetweenSidechainsNBackbone)
{
    rg_REAL sumXCircleAreaBetweenSidechainsNBackbone = 0.0;
    sumXCircleRadiusBetweenSidechainsNBackbone = 0.0;
    numXBetweenSidechainsNBackbone = 0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        rg_REAL xCircleAreaBetweenCurrentSidechainNBackbone = 0.0;
        rg_REAL xCircleRadiusBetweenCurrentSidechainNBackbone = 0.0;
        rg_INT  numXAtoms = 0;
        
        xCircleAreaBetweenCurrentSidechainNBackbone 
            = computeSumPairwiseXCircleAreaBetweenTwoAtomSets(sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms, xCircleRadiusBetweenCurrentSidechainNBackbone, numXAtoms);

        sumXCircleAreaBetweenSidechainsNBackbone += xCircleAreaBetweenCurrentSidechainNBackbone;
        xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = xCircleAreaBetweenCurrentSidechainNBackbone;
        xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = xCircleRadiusBetweenCurrentSidechainNBackbone;
        numXWithBackboneOfResidue_orderedByChainIDSeqNum[ i ] = numXAtoms;
        sumXCircleRadiusBetweenSidechainsNBackbone += xCircleRadiusBetweenCurrentSidechainNBackbone;
        numXBetweenSidechainsNBackbone += numXAtoms;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }    

    return sumXCircleAreaBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeXCircleAreaBetweenSidechainPairs(
                                                        V::GeometryTier::Residue** residues,
                                                        const rg_INT& numResidues,
                                                        rg_REAL*& xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                        rg_REAL*& xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                        rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                        rg_REAL& sumXCircleRadiusBetweenSidechainPairs,
                                                        rg_INT& numXBetweenSidechainPairs)
{
    rg_REAL sumXCircleAreaBetweenSidechainPairs = 0.0;
    sumXCircleRadiusBetweenSidechainPairs = 0.0;
    numXBetweenSidechainPairs = 0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

            rg_REAL xCircleAreaBetweenCurrentSidechainPairs = 0;
            rg_REAL xCircleRadiusBetweenCurrentSidechainPairs = 0;
            rg_INT  numXAtoms = 0;

            xCircleAreaBetweenCurrentSidechainPairs 
                = computeSumPairwiseXCircleAreaBetweenTwoAtomSets(sidechainAtoms_j, numSidechainAtoms_j, sidechainAtoms_i, numSidechainAtoms_i, xCircleRadiusBetweenCurrentSidechainPairs, numXAtoms);
            
            sumXCircleAreaBetweenSidechainPairs += xCircleAreaBetweenCurrentSidechainPairs;
            
            xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += xCircleAreaBetweenCurrentSidechainPairs / 2.0;
            xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += xCircleAreaBetweenCurrentSidechainPairs / 2.0;

            numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ i ] += numXAtoms / 2.0;
            numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum[ j ] += numXAtoms / 2.0;

            sumXCircleRadiusBetweenSidechainPairs += xCircleRadiusBetweenCurrentSidechainPairs;

            numXBetweenSidechainPairs += numXAtoms;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;
    }


    return sumXCircleAreaBetweenSidechainPairs;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyOfProtein(Molecule& protein)
{		
	// get backbone atoms
	Atom** backboneAtoms = rg_NULL;
	rg_INT numBackboneAtoms
		= protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

	// get target residues
	Residue** sortedResidues = rg_NULL;
	rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

	rg_REAL energyOfSidechainNBackbone = 0;
	energyOfSidechainNBackbone = computeVanDerWaalsEnergyBetweenSidechainsNBackbone(sortedResidues, numResidues, backboneAtoms, numBackboneAtoms);

	if(backboneAtoms != rg_NULL)
		delete[] backboneAtoms;

	rg_REAL energyOfSidechainPairs = 0;
	energyOfSidechainPairs = computeVanDerWaalsEnergyBetweenSidechainPairs(sortedResidues, numResidues);

	if(sortedResidues != rg_NULL)
		delete [] sortedResidues;

	rg_REAL energy = 0.0;
	energy = energyOfSidechainNBackbone + energyOfSidechainPairs;

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainsNBackbone(V::GeometryTier::Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms)
{
    rg_REAL energy = 0.0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;
        //numSidechainAtoms = residues[ i ]->getAtomsOnSidechain(sidechainAtoms);
        //This function should be called when we construct the energy of protein structure computed from other softwares
        //Hence, we have to use "FunctionsForRotamerLibrary::getAtomsOnSidechain()".
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        // test
        //if(residues[ i ]->getChain()->getChainIDFromInputFileInDecimal() == 'A' && residues[ i ]->getSequenceNumber() == 29)
        //	int here = 1;

        rg_REAL currEnergy = 0;
        currEnergy = computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms);

        // test
        //if(currEnergy > 100)
        //	cout << "Big energy: " << currEnergy << " " << residues[ i ]->getChain()->getChainIDFromInputFileInString() << " " << residues[ i ]->getSequenceNumber() << " " << residues[ i ]->getResidueName() << endl;

        energy += currEnergy;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }

    return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainPairs(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues)
{
    rg_REAL energy = 0.0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        //numSidechainAtoms_i = residues[ i ]->getAtomsOnSidechain(sidechainAtoms_i);
        //This function should be called when we construct the energy of protein structure computed from other softwares
        //Hence, we have to use "FunctionsForRotamerLibrary::getAtomsOnSidechain()".
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            //numSidechainAtoms_j = residues[ j ]->getAtomsOnSidechain(sidechainAtoms_j);
            //This function should be called when we construct the energy of protein structure computed from other softwares
            //Hence, we have to use "FunctionsForRotamerLibrary::getAtomsOnSidechain()".
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);


            // test
            //if(residues[ i ]->getChain()->getChainIDFromInputFileInDecimal() == 'A' && residues[ i ]->getSequenceNumber() == 29)
            //	int here = 1;

            rg_REAL currEnergy = 0;
            currEnergy = computeVanDerWaalsEnergyBetweenTwoAtomSets(sidechainAtoms_i, numSidechainAtoms_i, sidechainAtoms_j, numSidechainAtoms_j);

            // test
            //if(currEnergy > 100)
            //	cout << "Big energy: " << currEnergy << " " << residues[ i ]->getChain()->getChainIDFromInputFileInString() << " " << residues[ i ]->getSequenceNumber() << " " << residues[ i ]->getResidueName() << " "
            //	     << residues[ j ]->getChain()->getChainIDFromInputFileInString() << " " << residues[ j ]->getSequenceNumber() << " " << residues[ j ]->getResidueName() << endl;

            energy += currEnergy;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;
    }

    return energy;
}



rg_REAL     V::GeometryTier::computeElectrostaticEnergyOfProtein(Molecule& protein)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** sortedResidues = rg_NULL;
    rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

    rg_REAL energyOfSidechainNBackbone = 0;
    energyOfSidechainNBackbone = computeElectrostaticEnergyBetweenSidechainsNBackbone(sortedResidues, numResidues, backboneAtoms, numBackboneAtoms);

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    rg_REAL energyOfSidechainPairs = 0;
    energyOfSidechainPairs = computeElectrostaticEnergyBetweenSidechainPairs(sortedResidues, numResidues);

    if(sortedResidues != rg_NULL)
        delete [] sortedResidues;

    rg_REAL energy = 0.0;
    energy = energyOfSidechainNBackbone + energyOfSidechainPairs;

    return energy;
}



rg_REAL     V::GeometryTier::computeElectrostaticEnergyBetweenSidechainsNBackbone(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues, 
                                    V::GeometryTier::Atom** backboneAtoms, 
                                    const rg_INT& numBackboneAtoms)
{
    Atom* atomsOnBackbone = rg_NULL;
    atomsOnBackbone = new Atom[numBackboneAtoms];
    rg_INDEX i;
    for (i = 0;i < numBackboneAtoms;i++)
    {
        atomsOnBackbone[ i ] = *backboneAtoms[ i ];
    }

    PotentialEnergy PE;
    rg_REAL energyBetweenSidechainsNBackbone = 0.0;

    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        rg_REAL energyBetweenCurrentSidechainNBackbone = 0.0;
        energyBetweenCurrentSidechainNBackbone 
            = PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
        energyBetweenSidechainsNBackbone += energyBetweenCurrentSidechainNBackbone;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }

    if(atomsOnBackbone != rg_NULL)
        delete [] atomsOnBackbone;

    return energyBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeElectrostaticEnergyBetweenSidechainPairs(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues)
{
    PotentialEnergy PE;

    rg_REAL energyBetweenSidechainPairs = 0.0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        Atom* atomsOnSidechain_i = rg_NULL;
        atomsOnSidechain_i = new Atom[numSidechainAtoms_i];
        rg_INDEX k;
        for (k = 0;k < numSidechainAtoms_i;k++)
        {
            atomsOnSidechain_i[ k ] = *sidechainAtoms_i[ k ];
        }

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

            rg_REAL energyBetweenCurrentSidechainPairs = 0;
            energyBetweenCurrentSidechainPairs 
                = PE.computeElectrostaticBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
            energyBetweenSidechainPairs += energyBetweenCurrentSidechainPairs;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;

        if(atomsOnSidechain_i != rg_NULL)
            delete [] atomsOnSidechain_i;
    }

    return energyBetweenSidechainPairs;
}



rg_REAL     V::GeometryTier::computeHydrogenBondingEnergyOfProtein(Molecule& protein)
{
    // get backbone atoms
    Atom** backboneAtoms = rg_NULL;
    rg_INT numBackboneAtoms
        = protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

    // get target residues
    Residue** sortedResidues = rg_NULL;
    rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

    rg_REAL energyOfSidechainNBackbone = 0;
    energyOfSidechainNBackbone = computeHydrogenBondingEnergyBetweenSidechainsNBackbone(sortedResidues, numResidues, backboneAtoms, numBackboneAtoms);

    if(backboneAtoms != rg_NULL)
        delete[] backboneAtoms;

    rg_REAL energyOfSidechainPairs = 0;
    energyOfSidechainPairs = computeHydrogenBondingEnergyBetweenSidechainPairs(sortedResidues, numResidues);

    if(sortedResidues != rg_NULL)
        delete [] sortedResidues;

    rg_REAL energy = 0.0;
    energy = energyOfSidechainNBackbone + energyOfSidechainPairs;

    return energy;
}



rg_REAL     V::GeometryTier::computeHydrogenBondingEnergyBetweenSidechainsNBackbone(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues, 
                                    V::GeometryTier::Atom** backboneAtoms, 
                                    const rg_INT& numBackboneAtoms)
{
    Atom* atomsOnBackbone = rg_NULL;
    atomsOnBackbone = new Atom[numBackboneAtoms];
    rg_INDEX i;
    for (i = 0;i < numBackboneAtoms;i++)
    {
        atomsOnBackbone[ i ] = *backboneAtoms[ i ];
    }

    PotentialEnergy PE;
    rg_REAL energyBetweenSidechainsNBackbone = 0.0;

    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms = 0;
        Atom** sidechainAtoms = rg_NULL;		
        numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

        rg_REAL energyBetweenCurrentSidechainNBackbone = 0.0;
        energyBetweenCurrentSidechainNBackbone 
            = PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms, numSidechainAtoms, atomsOnBackbone, numBackboneAtoms);
        energyBetweenSidechainsNBackbone += energyBetweenCurrentSidechainNBackbone;

        if(sidechainAtoms != rg_NULL)
            delete [] sidechainAtoms;
    }

    if(atomsOnBackbone != rg_NULL)
        delete [] atomsOnBackbone;

    return energyBetweenSidechainsNBackbone;
}



rg_REAL     V::GeometryTier::computeHydrogenBondingEnergyBetweenSidechainPairs(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues)
{
    PotentialEnergy PE;

    rg_REAL energyBetweenSidechainPairs = 0.0;

    rg_INDEX i;
    for(i = 0;i < numResidues;i++)
    {
        rg_INT numSidechainAtoms_i = 0;
        Atom** sidechainAtoms_i = rg_NULL;
        numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

        Atom* atomsOnSidechain_i = rg_NULL;
        atomsOnSidechain_i = new Atom[numSidechainAtoms_i];
        rg_INDEX k;
        for (k = 0;k < numSidechainAtoms_i;k++)
        {
            atomsOnSidechain_i[ k ] = *sidechainAtoms_i[ k ];
        }

        rg_INDEX j;
        for(j = i + 1;j < numResidues;j++)
        {
            rg_INT numSidechainAtoms_j = 0;
            Atom** sidechainAtoms_j = rg_NULL;
            numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

            rg_REAL energyBetweenCurrentSidechainPairs = 0;
            energyBetweenCurrentSidechainPairs 
                = PE.computeHydrogenBondBetweenTwoMolecules(sidechainAtoms_j, numSidechainAtoms_j, atomsOnSidechain_i, numSidechainAtoms_i);
            energyBetweenSidechainPairs += energyBetweenCurrentSidechainPairs;

            if(sidechainAtoms_j != rg_NULL)
                delete [] sidechainAtoms_j;
        }

        if(sidechainAtoms_i != rg_NULL)
            delete [] sidechainAtoms_i;

        if(atomsOnSidechain_i != rg_NULL)
            delete [] atomsOnSidechain_i;
    }

    return energyBetweenSidechainPairs;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyOfProteinBySplineFitting(Molecule& protein)
{
	// get backbone atoms
	Atom** backboneAtoms = rg_NULL;
	rg_INT numBackboneAtoms
		= protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

	// get target residues
	Residue** sortedResidues = rg_NULL;
	rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

	rg_REAL energyOfSidechainNBackbone = 0;
	energyOfSidechainNBackbone = computeVanDerWaalsEnergyBetweenSidechainsNBackboneBySplineFitting(sortedResidues, numResidues, backboneAtoms, numBackboneAtoms);

	if(backboneAtoms != rg_NULL)
		delete[] backboneAtoms;

	rg_REAL energyOfSidechainPairs = 0;
	energyOfSidechainPairs = computeVanDerWaalsEnergyBetweenSidechainPairsBySplineFitting(sortedResidues, numResidues);

	if(sortedResidues != rg_NULL)
		delete [] sortedResidues;

	rg_REAL energy = 0.0;
	energy = energyOfSidechainNBackbone + energyOfSidechainPairs;

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainsNBackboneBySplineFitting(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues, 
                                    V::GeometryTier::Atom** backboneAtoms, 
                                    const rg_INT& numBackboneAtoms)
{
	rg_REAL energy = 0.0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_INT numSidechainAtoms = 0;
		Atom** sidechainAtoms = rg_NULL;

		numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

		rg_REAL currEnergy = 0;
		currEnergy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms);

		energy += currEnergy;

		if(sidechainAtoms != rg_NULL)
			delete [] sidechainAtoms;
	}

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainNBackbone(
                                    V::GeometryTier::Residue* residue, 
                                    V::GeometryTier::Atom** backboneAtoms, 
                                    const rg_INT& numBackboneAtoms)
{
	rg_INT numSidechainAtoms = 0;
	Atom** sidechainAtoms = rg_NULL;
	//numSidechainAtoms = sortedResidues[ i ]->getAtomsOnSidechain(sidechainAtoms);
	numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residue, sidechainAtoms);

	rg_REAL energy = 0.0;
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms);

	if(sidechainAtoms != rg_NULL)
		delete [] sidechainAtoms;

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainPairsBySplineFitting(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues)
{
	rg_REAL energy = 0.0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_INT numSidechainAtoms_i = 0;
		Atom** sidechainAtoms_i = rg_NULL;

		numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

		rg_INDEX j;
		for(j = i + 1;j < numResidues;j++)
		{
			rg_INT numSidechainAtoms_j = 0;
			Atom** sidechainAtoms_j = rg_NULL;

			numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

			rg_REAL currEnergy = 0;
			currEnergy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting(sidechainAtoms_i, numSidechainAtoms_i, sidechainAtoms_j, numSidechainAtoms_j);

			energy += currEnergy;

			if(sidechainAtoms_j != rg_NULL)
				delete [] sidechainAtoms_j;
		}

		if(sidechainAtoms_i != rg_NULL)
			delete [] sidechainAtoms_i;
	}

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenTwoSidechains(
                                    V::GeometryTier::Residue* residue1, 
                                    V::GeometryTier::Residue* residue2)
{
	Atom** sidechainAtomsOfResidue1 = rg_NULL;
	rg_INT numAtoms1 = 0;
	numAtoms1 = FunctionsForRotamerLibrary::getAtomsOnSidechain(residue1, sidechainAtomsOfResidue1);
	Atom** sidechainAtomsOfResidue2 = rg_NULL;
	rg_INT numAtoms2 = 0;
	numAtoms2 = FunctionsForRotamerLibrary::getAtomsOnSidechain(residue2, sidechainAtomsOfResidue2);

	rg_REAL energy = 0.0;
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(sidechainAtomsOfResidue1, numAtoms1, sidechainAtomsOfResidue2, numAtoms2);

	if(sidechainAtomsOfResidue1 != rg_NULL)
		delete [] sidechainAtomsOfResidue1;

	if(sidechainAtomsOfResidue2 != rg_NULL)
		delete [] sidechainAtomsOfResidue2;

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenTwoResidues(
                                    V::GeometryTier::Residue* residue1, 
                                    V::GeometryTier::Residue* residue2)
{
	Atom** atomsOfResidue1 = rg_NULL;
	rg_INT numAtoms1 = 0;
	numAtoms1 = FunctionsForRotamerLibrary::getAtomsOfResidue(residue1, atomsOfResidue1);
	Atom** atomsOfResidue2 = rg_NULL;
	rg_INT numAtoms2 = 0;
	numAtoms2 = FunctionsForRotamerLibrary::getAtomsOfResidue(residue2, atomsOfResidue2);

	rg_REAL energy = 0.0;
	energy = computeVanDerWaalsEnergyBetweenTwoAtomSets(atomsOfResidue1, numAtoms1, atomsOfResidue2, numAtoms2);

	if(atomsOfResidue1 != rg_NULL)
		delete [] atomsOfResidue1;

	if(atomsOfResidue2 != rg_NULL)
		delete [] atomsOfResidue2;

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenTwoAtomSets( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2 )
{
	rg_REAL totalEnergy = 0.0;

	for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
	{
		rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
		rg_INT atomTypeInAmberForAtom1  = (rg_INT)atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();

		for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
		{
			if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
				continue;

			rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
			rg_INT atomTypeInAmberForAtom2 = (rg_INT)atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();

			rg_REAL squaredDistance = centerOfAtom2.squaredDistance(centerOfAtom1);
			rg_REAL sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

			rg_REAL energy = 0.0;
			energy = (   ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][0] / (sixthPoweredDistance*sixthPoweredDistance) )
				       - ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][1] /  sixthPoweredDistance )                       );
			
			totalEnergy += energy;
		}
	}

	return totalEnergy; 
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2 )
{	
	// test
	//clock_t startTime = clock();

	rg_REAL totalEnergy = 0.0;

	for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
	{
		rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
		AmberAtomTypes atomTypeInAmberForAtom1 = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();

		for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
		{
			if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
				continue;

			rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
			AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();

			rg_REAL energy = 0.0;
			//rg_REAL distBTWCenters = centerOfAtom2.distance(centerOfAtom1);			
			//energy = FunctionsForSplineFittingLJPotEnergy::evaluateLJPotEnergyBySplineFitting(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters);
			//Currently, we use squared distance instead of distance
			rg_REAL squaredDistance = centerOfAtom2.squaredDistance(centerOfAtom1);
			energy = FunctionsForSplineFittingLJPotEnergy::evaluateLJPotEnergyBySplineFitting(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, squaredDistance);

			// test
			//rg_REAL sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;
			//rg_REAL LJEnergy = 0.0;
			//LJEnergy = (   ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][0] / (sixthPoweredDistance*sixthPoweredDistance) )
			//	- ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][1] /  sixthPoweredDistance )                       );
			//if(fabs(LJEnergy - energy) > 100)
			//{
			//	cout << centerOfAtom2.distance(centerOfAtom1)   << " " 
			//		<< squaredDistance                         << " "				       
			//		<< LJEnergy                                << " " 
			//		<< energy                                  << " " 
			//		<< LJEnergy - energy                       << " " 
			//		<< atomSet1[i]->getAtomNameFromInputFile() << " " 
			//		<< atomSet2[j]->getAtomNameFromInputFile() << " " 
			//		<< AmberAtomTypeStrings[atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber()] << " "
			//		<< AmberAtomTypeStrings[atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber()] << " "
			//		<< endl;

			//	energy = FunctionsForSplineFittingLJPotEnergy::evaluateLJPotEnergyBySplineFitting(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, squaredDistance);
			//}

			totalEnergy += energy;
		}
	}

	// test
	//clock_t endTime = clock();
	//rg_REAL computationTime = (rg_REAL)(endTime - startTime) / 1000000;
	//FunctionsForSplineFittingLJPotEnergy::m_time4EvalSplineOfLJEnergy += computationTime;

	return totalEnergy; 

}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenTwoAtoms( 
                                    V::GeometryTier::Atom* atom1, 
                                    V::GeometryTier::Atom* atom2 )
{
	rg_REAL totalEnergy = 0.0;
			
	if(atom1->hasChemicalBondWith( atom2 ))
		return totalEnergy;

	rg_Point3D centerOfAtom1 = atom1->getpAtomBall()->getCenter();
	rg_INT atomTypeInAmberForAtom1  = -1;
	atomTypeInAmberForAtom1 = (rg_INT)atom1->getpChemicalProperties()->getAtomTypeInAmber();

	rg_Point3D centerOfAtom2 = atom2->getpAtomBall()->getCenter();
	rg_INT atomTypeInAmberForAtom2 = -1;
	atomTypeInAmberForAtom2 = (rg_INT)atom2->getpChemicalProperties()->getAtomTypeInAmber();

	rg_REAL squaredDistance = 0.0;
	squaredDistance = ( centerOfAtom2 - centerOfAtom1 ).squaredMagnitude();

	rg_REAL sixthPoweredDistance = 0.0;
	sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

	totalEnergy = (   ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][0] / (sixthPoweredDistance*sixthPoweredDistance) )
		            - ( L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][1] /  sixthPoweredDistance )                       );

	return totalEnergy;
}



Vector3D    V::GeometryTier::evaluateGradientVectorForVanDerWaalsEnergyBetweenTwoAtomsWRTCoordinate(
                                    V::GeometryTier::Atom* atom1, 
                                    V::GeometryTier::Atom* atom2)
{
	//if (atom1->hasChemicalBondWith( atom2 ))
	//{
	//	Vector3D gradientVec(0.0, 0.0, 0.0);
	//	return gradientVec;
	//}
	//else
	{
		rg_Point3D center1 = atom1->getAtomBall().getCenter();
		rg_Point3D center2 = atom2->getAtomBall().getCenter();
		rg_REAL dist                   = center1.distance( center2 );
		rg_REAL squaredDist            = dist * dist;
		rg_REAL quarticPoweredDist     = squaredDist * squaredDist;
		rg_REAL seventhPoweredDist     = quarticPoweredDist * squaredDist * dist;
		rg_REAL threeteenthPoweredDist = seventhPoweredDist * quarticPoweredDist * squaredDist;
		
		rg_INT atomTypeInAmber1 = (rg_INT)atom1->getpChemicalProperties()->getAtomTypeInAmber();
		rg_INT atomTypeInAmber2 = (rg_INT)atom2->getpChemicalProperties()->getAtomTypeInAmber();
		
		rg_REAL derivativeForDist = 0.0;
		derivativeForDist = 
			-12.0 * L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmber1][atomTypeInAmber2][0] / threeteenthPoweredDist
			+ 6.0 * L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmber1][atomTypeInAmber2][1] / seventhPoweredDist    ;

		rg_REAL x1, y1, z1;
		center1.getCoordinates(x1, y1, z1);
		rg_REAL x2, y2, z2;
		center2.getCoordinates(x2, y2, z2);

		Vector3D gradientVec( derivativeForDist * (x1 - x2) / dist, derivativeForDist * (y1 - y2) / dist, derivativeForDist * (z1 - z2) / dist);
		return gradientVec;
	}
}



Vector3D    V::GeometryTier::evaluateGradientVectorForNewScoreFunctionBetweenTwoAtomsWRTCoordinate(
                                    V::GeometryTier::Atom* atom1, 
                                    V::GeometryTier::Atom* atom2)
{
    if (atom1->hasChemicalBondWith( atom2 ))
    {
        Vector3D gradientVec(0.0, 0.0, 0.0);
        return gradientVec;
    }
    else
    {
        rg_Point3D center1        = atom1->getAtomBall().getCenter();
        rg_Point3D center2        = atom2->getAtomBall().getCenter();
        rg_REAL    distBTWCenters = center1.distance( center2 );
        rg_REAL    radius1        = atom1->getAtomBall().getRadius();
        rg_REAL    radius2        = atom2->getAtomBall().getRadius();

        // inflate radii 10%
        radius1 *= 1.1;
        radius2 *= 1.1;

        rg_REAL distFromBoundary  = distBTWCenters - (radius1 + radius2);
        rg_REAL derivativeForDist = 0.0;

        if(rg_LT(distFromBoundary, 0.0))
        {
            rg_REAL sumOfTwoRadii  = radius1 + radius2;
            rg_REAL distDiff = fabs(radius1 - radius2);
            if( (distDiff < distBTWCenters) && (distBTWCenters < sumOfTwoRadii) )
            {
                rg_REAL dist_squared    = distBTWCenters * distBTWCenters;
                rg_REAL radius1_squared = radius1 * radius1;
                rg_REAL radius2_squared = radius2 * radius2;
                rg_REAL radius1_squared_minus_radius2_squared = radius1_squared - radius2_squared;

                derivativeForDist = rg_PI * dist_squared / 4.0 - rg_PI * (radius1_squared + radius2_squared) / 2.0 + rg_PI * radius1_squared_minus_radius2_squared * radius1_squared_minus_radius2_squared / (4.0 * dist_squared);
            }
            else
            {
                derivativeForDist = 0.0;
            }
        }
        else
        {
            derivativeForDist = 1.0;
        }

        rg_REAL x1, y1, z1;
        center1.getCoordinates(x1, y1, z1);
        rg_REAL x2, y2, z2;
        center2.getCoordinates(x2, y2, z2);

        Vector3D gradientVec( derivativeForDist * (x1 - x2) / distBTWCenters, derivativeForDist * (y1 - y2) / distBTWCenters, derivativeForDist * (z1 - z2) / distBTWCenters);
        return gradientVec;
    }
}



Vector3D    V::GeometryTier::evaluateGradientVectorForPositionBetweenTwoAtomsBySplineFitting(
                                    V::GeometryTier::Atom* atom1, 
                                    V::GeometryTier::Atom* atom2)
{
	if (atom1->hasChemicalBondWith( atom2 ))
	{
		Vector3D gradientVec(0.0, 0.0, 0.0);
		return gradientVec;
	}
	else
	{
		//rg_INT atomTypeInAmber1 = (rg_INT)atom1->getpChemicalProperties()->getAtomTypeInAmber();
		//rg_INT atomTypeInAmber2 = (rg_INT)atom2->getpChemicalProperties()->getAtomTypeInAmber();

		//rg_REAL squaredDist            = dist * dist;
		//rg_REAL quarticPoweredDist     = squaredDist * squaredDist;
		//rg_REAL seventhPoweredDist     = quarticPoweredDist * squaredDist * dist;
		//rg_REAL threeteenthPoweredDist = seventhPoweredDist * quarticPoweredDist * squaredDist;

		rg_REAL derivativeForDist = 0.0;
		AmberAtomTypes atomTypeInAmber1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
		AmberAtomTypes atomTypeInAmber2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();
		rg_Point3D center1 = atom1->getAtomBall().getCenter();
		rg_Point3D center2 = atom2->getAtomBall().getCenter();
		//rg_REAL dist                   = center1.distance( center2 );
		//derivativeForDist  = FunctionsForSplineFittingLJPotEnergy::evaluateDerivativeOfLJPotEnergyBySplineFitting(atomTypeInAmber1, atomTypeInAmber2, dist);
		//Currently, we use squared distance instead of distance
		rg_REAL squaredDistance = center1.squaredDistance(center2);
		derivativeForDist  = FunctionsForSplineFittingLJPotEnergy::evaluateDerivativeOfLJPotEnergyBySplineFitting(atomTypeInAmber1, atomTypeInAmber2, squaredDistance);

		//derivativeForDist = computeDerivativeOfVanDerWaalsEnergyBetweenTwoAtoms(atom1, atom2, dist);
		//derivativeForDist = 
		//	-12.0 * L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmber1][atomTypeInAmber2][0] / threeteenthPoweredDist
		//	+ 6.0 * L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmber1][atomTypeInAmber2][1] / seventhPoweredDist    ;

		rg_REAL dist = sqrt( squaredDistance );
		rg_REAL x1, y1, z1;
		center1.getCoordinates(x1, y1, z1);
		rg_REAL x2, y2, z2;
		center2.getCoordinates(x2, y2, z2);		

		Vector3D gradientVec( derivativeForDist * (x1 - x2) / dist, derivativeForDist * (y1 - y2) / dist, derivativeForDist * (z1 - z2) / dist);
		return gradientVec;
	}
}



rg_REAL     V::GeometryTier::computeDerivativeOfVanDerWaalsEnergyBetweenTwoAtoms(
                                    V::GeometryTier::Atom* atom1, 
                                    V::GeometryTier::Atom* atom2, 
                                    const rg_REAL& distCenters)
{
	AmberAtomTypes atomTypeInAmber1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomTypeInAmber2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();

	rg_REAL derivativeOfLJPotEnergy = 0.0;
	derivativeOfLJPotEnergy = FunctionsForSplineFittingLJPotEnergy::evaluateDerivativeOfLJPotEnergyBySplineFitting(atomTypeInAmber1, atomTypeInAmber2, distCenters);

	return derivativeOfLJPotEnergy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBySCWRL3(Molecule& protein)
{
	// get backbone atoms
	Atom** backboneAtoms = rg_NULL;
	rg_INT numBackboneAtoms
		= protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

	// get target residues
	Residue** sortedResidues = rg_NULL;
	rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues );

	rg_REAL energyBTWSidechainNBackbone = 0.0;
	energyBTWSidechainNBackbone = computeVanDerWaalsEnergyBetweenSidechainNBackboneBySCWRL3( sortedResidues, numResidues, backboneAtoms, numBackboneAtoms );

	if(backboneAtoms != rg_NULL)
		delete[] backboneAtoms;

	rg_REAL energyBTWSidechains = 0.0;
	energyBTWSidechains = computeVanDerWaalsEnergyBetweenSidechainsBySCWRL3( sortedResidues, numResidues );

	if(sortedResidues != rg_NULL)
		delete [] sortedResidues;

	rg_REAL energy = 0.0;
	energy = energyBTWSidechainNBackbone + energyBTWSidechains;

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainNBackboneBySCWRL3(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues, 
                                    V::GeometryTier::Atom** backboneAtoms, 
                                    const rg_INT& numBackboneAtoms)
{
	rg_REAL energy = 0.0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_INT numSidechainAtoms = 0;
		Atom** sidechainAtoms = rg_NULL;
		//numSidechainAtoms = sortedResidues[ i ]->getAtomsOnSidechain(sidechainAtoms);
		numSidechainAtoms = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms);

		rg_REAL currEnergy = 0;
        rg_INT numXAtoms = 0;
		currEnergy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3( sidechainAtoms, numSidechainAtoms, backboneAtoms, numBackboneAtoms, numXAtoms);
		energy += currEnergy;

		if(sidechainAtoms != rg_NULL)
			delete [] sidechainAtoms;
	}

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenSidechainsBySCWRL3(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues)
{
	rg_REAL energy = 0.0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_INT numSidechainAtoms_i = 0;
		Atom** sidechainAtoms_i = rg_NULL;
		//numSidechainAtoms_i = sortedResidues[ i ]->getAtomsOnSidechain(sidechainAtoms_i);
		numSidechainAtoms_i = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ i ], sidechainAtoms_i);

		rg_INDEX j;
		for(j = i + 1;j < numResidues;j++)
		{
			rg_INT numSidechainAtoms_j = 0;
			Atom** sidechainAtoms_j = rg_NULL;
			//numSidechainAtoms_j = sortedResidues[ j ]->getAtomsOnSidechain(sidechainAtoms_j);
			numSidechainAtoms_j = FunctionsForRotamerLibrary::getAtomsOnSidechain(residues[ j ], sidechainAtoms_j);

			rg_REAL currEnergy = 0;
            rg_INT numXAtoms = 0;
			currEnergy = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3(sidechainAtoms_i, numSidechainAtoms_i, sidechainAtoms_j, numSidechainAtoms_j, numXAtoms);
			energy += currEnergy;

			if(sidechainAtoms_j != rg_NULL)
				delete [] sidechainAtoms_j;
		}

		if(sidechainAtoms_i != rg_NULL)
			delete [] sidechainAtoms_i;
	}

	return energy;
}



rg_REAL     V::GeometryTier::computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_INT& numXAtoms )
{
    numXAtoms = 0;
	rg_REAL totalEnergy = 0.0;
	const rg_REAL k_rep = 5.882;
	const rg_REAL k_att = 0.08;

	for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
	{
		rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
		rg_REAL    radiusOfAtom1 = atomSet1[i]->getpAtomBall()->getRadius();

		for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
		{
			if (atomSet1[i]->hasChemicalBondWith( atomSet2[j]))
			{
				continue;
			}

			rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
			rg_REAL    radiusOfAtom2 = atomSet2[j]->getpAtomBall()->getRadius();

			rg_REAL distance = 0.0;
			distance = centerOfAtom1.distance( centerOfAtom2 );
			rg_REAL sumOfRadii = radiusOfAtom1 + radiusOfAtom2;
			rg_REAL twiceSumOfRadii = 2.0*sumOfRadii;

			rg_REAL energy = 0.0;

			if( rg_LE(distance, sumOfRadii) )
			{
				energy = k_rep * ( 1 - distance / sumOfRadii);
                numXAtoms++;
			}
			else if( rg_GT(distance, sumOfRadii) && rg_LT(distance, twiceSumOfRadii))
			{
				rg_REAL distDividedBySumOfRadii = distance / sumOfRadii;
				energy = k_att* ( distDividedBySumOfRadii*distDividedBySumOfRadii - 3.0 * distDividedBySumOfRadii + 2.0);
			}
			else
			{}

			totalEnergy += energy;
		}
	}

	return totalEnergy; 
}



Vector3D    V::GeometryTier::evaluateGradientVectorForPositionBetweenTwoAtomsBySCWRL3(
                                    V::GeometryTier::Atom* atom1, 
                                    V::GeometryTier::Atom* atom2)
{
	if (atom1->hasChemicalBondWith( atom2 ))
	{
		Vector3D gradientVec;
		gradientVec.setPoint(0.0, 0.0, 0.0);
		return gradientVec;
	}
	else
	{
		const rg_REAL k_rep             = 5.882;
		const rg_REAL two_times_k_att   = 0.16;
		const rg_REAL three_times_k_att = 0.24;

		rg_Point3D centerOfAtom1 = atom1->getpAtomBall()->getCenter();
		rg_REAL    radiusOfAtom1 = atom1->getpAtomBall()->getRadius();

		rg_Point3D centerOfAtom2 = atom2->getpAtomBall()->getCenter();
		rg_REAL    radiusOfAtom2 = atom2->getpAtomBall()->getRadius();

		rg_REAL distance = 0.0;
		distance = centerOfAtom1.distance( centerOfAtom2 );
		rg_REAL sumOfRadii = radiusOfAtom1 + radiusOfAtom2;
		rg_REAL twiceSumOfRadii = 2.0*sumOfRadii;

		rg_REAL derivativeForDist = 0.0;	
		Vector3D gradientVec;

		if( rg_LE(distance, sumOfRadii) )
		{
			derivativeForDist = -k_rep / sumOfRadii;

			rg_REAL x1, y1, z1;
			centerOfAtom1.getCoordinates(x1, y1, z1);
			rg_REAL x2, y2, z2;
			centerOfAtom2.getCoordinates(x2, y2, z2);

			gradientVec.setPoint( derivativeForDist * (x1 - x2) / distance, derivativeForDist * (y1 - y2) / distance, derivativeForDist * (z1 - z2) / distance);
		}
		else if( rg_GT(distance, sumOfRadii) && rg_LT(distance, twiceSumOfRadii))
		{
			derivativeForDist = two_times_k_att * distance / (sumOfRadii * sumOfRadii) - three_times_k_att / sumOfRadii;

			rg_REAL x1, y1, z1;
			centerOfAtom1.getCoordinates(x1, y1, z1);
			rg_REAL x2, y2, z2;
			centerOfAtom2.getCoordinates(x2, y2, z2);

			gradientVec.setPoint( derivativeForDist * (x1 - x2) / distance, derivativeForDist * (y1 - y2) / distance, derivativeForDist * (z1 - z2) / distance);
		}
		else
		{
			gradientVec.setPoint(0.0, 0.0, 0.0);
		}

		return gradientVec;
	}
}



rg_BOOL     V::GeometryTier::detectClashBetweenRotamersByEvaluatingVanDerWaalsEnergy( 
                                    Rotamer*  targetRotamer,
                                    const rg_INT& targetResidueIndex, 
                                    ManagerOfRotamers_Assigned_At_Residues* assignedRotamerSetOfProtein,  
                                    Backbone* backbone,
                                    rg_dList<rg_INDEX>& clashedResidueIndeices )
{
	rg_BOOL bClashed = rg_FALSE;
	Rotamer** rotamersOrderedByChainIDSeqNum  = assignedRotamerSetOfProtein->getRotamersOrderedByChainIDSeqNum();
	rg_INT    numAssignedRotamers = assignedRotamerSetOfProtein->getNumRotamers();

	//rg_REAL energyWithBackbone = 0.0;
	//energyWithBackbone = targetRotamer->computeEnergyWithBackboneBySCWRL3( backbone );

	rg_INDEX i;
	for (i = 0;i < numAssignedRotamers;i++)
	{
		if(i != targetResidueIndex)
		{			
			rg_REAL energyWithOtherRotamer = 0.0;
			energyWithOtherRotamer = targetRotamer->computeEnergyWithOtherRotamerBySCWRL3( *rotamersOrderedByChainIDSeqNum[ i ] );
			
			rg_REAL energy = 0.0;
			energy = /*energyWithBackbone + */ energyWithOtherRotamer;

			if( rg_GT(energy, 2) )
			{
				clashedResidueIndeices.add( i );
			}
		}
	}

	if(clashedResidueIndeices.getSize() > 0)
	{
		bClashed = rg_TRUE;
	}

	return bClashed;
}



void    V::GeometryTier::getRotamerSetsForResiduesSortedByChainIDSeqNumber(
                                Rotamer* sortedAssignedRotamersOfProtein, 
                                const rg_INT& numRotamers, 
                                RotamerSetOfResidue*& sortedRotamerSetsOfResidues)
{
	sortedRotamerSetsOfResidues = rg_NULL;
	rg_INT numRotamerSets = numRotamers;
	sortedRotamerSetsOfResidues = new RotamerSetOfResidue[ numRotamerSets ];

	rg_INDEX i;
	for (i = 0;i < numRotamers;i++)
	{
		sortedRotamerSetsOfResidues[ i ].set( &sortedAssignedRotamersOfProtein[ i ] );
	}
}



rg_REAL     V::GeometryTier::computeXVolumeOfProtein(Molecule& protein)
{
	// get backbone atoms
	Atom** backboneAtoms = rg_NULL;
	rg_INT numBackboneAtoms
		= protein.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms );

	// get target residues
	Residue** sortedResidues = rg_NULL;
	rg_INT numResidues = protein.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

	rg_REAL XVolumeOfSidechainNBackbone = 0;
	XVolumeOfSidechainNBackbone = computeXVolumeBetweenSidechainsNBackbone(sortedResidues, numResidues, backboneAtoms, numBackboneAtoms);

	if(backboneAtoms != rg_NULL)
		delete[] backboneAtoms;

	rg_REAL XVolumeOfSidechainPairs = 0;
	XVolumeOfSidechainPairs = computeXVolumeBetweenSidechainPairs(sortedResidues, numResidues);

	if(sortedResidues != rg_NULL)
		delete [] sortedResidues;

	rg_REAL XVolume = 0.0;
	XVolume = XVolumeOfSidechainNBackbone + XVolumeOfSidechainPairs;

	return XVolume;
}



rg_REAL     V::GeometryTier::computeXVolumeBetweenSidechainsNBackbone(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues, 
                                    V::GeometryTier::Atom** backboneAtoms, 
                                    const rg_INT& numBackboneAtoms)
{
	rg_REAL XVolumeOfSidechainNBackbone = 0.0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_REAL currXVolume = 0.0;
		currXVolume = residues[ i ]->computeSumOfPairwiseAtomXVolumesOfSidechainWithItsBackbone();
		XVolumeOfSidechainNBackbone += currXVolume;
	}

	return XVolumeOfSidechainNBackbone;
}



rg_REAL     V::GeometryTier::computeXVolumeBetweenSidechainPairs(
                                    V::GeometryTier::Residue** residues, 
                                    const rg_INT& numResidues)
{
	rg_REAL XVolumeOfSidechainPairs = 0;

	rg_INDEX i;
	for(i = 0;i < numResidues;i++)
	{
		rg_INDEX j;
		for(j = i + 1;j < numResidues;j++)
		{
			rg_REAL currXVolume = 0;
			currXVolume = residues[ i ]->computeSumOfPairwiseAtomXVolumesOfSidechainWithOtherResidue(residues[ j ]);

			XVolumeOfSidechainPairs += currXVolume;
		}
	}

	return XVolumeOfSidechainPairs;
}



rg_BOOL     V::GeometryTier::haveIntersectionBetweenTwoResidues(
                                    V::GeometryTier::Residue* residue1, 
                                    V::GeometryTier::Residue* residue2)
{
	rg_dList<Atom*> atomsOfResidue1;	
	FunctionsForRotamerLibrary::getAtomsOfResidue(residue1, atomsOfResidue1);

	rg_dList<Atom*> atomsOfResidue2;	
	FunctionsForRotamerLibrary::getAtomsOfResidue(residue1, atomsOfResidue2);

	rg_BOOL bIntersected = rg_FALSE;

	atomsOfResidue1.reset4Loop();
	while (atomsOfResidue1.setNext4Loop())
	{
		Atom* currAtom = atomsOfResidue1.getEntity();
		bIntersected = haveIntersectionBetweenAtomWithAtomSet(currAtom, atomsOfResidue2);
		if(bIntersected)
			break;
	}
	return bIntersected;
}



//void computeSumPairwiseXVolOfBallWithBallSet(const Sphere& ball1, rg_dList<Sphere>& ballSet2, rg_REAL& xVolume)
//{
//	EnclosingSphereOfSpheres enclosingSphere2;
//	enclosingSphere2.setSpheres(ballSet2);
//	enclosingSphere2.computeEnclosingSphere();
//
//	xVolume = 0.0;	
//	Circle3D xCircle;
//	if(ball1.intersect(enclosingSphere2, xCircle) != SI_NOT_INTERSECT)
//	{
//		ballSet2.reset4Loop();
//		while (ballSet2.setNext4Loop())
//		{
//			Sphere currBall2 = ballSet2.getEntity();
//			xVolume += ball1.computeIntersectingVolume( currBall2 );
//		}		
//	}
//}



rg_REAL     V::GeometryTier::computeSumPairwiseXVolOfAtomWithAtomSet(
                                    V::GeometryTier::Atom* atom1, 
                                    rg_dList<V::GeometryTier::Atom*>& atomSet2)
{
	// test
	//EnclosingSphereOfSpheres enclosingSphere2;
	//enclosingSphere2.setSpheres(atomSet2);
	//enclosingSphere2.computeEnclosingSphere();

	Sphere atomBall1 = atom1->getAtomBall();
	rg_REAL radiusOfAtomBall1 = atomBall1.getRadius();
	rg_Point3D centerOfAtomBall1 = atomBall1.getCenter();

	rg_REAL sumXVolume = 0.0;	
	//if(atomBall1.intersect(enclosingSphere2) != SI_NOT_INTERSECT)
	{
		atomSet2.reset4Loop();
		while (atomSet2.setNext4Loop())
		{
			Atom* currAtom2 = atomSet2.getEntity();
			if(atom1->hasChemicalBondWith( currAtom2 ))
				continue;

			Sphere currAtomBall2 = currAtom2->getAtomBall();			
			rg_REAL xVolume = 0.0;	
			//xVolume = atomBall1.computeIntersectingVolume( currAtomBall2 );
			//xVolume = atomBall1.computeIntersectingVolumeA( currAtomBall2 );
			
			// test
			rg_REAL radiusOfCurrAtomBall2 = currAtomBall2.getRadius();
			rg_REAL distBTWCenters = centerOfAtomBall1.distance(currAtomBall2.getCenter());
			xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfCurrAtomBall2, distBTWCenters);

			AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
			AmberAtomTypes atomType2 = currAtom2->getpChemicalProperties()->getAtomTypeInAmber();
			sumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, distBTWCenters, xVolume );
			//sumXVolume += xVolume;



			//sumXVolume += pow(xVolume, 6);
			//sumXVolume += computeVanDerWaalsEnergyBetweenTwoAtoms(atom1, currAtom2);
		}
	}

	return sumXVolume;
}



rg_REAL     V::GeometryTier::computeSumPairwiseXVolOfAtomWithAtomSet(
                                    V::GeometryTier::Atom* atom1, 
                                    rg_dList<V::GeometryTier::Atom*>& atomSet2, 
                                    rg_INT& numXAtoms)
{	
	Sphere atomBall1 = atom1->getAtomBall();
	rg_REAL radiusOfAtomBall1 = atomBall1.getRadius();
	rg_Point3D centerOfAtomBall1 = atomBall1.getCenter();

	numXAtoms = 0;
	rg_REAL sumXVolume = 0.0;

	atomSet2.reset4Loop();
	while (atomSet2.setNext4Loop())
	{
		Atom* currAtom2 = atomSet2.getEntity();
		if(atom1->hasChemicalBondWith( currAtom2 ))
			continue;

		Sphere currAtomBall2 = currAtom2->getAtomBall();
		rg_REAL xVolume = 0.0;	
		//xVolume = atomBall1.computeIntersectingVolume( currAtomBall2 );
		//xVolume = atomBall1.computeIntersectingVolumeA( currAtomBall2 );
		rg_REAL radiusOfCurrAtomBall2 = currAtomBall2.getRadius();
		rg_REAL distBTWCenters = centerOfAtomBall1.distance(currAtomBall2.getCenter());
		xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfCurrAtomBall2, distBTWCenters);

		if(rg_POS( xVolume ))
			numXAtoms++;
		
		// test
		//rg_REAL distBTWCenters = atomBall1.getCenter().distance(currAtomBall2.getCenter());
		AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
		AmberAtomTypes atomType2 = currAtom2->getpChemicalProperties()->getAtomTypeInAmber();
		sumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, distBTWCenters, xVolume );
		//sumXVolume += xVolume;

		//rg_REAL r1 = atomBall1.getRadius();
		//rg_REAL r2 = currAtomBall2.getRadius();
		//rg_REAL sumOfRadii = r1 + r2;
		//for(rg_INDEX z = 0;z < 99;z++)
		//{
		//	rg_REAL dist = sumOfRadii * (2. - (2. * z + 1.) / 100) / 2.;
		//	rg_REAL xvol = rg_GeoFunc::computeXVolumeOfTwoSpheres(r1, r2, dist);
		//	FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, dist, xvol);
		//}
		//rg_REAL xvol = rg_GeoFunc::computeXVolumeOfTwoSpheres(r1, r2, sumOfRadii);
		//FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, sumOfRadii, xvol);
		//xvol = rg_GeoFunc::computeXVolumeOfTwoSpheres(r1, r2, sumOfRadii / 100.);
		//FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, sumOfRadii / 100., xvol);

		//xvol = rg_GeoFunc::computeXVolumeOfTwoSpheres(r1, r2, sumOfRadii + 1.);
		//FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, sumOfRadii + 1., xvol);
		//xvol = rg_GeoFunc::computeXVolumeOfTwoSpheres(r1, r2, sumOfRadii / 101.);
		//FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomType1, atomType2, sumOfRadii / 101., xvol);

		//sumXVolume += xVolume;



		//sumXVolume += pow(xVolume, 6);
		//sumXVolume += computeVanDerWaalsEnergyBetweenTwoAtoms(atom1, currAtom2);
	}

	return sumXVolume;
}



rg_REAL     V::GeometryTier::computeSumPairwiseXVolBetweenTwoAtomSets( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_INT& numXAtoms )
{
	rg_REAL sumXVolume = 0.0;
    numXAtoms = 0;

	for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
	{
		rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
		rg_REAL radiusOfAtomBall1 = atomSet1[i]->getpAtomBall()->getRadius();
        AmberAtomTypes atomTypeInAmberForAtom1  = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();
        
		for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
		{
			if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
				continue;

			rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
            rg_REAL radiusOfAtomBall2 = atomSet2[j]->getpAtomBall()->getRadius();
			AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();
            
		    rg_REAL distBTWCenters = centerOfAtom1.distance(centerOfAtom2);
            rg_REAL xVolume = 0.0;
            xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);
            		
            if(rg_POS( xVolume ))
                numXAtoms++;

			sumXVolume += xVolume;
            //sumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters, xVolume );
		}
	}

    return sumXVolume; 
}



rg_REAL     V::GeometryTier::computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_INT& numXAtoms )
{
	// March 30, 2016 by Joonghyun
    // test
    return computeNewScoreFunction( atomSet1, numOfAtomSet1, atomSet2, numOfAtomSet2, numXAtoms );

	rg_REAL sumXVolume = 0.0;
    numXAtoms = 0;

	for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
	{
		rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
		rg_REAL radiusOfAtomBall1 = atomSet1[i]->getpAtomBall()->getRadius();
        AmberAtomTypes atomTypeInAmberForAtom1  = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();
        
        rg_REAL currSumXVolume = 0.0;

		for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
		{
			if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
				continue;

			rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
            rg_REAL radiusOfAtomBall2 = atomSet2[j]->getpAtomBall()->getRadius();
			AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();
            
#ifdef BETASCP_STAND_ALONE_VERSION
			// March 30, 2016 by Joonghyun
			radiusOfAtomBall1 = radiusOfAtomBall1 * (1.0 + RADIUS_INFLATION_RATE);
            radiusOfAtomBall2 = radiusOfAtomBall2 * (1.0 + RADIUS_INFLATION_RATE);
#endif

		    rg_REAL distBTWCenters = centerOfAtom1.distance(centerOfAtom2);            
            rg_REAL xVolume = 0.0;
            xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);

            if(rg_POS( xVolume ))
                numXAtoms++;

			// March 30, 2016 by Joonghyun
            //sumXVolume += xVolume;
            sumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters, xVolume );
		}
	}

    return sumXVolume; 
}



rg_REAL     V::GeometryTier::computeNewScoreFunction_old( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_INT& numXAtoms )
{
    rg_REAL sumXVolume = 0.0;
    numXAtoms = 0;
    // test
    //sumXVolume = computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3(atomSet1, numOfAtomSet1, atomSet2, numOfAtomSet2, numXAtoms);
    //sumXVolume = computeVanDerWaalsEnergyBetweenTwoAtomSets(atomSet1, numOfAtomSet1, atomSet2, numOfAtomSet2);
    //return sumXVolume;

    for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
    {
        rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
        rg_REAL radiusOfAtomBall1 = atomSet1[i]->getpAtomBall()->getRadius();
        AmberAtomTypes atomTypeInAmberForAtom1  = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();

        // test
        //rg_REAL minXVolume = rg_REAL_INFINITY;
        rg_REAL currSumXVolume = 0.0;

        for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
        {
            if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
                continue;

            rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
            rg_REAL radiusOfAtomBall2 = atomSet2[j]->getpAtomBall()->getRadius();
            AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();

            rg_REAL distBTWCenters = centerOfAtom1.distance(centerOfAtom2);

            // test
            //radiusOfAtomBall1 *= 1.1;
            //radiusOfAtomBall2 *= 1.1;
            //radiusOfAtomBall1 *= 0.9;
            //radiusOfAtomBall2 *= 0.9;

            rg_REAL dist = distBTWCenters - (radiusOfAtomBall1 + radiusOfAtomBall2);

            //xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);
            //sumXVolume += fabs(xVolume);
            //sumXVolume += xVolume;
            //sumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters, xVolume );
            //sumXVolume +=scale( xVolume );

            //if(rg_LE(xVolume, 5.0))
            //{
            //    currSumXVolume += fabs(xVolume);
            //    if(rg_NEG( xVolume ))
            //    {
            //        numXAtoms++;
            //    }
            //}

            rg_REAL xVolume = 0.0;
            if(rg_LT(dist, 0.0))
            {                
                xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);
                //xVolume = xVolume * xVolume * xVolume * xVolume;
                //xVolume = xVolume + fabs(dist);
                //xVolume = xVolume * fabs(dist);
                //if(rg_GT(xVolume, 0.0))
                numXAtoms++;
            }
            else
            {
                //if(rg_LE(xVolume, 4.0)) 이 조건 없는 경우, 4.0 or 5.0 모두 동일한 결과(HiQ45에 대해서)
                xVolume = dist;
            }
            currSumXVolume += xVolume;
        }

        //sumXVolume += currSumXVolume / numOfAtomSet2;
        //if(numOfAtomSet2 > 0)
        //    sumXVolume += minXVolume / numOfAtomSet2;
        //if(numOfAtomSet2 > 0)
        //    currSumXVolume = currSumXVolume / numOfAtomSet2;
        sumXVolume += currSumXVolume;
    }

    // test
    //if(numXAtoms > 1)
    //    sumXVolume = sumXVolume / numXAtoms;

    //if(rg_GT(sumXVolume, 1.0))
    //    sumXVolume = sqrt(sumXVolume);

    if(numOfAtomSet1 > 1)
        sumXVolume = sumXVolume / numOfAtomSet1;

    return sumXVolume; 

    // minXVolume
    // rg_LE(xVolume, 5.0), currSumXVolume += fabs(xVolume), numOfAtomSet1로만 나누어 준 경우 최적임.
}



rg_REAL     V::GeometryTier::computeNewScoreFunction( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_INT& numXAtoms )
{
    rg_REAL sumXVolume = 0.0;
    numXAtoms = 0;

    for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
    {
        rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
        rg_REAL radiusOfAtomBall1 = atomSet1[i]->getpAtomBall()->getRadius();
        AmberAtomTypes atomTypeInAmberForAtom1  = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();

        rg_REAL currSumXVolume = 0.0;

        for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
        {
            if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
                continue;

            rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
            rg_REAL radiusOfAtomBall2 = atomSet2[j]->getpAtomBall()->getRadius();
            AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();

            rg_REAL distBTWCenters = centerOfAtom1.distance(centerOfAtom2);
            
#ifdef BETASCP_STAND_ALONE_VERSION
			// March 28, 2016 Joonghyun
            // inflate radii 10%
            radiusOfAtomBall1 = radiusOfAtomBall1 * (1.0 + RADIUS_INFLATION_RATE);
            radiusOfAtomBall2 = radiusOfAtomBall2 * (1.0 + RADIUS_INFLATION_RATE);
#endif
            rg_REAL dist = distBTWCenters - (radiusOfAtomBall1 + radiusOfAtomBall2);

            rg_REAL xVolume = 0.0;
            if(rg_LT(dist, 0.0))
            {                
                xVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);
				currSumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters, xVolume );
                numXAtoms++;
            }
            else
            {
            	// March 28, 2016 Joonghyun
                //xVolume = dist;
                //currSumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters, xVolume );
                
				// March 25, 2016 Joonghyun
				rg_REAL squaredDist = distBTWCenters * distBTWCenters;
				rg_REAL sixthMultipleOfDist = squaredDist * squaredDist * squaredDist;	
				xVolume = -L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForAtom1][atomTypeInAmberForAtom2][1] / sixthMultipleOfDist;
#ifdef BETASCP_STAND_ALONE_VERSION				
				// May 11, 2016 by Joonghyun
				currSumXVolume += (LAMBDA_COEFF_COMBINATORIAL_SEARCH * xVolume);
				//currSumXVolume += (ALPHA_COEFF_COMBINATORIAL_SEARCH * xVolume);
#endif

#ifndef BETASCP_STAND_ALONE_VERSION
				currSumXVolume += (xVolume);
#endif
            }
            //currSumXVolume += xVolume;
			// March 25, 2016 Joonghyun
			// test
			//currSumXVolume += FunctionsForMappingFromXVolToLJPotEnergy::getLJPotEnergyValueCorrToXVol(atomTypeInAmberForAtom1, atomTypeInAmberForAtom2, distBTWCenters, xVolume );
        }

        sumXVolume += currSumXVolume;
    }
	// March 25, 2016 Joonghyun
    //if(numOfAtomSet1 > 1)
    //    sumXVolume = sumXVolume / numOfAtomSet1;

    return sumXVolume; 
}



rg_REAL     V::GeometryTier::scale(const rg_REAL& valueBeforeScaling)
{
    // logistics CDF
    rg_REAL m = 0.0;
    rg_REAL s = 0.35;
    rg_REAL mult = 10.0;

    rg_REAL scaledValue = mult / ( 1.0 + exp( ( valueBeforeScaling + m ) / s ) );

    return scaledValue;    
}



rg_REAL     V::GeometryTier::computeSumPairwiseXAreaBetweenTwoAtomSets( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_INT& numXAtoms )
{
    rg_REAL sumXArea = 0.0;
    numXAtoms = 0;

    for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
    {
        rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
        rg_REAL radiusOfAtomBall1 = atomSet1[i]->getpAtomBall()->getRadius();
        AmberAtomTypes atomTypeInAmberForAtom1  = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();

        for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
        {
            if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
                continue;

            rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
            rg_REAL radiusOfAtomBall2 = atomSet2[j]->getpAtomBall()->getRadius();
            AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();

            rg_REAL distBTWCenters = centerOfAtom1.distance(centerOfAtom2);
            rg_REAL xArea = 0.0;
            xArea = rg_GeoFunc::computeXAreaOfTwoSphere(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);

            if(rg_POS( xArea ))
                numXAtoms++;

            sumXArea += xArea;
        }
    }

    return sumXArea; 
}



rg_REAL     V::GeometryTier::computeSumPairwiseXCircleAreaBetweenTwoAtomSets( 
                                    V::GeometryTier::Atom** atomSet1, 
                                    const rg_INT& numOfAtomSet1, 
                                    V::GeometryTier::Atom** atomSet2, 
                                    const rg_INT& numOfAtomSet2, 
                                    rg_REAL& sumXCircleRadii, 
                                    rg_INT& numXAtoms )
{
    rg_REAL sumXCircleAreas = 0.0;
    sumXCircleRadii = 0.0;
    numXAtoms = 0;

    for( rg_INT i=0; i<numOfAtomSet1; i++ ) 
    {
        rg_Point3D centerOfAtom1 = atomSet1[i]->getpAtomBall()->getCenter();
        rg_REAL radiusOfAtomBall1 = atomSet1[i]->getpAtomBall()->getRadius();
        AmberAtomTypes atomTypeInAmberForAtom1  = atomSet1[i]->getpChemicalProperties()->getAtomTypeInAmber();

        for( rg_INT j=0; j<numOfAtomSet2; j++ ) 
        {
            if(atomSet1[i]->hasChemicalBondWith( atomSet2[j] )) // inline function
                continue;

            rg_Point3D centerOfAtom2 = atomSet2[j]->getpAtomBall()->getCenter();
            rg_REAL radiusOfAtomBall2 = atomSet2[j]->getpAtomBall()->getRadius();
            AmberAtomTypes atomTypeInAmberForAtom2 = atomSet2[j]->getpChemicalProperties()->getAtomTypeInAmber();

            rg_REAL distBTWCenters = centerOfAtom1.distance(centerOfAtom2);

            rg_REAL xCircleRadius = 0.0;
            xCircleRadius = rg_GeoFunc::computeXCircleRadiusOfTwoSphere(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);

            if(rg_POS( xCircleRadius ))
            {
                rg_REAL xCircleArea = 0.0;
                xCircleArea = rg_GeoFunc::computeXCircleAreaOfTwoSphere(radiusOfAtomBall1, radiusOfAtomBall2, distBTWCenters);

                sumXCircleRadii += xCircleRadius;
                sumXCircleAreas += xCircleArea;
                numXAtoms++;
            }
        }
    }

    return sumXCircleAreas; 
}



rg_BOOL     V::GeometryTier::haveIntersectionBetweenAtomWithAtomSet(
                                    V::GeometryTier::Atom* atom1, 
                                    rg_dList<V::GeometryTier::Atom*>& atomSet2)
{
	EnclosingSphereOfSpheres enclosingSphere2;
	enclosingSphere2.setSpheres(atomSet2);
	enclosingSphere2.computeEnclosingSphere();

	Sphere atomBall1 = atom1->getAtomBall();

	rg_BOOL bIntersected = rg_FALSE;
	if(atomBall1.intersect(enclosingSphere2) != SI_NOT_INTERSECT)
	{
		atomSet2.reset4Loop();
		while (atomSet2.setNext4Loop())
		{
			Atom* currAtom2 = atomSet2.getEntity();
			Sphere currAtomBall2 = currAtom2->getAtomBall();
			if(atomBall1.intersect( currAtomBall2 ) != SI_NOT_INTERSECT)
			{
				bIntersected = rg_TRUE;
				break;
			}			
		}
	}

	return bIntersected;
}



void    V::GeometryTier::computeSumPairwiseXVolOfAtomWithAtomSetWithoutChemBond(
                                V::GeometryTier::Atom* atom1, 
                                rg_dList<V::GeometryTier::Atom*>& atomSet2, 
                                rg_REAL& xVolume)
{
	EnclosingSphereOfSpheres enclosingSphere2;
	enclosingSphere2.setSpheres(atomSet2);
	enclosingSphere2.computeEnclosingSphere();

	Sphere atomBall1 = atom1->getAtomBall();

	xVolume = 0.0;	
	Circle3D xCircle;
	if(atomBall1.intersect(enclosingSphere2, xCircle) != SI_NOT_INTERSECT)
	{
		atomSet2.reset4Loop();
		while (atomSet2.setNext4Loop())
		{
			Atom* currAtom2 = atomSet2.getEntity();

			// check whether two atoms are connected by chemical bond
			if(!isThereBondBetween(*atom1, *currAtom2))
			{
				Sphere currAtomBall2 = currAtom2->getAtomBall();
				//xVolume += atomBall1.computeIntersectingVolume( currAtomBall2 );
                xVolume += rg_GeoFunc::computeXVolumeOfTwoSpheres( atomBall1, currAtomBall2 );
			}
		}
	}
}



rg_REAL     V::GeometryTier::computeRMSDBetweenTwoProteinsSharingCommonBackboneNHavingDifferentSidechains(
                                    Molecule& protein1, 
                                    Molecule& protein2)
{
	protein1.moveOxygenInCarboxylGroupToBackbone();
	protein2.moveOxygenInCarboxylGroupToBackbone();

	Residue** residuesOrderedByChainIDSeqNum1 = rg_NULL;
	rg_INT numResidues1 = 0;
	numResidues1 = protein1.getAminoResiduesOrderedByChainIDSequenceNumber(residuesOrderedByChainIDSeqNum1);

	rg_INT numResidues2 = 0;
	Residue** residuesOrderedByChainIDSeqNum2 = rg_NULL;
	numResidues2 = protein2.getAminoResiduesOrderedByChainIDSequenceNumber(residuesOrderedByChainIDSeqNum2);

	rg_REAL rmsd = 0.0;
	rg_INT  numAtoms = 0;
	rg_INDEX i;
	for (i = 0;i < numResidues1;i++)
	{
		rg_dList<Atom*> atomsOnSidechain1;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(residuesOrderedByChainIDSeqNum1[ i ], atomsOnSidechain1 );

		rg_dList<Atom*> atomsOnSidechain2;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(residuesOrderedByChainIDSeqNum2[ i ], atomsOnSidechain2 );

		rg_REAL squaredDiffOfResiduePair = 0.0;
		atomsOnSidechain1.reset4Loop();
		atomsOnSidechain2.reset4Loop();
		while (atomsOnSidechain1.setNext4Loop() && atomsOnSidechain2.setNext4Loop())
		{
			Atom* currAtom1 = rg_NULL;
			currAtom1 = atomsOnSidechain1.getEntity();
			rg_Point3D atomCenter1 = currAtom1->getAtomBall().getCenter();

			Atom* currAtom2 = rg_NULL;
			currAtom2 = atomsOnSidechain2.getEntity();
			rg_Point3D atomCenter2 = currAtom2->getAtomBall().getCenter();

			rg_REAL currSD = 0.0;
			currSD = atomCenter1.squaredDistance( atomCenter2 );
			squaredDiffOfResiduePair += currSD;
		}

		rmsd += squaredDiffOfResiduePair;
		numAtoms += atomsOnSidechain1.getSize();
	}

	rmsd = sqrt(rmsd / numAtoms);

	if(residuesOrderedByChainIDSeqNum1 != rg_NULL)
		delete [] residuesOrderedByChainIDSeqNum1;

	if(residuesOrderedByChainIDSeqNum2 != rg_NULL)
		delete [] residuesOrderedByChainIDSeqNum2;

	return rmsd;
}



rg_BOOL     V::GeometryTier::computeRMSDsBetweenAminoResiduePairsOfTwoProteinsSharingFixedCommonBackbone(
                                    Molecule& protein1,
                                    Molecule& protein2,
                                    rg_REAL& rmsdOfProtein,
                                    vector<rg_REAL>& RMSDsOrderedByChainIDSeqNumbers)
{
	protein1.moveOxygenInCarboxylGroupToBackbone();
	protein2.moveOxygenInCarboxylGroupToBackbone();

	rg_dList<Residue*> residuesOrderedByChainIDSeqNum1;
	protein1.getAminoResiduesOrderedByChainIDSequenceNumber(residuesOrderedByChainIDSeqNum1);

	rg_dList<Residue*> residuesOrderedByChainIDSeqNum2;
	protein2.getAminoResiduesOrderedByChainIDSequenceNumber(residuesOrderedByChainIDSeqNum2);

    rg_REAL sumOfSquaredDiff = 0.0;
    rg_INT  sumOfNumAtoms = 0;

	residuesOrderedByChainIDSeqNum1.reset4Loop();
	residuesOrderedByChainIDSeqNum2.reset4Loop();
	while(residuesOrderedByChainIDSeqNum1.setNext4Loop() &&
		  residuesOrderedByChainIDSeqNum2.setNext4Loop()   )
	{
		Residue* currResidue1 = residuesOrderedByChainIDSeqNum1.getEntity();
		Residue* currResidue2 = residuesOrderedByChainIDSeqNum2.getEntity();

        ResidueCode code1 = currResidue1->getResidueCode();
        ResidueCode code2 = currResidue2->getResidueCode();
        if(code1 != code2)
            return rg_FALSE;
        if( code1 == GLY_AMINO_RESIDUE || code1 == ALA_AMINO_RESIDUE )
        {
            RMSDsOrderedByChainIDSeqNumbers.push_back( 0.0 );
            continue;
        }            

		rg_dList<Atom*> atomsOnSidechain1;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(currResidue1, atomsOnSidechain1);

		rg_dList<Atom*> atomsOnSidechain2;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(currResidue2, atomsOnSidechain2);

        rg_INT numAtoms1 = atomsOnSidechain1.getSize();
        rg_INT numAtoms2 = atomsOnSidechain2.getSize();

		// April 2, 2016 by Joonghyun
		//if(numAtoms1 != numAtoms2)
		//	return rg_FALSE;
		
		rg_REAL squaredDiffOfResiduePair = 0.0;        

		atomsOnSidechain1.reset4Loop();
		atomsOnSidechain2.reset4Loop();
		while (atomsOnSidechain1.setNext4Loop() && atomsOnSidechain2.setNext4Loop())
		{
			Atom* currAtom1 = rg_NULL;
			currAtom1 = atomsOnSidechain1.getEntity();
			rg_Point3D atomCenter1 = currAtom1->getpAtomBall()->getCenter();

			Atom* currAtom2 = rg_NULL;
			currAtom2 = atomsOnSidechain2.getEntity();
			rg_Point3D atomCenter2 = currAtom2->getpAtomBall()->getCenter();

			rg_REAL currSD = 0.0;
			currSD = atomCenter1.squaredDistance( atomCenter2 );
			squaredDiffOfResiduePair += currSD;
		}
                
		rg_REAL rmsdOfResiduePair = sqrt(squaredDiffOfResiduePair / numAtoms1);
		RMSDsOrderedByChainIDSeqNumbers.push_back(rmsdOfResiduePair);

        sumOfSquaredDiff += squaredDiffOfResiduePair;
        sumOfNumAtoms    += numAtoms1;
	}

    rmsdOfProtein = sqrt(sumOfSquaredDiff / sumOfNumAtoms);

	return rg_TRUE;
}



void    V::GeometryTier::writeDihedralAnglesOfAminoResidueSidechains(
                                Molecule& molecule,     
                                const string& sidechainDihedralAngleFIle, 
                                const rg_BOOL& bAppend /* = rg_FALSE */)
{
    ofstream fout;
    if(bAppend)
    {
        fout.open(sidechainDihedralAngleFIle.c_str(), ios::app);        
    }
    else
    {
        fout.open(sidechainDihedralAngleFIle.c_str());
        writeHeaderForDihedralAnglesOfAminoResidueSidechains(fout);
    }

    writeBodyForDihedralAnglesOfAminoResidueSidechains(molecule, fout);
}



void    V::GeometryTier::writeHeaderForDihedralAnglesOfAminoResidueSidechains(ofstream& fout)
{
    fout << "  CODE  CHAIN    SEQ   TYPE       phi       psi      chi1      chi2      chi3    chi4 " << endl;
    fout << "--------------------------------------------------------------------------------------" << endl;
}



void    V::GeometryTier::writeBodyForDihedralAnglesOfAminoResidueSidechains(Molecule& molecule, ofstream& fout)
{
    // get PDB ID
    string fileNameWithoutPathAndExt = StringFunctions::getFileNameWithoutPathAndExtension( molecule.getMoleculeFileName() );
    rg_INT strLength = fileNameWithoutPathAndExt.size();
    rg_INDEX stPosOfPDBCode = strLength - 4;
    string PDBCODE = fileNameWithoutPathAndExt.substr(stPosOfPDBCode, 4);

    rg_dList<Residue*> aminoResiduesOrderedByChainIDSeqNum;
    molecule.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResiduesOrderedByChainIDSeqNum);

    aminoResiduesOrderedByChainIDSeqNum.reset4Loop();
    while (aminoResiduesOrderedByChainIDSeqNum.setNext4Loop())
    {
        Residue* currResidue = aminoResiduesOrderedByChainIDSeqNum.getEntity();

        string chainID = currResidue->getChain()->getChainIDFromInputFileInString();
        rg_INT seqNum  = currResidue->getSequenceNumber();
        ResidueCode code = currResidue->getResidueCode();
        rg_REAL phi, psi;
        currResidue->getBackBoneDihedralAngles(phi, psi);
        
        fout.setf(ios::fixed, ios::floatfield);
        fout.precision(3);

        fout.width( 6 );
        fout << PDBCODE;
        fout.width( 7 );
        fout << chainID;
        fout.width( 7 );
        fout << seqNum;
        fout.width( 7 );
        fout << RESIDUE_FEATURES[code].threeCodeName;
        fout.width( 10 );
        fout << phi;
        fout.width( 10 );
        fout << psi;
        
        rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(code);
        rg_REAL* dihedralAngles = new rg_REAL[numDihedralAngles];
        currResidue->getSideChainDihedralAngles(dihedralAngles);

        rg_INDEX i;
        for (i = 0;i < numDihedralAngles;i++)
        {
            fout.width( 10 );
            fout << dihedralAngles[ i ];
        }

        fout << endl;

        if(dihedralAngles != rg_NULL)
            delete [] dihedralAngles;
    }
}



void    V::GeometryTier::writeDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(
                                Molecule& referenceStructure,
                                Molecule& targetStructure, 
                                const string& sidechainDihedralAngleFIle, 
                                const rg_BOOL& bAppend /* = rg_FALSE */)
{
    ofstream fout;
    if(bAppend)
    {
        fout.open(sidechainDihedralAngleFIle.c_str(), ios::app);        
    }
    else
    {
        fout.open(sidechainDihedralAngleFIle.c_str());
        writeHeaderForDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(fout);
    }

    writeBodyForDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(referenceStructure, targetStructure, fout);
}



void    V::GeometryTier::writeHeaderForDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(ofstream& fout)
{
    fout << "  CODE  CHAIN    SEQ   TYPE       phi       psi      chi1      chi2      chi3      chi4     dchi1     dchi2     dchi3     dchi4" << endl;
    fout << "-------------------------------------------------------------------------------------------------------------------------------" << endl;
}



void    V::GeometryTier::writeBodyForDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(
                                Molecule& referenceStructure,
                                Molecule& targetStructure,
                                ofstream& fout)
{
    // get PDB ID
    string fileNameWithoutPathAndExt = StringFunctions::getFileNameWithoutPathAndExtension( referenceStructure.getMoleculeFileName() );
    rg_INT strLength = fileNameWithoutPathAndExt.size();
    rg_INDEX stPosOfPDBCode = strLength - 4;
    string PDBCODE = fileNameWithoutPathAndExt.substr(stPosOfPDBCode, 4);

    rg_dList<Residue*> aminoResiduesOfTargetMoleculeOrderedByChainIDSeqNum;
    targetStructure.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResiduesOfTargetMoleculeOrderedByChainIDSeqNum);

    rg_dList<Residue*> aminoResiduesOfReferenceMoleculeOrderedByChainIDSeqNum;
    referenceStructure.getAminoResiduesOrderedByChainIDSequenceNumber(aminoResiduesOfReferenceMoleculeOrderedByChainIDSeqNum);

    if(aminoResiduesOfTargetMoleculeOrderedByChainIDSeqNum.getSize() != aminoResiduesOfReferenceMoleculeOrderedByChainIDSeqNum.getSize())
        return;

    aminoResiduesOfTargetMoleculeOrderedByChainIDSeqNum.reset4Loop();
    aminoResiduesOfReferenceMoleculeOrderedByChainIDSeqNum.reset4Loop();

    while (aminoResiduesOfTargetMoleculeOrderedByChainIDSeqNum.setNext4Loop() && aminoResiduesOfReferenceMoleculeOrderedByChainIDSeqNum.setNext4Loop())
    {
        Residue* currTargetResidue = aminoResiduesOfTargetMoleculeOrderedByChainIDSeqNum.getEntity();
        Residue* currReferenceResidue = aminoResiduesOfReferenceMoleculeOrderedByChainIDSeqNum.getEntity();

        string chainID = currReferenceResidue->getChain()->getChainIDFromInputFileInString();
        rg_INT seqNum  = currReferenceResidue->getSequenceNumber();
        ResidueCode code = currReferenceResidue->getResidueCode();
        rg_REAL phi, psi;
        currReferenceResidue->getBackBoneDihedralAngles(phi, psi);

        fout.setf(ios::fixed, ios::floatfield);
        fout.precision(3);

        fout.width( 6 );
        fout << PDBCODE;
        fout.width( 7 );
        fout << chainID;
        fout.width( 7 );
        fout << seqNum;
        fout.width( 7 );
        fout << RESIDUE_FEATURES[code].threeCodeName;
        fout.width( 10 );
        fout << phi;
        fout.width( 10 );
        fout << psi;

        rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(code);
        rg_REAL* dihedralAnglesOfTargetResidue = new rg_REAL[numDihedralAngles];
        rg_REAL* dihedralAnglesOfReferenceResidue = new rg_REAL[numDihedralAngles];

        currTargetResidue->getSideChainDihedralAngles(dihedralAnglesOfTargetResidue);
        currReferenceResidue->getSideChainDihedralAngles(dihedralAnglesOfReferenceResidue);

        rg_INDEX i;
        for (i = 0;i < numDihedralAngles;i++)
        {
            fout.width( 10 );
            fout << dihedralAnglesOfTargetResidue[ i ];
        }

        for (i = 0;i < numDihedralAngles;i++)
        {
            fout.width( 10 );
            fout << (dihedralAnglesOfReferenceResidue[ i ] - dihedralAnglesOfTargetResidue[ i ]);
        }

        fout << endl;

        if(dihedralAnglesOfTargetResidue != rg_NULL)
            delete [] dihedralAnglesOfTargetResidue;

        if(dihedralAnglesOfReferenceResidue != rg_NULL)
            delete [] dihedralAnglesOfReferenceResidue;
    }
}



// CMKIM added on April 7th, 2011 for JHRYU
void   V::GeometryTier::getChemicalBondsInResidueStartFromBetaCarbon( 
    V::GeometryTier::Residue* aResidue, 
    rg_dList<V::GeometryTier::ChemicalBond*>& targetListOfChemicalBonds )
{
    rg_dList<Atom*> atomsOnSideChain;

    aResidue->getAtomsOnSideChain( &atomsOnSideChain );

    if ( atomsOnSideChain.getSize() == 0 )
        return;

    atomsOnSideChain.reset4Loop();
    while ( atomsOnSideChain.setNext4Loop() ) {
        Atom* currAtom = atomsOnSideChain.getEntity();
        
//         // If currAtom is a Beta carbon then continue.
//         // A chemical bond of Beta carbon will be considered with connected Gamma atom.
//         if ( currAtom->getAtomCode() == C_ATOM &&
//              currAtom->getpChemicalProperties()->getRemoteIndicator()  == BETA_REMOTE &&
//              currAtom->getpChemicalProperties()->getBrangeDesignator() == UNK_BRANCH   ) {
//             
//             continue;
//         }
        
        rg_dList<ChemicalBond*>* chemicalBonds = currAtom->getListChemicalBond();
        chemicalBonds->reset4Loop();
        while ( chemicalBonds->setNext4Loop() ) {
            ChemicalBond* currChemicalBond = chemicalBonds->getEntity();

            targetListOfChemicalBonds.addWithoutSame( currChemicalBond );
        }
    }
}



void    V::GeometryTier::getRotationalBonds( 
                                Molecule& aMolecule, 
                                rg_dList<ChemicalBond*>& rotationalBonds )
{
    rotationalBonds.removeAll();

    rg_dList<ChemicalBond>* allChemicalBonds = aMolecule.getChemicalBonds();
    
    allChemicalBonds->reset4Loop();    
    while(allChemicalBonds->setNext4Loop()) {
        ChemicalBond* currChemicalBond = allChemicalBonds->getpEntity();
        
        rotationalBonds.add(currChemicalBond);
    }

    removeBondOfHydrogenAtom(rotationalBonds);
    removeDoubleTripleAromaticAndAmideBond(rotationalBonds);
    removeThreeFoldSymmetryBond(rotationalBonds);
    removeBondOfAminoCarboxylAndHydroxylGroup(rotationalBonds);
    removeBondOfRing(rotationalBonds, allChemicalBonds);
}



void    V::GeometryTier::removeBondOfHydrogenAtom( rg_dList<V::GeometryTier::ChemicalBond*>& rotationalBonds )
{
    rotationalBonds.reset4Loop();
    while(rotationalBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = rotationalBonds.getEntity();

        Atom* firstAtom  = currChemicalBond->getFirstAtom();
        Atom* secondAtom = currChemicalBond->getSecondAtom();

        if(firstAtom->getAtomCode() == H_ATOM || secondAtom->getAtomCode() == H_ATOM) {
            rotationalBonds.killCurrent();
        }
    }
}



void    V::GeometryTier::removeDoubleTripleAromaticAndAmideBond( rg_dList<V::GeometryTier::ChemicalBond*>& rotationalBonds )
{
    rotationalBonds.reset4Loop();
    while(rotationalBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = rotationalBonds.getEntity();
        
        if( currChemicalBond->getTypeOfBond() == DOUBLE_BOND   ||
            currChemicalBond->getTypeOfBond() == TRIPLE_BOND   ||
            currChemicalBond->getTypeOfBond() == AROMATIC_BOND ||
            currChemicalBond->getTypeOfBond() == AMIDE_BOND       ) {
            
            rotationalBonds.killCurrent();
        }
    }
}



void    V::GeometryTier::removeThreeFoldSymmetryBond( rg_dList<V::GeometryTier::ChemicalBond*>& rotationalBonds )
{
    rotationalBonds.reset4Loop();
    while(rotationalBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = rotationalBonds.getEntity();

                        
        Atom* firstAtom  = currChemicalBond->getFirstAtom();
        Atom* secondAtom = currChemicalBond->getSecondAtom();


        rg_INT numOAtom  = 0;
        rg_INT numOHAtom = 0;        
        rg_INT numOSAtom = 0;
        rg_INT numHAtom  = 0;
        rg_INT numHCAtom = 0;

        
        if(firstAtom->getChemicalProperties().getAtomTypeInAmber() == P_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();

            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }

            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }

                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == O_ATM) {
                    numOAtom++;
                }
                else if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == OH_ATM) {
                    numOHAtom++;
                }
            }

            if(numOAtom > 0 && numOHAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }
        else if(firstAtom->getChemicalProperties().getAtomTypeInAmber() == S_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }

            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == OS_ATM) {
                    numOSAtom++;
                }                
            }
            
            if(numOSAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }        
        else if( firstAtom->getChemicalProperties().getAtomTypeInAmber() == CT_ATM ) {

            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == HC_ATM) {
                    numHCAtom++;
                }                
            }
            
            if(numHCAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }       
        else if( firstAtom->getChemicalProperties().getAtomTypeInAmber() == N3_ATM   ) {
            
            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == H_ATM) {
                    numHAtom++;
                }                
            }
            
            if(numHAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }   


        numOAtom  = 0;
        numOHAtom = 0;        
        numOSAtom = 0;
        numHAtom  = 0;
        numHCAtom = 0;
        
        
        if(secondAtom->getChemicalProperties().getAtomTypeInAmber() == P_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }

            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == O_ATM) {
                    numOAtom++;
                }
                else if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == OH_ATM) {
                    numOHAtom++;
                }
            }
            
            if(numOAtom > 0 && numOHAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }
        else if(secondAtom->getChemicalProperties().getAtomTypeInAmber() == S_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }

            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == OS_ATM) {
                    numOSAtom++;
                }                
            }
            
            if(numOSAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }    
        else if( secondAtom->getChemicalProperties().getAtomTypeInAmber() == CT_ATM ) {
            
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == HC_ATM) {
                    numHCAtom++;
                }                
            }
            
            if(numHCAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }       
        else if( secondAtom->getChemicalProperties().getAtomTypeInAmber() == N3_ATM   ) {
            
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 4) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == H_ATM) {
                    numHAtom++;
                }                
            }
            
            if(numHAtom > 2) {
                rotationalBonds.killCurrent();
            }
        }   
    }    
}



void    V::GeometryTier::removeBondOfAminoCarboxylAndHydroxylGroup( rg_dList<V::GeometryTier::ChemicalBond*>& rotationalBonds )
{
    rotationalBonds.reset4Loop();
    while(rotationalBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = rotationalBonds.getEntity();
                        
        Atom* firstAtom  = currChemicalBond->getFirstAtom();
        Atom* secondAtom = currChemicalBond->getSecondAtom();
        
        rg_INT numOAtom  = 0;
        rg_INT numOHAtom = 0;        
        rg_INT numO2Atom = 0;
        rg_INT numHAtom  = 0;
        rg_INT numHOAtom = 0;

        
        if(firstAtom->getChemicalProperties().getAtomTypeInAmber() == C_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();

            if(incidentChemicalBond->getSize() != 3) {
                continue;
            }

            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }

                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == O_ATM) {
                    numOAtom++;
                }
                else if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == OH_ATM) {
                    numOHAtom++;
                }
                else if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == O2_ATM) {
                    numO2Atom++;
                }
            }

            if( (numOAtom > 0 && numOHAtom > 0) ||
                (numOAtom > 0 && numO2Atom > 0)   ) {
               
                rotationalBonds.killCurrent();
            }
        }
        else if(firstAtom->getChemicalProperties().getAtomTypeInAmber() == N_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 3) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == H_ATM) {
                    numHAtom++;
                }
                
            }
            
            if( numHAtom > 1) {                
                rotationalBonds.killCurrent();
            }
        }
        else if(firstAtom->getChemicalProperties().getAtomTypeInAmber() == OH_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = firstAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 2) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(firstAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(firstAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == HO_ATM) {
                    numHOAtom++;
                }
                
            }
            
            if( numHOAtom > 0) {                
                rotationalBonds.killCurrent();
            }
        }
        

        numOAtom  = 0;
        numOHAtom = 0;        
        numO2Atom = 0;
        numHAtom  = 0;
        numHOAtom = 0;
        
        
        if(secondAtom->getChemicalProperties().getAtomTypeInAmber() == C_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();

            if(incidentChemicalBond->getSize() != 3) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == O_ATM) {
                    numOAtom++;
                }
                else if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == OH_ATM) {
                    numOHAtom++;
                }
                else if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == O2_ATM) {
                    numO2Atom++;
                }
            }
            
            if( (numOAtom > 0 && numOHAtom > 0) ||
                (numOAtom > 0 && numO2Atom > 0)    ) {
                
                rotationalBonds.killCurrent();
            }  
        }    
        else if(secondAtom->getChemicalProperties().getAtomTypeInAmber() == N_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 3) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == H_ATM) {
                    numHAtom++;
                }
                
            }
            
            if( numHAtom > 1) {                
                rotationalBonds.killCurrent();
            }
        }
        else if(secondAtom->getChemicalProperties().getAtomTypeInAmber() == OH_ATM) {
            rg_dList<ChemicalBond*>* incidentChemicalBond = secondAtom->getListChemicalBond();
            
            if(incidentChemicalBond->getSize() != 2) {
                continue;
            }
            
            incidentChemicalBond->reset4Loop();
            while(incidentChemicalBond->setNext4Loop()) {
                ChemicalBond* currIncidentChemicalBond = incidentChemicalBond->getEntity();
                
                Atom* oppositeAtom = rg_NULL;
                if(secondAtom == currIncidentChemicalBond->getFirstAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getSecondAtom();
                }
                else if(secondAtom == currIncidentChemicalBond->getSecondAtom()) {
                    oppositeAtom = currIncidentChemicalBond->getFirstAtom();
                }
                
                if(oppositeAtom->getChemicalProperties().getAtomTypeInAmber() == HO_ATM) {
                    numHOAtom++;
                }
                
            }
            
            if( numHOAtom > 0) {                
                rotationalBonds.killCurrent();
            }
        }
    }  
}



void    V::GeometryTier::removeBondOfRing( 
                                rg_dList<V::GeometryTier::ChemicalBond*>& rotationalBonds, 
                                rg_dList<V::GeometryTier::ChemicalBond>* allChemicalBonds)
{
    rotationalBonds.reset4Loop();
    while(rotationalBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = rotationalBonds.getEntity();

        rg_BOOL isOnRing = isChemicalBondOnRing(currChemicalBond, allChemicalBonds);

        if(isOnRing == rg_TRUE) {
            rotationalBonds.killCurrent();
        }
    }
}



rg_BOOL     V::GeometryTier::isChemicalBondOnRing( 
                                    V::GeometryTier::ChemicalBond* aChemicalBond, 
                                    rg_dList<V::GeometryTier::ChemicalBond>* allChemicalBonds )
{
    rg_BOOL isOnRing = rg_FALSE;
    rg_dList<ChemicalBond*> tempBonds;

    allChemicalBonds->reset4Loop();
    while(allChemicalBonds->setNext4Loop()) {
        ChemicalBond* currChemicalBond = allChemicalBonds->getpEntity();

        tempBonds.add(currChemicalBond);
    }

    tempBonds.reset4Loop();
    while(tempBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = tempBonds.getEntity();
        
        if(currChemicalBond == aChemicalBond) {
            tempBonds.killCurrent();
            break;
        }
    }
    
    Atom* firstAtom = aChemicalBond->getFirstAtom();

    rg_dList<Atom*> atomsQueue;
    
    atomsQueue.pushBack(firstAtom);
    
    while( atomsQueue.getSize() !=0 ) {
        Atom* currAtom = atomsQueue.popFront();

        rg_dList<ChemicalBond*> incidentChemicalBonds;

        getIncidentChemicalBondsAndKillIncidentBonds( currAtom, tempBonds, incidentChemicalBonds );

        incidentChemicalBonds.reset4Loop();
        while(incidentChemicalBonds.setNext4Loop()) {
            ChemicalBond* currBond = incidentChemicalBonds.getEntity();
            
            Atom* firstAtom  = currBond->getFirstAtom();
            Atom* secondAtom = currBond->getSecondAtom();

            if(firstAtom == currAtom) {
                atomsQueue.pushBack(secondAtom);
            }
            else if(secondAtom == currAtom) {
                atomsQueue.pushBack(firstAtom);
            }
        }
    }

    if(tempBonds.getSize() == 0 ) {
        isOnRing = rg_TRUE;
    }

    return isOnRing;

}



void    V::GeometryTier::getIncidentChemicalBondsAndKillIncidentBonds( 
                                V::GeometryTier::Atom* atom, 
                                rg_dList<V::GeometryTier::ChemicalBond*>& allChemicalBonds, 
                                rg_dList<V::GeometryTier::ChemicalBond*>& incidentChemicalBonds )
{
    incidentChemicalBonds.removeAll();

    allChemicalBonds.reset4Loop();
    while(allChemicalBonds.setNext4Loop()) {
        ChemicalBond* currBond = allChemicalBonds.getEntity();
        Atom* firstAtom  = currBond->getFirstAtom();
        Atom* secondAtom = currBond->getSecondAtom();
        
        if(firstAtom == atom) {
            incidentChemicalBonds.add(currBond);
            allChemicalBonds.killCurrent();   
        }
        else if(secondAtom == atom) {
            incidentChemicalBonds.add(currBond);
            allChemicalBonds.killCurrent();   
        }
    }    
}



void    V::GeometryTier::getAtomSetConnectedWithFirstAtomOnGivenBond( 
                                Molecule* aMolecule, 
                                ChemicalBond* givenBond, 
                                rg_dList<V::GeometryTier::Atom*>& atomSetConnectedWithFirstAtom)
{
    rg_BOOL isOnRing = rg_FALSE;
    rg_dList<ChemicalBond*> tempBonds;
    
    rg_dList<ChemicalBond>* allChemicalBonds = aMolecule->getChemicalBonds();

    allChemicalBonds->reset4Loop();
    while(allChemicalBonds->setNext4Loop()) {
        ChemicalBond* currChemicalBond = allChemicalBonds->getpEntity();
        
        tempBonds.add(currChemicalBond);
    }
    
    tempBonds.reset4Loop();
    while(tempBonds.setNext4Loop()) {
        ChemicalBond* currChemicalBond = tempBonds.getEntity();
        
        if(currChemicalBond == givenBond) {
            tempBonds.killCurrent();
            break;
        }
    }
    
    Atom* firstAtom = givenBond->getFirstAtom();
    
    rg_dList<Atom*> atomsQueue;
    
    atomsQueue.pushBack(firstAtom);
    
    while( atomsQueue.getSize() !=0 ) {
        Atom* currAtom = atomsQueue.popFront();

        atomSetConnectedWithFirstAtom.add(currAtom);
        
        rg_dList<ChemicalBond*> incidentChemicalBonds;
        
        getIncidentChemicalBondsAndKillIncidentBonds( currAtom, tempBonds, incidentChemicalBonds );
        
        incidentChemicalBonds.reset4Loop();
        while(incidentChemicalBonds.setNext4Loop()) {
            ChemicalBond* currBond = incidentChemicalBonds.getEntity();
            
            Atom* firstAtom  = currBond->getFirstAtom();
            Atom* secondAtom = currBond->getSecondAtom();
            
            if(firstAtom == currAtom) {
                atomsQueue.pushBack(secondAtom);
            }
            else if(secondAtom == currAtom) {
                atomsQueue.pushBack(firstAtom);
            }
        }
    }
}



void    V::GeometryTier::rotateAtomsByGivenBond( 
                                Molecule* aMolecule, 
                                ChemicalBond* givenBond, 
                                const double& angleByRadian )
{
    rg_dList<Atom*> atomSetConnectedWithFirstAtom;

    getAtomSetConnectedWithFirstAtomOnGivenBond( aMolecule, givenBond, atomSetConnectedWithFirstAtom);

    rg_Point3D startPt = givenBond->getSecondAtom()->getAtomBall().getCenter();
    rg_Point3D endPt   = givenBond->getFirstAtom()->getAtomBall().getCenter();

    rg_Point3D axis    = endPt - startPt;

    atomSetConnectedWithFirstAtom.reset4Loop();
    while(atomSetConnectedWithFirstAtom.setNext4Loop()) {
        Atom* currAtom = atomSetConnectedWithFirstAtom.getEntity();

        rg_Point3D currPt = currAtom->getAtomBall().getCenter();

        rg_Point3D rotatedPt = currPt - startPt;

        rg_TMatrix3D rotationMat;

        rotationMat.rotateArbitraryAxis(axis, angleByRadian);

        rotatedPt = rotationMat* rotatedPt;

        rotatedPt = rotatedPt + startPt;

        currAtom->setCenterOfAtomBall(rotatedPt);
    }
}


