#include "rg_TMatrix3D.h"
#include "rg_GeoFunc.h"
#include "FunctionsForRotamerLibrary.h"
#include "RotamerSetOfResidue.h"
#include "FunctionsForMolecule.h"
#include "rg_MathFunc.h"
#include <iostream>
#include <fstream>
typedef rg_Point3D Vector3D;

BBDepRotamer         FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[ SIZE_OF_BBDEP_ROTAMER_LIB ];
BBDepRotamer         FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB_2010[SIZE_OF_BBDEP_ROTAMER_LIB_2010];
//BBDepRotamer*         FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB = rg_NULL;
BBIndRotamer         FunctionsForRotamerLibrary::BBINDEP_ROTAMER_LIB[ SIZE_OF_BBINDEP_ROTAMER_LIB ];
TypeOfRotamerLibrary FunctionsForRotamerLibrary::m_currentlyUsedRotamerLibType = UNK_ROT_LIB;
rg_BOOL              FunctionsForRotamerLibrary::m_bBBINDLibLoaded  = rg_FALSE;
rg_BOOL              FunctionsForRotamerLibrary::m_bBBDEPLibLoaded  = rg_FALSE;


rg_BOOL FunctionsForRotamerLibrary::getAngleBinIndex_old(const rg_REAL& angle, rg_INDEX& index)
{
	if( (-180 <= angle) && (angle < -175) )
	{
		index = 0;
		return rg_TRUE;
	}
	else if( (-175 <= angle) && (angle < 175) )
	{
		rg_INT numBins   = 35;   // 37 - 1 [-180,-175) - 1 [175, 180)
		rg_REAL maxAngle =  175;
		rg_REAL minAngle = -175;

		const rg_INDEX offset = 1;
		index = offset + (rg_INDEX) ( numBins * (angle - minAngle) / (maxAngle - minAngle) );

		return rg_TRUE;
	}
	else if( (175 <= angle) && (angle < 180) )
	{
		index = 36;
		return rg_TRUE;
	}
	else
	{
		index = -1;
		return rg_FALSE;
	}
}

rg_BOOL FunctionsForRotamerLibrary::getAngleBinIndex(const rg_REAL& angle, rg_INDEX& index)
{
    rg_INT angle_dividedBy_GridSize = rg_MathFunc::round( angle / GRID_SIZE_OF_BACKBONE_DIHEDRAL_ANGLESPACE_FOR_BBDEPLIB );
    index = angle_dividedBy_GridSize + NUM_ANGLE_BINS_FOR_BBDEPLIB;
    if((0 <= index) && (index <= 36))
        return rg_TRUE;
    else
        return rg_FALSE;
}

rg_INDEX FunctionsForRotamerLibrary::getIndexOfResCode(const ResidueCode& code)
{
	rg_INDEX indexOfResidue = -1;
	switch(code)
	{
	case ARG_AMINO_RESIDUE:
		indexOfResidue	= 0;
		break;
	case ASN_AMINO_RESIDUE:
		indexOfResidue	= 1;
		break;
	case ASP_AMINO_RESIDUE:
		indexOfResidue	= 2;
		break;
	case CYS_AMINO_RESIDUE:
		indexOfResidue	= 3;
		break;	
	case GLN_AMINO_RESIDUE:
		indexOfResidue	= 4;
		break;	
	case GLU_AMINO_RESIDUE:
		indexOfResidue	= 5;
		break;	
	case HIS_AMINO_RESIDUE:
		indexOfResidue	= 6;
		break;	
	case ILE_AMINO_RESIDUE:
		indexOfResidue	= 7;
		break;	
	case LEU_AMINO_RESIDUE:
		indexOfResidue	= 8;
		break;	
	case LYS_AMINO_RESIDUE:
		indexOfResidue	= 9;
		break;	
	case MET_AMINO_RESIDUE:
		indexOfResidue	= 10;
		break;
	case PHE_AMINO_RESIDUE:
		indexOfResidue	= 11;
		break;	
	case PRO_AMINO_RESIDUE:
		indexOfResidue	= 12;
		break;	
	case SER_AMINO_RESIDUE:
		indexOfResidue	= 13;
		break;	
	case THR_AMINO_RESIDUE:
		indexOfResidue	= 14;
		break;	
	case TRP_AMINO_RESIDUE:
		indexOfResidue	= 15;
		break;	
	case TYR_AMINO_RESIDUE:
		indexOfResidue	= 16;
		break;	
	case VAL_AMINO_RESIDUE:
		indexOfResidue	= 17;
		break;
	//case ALA_AMINO_RESIDUE:
	//	indexOfResidue = 18;
	//	break;
	//case GLY_AMINO_RESIDUE:
	//	indexOfResidue = 19;
	//	break;
	default:
		break;
	}

	return indexOfResidue;
}

void FunctionsForRotamerLibrary::getRotamerSetsOfResidueSetFromBBDepRotLib(Backbone             & backbone,
	                                                                       Rotamer*               sortedAssignedRotamersOfProtein, 
	                                                                       RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
	                                                                       ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP             )
{
	// get backbone conformation
	BackboneConformation backboneConformation;
	backbone.getBackboneConformation(backboneConformation);

	rg_dList<BackboneDihedralAnglePair>* backboneDihedralAnglePairs 
		= backboneConformation.getBackboneDihedralAnglePairs();
	rg_INT numResiduePosition = backboneDihedralAnglePairs->getSize();

	// get sequence
	rg_dList<ResidueCode> sequence;
	backbone.getSequence(sequence);
		
	getRotamerSetsForResiduesSortedByChainIDSeqNumber(sortedAssignedRotamersOfProtein, 
		                                              numResiduePosition, 
		                                              sortedRotamerSetsOfResidues     );
	// reserve size for rotamer sets
	managerOfRotamerSetsForSCP.setNumResidues( numResiduePosition );

	// get rotamer sets considering the backbone dihedral angles
	rg_INDEX i = 0;
	backboneDihedralAnglePairs->reset4Loop();
	sequence.reset4Loop();
	while (backboneDihedralAnglePairs->setNext4Loop() && sequence.setNext4Loop())
	{
		rg_REAL dihedralAnglePair[ 2 ] = {DEFAULT_PHI_N_TERMINAL, DEFAULT_PSI_C_TERMINAL};
		BackboneDihedralAnglePair* currAnglePair = backboneDihedralAnglePairs->getpEntity();

#ifdef _DEBUG
        // test
        if(currAnglePair->getChainIDSeqNum() == 100187)
            int here = 1;
#endif

		currAnglePair->getDihedralAnglePair(dihedralAnglePair);
		ResidueCode currCode = sequence.getEntity();

		rg_INDEX start = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		rg_INDEX end   = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;

		getBBDepLibIndexOfResidue(currCode,
			                      dihedralAnglePair[ 0 ],
			                      dihedralAnglePair[ 1 ],
			                      start,
			                      end);

		sortedRotamerSetsOfResidues[ i ].setRotLibIDs(start, end);
		managerOfRotamerSetsForSCP.setRotamerSet(&sortedRotamerSetsOfResidues[ i ], i);
		i++;
	}
}

void FunctionsForRotamerLibrary::getRotamerSetsOfResidueSetFromBBDepRotLib_Dunbrack2010(Backbone             & backbone,
	                                                                                    Rotamer*               sortedAssignedRotamersOfProtein, 
	                                                                                    RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
	                                                                                    ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP             )
{
	// get backbone conformation
	BackboneConformation backboneConformation;
	backbone.getBackboneConformation(backboneConformation);

	rg_dList<BackboneDihedralAnglePair>* backboneDihedralAnglePairs 
		= backboneConformation.getBackboneDihedralAnglePairs();
	rg_INT numResiduePosition = backboneDihedralAnglePairs->getSize();

	// get sequence
	rg_dList<ResidueCode> sequence;
	backbone.getSequence(sequence);
		
	getRotamerSetsForResiduesSortedByChainIDSeqNumber(sortedAssignedRotamersOfProtein, 
		                                              numResiduePosition, 
		                                              sortedRotamerSetsOfResidues     );
	// reserve size for rotamer sets
	managerOfRotamerSetsForSCP.setNumResidues( numResiduePosition );

	// get rotamer sets considering the backbone dihedral angles
	rg_INDEX i = 0;
	backboneDihedralAnglePairs->reset4Loop();
	sequence.reset4Loop();
	while (backboneDihedralAnglePairs->setNext4Loop() && sequence.setNext4Loop())
	{
		rg_REAL dihedralAnglePair[ 2 ] = {DEFAULT_PHI_N_TERMINAL, DEFAULT_PSI_C_TERMINAL};
		BackboneDihedralAnglePair* currAnglePair = backboneDihedralAnglePairs->getpEntity();

		currAnglePair->getDihedralAnglePair(dihedralAnglePair);
		ResidueCode currCode = sequence.getEntity();

		rg_INDEX start = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		rg_INDEX end   = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;

		getBBDepLibIndexOfResidue_2010(currCode,
			                         dihedralAnglePair[ 0 ],
			                         dihedralAnglePair[ 1 ],
			                         start,
			                         end);

		sortedRotamerSetsOfResidues[ i ].setRotLibIDs(start, end);
		managerOfRotamerSetsForSCP.setRotamerSet(&sortedRotamerSetsOfResidues[ i ], i);
		i++;
	}
}

rg_BOOL FunctionsForRotamerLibrary::getBBDepLibIndexOfResidue(const ResidueCode& code, 
	                                                          rg_INDEX& start, 
															  rg_INDEX& end)
{
	rg_INDEX index = getIndexOfResCode( code );
	if( (0 <= index) && (index <= 17) )
	{
		start = INDEX_PAIR_BBDEP_ROTAMERS[ index ].m_startIndex;
		end   = INDEX_PAIR_BBDEP_ROTAMERS[ index ].m_endIndex  ;
		return rg_TRUE;
	}
	else
	{
		start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		return rg_FALSE;
	}
}

rg_BOOL FunctionsForRotamerLibrary::getBBDepLibIndexOfResidue(const ResidueCode& code, 
	                                                          const rg_REAL& phi, 
	                                                          const rg_REAL& psi, 
	                                                               rg_INDEX& start, 
	                                                               rg_INDEX& end)
{
	rg_INDEX residueIndex = getIndexOfResCode( code );
	if( (0 <= residueIndex) && (residueIndex <= 17) )
	{
		// find each index of bin which contains the backbone dihedral angle (phi, psi) using pre-defined bucket
		rg_INDEX phiIndex = -1;
		rg_INDEX psiIndex = -1;
		if(getAngleBinIndex(phi, phiIndex) && getAngleBinIndex(psi, psiIndex))
		{
			start = INDEX_PAIR_BUCKET_FOR_BBDEP_ROTAMERS[ residueIndex ][ phiIndex ][ psiIndex ].m_startIndex;
			end   = INDEX_PAIR_BUCKET_FOR_BBDEP_ROTAMERS[ residueIndex ][ phiIndex ][ psiIndex ].m_endIndex;
			return rg_TRUE;
		}
		else
		{
			start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
			return rg_FALSE;
		}		
	}
	else
	{
		start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		return rg_FALSE;
	}
}

rg_INT  FunctionsForRotamerLibrary::getNumBBDepRotamersOfResidue(const ResidueCode& code, 
	                                                             const rg_REAL& phi, 
	                                                             const rg_REAL& psi)
{
	rg_INDEX start = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
	rg_INDEX end   = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
	rg_BOOL bExist = getBBDepLibIndexOfResidue(code, phi, psi, start, end);

	if( bExist )
	{
		return (end - start + 1);
	}
	else
	{
		return 0;
	}	
}

rg_INT  FunctionsForRotamerLibrary::getNumBBDepRotamersOfResidue(const ResidueCode& code)
{
	rg_INDEX index = getIndexOfResCode( code );
	if( (0 <= index) && (index <= 17) )
	{
		return NUM_BBDEP_ROTAMERS[ index ];
	}
	else
	{
		return 0;
	}
}

rg_INDEX FunctionsForRotamerLibrary::getBBDepLibIndexOfRotamerWithHighestProbability(const ResidueCode& code, 
	                                                                                 const rg_REAL& phi, 
	                                                                                 const rg_REAL& psi)
{
	rg_INDEX start = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
	rg_INDEX end   = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
	rg_BOOL bExist = getBBDepLibIndexOfResidue(code, phi, psi, start, end);

	if( bExist )
	{
		return start;
	}
	else
	{
		return FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
	}
}

void  FunctionsForRotamerLibrary::getRotamerSetsOfResidueSetFromBBIndepRotLib(Backbone             & backbone,
	                                                                          Rotamer*               sortedAssignedRotamersOfProtein, 
	                                                                          RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
	                                                                          ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP   )
{
	// get sequence
	rg_dList<ResidueCode> sequence;
	backbone.getSequence(sequence);

	rg_INT numResiduePosition = -1;
	numResiduePosition = sequence.getSize();

	getRotamerSetsForResiduesSortedByChainIDSeqNumber(sortedAssignedRotamersOfProtein, 
		                                              numResiduePosition, 
		                                              sortedRotamerSetsOfResidues     );
	// reserve size for rotamer sets
	managerOfRotamerSetsForSCP.setNumResidues( numResiduePosition );

	// get rotamer sets considering the sequence
	rg_INDEX i = 0;
	sequence.reset4Loop();
	while (sequence.setNext4Loop())
	{
		ResidueCode currCode = sequence.getEntity();

		rg_INDEX start = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		rg_INDEX end   = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;

		getBBIndepLibIndexOfResidue(currCode,
			                        start,
			                        end);

		sortedRotamerSetsOfResidues[ i ].setRotLibIDs(start, end);
		managerOfRotamerSetsForSCP.setRotamerSet(&sortedRotamerSetsOfResidues[ i ], i);
		i++;
	}
}

rg_BOOL FunctionsForRotamerLibrary::getBBIndepLibIndexOfResidue(const ResidueCode& code, 
	                                                            rg_INDEX& start, 
	                                                            rg_INDEX& end)
{
	rg_INDEX residueIndex = getIndexOfResCode( code );
	if( (0 <= residueIndex) && (residueIndex <= 17) )
	{
		start = INDEX_PAIR_BBINDEP_ROTAMERS[ residueIndex ].m_startIndex;
		end   = INDEX_PAIR_BBINDEP_ROTAMERS[ residueIndex ].m_endIndex;
		return rg_TRUE;
	}
	else
	{
		start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		return rg_FALSE;
	}
}

//TypeOfRotamerLibrary FunctionsForRotamerLibrary::getTypeOfRotamerLibrary()
//{
//	return m_currentlyUsedRotamerLibType;
//}

rg_INT FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(const ResidueCode& code)
{
#ifdef _DEBUG
	if(code == PRO_AMINO_RESIDUE)
		int here = 1;
#endif

	if(code == ALA_AMINO_RESIDUE || code == GLY_AMINO_RESIDUE)
	{
		return 0;
	}
	else
	{
		return NUM_DIHEDRAL_ANGLES[getIndexOfResCode(code)];
	}

	//if(code == SER_AMINO_RESIDUE || code == THR_AMINO_RESIDUE ||
	//   code == VAL_AMINO_RESIDUE || code == CYS_AMINO_RESIDUE  )
	//{
	//	return 1;
	//}
	//else if(code == ASP_AMINO_RESIDUE || code == LEU_AMINO_RESIDUE || code == ILE_AMINO_RESIDUE ||
	//	    code == ASN_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE ||
	//		code == PRO_AMINO_RESIDUE || code == HIS_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE   )
	//{
	//		return 2;
	//}
	//else if(code == LYS_AMINO_RESIDUE || code == ARG_AMINO_RESIDUE)
	//{
	//	return 4;
	//}
	//else if(code == GLU_AMINO_RESIDUE ||
	//	    code == GLN_AMINO_RESIDUE ||
	//	    code == MET_AMINO_RESIDUE  )
	//{
	//	return 3;
	//}
	////if(code == ALA_AMINO_RESIDUE || code == GLY_AMINO_RESIDUE)
	//else
	//{
	//	return 0;
	//}
}

rg_INT  FunctionsForRotamerLibrary::getDihedralAngles(const rg_INDEX& rotLidIndex, rg_REAL*& dihedralAngles)
{
	rg_INT numDihedralAngles = 0;

	switch (m_currentlyUsedRotamerLibType)
	{
	case BACKBONE_INDEPENDENT_LIB_DUNBRACK2002:
		{
			numDihedralAngles = getNumDihedralAnglesOfResidue( BBINDEP_ROTAMER_LIB[ rotLidIndex ].code );
			rg_INDEX i;
			for(i = 0;i < numDihedralAngles;i++)
			{
				dihedralAngles[ i ] = BBINDEP_ROTAMER_LIB[ rotLidIndex ].chi[ i ];
			}
		}
		break;
	case BACKBONE_DEPENDENT_LIB_DUNBRACK2002:
		{
			numDihedralAngles = getNumDihedralAnglesOfResidue( BBDEP_ROTAMER_LIB[ rotLidIndex ].code );
			rg_INDEX i;
			for(i = 0;i < numDihedralAngles;i++)
			{
				dihedralAngles[ i ] = BBDEP_ROTAMER_LIB[ rotLidIndex ].chi[ i ];
			}
		}
		break;
	case BACKBONE_DEPENDENT_LIB_DUNBRACK2010:
		{
			//Coded by JHCHa - Apr. 04. 2016.
			numDihedralAngles = getNumDihedralAnglesOfResidue_2010(BBDEP_ROTAMER_LIB_2010[rotLidIndex].code);
			rg_INDEX i;
			for (i = 0; i < numDihedralAngles; i++)
			{
				dihedralAngles[i] = BBDEP_ROTAMER_LIB_2010[rotLidIndex].chi[i];
			}
		}
		break;
	default:
		break;
	}

	return numDihedralAngles;
}

void FunctionsForRotamerLibrary::transformAtomCentersOfSidechainUsingDihedralAngle(Rotamer& rotamer, const rg_INT& rotamerIndex)
{
	
	rg_INT numDihedralAngles = rotamer.getNumDihedralAngles();

	rg_REAL* dihedralAngles = rg_NULL;
	dihedralAngles = new rg_REAL[numDihedralAngles];
	FunctionsForRotamerLibrary::getDihedralAngles(rotamerIndex, dihedralAngles);

	Residue* targetResidue = rotamer.getResidue();
	transformAtomCentersOfSidechainUsingDihedralAngle(targetResidue, dihedralAngles, numDihedralAngles);
	rotamer.setDihedralAngles( dihedralAngles );

	if(dihedralAngles != rg_NULL)
		delete [] dihedralAngles;
}

void FunctionsForRotamerLibrary::transformAtomCentersOfSidechainUsingDihedralAngle(Rotamer& rotamer, const rg_REAL* dihedralAngles)
{
	Residue* targetResidue = rotamer.getResidue();
	transformAtomCentersOfSidechainUsingDihedralAngle(targetResidue, dihedralAngles);
	rotamer.setDihedralAngles( dihedralAngles );
}

void FunctionsForRotamerLibrary::transformAtomCentersOfSidechainUsingDihedralAngle(Residue* targetResidue, const rg_INT& rotamerIndex)
{
	rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(targetResidue->getResidueCode());

	rg_REAL* dihedralAngles = rg_NULL;
	dihedralAngles = new rg_REAL[numDihedralAngles];
	FunctionsForRotamerLibrary::getDihedralAngles(rotamerIndex, dihedralAngles);

	transformAtomCentersOfSidechainUsingDihedralAngle(targetResidue, dihedralAngles, numDihedralAngles);

	if(dihedralAngles != rg_NULL)
		delete [] dihedralAngles;
}

void FunctionsForRotamerLibrary::transformAtomCentersOfSidechainUsingDihedralAngle(Residue* targetResidue, const rg_REAL* dihedralAngles)
{
	rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(targetResidue->getResidueCode());
	transformAtomCentersOfSidechainUsingDihedralAngle(targetResidue, dihedralAngles, numDihedralAngles);
}

void FunctionsForRotamerLibrary::transformAtomCentersOfSidechainUsingDihedralAngle(Residue* targetResidue, const rg_REAL* dihedralAngles, const rg_INT& numDihedralAngles)
{
#ifdef _DEBUG
    ResidueCode code = targetResidue->getResidueCode();
    if(code == HIS_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE)
        int here = 1;
#endif


#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == ARG_AMINO_RESIDUE)
        int here = 1;
#endif

	// initialize the array of atom pointers which contains the consecutive four atom pointers on side chain
	Atom* consecutiveFourAtoms[ 4 ] 
	= {rg_NULL,
	   targetResidue->getNitrogenInAminoGroupOfAminoResidue(),
	   targetResidue->getAlphaCarbonOfAminoResidue(),
	   targetResidue->getBetaCarbonInSideChainOfAminoResidue()};

	rg_INT dihedralAngleIndex = -1;

	// get atoms on side chain in breath first order
	rg_dList<Atom*> atomsOnSideChainInBreathFirstOrder;
	RemoteIndicator  startRemote = GAMMA_REMOTE;
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	// getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsOnSideChainInBreathFirstOrder, startRemote);	
	targetResidue->getAtomsOnSidechain(atomsOnSideChainInBreathFirstOrder);
	atomsOnSideChainInBreathFirstOrder.reset4Loop();
	while (atomsOnSideChainInBreathFirstOrder.setNext4Loop())
	{
		Atom* currSidechainAtom = atomsOnSideChainInBreathFirstOrder.getEntity();
		RemoteIndicator remote = currSidechainAtom->getpChemicalProperties()->getRemoteIndicator();
		if(remote < startRemote)
			atomsOnSideChainInBreathFirstOrder.killCurrent();
		else
			break;
	}

	//rg_dList<Atom*> atomsOnSideChainInBreathFirstOrder;
	//RemoteIndicator  endRemote = UNK_REMOTE;
	//BranchDesignator endBranch = UNK_BRANCH;
	//FunctionsForRotamerLibrary::getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(targetResidue->getResidueCode(), endRemote, endBranch);
	//RemoteIndicator  startRemote = GAMMA_REMOTE;
	//getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsOnSideChainInBreathFirstOrder, startRemote, endRemote, endBranch);

	RemoteIndicator prevRemote  = BETA_REMOTE;
	atomsOnSideChainInBreathFirstOrder.reset4Loop();

	while(atomsOnSideChainInBreathFirstOrder.setNext4Loop())
	{
		Atom* currAtom = atomsOnSideChainInBreathFirstOrder.getEntity();
		RemoteIndicator currRemote = currAtom->getpChemicalProperties()->getRemoteIndicator();
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
				consecutiveFourAtoms[ i ] = consecutiveFourAtoms[i + 1];
			// construct the center of fourth atom
			consecutiveFourAtoms[ 3 ] = currAtom;
			// update atom coordinate one by one
			computeAtomCenterUsingDihedralAngle(consecutiveFourAtoms[ 0 ],
				                                consecutiveFourAtoms[ 1 ],
				                                consecutiveFourAtoms[ 2 ],
				                                consecutiveFourAtoms[ 3 ],
				                                dihedralAngles[dihedralAngleIndex]);
            //computeAtomCenterUsingDihedralAngle3(consecutiveFourAtoms[ 0 ],
            //                                     consecutiveFourAtoms[ 1 ],
            //                                     consecutiveFourAtoms[ 2 ],
            //                                     consecutiveFourAtoms[ 3 ],
            //                                     dihedralAngles[dihedralAngleIndex]);

			// test
			//computeAtomCenterUsingDihedralAngle2(consecutiveFourAtoms[ 0 ],
			//	consecutiveFourAtoms[ 1 ],
			//	consecutiveFourAtoms[ 2 ],
			//	consecutiveFourAtoms[ 3 ],
			//	dihedralAngles[dihedralAngleIndex]);

			prevRemote = currRemote;
		}
		// same remote yet different branch
		// The case (prevRemote == currRemote) is only possible
		else
		{
			// update second branch atom
			//computeSecondBranchAtomCenter(consecutiveFourAtoms[ 1 ], 
			//	                             consecutiveFourAtoms[ 2 ], 
			//	                             consecutiveFourAtoms[ 3 ], 
			//	                             currAtom);			
			computeSecondBranchAtomCenter2(consecutiveFourAtoms[ 0 ],
				                              consecutiveFourAtoms[ 1 ], 
				                              consecutiveFourAtoms[ 2 ], 
				                              consecutiveFourAtoms[ 3 ], 
				                              currAtom);
		}
	}

#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == HIS_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE)
    {
        rg_REAL* angles = new rg_REAL[numDihedralAngles];
        targetResidue->getSideChainDihedralAngles(angles);
        if(angles != rg_NULL)
            delete [] angles;
    }
#endif

#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == ARG_AMINO_RESIDUE)
    {
        rg_REAL* angles = new rg_REAL[numDihedralAngles];
        targetResidue->getSideChainDihedralAngles(angles);
        if(angles != rg_NULL)
            delete [] angles;
    }
#endif
}

void FunctionsForRotamerLibrary::computeAtomCenterUsingDihedralAngle(Atom* atom0, 
	                                                                 Atom* atom1, 
	                                                                 Atom* atom2, 
	                                                                 Atom* atom3, 
	                                                                 const rg_REAL& dihedralAngle)
{
	rg_Point3D atomCenter[ 3 ] 
	= {atom0->getpAtomBall()->getCenter(), 
	   atom1->getpAtomBall()->getCenter(),
	   atom2->getpAtomBall()->getCenter()};

	// rotation via dihedral angle
	Vector3D rotationAxisForDihedralAngle = atomCenter[ 2 ] - atomCenter[ 1 ];
	rotationAxisForDihedralAngle.normalize();
	rg_TMatrix3D rotMat4DihedralAngle;
	rotMat4DihedralAngle.rotateArbitraryAxis(rotationAxisForDihedralAngle, rg_PI * dihedralAngle / 180.);
	Vector3D vec10 = atomCenter[ 0 ] - atomCenter[ 1 ];
	vec10.normalize();
	Vector3D dirVecCoplanarWithAtomCenter012 = rotMat4DihedralAngle * vec10;

	// get bond lengthes and angles from AMBER
	AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomType2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomType3 = atom3->getpChemicalProperties()->getAtomTypeInAmber();

	rg_REAL bondLength23 = BOND_LENGTH[atomType2][atomType3];
	rg_REAL bondAngle123 = rg_PI * BOND_ANGLE[atomType1][atomType2][atomType3] / 180.;

	// rotation via bond angle	
	Vector3D vec21 = (-rotationAxisForDihedralAngle);
	Vector3D rotationAxisForBondAngle = vec21.crossProduct(dirVecCoplanarWithAtomCenter012);
	rotationAxisForBondAngle.normalize();

	rg_TMatrix3D rotMat4BondAngle;
	rotMat4BondAngle.rotateArbitraryAxis(rotationAxisForBondAngle, bondAngle123);
	vec21 = rotMat4BondAngle * vec21;

	rg_Point3D updatedCenter = atomCenter[ 2 ] + vec21 * bondLength23;	
	atom3->setCenterOfAtomBall(updatedCenter);
}

void FunctionsForRotamerLibrary::computeAtomCenterUsingDihedralAngle2(Atom* atom0, 
	                                                                  Atom* atom1, 
	                                                                  Atom* atom2, 
	                                                                  Atom* atom3, 
	                                                                  const rg_REAL& dihedralAngle)
{
	
	rg_Point3D atomCenter[ 3 ] 
	= {atom0->getpAtomBall()->getCenter(), 
	   atom1->getpAtomBall()->getCenter(),
	   atom2->getpAtomBall()->getCenter()};

	// get bond lengths and angles from AMBER
	AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomType2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomType3 = atom3->getpChemicalProperties()->getAtomTypeInAmber();
	rg_REAL bondLength23 = BOND_LENGTH[atomType2][atomType3];
	rg_REAL bondAngle123 = BOND_ANGLE[atomType1][atomType2][atomType3];

	// translate atom center2 
	Vector3D vec12 = atomCenter[ 2 ] - atomCenter[ 1 ];
	vec12.normalize();
	rg_Point3D atomCenter3 = atomCenter[ 2 ] + bondLength23 * vec12;

	// rotation via bond angle
	Vector3D vec10 = atomCenter[ 0 ] - atomCenter[ 1 ];
	vec10.normalize();
	Vector3D normalOfPlane012 = vec10.crossProduct( vec12 );
	rg_TMatrix3D tMat4RotByBondAngle;
	rg_REAL rotAngle = (180.0 - bondAngle123) * rg_PI / 180.0;
	tMat4RotByBondAngle.rotateArbitraryAxis(normalOfPlane012, rotAngle);

	Vector3D vec23 = atomCenter3 - atomCenter[ 2 ];
	vec23.normalize();
	vec23 = tMat4RotByBondAngle * vec23;
	atomCenter3 = atomCenter[ 2 ] + bondLength23 * vec23;

	// rotation via dihedral angle
	rg_TMatrix3D tMat4RotByDihedralAngle;
	rotAngle = rg_PI * dihedralAngle / 180.0;
	tMat4RotByDihedralAngle.rotateArbitraryAxis(vec12, rotAngle);
	vec23 = atomCenter3 - atomCenter[ 2 ];
	vec23.normalize();
	vec23 = tMat4RotByDihedralAngle * vec23;
	atomCenter3 = atomCenter[ 2 ] + bondLength23 * vec23;

	// set new coordinate
	atom3->setCenterOfAtomBall(atomCenter3);
}

void FunctionsForRotamerLibrary::computeAtomCenterUsingDihedralAngle3(Atom* atom0, 
                                                                      Atom* atom1, 
                                                                      Atom* atom2, 
                                                                      Atom* atom3, 
                                                                      const rg_REAL& dihedralAngle)
{
    // get bond lengthes and angles from AMBER
    AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
    AmberAtomTypes atomType2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();
    AmberAtomTypes atomType3 = atom3->getpChemicalProperties()->getAtomTypeInAmber();

    rg_REAL bondLength23 = BOND_LENGTH[atomType2][atomType3];
    rg_REAL bondAngle123 = rg_PI * BOND_ANGLE[atomType1][atomType2][atomType3] / 180.;

    rg_Point3D atomCenter[ 3 ] 
    = {atom0->getpAtomBall()->getCenter(), 
       atom1->getpAtomBall()->getCenter(),
       atom2->getpAtomBall()->getCenter()};

    // translate 
    Vector3D transVec = atomCenter[ 2 ] - atomCenter[ 1 ];
    transVec.normalize();
    rg_Point3D updatedCenter = atomCenter[ 2 ] + transVec * bondLength23;

    // rotation via bond angle	
    Vector3D vec10 = atomCenter[ 0 ] - atomCenter[ 1 ];
    vec10.normalize();
    Vector3D vec12 = atomCenter[ 2 ] - atomCenter[ 1 ];
    vec12.normalize();
    Vector3D rotationAxisForBondAngle = vec12.crossProduct(vec10);
    rotationAxisForBondAngle.normalize();

    rg_TMatrix3D rotMat4BondAngle;
    rotMat4BondAngle.rotateArbitraryAxis(rotationAxisForBondAngle, rg_PI * (180.0 - bondAngle123) / 180.0);
    Vector3D vec23 = updatedCenter - atomCenter[ 2 ];
    vec23.normalize();
    vec23 = rotMat4BondAngle * vec23;

    // rotation via dihedral angle
    Vector3D rotationAxisForDihedralAngle = vec12;
    rg_TMatrix3D rotMat4DihedralAngle;
    rotMat4DihedralAngle.rotateArbitraryAxis(rotationAxisForDihedralAngle, rg_PI * dihedralAngle / 180.);
    vec23 = rotMat4DihedralAngle * vec23;

    updatedCenter = atomCenter[ 2 ] + vec23 * bondLength23;	
    atom3->setCenterOfAtomBall(updatedCenter);
}


void FunctionsForRotamerLibrary::computePointUsingDihedralAngle(const rg_Point3D& point0, 
										                        const rg_Point3D& point1,
										                        const rg_Point3D& point2,
										                        const rg_REAL   & bondLength23,
										                        const rg_REAL   & bondAngle123,
										                        const rg_REAL   & dihedralAngle,
										                        rg_Point3D      & point3)
{
	// rotation via dihedral angle
	Vector3D rotationAxisForDihedralAngle = point2 - point1;
	rotationAxisForDihedralAngle.normalize();
	rg_TMatrix3D rotMat4DihedralAngle;
	rotMat4DihedralAngle.rotateArbitraryAxis(rotationAxisForDihedralAngle, rg_PI * dihedralAngle / 180.);
	Vector3D vec10 = point0 - point1;
	vec10.normalize();
	Vector3D dirVecCoplanarWithAtomCenter012 = rotMat4DihedralAngle * vec10;

	// rotation via bond angle	
	Vector3D vec21 = (-rotationAxisForDihedralAngle);
	Vector3D rotationAxisForBondAngle = vec21.crossProduct(dirVecCoplanarWithAtomCenter012);
	rotationAxisForBondAngle.normalize();

	rg_TMatrix3D rotMat4BondAngle;
	rotMat4BondAngle.rotateArbitraryAxis(rotationAxisForBondAngle, bondAngle123);
	vec21 = rotMat4BondAngle * vec21;

	point3 = point2 + vec21 * bondLength23;
}

void FunctionsForRotamerLibrary::computeSecondBranchAtomCenter(Atom* atom0, 
	                                                           Atom* atom1, 
															   Atom* atom2, 
															   Atom* secondBranchAtom)
{
	// get bond lengthes and angles from AMBER
	AmberAtomTypes atomType0 = atom0->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes atomType2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();
	AmberAtomTypes secondBranchAtomType = secondBranchAtom->getpChemicalProperties()->getAtomTypeInAmber();

	rg_REAL bondLength[ 3 ] = {BOND_LENGTH[atomType1][atomType0],
		                       BOND_LENGTH[atomType1][atomType2],
		                       BOND_LENGTH[atomType1][secondBranchAtomType]};

	rg_REAL bondAngle[ 2 ] = {rg_PI * BOND_ANGLE[atomType0][atomType1][secondBranchAtomType] / 180.,
		                      rg_PI * BOND_ANGLE[atomType2][atomType1][secondBranchAtomType] / 180.};

	// construct the center for second branch atom
	// via computing the fourth vertex of a tetrahedron given three vertices and their edges
	// This could be solved by computing the intersection among three spheres corresponding to three atoms
	rg_REAL edgeLength[ 3 ] = 
	{rg_GeoFunc::computeLengthOfTriangleGivenTwoLengthesNAngle(bondLength[ 0 ], bondLength[ 2 ], bondAngle[ 0 ]), 
	bondLength[ 2 ], 
	rg_GeoFunc::computeLengthOfTriangleGivenTwoLengthesNAngle(bondLength[ 1 ], bondLength[ 2 ], bondAngle[ 1 ])};

	rg_Point3D centerOfSecondBranch;
	rg_Point3D centerOfNearbyThreeAtoms[ 3 ] 
	= { atom0->getpAtomBall()->getCenter(),
		atom1->getpAtomBall()->getCenter(),
		atom2->getpAtomBall()->getCenter() };

	rg_Point3D* intersectionPoint = new rg_Point3D[ 2 ];
	rg_GeoFunc::computeIntersectionPointAmongThreeSpheres(Sphere(centerOfNearbyThreeAtoms[ 0 ], edgeLength[ 0 ]),
		                                                  Sphere(centerOfNearbyThreeAtoms[ 1 ], edgeLength[ 1 ]),
		                                                  Sphere(centerOfNearbyThreeAtoms[ 2 ], edgeLength[ 2 ]),
		                                                  intersectionPoint);

	rg_REAL dist0 = centerOfNearbyThreeAtoms[ 2 ].distance(intersectionPoint[ 0 ]);
	rg_REAL dist1 = centerOfNearbyThreeAtoms[ 2 ].distance(intersectionPoint[ 1 ]);

	// choose the intersection whose distance from FIRST BRANCH atom is shorter
	if(dist0 > dist1)
		secondBranchAtom->setCenterOfAtomBall(intersectionPoint[ 0 ]);
	else
		secondBranchAtom->setCenterOfAtomBall(intersectionPoint[ 1 ]);

	if(intersectionPoint != rg_NULL)
		delete[] intersectionPoint;

}

void FunctionsForRotamerLibrary::computeSecondBranchAtomCenter2(Atom* atom0, 
	                                                            Atom* atom1, 
																Atom* atom2, 
																Atom* atom3, 
																Atom* secondBranchAtom)
{
	rg_Point3D* candidatesForAtomCenter = rg_NULL;

	computeTwoCandidatesForSecondBranchAtomCenter(atom1, 
		                                          atom2, 
												  atom3, 
												  secondBranchAtom, 
												  candidatesForAtomCenter);

	rg_Point3D centerOfPrevAtom = atom0->getpAtomBall()->getCenter();

	rg_REAL dist0 = centerOfPrevAtom.distance(candidatesForAtomCenter[ 0 ]);
	rg_REAL dist1 = centerOfPrevAtom.distance(candidatesForAtomCenter[ 1 ]);

	// choose the intersection whose distance from FIRST BRANCH atom is shorter
	if(dist0 > dist1)
		secondBranchAtom->setCenterOfAtomBall(candidatesForAtomCenter[ 0 ]);
	else
		secondBranchAtom->setCenterOfAtomBall(candidatesForAtomCenter[ 1 ]);

	if(candidatesForAtomCenter != rg_NULL)
		delete[] candidatesForAtomCenter;
}

 void FunctionsForRotamerLibrary::computeTwoCandidatesForSecondBranchAtomCenter(Atom* atom1, 
	                                                                            Atom* atom2, 
																				Atom* atom3, 
																				Atom* secondBranchAtom, 
																				rg_Point3D*& candidatesForAtomCenter)
 {
	 // get bond lengthes and angles from AMBER
	 AmberAtomTypes atomType1 = atom1->getpChemicalProperties()->getAtomTypeInAmber();
	 AmberAtomTypes atomType2 = atom2->getpChemicalProperties()->getAtomTypeInAmber();
	 AmberAtomTypes atomType3 = atom3->getpChemicalProperties()->getAtomTypeInAmber();
	 AmberAtomTypes secondBranchAtomType = secondBranchAtom->getpChemicalProperties()->getAtomTypeInAmber();

	 rg_REAL bondLength[ 3 ] = {BOND_LENGTH[atomType2][atomType1],
		 BOND_LENGTH[atomType2][atomType3],
		 BOND_LENGTH[atomType2][secondBranchAtomType]};

	 rg_REAL bondAngle[ 2 ] = {rg_PI * BOND_ANGLE[atomType1][atomType2][secondBranchAtomType] / 180.,
		 rg_PI * BOND_ANGLE[atomType3][atomType2][secondBranchAtomType] / 180.};

	 // construct the center for second branch atom
	 // via computing the fourth vertex of a tetrahedron given three vertices and their edges
	 // This could be solved by computing the intersection among three spheres corresponding to three atoms
	 rg_REAL edgeLength[ 3 ] = 
	 {rg_GeoFunc::computeLengthOfTriangleGivenTwoLengthesNAngle(bondLength[ 0 ], bondLength[ 2 ], bondAngle[ 0 ]), 
	 bondLength[ 2 ], 
	 rg_GeoFunc::computeLengthOfTriangleGivenTwoLengthesNAngle(bondLength[ 1 ], bondLength[ 2 ], bondAngle[ 1 ])};

	 rg_Point3D centerOfSecondBranch;
	 rg_Point3D centerOfNearbyThreeAtoms[ 3 ] 
	 = { atom1->getpAtomBall()->getCenter(),
		 atom2->getpAtomBall()->getCenter(),
		 atom3->getpAtomBall()->getCenter() };

	 candidatesForAtomCenter = new rg_Point3D[ 2 ];
	 rg_GeoFunc::computeIntersectionPointAmongThreeSpheres(Sphere(centerOfNearbyThreeAtoms[ 0 ], edgeLength[ 0 ]),
		                                                   Sphere(centerOfNearbyThreeAtoms[ 1 ], edgeLength[ 1 ]),
		                                                   Sphere(centerOfNearbyThreeAtoms[ 2 ], edgeLength[ 2 ]),
		                                                   candidatesForAtomCenter);

 }

void FunctionsForRotamerLibrary::glueFixedStructureIntoResidue(Residue* targetResidue)
{
	if(! hasFixedSubstructure(targetResidue) )
		return;

	ResidueCode code = targetResidue->getResidueCode();
	RemoteIndicator gluePoint;
	if(code == ARG_AMINO_RESIDUE)
		gluePoint = DELTA_REMOTE;
	else
		gluePoint = GAMMA_REMOTE;

	// old origin, old X-Axis, and old Z-Axis vectors from fixed substructure
	rg_Point3D oldOrigin(0.0, 0.0, 0.0);
	Vector3D   oldXAxisVec(1.0, 0.0, 0.0);
	Vector3D   oldZAxisVec(0.0, 0.0, 1.0);

	// get three reference points from substructure updated by using dihedral angles
	Atom* referenceAtom[ 3 ] = {rg_NULL, rg_NULL, rg_NULL};
	rg_dList<Atom*> atomsWithThisRemote;
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	//getAtomsSortedByBranchDesignator(gluePoint, targetResidue, atomsWithThisRemote);
	targetResidue->getAtoms(gluePoint, atomsWithThisRemote);
	referenceAtom[ 0 ] = atomsWithThisRemote.getFirstEntity();
	atomsWithThisRemote.removeAll();

	if(code == ARG_AMINO_RESIDUE)
	{
		//RemoteIndicator nextRemote = getNextRemoteIndicator(gluePoint);
		RemoteIndicator nextRemote = (RemoteIndicator)( (rg_INT)gluePoint + 1 );
		// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
		//getAtomsSortedByBranchDesignator(nextRemote, targetResidue, atomsWithThisRemote);
		targetResidue->getAtoms(nextRemote, atomsWithThisRemote);
		referenceAtom[ 1 ] = atomsWithThisRemote.getFirstEntity(); // get epsilon nitrogen
		atomsWithThisRemote.removeAll();

		nextRemote = (RemoteIndicator)( (rg_INT)nextRemote + 1 );		
		//getAtomsSortedByBranchDesignator(getNextRemoteIndicator(nextRemote), targetResidue, atomsWithThisRemote);
		//getAtomsSortedByBranchDesignator(nextRemote, targetResidue, atomsWithThisRemote);
		targetResidue->getAtoms(nextRemote, atomsWithThisRemote);
		referenceAtom[ 2 ] = atomsWithThisRemote.getFirstEntity(); // get zeta carbon
	}
    // Other four residues: HIS, PHE, TYR, and TRP
	else
	{
		RemoteIndicator nextRemote = (RemoteIndicator)( (rg_INT)gluePoint + 1 );
		//getAtomsSortedByBranchDesignator(getNextRemoteIndicator(gluePoint), targetResidue, atomsWithThisRemote);
		//getAtomsSortedByBranchDesignator(nextRemote, targetResidue, atomsWithThisRemote);
		targetResidue->getAtoms(nextRemote, atomsWithThisRemote);
        // get both first and second delta branches
		referenceAtom[ 1 ] = atomsWithThisRemote.getFirstEntity();
		referenceAtom[ 2 ] = atomsWithThisRemote.getSecondEntity();
	}

	// three reference points on updated substructure
	rg_Point3D referencePnt[ 3 ] = {referenceAtom[ 0 ]->getpAtomBall()->getCenter(),
		                            referenceAtom[ 1 ]->getpAtomBall()->getCenter(),
		                            referenceAtom[ 2 ]->getpAtomBall()->getCenter()};

	// new origin, new X-Axis, and new Z-Axis vectors from three reference points
	Vector3D vec01 = referencePnt[ 1 ] - referencePnt[ 0 ];
	Vector3D vec02 = referencePnt[ 2 ] - referencePnt[ 0 ];
	Vector3D normalVecToPlane = vec01.crossProduct(vec02);		

	rg_Point3D newOrigin   = referencePnt[ 0 ];
	Vector3D   newXAxisVec = vec01.getUnitVector();
	Vector3D   newZAxisVec = normalVecToPlane.getUnitVector();

	// construct transformation matrix for transforming the old local frame to new local frame
	rg_TMatrix3D transMat;
	transMat.convertAxis(oldOrigin, oldXAxisVec, oldZAxisVec, 
		                 newOrigin, newXAxisVec, newZAxisVec);

	// glue fixed substructure into residue

	// get the points on the fixed substructure
	rg_dList<rg_Point3D> pointsOnFixedStructure;
	getCoordinatesOfFixedStructure(code, pointsOnFixedStructure);
	// apply the transformation to the fixed substructure 
	pointsOnFixedStructure.reset4Loop();
	while(pointsOnFixedStructure.setNext4Loop())
	{
		rg_Point3D* currPoint = pointsOnFixedStructure.getpEntity();
		currPoint->setPoint(transMat * (*currPoint));
	}

	// get atom pointers corresponding to fixed substructure
	//RemoteIndicator currRemote = gluePoint;
	rg_dList<Atom*> atomsToBeTransformed;
	// We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
	//getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsToBeTransformed, gluePoint);
	targetResidue->getAtomsOnSidechain(atomsToBeTransformed);
	atomsToBeTransformed.reset4Loop();
	while (atomsToBeTransformed.setNext4Loop())
	{
		Atom* currSidechainAtom = atomsToBeTransformed.getEntity();
		RemoteIndicator remote = currSidechainAtom->getpChemicalProperties()->getRemoteIndicator();
		if(remote < gluePoint)
			atomsToBeTransformed.killCurrent();
		else
			break;
	}
	//RemoteIndicator  endRemote = UNK_REMOTE;
	//BranchDesignator endBranch = UNK_BRANCH;
	//FunctionsForRotamerLibrary::getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(targetResidue->getResidueCode(), endRemote, endBranch);
	//RemoteIndicator  startRemote = gluePoint;
	//getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsToBeTransformed, startRemote, endRemote, endBranch);

#ifdef _DEBUG
    // test
    if(code == TYR_AMINO_RESIDUE)
        int here = 1;
#endif

	//// update atom coordinates
	atomsToBeTransformed.reset4Loop();
	pointsOnFixedStructure.reset4Loop();
	while(atomsToBeTransformed.setNext4Loop() &&
		pointsOnFixedStructure.setNext4Loop() )
	{
		Atom* currAtom = atomsToBeTransformed.getEntity();
		rg_Point3D currUpdatedCenter = pointsOnFixedStructure.getEntity();
		currAtom->getpAtomBall()->setCenter(currUpdatedCenter);
	}


#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == HIS_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE)
    {
        rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(code);
        rg_REAL* angles = new rg_REAL[numDihedralAngles];
        targetResidue->getSideChainDihedralAngles(angles);
        if(angles != rg_NULL)
            delete [] angles;
    }
#endif

#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == ARG_AMINO_RESIDUE)
    {
        rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(code);
        rg_REAL* angles = new rg_REAL[numDihedralAngles];
        targetResidue->getSideChainDihedralAngles(angles);
        if(angles != rg_NULL)
            delete [] angles;
    }
#endif
}




void FunctionsForRotamerLibrary::glueFixedStructureIntoResidue2(Residue* targetResidue)
{
    if(! hasFixedSubstructure(targetResidue) )
        return;

    ResidueCode code = targetResidue->getResidueCode();
    RemoteIndicator remoteOfReferencePt;
    if(code == ARG_AMINO_RESIDUE)
        remoteOfReferencePt = DELTA_REMOTE;
    else
        remoteOfReferencePt = GAMMA_REMOTE;

    // old origin, old X-Axis, and old Z-Axis vectors from fixed substructure
    rg_Point3D oldOrigin(0.0, 0.0, 0.0);
    Vector3D   oldXAxisVec(1.0, 0.0, 0.0);
    Vector3D   oldZAxisVec(0.0, 0.0, 1.0);

    // get three reference points from substructure updated by using dihedral angles
    Atom* referenceAtom[ 3 ] = {rg_NULL, rg_NULL, rg_NULL};
    rg_dList<Atom*> atomsWithThisRemote;
    // We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
    //getAtomsSortedByBranchDesignator(remoteOfReferencePt, targetResidue, atomsWithThisRemote);
#ifdef _DEBUG
	if(targetResidue->getSequenceNumber() == 77)
		int here = 0;
#endif
    targetResidue->getAtoms(remoteOfReferencePt, atomsWithThisRemote);
    referenceAtom[ 0 ] = atomsWithThisRemote.getFirstEntity();
    atomsWithThisRemote.removeAll();

    if(code == ARG_AMINO_RESIDUE)
    {
        //RemoteIndicator nextRemote = getNextRemoteIndicator(remoteOfReferencePt);
        RemoteIndicator nextRemote = (RemoteIndicator)( (rg_INT)remoteOfReferencePt + 1 );
        // We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
        //getAtomsSortedByBranchDesignator(nextRemote, targetResidue, atomsWithThisRemote);
        targetResidue->getAtoms(nextRemote, atomsWithThisRemote);
        referenceAtom[ 1 ] = atomsWithThisRemote.getFirstEntity(); // get epsilon nitrogen
        atomsWithThisRemote.removeAll();

        nextRemote = (RemoteIndicator)( (rg_INT)nextRemote + 1 );		
        //getAtomsSortedByBranchDesignator(getNextRemoteIndicator(nextRemote), targetResidue, atomsWithThisRemote);
        //getAtomsSortedByBranchDesignator(nextRemote, targetResidue, atomsWithThisRemote);
        targetResidue->getAtoms(nextRemote, atomsWithThisRemote);
        referenceAtom[ 2 ] = atomsWithThisRemote.getFirstEntity(); // get zeta carbon
    }
    // Other four residues: HIS, PHE, TYR, and TRP
    else
    {
        RemoteIndicator nextRemote = (RemoteIndicator)( (rg_INT)remoteOfReferencePt + 1 );
        //getAtomsSortedByBranchDesignator(getNextRemoteIndicator(remoteOfReferencePt), targetResidue, atomsWithThisRemote);
        //getAtomsSortedByBranchDesignator(nextRemote, targetResidue, atomsWithThisRemote);
        targetResidue->getAtoms(nextRemote, atomsWithThisRemote);
        // get both first and second delta branches
        referenceAtom[ 1 ] = atomsWithThisRemote.getFirstEntity();
        referenceAtom[ 2 ] = atomsWithThisRemote.getSecondEntity();
    }

    // three reference points on updated substructure
    rg_Point3D referencePnt[ 3 ] = {referenceAtom[ 0 ]->getpAtomBall()->getCenter(),
                                    referenceAtom[ 1 ]->getpAtomBall()->getCenter(),
                                    referenceAtom[ 2 ]->getpAtomBall()->getCenter()};

    // new origin, new X-Axis, and new Z-Axis vectors from three reference points
    Vector3D vec01 = referencePnt[ 1 ] - referencePnt[ 0 ];
    Vector3D vec02 = referencePnt[ 2 ] - referencePnt[ 0 ];
    Vector3D normalVecToPlane = vec01.crossProduct(vec02);		

    rg_Point3D newOrigin   = referencePnt[ 0 ];
    Vector3D   newXAxisVec = vec01.getUnitVector();
    Vector3D   newZAxisVec = normalVecToPlane.getUnitVector();

    // construct transformation matrix for transforming the old local frame to new local frame
    rg_TMatrix3D transMat;
    transMat.convertAxis(oldOrigin, oldXAxisVec, oldZAxisVec, 
                         newOrigin, newXAxisVec, newZAxisVec);

    // glue fixed substructure into residue

    // get the points on the fixed substructure
    rg_dList<rg_Point3D> pointsOnFixedStructure;
    getCoordinatesOfFixedStructure2(code, pointsOnFixedStructure);
    // apply the transformation to the fixed substructure 
    pointsOnFixedStructure.reset4Loop();
    while(pointsOnFixedStructure.setNext4Loop())
    {
        rg_Point3D* currPoint = pointsOnFixedStructure.getpEntity();
        currPoint->setPoint(transMat * (*currPoint));
    }

    // get atom pointers corresponding to fixed substructure
    //RemoteIndicator currRemote = remoteOfReferencePt;
    rg_dList<Atom*> atomsToBeTransformed;
    // We need not call this function any more because we have already reordered the atoms of sidechain in breath first order via "InputStructureAdaptorForSCP".
    //getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsToBeTransformed, remoteOfReferencePt);

    RemoteIndicator remoteOfGluePt = (RemoteIndicator)( (rg_INT)remoteOfReferencePt + 1 );
    if(code == ARG_AMINO_RESIDUE)
    {
        remoteOfGluePt = (RemoteIndicator)( (rg_INT)remoteOfGluePt + 1 );
    }

    targetResidue->getAtomsOnSidechain(atomsToBeTransformed);
    atomsToBeTransformed.reset4Loop();
    while (atomsToBeTransformed.setNext4Loop())
    {
        Atom* currSidechainAtom = atomsToBeTransformed.getEntity();
        RemoteIndicator remote = currSidechainAtom->getpChemicalProperties()->getRemoteIndicator();
        if(remote <= remoteOfGluePt)
            atomsToBeTransformed.killCurrent();
        else
            break;
    }
    //RemoteIndicator  endRemote = UNK_REMOTE;
    //BranchDesignator endBranch = UNK_BRANCH;
    //FunctionsForRotamerLibrary::getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(targetResidue->getResidueCode(), endRemote, endBranch);
    //RemoteIndicator  startRemote = remoteOfReferencePt;
    //getAtomsOnSideChainInBreathFirstOrder(targetResidue, atomsToBeTransformed, startRemote, endRemote, endBranch);

#ifdef _DEBUG
    // test
    if(code == TYR_AMINO_RESIDUE)
        int here = 1;
#endif

    //// update atom coordinates
    atomsToBeTransformed.reset4Loop();
    pointsOnFixedStructure.reset4Loop();
    while(atomsToBeTransformed.setNext4Loop() &&
        pointsOnFixedStructure.setNext4Loop() )
    {
        Atom* currAtom = atomsToBeTransformed.getEntity();
        rg_Point3D currUpdatedCenter = pointsOnFixedStructure.getEntity();
        currAtom->getpAtomBall()->setCenter(currUpdatedCenter);
    }


#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == HIS_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE)
    {
        rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(code);
        rg_REAL* angles = new rg_REAL[numDihedralAngles];
        targetResidue->getSideChainDihedralAngles(angles);
        if(angles != rg_NULL)
            delete [] angles;
    }
#endif

#ifdef _DEBUG
    code = targetResidue->getResidueCode();
    if(code == ARG_AMINO_RESIDUE)
    {
        rg_INT numDihedralAngles = FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue(code);
        rg_REAL* angles = new rg_REAL[numDihedralAngles];
        targetResidue->getSideChainDihedralAngles(angles);
        if(angles != rg_NULL)
            delete [] angles;
    }
#endif
}




rg_BOOL FunctionsForRotamerLibrary::hasFixedSubstructure(Residue* targetResidue)
{
	ResidueCode code = targetResidue->getResidueCode();
	if( code != HIS_AMINO_RESIDUE && 
		code != PHE_AMINO_RESIDUE && 
		code != TRP_AMINO_RESIDUE && 
		code != TYR_AMINO_RESIDUE &&
		code != ARG_AMINO_RESIDUE   )
	{
		return rg_FALSE;
	}
	else
	{
		return rg_TRUE;
	}
}

void FunctionsForRotamerLibrary::getCoordinatesOfFixedStructure(const rg_INT& residueCode, rg_dList<rg_Point3D>& pointsOnFixedStructure )
{
	rg_INT index;
	getIndexOfFixedStructureInResidue(residueCode, index);
	rg_INT i;
	rg_INT num = FIXED_STRUCTURE_IN_RESIDUE[index].numOfAtoms;
	for(i = 0;i < num;i++)
	{
		pointsOnFixedStructure.add(rg_Point3D(FIXED_STRUCTURE_IN_RESIDUE[index].pts[ i ][ 0 ], 
			                                  FIXED_STRUCTURE_IN_RESIDUE[index].pts[ i ][ 1 ],
			                                  FIXED_STRUCTURE_IN_RESIDUE[index].pts[ i ][ 2 ]));
	}
}

bool FunctionsForRotamerLibrary::getIndexOfFixedStructureInResidue(const rg_INT& residueCode, rg_INT& index)
{
	bool bFixedStructureExisted = true;
	switch(residueCode)
	{
	case HIS_AMINO_RESIDUE:
		index	= 0;
		break;
	case PHE_AMINO_RESIDUE:
		index	= 1;
		break;
	case TYR_AMINO_RESIDUE:	
		index	= 2;
		break;
    case TRP_AMINO_RESIDUE:
		index	= 3;
		break;	
	case ARG_AMINO_RESIDUE:
		index	= 4;
		break;	
	default:
		index	= -1;
		bFixedStructureExisted = false;
		break;
	}
	return bFixedStructureExisted;
}

void FunctionsForRotamerLibrary::getCoordinatesOfFixedStructure2(const rg_INT& residueCode, rg_dList<rg_Point3D>& pointsOnFixedStructure )
{
    rg_INT index;
    getIndexOfFixedStructureInResidue(residueCode, index);
    rg_INT i;
    rg_INT num = FIXED_STRUCTURE_IN_RESIDUE3[index].numOfAtoms;
    for(i = 0;i < num;i++)
    {
        pointsOnFixedStructure.add(rg_Point3D(FIXED_STRUCTURE_IN_RESIDUE3[index].pts[ i ][ 0 ], 
                                              FIXED_STRUCTURE_IN_RESIDUE3[index].pts[ i ][ 1 ],
                                              FIXED_STRUCTURE_IN_RESIDUE3[index].pts[ i ][ 2 ]));
    }
}

rg_BOOL FunctionsForRotamerLibrary::getStartNEndRemoteIndicatorsOfSidechainForResidue(const ResidueCode& code, RemoteIndicator& start, RemoteIndicator& end)
{
	rg_BOOL bExist = rg_TRUE;
	start = BETA_REMOTE;
	end   = BETA_REMOTE;

	switch(code)
	{
	case ARG_AMINO_RESIDUE:		
		end = ETA_REMOTE;
		break;
	case ASN_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		break;
	case ASP_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		break;
	case CYS_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		break;	
	case GLN_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		break;	
	case GLU_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		break;	
	case HIS_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		break;	
	case ILE_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		break;	
	case LEU_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		break;	
	case LYS_AMINO_RESIDUE:
		end = ZETA_REMOTE;
		break;	
	case MET_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		break;
	case PHE_AMINO_RESIDUE:
		end = ZETA_REMOTE;
		break;	
	case PRO_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		break;	
	case SER_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		break;	
	case THR_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		break;	
	case TRP_AMINO_RESIDUE:
		end = ETA_REMOTE;
		break;	
	case TYR_AMINO_RESIDUE:
		end = ETA_REMOTE;
		break;	
	case VAL_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		break;
	default:
		bExist = rg_FALSE; 
		break;
	}

	return bExist;
}

rg_BOOL FunctionsForRotamerLibrary::getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(const ResidueCode& code, RemoteIndicator& end, BranchDesignator& branchDesignator)
{
	rg_BOOL bExist = rg_TRUE;
	end = BETA_REMOTE;

	switch(code)
	{
	case ARG_AMINO_RESIDUE:		
		end = ETA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;
	case ASN_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;
	case ASP_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;
	case CYS_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;	
	case GLN_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;	
	case GLU_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;	
	case HIS_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;	
	case ILE_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		branchDesignator = FIRST_BRANCH;
		break;	
	case LEU_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;	
	case LYS_AMINO_RESIDUE:
		end = ZETA_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;	
	case MET_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;
	case PHE_AMINO_RESIDUE:
		end = ZETA_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;	
	case PRO_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;	
	case SER_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;	
	case THR_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;	
	case TRP_AMINO_RESIDUE:
		end = ETA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;	
	case TYR_AMINO_RESIDUE:
		end = ETA_REMOTE;
		branchDesignator = UNK_BRANCH;
		break;	
	case VAL_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		branchDesignator = SECOND_BRANCH;
		break;
	default:
		bExist = rg_FALSE; 
		break;
	}

	return bExist;

	/*{m_remoteIndicator=ETA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=DELTA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=DELTA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=GAMMA_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=EPSILON_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=EPSILON_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=EPSILON_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=DELTA_REMOTE m_brangeDesignator=FIRST_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=DELTA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=ZETA_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=EPSILON_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=ZETA_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=DELTA_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=GAMMA_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=GAMMA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=ETA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=ETA_REMOTE m_brangeDesignator=UNK_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	
	{m_remoteIndicator=GAMMA_REMOTE m_brangeDesignator=SECOND_BRANCH m_extraBrangeDesignator=UNK_BRANCH ...}	*/	
}

rg_BOOL FunctionsForRotamerLibrary::getEndRemoteIndicatorNMaximumBranchDesignatorOfSidechainForResidue(const ResidueCode& code, RemoteIndicator& end, BranchDesignator& maxBranchDesignator)
{
	rg_BOOL bExist = rg_TRUE;
	end = BETA_REMOTE;

	switch(code)
	{
	case ARG_AMINO_RESIDUE:		
		end = ETA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;
	case ASN_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;
	case ASP_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;
	case CYS_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		maxBranchDesignator = UNK_BRANCH;
		break;	
	case GLN_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case GLU_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case HIS_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case ILE_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case LEU_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case LYS_AMINO_RESIDUE:
		end = ZETA_REMOTE;
		maxBranchDesignator = UNK_BRANCH;
		break;	
	case MET_AMINO_RESIDUE:
		end = EPSILON_REMOTE;
		maxBranchDesignator = UNK_BRANCH;
		break;
	case PHE_AMINO_RESIDUE:
		end = ZETA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case PRO_AMINO_RESIDUE:
		end = DELTA_REMOTE;
		maxBranchDesignator = UNK_BRANCH;
		break;	
	case SER_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		maxBranchDesignator = UNK_BRANCH;
		break;	
	case THR_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case TRP_AMINO_RESIDUE:
		end = ETA_REMOTE;
		maxBranchDesignator = THIRD_BRANCH;
		break;	
	case TYR_AMINO_RESIDUE:
		end = ETA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;	
	case VAL_AMINO_RESIDUE:
		end = GAMMA_REMOTE;
		maxBranchDesignator = SECOND_BRANCH;
		break;
	default:
		bExist = rg_FALSE; 
		break;
	}

	return bExist;
}

rg_INT FunctionsForRotamerLibrary::getAtomsOfResidue(Residue* residue, rg_dList<Atom*>& atomsOfResidue)
{
	RemoteIndicator  endRemote = UNK_REMOTE;
	BranchDesignator maxBranch = UNK_BRANCH;
	getEndRemoteIndicatorNMaximumBranchDesignatorOfSidechainForResidue(residue->getResidueCode(), endRemote, maxBranch);

	getAtomsOfResidueWithGivenRemoteIndicatorNBranchDesignator(residue, atomsOfResidue, endRemote, maxBranch);

	return atomsOfResidue.getSize();
}

rg_INT FunctionsForRotamerLibrary::getAtomsOfResidue(Residue* residue, Atom**& atomsOfResidue)
{
	rg_dList<Atom*> residueAtomsInList;
	rg_INT numAtoms = getAtomsOfResidue(residue, residueAtomsInList);

	atomsOfResidue = new Atom* [numAtoms];

	rg_INDEX index = 0;
	residueAtomsInList.reset4Loop();
	while (residueAtomsInList.setNext4Loop())
	{
		atomsOfResidue[ index ] = residueAtomsInList.getEntity();
		index++;
	}

	return numAtoms;
}

rg_INT FunctionsForRotamerLibrary::getAtomsOnSidechain(Residue* residue, rg_dList<Atom*>& atomsOnSidechain )
{
	RemoteIndicator  endRemote = UNK_REMOTE;
	BranchDesignator maxBranch = UNK_BRANCH;
	getEndRemoteIndicatorNMaximumBranchDesignatorOfSidechainForResidue(residue->getResidueCode(), endRemote, maxBranch);
	RemoteIndicator  startRemote = BETA_REMOTE;

	getAtomsOnSidechainWithGivenRemoteIndicatorNBranchDesignator(residue, atomsOnSidechain, startRemote, endRemote, maxBranch);

	return atomsOnSidechain.getSize();
}

rg_INT FunctionsForRotamerLibrary::getAtomsOnSidechain(Residue* residue, Atom**& atomsOnSidechain )
{
	// get the necessary atoms by considering the remote indicator and branch designator

	rg_dList<Atom*> sidechainAtomsInList;
	rg_INT numAtoms = getAtomsOnSidechain(residue, sidechainAtomsInList);

    if(numAtoms > 0)
	    atomsOnSidechain = new Atom* [numAtoms];

	rg_INDEX index = 0;
	sidechainAtomsInList.reset4Loop();
	while (sidechainAtomsInList.setNext4Loop())
	{
		atomsOnSidechain[ index ] = sidechainAtomsInList.getEntity();
		index++;
	}

	return numAtoms;
}

/*
void FunctionsForRotamerLibrary::writeBBDepRotLibIntoBinaryFile(const string& fileNameWithPath)
{
	ofstream fout( fileNameWithPath.c_str(), ios::out | ios::binary );

	if(! fout)
	{
		cout << "cannot open file. \n";
	}

	// write BBDep Rotamer Library into file
	rg_INDEX i;
	for (i = 0;i < SIZE_OF_BBDEP_ROTAMER_LIB;i++)
	{
		
		// read from BBD_ROTAMER_LIB0
		BBDEP_ROTAMER_LIB[ i ].code = BBD_ROTAMER_LIB0[ i ].code;
		BBDEP_ROTAMER_LIB[ i ].phi  = BBD_ROTAMER_LIB0[ i ].phi;
		BBDEP_ROTAMER_LIB[ i ].psi  = BBD_ROTAMER_LIB0[ i ].psi;
		BBDEP_ROTAMER_LIB[ i ].numOfObservations = BBD_ROTAMER_LIB0[ i ].numOfObservations;
		rg_INDEX j;
		for (j = 0;j < 4;j++)
		{
			BBDEP_ROTAMER_LIB[ i ].chiRotamer[ j ]   = BBD_ROTAMER_LIB0[ i ].chiRotamer[ j ];
			BBDEP_ROTAMER_LIB[ i ].chi[ j ]          = BBD_ROTAMER_LIB0[ i ].chi[ j ];
			BBDEP_ROTAMER_LIB[ i ].stdDeviation[ j ] = BBD_ROTAMER_LIB0[ i ].stdDeviation[ j ];			
		}
		BBDEP_ROTAMER_LIB[ i ].probability = BBD_ROTAMER_LIB0[ i ].probability;
		
		fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ], sizeof(BBDepRotamer) );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].residueThreeCodeName, sizeof(char) * 4);
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].phi, sizeof(char) );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].psi, sizeof(char) );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].numOfObservations, sizeof(char) );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].chiRotamer, sizeof(char) * 4 );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].probability, sizeof(rg_REAL) );
		////fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].numOfChis, sizeof(char) );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].chi, sizeof(rg_REAL) * 4 );
		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].stdDeviation, sizeof(rg_REAL) * 4);
		
		fout.write((const char*) &BBD_ROTAMER_LIB0[ i ], sizeof(BBDepRotamer) );
	}

	fout.close();

	if(! fout.good() )
	{
		cout << "A file error occurred.";
	}

	//rg_INT maxNumOfObservations = - 10000000;
	//rg_INT minNumOfObservations = + 10000000;

	//rg_INT maxChiRotamer = - 10000000;
	//rg_INT minChiRotamer = + 10000000;

	//rg_INDEX i;
	//for (i = 0;i < SIZE_OF_BBDEP_ROTAMER_LIB;i++)
	//{
	//	rg_INT numObservations = BBDEP_ROTAMER_LIB[ i ].numOfObservations;
	//	
	//	if(numObservations > maxNumOfObservations)
	//		maxNumOfObservations = numObservations;

	//	if(numObservations < minNumOfObservations)
	//		minNumOfObservations = numObservations;

	//	rg_INDEX j;
	//	for (j = 0;j < 4;j++)
	//	{
	//		rg_INT chiRotamer = BBDEP_ROTAMER_LIB[ i ].chiRotamer[ j ];

	//		if(chiRotamer > maxChiRotamer)
	//			maxChiRotamer = chiRotamer;

	//		if(chiRotamer < minChiRotamer)
	//			minChiRotamer = chiRotamer;
	//	}		
	//}

	//cout << "MAX. Number of Observations: " << maxNumOfObservations << endl;
	//cout << "MIN. Number of Observations: " << minNumOfObservations << endl;

	//cout << "MAX. Chi Rotamer: " << maxChiRotamer << endl;
	//cout << "MIN. Chi Rotamer: " << minChiRotamer << endl;
}
*/



void FunctionsForRotamerLibrary::writeBBIndepRotLibIntoBinaryFile(const string& fileNameWithPath)
{
	ofstream fout( fileNameWithPath.c_str(), ios::out | ios::binary );

	if(! fout)
	{
		cout << "cannot open file. \n";
	}

	// write BBDep Rotamer Library into file
	rg_INDEX i;
	for (i = 0;i < SIZE_OF_BBINDEP_ROTAMER_LIB;i++)
	{
		// read from BBINDEP_ROTAMER_LIB0
//		BBINDEP_ROTAMER_LIB[ i ].code = BBINDEP_ROTAMER_LIB0[ i ].code;
//		BBINDEP_ROTAMER_LIB[ i ].numOfObservationsForChiOneRotamers  = BBINDEP_ROTAMER_LIB0[ i ].numOfObservationsForChiOneRotamers;
//		BBINDEP_ROTAMER_LIB[ i ].numOfObservations  = BBINDEP_ROTAMER_LIB0[ i ].numOfObservations;
//		BBINDEP_ROTAMER_LIB[ i ].probability = BBINDEP_ROTAMER_LIB0[ i ].probability;		
//		BBINDEP_ROTAMER_LIB[ i ].stdDeviation = BBINDEP_ROTAMER_LIB0[ i ].stdDeviation;
//		BBINDEP_ROTAMER_LIB[ i ].condProbabilityGivenChiOneRotamer = BBINDEP_ROTAMER_LIB0[ i ].condProbabilityGivenChiOneRotamer;
//		BBINDEP_ROTAMER_LIB[ i ].stdDeviationGivenChiOneRotamer = BBINDEP_ROTAMER_LIB0[ i ].stdDeviationGivenChiOneRotamer;
//		
//		rg_INDEX j;
//		for (j = 0;j < 4;j++)
//		{
//			BBINDEP_ROTAMER_LIB[ i ].chiRotamer[ j ]   = BBINDEP_ROTAMER_LIB0[ i ].chiRotamer[ j ];
//			BBINDEP_ROTAMER_LIB[ i ].chi[ j ]          = BBINDEP_ROTAMER_LIB0[ i ].chi[ j ];
//			BBINDEP_ROTAMER_LIB[ i ].stdDeviationOfChi[ j ] = BBINDEP_ROTAMER_LIB0[ i ].stdDeviationOfChi[ j ];			
//		}
		
		fout.write((const char*) &BBINDEP_ROTAMER_LIB0[ i ], sizeof(BBIndRotamer0) );
	}

	fout.close();

	if(! fout.good() )
	{
		cout << "A file error occurred.";
	}
}


void FunctionsForRotamerLibrary::loadBBDepRotLib(const string& fileNameWithPath)
{
    m_currentlyUsedRotamerLibType = BACKBONE_DEPENDENT_LIB_DUNBRACK2002;

    if(m_bBBDEPLibLoaded)
        return;
    else
        m_bBBDEPLibLoaded = rg_TRUE;

	ifstream fin( fileNameWithPath.c_str(), ios::in | ios::binary );

	if(! fin)
	{
		cout << "cannot open rotamer library file. \n";
		// May 26, 2016 by Joonghyun
		exit( 1 );
	}

	// read BBDep Rotamer Library from file

    //BBDEP_ROTAMER_LIB = new BBDepRotamer[SIZE_OF_BBDEP_ROTAMER_LIB];

	rg_INDEX i;
	for (i = 0;i < SIZE_OF_BBDEP_ROTAMER_LIB;i++)
	{
		fin.read((char*) &BBDEP_ROTAMER_LIB[ i ], sizeof(BBDepRotamer) );
		//char residueThreeCodeName[ 4 ];
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].residueThreeCodeName, sizeof(char) * 4);
		//char phi, psi;
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].phi, sizeof(char) );
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].psi, sizeof(char) );
		//char numOfObservations;
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].numOfObservations, sizeof(char) );
		//char chiRotamer[ 4 ];
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].chiRotamer, sizeof(char) * 4 );
		//rg_REAL probability;
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].probability, sizeof(rg_REAL) );
		////fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].numOfChis, sizeof(char) );
		//rg_REAL chi[ 4 ];
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].chi, sizeof(rg_REAL) * 4 );
		//rg_REAL stdDeviation[ 4 ];
		//fin.read((const char*) &BBDEP_ROTAMER_LIB[ i ].stdDeviation, sizeof(rg_REAL) * 4);
	}

	fin.close();

	if(! fin.good() )
	{
		// cout << "A file error occurred.";
		// May 26, 2016 by Joonghyun
		cout << "An error occurred during loading the backbone dependnet rotamer library.";
		exit( 1 );
	}
}

void FunctionsForRotamerLibrary::loadBBIndepRotLib(const string& fileNameWithPath)
{
    m_currentlyUsedRotamerLibType = BACKBONE_INDEPENDENT_LIB_DUNBRACK2002;

    if(m_bBBINDLibLoaded)
        return;
    else
        m_bBBINDLibLoaded = rg_TRUE;

	ifstream fin( fileNameWithPath.c_str(), ios::in | ios::binary );

	if(! fin)
	{
		cout << "cannot open rotamer library file. \n";
		// May 26, 2016 by Joonghyun
		exit(1);
	}

	// read BBIndep Rotamer Library from file
	rg_INDEX i;
	for (i = 0;i < SIZE_OF_BBINDEP_ROTAMER_LIB;i++)
	{
		fin.read((char*) &BBINDEP_ROTAMER_LIB[ i ], sizeof(BBIndRotamer) );
	}

	fin.close();

	if(! fin.good() )
	{
		//cout << "A file error occurred.";
		// May 26, 2016 by Joonghyun
		cout << "An error occurred during loading the backbone independent rotamer library.";
		exit(1);
	}
}

void FunctionsForRotamerLibrary::loadBBDepRotLib_Dunbrack10(const string& fileNameWithPath)
{
    m_currentlyUsedRotamerLibType = BACKBONE_DEPENDENT_LIB_DUNBRACK2010;

    if(m_bBBDEPLibLoaded)
        return;
    else
        m_bBBDEPLibLoaded = rg_TRUE;


}

//void FunctionsForRotamerLibrary::loadBBDepRotLib(const string& fileNameWithPath, BBDepRotamer rotLib[SIZE_OF_BBDEP_ROTAMER_LIB])
//{
//	ifstream fin( fileNameWithPath.c_str(), ios::in | ios::binary );
//
//	if(! fin)
//	{
//		cout << "cannot open rotamer library file. \n";
//	}
//
//	// read BBDep Rotamer Library into file
//	rg_INDEX i;
//	for (i = 0;i < SIZE_OF_BBDEP_ROTAMER_LIB;i++)
//	{
//		fin.read((char*) &rotLib[ i ], sizeof(BBDepRotamer) );
//		cout << rotLib[ i ].code << " " 
//			 << rotLib[ i ].phi << " " << rotLib[ i ].psi << " " 
//			 << rotLib[ i ].numOfObservations << " " 
//			 << rotLib[ i ].chiRotamer[ 0 ] << rotLib[ i ].chiRotamer[ 1 ] << rotLib[ i ].chiRotamer[ 2 ] << rotLib[ i ].chiRotamer[ 3 ] << " "
//			 << rotLib[ i ].probability << " "
//			 << rotLib[ i ].chi[ 0 ] << " " 
//			 << rotLib[ i ].chi[ 1 ] << " "
//			 << rotLib[ i ].chi[ 2 ] << " "
//			 << rotLib[ i ].chi[ 3 ] << " "
//			 << rotLib[ i ].stdDeviation[ 0 ] << " "
//			 << rotLib[ i ].stdDeviation[ 1 ] << " "
//			 << rotLib[ i ].stdDeviation[ 2 ] << " "
//			 << rotLib[ i ].stdDeviation[ 3 ] << endl;
//	}
//
//	fin.close();
//
//	if(! fin.good() )
//	{
//		cout << "A file error occurred.";
//	}
//}





//////////////////////////////////////////////////////////////////////////
// For Dunbrack 2010 - by JHCha
//////////////////////////////////////////////////////////////////////////

rg_BOOL FunctionsForRotamerLibrary::getBBDepLibIndexOfResidue_2010(const ResidueCode& code,
	rg_INDEX& start,
	rg_INDEX& end)
{
	rg_INDEX index = getIndexOfResCode(code);
	if ((0 <= index) && (index <= 17))
	{
		start = INDEX_PAIR_BBDEP_ROTAMERS[index].m_startIndex;
		end = INDEX_PAIR_BBDEP_ROTAMERS[index].m_endIndex;
		return rg_TRUE;
	}
	else
	{
		start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		return rg_FALSE;
	}
}

rg_BOOL FunctionsForRotamerLibrary::getBBDepLibIndexOfResidue_2010(const ResidueCode& code,
	const rg_REAL& phi,
	const rg_REAL& psi,
	rg_INDEX& start,
	rg_INDEX& end)
{
	rg_INDEX residueIndex = getIndexOfResCode_2010(code);
	if ((0 <= residueIndex) && (residueIndex <= 21))
	{
		// find each index of bin which contains the backbone dihedral angle (phi, psi) using pre-defined bucket
		rg_INDEX phiIndex = -1;
		rg_INDEX psiIndex = -1;
		if (getAngleBinIndex(phi, phiIndex) && getAngleBinIndex(psi, psiIndex))
		{
			start = INDEX_PAIR_BUCKET_FOR_BBDEP_ROTAMERS_2010[residueIndex][phiIndex][psiIndex].m_startIndex;
			end = INDEX_PAIR_BUCKET_FOR_BBDEP_ROTAMERS_2010[residueIndex][phiIndex][psiIndex].m_endIndex;

			//cout << "Curr code: " << residueIndex << ", Phi: " << phi << ", Psi: " << psi << ", start: " << start << ", end: " << end << endl;
			return rg_TRUE;
		}
		else
		{
			start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
			return rg_FALSE;
		}
	}
	else
	{
		start = end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		return rg_FALSE;
	}
}

rg_INDEX FunctionsForRotamerLibrary::getIndexOfResCode_2010(const ResidueCode& code)
{
	rg_INDEX indexOfResidue = -1;
	switch (code)
	{
	case ARG_AMINO_RESIDUE:
		indexOfResidue = 0;
		break;
	case ASN_AMINO_RESIDUE:
		indexOfResidue = 1;
		break;
	case ASP_AMINO_RESIDUE:
		indexOfResidue = 2;
		break;
/*
	case CPR_AMINO_RESIDUE:
		indexOfResidue = 3;
		break;
	case CYD_AMINO_RESIDUE:
		indexOfResidue = 4;
		break;
	case CYH_AMINO_RESIDUE:
		indexOfResidue = 5;
		break;
	*/
	case CYS_AMINO_RESIDUE:
		indexOfResidue = 6;
		break;
	case GLN_AMINO_RESIDUE:
		indexOfResidue = 7;
		break;
	case GLU_AMINO_RESIDUE:
		indexOfResidue = 8;
		break;
	case HIS_AMINO_RESIDUE:
		indexOfResidue = 9;
		break;
	case ILE_AMINO_RESIDUE:
		indexOfResidue = 10;
		break;
	case LEU_AMINO_RESIDUE:
		indexOfResidue = 11;
		break;
	case LYS_AMINO_RESIDUE:
		indexOfResidue = 12;
		break;
	case MET_AMINO_RESIDUE:
		indexOfResidue = 13;
		break;
	case PHE_AMINO_RESIDUE:
		indexOfResidue = 14;
		break;
	case PRO_AMINO_RESIDUE:
		indexOfResidue = 15;
		break;
	case SER_AMINO_RESIDUE:
		indexOfResidue = 16;
		break;
	case THR_AMINO_RESIDUE:
		indexOfResidue = 17;
		break;
	/*case TPR_AMINO_RESIDUE:
		indexOfResidue = 18;
		break;
	*/
	case TRP_AMINO_RESIDUE:
		indexOfResidue = 19;
		break;
	case TYR_AMINO_RESIDUE:
		indexOfResidue = 20;
		break;
	case VAL_AMINO_RESIDUE:
		indexOfResidue = 21;
		break;
		//case ALA_AMINO_RESIDUE:
		//	indexOfResidue = 18;
		//	break;
		//case GLY_AMINO_RESIDUE:
		//	indexOfResidue = 19;
		//	break;
	default:
		break;
	}

	return indexOfResidue;
}


rg_INT FunctionsForRotamerLibrary::getNumDihedralAnglesOfResidue_2010(const ResidueCode& code)
{
#ifdef _DEBUG
	if (code == PRO_AMINO_RESIDUE)
		int here = 1;
#endif

	if (code == ALA_AMINO_RESIDUE || code == GLY_AMINO_RESIDUE)
	{
		return 0;
	}
	else
	{
		// April 5, 2016 by Joonghyun
		//return NUM_DIHEDRAL_ANGLES_2010[getIndexOfResCode(code)];
		return NUM_DIHEDRAL_ANGLES_2010[getIndexOfResCode_2010(code)];
	}

	//if(code == SER_AMINO_RESIDUE || code == THR_AMINO_RESIDUE ||
	//   code == VAL_AMINO_RESIDUE || code == CYS_AMINO_RESIDUE  )
	//{
	//	return 1;
	//}
	//else if(code == ASP_AMINO_RESIDUE || code == LEU_AMINO_RESIDUE || code == ILE_AMINO_RESIDUE ||
	//	    code == ASN_AMINO_RESIDUE || code == PHE_AMINO_RESIDUE || code == TYR_AMINO_RESIDUE ||
	//		code == PRO_AMINO_RESIDUE || code == HIS_AMINO_RESIDUE || code == TRP_AMINO_RESIDUE   )
	//{
	//		return 2;
	//}
	//else if(code == LYS_AMINO_RESIDUE || code == ARG_AMINO_RESIDUE)
	//{
	//	return 4;
	//}
	//else if(code == GLU_AMINO_RESIDUE ||
	//	    code == GLN_AMINO_RESIDUE ||
	//	    code == MET_AMINO_RESIDUE  )
	//{
	//	return 3;
	//}
	////if(code == ALA_AMINO_RESIDUE || code == GLY_AMINO_RESIDUE)
	//else
	//{
	//	return 0;
	//}
}



void FunctionsForRotamerLibrary::getRotamerSetsOfResidueSetFromBBDepRotLib_2010(Backbone & backbone, Rotamer* sortedAssignedRotamersOfProtein, RotamerSetOfResidue* & sortedRotamerSetsOfResidues, ManagerOfRotamerSetsForSCP & managerOfRotamerSetsForSCP)
{
	// get backbone conformation
	BackboneConformation backboneConformation;
	backbone.getBackboneConformation(backboneConformation);

	rg_dList<BackboneDihedralAnglePair>* backboneDihedralAnglePairs
		= backboneConformation.getBackboneDihedralAnglePairs();
	rg_INT numResiduePosition = backboneDihedralAnglePairs->getSize();

	// get sequence
	rg_dList<ResidueCode> sequence;
	backbone.getSequence(sequence);

	getRotamerSetsForResiduesSortedByChainIDSeqNumber(sortedAssignedRotamersOfProtein,
		numResiduePosition,
		sortedRotamerSetsOfResidues);
	// reserve size for rotamer sets
	managerOfRotamerSetsForSCP.setNumResidues(numResiduePosition);

	// get rotamer sets considering the backbone dihedral angles
	rg_INDEX i = 0;
	backboneDihedralAnglePairs->reset4Loop();
	sequence.reset4Loop();
	while (backboneDihedralAnglePairs->setNext4Loop() && sequence.setNext4Loop())
	{
		rg_REAL dihedralAnglePair[2] = { DEFAULT_PHI_N_TERMINAL, DEFAULT_PSI_C_TERMINAL };
		BackboneDihedralAnglePair* currAnglePair = backboneDihedralAnglePairs->getpEntity();

#ifdef _DEBUG
		// test
		if (currAnglePair->getChainIDSeqNum() == 100187)
			int here = 1;
#endif

		currAnglePair->getDihedralAnglePair(dihedralAnglePair);
		ResidueCode currCode = sequence.getEntity();

		rg_INDEX start = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
		rg_INDEX end = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;

		getBBDepLibIndexOfResidue_2010(currCode,
			dihedralAnglePair[0],
			dihedralAnglePair[1],
			start,
			end);

		sortedRotamerSetsOfResidues[i].setRotLibIDs(start, end);
		managerOfRotamerSetsForSCP.setRotamerSet(&sortedRotamerSetsOfResidues[i], i);
		i++;
	}
}



void FunctionsForRotamerLibrary::loadBBDepRotLib_Dunbrack2010_ASCII(const string& fileNameWithPath)
{
	m_currentlyUsedRotamerLibType = BACKBONE_DEPENDENT_LIB_DUNBRACK2010;
	const int rotamerLibSize = 741216;

	if (m_bBBDEPLibLoaded)
		return;
	else
		m_bBBDEPLibLoaded = rg_TRUE;

	char	    lineOfData[1000];
	const char* seps = " \t\n";
	char*	    token = NULL;
	char*       context = NULL;

	int resID = 0;

	ifstream    fin(fileNameWithPath.c_str(), ios::in);
	for (int i = 0; i < rotamerLibSize; i++)
	{
		fin.getline(lineOfData, 1000);
		//token = strtok_s(lineOfData, seps, &context); //# or ResCode
		token = strtok(lineOfData, seps); //# or ResCode
		if (token[0] != '#')
		{
			BBDEP_ROTAMER_LIB_2010[resID].code = translate_three_character_residue_code(token);

			//token = strtok_s(NULL, seps, &context); // phi
			token = strtok(NULL, seps); // phi
			BBDEP_ROTAMER_LIB_2010[resID].phi = atoi(token);

			//token = strtok_s(NULL, seps, &context); // psi
			token = strtok(NULL, seps); // psi
			BBDEP_ROTAMER_LIB_2010[resID].psi = atoi(token);

			//token = strtok_s(NULL, seps, &context); // count
			token = strtok(NULL, seps); // count
			BBDEP_ROTAMER_LIB_2010[resID].numOfObservations = atoi(token);

			for (int chiIndex = 0; chiIndex < 4; chiIndex++)
			{
				//token = strtok_s(NULL, seps, &context); // chi rotamer
				token = strtok(NULL, seps); // chi rotamer
				BBDEP_ROTAMER_LIB_2010[resID].chiRotamer[chiIndex] = atoi(token);
			}

			//token = strtok_s(NULL, seps, &context); // probability
			token = strtok(NULL, seps); // probability
			BBDEP_ROTAMER_LIB_2010[resID].probability = atof(token);

			for (int chiIndex = 0; chiIndex < 4; chiIndex++)
			{
				//token = strtok_s(NULL, seps, &context); // chi value
				token = strtok(NULL, seps); // chi value
				BBDEP_ROTAMER_LIB_2010[resID].chi[chiIndex] = atof(token);
			}

			for (int chiIndex = 0; chiIndex < 4; chiIndex++)
			{
				//token = strtok_s(NULL, seps, &context); // stdDeviation
				token = strtok(NULL, seps); // stdDeviation
				BBDEP_ROTAMER_LIB_2010[resID].stdDeviation[chiIndex] = atof(token);
			}
			resID++;
		}
	}
	fin.close();
}

ResidueCode FunctionsForRotamerLibrary::translate_three_character_residue_code(char* threeCharResCode)
{
	string s_resCode(threeCharResCode);

	if (s_resCode == "ARG")
	{
		return ARG_AMINO_RESIDUE;
	}
	if (s_resCode == "ASN")
	{
		return ASN_AMINO_RESIDUE;
	}
	if (s_resCode == "ASP")
	{
		return ASP_AMINO_RESIDUE;
	}
	if (s_resCode == "CPR")
	{
		//return CPR_AMINO_RESIDUE;
		return PRO_AMINO_RESIDUE;
	}
	if (s_resCode == "CYD")
	{
		//return CYD_AMINO_RESIDUE;
		return CYS_AMINO_RESIDUE;
	}
	if (s_resCode == "CYH")
	{
		//return CYH_AMINO_RESIDUE;
		return CYS_AMINO_RESIDUE;
	}
	if (s_resCode == "CYS")
	{
		return CYS_AMINO_RESIDUE;
	}
	if (s_resCode == "GLN")
	{
		return GLN_AMINO_RESIDUE;
	}
	if (s_resCode == "GLU")
	{
		return GLU_AMINO_RESIDUE;
	}
	if (s_resCode == "HIS")
	{
		return HIS_AMINO_RESIDUE;
	}
	if (s_resCode == "ILE")
	{
		return ILE_AMINO_RESIDUE;
	}
	if (s_resCode == "LEU")
	{
		return LEU_AMINO_RESIDUE;
	}
	if (s_resCode == "LYS")
	{
		return LYS_AMINO_RESIDUE;
	}
	if (s_resCode == "MET")
	{
		return MET_AMINO_RESIDUE;
	}
	if (s_resCode == "PHE")
	{
		return PHE_AMINO_RESIDUE;
	}
	if (s_resCode == "PRO")
	{
		return PRO_AMINO_RESIDUE;
	}
	if (s_resCode == "SER")
	{
		return SER_AMINO_RESIDUE;
	}
	if (s_resCode == "THR")
	{
		return THR_AMINO_RESIDUE;
	}
	if (s_resCode == "TPR")
	{
		//return TPR_AMINO_RESIDUE;
		return PRO_AMINO_RESIDUE;
	}
	if (s_resCode == "TRP")
	{
		return TRP_AMINO_RESIDUE;
	}
	if (s_resCode == "TYR")
	{
		return TYR_AMINO_RESIDUE;
	}
	if (s_resCode == "VAL")
	{
		return VAL_AMINO_RESIDUE;
	}
}

void FunctionsForRotamerLibrary::loadBBDepRotLib_Dunbrack2010_BIN(const string& fileNameWithPath)
{
	m_currentlyUsedRotamerLibType = BACKBONE_DEPENDENT_LIB_DUNBRACK2010;

	if (m_bBBDEPLibLoaded)
		return;
	else
		m_bBBDEPLibLoaded = rg_TRUE;

	ifstream fin(fileNameWithPath.c_str(), ios::in | ios::binary);

	if (!fin)
	{
		cout << "cannot open rotamer library file. \n";
		// May 26, 2016 by Joonghyun
		exit( 1 );
	}

	// read BBDep Rotamer Library from file

	//BBDEP_ROTAMER_LIB = new BBDepRotamer[SIZE_OF_BBDEP_ROTAMER_LIB];

	rg_INDEX i;
	for (i = 0; i < SIZE_OF_BBDEP_ROTAMER_LIB_2010; i++)
	{
		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i], sizeof(BBDepRotamer));
		// 		char residueThreeCodeName[ 4 ];
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].code, sizeof(char) * 4);
		// 		char phi, psi;
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].phi, sizeof(char));
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].psi, sizeof(char));
		// 		char numOfObservations;
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].numOfObservations, sizeof(char));
		// 		char chiRotamer[ 4 ];
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].chiRotamer, sizeof(char) * 4);
		// 		rg_REAL probability;
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].probability, sizeof(rg_REAL));
		// 		//fin.read((char*) &BBDEP_ROTAMER_LIB[ i ].numOfChis, sizeof(char) );
		// 		rg_REAL chi[ 4 ];
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].chi, sizeof(rg_REAL) * 4);
		// 		rg_REAL stdDeviation[ 4 ];
		// 		fin.read((char*)&BBDEP_ROTAMER_LIB_2010[i].stdDeviation, sizeof(rg_REAL) * 4);
	}

	fin.close();

	if (!fin.good())
	{
		//cout << "A file error occurred.";
		// May 26, 2016 by Joonghyun
		cout << "An error occurred during loading the backbone dependent rotamer library 2010.";
		exit( 1 );
	}
}



void FunctionsForRotamerLibrary::writeBBDepRotLibIntoBinaryFile_2010(const string& fileNameWithPath)
{
	ofstream fout(fileNameWithPath.c_str(), ios::out | ios::binary);

	if (!fout)
	{
		cout << "cannot open file. \n";
	}

	// write BBDep Rotamer Library into file
	rg_INDEX i;
	for (i = 0; i < SIZE_OF_BBDEP_ROTAMER_LIB_2010; i++)
	{
		/*
		// read from BBD_ROTAMER_LIB0
		BBDEP_ROTAMER_LIB[i].code = BBD_ROTAMER_LIB0[i].code;
		BBDEP_ROTAMER_LIB[i].phi = BBD_ROTAMER_LIB0[i].phi;
		BBDEP_ROTAMER_LIB[i].psi = BBD_ROTAMER_LIB0[i].psi;
		BBDEP_ROTAMER_LIB[i].numOfObservations = BBD_ROTAMER_LIB0[i].numOfObservations;
		rg_INDEX j;
		for (j = 0; j < 4; j++)
		{
		BBDEP_ROTAMER_LIB[i].chiRotamer[j] = BBD_ROTAMER_LIB0[i].chiRotamer[j];
		BBDEP_ROTAMER_LIB[i].chi[j] = BBD_ROTAMER_LIB0[i].chi[j];
		BBDEP_ROTAMER_LIB[i].stdDeviation[j] = BBD_ROTAMER_LIB0[i].stdDeviation[j];
		}
		BBDEP_ROTAMER_LIB[i].probability = BBD_ROTAMER_LIB0[i].probability;
		*/

		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i], sizeof(BBDepRotamer));
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].code, sizeof(char) * 4);
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].phi, sizeof(char));
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].psi, sizeof(char));
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].numOfObservations, sizeof(char));
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].chiRotamer, sizeof(char) * 4);
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].probability, sizeof(rg_REAL));
		// 		//fout.write((const char*) &BBDEP_ROTAMER_LIB[ i ].numOfChis, sizeof(char) );
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].chi, sizeof(rg_REAL) * 4);
		// 		fout.write((const char*)&BBDEP_ROTAMER_LIB_2010[i].stdDeviation, sizeof(rg_REAL) * 4);
	}

	fout.close();

	if (!fout.good())
	{
		cout << "A file error occurred.";
	}

	//rg_INT maxNumOfObservations = - 10000000;
	//rg_INT minNumOfObservations = + 10000000;

	//rg_INT maxChiRotamer = - 10000000;
	//rg_INT minChiRotamer = + 10000000;

	//rg_INDEX i;
	//for (i = 0;i < SIZE_OF_BBDEP_ROTAMER_LIB;i++)
	//{
	//	rg_INT numObservations = BBDEP_ROTAMER_LIB[ i ].numOfObservations;
	//
	//	if(numObservations > maxNumOfObservations)
	//		maxNumOfObservations = numObservations;

	//	if(numObservations < minNumOfObservations)
	//		minNumOfObservations = numObservations;

	//	rg_INDEX j;
	//	for (j = 0;j < 4;j++)
	//	{
	//		rg_INT chiRotamer = BBDEP_ROTAMER_LIB[ i ].chiRotamer[ j ];

	//		if(chiRotamer > maxChiRotamer)
	//			maxChiRotamer = chiRotamer;

	//		if(chiRotamer < minChiRotamer)
	//			minChiRotamer = chiRotamer;
	//	}
	//}

	//cout << "MAX. Number of Observations: " << maxNumOfObservations << endl;
	//cout << "MIN. Number of Observations: " << minNumOfObservations << endl;

	//cout << "MAX. Chi Rotamer: " << maxChiRotamer << endl;
	//cout << "MIN. Chi Rotamer: " << minChiRotamer << endl;
}



bool FunctionsForRotamerLibrary::compare_BBDepRotLib_Dunbrack2010_BIN_to_ASCII(const string& fileNameWithPath)
{
	bool isTwoRotLibSame = true;

	const int rotamerLibSize = 741216;

	char	    lineOfData[1000];
	const char* seps = " \t\n";
	char*	    token = NULL;
	char*       context = NULL;

	int resID = 0;

	ifstream    fin(fileNameWithPath.c_str(), ios::in);
	for (int i = 0; i < rotamerLibSize; i++)
	{
		fin.getline(lineOfData, 1000);
		//token = strtok_s(lineOfData, seps, &context); //# or ResCode
		token = strtok(lineOfData, seps); //# or ResCode
		if (token[0] != '#')
		{
			if (resID % 1000 == 0)
			{
				cout << "Compare index: " << resID << endl;
			}

			BBDepRotamer tempRotamer;

			tempRotamer.code = translate_three_character_residue_code(token);

			//token = strtok_s(NULL, seps, &context); // phi
			token = strtok(NULL, seps); // phi
			tempRotamer.phi = atoi(token);

			//token = strtok_s(NULL, seps, &context); // psi
			token = strtok(NULL, seps); // psi
			tempRotamer.psi = atoi(token);

			//token = strtok_s(NULL, seps, &context); // count
			token = strtok(NULL, seps); // count
			tempRotamer.numOfObservations = atoi(token);

			for (int chiIndex = 0; chiIndex < 4; chiIndex++)
			{
				//token = strtok_s(NULL, seps, &context); // chi rotamer
				token = strtok(NULL, seps); // chi rotamer
				tempRotamer.chiRotamer[chiIndex] = atoi(token);
			}

			//token = strtok_s(NULL, seps, &context); // probability
			token = strtok(NULL, seps); // probability
			tempRotamer.probability = atof(token);

			for (int chiIndex = 0; chiIndex < 4; chiIndex++)
			{
				//token = strtok_s(NULL, seps, &context); // chi value
				token = strtok(NULL, seps); // chi value
				tempRotamer.chi[chiIndex] = atof(token);
			}

			for (int chiIndex = 0; chiIndex < 4; chiIndex++)
			{
				//token = strtok_s(NULL, seps, &context); // stdDeviation
				token = strtok(NULL, seps); // stdDeviation
				tempRotamer.stdDeviation[chiIndex] = atof(token);
			}

			bool isCurrRotamerSame = compare_BBDepRotLib_Data(tempRotamer, BBDEP_ROTAMER_LIB_2010[resID]);
			if (isCurrRotamerSame == false)
			{
				isTwoRotLibSame = false;;
			}

			resID++;
		}
	}
	fin.close();

	return isTwoRotLibSame;
}

bool FunctionsForRotamerLibrary::compare_BBDepRotLib_Data(const BBDepRotamer& rotamer1, const BBDepRotamer& rotamer2)
{
	bool isTwoRotamerSame = true;
	const float floatCompareCutoff = 0.1;

	if (rotamer1.code != rotamer2.code)
	{
		cout << "Code is different" << endl;
		isTwoRotamerSame = false;
	}

	if (rotamer1.phi != rotamer2.phi)
	{
		cout << "Phi is different" << endl;
		isTwoRotamerSame = false;
	}

	if (rotamer1.psi != rotamer2.psi)
	{
		cout << "Psi is different" << endl;
		isTwoRotamerSame = false;
	}

	if (rotamer1.numOfObservations != rotamer2.numOfObservations)
	{
		cout << "numOfObservations is different" << endl;
		isTwoRotamerSame = false;
	}

	if (abs(rotamer1.probability - rotamer2.probability) > floatCompareCutoff)
	{
		cout << "probability is different" << endl;
		isTwoRotamerSame = false;
	}

	for (int chiIndex = 0; chiIndex < 4; chiIndex++)
	{
		if (rotamer1.chiRotamer[chiIndex] != rotamer2.chiRotamer[chiIndex])
		{
			cout << "chirotamer is different" << endl;
			isTwoRotamerSame = false;
		}

		if (abs(rotamer1.chi[chiIndex] - rotamer2.chi[chiIndex]) > floatCompareCutoff)
		{
			cout << "chi is different" << endl;
			isTwoRotamerSame = false;
		}

		if (abs(rotamer1.stdDeviation[chiIndex] - rotamer2.stdDeviation[chiIndex]) > floatCompareCutoff)
		{
			cout << "stddev is different" << endl;
			isTwoRotamerSame = false;
		}
	}

	return isTwoRotamerSame;
}

bool FunctionsForRotamerLibrary::test_BBDepLib_index_of_Residue_2010()
{
	bool isBBDepLibIndexRight = true;

	/*ResidueCode resCodeSetForTest[22] =
	{ ARG_AMINO_RESIDUE, ASN_AMINO_RESIDUE, ASP_AMINO_RESIDUE, CPR_AMINO_RESIDUE, CYD_AMINO_RESIDUE, CYH_AMINO_RESIDUE,
	CYS_AMINO_RESIDUE, GLN_AMINO_RESIDUE, GLU_AMINO_RESIDUE, HIS_AMINO_RESIDUE, ILE_AMINO_RESIDUE, LEU_AMINO_RESIDUE,
	LYS_AMINO_RESIDUE, MET_AMINO_RESIDUE, PHE_AMINO_RESIDUE, PRO_AMINO_RESIDUE, SER_AMINO_RESIDUE, THR_AMINO_RESIDUE,
	TPR_AMINO_RESIDUE, TRP_AMINO_RESIDUE, TYR_AMINO_RESIDUE, VAL_AMINO_RESIDUE };*/

	ResidueCode resCodeSetForTest[18] =
	{ ARG_AMINO_RESIDUE, ASN_AMINO_RESIDUE, ASP_AMINO_RESIDUE, CYS_AMINO_RESIDUE, GLN_AMINO_RESIDUE, GLU_AMINO_RESIDUE, HIS_AMINO_RESIDUE, ILE_AMINO_RESIDUE, LEU_AMINO_RESIDUE,
	LYS_AMINO_RESIDUE, MET_AMINO_RESIDUE, PHE_AMINO_RESIDUE, PRO_AMINO_RESIDUE, SER_AMINO_RESIDUE, THR_AMINO_RESIDUE,
	TRP_AMINO_RESIDUE, TYR_AMINO_RESIDUE, VAL_AMINO_RESIDUE };

	for (int resIndex = 0; resIndex < 18; resIndex++)
	{
		cout << "Curr res: " << resIndex << endl;
		for (int phiIndex = 0; phiIndex <= 36; phiIndex++)
		{
			for (int psiIndex = 0; psiIndex <= 36; psiIndex++)
			{
				int phi = -180 + 10 * phiIndex;
				int psi = -180 + 10 * psiIndex;

				int startIndex = 0;
				int endIndex = 0;
				getBBDepLibIndexOfResidue_2010(resCodeSetForTest[resIndex], phi, psi, startIndex, endIndex);

				for (int rotamerIndex = startIndex; rotamerIndex <= endIndex; rotamerIndex++)
				{
					if (BBDEP_ROTAMER_LIB_2010[rotamerIndex].phi != phi)
					{
						cout << "Phi is different -  resCode: " << resIndex << ", phi: " << phi << ", psi: " << psi << ", rotamerIndex: " << rotamerIndex << endl;
						isBBDepLibIndexRight = false;
					}

					if (BBDEP_ROTAMER_LIB_2010[rotamerIndex].psi != psi)
					{
						cout << "Psi is different -  resCode: " << resIndex << ", phi: " << phi << ", psi: " << psi << ", rotamerIndex: " << rotamerIndex << endl;
						isBBDepLibIndexRight = false;
					}
				}
			}
		}
	}

	return isBBDepLibIndexRight;
}


//////////////////////////////////////////////////////////////////////////
// For Dunbrack 2010 - by JHCha END
//////////////////////////////////////////////////////////////////////////

