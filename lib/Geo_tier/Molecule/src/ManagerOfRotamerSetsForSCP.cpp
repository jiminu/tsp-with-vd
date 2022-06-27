#include "ManagerOfRotamerSetsForSCP.h"
#include "Rotamer.h"
#include "FunctionsForMolecule.h"
using namespace V::GeometryTier;

ManagerOfRotamerSetsForSCP::ManagerOfRotamerSetsForSCP()
{
	setInitialValues();
}

ManagerOfRotamerSetsForSCP::ManagerOfRotamerSetsForSCP(const rg_INT& numResidues)
{
	setInitialValues();

	setNumResidues( numResidues );
}

ManagerOfRotamerSetsForSCP::ManagerOfRotamerSetsForSCP(const ManagerOfRotamerSetsForSCP& managerOfRotamerSetsForSCP)
{
	setInitialValues();

	setRotamerSets(managerOfRotamerSetsForSCP.m_rotamerSetsOfResidues, managerOfRotamerSetsForSCP.m_numResidues);
}

ManagerOfRotamerSetsForSCP::~ManagerOfRotamerSetsForSCP()
{
	destroy();
}

void ManagerOfRotamerSetsForSCP::getRotamerSets(RotamerSetOfResidue**& rotamerSetsOfResidues) const
{
	rotamerSetsOfResidues = new RotamerSetOfResidue*[m_numResidues];
	rg_INDEX i;
	for (i = 0;i < m_numResidues;i++)
	{
		rotamerSetsOfResidues[ i ] = m_rotamerSetsOfResidues[ i ];
	}
}

rg_INT ManagerOfRotamerSetsForSCP::getAssignedRotamerSet( Rotamer**& assignedRotamerSet )
{
	//assignedRotamerSet = new Rotamer* [m_numResidues];
	rg_INDEX i;
	for (i = 0;i < m_numResidues;i++)
	{
		assignedRotamerSet[ i ] = m_rotamerSetsOfResidues[ i ]->getRotamer();
	}

	return m_numResidues;
}

rg_BOOL ManagerOfRotamerSetsForSCP::isOneAndOnlyOneRotamerAssignedToEachResidue() const
{
	rg_BOOL bOneAndOnlyOneRotamerAssigned = rg_TRUE;

	rg_INDEX i;
	for (i = 0;i < m_numResidues;i++)
	{
		rg_INDEX assignedRotLibID = UNKNOWN_ROT_INDEX;
		assignedRotLibID = m_rotamerSetsOfResidues[ i ]->getAssignedRotLibID();
		//if(assignedRotLibID == UNKNOWN_ROT_INDEX || assignedRotLibID == FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX)
        if(assignedRotLibID == UNKNOWN_ROT_INDEX)
		{
			bOneAndOnlyOneRotamerAssigned = rg_FALSE;
			break;
		}
	}

	return bOneAndOnlyOneRotamerAssigned;
}

void ManagerOfRotamerSetsForSCP::setRotamerSets(RotamerSetOfResidue** rotamerSetsOfResidues, const rg_INT& numResidues)
{
	destroy();

	m_numResidues = numResidues;
	m_rotamerSetsOfResidues = new RotamerSetOfResidue* [m_numResidues];
	rg_INDEX i;
	for (i = 0;i < m_numResidues;i++)
	{
		m_rotamerSetsOfResidues[ i ] = rotamerSetsOfResidues[ i ];
	}
}

void ManagerOfRotamerSetsForSCP::setRotamerSet(RotamerSetOfResidue* rotamerSetOfResidue, const rg_INDEX index)
{
	if(index >= m_numResidues)
		return;

	m_rotamerSetsOfResidues[ index ] = rotamerSetOfResidue;
}

void ManagerOfRotamerSetsForSCP::setNumResidues(const rg_INT& numResidues)
{
	destroy();
	m_numResidues = numResidues;
	m_rotamerSetsOfResidues = new RotamerSetOfResidue* [m_numResidues];
}

void ManagerOfRotamerSetsForSCP::setAssignedRotLibIDs(ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
{
	chainIDSeqNumRotLibIDPairsForSidechainConformation.sortChainIDSeqNumRotLibIDPairs();

	rg_dList<ChainIDSeqNumRotLibIDPair>* pairs =
		chainIDSeqNumRotLibIDPairsForSidechainConformation.getChainIDSeqNumRotLibIDPairs();

	rg_INDEX i = 0;
	pairs->reset4Loop();
	while (pairs->setNext4Loop() && i < m_numResidues)
	{
		ChainIDSeqNumRotLibIDPair* currPair = pairs->getpEntity();
		rg_INDEX rotLibID = UNKNOWN_ROT_INDEX;
		rotLibID = currPair->getRotLibID();

		m_rotamerSetsOfResidues[ i ]->setAssignedRotLibID( rotLibID );
		i++;
	}
}

void ManagerOfRotamerSetsForSCP::fixSidechainsOfResidueInstancesWithRotamerID(rg_dList<Residue*>& residuesWithFixedSidechainConformation, const rg_INDEX& rotID /*= FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX*/)
{
    rg_INDEX i;
    for (i = 0;i < m_numResidues;i++)
    {
        Residue* currResidue = m_rotamerSetsOfResidues[ i ]->getRotamer()->getResidue();
        rg_INT chainIDSeqNumOfCurrResidue = getChainIDSeqNumOfResidue( currResidue );

        residuesWithFixedSidechainConformation.reset4Loop();
        while (residuesWithFixedSidechainConformation.setNext4Loop())
        {
            Residue* currResidueToBeFixed = residuesWithFixedSidechainConformation.getEntity();
            rg_INT chainIDSeqNumOfCurrResidueToBeFixed = getChainIDSeqNumOfResidue( currResidueToBeFixed );
            if(chainIDSeqNumOfCurrResidue == chainIDSeqNumOfCurrResidueToBeFixed)
            {
                m_rotamerSetsOfResidues[ i ]->fixSidechainConformationWithRotamerID( rotID );
                //m_rotamerSetsOfResidues[ i ]->setAssignedRotLibID( rotID );
                //m_rotamerSetsOfResidues[ i ]->setFixSidechainFlag( rg_TRUE );                
                //m_rotamerSetsOfResidues[ i ]->setRotLibIDs(FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX, FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX);
                //residuesToBeFixedWithCurrentSidechainConformation.killCurrent();
                break;
            }
        }
    }
}


void ManagerOfRotamerSetsForSCP::fixSidechainsOfResidueInstancesWithRotamerID(rg_dList<ResidueCode>& residueTypesWithFixedSidechainConformation, const rg_INDEX& rotID /*= FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX*/)
{
    rg_INDEX i;
    for (i = 0;i < m_numResidues;i++)
    {
        Residue* currResidue = m_rotamerSetsOfResidues[ i ]->getRotamer()->getResidue();
        ResidueCode currCode = currResidue->getResidueCode();

        residueTypesWithFixedSidechainConformation.reset4Loop();
        while (residueTypesWithFixedSidechainConformation.setNext4Loop())
        {
            ResidueCode currResidueCodeWithFixedSidechainConformation = residueTypesWithFixedSidechainConformation.getEntity();
            if (currCode == currResidueCodeWithFixedSidechainConformation)
            {
                m_rotamerSetsOfResidues[ i ]->fixSidechainConformationWithRotamerID( rotID );
                break;
            }
        }
    }
}


//void ManagerOfRotamerSetsForSCP::fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<Residue*>& residuesWithFixedSidechainConformation)
//{
//    rg_INDEX i;
//    for (i = 0;i < m_numResidues;i++)
//    {
//        Residue* currResidue = m_rotamerSetsOfResidues[ i ]->getRotamer()->getResidue();
//        rg_INT chainIDSeqNumOfCurrResidue = getChainIDSeqNumOfResidue( currResidue );
//
//        residuesWithFixedSidechainConformation.reset4Loop();
//        while (residuesWithFixedSidechainConformation.setNext4Loop())
//        {
//            Residue* currResidueToBeFixed = residuesWithFixedSidechainConformation.getEntity();
//            rg_INT chainIDSeqNumOfCurrResidueToBeFixed = getChainIDSeqNumOfResidue( currResidueToBeFixed );
//            if(chainIDSeqNumOfCurrResidue == chainIDSeqNumOfCurrResidueToBeFixed)
//            {
//                m_rotamerSetsOfResidues[ i ]->setFixSidechainFlag( rg_TRUE );
//                break;
//            }
//        }
//    }
//}


//void ManagerOfRotamerSetsForSCP::fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<ResidueCode>& residueTypesWithFixedSidechainConformation)
//{
//    rg_INDEX i;
//    for (i = 0;i < m_numResidues;i++)
//    {
//        Residue* currResidue = m_rotamerSetsOfResidues[ i ]->getRotamer()->getResidue();
//        ResidueCode currCode = currResidue->getResidueCode();
//
//        residueTypesWithFixedSidechainConformation.reset4Loop();
//        while (residueTypesWithFixedSidechainConformation.setNext4Loop())
//        {
//            ResidueCode currResidueCodeWithFixedSidechainConformation = residueTypesWithFixedSidechainConformation.getEntity();
//            if (currCode == currResidueCodeWithFixedSidechainConformation)
//            {
//                m_rotamerSetsOfResidues[ i ]->setFixSidechainFlag( rg_TRUE );
//                break;
//            }
//        }
//    }
//}


void ManagerOfRotamerSetsForSCP::removeRotamersWithProbabilityLessThan(const rg_REAL& thresholdOfProbability)
{
	rg_INDEX i;
	for (i = 0;i < m_numResidues;i++)
	{
		m_rotamerSetsOfResidues[ i ]->removeRotamersWithProbabilityLessThan( thresholdOfProbability );
	}
}

void ManagerOfRotamerSetsForSCP::resetAssignedRotLibIDs()
{
	rg_INDEX i;
	for (i = 0;i < m_numResidues;i++)
	{
		m_rotamerSetsOfResidues[ i ]->resetAssignedRotLibID();
	}
}

rg_REAL ManagerOfRotamerSetsForSCP::computeRotamerPreferenceEnergyBySCWRL3(ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
{
	rg_dList<ChainIDSeqNumRotLibIDPair>* pairs =
		chainIDSeqNumRotLibIDPairsForSidechainConformation.getChainIDSeqNumRotLibIDPairs();

	rg_REAL energy = 0.0;

	rg_INDEX i = 0;
	pairs->reset4Loop();
	while (pairs->setNext4Loop() && i < m_numResidues)
	{
		ChainIDSeqNumRotLibIDPair* currPair = pairs->getpEntity();

		rg_INDEX rotLibID = UNKNOWN_ROT_INDEX;
		rotLibID = currPair->getRotLibID();
		rg_INDEX highestProbabilityRotLibID = UNKNOWN_ROT_INDEX;
		highestProbabilityRotLibID = m_rotamerSetsOfResidues[ i ]->getRotLibIDWithHighestProbability();
		Rotamer* currRotamer = m_rotamerSetsOfResidues[ i ]->getRotamer();

		rg_REAL currEnergy = 0.0;
		currEnergy = currRotamer->computeRotamerPreferenceEnergyBySCWRL3(rotLibID, highestProbabilityRotLibID);
		energy += currEnergy;
	}

	return energy;
}

ManagerOfRotamerSetsForSCP& ManagerOfRotamerSetsForSCP::operator=(const ManagerOfRotamerSetsForSCP& managerOfRotamerSetsForSCP)
{
	if(this == &managerOfRotamerSetsForSCP)
		return *this;

	setRotamerSets(managerOfRotamerSetsForSCP.m_rotamerSetsOfResidues, managerOfRotamerSetsForSCP.m_numResidues);

	return *this;
}

void ManagerOfRotamerSetsForSCP::destroy()
{
	if(m_numResidues > 0 && m_rotamerSetsOfResidues != rg_NULL)
	{
		delete [] m_rotamerSetsOfResidues;
		m_rotamerSetsOfResidues = rg_NULL;
		m_numResidues = 0;
	}
}