#include "ManagerOfRotamers_Assigned_At_Residues.h"
#include "Rotamer.h"
#include "ConstForRotamerLibrary.h"
#include "FunctionsForMolecule.h"
using namespace V::GeometryTier;

ManagerOfRotamers_Assigned_At_Residues::ManagerOfRotamers_Assigned_At_Residues()
{
	setInitialValues();
}

ManagerOfRotamers_Assigned_At_Residues::ManagerOfRotamers_Assigned_At_Residues(Rotamer* rotamers, const rg_INT& numRotamers)
{
	setInitialValues();
	setOrderedRotamers(rotamers, numRotamers);    
}

ManagerOfRotamers_Assigned_At_Residues::ManagerOfRotamers_Assigned_At_Residues(const ManagerOfRotamers_Assigned_At_Residues& managerOfRotamers_Assigned_At_Residues)
{
	setInitialValues();
	setOrderedRotamers(managerOfRotamers_Assigned_At_Residues.m_rotamersOrderedByChainIDSeqNum, managerOfRotamers_Assigned_At_Residues.m_numRotamers);
	m_chainIDSeqNumRotLibIDPairsForSidechainConformation = managerOfRotamers_Assigned_At_Residues.m_chainIDSeqNumRotLibIDPairsForSidechainConformation;
}

ManagerOfRotamers_Assigned_At_Residues::~ManagerOfRotamers_Assigned_At_Residues()
{
	destroy();
}

rg_INT ManagerOfRotamers_Assigned_At_Residues::getRotamersOrderedByChainIDSeqNum(Rotamer**& rotamers_in_ordered) 
{
    rg_INDEX i;
    for (i = 0;i < m_numRotamers;i++)
    {
        rotamers_in_ordered[ i ] = m_rotamersOrderedByChainIDSeqNum[ i ];
    }

    return m_numRotamers;
}

void ManagerOfRotamers_Assigned_At_Residues::getRecordKeepingForFixingDihedralAnglesInOderOfRotamers(Rotamer** rotamers_in_ordered, rg_BOOL*& recordKeepingForFixingDihedralAngles_in_ordered)
{
    rg_INDEX i;
    for(i = 0;i < m_numRotamers;i++)
        recordKeepingForFixingDihedralAngles_in_ordered[ i ] = rg_FALSE;

    for (i = 0;i < m_numRotamers;i++)
    {
        if(m_recordKeepingForFixingDihedralAngles[ i ])
        {
            rg_INT chainIDSeqNumOfCurrResidueWithFixedConformation = getChainIDSeqNumOfResidue( m_rotamersOrderedByChainIDSeqNum[ i ]->getResidue() );
            rg_INDEX j;
            for (j = 0;j < m_numRotamers;j++)
            {
                rg_INT currChainIDSeqNum = getChainIDSeqNumOfResidue( rotamers_in_ordered[ j ]->getResidue() );
                if(chainIDSeqNumOfCurrResidueWithFixedConformation == currChainIDSeqNum)
                {
                    recordKeepingForFixingDihedralAngles_in_ordered[ j ] = rg_TRUE;
                    break;
                }
            }
        }        
    }
}

void     ManagerOfRotamers_Assigned_At_Residues::getSidechainConformation(SidechainConformation& sidechainConformation)
{
	sidechainConformation.set(m_rotamersOrderedByChainIDSeqNum, m_numRotamers);
}

void     ManagerOfRotamers_Assigned_At_Residues::setOrderedRotamers(Rotamer** rotamers, const rg_INT& numRotamers)
{
	destroy();

	m_numRotamers = numRotamers;
	m_rotamersOrderedByChainIDSeqNum = new Rotamer* [ m_numRotamers ];
	rg_INDEX i;
	for (i = 0;i < m_numRotamers;i++)
	{
		m_rotamersOrderedByChainIDSeqNum[ i ] = rotamers[ i ];
	}
    setRecordKeepingForFixingDihedralAngles();
}

void     ManagerOfRotamers_Assigned_At_Residues::setOrderedRotamers(Rotamer* rotamers, const rg_INT& numRotamers)
{
	destroy();

	m_numRotamers = numRotamers;
	m_rotamersOrderedByChainIDSeqNum = new Rotamer* [ m_numRotamers ];
	rg_INDEX i;
	for (i = 0;i < m_numRotamers;i++)
	{
		m_rotamersOrderedByChainIDSeqNum[ i ] = &rotamers[ i ];
	}
    setRecordKeepingForFixingDihedralAngles();
}

void ManagerOfRotamers_Assigned_At_Residues::setChainIDSeqNumRotLibPairsForSidechainConformation(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
{ 
	m_chainIDSeqNumRotLibIDPairsForSidechainConformation = chainIDSeqNumRotLibIDPairsForSidechainConformation;
	m_chainIDSeqNumRotLibIDPairsForSidechainConformation.sortChainIDSeqNumRotLibIDPairs();
}

void ManagerOfRotamers_Assigned_At_Residues::fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<Residue*>& residuesWithFixedSidechainConformation)
{
    rg_INDEX i;
    for (i = 0;i < m_numRotamers;i++)
    {
        Residue* currResidue = m_rotamersOrderedByChainIDSeqNum[ i ]->getResidue();
        rg_INT chainIDSeqNumOfCurrResidue = getChainIDSeqNumOfResidue( currResidue );

        residuesWithFixedSidechainConformation.reset4Loop();
        while (residuesWithFixedSidechainConformation.setNext4Loop())
        {
            Residue* currResidueToBeFixed = residuesWithFixedSidechainConformation.getEntity();
            rg_INT chainIDSeqNumOfCurrResidueToBeFixed = getChainIDSeqNumOfResidue( currResidueToBeFixed );

            if(chainIDSeqNumOfCurrResidue == chainIDSeqNumOfCurrResidueToBeFixed)
            {
                m_recordKeepingForFixingDihedralAngles[ i ] = rg_TRUE;
                //residuesToBeFixedWithCurrentSidechainConformation.killCurrent();
                break;
            }
        }
    }
}

void ManagerOfRotamers_Assigned_At_Residues::fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<ResidueCode>& residueTypesWithFixedSidechainConformation)
{
    rg_INDEX i;
    for (i = 0;i < m_numRotamers;i++)
    {
        Residue* currResidue = m_rotamersOrderedByChainIDSeqNum[ i ]->getResidue();
        ResidueCode currCode = currResidue->getResidueCode();

        residueTypesWithFixedSidechainConformation.reset4Loop();
        while (residueTypesWithFixedSidechainConformation.setNext4Loop())
        {
            ResidueCode currResidueCodeWithFixedSidechainConformation = residueTypesWithFixedSidechainConformation.getEntity();
            if (currCode == currResidueCodeWithFixedSidechainConformation)
            {
                m_recordKeepingForFixingDihedralAngles[ i ] = rg_TRUE;
                //residuesToBeFixedWithCurrentSidechainConformation.killCurrent();
                break;
            }
        }
    }
}


void ManagerOfRotamers_Assigned_At_Residues::transformCoordinatesOfAtomsByApplyingAssignedRotLibIDs(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
{
	setChainIDSeqNumRotLibPairsForSidechainConformation( chainIDSeqNumRotLibIDPairsForSidechainConformation );

	rg_dList<ChainIDSeqNumRotLibIDPair>* chainIDSeqNumRotLibIDPairs = rg_NULL;
	chainIDSeqNumRotLibIDPairs = m_chainIDSeqNumRotLibIDPairsForSidechainConformation.getChainIDSeqNumRotLibIDPairs();

	rg_INDEX index = 0;
	chainIDSeqNumRotLibIDPairs->reset4Loop();
	while (chainIDSeqNumRotLibIDPairs->setNext4Loop())
	{
		ChainIDSeqNumRotLibIDPair* currPair = chainIDSeqNumRotLibIDPairs->getpEntity();
		rg_INDEX rotLibID = UNKNOWN_ROT_INDEX;
		rotLibID = currPair->getRotLibID();
		m_rotamersOrderedByChainIDSeqNum[ index ]->updateSidechainAtomCoordinatesByApplyingRotLibID( rotLibID );
		index++;
	}
}

ManagerOfRotamers_Assigned_At_Residues& ManagerOfRotamers_Assigned_At_Residues::operator=(const ManagerOfRotamers_Assigned_At_Residues& managerOfRotamers_Assigned_At_Residues)
{
	if(this == &managerOfRotamers_Assigned_At_Residues)
		return *this;

	setOrderedRotamers(managerOfRotamers_Assigned_At_Residues.m_rotamersOrderedByChainIDSeqNum, managerOfRotamers_Assigned_At_Residues.m_numRotamers);
	m_chainIDSeqNumRotLibIDPairsForSidechainConformation = managerOfRotamers_Assigned_At_Residues.m_chainIDSeqNumRotLibIDPairsForSidechainConformation;

	return *this;
}

void ManagerOfRotamers_Assigned_At_Residues::destroy()
{
	if(m_numRotamers > 0 && m_rotamersOrderedByChainIDSeqNum != rg_NULL && m_recordKeepingForFixingDihedralAngles != rg_NULL)
	{
		delete [] m_rotamersOrderedByChainIDSeqNum;
		m_rotamersOrderedByChainIDSeqNum = rg_NULL;
        delete [] m_recordKeepingForFixingDihedralAngles;
        m_recordKeepingForFixingDihedralAngles = rg_NULL;
		m_numRotamers = 0;
	}
}