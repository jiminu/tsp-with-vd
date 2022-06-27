#ifndef _MANAGER_OF_ROTAMERS_ASSINGNED_AT_RESIDUES_H_
#define _MANAGER_OF_ROTAMERS_ASSINGNED_AT_RESIDUES_H_

class Rotamer;
#include "SidechainConformation.h"
#include "ChainIDSeqNumRotLibIDPairsForSidechainConformation.h"

// This class represents a set of rotamers where each one is assigned to each residue in a protein.

class ManagerOfRotamers_Assigned_At_Residues
{
private:
	Rotamer**                                          m_rotamersOrderedByChainIDSeqNum                    ;  // rotamers ordered by ID (chain ID + sequence number)
	rg_INT                                             m_numRotamers                                       ;
	ChainIDSeqNumRotLibIDPairsForSidechainConformation m_chainIDSeqNumRotLibIDPairsForSidechainConformation;  // lastly updated rotamer IDs
    rg_BOOL*                                           m_recordKeepingForFixingDihedralAngles              ;  // recordKeeping for fixing the dihedral angles of sidechains

	inline void setInitialValues()
	{
		m_rotamersOrderedByChainIDSeqNum = rg_NULL;
		m_numRotamers = 0;
        m_recordKeepingForFixingDihedralAngles = rg_NULL;
	}
public:
	ManagerOfRotamers_Assigned_At_Residues();
	ManagerOfRotamers_Assigned_At_Residues(Rotamer* rotamers, const rg_INT& numRotamers);
	ManagerOfRotamers_Assigned_At_Residues(const ManagerOfRotamers_Assigned_At_Residues& managerOfRotamers_Assigned_At_Residues);
	~ManagerOfRotamers_Assigned_At_Residues();

    rg_INT getRotamersOrderedByChainIDSeqNum(Rotamer**& rotamers_in_ordered);

	inline Rotamer** getRotamersOrderedByChainIDSeqNum() 
    { 
        return m_rotamersOrderedByChainIDSeqNum; 
    }

	inline Rotamer*  getRotamer(const rg_INDEX& index) 
    { 
        return m_rotamersOrderedByChainIDSeqNum[ index ]; 
    }

	inline rg_INT    getNumRotamers() 
    { 
        return m_numRotamers; 
    }

	inline ChainIDSeqNumRotLibIDPairsForSidechainConformation*
		getChainIDSeqNumRotLibIDPairsForSidechainConformation() 
    { 
        return &m_chainIDSeqNumRotLibIDPairsForSidechainConformation; 
    }

    inline rg_BOOL* getRecordKeepingForFixingDihedralAngles()
    {
        return m_recordKeepingForFixingDihedralAngles;
    }

    void getRecordKeepingForFixingDihedralAnglesInOderOfRotamers(Rotamer** rotamers_in_ordered, rg_BOOL*& recordKeepingForFixingDihedralAngles_in_ordered);

	inline void getChainIDSeqNumRotLibIDPairsForSidechainConformation(ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
	{
		chainIDSeqNumRotLibIDPairsForSidechainConformation = m_chainIDSeqNumRotLibIDPairsForSidechainConformation;
	}

	void             getSidechainConformation(SidechainConformation& sidechainConformation);

	void             setOrderedRotamers(Rotamer** rotamers, const rg_INT& numRotamers);
	void             setOrderedRotamers(Rotamer*  rotamers, const rg_INT& numRotamers);
	void             setChainIDSeqNumRotLibPairsForSidechainConformation(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation);
    inline void      setRecordKeepingForFixingDihedralAngles()
    {
        if(m_numRotamers > 0)
        {
            m_recordKeepingForFixingDihedralAngles = new rg_BOOL[ m_numRotamers ];
            rg_INDEX i;
            for (i = 0;i < m_numRotamers;i++)
                m_recordKeepingForFixingDihedralAngles[ i ] = rg_FALSE;
        }
    }

    void             fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<V::GeometryTier::Residue*>& residuesWithFixedSidechainConformation);

    void             fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<V::GeometryTier::ResidueCode>& residueTypesWithFixedSidechainConformation);

	void             transformCoordinatesOfAtomsByApplyingAssignedRotLibIDs(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation);	

	ManagerOfRotamers_Assigned_At_Residues& operator=(const ManagerOfRotamers_Assigned_At_Residues& managerOfRotamers_Assigned_At_Residues);

private:
	void destroy();
};

#endif