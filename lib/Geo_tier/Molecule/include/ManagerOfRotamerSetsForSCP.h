#ifndef _MANAGER_OF_ROTAMERSETS_FOR_SCP_H_
#define _MANAGER_OF_ROTAMERSETS_FOR_SCP_H_

#include "RotamerSetOfResidue.h"
#include "ChainIDSeqNumRotLibIDPairsForSidechainConformation.h"

// This class represents for a set of rotamers sets where each set corresponds to a residue in a protein
// The same set of rotamer sets for protein design problem can be represented by using the derived class of this class.
// For the design problem, we need "Allowable residue types" for each residue position.


class ManagerOfRotamerSetsForSCP
{
private:
	RotamerSetOfResidue** m_rotamerSetsOfResidues;
	rg_INT                m_numResidues          ; // # residue positions

	inline void setInitialValues()
	{
		m_rotamerSetsOfResidues = rg_NULL;
		m_numResidues = 0;
	}

public:
	ManagerOfRotamerSetsForSCP();
	ManagerOfRotamerSetsForSCP(const rg_INT& numResidues);
	ManagerOfRotamerSetsForSCP(const ManagerOfRotamerSetsForSCP& managerOfRotamerSetsForSCP);
	~ManagerOfRotamerSetsForSCP();

	inline RotamerSetOfResidue** getRotamerSets() { return m_rotamerSetsOfResidues; }
	void                         getRotamerSets(RotamerSetOfResidue**& rotamerSetsOfResidues) const;
	inline rg_INT                getNumResidues() const { return m_numResidues;           }

	rg_INT                       getAssignedRotamerSet( Rotamer**& assignedRotamerSet );

	rg_BOOL                      isOneAndOnlyOneRotamerAssignedToEachResidue() const;

	void                         setRotamerSets(RotamerSetOfResidue** rotamerSetsOfResidues, const rg_INT& numResidues);
	void                         setRotamerSet(RotamerSetOfResidue* rotamerSetOfResidue, const rg_INDEX index);
	void                         setNumResidues(const rg_INT& numResidues);
	void                         setAssignedRotLibIDs(ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation);

    void                         fixSidechainsOfResidueInstancesWithRotamerID(
                                        rg_dList<V::GeometryTier::Residue*>& residuesWithFixedSidechainConformation, 
                                        const rg_INDEX& rotID = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX);

    void                         fixSidechainsOfResidueInstancesWithRotamerID(
                                        rg_dList<V::GeometryTier::ResidueCode>& residueTypesWithFixedSidechainConformation, 
                                        const rg_INDEX& rotID = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX);

    //void                         fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<Residue*>& residuesWithFixedSidechainConformation);

    //void                         fixSidechainsOfResidueInstancesWithCurrentConformation(rg_dList<ResidueCode>& residueTypesWithFixedSidechainConformation);

	void                         removeRotamersWithProbabilityLessThan(const rg_REAL& thresholdOfProbability);

	void                         resetAssignedRotLibIDs();

	rg_REAL                      computeRotamerPreferenceEnergyBySCWRL3(ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation);

	// energy computation for arbitrary dihedral angle values

	// dihedral angle optimization

	ManagerOfRotamerSetsForSCP& operator=(const ManagerOfRotamerSetsForSCP& managerOfRotamerSetsForSCP);

private:
	void destroy();

};

#endif