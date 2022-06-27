#include "BackboneConformation.h"
#include "rg_Atom.h"
#include "Residue.h"
#include "ConstForRotamerLibrary.h"
#include "FunctionsForMolecule.h"
using namespace V::GeometryTier;

BackboneConformation::BackboneConformation()
{
}

BackboneConformation::BackboneConformation(const BackboneConformation& backboneConformation)
{
	m_backboneDihedralAnglePairs = backboneConformation.m_backboneDihedralAnglePairs;
}

BackboneConformation::~BackboneConformation()
{
}

rg_dList<BackboneDihedralAnglePair>* BackboneConformation::getBackboneDihedralAnglePairs()
{
	// sort backbone dihedral angle pairs w. r. t. their chain-sequence IDs.
	sortBackboneDihedralAnglePairs(m_backboneDihedralAnglePairs);
	return &m_backboneDihedralAnglePairs;
}

void BackboneConformation::set(Atom** backboneAtoms, const rg_INT& numBackboneAtoms)
{
	rg_INDEX i;
	for(i = 0;i < numBackboneAtoms;i++)
	{
		if(backboneAtoms[ i ]->getAtomCode() == C_ATOM &&
		   backboneAtoms[ i ]->getpChemicalProperties()->getRemoteIndicator() == ALPHA_REMOTE)
		{
			Residue* residue = backboneAtoms[ i ]->getResidue();
			rg_REAL dihedralAnglePair[ 2 ] = {DEFAULT_PHI_N_TERMINAL, DEFAULT_PSI_C_TERMINAL};
			residue->getBackBoneDihedralAngles(dihedralAnglePair[ 0 ], dihedralAnglePair[ 1 ]);
			rg_INT chainSeqID = getChainIDSeqNumOfResidue( residue );

			BackboneDihedralAnglePair backboneDihedralAnglePair(dihedralAnglePair, chainSeqID);
			add( backboneDihedralAnglePair );
		}
	}
}

void BackboneConformation::getConformationOfResidue(const rg_INT& chainIDSeqNum, BackboneDihedralAnglePair& backboneDihedralAnglePair)
{
    m_backboneDihedralAnglePairs.reset4Loop();
    while (m_backboneDihedralAnglePairs.setNext4Loop())
    {
        BackboneDihedralAnglePair* currPair          = m_backboneDihedralAnglePairs.getpEntity();
        rg_INT                     currChainIDSeqNum = currPair->getChainIDSeqNum();
        if(currChainIDSeqNum == chainIDSeqNum)
        {
            backboneDihedralAnglePair = *currPair;
            return;
        }
    }
}

void BackboneConformation::getConformationOfResidue(Residue* residue, BackboneDihedralAnglePair& backboneDihedralAnglePair)
{
    rg_INT chainIDSeqNum = getChainIDSeqNumOfResidue(residue);
    getConformationOfResidue(chainIDSeqNum, backboneDihedralAnglePair);
}

void BackboneConformation::add(const BackboneDihedralAnglePair& backboneDihedralAnglePair)
{
	m_backboneDihedralAnglePairs.add(backboneDihedralAnglePair);
}

BackboneConformation& BackboneConformation::operator=(const BackboneConformation& backboneConformation)
{
	if(this == &backboneConformation)
		return *this;

	m_backboneDihedralAnglePairs = backboneConformation.m_backboneDihedralAnglePairs;

	return *this;	
}