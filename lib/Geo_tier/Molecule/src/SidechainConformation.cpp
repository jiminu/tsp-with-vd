#include "SidechainConformation.h"
#include "FunctionsForMolecule.h"
using namespace V::GeometryTier;


SidechainConformation::SidechainConformation()
{
}

SidechainConformation::SidechainConformation(const SidechainConformation& sidechainConformation)
{
	m_sidechainDihedralAngles = sidechainConformation.m_sidechainDihedralAngles;
}

SidechainConformation::SidechainConformation(Rotamer** rotamers, const rg_INT& numRotamers)
{
    set(rotamers, numRotamers);
}

SidechainConformation::SidechainConformation(Rotamer* rotamers, const rg_INT& numRotamers)
{
    set(rotamers, numRotamers);
}

SidechainConformation::~SidechainConformation()
{
}

rg_dList<SidechainDihedralAngleSet>* SidechainConformation::getSidechainDihedralAngleSets() 
{ 
	// sort sidechain dihedral angle sets w. r. t. their chain-sequence IDs.
	sortSidechainDihedralAngleSets(m_sidechainDihedralAngles);
	return &m_sidechainDihedralAngles; 
}


void SidechainConformation::getConformationOfResidue(const rg_INT& chainIDSeqNum, SidechainDihedralAngleSet& sideChainDihedralAngleSet)
{
    m_sidechainDihedralAngles.reset4Loop();
    while (m_sidechainDihedralAngles.setNext4Loop())
    {
        SidechainDihedralAngleSet* currAngleSet = m_sidechainDihedralAngles.getpEntity();
        rg_INT currChainIDSeqNum = currAngleSet->getChainIDSeqNum();
        if(currChainIDSeqNum == chainIDSeqNum)
        {
            sideChainDihedralAngleSet = *currAngleSet;
            return;
        }
    }
}

void SidechainConformation::getConformationOfResidue(Residue* residue,            SidechainDihedralAngleSet& sideChainDihedralAngleSet)
{
    rg_INT chainIDSeqNum = getChainIDSeqNumOfResidue(residue);
    getConformationOfResidue(chainIDSeqNum, sideChainDihedralAngleSet);
}

void SidechainConformation::set(Rotamer** rotamers, const rg_INT& numRotamers)
{
	rg_INDEX i;
	for(i = 0;i < numRotamers;i++)
	{
		SidechainDihedralAngleSet sidechainDihedralAngleSet( rotamers[ i ] );
		add(sidechainDihedralAngleSet);
	}
}


void SidechainConformation::set(Rotamer* rotamers, const rg_INT& numRotamers)
{
    rg_INDEX i;
    for(i = 0;i < numRotamers;i++)
    {
        SidechainDihedralAngleSet sidechainDihedralAngleSet( rotamers[ i ] );
        add(sidechainDihedralAngleSet);
    }
}

void SidechainConformation::add(const SidechainDihedralAngleSet& sidechainDihedralAngleSet)
{
	m_sidechainDihedralAngles.add( sidechainDihedralAngleSet );
}

SidechainConformation& SidechainConformation::operator=(const SidechainConformation& sidechainConformation)
{
	if(this == &sidechainConformation)
		return *this;

	m_sidechainDihedralAngles = sidechainConformation.m_sidechainDihedralAngles;

	return *this;
}