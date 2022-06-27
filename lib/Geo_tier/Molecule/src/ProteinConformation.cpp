#include "ProteinConformation.h"
#include "FunctionsForMolecule.h"
using namespace V::GeometryTier;


ProteinConformation::ProteinConformation()
{
	setInitialValues();
}

ProteinConformation::ProteinConformation(ProteinConformation& proteinConformation)
{
	setInitialValues();

	set(proteinConformation.m_backboneConformation, proteinConformation.m_sidechainConformation/*, proteinConformation.m_type*/);
}

//ProteinConformation::ProteinConformation( const ProteinConformationType& type )
//{
//	m_type = type;
//}

ProteinConformation::~ProteinConformation()
{
	//m_type = UNKNOWN_TYPE;
}


void ProteinConformation::getConformationOfResidue(const rg_INT& chainIDSeqNum, 
                                                   BackboneDihedralAnglePair& backboneDihedralAnglePair, 
                                                   SidechainDihedralAngleSet& sideChainDihedralAngleSet)
{
    m_backboneConformation.getConformationOfResidue(chainIDSeqNum, backboneDihedralAnglePair);
    m_sidechainConformation.getConformationOfResidue(chainIDSeqNum, sideChainDihedralAngleSet);
}

void ProteinConformation::getConformationOfResidue(Residue* residue,  
                                                   BackboneDihedralAnglePair& backboneDihedralAnglePair, 
                                                   SidechainDihedralAngleSet& sideChainDihedralAngleSet)
{
    rg_INT chainIDSeqNum = getChainIDSeqNumOfResidue(residue);
    getConformationOfResidue(chainIDSeqNum,
                             backboneDihedralAnglePair,
                             sideChainDihedralAngleSet);
}


//rg_REAL ProteinConformation::evaluateEnergyFunction()
//{
//	rg_REAL energyVal = 0.0;
//	return energyVal;
//}
//
//void ProteinConformation::evaluateDerivativeOfEnergyFunction()
//{
//
//}



ProteinConformation& ProteinConformation::operator= (const ProteinConformation& proteinConformation)
{
	if(this == &proteinConformation)
		return *this;

	set(proteinConformation.m_backboneConformation, proteinConformation.m_sidechainConformation/*, proteinConformation.m_type*/);

	return *this;
}