#include "BackboneDihedralAnglePair.h"
#include "ConstForRotamerLibrary.h"

BackboneDihedralAnglePair::BackboneDihedralAnglePair()
{
	setInitialValues();
}

BackboneDihedralAnglePair::BackboneDihedralAnglePair(const rg_REAL dihedralAnglePair[], const rg_INT& chainSeqID)
{
	setInitialValues();
	set(dihedralAnglePair, chainSeqID);
}

BackboneDihedralAnglePair::BackboneDihedralAnglePair(const BackboneDihedralAnglePair& backboneDihedralAnglePair)
{
	setInitialValues();
	set(backboneDihedralAnglePair.m_dihedralAnglePair, backboneDihedralAnglePair.m_chainIDSeqNum);
}

BackboneDihedralAnglePair::~BackboneDihedralAnglePair()
{
}

BackboneDihedralAnglePair& BackboneDihedralAnglePair::operator=(const BackboneDihedralAnglePair& backboneDihedralAnglePair)
{
	if(this == &backboneDihedralAnglePair)
		return *this;

	set(backboneDihedralAnglePair.m_dihedralAnglePair, backboneDihedralAnglePair.m_chainIDSeqNum);

	return *this;
}