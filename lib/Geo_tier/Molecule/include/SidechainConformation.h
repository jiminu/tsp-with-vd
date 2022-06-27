#ifndef _SIDECHAINCONFORMATION_H_
#define _SIDECHAINCONFORMATION_H_

#include "SidechainDihedralAngleSet.h"
#include "rg_dList.h"

class V::GeometryTier::Residue;

class SidechainConformation
{
private:
	rg_dList<SidechainDihedralAngleSet> m_sidechainDihedralAngles;

public:
	SidechainConformation();
	SidechainConformation(const SidechainConformation& sidechainConformation);
    SidechainConformation(Rotamer** rotamers, const rg_INT& numRotamers);
    SidechainConformation(Rotamer* rotamers, const rg_INT& numRotamers);
	~SidechainConformation();

	rg_dList<SidechainDihedralAngleSet>* getSidechainDihedralAngleSets();
    void getConformationOfResidue(const rg_INT& chainIDSeqNum, SidechainDihedralAngleSet& sideChainDihedralAngleSet);
    void getConformationOfResidue(V::GeometryTier::Residue* residue,            SidechainDihedralAngleSet& sideChainDihedralAngleSet);
	void set(Rotamer** rotamers, const rg_INT& numRotamers);
    void set(Rotamer* rotamers, const rg_INT& numRotamers);
	void add(const SidechainDihedralAngleSet& sidechainDihedralAngleSet);

	SidechainConformation& operator=(const SidechainConformation& sidechainConformation);
};

#endif