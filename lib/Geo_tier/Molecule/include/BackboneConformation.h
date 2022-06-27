#ifndef _BACKBONECONFORMATION_H_
#define _BACKBONECONFORMATION_H_

#include "ConstForMolecule.h"
#include "BackboneDihedralAnglePair.h"
#include "rg_dList.h"

#include "rg_Atom.h"
#include "Residue.h"
//class V::GeometryTier::Residue;
//class V::GeometryTier::Atom;


class BackboneConformation
{
private:
	rg_dList<BackboneDihedralAnglePair> m_backboneDihedralAnglePairs;

public:
	BackboneConformation();
	BackboneConformation(const BackboneConformation& backboneConformation);
	~BackboneConformation();

	rg_dList<BackboneDihedralAnglePair>* getBackboneDihedralAnglePairs();
    void getConformationOfResidue(const rg_INT& chainIDSeqNum, BackboneDihedralAnglePair& backboneDihedralAnglePair);
    void getConformationOfResidue(V::GeometryTier::Residue* residue, BackboneDihedralAnglePair& backboneDihedralAnglePair);
	void set(V::GeometryTier::Atom** backboneAtoms, const rg_INT& numBackboneAtoms);
	void add(const BackboneDihedralAnglePair& backboneDihedralAnglePair);
    
	BackboneConformation& operator=(const BackboneConformation& backboneConformation);
};

#endif