#ifndef _SIDECHAINDIHEDRALANGLESET_H_
#define _SIDECHAINDIHEDRALANGLESET_H_

// This class represents a set of dihedral angles 
// which define the conformation of a sidechain.

#include "Rotamer.h"

class SidechainDihedralAngleSet
{
private:
	rg_REAL* m_dihedralAngles;
	rg_INT   m_numDihedralAngles;
	rg_INT   m_chainIDSeqNum;

	inline void setInitialValues()
	{
		m_dihedralAngles    = rg_NULL;
		m_numDihedralAngles = 0;
		m_chainIDSeqNum     = -1;
	}

public:
	SidechainDihedralAngleSet();
	SidechainDihedralAngleSet(Rotamer& rotamer);
	SidechainDihedralAngleSet(Rotamer* rotamer);
	SidechainDihedralAngleSet(const SidechainDihedralAngleSet& sidechainDihedralAngleSet);
	~SidechainDihedralAngleSet();

	inline rg_REAL* getDihedralAngles() { return m_dihedralAngles; }
	inline rg_INT   getNumDihedralAngles() const { return m_numDihedralAngles; }
	inline rg_INT   getChainIDSeqNum() const { return m_chainIDSeqNum; }

	void            setDihedralAngles( rg_REAL* dihedralAngles, const rg_INT& numDihedralAngles );
	inline void     setNumDihedralAngles( const rg_INT& numDIhedralAngles ) { m_numDihedralAngles = numDIhedralAngles; }
	inline void     setChainSeqID( const rg_INT& chainSeqID ) { m_chainIDSeqNum = chainSeqID; }

	SidechainDihedralAngleSet& operator=(const SidechainDihedralAngleSet& sidechainDihedralAngleSet);

private:
	void destroyDihedralAngles();
};

#endif