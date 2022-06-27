#ifndef _BACKBONEDIHEDRALANGLEPAIR_H_
#define _BACKBONEDIHEDRALANGLEPAIR_H_

#include "ConstForRotamerLibrary.h"

class BackboneDihedralAnglePair
{
private:
	rg_REAL m_dihedralAnglePair[ 2 ]; // phi and psi angles
	rg_INT  m_chainIDSeqNum         ; // chain ID + sequence number
	inline void setInitialValues()
	{
		m_dihedralAnglePair[ 0 ] = DEFAULT_PHI_N_TERMINAL;
		m_dihedralAnglePair[ 1 ] = DEFAULT_PSI_C_TERMINAL;
		m_chainIDSeqNum             = -1;
	}
	inline void set(const rg_REAL dihedralAnglePair[], const rg_INT& chainSeqID)
	{
		m_dihedralAnglePair[ 0 ] = dihedralAnglePair[ 0 ];
		m_dihedralAnglePair[ 1 ] = dihedralAnglePair[ 1 ];
		m_chainIDSeqNum = chainSeqID;
	}

public:
	BackboneDihedralAnglePair();
	BackboneDihedralAnglePair(const rg_REAL dihedralAnglePair[], const rg_INT& chainSeqID);
	BackboneDihedralAnglePair(const BackboneDihedralAnglePair& backboneDihedralAnglePair);
	~BackboneDihedralAnglePair();

	inline void   getDihedralAnglePair(rg_REAL dihedralAnglePair[])
	{
		dihedralAnglePair[ 0 ] = m_dihedralAnglePair[ 0 ];
		dihedralAnglePair[ 1 ] = m_dihedralAnglePair[ 1 ];
	}
	inline rg_INT getChainIDSeqNum() const { return m_chainIDSeqNum;}
	inline void   setDihedralAnglePair(const rg_REAL dihedralAnglePair[])
	{
		m_dihedralAnglePair[ 0 ] = dihedralAnglePair[ 0 ];
		m_dihedralAnglePair[ 1 ] = dihedralAnglePair[ 1 ];
	}
	inline void   setChainSeqID(const rg_INT& chainSeqID) { m_chainIDSeqNum = chainSeqID; }

	BackboneDihedralAnglePair& operator=(const BackboneDihedralAnglePair& backboneDihedralAnglePair);

};

#endif 