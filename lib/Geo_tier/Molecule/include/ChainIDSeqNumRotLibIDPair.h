#ifndef _CHAINIDSEQNUMROTLIBIDPAIR_H_
#define _CHAINIDSEQNUMROTLIBIDPAIR_H_

#include "ConstForRotamerLibrary.h"

class ChainIDSeqNumRotLibIDPair
{
private:
	rg_INT m_chainIDSeqNum;
	rg_INT m_rotLibID     ;
	inline void setInitialValues()
	{
		m_chainIDSeqNum = -1;
		m_rotLibID      = UNKNOWN_ROT_INDEX;
	}
	inline void set(const rg_INT& chainIDSeqNum, const rg_INT& rotLibID)
	{
		m_chainIDSeqNum = chainIDSeqNum;
		m_rotLibID      = rotLibID;
	}

public:
	ChainIDSeqNumRotLibIDPair();
	ChainIDSeqNumRotLibIDPair(const rg_INT& chainIDSeqNum, const rg_INT& rotLibID);
	ChainIDSeqNumRotLibIDPair(const ChainIDSeqNumRotLibIDPair& chainIDSeqNumRotLibIDPair);
	~ChainIDSeqNumRotLibIDPair();

	inline rg_INT getChainIDSeqNum() const { return m_chainIDSeqNum; }
	inline rg_INT getRotLibID()      const { return m_rotLibID     ; }

	inline void   setChainIDSeqNum( const rg_INT& chainIDSeqNum ) { m_chainIDSeqNum = chainIDSeqNum; }
	inline void   setRotLibID     ( const rg_INT& rotLibID      ) { m_rotLibID      = rotLibID     ; }

	ChainIDSeqNumRotLibIDPair& operator=(const ChainIDSeqNumRotLibIDPair& chainIDSeqNumRotLibIDPair);
};

#endif