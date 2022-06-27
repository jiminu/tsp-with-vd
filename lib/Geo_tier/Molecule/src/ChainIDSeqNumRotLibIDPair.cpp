#include "ChainIDSeqNumRotLibIDPair.h"
#include "ConstForRotamerLibrary.h"

ChainIDSeqNumRotLibIDPair::ChainIDSeqNumRotLibIDPair()
{
	setInitialValues();
}

ChainIDSeqNumRotLibIDPair::ChainIDSeqNumRotLibIDPair(const rg_INT& chainIDSeqNum, const rg_INT& rotLibID)
{
	setInitialValues();
	set(chainIDSeqNum, rotLibID);
}

ChainIDSeqNumRotLibIDPair::ChainIDSeqNumRotLibIDPair(const ChainIDSeqNumRotLibIDPair& chainIDSeqNumRotLibIDPair)
{
	setInitialValues();
	set(chainIDSeqNumRotLibIDPair.m_chainIDSeqNum, chainIDSeqNumRotLibIDPair.m_rotLibID);
}

ChainIDSeqNumRotLibIDPair::~ChainIDSeqNumRotLibIDPair()
{
}

ChainIDSeqNumRotLibIDPair& ChainIDSeqNumRotLibIDPair::operator=(const ChainIDSeqNumRotLibIDPair& chainIDSeqNumRotLibIDPair)
{
	if (this == &chainIDSeqNumRotLibIDPair)
	{
		return *this;
	}

	set(chainIDSeqNumRotLibIDPair.m_chainIDSeqNum, chainIDSeqNumRotLibIDPair.m_rotLibID);

	return *this;
}