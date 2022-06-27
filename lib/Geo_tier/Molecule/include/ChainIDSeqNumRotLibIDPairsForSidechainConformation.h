#ifndef _ROTAMERLIBINDICESFORSIDECHAINCONFORMATION_H_
#define _ROTAMERLIBINDICESFORSIDECHAINCONFORMATION_H_

#include "ChainIDSeqNumRotLibIDPair.h"
#include "rg_dList.h"

class ChainIDSeqNumRotLibIDPairsForSidechainConformation
{
private:
	rg_dList<ChainIDSeqNumRotLibIDPair> m_chainIDSeqNumRotLibIDPair;

public:
	ChainIDSeqNumRotLibIDPairsForSidechainConformation();
	ChainIDSeqNumRotLibIDPairsForSidechainConformation(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation);
	~ChainIDSeqNumRotLibIDPairsForSidechainConformation();

	ChainIDSeqNumRotLibIDPair* getChainIDSeqNumRotLibIDPairCorrTo(const rg_INT& chainIDSeqNum);
	rg_INT                     getRotIDCorrToChainIDSeqNum(const rg_INT& chainIDSeqNum);
	inline rg_dList<ChainIDSeqNumRotLibIDPair>* getChainIDSeqNumRotLibIDPairs()
    {
        return &m_chainIDSeqNumRotLibIDPair;
    }
	rg_INT                     getChainIDSeqNumRotLibIDPairs(ChainIDSeqNumRotLibIDPair*& chainIDSeqNumRotLibIDPair);
	void add(const ChainIDSeqNumRotLibIDPair& chainIDSeqNumRotLibIDPair);

	void setRotLibID(const rg_INT& chainIDSeqNum, const rg_INT& rotLibID);
	void sortChainIDSeqNumRotLibIDPairs();

	ChainIDSeqNumRotLibIDPairsForSidechainConformation& operator=(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation);

};

#endif