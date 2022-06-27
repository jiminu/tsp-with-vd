#include "ChainIDSeqNumRotLibIDPairsForSidechainConformation.h"
#include "FunctionsForMolecule.h"
#include "ConstForRotamerLibrary.h"
using namespace V::GeometryTier;

ChainIDSeqNumRotLibIDPairsForSidechainConformation::ChainIDSeqNumRotLibIDPairsForSidechainConformation()
{
}

ChainIDSeqNumRotLibIDPairsForSidechainConformation::ChainIDSeqNumRotLibIDPairsForSidechainConformation(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
{
	m_chainIDSeqNumRotLibIDPair = chainIDSeqNumRotLibIDPairsForSidechainConformation.m_chainIDSeqNumRotLibIDPair;
}

ChainIDSeqNumRotLibIDPairsForSidechainConformation::~ChainIDSeqNumRotLibIDPairsForSidechainConformation()
{
}

ChainIDSeqNumRotLibIDPair* ChainIDSeqNumRotLibIDPairsForSidechainConformation::getChainIDSeqNumRotLibIDPairCorrTo(const rg_INT& chainIDSeqNum)
{
	ChainIDSeqNumRotLibIDPair* pair = rg_NULL;
	// linear search
	m_chainIDSeqNumRotLibIDPair.reset4Loop();
	while (m_chainIDSeqNumRotLibIDPair.setNext4Loop())
	{
		ChainIDSeqNumRotLibIDPair* currPair = m_chainIDSeqNumRotLibIDPair.getpEntity();
		if(currPair->getChainIDSeqNum() == chainIDSeqNum)
		{
			pair = currPair;
			break;
		}
	}
	return pair;
}

rg_INT ChainIDSeqNumRotLibIDPairsForSidechainConformation::getRotIDCorrToChainIDSeqNum(const rg_INT& chainIDSeqNum)
{
	rg_INT rotID = UNKNOWN_ROT_INDEX;
	ChainIDSeqNumRotLibIDPair* pair = getChainIDSeqNumRotLibIDPairCorrTo( chainIDSeqNum );
	if(pair != rg_NULL)
		rotID = pair->getRotLibID();

	return rotID;
}

//rg_dList<ChainIDSeqNumRotLibIDPair>* ChainIDSeqNumRotLibIDPairsForSidechainConformation::getChainIDSeqNumRotLibIDPairs() 
//{ 
//	//sortChainIDSeqNumRotLibIDPairs();
//	return &m_chainIDSeqNumRotLibIDPair; 
//}

rg_INT ChainIDSeqNumRotLibIDPairsForSidechainConformation::getChainIDSeqNumRotLibIDPairs(ChainIDSeqNumRotLibIDPair*& chainIDSeqNumRotLibIDPair)
{
	//sortChainIDSeqNumRotLibIDPairs();

	rg_INT numPairs = m_chainIDSeqNumRotLibIDPair.getSize();
	chainIDSeqNumRotLibIDPair = new ChainIDSeqNumRotLibIDPair[ numPairs ];

	rg_INDEX i = 0;
	m_chainIDSeqNumRotLibIDPair.reset4Loop();
	while (m_chainIDSeqNumRotLibIDPair.setNext4Loop())
	{
		chainIDSeqNumRotLibIDPair[ i ] = m_chainIDSeqNumRotLibIDPair.getEntity();
		i++;
	}

	return numPairs;
}

void ChainIDSeqNumRotLibIDPairsForSidechainConformation::add(const ChainIDSeqNumRotLibIDPair& chainIDSeqNumRotLibIDPair)
{
	m_chainIDSeqNumRotLibIDPair.add( chainIDSeqNumRotLibIDPair );
}

void ChainIDSeqNumRotLibIDPairsForSidechainConformation::setRotLibID(const rg_INT& chainIDSeqNum, const rg_INT& rotLibID)
{
	// linear search
	m_chainIDSeqNumRotLibIDPair.reset4Loop();
	while (m_chainIDSeqNumRotLibIDPair.setNext4Loop())
	{
		ChainIDSeqNumRotLibIDPair* currPair = m_chainIDSeqNumRotLibIDPair.getpEntity();
		if(currPair->getChainIDSeqNum() == chainIDSeqNum)
		{
			currPair->setRotLibID( rotLibID );
			return;
		}
	}
}

void ChainIDSeqNumRotLibIDPairsForSidechainConformation::sortChainIDSeqNumRotLibIDPairs()
{
	// sort chainIDSeqNum rotamer library ID pairs w. r. t. their chain-sequence IDs.
	sortChainIDSeqNumRotLibIDPairSet(m_chainIDSeqNumRotLibIDPair);
}

ChainIDSeqNumRotLibIDPairsForSidechainConformation& ChainIDSeqNumRotLibIDPairsForSidechainConformation::operator=(const ChainIDSeqNumRotLibIDPairsForSidechainConformation& chainIDSeqNumRotLibIDPairsForSidechainConformation)
{
	if(this == &chainIDSeqNumRotLibIDPairsForSidechainConformation)
		return *this;

	m_chainIDSeqNumRotLibIDPair = chainIDSeqNumRotLibIDPairsForSidechainConformation.m_chainIDSeqNumRotLibIDPair;

	return *this;
}