#ifndef _CHAINIDWITHSEQUENCENUMBERSOFMISSINGRESIDUES_
#define _CHAINIDWITHSEQUENCENUMBERSOFMISSINGRESIDUES_

#include "rg_dList.h"
#include <string>
using namespace std;

class ChainIDWithSequenceNumbersOfMissingResidues
{
private:
	string           m_chainIDFromInputPDBFile;
	rg_dList<rg_INT> m_seqNumbersOfMissingResidues;

	inline void set(const string& chainIDFromInputPDBFile, const rg_dList<rg_INT>& seqNumbersOfMissingResidues)
	{
		m_chainIDFromInputPDBFile = chainIDFromInputPDBFile;
		m_seqNumbersOfMissingResidues = seqNumbersOfMissingResidues;
	}

public:
	ChainIDWithSequenceNumbersOfMissingResidues();
	ChainIDWithSequenceNumbersOfMissingResidues(const ChainIDWithSequenceNumbersOfMissingResidues& chainIDWithSequenceNumbersOfMissingResidues);
	ChainIDWithSequenceNumbersOfMissingResidues(string& chainIDFromInputPDBFile, rg_dList<rg_INT>& seqNumbersOfMissingResidues);
	~ChainIDWithSequenceNumbersOfMissingResidues();
	
	inline string getChainIDFromInputPDBFile() const { return m_chainIDFromInputPDBFile; }
	inline rg_dList<rg_INT>* getSeqNumbersOfMissingResidues() { return &m_seqNumbersOfMissingResidues; }
	inline void   getSeqNumbersOfMissingResidues( rg_dList<rg_INT>& seqNumbersOfMissingResidues ) { seqNumbersOfMissingResidues = m_seqNumbersOfMissingResidues; }

	ChainIDWithSequenceNumbersOfMissingResidues& operator=(const ChainIDWithSequenceNumbersOfMissingResidues& chainIDWithSequenceNumbersOfMissingResidues);
};

#endif
