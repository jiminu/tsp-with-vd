#include "FunctionsForMolecule.h"
#include "SidechainDihedralAngleSet.h"
using namespace V::GeometryTier;


SidechainDihedralAngleSet::SidechainDihedralAngleSet()
{
	setInitialValues();
}

SidechainDihedralAngleSet::SidechainDihedralAngleSet(const SidechainDihedralAngleSet& sidechainDihedralAngleSet)
{
	setInitialValues();

	setDihedralAngles(sidechainDihedralAngleSet.m_dihedralAngles, sidechainDihedralAngleSet.m_numDihedralAngles);
	m_chainIDSeqNum = sidechainDihedralAngleSet.m_chainIDSeqNum;
}

SidechainDihedralAngleSet::SidechainDihedralAngleSet(Rotamer& rotamer)
{
	setInitialValues();

	rg_INT numDihedralAngles = rotamer.getNumDihedralAngles();
	rg_REAL* dihedralAngles = rotamer.getDihedralAngles();
	setDihedralAngles(dihedralAngles, numDihedralAngles);
	Residue* residue = rotamer.getResidue();
	m_chainIDSeqNum = getChainIDSeqNumOfResidue( residue );
}

SidechainDihedralAngleSet::SidechainDihedralAngleSet(Rotamer* rotamer)
{
    if(rotamer != rg_NULL)
    {        
	    setInitialValues();

	    rg_INT numDihedralAngles = rotamer->getNumDihedralAngles();
	    rg_REAL* dihedralAngles = rotamer->getDihedralAngles();
	    setDihedralAngles(dihedralAngles, numDihedralAngles);
	    Residue* residue = rotamer->getResidue();
	    m_chainIDSeqNum = getChainIDSeqNumOfResidue( residue );
    }
}

SidechainDihedralAngleSet::~SidechainDihedralAngleSet()
{
	destroyDihedralAngles() ;
	m_chainIDSeqNum   = -1;
}

void     SidechainDihedralAngleSet::setDihedralAngles( rg_REAL* dihedralAngles, const rg_INT& numDihedralAngles )
{
	destroyDihedralAngles();
	
	m_numDihedralAngles  = numDihedralAngles;

	if(m_numDihedralAngles > 0)
	{
		m_dihedralAngles = new rg_REAL[ m_numDihedralAngles ];
		rg_INDEX i;
		for (i = 0;i < m_numDihedralAngles;i++)
		{
			m_dihedralAngles[ i ] = dihedralAngles[ i ];
		}
	}
}

SidechainDihedralAngleSet& SidechainDihedralAngleSet::operator=(const SidechainDihedralAngleSet& sidechainDihedralAngleSet)
{
	if(this == &sidechainDihedralAngleSet)
		return *this;

	setDihedralAngles(sidechainDihedralAngleSet.m_dihedralAngles, sidechainDihedralAngleSet.m_numDihedralAngles);
	m_chainIDSeqNum = sidechainDihedralAngleSet.m_chainIDSeqNum;

	return *this;
}

void SidechainDihedralAngleSet::destroyDihedralAngles()
{
	if(m_numDihedralAngles > 0 && m_dihedralAngles != rg_NULL)
	{
		delete [] m_dihedralAngles;		
		m_dihedralAngles = rg_NULL;
		m_numDihedralAngles = 0;
	}
}