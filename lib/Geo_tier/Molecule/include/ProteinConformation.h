#ifndef _PROTEINCONFORMATION_H_
#define _PROTEINCONFORMATION_H_

#include "BackboneConformation.h"
#include "SidechainConformation.h"

// This class represents the conformation for a protein structure 
// which consists of backbone conformation and sidechain conformation.

//enum ProteinConformationType {UNKNOWN_TYPE, FIXED_BACKBONE, FLEXIBLE_BACKBONE};

class ProteinConformation
{
private:
	BackboneConformation     m_backboneConformation ;
	SidechainConformation    m_sidechainConformation;
	//ProteinConformationType  m_type                 ;
	inline void setInitialValues()
	{
		//m_type = UNKNOWN_TYPE;
	}

public:
	ProteinConformation();
	ProteinConformation(ProteinConformation& proteinConformation);
	//ProteinConformation( const ProteinConformationType& type );
	~ProteinConformation();

	inline BackboneConformation * getBackboneConformation() { return &m_backboneConformation; }
	inline SidechainConformation* getSidechainConformation() { return &m_sidechainConformation; }

    void getConformationOfResidue(const rg_INT& chainIDSeqNum, 
                                  BackboneDihedralAnglePair& backboneDihedralAnglePair, 
                                  SidechainDihedralAngleSet& sideChainDihedralAngleSet);
    void getConformationOfResidue(V::GeometryTier::Residue* residue,
                                  BackboneDihedralAnglePair& backboneDihedralAnglePair, 
                                  SidechainDihedralAngleSet& sideChainDihedralAngleSet);

    inline void set(const BackboneConformation& backboneConformation, 
                    const SidechainConformation& sidechainConformation 
                    /*const ProteinConformationType& type = FIXED_BACKBONE*/)
    {
        m_backboneConformation = backboneConformation;
        m_sidechainConformation = sidechainConformation;
        //m_type = type;
    }

	//inline void setType(const ProteinConformationType& type) {m_type = type;}
	
	//rg_REAL evaluateEnergyFunction();
	//void    evaluateDerivativeOfEnergyFunction();

	ProteinConformation& operator= (const ProteinConformation& proteinConformation);
};

#endif
