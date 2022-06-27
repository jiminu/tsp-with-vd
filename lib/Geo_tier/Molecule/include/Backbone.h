#ifndef _BACKBONE_H_
#define _BACKBONE_H_


#include "BackboneConformation.h"
#include "ConstForMolecule.h"

//using namespace V::GeometryTier;

class V::GeometryTier::Atom;



// This class represents an object for representing the atoms on the protein backbone.
// The atoms are stored in the order of the protein sequence.
// In the case of the multiple chains, we concatenate the backbone atoms of each chain with each other in the order of chain ID.

class Backbone
{
private:
	// N, Ca, C, O atoms of each residue
    V::GeometryTier::Atom**    m_backboneAtoms;
	rg_INT    m_numAtoms     ;
	inline void setInitialValues()
	{
		m_backboneAtoms = rg_NULL;
		m_numAtoms      = 0;
	}

public:
	Backbone();
	Backbone(V::GeometryTier::Atom** backboneAtoms, const rg_INT& numAtoms );
	Backbone( Backbone& backbone );
	~Backbone();

	inline rg_INT  getBackboneAtoms(V::GeometryTier::Atom**& backboneAtoms )
	{
		backboneAtoms = m_backboneAtoms;
		return m_numAtoms;
	}

	rg_INT getBackboneAtoms( rg_dList<V::GeometryTier::Atom*>& backboneAtoms );

	inline V::GeometryTier::Atom**  getBackboneAtoms( rg_INT& numAtoms      )
	{
		numAtoms = m_numAtoms;
		return m_backboneAtoms;
	}

	void    getBackboneConformation(BackboneConformation& backboneConformation);
	void    getSequence(rg_dList<V::GeometryTier::ResidueCode>& sequence);

    void    setBackboneAtoms(V::GeometryTier::Atom**& backboneAtoms, const rg_INT& numAtoms );
    void    setBackbonConformation(BackboneConformation& backboneConformation);

    //void    perturbAlphaCarbonByRotation(Atom* targetAlphaCarbon, const rg_REAL& angle);

	Backbone& operator=( Backbone& backbone );

private:
	void destroy();
};

#endif