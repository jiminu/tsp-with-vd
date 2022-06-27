#include "Backbone.h"
#include "Residue.h"
#include "rg_Atom.h"
#include "rg_TMatrix3D.h"

typedef rg_Point3D Vector3D;


using namespace V::GeometryTier;




Backbone::Backbone()
{
	setInitialValues();
}

Backbone::Backbone(V::GeometryTier::Atom** backboneAtoms, const rg_INT& numAtoms )
{
	setInitialValues();
	setBackboneAtoms(backboneAtoms, numAtoms);
}

Backbone::Backbone( Backbone& backbone )
{
	setInitialValues();
	setBackboneAtoms(backbone.m_backboneAtoms, backbone.m_numAtoms);
}

Backbone::~Backbone()
{
	destroy();
}

rg_INT Backbone::getBackboneAtoms( rg_dList<V::GeometryTier::Atom*>& backboneAtoms )
{
	rg_INDEX i;
	for (i = 0;i < m_numAtoms;i++)
	{
		backboneAtoms.add( m_backboneAtoms[ i ] );
	}
	return backboneAtoms.getSize();
}

//rg_INT  Backbone::getBackboneAtoms( Atom**& backboneAtoms )
//{
//}

//Atom**  Backbone::getBackboneAtoms( rg_INT& numAtoms      )
//{
//}

void    Backbone::getBackboneConformation(BackboneConformation& backboneConformation)
{
	backboneConformation.set(m_backboneAtoms, m_numAtoms);
}

void    Backbone::getSequence(rg_dList<V::GeometryTier::ResidueCode>& sequence)
{
	rg_INDEX i;
	for(i = 0;i < m_numAtoms;i++)
	{
		if(m_backboneAtoms[ i ]->getAtomCode() == C_ATOM &&
			m_backboneAtoms[ i ]->getpChemicalProperties()->getRemoteIndicator() == ALPHA_REMOTE)
		{
			ResidueCode code = m_backboneAtoms[ i ]->getResidue()->getResidueCode();
			sequence.add( code );
		}
	}
}

void    Backbone::setBackboneAtoms(V::GeometryTier::Atom**& backboneAtoms, const rg_INT& numAtoms )
{
	destroy();

	m_numAtoms = numAtoms;
	m_backboneAtoms = new V::GeometryTier::Atom*[ m_numAtoms ];
	
	rg_INDEX i;
	for (i = 0;i < m_numAtoms;i++)
	{
		m_backboneAtoms[ i ] = backboneAtoms[ i ];
	}
}

void    Backbone::setBackbonConformation(BackboneConformation& backboneConformation)
{

}

//void    Backbone::perturbAlphaCarbonByRotation(Atom* targetAlphaCarbon, const rg_REAL& angle)
//{
//    Residue*   targetResidue         = targetAlphaCarbon->getResidue();
//    Atom*      nitrogenInAminoGroup  = targetResidue->getNitrogenInAminoGroupOfAminoResidue();
//    Atom*      carbonInCarboxylGroup = targetResidue->getCarbonInCarboxylGroupOfAminoResidue();
//
//    rg_Point3D centerOfAlphaCarbon   = targetAlphaCarbon->getpAtomBall()->getCenter();
//    rg_Point3D centerOfNitrogen      = nitrogenInAminoGroup->getpAtomBall()->getCenter();
//    rg_Point3D centerOfCarbon        = carbonInCarboxylGroup->getpAtomBall()->getCenter();
//
//    Vector3D rotationAxis = centerOfCarbon - centerOfNitrogen;
//    rotationAxis.normalize();
//    rg_TMatrix3D rotMat;
//    rotMat.rotateArbitraryAxis(rotationAxis, rg_PI * angle / 180.);
//
//    targetAlphaCarbon->getpAtomBall()->setCenter( rotMat * centerOfAlphaCarbon );
//}

Backbone&  Backbone::operator=(Backbone& backbone)
{
	if(this == &backbone)
		return *this;

	setBackboneAtoms(backbone.m_backboneAtoms, backbone.m_numAtoms);
	
	return *this;
}

void    Backbone::destroy()
{
	if (m_numAtoms > 0 && m_backboneAtoms != rg_NULL)
	{
		delete [] m_backboneAtoms;
		m_backboneAtoms = rg_NULL;
		m_numAtoms = 0;
	}
}

