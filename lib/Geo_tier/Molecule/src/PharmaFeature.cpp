#include "PharmaFeature.h"
#include "EnclosingSphereOfSpheres.h"
#include "rg_Atom.h"
using namespace V::GeometryTier;



PharmaFeature::PharmaFeature()
{
    m_ID = -1;
    m_flagsOfPharmaFeatureType = 0;
    m_flagsOfChemFuncGroupType = 0;
    
}



PharmaFeature::PharmaFeature( const rg_INT& ID )
{
    m_ID = ID;
    m_flagsOfPharmaFeatureType = 0;
    m_flagsOfChemFuncGroupType = 0;    
}



PharmaFeature::PharmaFeature( const PharmaFeature& aPharmaFeature )
{
    m_ID                        = aPharmaFeature.m_ID ;

    m_flagsOfPharmaFeatureType  = aPharmaFeature.m_flagsOfPharmaFeatureType;
    m_flagsOfChemFuncGroupType  = aPharmaFeature.m_flagsOfChemFuncGroupType;

    m_listpAtoms.removeAll();
    m_listpAtoms.append( aPharmaFeature.m_listpAtoms );

    m_centerOfMass              = aPharmaFeature.m_centerOfMass;
    m_minEnclosingSphere        = aPharmaFeature.m_minEnclosingSphere;
}



PharmaFeature::~PharmaFeature()
{   
}



rg_INT PharmaFeature::getID() const
{
    return m_ID;
}



rg_dList<Atom*>* PharmaFeature::getAtoms()
{
    return &m_listpAtoms;
}



rg_Point3D PharmaFeature::getCenterOfMass() const
{
    return m_centerOfMass;
}



Sphere PharmaFeature::getMinEnclosingSphere() const
{
    return m_minEnclosingSphere;
}



rg_INT PharmaFeature::getPharmaFeatureType() const
{
    return m_flagsOfPharmaFeatureType;
}



rg_INT PharmaFeature::getChemFuncGroupType() const
{
    return m_flagsOfChemFuncGroupType;
}



rg_FLAG PharmaFeature::isPharmaFeatureType( const rg_INT& aPharmaFeatureType ) const
{
    if( m_flagsOfPharmaFeatureType&(1<<aPharmaFeatureType) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG PharmaFeature::isChemFuncGroupType( const rg_INT& aChemFuncGroupType ) const
{
    if( m_flagsOfChemFuncGroupType&(1<<aChemFuncGroupType) )
        return rg_TRUE;
    else
        return rg_FALSE;    
}



void PharmaFeature::setID( const rg_INT& ID )
{
    m_ID = ID;
}



void PharmaFeature::addPharmaFeatureType( const rg_INT& aPharmaFeatureType )
{
    m_flagsOfPharmaFeatureType = m_flagsOfPharmaFeatureType|(1<<aPharmaFeatureType);
}




void PharmaFeature::setFlagsOfPharmaFeatureType( const rg_INT& flagsOfPharmaFeatureType )
{
    m_flagsOfPharmaFeatureType = flagsOfPharmaFeatureType;
}


void PharmaFeature::mergePharmaFeatureType( const rg_INT& flagsOfPharmaFeatureType )
{
    m_flagsOfPharmaFeatureType = m_flagsOfPharmaFeatureType|flagsOfPharmaFeatureType;
}



void PharmaFeature::addChemFuncGroupType( const rg_INT& aChemFuncGroupType )
{
    m_flagsOfChemFuncGroupType = m_flagsOfChemFuncGroupType|(1<<aChemFuncGroupType);
}



void PharmaFeature::setFlagsOfChemFuncGroupType( const rg_INT& flagsOfChemFuncGroupType )
{
    m_flagsOfChemFuncGroupType = flagsOfChemFuncGroupType;
}



void PharmaFeature::mergeChemFuncGroupType( const rg_INT& flagsOfChemFuncGroupType )
{
    m_flagsOfChemFuncGroupType = m_flagsOfChemFuncGroupType|flagsOfChemFuncGroupType;
}



void PharmaFeature::addAtom( Atom* aFeatureAtom )
{
    m_listpAtoms.addTail( aFeatureAtom );
}



void PharmaFeature::setTypesAndAtoms( const rg_INT& aPharmaFeatureType, const rg_INT& aChemFuncGroupType, Atom* atomA, Atom* atomB /*= rg_NULL*/, Atom* atomC /*= rg_NULL*/, Atom* atomD /*= rg_NULL*/, Atom* atomE /*= rg_NULL*/, Atom* atomF /*= rg_NULL */ )
{
    m_flagsOfPharmaFeatureType = m_flagsOfPharmaFeatureType|(1<<aPharmaFeatureType);
    m_flagsOfChemFuncGroupType = m_flagsOfChemFuncGroupType|(1<<aChemFuncGroupType);

    m_listpAtoms.addTail( atomA );

	if( atomB != rg_NULL )
		m_listpAtoms.addTail( atomB );

	if( atomC != rg_NULL )
		m_listpAtoms.addTail( atomC );

	if( atomD != rg_NULL )
		m_listpAtoms.addTail( atomD );

	if( atomE != rg_NULL )
		m_listpAtoms.addTail( atomE );

	if( atomF != rg_NULL )
		m_listpAtoms.addTail( atomF );
}



void PharmaFeature::setTypesAndAtoms( const rg_INT& aPharmaFeatureType, const rg_INT& aChemFuncGroupType, list<Atom*>* atoms )
{
    m_flagsOfPharmaFeatureType = m_flagsOfPharmaFeatureType|(1<<aPharmaFeatureType);
    m_flagsOfChemFuncGroupType = m_flagsOfChemFuncGroupType|(1<<aChemFuncGroupType);

    list<Atom*>::iterator i_atoms = atoms->begin();

    while( i_atoms != atoms->end() ) {
        Atom* currAtom = *i_atoms;
        
        m_listpAtoms.addTail( currAtom );
        i_atoms++;
    }
}



PharmaFeature& PharmaFeature::operator=( const PharmaFeature& aPharmaFeature )
{
    if( this == &aPharmaFeature)
        return *this;

    m_ID                        = aPharmaFeature.m_ID ;

    m_flagsOfPharmaFeatureType  = aPharmaFeature.m_flagsOfPharmaFeatureType;
    m_flagsOfChemFuncGroupType  = aPharmaFeature.m_flagsOfChemFuncGroupType;

    m_listpAtoms.removeAll();
    m_listpAtoms.append( aPharmaFeature.m_listpAtoms );

    m_centerOfMass              = aPharmaFeature.m_centerOfMass;
    m_minEnclosingSphere        = aPharmaFeature.m_minEnclosingSphere;

    return *this;    
}



void PharmaFeature::setCentersOfPharmaFeature()
{
    computeAndSetCenterOfMass();
    computeAndSetMinEnclosingSphere();
}



void PharmaFeature::computeAndSetCenterOfMass()
{
    rg_Point3D sumOfAllAtomCenters;
    m_listpAtoms.reset4Loop();
    while( m_listpAtoms.setNext4Loop() ) {
        sumOfAllAtomCenters += m_listpAtoms.getEntity()->getAtomBall().getCenter();
    }
    m_centerOfMass = sumOfAllAtomCenters/m_listpAtoms.getSize();
}



void PharmaFeature::computeAndSetMinEnclosingSphere()
{
    EnclosingSphereOfSpheres enclosingSphereOfMolecule;

    list<Sphere> listOfAtomBalls;
    m_listpAtoms.reset4Loop();
    while( m_listpAtoms.setNext4Loop() ) {
        Sphere atomBall = m_listpAtoms.getEntity()->getAtomBall();
        listOfAtomBalls.push_back( atomBall );
    }
        
    enclosingSphereOfMolecule.setSpheres( listOfAtomBalls );
    enclosingSphereOfMolecule.computeEnclosingSphere();
    
    m_minEnclosingSphere = enclosingSphereOfMolecule;  
}
