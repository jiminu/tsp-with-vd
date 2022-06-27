#ifndef _PHARMAFEATURE_H
#define _PHARMAFEATURE_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "Sphere.h"
#include "ConstForMolecule.h"

#include <list>
using namespace std;



namespace V {
namespace GeometryTier {



class Atom;

class PharmaFeature
{
private:
    rg_INT              m_ID;
    
    rg_INT              m_flagsOfPharmaFeatureType; // BIT OPERATORS : CONTAINS PharmaFeatureType
    rg_INT              m_flagsOfChemFuncGroupType; // BIT OPERATORS : CONTAINS ChemFuncGroupType

    rg_dList<Atom*>     m_listpAtoms;

    rg_Point3D          m_centerOfMass;
    Sphere              m_minEnclosingSphere;
    
    
public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    PharmaFeature();
    PharmaFeature( const rg_INT& ID );
    PharmaFeature( const PharmaFeature& aPharmaFeature );
    ~PharmaFeature();


    //  GET FUNCTION
    rg_INT              getID() const;
    rg_dList<Atom*>*    getAtoms();

    rg_Point3D          getCenterOfMass() const;
    Sphere              getMinEnclosingSphere() const;

    rg_INT              getPharmaFeatureType() const;
    rg_INT              getChemFuncGroupType() const;


    rg_FLAG             isPharmaFeatureType( const rg_INT& aPharmaFeatureType ) const;
    rg_FLAG             isChemFuncGroupType( const rg_INT& aChemFuncGroupType ) const;



    //  SET FUNCTION
    void                setID( const rg_INT& ID );

    void                addPharmaFeatureType( const rg_INT& aPharmaFeatureType );
    void                setFlagsOfPharmaFeatureType( const rg_INT& flagsOfPharmaFeatureType );

    void                addChemFuncGroupType( const rg_INT& aChemFuncGroupType );
    void                setFlagsOfChemFuncGroupType( const rg_INT& flagsOfChemFuncGroupType );

    void                addAtom( Atom* aFeatureAtom );

    void                setTypesAndAtoms( const rg_INT& aPharmaFeatureType, const rg_INT& aChemFuncGroupType, Atom* atomA, Atom* atomB = rg_NULL, Atom* atomC = rg_NULL, Atom* atomD = rg_NULL, Atom* atomE = rg_NULL, Atom* atomF = rg_NULL );
    void                setTypesAndAtoms( const rg_INT& aPharmaFeatureType, const rg_INT& aChemFuncGroupType, list<Atom*>* atoms );
    void                setCentersOfPharmaFeature();

    void                mergePharmaFeatureType( const rg_INT& flagsOfPharmaFeatureType );
    void                mergeChemFuncGroupType( const rg_INT& flagsOfChemFuncGroupType );


    //  OPERATOR OVERLOADING
    PharmaFeature& operator = ( const PharmaFeature& aPharmaFeature );


private:

    void                computeAndSetCenterOfMass();
    void                computeAndSetMinEnclosingSphere();
};



} // namespace GeometryTier
} // namespace V


#endif

