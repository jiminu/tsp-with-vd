#ifndef _ATOM_H
#define _ATOM_H

#include "rg_Const.h"
#include "ConstForMolecule.h"
#include "Sphere.h"
#include "ChemicalPropertiesOfAtom.h"
#include "rg_dList.h"


namespace V {
namespace GeometryTier {


class ChemicalBond;
class Residue;
class Chain;

class Atom
{
private:
    rg_INT                    m_ID;
    AtomCode                  m_atomCode;
    Sphere                    m_atomBall;
    ChemicalPropertiesOfAtom  m_chemProperties;

    rg_dList<ChemicalBond*>   m_chemicalBonds;
    Residue*                  m_residue;

//////////////////////////////////////////////////////
    rg_INT                    m_serialFromInputFile;
    string                    m_atomNameFromInputFile;
//////////////////////////////////////////////////////

public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Atom();
    Atom( const rg_INT& ID );
    Atom( const Atom& atom );
    ~Atom();

    //  GET FUNCTION
    rg_INT                      getID() const;
    AtomCode                    getAtomCode() const;
    Sphere                      getAtomBall() const;
    Sphere*                     getpAtomBall();
    ChemicalPropertiesOfAtom    getChemicalProperties() const;
    ChemicalPropertiesOfAtom*   getpChemicalProperties();
    rg_dList<ChemicalBond*>*    getListChemicalBond();
    Residue*                    getResidue();
    Chain*                      getChain();
    rg_INT                      getSerialFromInputFile() const;
    string                      getAtomNameFromInputFile() const;
    string                      getAtomNameInPDBFormat() const;


    //  SET FUNCTION
    void    setID( const rg_INT& ID );
    void    setAtomCode( const AtomCode& atomCode );
    void    setAtomBall( const Sphere& atomBall );
    void    setCenterOfAtomBall( const rg_Point3D& center );
    void    setChemicalProperties( const ChemicalPropertiesOfAtom& chemicalProperties );
    void    addChemicalBond( ChemicalBond* aChemicalBond );
    void    setResidue( Residue* residue );
    void    setSerialFromInputFile( const rg_INT& serial );
    void    setAtomNameFromInputFile( const string& atomName );


    //  OPERATOR OVERLOADING
    Atom& operator =(const Atom& atom);

    rg_BOOL hasChemicalBondWith(Atom* atom);   
    rg_REAL compute_Lennard_Jones_potential(const Atom & atom);

    static bool AtomIDLess(Atom* atom1, Atom* atom2) {
        return atom1->getSerialFromInputFile() < atom2->getSerialFromInputFile();
    }
};


inline  rg_INT                                          Atom::getID() const                     { return m_ID; }
inline  V::GeometryTier::AtomCode                    Atom::getAtomCode() const               { return m_atomCode; }
inline  Sphere                                          Atom::getAtomBall() const               { return m_atomBall; }
inline  Sphere*                                         Atom::getpAtomBall()                    { return &m_atomBall; }
inline  V::GeometryTier::ChemicalPropertiesOfAtom    Atom::getChemicalProperties() const     { return m_chemProperties; }
inline  V::GeometryTier::ChemicalPropertiesOfAtom*   Atom::getpChemicalProperties()          { return &m_chemProperties; }
inline  rg_dList<V::GeometryTier::ChemicalBond*>*    Atom::getListChemicalBond()             { return &m_chemicalBonds; }
inline  V::GeometryTier::Residue*                    Atom::getResidue()                      { return m_residue; }
inline  rg_INT                                          Atom::getSerialFromInputFile() const    { return m_serialFromInputFile; }
inline  string                                          Atom::getAtomNameFromInputFile() const  { return m_atomNameFromInputFile; }

inline  void    Atom::setID(const rg_INT& ID)                                                       { m_ID = ID; }
inline  void    Atom::setAtomCode(const AtomCode& atomCode)                                         { m_atomCode = atomCode; }
inline  void    Atom::setAtomBall(const Sphere& atomBall)                                           { m_atomBall = atomBall; }
inline  void    Atom::setCenterOfAtomBall(const rg_Point3D& center)                                 { m_atomBall.setCenter(center); }
inline  void    Atom::setChemicalProperties(const ChemicalPropertiesOfAtom& chemicalProperties)     { m_chemProperties = chemicalProperties; }
inline  void    Atom::addChemicalBond(ChemicalBond* aChemicalBond)                                  { m_chemicalBonds.addTail(aChemicalBond); }
inline  void    Atom::setResidue(Residue* residue)                                                  { m_residue = residue; }
inline  void    Atom::setSerialFromInputFile(const rg_INT& serial)                                  { m_serialFromInputFile = serial; }
inline  void    Atom::setAtomNameFromInputFile(const string& atomName)                              { m_atomNameFromInputFile = atomName; }



} // namespace GeometryTier
} // namespace V


#endif

