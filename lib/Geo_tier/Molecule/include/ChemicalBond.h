#ifndef _CHEMICALBOND_H
#define _CHEMICALBOND_H

#include "rg_Const.h"
#include "ConstForMolecule.h"



namespace V {
namespace GeometryTier {



class Atom;

class ChemicalBond
{
private:
    rg_INT          m_ID;
    Atom*           m_firstAtom;
    Atom*           m_secondAtom;
    BondType        m_typeOfBond;


//////////////////////////////////////////////////////
    rg_INT          m_serialFromInputFile;
//////////////////////////////////////////////////////


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    ChemicalBond();
    ChemicalBond( const rg_INT& ID );
    ChemicalBond( const rg_INT& ID, Atom* firstAtom, Atom* secondAtom );
    ChemicalBond( const rg_INT& ID, Atom* firstAtom, Atom* secondAtom, const BondType& typeOfBond );
    ChemicalBond( const ChemicalBond& chemicalBond );
    ~ChemicalBond();

    //  GET FUNCTION
    rg_INT                   getID() const;
    inline Atom*             getFirstAtom() { return m_firstAtom; };
    inline Atom*             getSecondAtom() { return m_secondAtom; };
    inline void              getAtoms( Atom*& firstAtom, Atom*& secondAtom ) { firstAtom = m_firstAtom;  secondAtom = m_secondAtom; };
    BondType                 getTypeOfBond() const;
    Atom*                    getBondedAtom( Atom* baseAtom );

    rg_FLAG	                 isOnBackBone() const;

    //  SET FUNCTION
    void                     setChemicalBond( const rg_INT& ID, Atom* firstAtom, Atom* secondAtom, const BondType& typeOfBond );
    void                     setID( const rg_INT& ID );
    void                     setFirstAtom( Atom* firstAtom );
    void                     setSecondAtom( Atom* secondAtom );
    void                     setAtoms( Atom* firstAtom, Atom* secondAtom );
    void                     setTypeOfBond( const BondType& typeOfBond );

//////////////////////////////////////////////////////
    rg_INT  getSerialFromInputFile() const;
    void    setSerialFromInputFile( const rg_INT& serial );
//////////////////////////////////////////////////////

    //  OPERATOR OVERLOADING
    ChemicalBond& operator =(const ChemicalBond& chemicalBond);
    rg_FLAG       operator==(const ChemicalBond& chemicalBond);

};



} // namespace GeometryTier
} // namespace V


#endif

