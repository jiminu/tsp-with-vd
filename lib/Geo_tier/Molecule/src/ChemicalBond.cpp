#include "ChemicalBond.h"
#include "rg_Atom.h"
using namespace V::GeometryTier;



ChemicalBond::ChemicalBond()
: m_ID(0), m_firstAtom(rg_NULL), m_secondAtom(rg_NULL), m_typeOfBond(UNK_BOND)
{
    m_serialFromInputFile = 0;
}



ChemicalBond::ChemicalBond( const rg_INT& ID )
: m_ID(ID), m_firstAtom(rg_NULL), m_secondAtom(rg_NULL), m_typeOfBond(UNK_BOND)
{
    m_serialFromInputFile = 0;
}



ChemicalBond::ChemicalBond( const rg_INT& ID, Atom* firstAtom, Atom* secondAtom )
: m_ID(ID), m_firstAtom(firstAtom), m_secondAtom(secondAtom), m_typeOfBond(UNK_BOND)
{
    m_serialFromInputFile = 0;
}



ChemicalBond::ChemicalBond( const rg_INT& ID, Atom* firstAtom, Atom* secondAtom, const BondType& typeOfBond )
: m_ID(ID), m_firstAtom(firstAtom), m_secondAtom(secondAtom), m_typeOfBond(typeOfBond)
{
    m_serialFromInputFile = 0;
}



ChemicalBond::ChemicalBond( const ChemicalBond& chemicalBond )
: m_ID(chemicalBond.m_ID), m_firstAtom(chemicalBond.m_firstAtom), m_secondAtom(chemicalBond.m_secondAtom), m_typeOfBond(chemicalBond.m_typeOfBond)
{
    m_serialFromInputFile = chemicalBond.m_serialFromInputFile;
}



ChemicalBond::~ChemicalBond()
{
}



rg_INT ChemicalBond::getID() const
{
    return m_ID;
}



// Atom* ChemicalBond::getFirstAtom()
// {
//     return m_firstAtom;
// }
// 
// 
// 
// Atom* ChemicalBond::getSecondAtom()
// {
//     return m_secondAtom;
// }



// void ChemicalBond::getAtoms( Atom*& firstAtom, Atom*& secondAtom )
// {
//     firstAtom  = m_firstAtom;
//     secondAtom = m_secondAtom;
// }



BondType ChemicalBond::getTypeOfBond() const
{
    return m_typeOfBond;
}



Atom* ChemicalBond::getBondedAtom( Atom* baseAtom )
{
    if( m_firstAtom == baseAtom ) {
        return m_secondAtom;
    }
    else if( m_secondAtom == baseAtom ) {
        return m_firstAtom;
    }
    else {
        return rg_NULL;
    }
}



rg_FLAG ChemicalBond::isOnBackBone() const
{
	if(  m_firstAtom->getpChemicalProperties()->isOnBackBone()  == rg_TRUE && 
        m_secondAtom->getpChemicalProperties()->isOnBackBone() == rg_TRUE    ) {
		return rg_TRUE;
    }
    else {
		return rg_FALSE;   
    }
}



void ChemicalBond::setChemicalBond( const rg_INT& ID, Atom* firstAtom, Atom* secondAtom, const BondType& typeOfBond )
{
    m_ID         = ID;
    m_firstAtom  = firstAtom;    
    m_secondAtom = secondAtom;
    m_typeOfBond  = typeOfBond;    
}



void ChemicalBond::setID( const rg_INT& ID )
{
    m_ID = ID;    
}



void ChemicalBond::setFirstAtom( Atom* firstAtom )
{
    m_firstAtom = firstAtom;
}



void ChemicalBond::setSecondAtom( Atom* secondAtom )
{
    m_secondAtom = secondAtom;
}



void ChemicalBond::setAtoms( Atom* firstAtom, Atom* secondAtom )
{
    m_firstAtom  = firstAtom;
    m_secondAtom = secondAtom;
}



void ChemicalBond::setTypeOfBond( const BondType& typeOfBond )
{
    m_typeOfBond = typeOfBond;    
}



//////////////////////////////////////////////////////
rg_INT  ChemicalBond::getSerialFromInputFile() const
{
    return m_serialFromInputFile;
}



void    ChemicalBond::setSerialFromInputFile( const rg_INT& serial )
{
    m_serialFromInputFile = serial;
}
//////////////////////////////////////////////////////



ChemicalBond& ChemicalBond::operator=( const ChemicalBond& chemicalBond )
{
    if( this == &chemicalBond)
        return *this;

    m_ID         = chemicalBond.m_ID;
    m_firstAtom  = chemicalBond.m_firstAtom;
    m_secondAtom = chemicalBond.m_secondAtom;
    m_typeOfBond = chemicalBond.m_typeOfBond;
    m_serialFromInputFile = chemicalBond.m_serialFromInputFile;

    return *this;
}



rg_FLAG ChemicalBond::operator==( const ChemicalBond& chemicalBond )
{
    if ( ( m_firstAtom==chemicalBond.m_firstAtom  && m_secondAtom==chemicalBond.m_secondAtom) || 
         ( m_firstAtom==chemicalBond.m_secondAtom && m_secondAtom==chemicalBond.m_firstAtom )   ) {
        return rg_TRUE;
    }   
    else {
        return rg_FALSE;
    }
}
