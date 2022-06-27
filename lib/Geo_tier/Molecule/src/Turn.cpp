#include "Turn.h"
#include "Residue.h"
#include "Chain.h"
//using namespace V::GeometryTier;


V::GeometryTier::Turn::Turn()
{
}



V::GeometryTier::Turn::Turn( const rg_INT& serial )
{
    m_serial = serial;
}



V::GeometryTier::Turn::Turn( const rg_INT& serial, const string& turnID )
{
    m_serial = serial;
    m_turnID = turnID;
}



V::GeometryTier::Turn::Turn( const rg_INT& serial, const string& turnID, rg_dList<V::GeometryTier::Residue*>* residues, const string& comments )
{
    m_serial     = serial;
    m_turnID    = turnID;

    m_residues.append( *residues );
    m_comments = comments;
}



V::GeometryTier::Turn::Turn( const V::GeometryTier::Turn& aTurn)
{
    m_serial         = aTurn.m_serial;
    m_turnID         = aTurn.m_turnID;
    m_residues       = aTurn.m_residues;
    m_comments       = aTurn.m_comments;
}



V::GeometryTier::Turn::~Turn()
{
}


//  GET FUNCTION

V::GeometryTier::Chain* V::GeometryTier::Turn::getChain()
{
    return m_residues.getFirstEntity()->getChain();
}



rg_BOOL V::GeometryTier::Turn::isResidueInTurn(V::GeometryTier::Residue* aResidue ) const
{
    rg_BOOL isInTurn = rg_FALSE;

    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        if( aResidue == m_residues.getEntity() ) {
            isInTurn = rg_TRUE;
            break;
        }
    }

    return isInTurn;
}



//  SET FUNCTION

//  OPERATOR OVERLOADING
V::GeometryTier::Turn& V::GeometryTier::Turn::operator =( const V::GeometryTier::Turn& aTurn )
{
    if( this == &aTurn )
        return *this;

    m_serial         = aTurn.m_serial;
    m_turnID         = aTurn.m_turnID;
    m_residues = aTurn.m_residues;
    m_comments       = aTurn.m_comments;

    return *this;
}
