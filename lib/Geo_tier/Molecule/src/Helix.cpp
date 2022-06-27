#include "Helix.h"
#include "Chain.h"


V::GeometryTier::Helix::Helix()
{
}



V::GeometryTier::Helix::Helix( const rg_INT& serial, const string& helixID, const V::GeometryTier::HelixClass& helixClass )
{
    m_serial     = serial;
    m_helixID    = helixID;
    m_helixClass = helixClass;
}



V::GeometryTier::Helix::Helix( const rg_INT& serial, const string& helixID, const rg_INT& helixClass )
{
    m_serial     = serial;
    m_helixID    = helixID;

    if( helixClass > 0 && helixClass < numberOfHelixClass ) {
        m_helixClass = (HelixClass)helixClass;
    }
    else {
        m_helixClass = UNK_HELIX_CLASS;
    }
}



V::GeometryTier::Helix::Helix( 
                                const rg_INT& serial, 
                                const string& helixID, 
                                const rg_INT& helixClass, 
                                rg_dList<V::GeometryTier::Residue*>* residues, 
                                const string& comments )
{
    m_serial     = serial;
    m_helixID    = helixID;

    if( helixClass > 0 && helixClass < numberOfHelixClass ) {
        m_helixClass = (HelixClass)helixClass;
    }
    else {
        m_helixClass = UNK_HELIX_CLASS;
    }

    m_residues.append( *residues );
    m_comments = comments;
}



V::GeometryTier::Helix::Helix( const V::GeometryTier::Helix& aHelix )
{
    m_serial         = aHelix.m_serial;
    m_helixID        = aHelix.m_helixID;
    m_residues       = aHelix.m_residues;
    m_helixClass     = aHelix.m_helixClass;
    m_comments       = aHelix.m_comments;
}



V::GeometryTier::Helix::~Helix()
{
}



//  GET FUNCTION

V::GeometryTier::Chain* V::GeometryTier::Helix::getChain()
{
    return m_residues.getFirstEntity()->getChain();
}



rg_BOOL V::GeometryTier::Helix::isResidueInHelix(V::GeometryTier::Residue* aResidue ) const
{
    rg_BOOL isInHelix = rg_FALSE;

    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        if( aResidue == m_residues.getEntity() ) {
            isInHelix = rg_TRUE;
            break;
        }
    }

    return isInHelix;
}



//  SET FUNCTION

//  OPERATOR OVERLOADING
V::GeometryTier::Helix& V::GeometryTier::Helix::operator =( const V::GeometryTier::Helix& aHelix )
{
    if( this == &aHelix )
        return *this;
    
    m_serial         = aHelix.m_serial;
    m_helixID        = aHelix.m_helixID;
    m_residues = aHelix.m_residues;
    m_helixClass     = aHelix.m_helixClass;
    m_comments       = aHelix.m_comments;
    
    return *this;
}
