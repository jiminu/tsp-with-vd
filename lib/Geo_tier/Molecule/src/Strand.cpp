#include "Strand.h"
#include "Sheet.h"
#include "Chain.h"
//using namespace V::GeometryTier;


V::GeometryTier::Strand::Strand()
{
    m_serial = 0;
    m_sheet  = rg_NULL;
}



V::GeometryTier::Strand::Strand( const rg_INT& strandSerial )
{
    m_serial = strandSerial;
    m_sheet  = rg_NULL;
}



V::GeometryTier::Strand::Strand( const rg_INT& strandSerial, V::GeometryTier::Sheet* aSheet )
{
    m_serial = strandSerial;
    m_sheet  = aSheet;
}



V::GeometryTier::Strand::Strand( const V::GeometryTier::Strand& aStrand)
{
    m_serial   = aStrand.m_serial;
    m_sheet    = aStrand.m_sheet;
    m_residues = aStrand.m_residues;
}



V::GeometryTier::Strand::~Strand()
{
}



//  GET FUNCTION

V::GeometryTier::Chain*  V::GeometryTier::Strand::getChain()
{
    return m_residues.getFirstEntity()->getChain();
}



rg_BOOL     V::GeometryTier::Strand::isResidueInStrand(V::GeometryTier::Residue* aResidue ) const
{
    rg_BOOL isInStrand = rg_FALSE;

    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        if( aResidue == m_residues.getEntity() ) {
            isInStrand = rg_TRUE;
            break;
        }
    }

    return isInStrand;    
}



//  SET FUNCTION



//  OPERATOR OVERLOADING
V::GeometryTier::Strand& V::GeometryTier::Strand::operator =( const V::GeometryTier::Strand& aStrand )
{
    if( this == &aStrand )
        return *this;

    m_serial   = aStrand.m_serial;
    m_sheet    = aStrand.m_sheet;
    m_residues = aStrand.m_residues;

    return *this;
}
