#include "Sheet.h"
#include "Chain.h"
//using namespace V::GeometryTier;


V::GeometryTier::Sheet::Sheet()
: m_numOfStrands(0)
, m_strands(rg_NULL)
, m_sensesOfStrands(rg_NULL)
, m_HBondsBetweenStrands(rg_NULL)
{
}



V::GeometryTier::Sheet::Sheet( const string& sheetID )
: m_sheetID(sheetID)
, m_numOfStrands(0)
, m_strands(rg_NULL)
, m_sensesOfStrands(rg_NULL)
, m_HBondsBetweenStrands(rg_NULL)
{
}



V::GeometryTier::Sheet::Sheet( const string& sheetID, const rg_INT& numOfStrands )
: m_sheetID(sheetID)
, m_numOfStrands(numOfStrands)
{
    m_strands              = new Strand[numOfStrands];
    m_sensesOfStrands      = new rg_FLAG[numOfStrands];
    m_HBondsBetweenStrands = new HydrogenBond[numOfStrands];

    for ( rg_INT i_strand=0; i_strand<numOfStrands; i_strand++ ) {
        m_sensesOfStrands[i_strand] = 0;
    }
}



V::GeometryTier::Sheet::Sheet( const V::GeometryTier::Sheet& aSheet )
{
    m_sheetID      = aSheet.m_sheetID;
    m_numOfStrands = aSheet.m_numOfStrands;
    
    m_strands              = new Strand[m_numOfStrands];
    m_sensesOfStrands      = new rg_FLAG[m_numOfStrands];
    m_HBondsBetweenStrands = new HydrogenBond[m_numOfStrands];

    for ( rg_INT i_strand=0; i_strand<m_numOfStrands; i_strand++ ) {
        m_strands[i_strand]              = aSheet.m_strands[i_strand];
        m_sensesOfStrands[i_strand]      = aSheet.m_sensesOfStrands[i_strand];
        m_HBondsBetweenStrands[i_strand] = aSheet.m_HBondsBetweenStrands[i_strand];
    }
}



V::GeometryTier::Sheet::~Sheet()
{
    if ( m_strands != rg_NULL )
        delete [] m_strands;

    if ( m_sensesOfStrands != rg_NULL )
        delete [] m_sensesOfStrands;

    if ( m_HBondsBetweenStrands != rg_NULL )
        delete [] m_HBondsBetweenStrands;
}



//  GET FUNCTION
V::GeometryTier::Strand* V::GeometryTier::Sheet::getStrand( const rg_INT& strandArrayID )
{
    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        return &m_strands[strandArrayID];
    }
    else {
        return rg_NULL;
    }
}



V::GeometryTier::Strand* V::GeometryTier::Sheet::getStrandBySerial( const rg_INT& strandSerial )
{
    rg_INT strandArrayID = strandSerial-1;

    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        return &m_strands[strandArrayID];
    }
    else {
        return rg_NULL;
    }
}



rg_FLAG V::GeometryTier::Sheet::getSenseOfStrand( const rg_INT& strandArrayID )
{
    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        return m_sensesOfStrands[strandArrayID];
    }
    else {
        return 0;
    }
}



V::GeometryTier::HydrogenBond*   V::GeometryTier::Sheet::getHBondsBetweenStrands()
{
    return m_HBondsBetweenStrands;
}



V::GeometryTier::HydrogenBond*   V::GeometryTier::Sheet::getHBondBetweenStrands( const rg_INT& strandArrayID )
{
    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        return &m_HBondsBetweenStrands[strandArrayID];
    }
    else {
        return rg_NULL;
    }

}



V::GeometryTier::Chain*  V::GeometryTier::Sheet::getChain()
{
    return m_strands[0].getResidues()->getFirstEntity()->getChain();
}



rg_BOOL     V::GeometryTier::Sheet::isResidueInSheet( Residue* aResidue ) const
{
    rg_BOOL isInSheet = rg_FALSE;

    rg_INT i_strand=0;
    for ( i_strand=0; i_strand<m_numOfStrands; i_strand++ ) {
        isInSheet = m_strands[i_strand].isResidueInStrand( aResidue );

        if ( isInSheet == rg_TRUE )
            break;
    }
    
    return isInSheet;
}



//  SET FUNCTION

void    V::GeometryTier::Sheet::setStrand( const rg_INT& strandArrayID, Strand& aStrand )
{
    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        m_strands[strandArrayID] = aStrand;
    }
}



void    V::GeometryTier::Sheet::setSensesOfStrands( rg_FLAG* sensesOfStrands )
{
    m_sensesOfStrands = sensesOfStrands;
}



void    V::GeometryTier::Sheet::setSenseOfStrand( const rg_INT& strandArrayID, const rg_FLAG& isParallel )
{
    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        m_sensesOfStrands[strandArrayID] = isParallel;
    }
}



void    V::GeometryTier::Sheet::setSenseOfStrandBySerial( const rg_INT& strandSerial, const rg_FLAG& isParallel )
{
    rg_INT strandArrayID = strandSerial-1;

    if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        m_sensesOfStrands[strandArrayID] = isParallel;
    }
}



void V::GeometryTier::Sheet::setHBondsBetweenStrands( HydrogenBond* hBondsBetweenStrands )
{
    m_HBondsBetweenStrands = hBondsBetweenStrands;
}



void V::GeometryTier::Sheet::setHBondsBetweenStrands( const rg_INT& strandArrayID, const HydrogenBond& hBondBetweenStrands )
{
   if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        m_HBondsBetweenStrands[strandArrayID] = hBondBetweenStrands;
    }
}



void V::GeometryTier::Sheet::setHBondsBetweenStrandsBySerial( const rg_INT& strandSerial, const HydrogenBond& hBondBetweenStrands )
{
    rg_INT strandArrayID = strandSerial-1;

   if ( strandArrayID >= 0 || strandArrayID < m_numOfStrands ) {
        m_HBondsBetweenStrands[strandArrayID] = hBondBetweenStrands;
    }
}



//  OPERATOR OVERLOADING
V::GeometryTier::Sheet& V::GeometryTier::Sheet::operator =( const V::GeometryTier::Sheet& aSheet )
{
    if( this == &aSheet )
        return *this;

    m_sheetID      = aSheet.m_sheetID;
    m_numOfStrands = aSheet.m_numOfStrands;
    
    m_strands              = new Strand[m_numOfStrands];
    m_sensesOfStrands      = new rg_FLAG[m_numOfStrands];
    m_HBondsBetweenStrands = new HydrogenBond[m_numOfStrands];

    for ( rg_INT i_strand=0; i_strand<m_numOfStrands; i_strand++ ) {
        m_strands[i_strand]              = aSheet.m_strands[i_strand];
        m_sensesOfStrands[i_strand]      = aSheet.m_sensesOfStrands[i_strand];
        m_HBondsBetweenStrands[i_strand] = aSheet.m_HBondsBetweenStrands[i_strand];
    }

    return *this;
}

