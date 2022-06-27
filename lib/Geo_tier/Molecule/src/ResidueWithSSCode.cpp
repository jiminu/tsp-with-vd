#include "ResidueWithSSCode.h"
#include "Residue.h"
#include "Helix.h"
#include "Sheet.h"
#include "Turn.h"
using namespace V::GeometryTier;


ResidueWithSSCode::ResidueWithSSCode()
{
}



ResidueWithSSCode::ResidueWithSSCode( const ResidueWithSSCode& aResidueWithSSCode )
{
    m_code                      = aResidueWithSSCode.m_code;
    m_secondaryStructureElement = aResidueWithSSCode.m_secondaryStructureElement;
    m_residue                   = aResidueWithSSCode.m_residue;
}



ResidueWithSSCode::ResidueWithSSCode( const SecondaryStructureCode& aCode, void* secondaryStructureElement, Residue* aResidue )
{
    m_code                      = aCode;
    m_secondaryStructureElement = secondaryStructureElement;
    m_residue                   = aResidue;
}



ResidueWithSSCode::~ResidueWithSSCode()
{
}



//  GET FUNCTION
SecondaryStructureCode ResidueWithSSCode::getCode() const
{
    return m_code;
}



void* ResidueWithSSCode::getSecondaryStructureElement()
{
    return m_secondaryStructureElement;
}



Residue* ResidueWithSSCode::getResidue()
{
    return m_residue;
}



rg_BOOL  ResidueWithSSCode::isResidueInSecondaryStructure() const
{
    if ( m_code != SSC_UNK )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_BOOL  ResidueWithSSCode::isResidueInHelix() const
{
    if ( m_code == SSC_HELIX )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_BOOL  ResidueWithSSCode::isResidueInSheet() const
{
    if ( m_code == SSC_SHEET )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_BOOL  ResidueWithSSCode::isResidueInTurn() const
{
    if ( m_code == SSC_TURN )
        return rg_TRUE;
    else
        return rg_FALSE;
}



Helix*  ResidueWithSSCode::getHelixOfResidue() const
{
    if ( m_code == SSC_HELIX )
        return (Helix*)m_secondaryStructureElement;
    else
        return rg_NULL;
}



Sheet*  ResidueWithSSCode::getSheetOfResidue() const
{
    if ( m_code == SSC_SHEET )
        return (Sheet*)m_secondaryStructureElement;
    else
        return rg_NULL;
}



Turn*  ResidueWithSSCode::getTurnOfResidue() const
{
    if ( m_code == SSC_TURN )
        return (Turn*)m_secondaryStructureElement;
    else
        return rg_NULL;
}



//  SET FUNCTION
void ResidueWithSSCode::setCode( const SecondaryStructureCode& aCode )
{
    m_code = aCode;
}



void ResidueWithSSCode::setSecondaryStructureElement( void* aSecondaryStructureElement )
{
    m_secondaryStructureElement = aSecondaryStructureElement;
}



void ResidueWithSSCode::setResidue( Residue* aResidue )
{
    m_residue = aResidue;
}



//  OPERATOR OVERLOADING
ResidueWithSSCode& ResidueWithSSCode::operator =(const ResidueWithSSCode& aResidueWithSSCode)
{
    if( this == &aResidueWithSSCode )
        return *this;

    m_code                      = aResidueWithSSCode.m_code;
    m_secondaryStructureElement = aResidueWithSSCode.m_secondaryStructureElement;
    m_residue                   = aResidueWithSSCode.m_residue;

    return *this;
}
