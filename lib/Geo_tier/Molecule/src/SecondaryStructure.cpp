#include "SecondaryStructure.h"
using namespace V::GeometryTier;


SecondaryStructure::SecondaryStructure()
{
}



SecondaryStructure::SecondaryStructure( const SecondaryStructure& aSecondaryStructure )
{
    m_listOfHelices = aSecondaryStructure.m_listOfHelices;
    m_listOfSheets  = aSecondaryStructure.m_listOfSheets;
    m_listOfTurns   = aSecondaryStructure.m_listOfTurns;
}



SecondaryStructure::~SecondaryStructure()
{
}



//  GET FUNCTION
rg_dList<Helix>* SecondaryStructure::getHelices()
{
    return &m_listOfHelices;
}



rg_dList<Sheet>* SecondaryStructure::getSheets()
{
    return &m_listOfSheets;
}



rg_dList<Turn>* SecondaryStructure::getTurns()
{
    return &m_listOfTurns;
}



Sheet* SecondaryStructure::getSheet( const string& sheetID )
{
    Sheet* targetSheet = rg_NULL;

    m_listOfSheets.reset4Loop();
    while ( m_listOfSheets.setNext4Loop() ) {
        Sheet* currSheet = m_listOfSheets.getpEntity();
        if ( currSheet->getSheetID().compare( sheetID ) == 0 ) {
            targetSheet = currSheet;
            break;
        }
    }
    return targetSheet;
}



//  SET FUNCTION
Helix* SecondaryStructure::addHelix( const Helix& aHelix )
{
    return m_listOfHelices.addTail( aHelix );
}



Sheet* SecondaryStructure::addSheet( const Sheet& aSheet )
{
    return m_listOfSheets.addTail( aSheet );
}



Turn* SecondaryStructure::addTurn( const Turn& aTurn )
{
    return m_listOfTurns.addTail( aTurn );
}



// QUERY FUNCTION
rg_BOOL SecondaryStructure::isResidueInSecondaryStructure( Residue* aResidue ) const
{
    if ( isResidueInHelix( aResidue ) == rg_TRUE )
        return rg_TRUE;
    
    if ( isResidueInSheet( aResidue ) == rg_TRUE )
        return rg_TRUE;

    if ( isResidueInTurn( aResidue ) == rg_TRUE )
        return rg_TRUE;

    return rg_FALSE;
}



rg_BOOL SecondaryStructure::isResidueInHelix( Residue * aResidue ) const
{
    rg_BOOL isInHelix = rg_FALSE;

    m_listOfHelices.reset4Loop();
    while ( m_listOfHelices.setNext4Loop() ) {
        isInHelix = m_listOfHelices.getpEntity()->isResidueInHelix( aResidue );

        if ( isInHelix == rg_TRUE )
            break;
    }

    return isInHelix;
}



rg_BOOL SecondaryStructure::isResidueInSheet( Residue * aResidue ) const
{
    rg_BOOL isInSheet = rg_FALSE;

    m_listOfSheets.reset4Loop();
    while ( m_listOfSheets.setNext4Loop() ) {
        isInSheet = m_listOfSheets.getpEntity()->isResidueInSheet( aResidue );
        
        if ( isInSheet == rg_TRUE )
            break;
    }

    return isInSheet;
}



rg_BOOL SecondaryStructure::isResidueInTurn( Residue * aResidue ) const	
{
    rg_BOOL isInTurn = rg_FALSE;

    m_listOfTurns.reset4Loop();
    while ( m_listOfTurns.setNext4Loop() ) {
        isInTurn = m_listOfTurns.getpEntity()->isResidueInTurn( aResidue );

        if ( isInTurn == rg_TRUE )
            break;
    }

    return isInTurn;
}



Helix* SecondaryStructure::getHelixOfResidue( Residue* aResidue )
{
    Helix*  targetHelix = rg_NULL;
    rg_BOOL isInHelix   = rg_FALSE;

    m_listOfHelices.reset4Loop();
    while ( m_listOfHelices.setNext4Loop() ) {
        Helix* currHelix = m_listOfHelices.getpEntity();
        isInHelix = currHelix->isResidueInHelix( aResidue );

        if ( isInHelix == rg_TRUE ) {
            targetHelix = currHelix;
            break;
        }
    }

    return targetHelix;
}



Sheet* SecondaryStructure::getSheetOfResidue( Residue* aResidue )
{
    Sheet*  targetSheet = rg_NULL;
    rg_BOOL isInSheet   = rg_FALSE;

    m_listOfSheets.reset4Loop();
    while ( m_listOfSheets.setNext4Loop() ) {
        Sheet* currSheet = m_listOfSheets.getpEntity();
        isInSheet = currSheet->isResidueInSheet( aResidue );

        if ( isInSheet == rg_TRUE ) {
            targetSheet = currSheet;
            break;
        }
    }

    return targetSheet;
}



Turn* SecondaryStructure::getTurnOfResidue( Residue* aResidue )
{
    Turn*   targetTurn = rg_NULL;
    rg_BOOL isInTurn   = rg_FALSE;

    m_listOfTurns.reset4Loop();
    while ( m_listOfTurns.setNext4Loop() ) {
        Turn* currTurn = m_listOfTurns.getpEntity();
        isInTurn = currTurn->isResidueInTurn( aResidue );

        if ( isInTurn == rg_TRUE ) {
            targetTurn = currTurn;
            break;
        }
    }

    return targetTurn;
}



rg_BOOL SecondaryStructure::isAtomInSecondaryStructure( Atom* anAtom ) const
{
    return isResidueInSecondaryStructure( anAtom->getResidue() );
}



rg_BOOL SecondaryStructure::isAtomInHelix( Atom * anAtom ) const
{
    return isResidueInHelix( anAtom->getResidue() );
}



rg_BOOL SecondaryStructure::isAtomInSheet( Atom * anAtom ) const
{
    return isResidueInSheet( anAtom->getResidue() );
}



rg_BOOL SecondaryStructure::isAtomInTurn( Atom * anAtom ) const
{
    return isResidueInTurn( anAtom->getResidue() );
}



Helix* SecondaryStructure::getHelixOfAtom( Atom* anAtom )
{
    return getHelixOfResidue( anAtom->getResidue() );
}



Sheet* SecondaryStructure::getSheetOfAtom( Atom* anAtom )
{
    return getSheetOfResidue( anAtom->getResidue() );
}



Turn* SecondaryStructure::getTurnOfAtom( Atom* anAtom )
{
    return getTurnOfResidue( anAtom->getResidue() );
}



//  OPERATOR OVERLOADING
SecondaryStructure& SecondaryStructure::operator =( const SecondaryStructure& aSecondaryStructure )
{
    if( this == &aSecondaryStructure )
        return *this;

    m_listOfHelices = aSecondaryStructure.m_listOfHelices;
    m_listOfSheets  = aSecondaryStructure.m_listOfSheets;
    m_listOfTurns   = aSecondaryStructure.m_listOfTurns;

    return *this;
}
