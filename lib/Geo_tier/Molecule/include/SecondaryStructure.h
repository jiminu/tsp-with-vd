#ifndef _SECONDARYSTRUCTURE_H
#define _SECONDARYSTRUCTURE_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "ConstForMolecule.h"
#include "Helix.h"
#include "Sheet.h"
#include "Turn.h"




namespace V {
namespace GeometryTier {



class SecondaryStructure
{
private:
    rg_dList<Helix> m_listOfHelices;
    rg_dList<Sheet> m_listOfSheets;
    rg_dList<Turn>  m_listOfTurns;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    SecondaryStructure();
    SecondaryStructure( const SecondaryStructure& aSecondaryStructure );
    ~SecondaryStructure();


    //  GET FUNCTION
    rg_dList<Helix>*    getHelices();
    rg_dList<Sheet>*    getSheets();
    rg_dList<Turn>*     getTurns();

    Sheet*              getSheet( const string& sheetID );


    //  SET FUNCTION
    Helix*  addHelix( const Helix& aHelix );
    Sheet*  addSheet( const Sheet& aSheet );
    Turn*   addTurn( const Turn& aTurn );


    // QUERY FUNCTION
    rg_BOOL		isResidueInSecondaryStructure( Residue* aResidue ) const;

	rg_BOOL		isResidueInHelix( Residue * aResidue ) const;
	rg_BOOL		isResidueInSheet( Residue * aResidue ) const;
	rg_BOOL		isResidueInTurn( Residue * aResidue ) const;	

    Helix*		getHelixOfResidue( Residue* aResidue );
	Sheet*		getSheetOfResidue( Residue* aResidue );
	Turn*		getTurnOfResidue( Residue* aResidue );

    
    rg_BOOL		isAtomInSecondaryStructure( Atom* anAtom ) const;

	rg_BOOL		isAtomInHelix( Atom * anAtom ) const;
	rg_BOOL		isAtomInSheet( Atom * anAtom ) const;
	rg_BOOL		isAtomInTurn( Atom * anAtom ) const;	

    Helix*		getHelixOfAtom( Atom* anAtom );
	Sheet*		getSheetOfAtom( Atom* anAtom );
	Turn*		getTurnOfAtom( Atom* anAtom );
    



    //  OPERATOR OVERLOADING
    SecondaryStructure& operator =( const SecondaryStructure& aSecondaryStructure );
};



} // namespace GeometryTier
} // namespace V


#endif

