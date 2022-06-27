#ifndef _RESIDUEWITHSSCODE_H
#define _RESIDUEWITHSSCODE_H

#include "rg_Const.h"
#include "ConstForMolecule.h"



namespace V {
namespace GeometryTier {



class Residue;
class Helix;
class Sheet;
class Turn;

class ResidueWithSSCode
{
private:
    SecondaryStructureCode    m_code;
    void*                     m_secondaryStructureElement;
    Residue*                  m_residue;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    ResidueWithSSCode();
    ResidueWithSSCode( const SecondaryStructureCode& aCode, void* secondaryStructureElement, Residue* aResidue );
    ResidueWithSSCode( const ResidueWithSSCode& aResidueWithSSCode );
    ~ResidueWithSSCode();


    //  GET FUNCTION
    SecondaryStructureCode    getCode() const;
    void*                     getSecondaryStructureElement();
    Residue*                  getResidue();

    rg_BOOL		isResidueInSecondaryStructure() const;

	rg_BOOL		isResidueInHelix() const;
	rg_BOOL		isResidueInSheet() const;
	rg_BOOL		isResidueInTurn() const;	

    Helix*		getHelixOfResidue() const;
	Sheet*		getSheetOfResidue() const;
	Turn*		getTurnOfResidue() const;


    //  SET FUNCTION
    void    setCode( const SecondaryStructureCode& aCode );
    void    setSecondaryStructureElement( void* aSecondaryStructureElement );
    void    setResidue( Residue* aResidue );


    //  OPERATOR OVERLOADING
    ResidueWithSSCode& operator =(const ResidueWithSSCode& aResidueWithSSCode);
};



} // namespace GeometryTier
} // namespace V


#endif

