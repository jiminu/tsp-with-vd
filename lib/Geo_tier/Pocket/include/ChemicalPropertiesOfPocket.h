#ifndef CHEMICAL_PROPERTIES_OF_POCKET_H_
#define CHEMICAL_PROPERTIES_OF_POCKET_H_

#include "rg_Const.h"

class ChemicalPropertiesOfPocket
{
private:
	rg_FLAG		m_isThisComputed;

	rg_INT		m_numOfTotalAtoms;
	
	rg_INT		m_numOfCarbonAtoms;
	rg_INT		m_numOfNitrogenAtoms;
	rg_INT		m_numOfOxygenAtoms;
	rg_INT		m_numOfSulfurAtoms;

	rg_INT		m_numOfBackboneAtoms;
	rg_INT		m_numOfSideChainAtoms;

	rg_INT		m_numOfHydrophobicResidues;
	rg_INT		m_numOfPolarResidues;
	rg_INT		m_numOfChargedResidues;

	rg_INT		m_numOfAlphaHelices;
	rg_INT		m_numOfBetaSheets;
	rg_INT		m_numOfLoops;

	rg_INT		m_numOfHBondDonors;
	rg_INT		m_numOfHBondAcceptor;

public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    ChemicalPropertiesOfPocket();
    ChemicalPropertiesOfPocket( const ChemicalPropertiesOfPocket& chemicalProperties );
    ~ChemicalPropertiesOfPocket();

    //  GET FUNCTION
	rg_FLAG		isThisComputed() const;

	rg_INT		getNumOfTotalAtoms() const;
	
	rg_INT		getNumOfCarbonAtoms() const;
	rg_INT		getNumOfNitrogenAtoms() const;
	rg_INT		getNumOfOxygenAtoms() const;
	rg_INT		getNumOfSulfurAtoms() const;

	rg_INT		getNumOfBackboneAtoms() const;
	rg_INT		getNumOfSideChainAtoms() const;

	rg_INT		getNumOfHydrophobicResidues() const;
	rg_INT		getNumOfPolarResidues() const;
	rg_INT		getNumOfChargedResidues() const;

	rg_INT		getNumOfAlphaHelices() const;
	rg_INT		getNumOfBetaSheets() const;
	rg_INT		getNumOfLoops() const;

	rg_INT		getNumOfHBondDonors() const;
	rg_INT		getNumOfHBondAcceptor() const;

    
    //  SET FUNCTION
	void		setIsThisComputed(const rg_FLAG& flag);

	void		setNumOfTotalAtoms(const rg_INT& num);

	void		setNumOfCarbonAtoms(const rg_INT& num);
	void		setNumOfNitrogenAtoms(const rg_INT& num);
	void		setNumOfOxygenAtoms(const rg_INT& num);
	void		setNumOfSulfurAtoms(const rg_INT& num);

	void		setNumOfBackboneAtoms(const rg_INT& num);
	void		setNumOfSideChainAtoms(const rg_INT& num);

	void		setNumOfHydrophobicResidues(const rg_INT& num);
	void		setNumOfPolarResidues(const rg_INT& num);
	void		setNumOfChargedResidues(const rg_INT& num);

	void		setNumOfAlphaHelices(const rg_INT& num);
	void		setNumOfBetaSheets(const rg_INT& num);
	void		setNumOfLoops(const rg_INT& num);

	void		setNumOfHBondDonors(const rg_INT& num);
	void		setNumOfHBondAcceptor(const rg_INT& num);

    
    //  OPERATOR OVERLOADING
    ChemicalPropertiesOfPocket& operator =(const ChemicalPropertiesOfPocket& chemicalProperties);
};


#endif