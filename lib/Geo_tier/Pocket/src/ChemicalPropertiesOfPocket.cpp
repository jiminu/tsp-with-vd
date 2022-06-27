#include "ChemicalPropertiesOfPocket.h"


ChemicalPropertiesOfPocket::ChemicalPropertiesOfPocket()
{
	m_isThisComputed = rg_FALSE;
}

ChemicalPropertiesOfPocket::ChemicalPropertiesOfPocket( const ChemicalPropertiesOfPocket& chemicalProperties )
{
	m_isThisComputed		= chemicalProperties.m_isThisComputed;

	m_numOfTotalAtoms		= chemicalProperties.m_numOfTotalAtoms;
	
	m_numOfCarbonAtoms		= chemicalProperties.m_numOfCarbonAtoms;
	m_numOfNitrogenAtoms	= chemicalProperties.m_numOfNitrogenAtoms;
	m_numOfOxygenAtoms		= chemicalProperties.m_numOfOxygenAtoms;
	m_numOfSulfurAtoms		= chemicalProperties.m_numOfSulfurAtoms;

	m_numOfBackboneAtoms	= chemicalProperties.m_numOfBackboneAtoms;
	m_numOfSideChainAtoms	= chemicalProperties.m_numOfSideChainAtoms;

	m_numOfHydrophobicResidues = chemicalProperties.m_numOfHydrophobicResidues;
	m_numOfPolarResidues	= chemicalProperties.m_numOfPolarResidues;
	m_numOfChargedResidues	= chemicalProperties.m_numOfChargedResidues;

	m_numOfAlphaHelices		= chemicalProperties.m_numOfAlphaHelices;
	m_numOfBetaSheets		= chemicalProperties.m_numOfBetaSheets;
	m_numOfLoops			= chemicalProperties.m_numOfLoops;

	m_numOfHBondDonors		= chemicalProperties.m_numOfHBondDonors;
	m_numOfHBondAcceptor	= chemicalProperties.m_numOfHBondAcceptor;
}



ChemicalPropertiesOfPocket::~ChemicalPropertiesOfPocket()
{	
}

//  GET FUNCTION
rg_FLAG ChemicalPropertiesOfPocket::isThisComputed() const
{
	return m_isThisComputed;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfTotalAtoms() const
{
	return m_numOfTotalAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfCarbonAtoms() const
{
	return m_numOfCarbonAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfNitrogenAtoms() const
{
	return m_numOfNitrogenAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfOxygenAtoms() const
{
	return m_numOfOxygenAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfSulfurAtoms() const
{
	return m_numOfSulfurAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfBackboneAtoms() const
{
	return m_numOfBackboneAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfSideChainAtoms() const
{
	return m_numOfSideChainAtoms;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfHydrophobicResidues() const
{
	return m_numOfHydrophobicResidues;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfPolarResidues() const
{
	return m_numOfPolarResidues;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfChargedResidues() const
{
	return m_numOfChargedResidues;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfAlphaHelices() const
{
	return m_numOfAlphaHelices;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfBetaSheets() const
{
	return m_numOfBetaSheets;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfLoops() const
{
	return m_numOfLoops;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfHBondDonors() const
{
	return m_numOfHBondDonors;
}

rg_INT ChemicalPropertiesOfPocket::getNumOfHBondAcceptor() const
{
	return m_numOfHBondAcceptor;
}




//  SET FUNCTION
void ChemicalPropertiesOfPocket::setIsThisComputed(const rg_FLAG& flag)
{
	m_isThisComputed = flag;
}

void ChemicalPropertiesOfPocket::setNumOfTotalAtoms(const rg_INT& num)
{
	m_numOfTotalAtoms = num;
}

void ChemicalPropertiesOfPocket::setNumOfCarbonAtoms(const rg_INT& num)
{
	m_numOfCarbonAtoms = num;
}
void ChemicalPropertiesOfPocket::setNumOfNitrogenAtoms(const rg_INT& num)
{
	m_numOfNitrogenAtoms = num;
}

void ChemicalPropertiesOfPocket::setNumOfOxygenAtoms(const rg_INT& num)
{
	m_numOfOxygenAtoms = num;
}

void ChemicalPropertiesOfPocket::setNumOfSulfurAtoms(const rg_INT& num)
{
	m_numOfSulfurAtoms = num;
}

void ChemicalPropertiesOfPocket::setNumOfBackboneAtoms(const rg_INT& num)
{
	m_numOfBackboneAtoms = num;
}

void ChemicalPropertiesOfPocket::setNumOfSideChainAtoms(const rg_INT& num)
{
	m_numOfSideChainAtoms = num;
}

void ChemicalPropertiesOfPocket::setNumOfHydrophobicResidues(const rg_INT& num)
{
	m_numOfHydrophobicResidues = num;
}

void ChemicalPropertiesOfPocket::setNumOfPolarResidues(const rg_INT& num)
{
	m_numOfPolarResidues = num;
}

void ChemicalPropertiesOfPocket::setNumOfChargedResidues(const rg_INT& num)
{
	m_numOfChargedResidues = num;
}

void ChemicalPropertiesOfPocket::setNumOfAlphaHelices(const rg_INT& num)
{
	m_numOfAlphaHelices = num;
}

void ChemicalPropertiesOfPocket::setNumOfBetaSheets(const rg_INT& num)
{
	m_numOfBetaSheets = num;
}

void ChemicalPropertiesOfPocket::setNumOfLoops(const rg_INT& num)
{
	m_numOfLoops = num;
}

void ChemicalPropertiesOfPocket::setNumOfHBondDonors(const rg_INT& num)
{
	m_numOfHBondDonors = num;
}

void ChemicalPropertiesOfPocket::setNumOfHBondAcceptor(const rg_INT& num)
{
	m_numOfHBondAcceptor = num;
}


//  OPERATOR OVERLOADING
ChemicalPropertiesOfPocket& ChemicalPropertiesOfPocket::operator =(const ChemicalPropertiesOfPocket& chemicalProperties)
{
	if ( this == &chemicalProperties)
		return *this;

	m_isThisComputed		= chemicalProperties.m_isThisComputed;

	m_numOfTotalAtoms		= chemicalProperties.m_numOfTotalAtoms;
	
	m_numOfCarbonAtoms		= chemicalProperties.m_numOfCarbonAtoms;
	m_numOfNitrogenAtoms	= chemicalProperties.m_numOfNitrogenAtoms;
	m_numOfOxygenAtoms		= chemicalProperties.m_numOfOxygenAtoms;
	m_numOfSulfurAtoms		= chemicalProperties.m_numOfSulfurAtoms;

	m_numOfBackboneAtoms	= chemicalProperties.m_numOfBackboneAtoms;
	m_numOfSideChainAtoms	= chemicalProperties.m_numOfSideChainAtoms;

	m_numOfHydrophobicResidues = chemicalProperties.m_numOfHydrophobicResidues;
	m_numOfPolarResidues	= chemicalProperties.m_numOfPolarResidues;
	m_numOfChargedResidues	= chemicalProperties.m_numOfChargedResidues;

	m_numOfAlphaHelices		= chemicalProperties.m_numOfAlphaHelices;
	m_numOfBetaSheets		= chemicalProperties.m_numOfBetaSheets;
	m_numOfLoops			= chemicalProperties.m_numOfLoops;

	m_numOfHBondDonors		= chemicalProperties.m_numOfHBondDonors;
	m_numOfHBondAcceptor	= chemicalProperties.m_numOfHBondAcceptor;

	return *this;
}

