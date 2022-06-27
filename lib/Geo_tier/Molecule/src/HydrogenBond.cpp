#include "HydrogenBond.h"
#include "ForceField.h"
#include "ChemicalBond.h"



V::GeometryTier::HydrogenBond::HydrogenBond()
: m_donor(rg_NULL)
, m_hydrogen(rg_NULL)
, m_acceptor(rg_NULL)
, m_distDonAndAcc(0.)
, m_distHydAndAcc(0.)
, m_angleDonHydAcc(0.)
{
}



V::GeometryTier::HydrogenBond::HydrogenBond(V::GeometryTier::Atom* hydrogen, V::GeometryTier::Atom* acceptor )
: m_hydrogen(hydrogen)
, m_acceptor(acceptor)
, m_distDonAndAcc(0.)
, m_distHydAndAcc(0.)
, m_angleDonHydAcc(0.)
{
    findAndSetDonorFromHydrogen();
    computeAndSetDistanceBetweenDonorAndAcceptor();
    computeAndSetDistanceBetweenHydrogenAndAcceptor();
    computeAndSetAngleBetweenDonorHydrogenAcceptor();
}



V::GeometryTier::HydrogenBond::HydrogenBond(V::GeometryTier::Atom* donor, V::GeometryTier::Atom* hydrogen, V::GeometryTier::Atom* acceptor  )
: m_donor(donor)
, m_hydrogen(hydrogen)
, m_acceptor(acceptor)
, m_distDonAndAcc(0.)
, m_distHydAndAcc(0.)
, m_angleDonHydAcc(0.)
{
    computeAndSetDistanceBetweenDonorAndAcceptor();
    computeAndSetDistanceBetweenHydrogenAndAcceptor();
    computeAndSetAngleBetweenDonorHydrogenAcceptor();
}



V::GeometryTier::HydrogenBond::HydrogenBond( const V::GeometryTier::HydrogenBond& hBond )
: m_donor(hBond.m_donor)
, m_hydrogen(hBond.m_hydrogen)
, m_acceptor(hBond.m_acceptor)
, m_distDonAndAcc(hBond.m_distDonAndAcc)
, m_distHydAndAcc(hBond.m_distHydAndAcc)
, m_angleDonHydAcc(hBond.m_angleDonHydAcc)
{
}



V::GeometryTier::HydrogenBond::~HydrogenBond()
{
}



rg_REAL V::GeometryTier::HydrogenBond::getHydrogenBondEnergy() const
{
    rg_REAL energy = 0.0;

    rg_REAL distHydAndAcc    = m_distHydAndAcc;
    rg_INT  idForHBondCoeffs = m_acceptor->getAtomCode() + NEG_VALUE_FOR_HBOND_COEFFS;
    rg_REAL hBondCoeff_A     = HYDROGEN_BOND_COEFFS[idForHBondCoeffs][0];
    rg_REAL hBondCoeff_B     = HYDROGEN_BOND_COEFFS[idForHBondCoeffs][1];

    rg_REAL denom = pow( m_distHydAndAcc, 10.0 );

    energy = (   HYDROGEN_BOND_COEFFS[idForHBondCoeffs][0] / ( denom * pow( m_distHydAndAcc, 2.0 ) )
               - HYDROGEN_BOND_COEFFS[idForHBondCoeffs][1] / ( denom ) );

    return energy;
}





void V::GeometryTier::HydrogenBond::setDonorAndAcceptor(V::GeometryTier::Atom* donor, V::GeometryTier::Atom* acceptor )
{
    m_donor    = donor;
    m_acceptor = acceptor;

    computeAndSetDistanceBetweenDonorAndAcceptor();

    if( m_hydrogen != rg_NULL ) {
        computeAndSetDistanceBetweenHydrogenAndAcceptor();
        computeAndSetAngleBetweenDonorHydrogenAcceptor();
    }
}



void V::GeometryTier::HydrogenBond::setAtoms(V::GeometryTier::Atom* donor, V::GeometryTier::Atom* hydrogen, V::GeometryTier::Atom* acceptor )
{
    m_donor    = donor;
    m_hydrogen = hydrogen;
    m_acceptor = acceptor;

    computeAndSetDistanceBetweenDonorAndAcceptor();
    computeAndSetDistanceBetweenHydrogenAndAcceptor();
    computeAndSetAngleBetweenDonorHydrogenAcceptor();
}



rg_BOOL V::GeometryTier::HydrogenBond::evaluateAndSetDonorAndAcceptorForProteinAtoms(V::GeometryTier::Atom* firstAtom, V::GeometryTier::Atom* secondAtom )
{
    rg_BOOL isOrdinaryHBond = rg_FALSE;
    if ( firstAtom->getAtomCode() == N_ATOM && secondAtom->getAtomCode() == O_ATOM ) {
        m_donor    = firstAtom;
        m_acceptor = secondAtom;
        computeAndSetDistanceBetweenDonorAndAcceptor();
        
        if( m_hydrogen != rg_NULL ) {
            computeAndSetDistanceBetweenHydrogenAndAcceptor();
            computeAndSetAngleBetweenDonorHydrogenAcceptor();
        }

        isOrdinaryHBond = rg_TRUE;
    }
    else if ( firstAtom->getAtomCode() == O_ATOM && secondAtom->getAtomCode() == N_ATOM ) {
        m_donor    = secondAtom;
        m_acceptor = firstAtom;
        computeAndSetDistanceBetweenDonorAndAcceptor();
        
        if( m_hydrogen != rg_NULL ) {
            computeAndSetDistanceBetweenHydrogenAndAcceptor();
            computeAndSetAngleBetweenDonorHydrogenAcceptor();
        }

        isOrdinaryHBond = rg_TRUE;
    }
    else {
        isOrdinaryHBond = rg_FALSE;
    }

    return isOrdinaryHBond;
}



V::GeometryTier::HydrogenBond& V::GeometryTier::HydrogenBond::operator=( const V::GeometryTier::HydrogenBond& hBond )
{
    if( this == &hBond )
        return *this;

    m_donor             = hBond.m_donor;
    m_hydrogen          = hBond.m_hydrogen;
    m_acceptor          = hBond.m_acceptor;
    m_distDonAndAcc     = hBond.m_distDonAndAcc;
    m_distHydAndAcc     = hBond.m_distHydAndAcc;
    m_angleDonHydAcc    = hBond.m_angleDonHydAcc;

    return *this;
}



void V::GeometryTier::HydrogenBond::findAndSetDonorFromHydrogen()
{
    ChemicalBond* chemBond  = m_hydrogen->getListChemicalBond()->getFirstEntity();

    if( chemBond->getFirstAtom()->getAtomCode() == H_ATOM )
        m_donor = chemBond->getSecondAtom();
    else if( chemBond->getSecondAtom()->getAtomCode() == H_ATOM )
        m_donor = chemBond->getFirstAtom();
    else{}
}



void V::GeometryTier::HydrogenBond::computeAndSetDistanceBetweenDonorAndAcceptor()
{
    m_distDonAndAcc = ( m_donor->getAtomBall().getCenter() - m_acceptor->getAtomBall().getCenter() ).magnitude();
}



void V::GeometryTier::HydrogenBond::computeAndSetDistanceBetweenHydrogenAndAcceptor()
{
    m_distHydAndAcc = ( m_hydrogen->getAtomBall().getCenter() - m_acceptor->getAtomBall().getCenter() ).magnitude();
}



void V::GeometryTier::HydrogenBond::computeAndSetAngleBetweenDonorHydrogenAcceptor()
{
    rg_Point3D firstVector  = m_donor->getAtomBall().getCenter() - m_hydrogen->getAtomBall().getCenter();
	rg_Point3D secondVector = m_acceptor->getAtomBall().getCenter() - m_hydrogen->getAtomBall().getCenter();
	
	m_angleDonHydAcc = firstVector.angle(secondVector) * 180 / rg_PI;
}



