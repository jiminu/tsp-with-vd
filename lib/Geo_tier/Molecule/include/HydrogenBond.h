#ifndef _HYDROGENBOND_H
#define _HYDROGENBOND_H

#include "ConstForMolecule.h"
#include "rg_Atom.h"



namespace V {
namespace GeometryTier {



class HydrogenBond
{
private:
    Atom*           m_donor;
    Atom*           m_hydrogen;
    Atom*           m_acceptor;
    rg_REAL         m_distDonAndAcc;
    rg_REAL         m_distHydAndAcc;
    rg_REAL         m_angleDonHydAcc;

public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    HydrogenBond();
    HydrogenBond( Atom* hydrogen, Atom* acceptor );
    HydrogenBond( Atom* donor, Atom* hydrogen, Atom* acceptor );
    HydrogenBond( const HydrogenBond& hBond );
    ~HydrogenBond();

    //  GET FUNCTION
    Atom*                    getDonor();
    Atom*                    getHydrogen();
    Atom*                    getAcceptor();
    rg_REAL                  getDistanceBetweenDonorAndAcceptor() const;
    rg_REAL                  getDistanceBetweenHydrogenAndAcceptor() const;
    rg_REAL                  getAngleBetweenDonorHydrogenAcceptor() const;

    rg_REAL                  getHydrogenBondEnergy() const;


    //  SET FUNCTION
    void                     setDonor( Atom* donor );
    void                     setHydrogen( Atom* hydrogen );
    void                     setAcceptor( Atom* acceptor );
    void                     setDonorAndAcceptor( Atom* donor, Atom* acceptor );
    void                     setAtoms( Atom* donor, Atom* hydrogen, Atom* acceptor );
    
    // Automatically evaluate donor and acceptor atoms. N = donor, O = acceptor
    rg_BOOL                  evaluateAndSetDonorAndAcceptorForProteinAtoms( Atom* firstAtom, Atom* secondAtom );


    //  OPERATOR OVERLOADING
    HydrogenBond& operator =(const HydrogenBond& hBond);

private:
    void                    findAndSetDonorFromHydrogen();
    void                    computeAndSetDistanceBetweenDonorAndAcceptor();
    void                    computeAndSetDistanceBetweenHydrogenAndAcceptor();
    void                    computeAndSetAngleBetweenDonorHydrogenAcceptor();
};


inline  V::GeometryTier::Atom*   V::GeometryTier::HydrogenBond::getHydrogen() { return m_hydrogen; };
inline  V::GeometryTier::Atom*   V::GeometryTier::HydrogenBond::getAcceptor() { return m_acceptor; };
inline  V::GeometryTier::Atom*   V::GeometryTier::HydrogenBond::getDonor()    { return m_donor; }
inline  rg_REAL V::GeometryTier::HydrogenBond::getDistanceBetweenDonorAndAcceptor() const    { return m_distDonAndAcc; };
inline  rg_REAL V::GeometryTier::HydrogenBond::getDistanceBetweenHydrogenAndAcceptor() const { return m_distHydAndAcc; };
inline  rg_REAL V::GeometryTier::HydrogenBond::getAngleBetweenDonorHydrogenAcceptor() const  { return m_angleDonHydAcc; };


inline  void    V::GeometryTier::HydrogenBond::setDonor(V::GeometryTier::Atom* donor)         { m_donor = donor; }
inline  void    V::GeometryTier::HydrogenBond::setHydrogen(V::GeometryTier::Atom* hydrogen)   { m_hydrogen = hydrogen; }
inline  void    V::GeometryTier::HydrogenBond::setAcceptor(V::GeometryTier::Atom* acceptor)   { m_acceptor = acceptor; }







} // namespace GeometryTier
} // namespace V


#endif

