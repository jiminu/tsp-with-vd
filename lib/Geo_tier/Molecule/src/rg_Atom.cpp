#include "rg_Atom.h"

#include "ChemicalBond.h"
#include "Residue.h"

#include "StringFunctions.h"



V::GeometryTier::Atom::Atom()
: m_ID(0), m_atomCode(UNK_ATOM), m_residue(rg_NULL)
{
    m_serialFromInputFile = 0;
    m_chemicalBonds.removeAll();
}




V::GeometryTier::Atom::Atom( const rg_INT& ID )
: m_ID(ID), m_atomCode(UNK_ATOM), m_residue(rg_NULL)
{
    m_serialFromInputFile = 0;
    m_chemicalBonds.removeAll();
}



V::GeometryTier::Atom::Atom( const V::GeometryTier::Atom& atom )
: m_ID(atom.m_ID), m_atomCode(atom.m_atomCode), m_atomBall(atom.m_atomBall), m_chemProperties(atom.m_chemProperties), m_chemicalBonds(atom.m_chemicalBonds), m_residue(atom.m_residue)
{
    m_serialFromInputFile   = atom.m_serialFromInputFile;
    m_atomNameFromInputFile = atom.m_atomNameFromInputFile;
}



V::GeometryTier::Atom::~Atom()
{
}



V::GeometryTier::Chain* V::GeometryTier::Atom::getChain() 
{ 
    return m_residue->getChain(); 
}

string  V::GeometryTier::Atom::getAtomNameInPDBFormat() const
{
    string atomName;
    if (m_residue == rg_NULL) {
        atomName = m_atomNameFromInputFile;
    }
    else {
        if (!m_residue->isStandardResidue()) {
            string atomSymbol = (string)ATOM_FEATURES[(int)m_atomCode].symbol;
            atomName = atomSymbol + "  ";
        }
        else {
            string extraBranch  = m_chemProperties.getExtraBrangeDesignatorInString();
            string atomName     = StringFunctions::strTrim((string)ATOM_FEATURES[(int)m_atomCode].symbol);
            string remoteInd    = m_chemProperties.getRemoteIndicatorInString();
            string branchDes    = m_chemProperties.getBrangeDesignatorInString();

            if (atomName.length() == 2)
                extraBranch = "";

            atomName = extraBranch + atomName + remoteInd + branchDes;
        }
    }

    return atomName;
}



V::GeometryTier::Atom& V::GeometryTier::Atom::operator=( const V::GeometryTier::Atom& atom )
{
    if( this == &atom )
        return *this;

    m_ID                = atom.m_ID;
    m_atomCode          = atom.m_atomCode;
    m_atomBall          = atom.m_atomBall;
    m_chemProperties    = atom.m_chemProperties;

    m_chemicalBonds     = atom.m_chemicalBonds;
    m_residue           = atom.m_residue;

    m_serialFromInputFile   = atom.m_serialFromInputFile;
    m_atomNameFromInputFile = atom.m_atomNameFromInputFile;

    return *this;
}



rg_BOOL V::GeometryTier::Atom::hasChemicalBondWith(V::GeometryTier::Atom* atom)
{
    rg_BOOL bHasChemicalBond = rg_FALSE;

    m_chemicalBonds.reset4Loop();
    while (m_chemicalBonds.setNext4Loop())
    {
        ChemicalBond* currBond = m_chemicalBonds.getEntity();
        if (atom == currBond->getBondedAtom(this))
        {
            bHasChemicalBond = rg_TRUE;
            break;
        }
    }

    return bHasChemicalBond;
}



rg_REAL V::GeometryTier::Atom::compute_Lennard_Jones_potential(const V::GeometryTier::Atom & atom)
{
    rg_Point3D centerOfThisAtom = m_atomBall.getCenter();
    rg_INT     atomTypeInAmberForThisAtom = (rg_INT)(m_chemProperties.getAtomTypeInAmber());

    rg_Point3D centerOfAtom = atom.m_atomBall.getCenter();
    rg_INT     atomTypeInAmberForAtom = (rg_INT)(atom.m_chemProperties.getAtomTypeInAmber());

    rg_REAL squaredDistance = centerOfAtom.squaredDistance(centerOfThisAtom);
    rg_REAL sixthPoweredDistance = squaredDistance*squaredDistance*squaredDistance;

    rg_REAL potential = 0.0;
    potential = ((L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForThisAtom][atomTypeInAmberForAtom][0] / (sixthPoweredDistance*sixthPoweredDistance))
        - (L_J_POTENTIAL_COEFFS_NON_BONDED_ATOM_PAIR[atomTypeInAmberForThisAtom][atomTypeInAmberForAtom][1] / sixthPoweredDistance));

    return potential;
}
