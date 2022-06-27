#ifndef _RESIDUE_H
#define _RESIDUE_H

#include "rg_Const.h"
#include "ConstForMolecule.h"
#include "rg_dList.h"
#include <string>

//#include "BackboneDihedralAnglePair.h"
//#include "SidechainDihedralAngleSet.h"

using namespace std;


class BackboneDihedralAnglePair;
class SidechainDihedralAngleSet;


namespace V {
namespace GeometryTier {



class Atom;
class Chain;

class Residue
{
private:
    rg_INT            m_ID;
    ResidueCode       m_residueCode;
    rg_dList<Atom*>   m_atoms;
    Chain*            m_chain;

    rg_INT            m_sequenceNumber;
    string            m_residueName;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Residue();
    Residue( const rg_INT& ID );
    Residue( const rg_INT& ID, const ResidueCode& aResidueCode );
    Residue( const rg_INT& ID, const ResidueCode& aResidueCode, const rg_INT& sequenceNumber );
    Residue( const Residue& aResidue );
    ~Residue();


    //  GET FUNCTION
    rg_INT              getID() const;
    ResidueCode         getResidueCode() const;
    rg_dList<Atom*>*    getAtoms();
	void                getAtoms(rg_dList<Atom*>& atomsOfResidue);
	rg_INT              getAtoms( Atom**& atomsOfResidue );
    Atom*               getAtom( const string& atomNameFromInputFile );
    Chain*              getChain();
    rg_INT              getSequenceNumber() const;
    string              getResidueName() const;
    string              getResidueNameInPDBFormat() const;


    Atom*               getAlphaCarbonOfAminoResidue() const;
    Atom*               getNitrogenInAminoGroupOfAminoResidue() const;
    Atom*               getCarbonInCarboxylGroupOfAminoResidue() const;
    Atom*               getOxygenInCarboxylGroupOfAminoResidue() const;
    Atom*               getBetaCarbonInSideChainOfAminoResidue() const;
    Atom*               getGammaCarbonInSideChainOfAminoResidue() const;

    Atom*               getAtom( const AtomCode& typeOfAtom, const RemoteIndicator& remote, const BranchDesignator& branch ) const;
	Atom*               getAtom( const RemoteIndicator& remote, const BranchDesignator& branch ) const;
	void                getAtoms( const RemoteIndicator& remote, rg_dList<Atom*>& atomsWithThisRemote ) const;    
    
	void                getAtomsOnBackbone( rg_dList<Atom*>& atomsOnBackbone );
    void                getAtomsOnBackbone( rg_dList<Atom*>* listOfAtomsOnBackbone );
	rg_INT              getAtomsOnBackbone( Atom**& atomsOnBackbone );
    void                getAtomsOnSideChain( rg_dList<Atom*>* listOfAtomsOnSideChain );
	rg_INT              getAtomsOnSidechain( rg_dList<Atom*>& atomsOnSidechain );
	rg_INT              getAtomsOnSidechain( Atom**& atomsOnSidechain );
    void                getAtomsNotOnSideChain( rg_dList<Atom*>* listOfAtomsNotOnSideChain ); // Difference with "getAtomsOnBackbone" : O atom.

	Residue*            getResidueConnectedByNitrogenInAminoGroupOfAminoResidue() const;
	Residue*            getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue() const;

    rg_BOOL             isResidueType( const ResidueType& typeOfResidue ) const;
    rg_BOOL             isAminoResidue() const;
    rg_BOOL             isNeucleicAcidResidue() const;
    rg_BOOL             isDNAResidue() const;
    rg_BOOL             isRNAResidue() const;
    rg_BOOL             isStandardResidue() const;

	rg_BOOL             getBackBoneDihedralAngles(rg_REAL& phi, rg_REAL& psi) const;
	rg_INT              getSideChainDihedralAngles(rg_REAL*& dihedralAngles)  ;
	rg_BOOL             hasSideChainDihedralAngle() const;
	rg_BOOL             isN_Terminal() const;
	rg_BOOL             isC_Terminal() const;
    void                changeConformationWith(const BackboneDihedralAnglePair& backboneDihedralAnglePair, const SidechainDihedralAngleSet& sideChainDihedralAngleSet);

    //  SET FUNCTION
    void    setID( const rg_INT& ID );
    void    setResidueCode( const ResidueCode& aResidueCode );
    void    setChain( Chain* aChain);
    void    setSequenceNumber( const rg_INT& sequenceNumber );
    void    setResidueName( const string& residueName );
    
    Atom*   addAtom( Atom* atom );

	void    replaceAtomsWith(rg_dList<Atom*>& atomsToBeReplaced);

	rg_REAL computeSumOfPairwiseAtomXVolumesWithOtherResidue(Residue* residue);
	rg_REAL computeSumOfPairwiseAtomXVolumesOfSidechainWithOtherResidue(Residue* residue);
	rg_REAL computeSumOfPairwiseAtomXVolumesOfSidechainWithItsBackbone();

    //  OPERATOR OVERLOADING
    Residue& operator =( const Residue& residue );

    static bool SerialLess(Residue* res1, Residue* res2) {
        return res1->m_sequenceNumber < res2->m_sequenceNumber;
    }

};


inline  rg_INT                                  V::GeometryTier::Residue::getID() const              { return m_ID; }
inline  V::GeometryTier::ResidueCode         V::GeometryTier::Residue::getResidueCode() const     { return m_residueCode; }
inline  rg_dList<V::GeometryTier::Atom*>*    V::GeometryTier::Residue::getAtoms()                 { return &m_atoms; }
inline  V::GeometryTier::Chain*              V::GeometryTier::Residue::getChain()                 { return m_chain; }
inline  rg_INT                                  V::GeometryTier::Residue::getSequenceNumber() const  { return m_sequenceNumber; }
inline  string                                  V::GeometryTier::Residue::getResidueName() const     { return m_residueName; }

inline  V::GeometryTier::Atom*   V::GeometryTier::Residue::getAlphaCarbonOfAminoResidue() const            { return getAtom(C_ATOM, ALPHA_REMOTE, UNK_BRANCH); }
inline  V::GeometryTier::Atom*   V::GeometryTier::Residue::getNitrogenInAminoGroupOfAminoResidue() const   { return getAtom(N_ATOM, UNK_REMOTE, UNK_BRANCH); }
inline  V::GeometryTier::Atom*   V::GeometryTier::Residue::getCarbonInCarboxylGroupOfAminoResidue() const  { return getAtom(C_ATOM, UNK_REMOTE, UNK_BRANCH); }
inline  V::GeometryTier::Atom*   V::GeometryTier::Residue::getBetaCarbonInSideChainOfAminoResidue() const  { return getAtom(C_ATOM, BETA_REMOTE, UNK_BRANCH); }
inline  V::GeometryTier::Atom*   V::GeometryTier::Residue::getGammaCarbonInSideChainOfAminoResidue() const { return getAtom(C_ATOM, GAMMA_REMOTE, UNK_BRANCH); }

inline  void    V::GeometryTier::Residue::setID(const rg_INT& ID)                            { m_ID = ID; }
inline  void    V::GeometryTier::Residue::setResidueCode(const ResidueCode& aResidueCode)    { m_residueCode = aResidueCode; }
inline  void    V::GeometryTier::Residue::setChain(Chain* aChain)                            { m_chain = aChain; }
inline  void    V::GeometryTier::Residue::setSequenceNumber(const rg_INT& sequenceNumber)    { m_sequenceNumber = sequenceNumber; }
inline  void    V::GeometryTier::Residue::setResidueName(const string& residueName)          { m_residueName = residueName; }

inline  V::GeometryTier::Atom*   V::GeometryTier::Residue::addAtom(V::GeometryTier::Atom* atom) { return *m_atoms.addTail(atom); }
inline  void    V::GeometryTier::Residue::replaceAtomsWith(rg_dList<V::GeometryTier::Atom*>& atomsToBeReplaced) { m_atoms = atomsToBeReplaced; }


} // namespace GeometryTier
} // namespace V


#endif

