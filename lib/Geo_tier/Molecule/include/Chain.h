#ifndef _CHAIN_H
#define _CHAIN_H

#include "rg_Const.h"
#include "ConstForMolecule.h"
#include "rg_dList.h"
#include "SecondaryStructure.h"
#include "ResidueWithSSCode.h"
#include <string>
using namespace std;



namespace V {
namespace GeometryTier {



class Residue;
class Molecule;

class Chain
{
private:
    rg_INT              m_ID;
    rg_INT              m_chainIDFromInputFileInDecimal;

    Molecule*           m_molecule;
    rg_dList<Residue*>  m_residues;
    
    ////////////////////////////////////////////////////

    SecondaryStructure  m_secondaryStructure;
    ChainCode           m_chainCode;    /* DNA or RNA ? modification required for query function ! */
    
    

public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Chain();
    Chain( const rg_INT& ID );
    Chain( const rg_INT& ID, Molecule* aMolecule );
    Chain( const rg_INT& ID, const ChainCode& aChainCode );
    Chain( const rg_INT& ID, const ChainCode& aChainCode, const rg_INT& chainIDFromInput );
    Chain( const Chain& aChain );
    ~Chain();


    //  GET FUNCTION
    rg_INT               getID() const;
    rg_dList<Residue*>*  getResidues();
    Molecule*            getMolecule();
    rg_INT               getChainIDFromInputFileInDecimal() const;
    string               getChainIDFromInputFileInString() const;
    
    //  by Youngsong Cho, 2012.03.14.
    rg_INT               getNumStandardResidues() const;

    ChainCode            getChainCode() const;
    SecondaryStructure*  getSecondaryStructure();
    
    rg_BOOL				 isProteinChain() const;
	rg_BOOL				 isDNAChain() const;
	rg_BOOL				 isRNAChain() const;
    rg_BOOL				 isLigandChain() const; // NOT WORKING.
	rg_BOOL				 isHOHChain() const;
	rg_BOOL				 isIONChain() const;

	rg_BOOL              haveMissingResidues();
	rg_INT               findMissingResidues(rg_dList<rg_INT>& seqNumbersOfMissingResidues);

	void			     getResiduesSortedBySequenceNumber( rg_dList<Residue*>& targetResidues );
	void			     getAminoResiduesSortedBySequenceNumber( rg_dList<Residue*>& aminoResiduesSortedBySeqNumbers );
    rg_INT			     getAminoResiduesSortedBySequenceNumber( Residue**& aminoResiduesSortedBySeqNumbers );
    void			     getResidueWithSSCsSortedBySequenceNumber( rg_dList<ResidueWithSSCode>& targetResidueWithSSCs );
    void			     getAtomsOnBackboneSortedBySequenceNumber( rg_dList<Atom*>& targetAtoms );

    void                 getResiduesInSecondaryStructureBySequenceNumber( rg_dList<Residue*>& targetResidues );
    void                 getResiduesNotInSecondaryStructureBySequenceNumber( rg_dList<Residue*>& targetResidues );


    //  SET FUNCTION
    void     setID( const rg_INT& ID );
    void     setChainCode( const ChainCode& aChainCode );
    void     setMolecule( Molecule* aMolecule );
    void     setChainIDFromInputFileInDecimal( const rg_INT& chainIDFromInput );

    void     evaluateAndSetChainCode();
    
    Residue* addResidue( Residue* aResidue );




    //  OPERATOR OVERLOADING
    Chain& operator =( const Chain& aChain );

};


inline  rg_INT                                  V::GeometryTier::Chain::getID() const { return m_ID; }
inline  rg_dList<V::GeometryTier::Residue*>* V::GeometryTier::Chain::getResidues() { return &m_residues; }
inline  Molecule*                               V::GeometryTier::Chain::getMolecule() { return m_molecule; }
inline  rg_INT                                  V::GeometryTier::Chain::getChainIDFromInputFileInDecimal() const { return m_chainIDFromInputFileInDecimal; }
inline  ChainCode                               V::GeometryTier::Chain::getChainCode() const { return m_chainCode; }
inline  SecondaryStructure*                     V::GeometryTier::Chain::getSecondaryStructure() { return &m_secondaryStructure; }

inline  rg_BOOL     V::GeometryTier::Chain::isProteinChain() const   { return (m_chainCode == PROTEIN_CHAIN) ? rg_TRUE : rg_FALSE; }
inline  rg_BOOL     V::GeometryTier::Chain::isDNAChain() const       { return (m_chainCode == DNA_CHAIN) ? rg_TRUE : rg_FALSE; }
inline  rg_BOOL     V::GeometryTier::Chain::isRNAChain() const       { return (m_chainCode == RNA_CHAIN) ? rg_TRUE : rg_FALSE; }
inline  rg_BOOL     V::GeometryTier::Chain::isHOHChain() const       { return (m_chainCode == HOH_CHAIN) ? rg_TRUE : rg_FALSE; }
inline  rg_BOOL     V::GeometryTier::Chain::isIONChain() const       { return (m_chainCode == ION_CHAIN) ? rg_TRUE : rg_FALSE; }


inline  void        V::GeometryTier::Chain::setID(const rg_INT& ID) { m_ID = ID; }
inline  void        V::GeometryTier::Chain::setChainCode(const ChainCode& aChainCode) { m_chainCode = aChainCode; }
inline  void        V::GeometryTier::Chain::setMolecule(Molecule* aMolecule) { m_molecule = aMolecule; }
inline  void        V::GeometryTier::Chain::setChainIDFromInputFileInDecimal(const rg_INT& chainIDFromInput) { m_chainIDFromInputFileInDecimal = chainIDFromInput; }




} // namespace GeometryTier
} // namespace V


#endif

