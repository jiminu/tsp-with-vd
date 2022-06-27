#ifndef _STRAND_H
#define _STRAND_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "Residue.h"



namespace V {
namespace GeometryTier {



class Sheet;

class Strand
{
private:
    rg_INT               m_serial;
    Sheet*               m_sheet;
    rg_dList<Residue*>   m_residues;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Strand();
    Strand( const rg_INT& serial );
    Strand( const rg_INT& serial, Sheet* aSheet );
    Strand( const Strand& aStrand);
    ~Strand();


    //  GET FUNCTION
    rg_INT                 getSerial() const;
    Sheet*                 getSheet();
    rg_dList<Residue*>*    getResidues();

    Chain*                 getChain();

    rg_BOOL                isResidueInStrand( Residue* aResidue ) const;

    //  SET FUNCTION
    void       setSerial( const rg_INT& serial );
    void       setSheet( Sheet* aSheet );
    
    void       addResidue( Residue* aResidue );
    void       setResidues( rg_dList<Residue*>* residues );

    //  OPERATOR OVERLOADING
    Strand& operator =( const Strand& aStrand );

};



inline  rg_INT                                   V::GeometryTier::Strand::getSerial() const  { return m_serial; }
inline  V::GeometryTier::Sheet*               V::GeometryTier::Strand::getSheet()         { return m_sheet; }
inline  rg_dList<V::GeometryTier::Residue*>*  V::GeometryTier::Strand::getResidues()      { return &m_residues; }

inline  void V::GeometryTier::Strand::setSerial(const rg_INT& serial)            { m_serial = serial; }
inline  void V::GeometryTier::Strand::setSheet(Sheet* aSheet)                    { m_sheet = aSheet; }
inline  void V::GeometryTier::Strand::addResidue(Residue* aResidue)              { m_residues.addTail(aResidue); }
inline  void V::GeometryTier::Strand::setResidues(rg_dList<Residue*>* residues)  { m_residues.append(*residues); }




} // namespace GeometryTier
} // namespace V


#endif

