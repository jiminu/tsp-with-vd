#ifndef _SHEET_H
#define _SHEET_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "Strand.h"
#include "HydrogenBond.h"

#include <string>
using namespace std;



namespace V {
namespace GeometryTier {



class Chain;

class Sheet
{
private:

    string              m_sheetID;
    rg_INT              m_numOfStrands;
    Strand*             m_strands;               // size of array = m_numOfStrands. 

    rg_FLAG*            m_sensesOfStrands;       // size of array = m_numOfStrands. 
                                                 // first strand = 0, isParallel = 1 or -1.

    HydrogenBond*       m_HBondsBetweenStrands;  // size of array = m_numOfStrands. 
                                                 // m_HBondsBetweenStrands[0]= empty.
                                                 // m_HBondsBetweenStrands[1]= HBonds between m_strands[0] and m_strands[1]

    
public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Sheet();
    Sheet( const string& sheetID );
    Sheet( const string& sheetID, const rg_INT& numOfStrands );
    Sheet( const Sheet& aSheet );
    ~Sheet();


    //  GET FUNCTION
    string              getSheetID() const;
    rg_INT              getNumOfStrands() const;
    Strand*             getStrands();
    Strand*             getStrand( const rg_INT& strandArrayID );
    Strand*             getStrandBySerial( const rg_INT& strandSerial );

    rg_FLAG*            getSensesOfStrands();
    rg_FLAG             getSenseOfStrand( const rg_INT& strandArrayID );
    
    HydrogenBond*       getHBondsBetweenStrands();
    HydrogenBond*       getHBondBetweenStrands( const rg_INT& strandArrayID );

    Chain*              getChain();

    rg_BOOL             isResidueInSheet( Residue* aResidue ) const;


    //  SET FUNCTION
    void    setSheetID( const string& sheetID );
    void    setNumOfStrands( const rg_INT& numOfStrands );
    void    setStrands( Strand* strands );
    void    setStrand( const rg_INT& strandArrayID, Strand& aStrand );
    
    void    setSensesOfStrands( rg_FLAG* sensesOfStrands );
    void    setSenseOfStrand( const rg_INT& strandArrayID, const rg_FLAG& isParallel );
    void    setSenseOfStrandBySerial( const rg_INT& strandSerial, const rg_FLAG& isParallel );

    void    setHBondsBetweenStrands( HydrogenBond* hBondsBetweenStrands );
    void    setHBondsBetweenStrands( const rg_INT& strandArrayID, const HydrogenBond& hBondBetweenStrands );
    void    setHBondsBetweenStrandsBySerial( const rg_INT& strandSerial, const HydrogenBond& hBondBetweenStrands );

   
    //  OPERATOR OVERLOADING
    Sheet& operator =( const Sheet& aSheet );
};


inline  string                      V::GeometryTier::Sheet::getSheetID() const       { return m_sheetID; }
inline  rg_INT                      V::GeometryTier::Sheet::getNumOfStrands() const  { return m_numOfStrands; }
inline  V::GeometryTier::Strand* V::GeometryTier::Sheet::getStrands()             { return m_strands; }
inline  rg_FLAG*                    V::GeometryTier::Sheet::getSensesOfStrands()     { return m_sensesOfStrands; }

inline  void V::GeometryTier::Sheet::setSheetID(const string& sheetID)           { m_sheetID = sheetID; }
inline  void V::GeometryTier::Sheet::setNumOfStrands(const rg_INT& numOfStrands) { m_numOfStrands = numOfStrands; }
inline  void V::GeometryTier::Sheet::setStrands(Strand* strands)                 { m_strands = strands; }


} // namespace GeometryTier
} // namespace V


#endif

