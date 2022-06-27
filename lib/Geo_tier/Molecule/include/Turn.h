#ifndef _TURN_H
#define _TURN_H

#include "rg_Const.h"
#include "rg_dList.h"
#include <string>

using namespace std;



namespace V {
namespace GeometryTier {


class Residue;
class Chain;

class Turn
{
private:
    rg_INT               m_serial;
    string               m_turnID;
    rg_dList<Residue*>   m_residues;

    string               m_comments;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Turn();
    Turn( const rg_INT& serial );
    Turn( const rg_INT& serial, const string& turnID );
    Turn( const rg_INT& serial, const string& turnID, rg_dList<Residue*>* residues, const string& comments );
    Turn( const Turn& aTurn);
    ~Turn();


    //  GET FUNCTION
    rg_INT                 getSerial() const;
    string                 getTurnID() const;
    rg_dList<Residue*>*    getResidues();
    string                 getComments() const;

    Chain*                 getChain();

    rg_BOOL                isResidueInTurn( Residue* aResidue ) const;

    //  SET FUNCTION
    void       setSerial( const rg_INT& serial );
    void       setTurnID( const string& turnID );
    
    void       addResidue( Residue* aResidue );

    void       setComments( const string& comments );


    //  OPERATOR OVERLOADING
    Turn& operator =( const Turn& aTurn );

};


inline  rg_INT                                  V::GeometryTier::Turn::getSerial() const     { return m_serial; }
inline  string                                  V::GeometryTier::Turn::getTurnID() const     { return m_turnID; }
inline  rg_dList<V::GeometryTier::Residue*>* V::GeometryTier::Turn::getResidues()         { return &m_residues; }
inline  string                                  V::GeometryTier::Turn::getComments() const   { return m_comments; }


inline  void    V::GeometryTier::Turn::setSerial(const rg_INT& serial) { m_serial = serial; }
inline  void    V::GeometryTier::Turn::setTurnID(const string& turnID) { m_turnID = turnID; }
inline  void    V::GeometryTier::Turn::addResidue(V::GeometryTier::Residue* aResidue) { m_residues.addTail(aResidue); }
inline  void    V::GeometryTier::Turn::setComments(const string& comments) { m_comments = comments; }




} // namespace GeometryTier
} // namespace V


#endif

