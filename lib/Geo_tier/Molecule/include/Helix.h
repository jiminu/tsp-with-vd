#ifndef _Helix_H
#define _Helix_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "ConstForMolecule.h"
#include "Residue.h"
#include <string>

using namespace std;



namespace V {
namespace GeometryTier {



class Chain;

class Helix
{
private:
    rg_INT               m_serial;
    string               m_helixID;
    rg_dList<Residue*>   m_residues;

    HelixClass           m_helixClass;
    string               m_comments;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    Helix();
    Helix( const rg_INT& serial, const string& helixID, const HelixClass& helixClass );
    Helix( const rg_INT& serial, const string& helixID, const rg_INT& helixClass );
    Helix( const rg_INT& serial, const string& helixID, const rg_INT& helixClass, rg_dList<Residue*>* residues, const string& comments );
    Helix( const Helix& aHelix);
    ~Helix();


    //  GET FUNCTION
    rg_INT                 getSerial() const;
    string                 getHelixID() const;
    rg_dList<Residue*>*    getResidues();
    HelixClass             getHelixClass() const;
    string                 getComments() const;

    Chain*                 getChain();
    rg_BOOL                isResidueInHelix( Residue* aResidue ) const;


    //  SET FUNCTION
    void       setSerial( const rg_INT& serial );
    void       setHelixID( const string& helixID );
    
    void       addResidue( Residue* aResidue );
    void       setResidues( rg_dList<Residue*>* residues );

    void       setHelixClass( const HelixClass& helixClass );
    void       setComments( const string& comments );


    //  OPERATOR OVERLOADING
    Helix& operator =( const Helix& aHelix );

};


inline  rg_INT                                  V::GeometryTier::Helix::getSerial() const        { return m_serial; }
inline  string                                  V::GeometryTier::Helix::getHelixID() const       { return m_helixID; }
inline  rg_dList<V::GeometryTier::Residue*>* V::GeometryTier::Helix::getResidues()            { return &m_residues; }
inline  V::GeometryTier::HelixClass          V::GeometryTier::Helix::getHelixClass() const    { return m_helixClass; }
inline  string                                  V::GeometryTier::Helix::getComments() const      { return m_comments; }


inline  void    V::GeometryTier::Helix::setSerial(const rg_INT& serial)              { m_serial = serial; }
inline  void    V::GeometryTier::Helix::setHelixID(const string& helixID)            { m_helixID = helixID; }
inline  void    V::GeometryTier::Helix::addResidue(Residue* aResidue)                { m_residues.addTail(aResidue); }
inline  void    V::GeometryTier::Helix::setResidues(rg_dList<Residue*>* residues)    { m_residues.append(*residues); }
inline  void    V::GeometryTier::Helix::setHelixClass(const HelixClass& helixClass)  { m_helixClass = helixClass; }
inline  void    V::GeometryTier::Helix::setComments(const string& comments)          { m_comments = comments; }



} // namespace GeometryTier
} // namespace V


#endif

