#ifndef _CHEMICALPROPERTIESOFATOM_H
#define _CHEMICALPROPERTIESOFATOM_H

#include "rg_Const.h"
#include "ConstForMolecule.h"
#include "ForceField.h"
#include "PharmaFeature.h"
#include <string>
using namespace std;


namespace V {
namespace GeometryTier {



class ChemicalPropertiesOfAtom
{
private:

    RemoteIndicator          m_remoteIndicator;
    BranchDesignator         m_brangeDesignator;
    BranchDesignator         m_extraBrangeDesignator;
    AmberAtomTypes           m_atomTypeInAmber;
    
    // JKIM ADDED///////////////////////////////
    string                   m_SYBYLAtomType;
    ////////////////////////////////////////////


    rg_REAL                  m_charge;
    rg_REAL                  m_chargeInAmber;
    rg_REAL                  m_occupancy;
    rg_REAL                  m_tempFactor;
    rg_INT                   m_numOfConnectedHydrogen;
    rg_INT                   m_numOfAcceptableHydrogen;
    
    rg_dList<PharmaFeature*> m_pharmaFeatures;
    
    rg_FLAG                  m_isOnBackBone;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    ChemicalPropertiesOfAtom();
    ChemicalPropertiesOfAtom( const RemoteIndicator& remoteIndicator, const BranchDesignator& brangeDesignator, 
                              const AmberAtomTypes& atomTypeInAmber, const rg_REAL& charge, const rg_REAL& occupancy, 
                              const rg_REAL& tempFactor, const rg_INT& numOfConnectedHydrogen, const rg_INT& numOfAcceptableHydrogen );
    ChemicalPropertiesOfAtom( const ChemicalPropertiesOfAtom& chemicalProperties );
    ~ChemicalPropertiesOfAtom();

    //  GET FUNCTION
    inline RemoteIndicator     getRemoteIndicator() const { return m_remoteIndicator; };
    string                     getRemoteIndicatorInString() const;
    inline BranchDesignator    getBrangeDesignator() const { return m_brangeDesignator; };
    string                     getBrangeDesignatorInString() const;
    inline BranchDesignator    getExtraBrangeDesignator() const { return m_extraBrangeDesignator; };
    string                     getExtraBrangeDesignatorInString() const;
    inline AmberAtomTypes      getAtomTypeInAmber() const { return m_atomTypeInAmber; };
    inline rg_REAL             getCharge() const { return m_charge; };
    inline rg_REAL             getChargeInAmber() const { return m_chargeInAmber; };
    inline rg_REAL             getOccupancy() const { return m_occupancy; };
    inline rg_REAL             getTempFactor() const { return m_tempFactor; };
    inline rg_INT              getNumOfConnectedHydrogen() const { return m_numOfConnectedHydrogen; };
    inline rg_INT              getNumOfAcceptableHydrogen() const { return m_numOfAcceptableHydrogen; };
    rg_dList<PharmaFeature*>*  getListOfPharmaFeatures();

    inline string              getSYBYLAtomType() const { return m_SYBYLAtomType;}

    rg_FLAG                    isOnBackBone() const;

    //  SET FUNCTION
    void        setAllProperties( const RemoteIndicator& remoteIndicator, const BranchDesignator& brangeDesignator, const BranchDesignator& extraBrangeDesignator,
                                  const AmberAtomTypes& atomTypeInAmber, const rg_REAL& charge, const rg_REAL& occupancy, 
                                  const rg_REAL& tempFactor, const rg_INT& numOfConnectedHydrogen, const rg_INT& numOfAcceptableHydrogen, const rg_FLAG isOnBackBone );
    void        setRemoteIndicator( const RemoteIndicator& remoteIndicator );
    void        setBrangeDesignator( const BranchDesignator& brangeDesignator );
    void        setExtraBrangeDesignator( const BranchDesignator& extraBrangeDesignator );
    void        setAtomTypeInAmber( const AmberAtomTypes& atomTypeInAmber );
    void        setCharge( const rg_REAL& charge );
    void        setChargeInAmber( const rg_REAL& chargeInAmber );
    void        setOccupancy( const rg_REAL& occupancy );
    void        setTempFactor( const rg_REAL& tempFactor );
    void        setNumOfConnectedHydrogen( const rg_INT& numOfConnectedHydrogen );
    void        setNumOfAcceptableHydrogen( const rg_INT& numOfAcceptableHydrogen );
    
    void        setIsOnBackBone( const rg_FLAG isOnBackBone );

    void        addPharmaFeature( PharmaFeature* aPharmaFeature );


    inline void setSYBYLAtomType(const string& sybylAtomType) { m_SYBYLAtomType = sybylAtomType;}

    //  OPERATOR OVERLOADING
    ChemicalPropertiesOfAtom& operator =(const ChemicalPropertiesOfAtom& chemicalProperties);
};



} // namespace GeometryTier
} // namespace V


#endif

