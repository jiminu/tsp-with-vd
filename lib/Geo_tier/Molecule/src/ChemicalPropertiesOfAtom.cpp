#include "ChemicalPropertiesOfAtom.h"
using namespace V::GeometryTier;


ChemicalPropertiesOfAtom::ChemicalPropertiesOfAtom()
: m_remoteIndicator(UNK_REMOTE), m_brangeDesignator(UNK_BRANCH), m_extraBrangeDesignator(UNK_BRANCH), m_atomTypeInAmber(UNK_ATM), m_charge(.0), m_chargeInAmber(.0), m_occupancy(.0), m_tempFactor(.0), m_numOfConnectedHydrogen(0), m_numOfAcceptableHydrogen(0), m_isOnBackBone(rg_FALSE)
{
    m_SYBYLAtomType           = "";
}



ChemicalPropertiesOfAtom::ChemicalPropertiesOfAtom( const RemoteIndicator& remoteIndicator, const BranchDesignator& brangeDesignator, 
                                                    const AmberAtomTypes& atomTypeInAmber, const rg_REAL& charge, const rg_REAL& occupancy,       
                                                    const rg_REAL& tempFactor, const rg_INT& numOfConnectedHydrogen, const rg_INT& numOfAcceptableHydrogen )
{
    m_remoteIndicator         = remoteIndicator;
    m_brangeDesignator        = brangeDesignator;
    m_atomTypeInAmber         = atomTypeInAmber;
    m_charge                  = charge;
    m_occupancy               = occupancy;
    m_tempFactor              = tempFactor;
    m_numOfConnectedHydrogen  = numOfConnectedHydrogen;
    m_numOfAcceptableHydrogen = numOfAcceptableHydrogen;
    
    m_extraBrangeDesignator   = UNK_BRANCH;

    m_SYBYLAtomType           = "";
}



ChemicalPropertiesOfAtom::ChemicalPropertiesOfAtom( const ChemicalPropertiesOfAtom& chemicalProperties )
{
    m_remoteIndicator         = chemicalProperties.m_remoteIndicator;
    m_brangeDesignator        = chemicalProperties.m_brangeDesignator;
    m_extraBrangeDesignator   = chemicalProperties.m_extraBrangeDesignator;
    m_atomTypeInAmber         = chemicalProperties.m_atomTypeInAmber;
    m_charge                  = chemicalProperties.m_charge;
    m_chargeInAmber           = chemicalProperties.m_chargeInAmber;
    m_occupancy               = chemicalProperties.m_occupancy;
    m_tempFactor              = chemicalProperties.m_tempFactor;
    m_numOfConnectedHydrogen  = chemicalProperties.m_numOfConnectedHydrogen;
    m_numOfAcceptableHydrogen = chemicalProperties.m_numOfAcceptableHydrogen;
    m_isOnBackBone            = chemicalProperties.m_isOnBackBone;
    m_pharmaFeatures          = chemicalProperties.m_pharmaFeatures;

    m_SYBYLAtomType           = chemicalProperties.m_SYBYLAtomType;
}



ChemicalPropertiesOfAtom::~ChemicalPropertiesOfAtom()
{
}


// RemoteIndicator ChemicalPropertiesOfAtom::getRemoteIndicator() const
// {
//     return m_remoteIndicator;
// }
// 
// 
// BranchDesignator ChemicalPropertiesOfAtom::getBrangeDesignator() const
// {
//     return m_brangeDesignator;
// }
// 
// 
// AmberAtomTypes ChemicalPropertiesOfAtom::getAtomTypeInAmber() const
// {
//     return m_atomTypeInAmber;
// }
// 
// 
// rg_REAL ChemicalPropertiesOfAtom::getCharge() const
// {
//     return m_charge;
// }
// 
// 
// rg_REAL ChemicalPropertiesOfAtom::getOccupancy() const
// {
//     return m_occupancy;
// }
// 
// 
// rg_REAL ChemicalPropertiesOfAtom::getTempFactor() const
// {
//     return m_tempFactor;
// }
// 
// 
// rg_INT ChemicalPropertiesOfAtom::getNumOfConnectedHydrogen() const
// {
//     return m_numOfConnectedHydrogen;
// }
// 
// 
// rg_INT ChemicalPropertiesOfAtom::getNumOfAcceptableHydrogen() const
// {
//     return m_numOfAcceptableHydrogen;
// }



string ChemicalPropertiesOfAtom::getRemoteIndicatorInString() const
{
    string strRemoteIndicator;

    switch( m_remoteIndicator ) {
        case ALPHA_REMOTE :
            strRemoteIndicator = "A";
    	    break;
        case BETA_REMOTE :
            strRemoteIndicator = "B";   
    	    break;
        case GAMMA_REMOTE :
            strRemoteIndicator = "G";
            break;
        case DELTA_REMOTE :
            strRemoteIndicator = "D";   
            break;
        case EPSILON_REMOTE :
            strRemoteIndicator = "E";
    	    break;
        case ZETA_REMOTE :
            strRemoteIndicator = "Z";
    	    break;
        case ETA_REMOTE :
            strRemoteIndicator = "H";   
            break;
        case TERMINATE_REMOTE :
            strRemoteIndicator = "X";   
            break;
        default:
            strRemoteIndicator = " ";
            break;
    }

    return strRemoteIndicator;
}



string ChemicalPropertiesOfAtom::getBrangeDesignatorInString() const
{
    string strBranchDesignator;

    switch( m_brangeDesignator ) {
        case FIRST_BRANCH :
            strBranchDesignator = "1";
    	    break;
        case SECOND_BRANCH :
            strBranchDesignator = "2";   
    	    break;
        case THIRD_BRANCH :
            strBranchDesignator = "3";
            break;
        case FOURTH_BRANCH :
            strBranchDesignator = "T";   
            break;
        default:
            strBranchDesignator = " ";
            break;
    }

    return strBranchDesignator;
}

string ChemicalPropertiesOfAtom::getExtraBrangeDesignatorInString() const
{
    string strExtraBranchDesignator;

    switch( m_extraBrangeDesignator ) {
        case FIRST_BRANCH :
            strExtraBranchDesignator = "1";
    	    break;
        case SECOND_BRANCH :
            strExtraBranchDesignator = "2";   
    	    break;
        case THIRD_BRANCH :
            strExtraBranchDesignator = "3";
            break;
        case FOURTH_BRANCH :
            strExtraBranchDesignator = "4";   
            break;
        default:
            strExtraBranchDesignator = " ";
            break;
    }

    return strExtraBranchDesignator;    
}



rg_dList<PharmaFeature*>* ChemicalPropertiesOfAtom::getListOfPharmaFeatures()
{
    return &m_pharmaFeatures;
}



rg_FLAG ChemicalPropertiesOfAtom::isOnBackBone() const
{
    return m_isOnBackBone;
}



void ChemicalPropertiesOfAtom::setAllProperties( const RemoteIndicator& remoteIndicator, const BranchDesignator& brangeDesignator, const BranchDesignator& extraBrangeDesignator, const AmberAtomTypes& atomTypeInAmber, const rg_REAL& charge,  const rg_REAL& occupancy,  const rg_REAL& tempFactor, const rg_INT& numOfConnectedHydrogen, const rg_INT& numOfAcceptableHydrogen, const rg_FLAG isOnBackBone )
{
    m_remoteIndicator         = remoteIndicator;
    m_brangeDesignator        = brangeDesignator;
    m_extraBrangeDesignator   = extraBrangeDesignator;
    m_atomTypeInAmber         = atomTypeInAmber;
    m_charge                  = charge;
    m_occupancy               = occupancy;
    m_tempFactor              = tempFactor;
    m_numOfConnectedHydrogen  = numOfConnectedHydrogen;
    m_numOfAcceptableHydrogen = numOfAcceptableHydrogen;
    m_isOnBackBone            = isOnBackBone;
}



void ChemicalPropertiesOfAtom::setRemoteIndicator( const RemoteIndicator& remoteIndicator )
{
    m_remoteIndicator = remoteIndicator;    
}



void ChemicalPropertiesOfAtom::setBrangeDesignator( const BranchDesignator& brangeDesignator )
{
    m_brangeDesignator = brangeDesignator;
}



void ChemicalPropertiesOfAtom::setExtraBrangeDesignator( const BranchDesignator& extraBrangeDesignator )
{
    m_extraBrangeDesignator = extraBrangeDesignator;
}



void ChemicalPropertiesOfAtom::setAtomTypeInAmber( const AmberAtomTypes& atomTypeInAmber )
{
    m_atomTypeInAmber = atomTypeInAmber;    
}



void ChemicalPropertiesOfAtom::setCharge( const rg_REAL& charge )
{
    m_charge = charge;
}



void ChemicalPropertiesOfAtom::setChargeInAmber( const rg_REAL& chargeInAmber )
{
    m_chargeInAmber = chargeInAmber;
}



void ChemicalPropertiesOfAtom::setOccupancy( const rg_REAL& occupancy )
{
    m_occupancy = occupancy;    
}



void ChemicalPropertiesOfAtom::setTempFactor( const rg_REAL& tempFactor )
{
    m_tempFactor = tempFactor;    
}



void ChemicalPropertiesOfAtom::setNumOfConnectedHydrogen( const rg_INT& numOfConnectedHydrogen )
{
    m_numOfConnectedHydrogen = numOfConnectedHydrogen;
}



void ChemicalPropertiesOfAtom::setNumOfAcceptableHydrogen( const rg_INT& numOfAcceptableHydrogen )
{
    m_numOfAcceptableHydrogen = numOfAcceptableHydrogen;
}



void ChemicalPropertiesOfAtom::setIsOnBackBone( const rg_FLAG isOnBackBone )
{
    m_isOnBackBone = isOnBackBone;
}



void ChemicalPropertiesOfAtom::addPharmaFeature( PharmaFeature* aPharmaFeature )
{
    m_pharmaFeatures.addTail( aPharmaFeature );
}



ChemicalPropertiesOfAtom& ChemicalPropertiesOfAtom::operator=( const ChemicalPropertiesOfAtom& chemicalProperties )
{
    if( this == &chemicalProperties )
        return *this;

    m_remoteIndicator         = chemicalProperties.m_remoteIndicator;
    m_brangeDesignator        = chemicalProperties.m_brangeDesignator;
    m_extraBrangeDesignator   = chemicalProperties.m_extraBrangeDesignator;
    m_atomTypeInAmber         = chemicalProperties.m_atomTypeInAmber;
    m_charge                  = chemicalProperties.m_charge;
    m_chargeInAmber           = chemicalProperties.m_chargeInAmber;
    m_occupancy               = chemicalProperties.m_occupancy;
    m_tempFactor              = chemicalProperties.m_tempFactor;
    m_numOfConnectedHydrogen  = chemicalProperties.m_numOfConnectedHydrogen;
    m_numOfAcceptableHydrogen = chemicalProperties.m_numOfAcceptableHydrogen;
    m_isOnBackBone            = chemicalProperties.m_isOnBackBone;
    m_pharmaFeatures          = chemicalProperties.m_pharmaFeatures;

    m_SYBYLAtomType           = chemicalProperties.m_SYBYLAtomType;


    return *this;
}
