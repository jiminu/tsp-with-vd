#ifndef _CONSTFORPHARMACOPHORE_H
#define _CONSTFORPHARMACOPHORE_H

#include "rg_Const.h"
#include <bitset>
#include <map>
using namespace std;


namespace V {
namespace GeometryTier {



const rg_INT    MAX_CHEM_ATOM_NUM = 200;

const rg_INT    NUM_OF_FEATRURE_TYPE = 6;

const rg_INT    UNK_FEATRURE_TYPE = 0;
const rg_INT     PI_FEATRURE_TYPE = 1;
const rg_INT    NI_FEATRURE_TYPE  = 2; 
const rg_INT    HBA_FEATRURE_TYPE = 3;
const rg_INT    HBD_FEATRURE_TYPE = 4;
const rg_INT     HY_FEATRURE_TYPE = 5;

// const rg_INT UNK_FEATRURE_TYPE_BIT  = 1<< UNK_FEATRURE_TYPE;
// const rg_INT  PI_FEATRURE_TYPE_BIT  = 1<<  PI_FEATRURE_TYPE;
// const rg_INT  NI_FEATRURE_TYPE_BIT  = 1<<  NI_FEATRURE_TYPE; 
// const rg_INT HBA_FEATRURE_TYPE_BIT  = 1<< HBA_FEATRURE_TYPE;
// const rg_INT HBD_FEATRURE_TYPE_BIT  = 1<< HBD_FEATRURE_TYPE;
// const rg_INT  HY_FEATRURE_TYPE_BIT  = 1<<  HY_FEATRURE_TYPE;



const rg_INT              UNK_CHEM_TYPE  = 0;

const rg_INT       POS_CHARGE_CHEM_TYPE  = 1;  // PI START
const rg_INT            AMINE_CHEM_TYPE  = 2;
const rg_INT          AMIDINE_CHEM_TYPE  = 3;
const rg_INT        GUANIDINE_CHEM_TYPE  = 4;  // PI END
                          
const rg_INT       NEG_CHARGE_CHEM_TYPE  = 5;  // NI START
const rg_INT         CARBOXYL_CHEM_TYPE  = 6;
const rg_INT           TRIFLU_CHEM_TYPE  = 7;
const rg_INT         SULFINIC_CHEM_TYPE  = 8;
const rg_INT         SULFONIC_CHEM_TYPE  = 9;
const rg_INT         SULFURIC_CHEM_TYPE  = 10;
const rg_INT       PHOSPHINIC_CHEM_TYPE  = 11;
const rg_INT       PHOSPHONIC_CHEM_TYPE  = 12;
const rg_INT PHOSPHORIC_ESTER_CHEM_TYPE  = 13;
const rg_INT        TETRAZOLE_CHEM_TYPE  = 14; // NI END
                          
const rg_INT       N_LONEPAIR_CHEM_TYPE  = 15; // HBA START
const rg_INT       O_LONEPAIR_CHEM_TYPE  = 16;
const rg_INT       S_LONEPAIR_CHEM_TYPE  = 17; // HBA END
                            
const rg_INT         HYDROXYL_CHEM_TYPE  = 18; // HBD START
const rg_INT            THIOL_CHEM_TYPE  = 19;
const rg_INT        ACETYLENE_CHEM_TYPE  = 20;
const rg_INT               NH_CHEM_TYPE  = 21; // HBD END

const rg_INT           C_RING_CHEM_TYPE  = 22; // HY START
const rg_INT          HALOGEN_CHEM_TYPE  = 23;
const rg_INT          C_GROUP_CHEM_TYPE  = 24; // HY END


// const rg_INT              UNK_CHEM_TYPE_BIT  = 1<<              UNK_CHEM_TYPE;
// 
// const rg_INT       POS_CHARGE_CHEM_TYPE_BIT  = 1<<       POS_CHARGE_CHEM_TYPE; // PI START
// const rg_INT            AMINE_CHEM_TYPE_BIT  = 1<<            AMINE_CHEM_TYPE;
// const rg_INT          AMIDINE_CHEM_TYPE_BIT  = 1<<          AMIDINE_CHEM_TYPE;
// const rg_INT      GUANIDINIUM_CHEM_TYPE_BIT  = 1<<      GUANIDINIUM_CHEM_TYPE; // PI END
//                                        
// const rg_INT       NEG_CHARGE_CHEM_TYPE_BIT  = 1<<       NEG_CHARGE_CHEM_TYPE; // NI START
// const rg_INT         CARBOXYL_CHEM_TYPE_BIT  = 1<<         CARBOXYL_CHEM_TYPE;
// const rg_INT            TFMSA_CHEM_TYPE_BIT  = 1<<            TFMSA_CHEM_TYPE;
// const rg_INT         SULFINIC_CHEM_TYPE_BIT  = 1<<         SULFINIC_CHEM_TYPE;
// const rg_INT         SULFONIC_CHEM_TYPE_BIT  = 1<<         SULFONIC_CHEM_TYPE;
// const rg_INT         SULFURIC_CHEM_TYPE_BIT  = 1<<         SULFURIC_CHEM_TYPE;
// const rg_INT       PHOSPHINIC_CHEM_TYPE_BIT  = 1<<       PHOSPHINIC_CHEM_TYPE;
// const rg_INT       PHOSPHONIC_CHEM_TYPE_BIT  = 1<<       PHOSPHONIC_CHEM_TYPE;
// const rg_INT PHOSPHORIC_ESTER_CHEM_TYPE_BIT  = 1<< PHOSPHORIC_ESTER_CHEM_TYPE;
// const rg_INT        TETRAZOLE_CHEM_TYPE_BIT  = 1<<        TETRAZOLE_CHEM_TYPE; // NI END
//                                        
// const rg_INT       N_LONEPAIR_CHEM_TYPE_BIT  = 1<<       N_LONEPAIR_CHEM_TYPE; // HBA START
// const rg_INT       O_LONEPAIR_CHEM_TYPE_BIT  = 1<<       O_LONEPAIR_CHEM_TYPE;
// const rg_INT       S_LONEPAIR_CHEM_TYPE_BIT  = 1<<       S_LONEPAIR_CHEM_TYPE; // HBA END
//                                            
// const rg_INT         HYDROXYL_CHEM_TYPE_BIT  = 1<<         HYDROXYL_CHEM_TYPE; // HBD START
// const rg_INT            THIOL_CHEM_TYPE_BIT  = 1<<            THIOL_CHEM_TYPE;
// const rg_INT        ACETYLENE_CHEM_TYPE_BIT  = 1<<        ACETYLENE_CHEM_TYPE;
// const rg_INT               NH_CHEM_TYPE_BIT  = 1<<               NH_CHEM_TYPE; // HBD END
// 
// const rg_INT        FIVE_RING_CHEM_TYPE_BIT  = 1<<        FIVE_RING_CHEM_TYPE; // HY START
// const rg_INT         SIX_RING_CHEM_TYPE_BIT  = 1<<         SIX_RING_CHEM_TYPE;
// const rg_INT       TERT_BUTYL_CHEM_TYPE_BIT  = 1<<       TERT_BUTYL_CHEM_TYPE;
// const rg_INT          HALOGEN_CHEM_TYPE_BIT  = 1<<          HALOGEN_CHEM_TYPE; // HY END

const float PHARMA_DEFAULT_COLOR[3] = {0.5f, 0.5f, 0.5f};
const float PHARMA_COMPLEX_COLOR[3] = {0.5f, 0.5f, 0.5f};
const float PHARMA_FEATURE_COLOR[NUM_OF_FEATRURE_TYPE][3] =                
                   {{1.0f,  1.0f,  1.0f},
                    {0.9f,  0.4f,  0.6f},
                    {0.4f,  0.6f,  0.9f},
                    {0.4f,  0.3f,  0.6f},
                    {1.0f,  0.6f,  0.1f},
                    {0.6f,  0.6f,  0.4f}};
                    
                                             




} // namespace GeometryTier
} // namespace V


#endif

