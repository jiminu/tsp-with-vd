#ifndef _CONSTFORPOTENTIALENERGY_H
#define _CONSTFORPOTENTIALENERGY_H

#include "rg_Const.h"



namespace V {
namespace GeometryTier {



const rg_INT    VDW_TERM_ON   = 1<<0;
const rg_INT    ELEC_TERM_ON  = 1<<1;
const rg_INT    HBOND_TERM_ON = 1<<2;

const rg_INT    ID_VDW_TERM   = 0;
const rg_INT    ID_ELEC_TERM  = 1;
const rg_INT    ID_HBOND_TERM = 2;




} // namespace GeometryTier
} // namespace V


#endif

