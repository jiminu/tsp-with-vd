#ifndef _CONSTFORXVOLLJPOTENERGY_H_
#define _CONSTFORXVOLLJPOTENERGY_H_


#include "ForceField.h"
using namespace V::GeometryTier;


// Conversion table from Amber atom type to Atom code
// Note: 
// Atom types which do not appear in "AMBER_ATOM_TYPE"
// are substituted with UNK_ATOM.
const AtomCode ATOMCODE_FOR_AMBER_ATOM_TYPE[NUM_ATOM_TYPE_IN_AMBER] =
{
	C_ATOM             , C_ATOM             , UNK_ATOM /*CM_ATM*/, UNK_ATOM /*Cs_ATM*/, C_ATOM             , UNK_ATOM /*F_ATM*/ , H_ATOM, H_ATOM, 
	UNK_ATOM /*H2_ATM*/, UNK_ATOM /*H3_ATM*/,              H_ATOM, H_ATOM             , H_ATOM             , H_ATOM             , H_ATOM, H_ATOM, 
	H_ATOM             , UNK_ATOM /*HW_ATM*/, UNK_ATOM /*IP_ATM*/, UNK_ATOM /*K_ATM*/ , UNK_ATOM /*Li_ATM*/, N_ATOM             , N_ATOM, O_ATOM, 
	O_ATOM             , O_ATOM             , UNK_ATOM /*OS_ATM*/, UNK_ATOM /*OW_ATM*/, UNK_ATOM /*P_ATM*/ , UNK_ATOM /*Rb_ATM*/, S_ATOM, S_ATOM, UNK_ATOM
};

const char ONE_CHARACTER_OF_ATOMCODE[NUM_ATOM_TYPE_IN_AMBER] =
{
	'C'           , 'C'           , 'X' /*CM_ATM*/, 'X' /*Cs_ATM*/, 'C'           , 'X' /*F_ATM*/ , 'H', 'H', 
	'X' /*H2_ATM*/, 'X' /*H3_ATM*/, 'H'           , 'H'           , 'H'           , 'H'           , 'H', 'H', 
	'H'           , 'X' /*HW_ATM*/, 'X' /*IP_ATM*/, 'X' /*K_ATM*/ , 'X' /*Li_ATM*/, 'N'           , 'N', 'O', 
	'O'           , 'O'           , 'X' /*OS_ATM*/, 'X' /*OW_ATM*/, 'X' /*P_ATM*/ , 'X' /*Rb_ATM*/, 'S', 'S', 'X'
};

// test
const char AmberAtomTypeStrings[33][8] = {  "C_ATM",  "CA_ATM", "CM_ATM", "Cs_ATM", "CT_ATM", "F_ATM", " H_ATM", "H1_ATM", 
	                                        "H2_ATM", "H3_ATM", "H4_ATM", "H5_ATM", "HA_ATM", "HC_ATM","HO_ATM", "HP_ATM", 
	                                        "HS_ATM", "HW_ATM", "IP_ATM", "K_ATM",  "Li_ATM", "N_ATM", "N3_ATM", "O_ATM" , 
	                                        "O2_ATM", "OH_ATM", "OS_ATM", "OW_ATM", "P_ATM",  "Rb_ATM","S_ATM",  "SH_ATM", "UNK_ATM" };

#endif