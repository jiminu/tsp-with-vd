#ifndef _ROTAMERSETOFRESIDUE_H_
#define _ROTAMERSETOFRESIDUE_H_

#include "Backbone.h"
#include "EnclosingSphereOfSpheres.h"
#include "ConstForRotamerLibrary.h"
class Rotamer;

// This class implicitly represents a set of rotamers ( in a rotamer library ) corresponding to a particular residue.
// Each rotamer can be obtained by applying the dihedral angle set to a residue (a sidechain).
// The instance of this class is named "rotamer set" or "UROAR (Union of Rotamer On A Residue".

// Refer to the following article for the details of the UROAR and related concepts.
// Joonghyun Ryu and Deok-Soo Kim, Protein Structure Optimization by Side-chain Positioning via Beta-complex, 
// Journal of Global Optimization, (DOI: 10.1007/s10898-012-9886-3), 2012.

class RotamerSetOfResidue
{
private:
	Rotamer*  m_rotamer         ;  // a rotamer corresponding to a given rotamer library index
	rg_INDEX  m_startRotLibID   ;  // start index of the rotamer library
	rg_INDEX  m_endRotLibLID    ;  // end   index of the rotamer library
	rg_INDEX  m_assignedRotLibID;  //       index of assigned (fixed) rotamer which could represent the native conformation (corresponding to FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX).
    rg_BOOL   m_bFixSidechainWithCurrentlyAssignedRotamer
                                ;  // a flag for determining whether the sidechain conformation is fixed or not

public:
	RotamerSetOfResidue();
	RotamerSetOfResidue(const RotamerSetOfResidue& rotamerSetOfResidue);
	RotamerSetOfResidue(Rotamer* rotamer, const rg_INDEX& startRotLibID, const rg_INDEX& endRotLibID);
	~RotamerSetOfResidue();

	inline Rotamer* getRotamer() { return m_rotamer; }
	Rotamer* getRotamerCorrToRotLibID( const rg_INDEX& rotLibID );

	inline rg_INDEX getStartRotLibID() const 
    { 
        return m_startRotLibID; 
    }

	inline rg_INDEX getRotLibIDWithHighestProbability() const 
    { 
        return m_startRotLibID; 
    }

	inline void getRotLibIDs(rg_INDEX& startRotLibID, rg_INDEX& endRotLibID) const
	{
		startRotLibID = m_startRotLibID;
		endRotLibID   = m_endRotLibLID;
	}

	inline rg_INDEX getAssignedRotLibID() const
	{
		return m_assignedRotLibID;
	}

    inline rg_BOOL isSidechainConformationFixed() const 
    {
        return m_bFixSidechainWithCurrentlyAssignedRotamer;
    }

	inline rg_BOOL isRotamerAssigned() const
	{
		//if(m_assignedRotLibID != UNKNOWN_ROT_INDEX && m_assignedRotLibID != FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX)
        if(m_assignedRotLibID != UNKNOWN_ROT_INDEX)
			return rg_TRUE;
		else
			return rg_FALSE;
	}

	inline rg_INT getNumRotamers() const
	{
		rg_INDEX startRotLibID = UNKNOWN_ROT_INDEX;
		rg_INDEX endRotLibID   = UNKNOWN_ROT_INDEX;
		getRotLibIDs(startRotLibID, endRotLibID);

		rg_INT numRotamers = (endRotLibID - startRotLibID) + 1;
		return numRotamers;
	}

	inline void  set(Rotamer* rotamer, 
		             const rg_INDEX& startRotLibID = UNKNOWN_ROT_INDEX, 
		             const rg_INDEX& endRotLibID   = UNKNOWN_ROT_INDEX, 
		             const rg_INDEX& assignedRotLibID = UNKNOWN_ROT_INDEX,
                     const rg_BOOL&  bFixSidechainWithCurrentRotamer = rg_FALSE)
	{
		m_rotamer       = rotamer;
		m_startRotLibID = startRotLibID;
		m_endRotLibLID  = endRotLibID;
		m_assignedRotLibID = assignedRotLibID;
        m_bFixSidechainWithCurrentlyAssignedRotamer = bFixSidechainWithCurrentRotamer;
	}

	inline void setRotLibIDs(const rg_INDEX& startRotLibID, const rg_INDEX& endRotLibID)
	{
		m_startRotLibID = startRotLibID;
		m_endRotLibLID  = endRotLibID  ;
	}

	inline void setAssignedRotLibID(const rg_INDEX& assignedRotLibID)
	{
		m_assignedRotLibID = assignedRotLibID;
	}

    inline void fixSidechainConformationWithRotamerID(const rg_INDEX& rotamerID)
    {
        m_bFixSidechainWithCurrentlyAssignedRotamer = rg_TRUE;
        m_assignedRotLibID = rotamerID;
    }

    inline void releaseSidechainConformation()
    {
        m_bFixSidechainWithCurrentlyAssignedRotamer = rg_FALSE;
        m_assignedRotLibID = UNKNOWN_ROT_INDEX;
    }

	inline void resetAssignedRotLibID()
	{
        if(! m_bFixSidechainWithCurrentlyAssignedRotamer)
		    m_assignedRotLibID = UNKNOWN_ROT_INDEX;
	}

	void  removeRotamersWithProbabilityLessThan(const rg_REAL& thresholdOfProbability);

	void  updateRotamer(const rg_INT& rotLibID);
	void  updateRotamer(const rg_REAL* dihedralAngles);

	void computeEnergyWithOtherRotamerSetOfResidue( const RotamerSetOfResidue& otherRotamerSetOfResidue );
	void computeEnergyWithBackbone( Backbone& backbone );
	void computeMinimumEnclosingSphere(V::GeometryTier::EnclosingSphereOfSpheres& MES_U, const rg_BOOL& bBackboneAtomsIncluded = rg_FALSE);

	RotamerSetOfResidue& operator=(const RotamerSetOfResidue& rotamerSetOfResidue);

};

#endif