#include "RotamerSetOfResidue.h"
#include "FunctionsForRotamerLibrary.h"
#include "Rotamer.h"
#include "Residue.h"

RotamerSetOfResidue::RotamerSetOfResidue()
{
	set(rg_NULL,
		UNKNOWN_ROT_INDEX,
		UNKNOWN_ROT_INDEX,
		UNKNOWN_ROT_INDEX,
        rg_FALSE);
    //set(rg_NULL,
    //    UNKNOWN_ROT_INDEX,
    //    UNKNOWN_ROT_INDEX,
    //    UNKNOWN_ROT_INDEX);
}

RotamerSetOfResidue::RotamerSetOfResidue(const RotamerSetOfResidue& rotamerSetOfResidue)
{
	set(rotamerSetOfResidue.m_rotamer, 
		rotamerSetOfResidue.m_startRotLibID, 
		rotamerSetOfResidue.m_endRotLibLID, 
		rotamerSetOfResidue.m_assignedRotLibID,
        rotamerSetOfResidue.m_bFixSidechainWithCurrentlyAssignedRotamer);
    //set(rotamerSetOfResidue.m_rotamer, 
    //    rotamerSetOfResidue.m_startRotLibID, 
    //    rotamerSetOfResidue.m_endRotLibLID, 
    //    rotamerSetOfResidue.m_assignedRotLibID);
}

RotamerSetOfResidue::RotamerSetOfResidue(Rotamer* rotamer, const rg_INDEX& startRotLibID, const rg_INDEX& endRotLibID)
{
	set(rotamer, 
		startRotLibID, 
		endRotLibID, 
		UNKNOWN_ROT_INDEX,
        rg_FALSE);
    //set(rotamer, 
    //    startRotLibID, 
    //    endRotLibID, 
    //    UNKNOWN_ROT_INDEX);
}

RotamerSetOfResidue::~RotamerSetOfResidue()
{
}

Rotamer* RotamerSetOfResidue::getRotamerCorrToRotLibID( const rg_INDEX& rotLibID )
{
	if (rotLibID != UNKNOWN_ROT_INDEX && rotLibID != FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX)
	{
		updateRotamer( rotLibID );
	}	
	return m_rotamer;
}

void RotamerSetOfResidue::removeRotamersWithProbabilityLessThan(const rg_REAL& thresholdOfProbability)
{
	rg_INDEX i = UNKNOWN_ROT_INDEX;
	for (i = m_startRotLibID;i <= m_endRotLibLID;i++)
	{
		rg_REAL probability = -1;
		TypeOfRotamerLibrary type = FunctionsForRotamerLibrary::getCurrentlyUsedRotamerLibType();
		if(type == BACKBONE_INDEPENDENT_LIB_DUNBRACK2002)
			probability = FunctionsForRotamerLibrary::BBINDEP_ROTAMER_LIB[ i ].probability;
		else if(type == BACKBONE_DEPENDENT_LIB_DUNBRACK2002)
			probability = FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB[ i ].probability;
		else if(type == BACKBONE_DEPENDENT_LIB_DUNBRACK2010)
			probability = FunctionsForRotamerLibrary::BBDEP_ROTAMER_LIB_2010[ i ].probability;

		if(probability < thresholdOfProbability)
		{
			if(m_startRotLibID <= (i - 1) && (i - 1) < m_endRotLibLID)
			{
				m_endRotLibLID = i - 1;
			}
			// i == m_startRotLibID 
			else
			{
				m_startRotLibID = m_endRotLibLID = FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX;
			}
			break;
		}
	}
}

void RotamerSetOfResidue::updateRotamer(const rg_INT& rotLibID)
{
	//FunctionsForRotamerLibrary::computeAtomCenterUsingDihedralAngle2(*m_rotamer, rotLibID);
	m_rotamer->updateSidechainAtomCoordinatesByApplyingRotLibID(rotLibID);
}

void RotamerSetOfResidue::updateRotamer(const rg_REAL* dihedralAngles)
{
	//FunctionsForRotamerLibrary::computeAtomCenterUsingDihedralAngle2(*m_rotamer, dihedralAngles);
	m_rotamer->updateSidechainAtomCoordinatesByApplyingDihedralAngles( dihedralAngles );
}

void RotamerSetOfResidue::computeEnergyWithOtherRotamerSetOfResidue( const RotamerSetOfResidue& otherRotamerSetOfResidue )
{

}

void RotamerSetOfResidue::computeEnergyWithBackbone( Backbone& backbone )
{

}

void RotamerSetOfResidue::computeMinimumEnclosingSphere(EnclosingSphereOfSpheres& MES_U, const rg_BOOL& bBackboneAtomsIncluded /*= rg_FALSE*/)
{
	rg_dList<Atom*> atomsOfUROAR;

	//if(m_assignedRotLibID != UNKNOWN_ROT_INDEX && m_assignedRotLibID != FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX)
	//{
	//	updateRotamer( m_assignedRotLibID );
	//	rg_dList<Atom*> atomsOnSidechain;
	//	m_rotamer->getAtomsOnSidechain( atomsOnSidechain );
	//	atomsOfUROAR.append(atomsOnSidechain);
	//}
	//else
	//{
	Residue* residue = m_rotamer->getResidue();
	ResidueCode code = residue->getResidueCode();

    switch (code)
    {
    case ALA_AMINO_RESIDUE:
        {
            Atom* betaCarbon = residue->getBetaCarbonInSideChainOfAminoResidue();
            atomsOfUROAR.add( betaCarbon );
        }
    	break;
    case GLY_AMINO_RESIDUE:
        {
            Atom* alphaCarbon = residue->getAlphaCarbonOfAminoResidue();
            atomsOfUROAR.add( alphaCarbon );
        }
        break;
    default:
        if(m_assignedRotLibID != UNKNOWN_ROT_INDEX)
        {
            if(m_assignedRotLibID != FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX)
                updateRotamer( m_assignedRotLibID );
            rg_dList<Atom*> atomsOnSidechain;
            //m_rotamer->getAtomsOnSidechain( atomsOnSidechain );
            residue->getAtomsOnSidechain(atomsOnSidechain);
            atomsOfUROAR.append(atomsOnSidechain);
        }
        else
        {
            rg_INDEX i;
            for (i = m_startRotLibID;i <= m_endRotLibLID;i++)
            {
                updateRotamer( i );
                rg_dList<Atom*> atomsOnSidechain;
                //m_rotamer->getAtomsOnSidechain( atomsOnSidechain );
                residue->getAtomsOnSidechain(atomsOnSidechain);
                atomsOfUROAR.append(atomsOnSidechain);
            }
        }
        break;
    }

    if(bBackboneAtomsIncluded)
    {
        switch (code)
        {
        case GLY_AMINO_RESIDUE:
            {
                atomsOfUROAR.removeAll();

                rg_dList<Atom*> atomsOnBackbone;
                residue->getAtomsOnBackbone(atomsOnBackbone);
                atomsOfUROAR.append(atomsOnBackbone);
            }
            break;
        default:
            {
                rg_dList<Atom*> atomsOnBackbone;
                residue->getAtomsOnBackbone(atomsOnBackbone);
                atomsOfUROAR.append(atomsOnBackbone);
            }
            break;
        }
    }

	//if (code != ALA_AMINO_RESIDUE && code != GLY_AMINO_RESIDUE)
	//{
 //       //if(m_bFixSidechainWithCurrentlyAssignedRotamer)
 //       if(m_assignedRotLibID != UNKNOWN_ROT_INDEX)
 //       {
 //           if(m_assignedRotLibID != FIX_CURRENT_CONFORMATION_WITHOUT_ROT_INDEX)
 //               updateRotamer( m_assignedRotLibID );
 //           rg_dList<Atom*> atomsOnSidechain;
 //           //m_rotamer->getAtomsOnSidechain( atomsOnSidechain );
 //           residue->getAtomsOnSidechain(atomsOnSidechain);
 //           atomsOfUROAR.append(atomsOnSidechain);
 //       }
 //       else
 //       {
 //           rg_INDEX i;
 //           for (i = m_startRotLibID;i <= m_endRotLibLID;i++)
 //           {
 //               updateRotamer( i );
 //               rg_dList<Atom*> atomsOnSidechain;
 //               //m_rotamer->getAtomsOnSidechain( atomsOnSidechain );
 //               residue->getAtomsOnSidechain(atomsOnSidechain);
 //               atomsOfUROAR.append(atomsOnSidechain);
 //           }
 //       }
	//}
	//else if(code == ALA_AMINO_RESIDUE)
	//{
	//	Atom* betaCarbon = residue->getBetaCarbonInSideChainOfAminoResidue();
	//	atomsOfUROAR.add( betaCarbon );
	//}
	//else if(code == GLY_AMINO_RESIDUE)
	//{
	//	Atom* alphaCarbon = residue->getAlphaCarbonOfAminoResidue();
	//	atomsOfUROAR.add( alphaCarbon );
	//}
	//else{}

	MES_U.setSpheres(atomsOfUROAR);
	MES_U.computeEnclosingSphere();
}

RotamerSetOfResidue& RotamerSetOfResidue::operator=(const RotamerSetOfResidue& rotamerSetOfResidue)
{
	if(this == &rotamerSetOfResidue)
		return *this;

	set(rotamerSetOfResidue.m_rotamer, 
		rotamerSetOfResidue.m_startRotLibID, 
		rotamerSetOfResidue.m_endRotLibLID, 
		rotamerSetOfResidue.m_assignedRotLibID,
        rotamerSetOfResidue.m_bFixSidechainWithCurrentlyAssignedRotamer);

    //set(rotamerSetOfResidue.m_rotamer, 
    //    rotamerSetOfResidue.m_startRotLibID, 
    //    rotamerSetOfResidue.m_endRotLibLID, 
    //    rotamerSetOfResidue.m_assignedRotLibID);

	return *this;
}