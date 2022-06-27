#ifndef _FUNCTIONSFORROTAMERLIBRARY_
#define _FUNCTIONSFORROTAMERLIBRARY_

#include "rg_Point3D.h"
#include "Backbone.h"
#include "ManagerOfRotamers_Assigned_At_Residues.h"
#include "ManagerOfRotamerSetsForSCP.h"
#include "BackboneDependentRotamerLibrary.h"
#include "BackboneIndependentRotamerLibrary.h"

typedef rg_Point3D Vector3D;

class V::GeometryTier::Atom;
class V::GeometryTier::Residue;
class RotamerSetOfResidue;

// For detailed description, refer to the included "ConstForRotamerLibrary.h" and the following article
// Dunbrack Jr., R.L., Cohen, F.E.: Bayesian statistical analysis of protein side-chain rotamer preferences.
// Protein Science 6(8), 1661{1681 (1997)

class FunctionsForRotamerLibrary
{
public:
	static BBDepRotamer BBDEP_ROTAMER_LIB[SIZE_OF_BBDEP_ROTAMER_LIB];
	static BBDepRotamer BBDEP_ROTAMER_LIB_2010[SIZE_OF_BBDEP_ROTAMER_LIB_2010];
    //static BBDepRotamer* BBDEP_ROTAMER_LIB;
	static BBIndRotamer BBINDEP_ROTAMER_LIB[SIZE_OF_BBINDEP_ROTAMER_LIB];
	
private:
	static TypeOfRotamerLibrary m_currentlyUsedRotamerLibType                         ;
    static rg_BOOL              m_bBBINDLibLoaded                          ;
    static rg_BOOL              m_bBBDEPLibLoaded                          ;
    static rg_BOOL  getAngleBinIndex_old(const rg_REAL& angle, rg_INDEX& index);
	static rg_BOOL  getAngleBinIndex(const rg_REAL& angle, rg_INDEX& index);
	static rg_INDEX getIndexOfResCode(const ResidueCode& code);	

public:
	// functions for backbone-dependent library
	static void  getRotamerSetsOfResidueSetFromBBDepRotLib(Backbone             & backbone,
		                                                   Rotamer*               sortedAssignedRotamersOfProtein, 
														   RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
		                                                   ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP   );

	static void  getRotamerSetsOfResidueSetFromBBDepRotLib_Dunbrack2010(Backbone             & backbone,
		                                                                Rotamer*               sortedAssignedRotamersOfProtein, 
														                RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
		                                                                ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP   );

	static rg_BOOL getBBDepLibIndexOfResidue(const ResidueCode& code, 
		                                     rg_INDEX& start, 
											 rg_INDEX& end);

	static rg_BOOL getBBDepLibIndexOfResidue(const ResidueCode& code, 
		                                     const rg_REAL& phi, 
											 const rg_REAL& psi, 
											 rg_INDEX& start, 
											 rg_INDEX& end);
	
	static rg_INT  getNumBBDepRotamersOfResidue(const ResidueCode& code, 
		                                        const rg_REAL& phi, 
												const rg_REAL& psi);

	static rg_INT  getNumBBDepRotamersOfResidue(const ResidueCode& code);

	static rg_INDEX getBBDepLibIndexOfRotamerWithHighestProbability(const ResidueCode& code, 
		                                                            const rg_REAL& phi, 
		                                                            const rg_REAL& psi);

	// functions for backbone-independent library
	static void  getRotamerSetsOfResidueSetFromBBIndepRotLib(Backbone             & backbone,
		                                                     Rotamer*               sortedAssignedRotamersOfProtein, 
		                                                     RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
		                                                     ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP   );

	static rg_BOOL getBBIndepLibIndexOfResidue(const ResidueCode& code, 
		                                       rg_INDEX& start, 
		                                       rg_INDEX& end);

	// common functions for backbone-dependent and -independent libraries

	static inline TypeOfRotamerLibrary getCurrentlyUsedRotamerLibType() { return m_currentlyUsedRotamerLibType; }
    static inline rg_BOOL              isBBINDLibLoaded()               { return m_bBBINDLibLoaded;  }
    static inline rg_BOOL              isBBDEPLibLoaded()               { return m_bBBDEPLibLoaded;  }

	static rg_INT  getNumDihedralAnglesOfResidue(const ResidueCode& code);

	static rg_INT  getDihedralAngles(const rg_INDEX& rotLidIndex, rg_REAL*& dihedralAngles);

	static void transformAtomCentersOfSidechainUsingDihedralAngle(Rotamer& rotamer, const rg_INT& rotamerIndex);
	static void transformAtomCentersOfSidechainUsingDihedralAngle(Rotamer& rotamer, const rg_REAL* dihedralAngles);
	static void transformAtomCentersOfSidechainUsingDihedralAngle(V::GeometryTier::Residue* targetResidue, const rg_INT& rotamerIndex);
	static void transformAtomCentersOfSidechainUsingDihedralAngle(V::GeometryTier::Residue* targetResidue, const rg_REAL* dihedralAngles);
	static void transformAtomCentersOfSidechainUsingDihedralAngle(V::GeometryTier::Residue* targetResidue, const rg_REAL* dihedralAngles, const rg_INT& numDihedralAngles);

		static void computeAtomCenterUsingDihedralAngle(V::GeometryTier::Atom* atom0, 
			                                            V::GeometryTier::Atom* atom1, 
														V::GeometryTier::Atom* atom2, 
														V::GeometryTier::Atom* atom3, 
														const rg_REAL& dihedralAngle);
		static void computeAtomCenterUsingDihedralAngle2(V::GeometryTier::Atom* atom0, 
			                                             V::GeometryTier::Atom* atom1, 
														 V::GeometryTier::Atom* atom2, 
														 V::GeometryTier::Atom* atom3, 
														 const rg_REAL& dihedralAngle);
        static void computeAtomCenterUsingDihedralAngle3(V::GeometryTier::Atom* atom0, 
                                                         V::GeometryTier::Atom* atom1, 
                                                         V::GeometryTier::Atom* atom2, 
                                                         V::GeometryTier::Atom* atom3, 
                                                         const rg_REAL& dihedralAngle);
		static void computePointUsingDihedralAngle(const rg_Point3D& point0, 
			                                       const rg_Point3D& point1,
			                                       const rg_Point3D& point2,
			                                       const rg_REAL   & bondLength23,
			                                       const rg_REAL   & bondAngle123,
			                                       const rg_REAL   & dihedralAngle,
			                                       rg_Point3D      & point3);

		static void computeSecondBranchAtomCenter(V::GeometryTier::Atom* atom0, 
			                                      V::GeometryTier::Atom* atom1, 
												  V::GeometryTier::Atom* atom2, 
												  V::GeometryTier::Atom* secondBranchAtom);
		static void computeSecondBranchAtomCenter2(V::GeometryTier::Atom* atom0, 
			                                       V::GeometryTier::Atom* atom1, 
												   V::GeometryTier::Atom* atom2, 
												   V::GeometryTier::Atom* atom3, 
												   V::GeometryTier::Atom* secondBranchAtom);
		static void computeTwoCandidatesForSecondBranchAtomCenter(V::GeometryTier::Atom* atom1, 
			                                                      V::GeometryTier::Atom* atom2, 
																  V::GeometryTier::Atom* atom3, 
																  V::GeometryTier::Atom* secondBranchAtom, 
																  rg_Point3D*& candidatesForAtomCenter);
	
		static void glueFixedStructureIntoResidue(V::GeometryTier::Residue* targetResidue);
        static void glueFixedStructureIntoResidue2(V::GeometryTier::Residue* targetResidue);
		static rg_BOOL hasFixedSubstructure(V::GeometryTier::Residue* targetResidue);
		static void getCoordinatesOfFixedStructure(const rg_INT& residueCode, rg_dList<rg_Point3D>& pointsOnFixedStructure );
			static bool getIndexOfFixedStructureInResidue(const rg_INT& residueCode, rg_INT& index);
        static void getCoordinatesOfFixedStructure2(const rg_INT& residueCode, rg_dList<rg_Point3D>& pointsOnFixedStructure );

	// Utility functions
	static rg_BOOL getStartNEndRemoteIndicatorsOfSidechainForResidue(const ResidueCode& code, RemoteIndicator& start, RemoteIndicator& end);
	static rg_BOOL getEndRemoteIndicatorNBranchDesignatorOfSidechainForResidue(const ResidueCode& code, RemoteIndicator& end, BranchDesignator& branchDesignator);
	static rg_BOOL getEndRemoteIndicatorNMaximumBranchDesignatorOfSidechainForResidue(const ResidueCode& code, RemoteIndicator& end, BranchDesignator& maxBranchDesignator);
	static rg_INT  getAtomsOfResidue(V::GeometryTier::Residue* residue, rg_dList<V::GeometryTier::Atom*>& atomsOfResidue);
	static rg_INT  getAtomsOfResidue(V::GeometryTier::Residue* residue, V::GeometryTier::Atom**& atomsOfResidue);
	static rg_INT  getAtomsOnSidechain(V::GeometryTier::Residue* residue, rg_dList<V::GeometryTier::Atom*>& atomsOnSidechain );
	static rg_INT  getAtomsOnSidechain(V::GeometryTier::Residue* residue, V::GeometryTier::Atom**& atomsOnSidechain );
	//static void  writeBBDepRotLibIntoBinaryFile(const string& fileNameWithPath);
	static void  writeBBIndepRotLibIntoBinaryFile(const string& fileNameWithPath);
	static void    loadBBDepRotLib(const string& fileNameWithPath);
	static void    loadBBIndepRotLib(const string& fileNameWithPath);
    static void    loadBBDepRotLib_Dunbrack10(const string& fileNameWithPath);
	//static void    loadBBDepRotLib(const string& fileNameWithPath, BBDepRotamer rotLib[SIZE_OF_BBDEP_ROTAMER_LIB]);

    //////////////////////////////////////////////////////////////////////////
	// For Dunbrack 2010 - by JHCha
	//////////////////////////////////////////////////////////////////////////

	static rg_BOOL getBBDepLibIndexOfResidue_2010(const ResidueCode& code,
		rg_INDEX& start,
		rg_INDEX& end);

	static rg_BOOL getBBDepLibIndexOfResidue_2010(const ResidueCode& code,
		const rg_REAL& phi,
		const rg_REAL& psi,
		rg_INDEX& start,
		rg_INDEX& end);

	static rg_INDEX getIndexOfResCode_2010(const ResidueCode& code);
	static rg_INT  getNumDihedralAnglesOfResidue_2010(const ResidueCode& code);
	static void  getRotamerSetsOfResidueSetFromBBDepRotLib_2010(Backbone             & backbone,
		Rotamer*               sortedAssignedRotamersOfProtein,
		RotamerSetOfResidue* & sortedRotamerSetsOfResidues,
		ManagerOfRotamerSetsForSCP    & managerOfRotamerSetsForSCP);

    static void    loadBBDepRotLib_Dunbrack2010_ASCII(const string& fileNameWithPath);
	static void    loadBBDepRotLib_Dunbrack2010_BIN(const string& fileNameWithPath);
	static void  writeBBDepRotLibIntoBinaryFile_2010(const string& fileNameWithPath);
	
	static ResidueCode translate_three_character_residue_code(char* threeCharResCode);
	static bool compare_BBDepRotLib_Dunbrack2010_BIN_to_ASCII(const string& fileNameWithPath);
	static bool compare_BBDepRotLib_Data(const BBDepRotamer& rotamer1, const BBDepRotamer& rotamer2);

	static bool test_BBDepLib_index_of_Residue_2010();

	//////////////////////////////////////////////////////////////////////////
	// For Dunbrack 2010 - by JHCha END
	//////////////////////////////////////////////////////////////////////////
};

#endif

