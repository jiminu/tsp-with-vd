#ifndef _MOLECULE_H
#define _MOLECULE_H

// #pragma warning (once: 1179)
// #pragma warning (once: 1266)
// #pragma warning (once: 2138)
#pragma warning (once: 4081)
#pragma warning (disable: 4786)

#include "ConstForMolecule.h"
#include "rg_Atom.h"
#include "ChemicalBond.h"
#include "Residue.h"
#include "Chain.h"
#include "rg_Point3D.h"
#include "Sphere.h"
#include "EnclosingSphereOfSpheres.h"
#include "PharmaFeature.h"
#include "ConstForPharmacophore.h"
#include "Backbone.h"
#include "ManagerOfRotamers_Assigned_At_Residues.h"
#include "ChainIDWithSequenceNumbersOfMissingResidues.h"
#include "ProteinConformation.h"

#include <iostream>
#include <string>
#include <bitset>
#include <map>
using namespace std;



namespace V {
namespace GeometryTier {



typedef map<int, Atom*>          AtomMap;
typedef map<int, Residue*>       ResidueMap;
typedef map<int, Chain*>         ChainMap;
typedef map<int, ChemicalBond*>  ChemBondMap;
typedef map<int, PharmaFeature*> PharmaFeatureMap;

typedef map<string, int>      AtomSymbolMap;

class Molecule
{
private:

    string                   m_moleculeFileName;
    rg_dList<string>         m_headerRecords;
    AtomRadiusType           m_atomRadiusType;

    rg_dList<Atom>           m_atoms;
    rg_dList<ChemicalBond>   m_chemicalBonds;
    rg_dList<Residue>        m_residues;
    rg_dList<Chain>          m_chains;

    list<PharmaFeature>      m_pharmaFeatures;

    rg_Point3D               m_centerOfMass;
    Sphere                   m_minEnclosingSphere;

    rg_INT                   m_modelSerialFromInputFile;


    // JKKIM ADDED //////////
    string                   m_moleculeName;

	bool					 m_isCryst;
	float					 m_cryst[6];


    // by Y.Cho at 2011-10-12
    string                   m_timeStamp;
    rg_INT                   m_fileSize;    // byte.

	string					 m_method;
	float					 m_resolution;

private:
    void duplicate(const Molecule& aMolecule);
        void    makeMapsOfEntitiesInMolecule(
                    const Molecule& aMolecule,
                    map<Atom*, Atom*>&           mapAtom2Atom,
                    map<Residue*, Residue*>&     mapResidue2Residue,
                    map<Chain*, Chain*>&         mapChain2Chain,
                    map<ChemicalBond*, ChemicalBond*>&   mapBond2Bond );

        void    setEntitiesAndRelationshipAmongEntities(
                    const map<Atom*, Atom*>&           mapAtom2Atom,
                    const map<Residue*, Residue*>&     mapResidue2Residue,
                    const map<Chain*, Chain*>&         mapChain2Chain,
                    const map<ChemicalBond*, ChemicalBond*>&   mapBond2Bond);

public:
    //  CONSTRUCTOR & DECONSTRUCTORs
    Molecule();
//    Molecule( const rg_INT& ID );
    Molecule( const Molecule& aMolecule );
    ~Molecule();

    //  GET FUNCTION
    string                   getMoleculeFileName() const;
    const rg_dList<string>&  getHeaderRecords() const { return m_headerRecords; }

    AtomRadiusType           getAtomRadiusType() const;
    string                   getDescriptionOfAtomRadiusType() const;

    inline string            getMoleculeName() const { return m_moleculeName; }

    const rg_dList<Atom>&    getAtomsInMolecule() const { return m_atoms; }
    rg_dList<Atom>*          getAtoms();
    void                     getPtrAtoms(rg_dList<Atom*>& ptrAtoms) const;
    rg_dList<ChemicalBond>*  getChemicalBonds();
    rg_dList<Residue>*       getResidues();
    rg_INT                   getNonAminoResidues(Residue**& residueArr);

	rg_INT                   getAminoResidues(Residue**& residueArr);
    rg_INT                   getNumNonAminoResidues();
    rg_INT                   getNumAminoResidues();
    rg_dList<Chain>*         getChains();
    Chain*                   getChain( const rg_INT& chainIDInDecimal );



	inline	float	getResolution()	{return m_resolution;}
	inline	string	getMethod()		{return m_method;}

	inline	void	setResolution(const float& resolution)	{m_resolution = resolution;}
	inline	void	setMethod(const string& method)			{m_method = method;}


	inline	bool	isCryst()					{return m_isCryst;}
	inline	float	getCrystInfo(const int& index)	{return m_cryst[index];}
	
	inline	void	setCryst(const bool& isCryst)						    {m_isCryst = isCryst;}
	inline	void	setCrystInfo(const int& index, const float& crystInfo)	{m_cryst[index] = crystInfo;}

    
//  void                     getResiduesSortedBySequenceNumber( rg_dList<Residue*>* targetListOfResidues );         // REPLACED TO: getResiduesSortedBySequenceNumber( rg_dList<Residue*>& targetListOfResidues )
//  void                     getAtomsOnBackboneOrderedBySequenceNumber( rg_dList<Atom*>* listOfAtomsOnBackbone );   // REPLACED TO: getAtomsOnBackboneSorteddBySequenceNumber( rg_dList<Atom*>& listOfAtomsOnBackbone )
    void                     getResiduesSortedBySequenceNumber( rg_dList<Residue*>& targetListOfResidues ); 
	void                     getResiduesSortedBySequenceNumberForEachChain(rg_dList<Residue*>& targetListOfResidues);
	rg_INT                   getResiduesSortedBySequenceNumberForEachChain(Residue**& sortedResidues);
	// Dec. 30, 2017 by Joonghyun
	rg_INT                   getAllResiduesOrderedByChainIDSequenceNumber(Residue**& sortedResidues);
	rg_INT                   getAminoResiduesOrderedByChainIDSequenceNumber( Residue**& sortedResidues );
    rg_INT                   getNonAminoResiduesOrderedByChainIDSequenceNumber(Residue **& sortedResidues);
	rg_INT                   getAminoResiduesOrderedByChainIDSequenceNumber( rg_dList<Residue*>& residuesOrderedByChainIDSeqNum );
    void                     getAtomsOnBackboneSortedBySequenceNumber( rg_dList<Atom*>& listOfAtomsOnBackbone );
	void                     getAtomsOnBackboneSortedByChainIDSequenceNumber( rg_dList<Atom*>& listOfAtomsOnBackbone );
	rg_INT                   getAtomsOnBackboneSortedByChainIDSequenceNumber( Atom**& backboneAtoms );

	rg_INT                   getSumAtomNumsOfAminoResidues() const;

	rg_BOOL                  haveMissingResidues() const;
	rg_INT                   findMissingResidues(rg_dList<ChainIDWithSequenceNumbersOfMissingResidues>& chainIDsWithSeqNumbersOfMissingResidues) const;

    void                     getListOfHAtomsFromHDonor( rg_dList<Atom*>* targetListOfHydrogenAtoms );
    void                     getListOfHAcceptorAtoms( rg_dList<Atom*>* targetListOfHAcceptorAtoms );

    rg_Point3D               getCenterOfMass();
    Sphere                   getMinEnclosingSphere();
    rg_REAL                  getMeanRadiusOfAtoms();
    rg_REAL                  getMaxRadiusOfAtoms();

    rg_REAL                  getMolecularWeight();

    //  by Youngsong Cho, 2012.03.14.
    rg_INT                   getNumStandardResidues() const;

	rg_INT                   getNumProteinChains() const;
    rg_INT                   getNumberOfChainsWithNonStandardResidues();
    rg_INT                   getNumberOfChainsWithNonStandardResiduesAndHOH();

    rg_INT                   getModelSerialFromInputFile() const;

    //  by Youngsong Cho, 2018.10.30.
    rg_INT                   getNumberOfAtoms() const { return m_atoms.getSize(); }
    rg_INT                   getNumberOfAtomsExcludingHOH() const;
    rg_INT                   getNumberOfAtomsInStandardResidues() const;
    rg_INT                   getNumberOfAtomsInNonStandardResidues() const;
    rg_INT                   getNumberOfResidues() const { return m_residues.getSize(); }
    rg_INT                   getNumberOfStandardResidues() const;
    rg_INT                   getNumberOfNonStandardResidues() const;
    rg_INT                   getNumberOfChains() const { return m_chains.getSize(); }

    rg_INT                   getAtomsInStandardResidues(list<Atom*>& atomsInStandardResidues) const;
    rg_INT                   getAtomsInNonStandardResidues(list<Atom*>& atomsInNonStandardResidues) const;


    //  by Youngsong Cho, 2019.01.15.
    rg_INT                   evaluateAtomFrequency(map<AtomCode, int>& atomFrequency) const;
    rg_INT                   evaluateAtomFrequency(map<string, int>& atomFrequency) const;

    // by Y.Cho at 2011-10-12
    string                   getTimeStamp() const;
    rg_INT                   getFileSize() const;

    void                     getConformation(ProteinConformation& proteinConformation);
	void                     getBackboneSortedByChainIDSeqNumber(Backbone& backbone);
	void                     getRotamersCorrToSidechainsOfAllResidueInstancesOrderedByChainIDSeqNumber(ManagerOfRotamers_Assigned_At_Residues& managerOfRotamers_Assigned_At_Residues);
	// Dec. 30, 2017 renamed by Joonghyun
    rg_INT                   getVirtualRotamersCorrToSidechainsOfAllNonAminoResidueInstancesOrderedByChainIDSeqNumber(Rotamer*& rotamersOderedByChainIDSeqNumber);
	rg_INT                   getRotamersCorrToSidechainsOfAllAminoResidueInstancesOrderedByChainIDSeqNumber(Rotamer*& sortedRotamerSet);
	// Dec. 30, 2017 by Joonghyun
	rg_INT                   getRotamersCorrToSidechainsOfAllResidueInstancesOrderedByChainIDSeqNumber(Rotamer*& sortedRotamerSet);

    //  SET FUNCTION
    void                     setMoleculeFileName( const string& moleculeFileName );

    inline void              setMoleculeName( const string& moleculeName ) {m_moleculeName = moleculeName;}

    void                     setAtomRadiusType( const AtomRadiusType& atomRadiusType );
    string*                  addHeaderRecords( const string& recLine );
    Atom*                    addAtom( const Atom& atom );
    ChemicalBond*            addChemicalBond( const ChemicalBond& chemicalBond );
    Residue*                 addResidue( const Residue& residue );
    Chain*                   addChain( const Chain& aChain);

    void                     setCenterOfMass( const rg_Point3D& centerOfMass );
    void                     setMinEnclosingSphere( const Sphere& minEnclosingSphere );

    void                     setModelSerialFromInputFile( const rg_INT& modelSerialFromInputFile );
    
    // by Y.Cho at 2011-10-12
    void                     setTimeStamp( const string& timeStamp );
    void                     setFileSize(  const rg_INT& filesize );

    // FUNCTIONS FOR MODIFICATION
    void                     moveOxygenInCarboxylGroupToBackbone(); // Changed from: addOxygenInCarboxylGroupToBackbone()

    void                     deleteAtomsInSideChains( rg_BOOL replaceResidueCodeToUnknown, const RemoteIndicator& startRemoteIndicator, const ResidueType& typeOfResidue );
    void                     deleteAtomsInSideChain( rg_BOOL replaceResidueCodeToUnknown, const RemoteIndicator& startRemoteIndicator, Residue* targetResidue );
    void                     deleteAtom( rg_BOOL replaceResidueCodeToUnknown, Atom* atomToDelete );

    void                     recoverMissingAtomsInResiduesWithoutCoordinates( const ResidueType& typeOfResidue );
    void                     recoverMissingAtomsInResidueWithoutCoordinates( Residue* targetResidue );
        rg_BOOL                 isHydrogenAtomName( const string& atomName );
        rg_INT                  getBiggestAtomSerialFromInputFile();

	void                     recoverMissingAtomsInResidueWithGivenCoordinates( Residue* targetResidue, 
		                                                                       const rg_Point3D& defaultAtomCenter, 
																			   AtomSymbolMap& mapOfAtomSymbols );
	
	void                     recoverMissingAtomsOnBackboneOfResidueWithGivenCoordinates( Residue* targetResidue, 
		                                                                                 const rg_Point3D& defaultAtomCenter, 
		                                                                                 AtomSymbolMap& mapOfAtomSymbols );

	void                     recoverMissingAtomsInSidechainOfResidueWithGivenCoordinates( Residue* targetResidue, 
		                                                                                  const rg_Point3D& defaultAtomCenter, 
		                                                                                  AtomSymbolMap& mapOfAtomSymbols );

	void                     recoverMissingAtomsInResidueWithGivenCoordinates( Residue* targetResidue, 
		                                                                       const rg_Point3D& defaultAtomCenter,
																			   const rg_INT& startIndexForPairOfBondedAtom,
																			   const rg_INT& endIndexForPairOfBondedAtom,
		                                                                       AtomSymbolMap& mapOfAtomSymbols );

	void                     resetIDs();
		void resetSerialsFromInputFileOfAtoms();

	// Modeling and designing
	void                     addOxygenOfCarboxylGroupWithCenter(Atom* carbonInCarboxylGroup, 
		                                                        const rg_Point3D& oxygenCenter, 
																Residue* targetResidue,
																AtomSymbolMap& atomSymbolMap);
	Residue*                 addNewRotamerWithFreezedBetaCarbon(Residue* prevResidue,
									                            Residue* nextResidue,
										                        const ResidueCode& residueCode,
										                        Chain*   chain,
										                        Residue* refResidue);

	Residue*                 addNewRotamerWithFreezedBetaCarbon(Residue* prevResidue,
									                            Residue* nextResidue,
										                        const ResidueCode& residueCode,
										                        Chain*   chain,
										                        rg_dList<Atom*>& freezedAtoms);

	Residue*                 addNewResidueWithoutSpecificationOfAtomCoordinates(Residue* prevResidue,
									                                            Residue* nextResidue,
										                                        const ResidueCode& residueCode,
										                                        Chain*   chain);

	Residue*                 addNewResidue(Residue* prevResidue,
									       Residue* nextResidue,
										   Residue* newResidue);

	Residue* makeAminoResidueNConstituentAtoms(const ResidueCode& residueCode,
		                                       Chain* chain);
		void makeAndAddAtomsOfResidue(Residue* targetResidue, 
			                          const rg_INT& startAtomID, 
								      const rg_INT& startAtomSerial);
		void makeChemicalBondsOfResidue(Residue* targetResidue, 
			                            rg_dList<Atom>& atomsOfResidue,
										const rg_INT& startChemBondID,
										rg_dList<ChemicalBond>& chemicalBondsOfResidue);

		void disconnectBondBetween(Residue* firstResidue,
		                           Residue* secondResidue);


	void removeAminoResidueNConstituentAtoms(Residue* prevResidue,
		                                     Residue* nextResidue,
								             Residue* targetResidue);

	void removeAminoResidueNConstituentAtoms(Residue* targetResidue);
		void removeAtomsOfResidue(Residue* targetResidue);
		void disconnectAndRemoveBondsOfResidue(Residue* targetResidue);
		void disconnectAndRemoveBondsOfAtom(Atom* targetAtom);
		void disconnectAndRemoveBondBetween(Atom* firstAtom, Atom* secondAtom);

	Residue* replaceAminoResidue(Residue* targetResidue, 
							     Residue* substitue);

	Residue* replaceAminoResidue(Residue* prevResidue,
		                         Residue* nextResidue,
		                         Residue* targetResidue, 
							     Residue* substitue);

	//void replaceAminoResidue(Residue* targetResidue, const ResidueCode& newResidueCode, const rg_BOOL& bBetaCarbonComputed = rg_FALSE);
    void replaceAminoResidue(Residue* targetResidue, const ResidueCode& newResidueCode);
		void removeAtomsNTheirChemicalBondsOfSidechain(Residue* targetResidue);
		void removeHydrongenInSidechain(Residue* targetResidue);

		void removeChemicalBondsOnSidechaOfResidue(Residue* targetResidue);
		void removeAtomsOnSidechainOfResidue(Residue* targetResidue);

    // TO BE MOVED TO NEW CLASS
    list<PharmaFeature>*     getListPharmaFeatures();

   void                     defineAndSetPharmaFeatures();
       void                     definePIPharmaFeatures();            
           rg_FLAG                  isPiPositiveCharge( Atom* currAtom, PharmaFeature& piPositiveChargeFeature );
           rg_FLAG                  isPiAmine( Atom* currAtom, PharmaFeature& piAmineFeature  );
           rg_FLAG                  isPiAmidine( Atom* currAtom, PharmaFeature& piAmidineFeature );
			rg_FLAG                  isPiGuanidine( Atom* currAtom, PharmaFeature& piGuanidineFeature );

       void                     defineNIPharmaFeatures();
           rg_FLAG                 isNiNegativeCharge( Atom* currAtom, PharmaFeature& niNegativeChargeFeature );
           rg_FLAG                 isNiCarboxyl( Atom* currAtom, PharmaFeature& niCarboxylFeature );
           rg_FLAG                 isNiTriflu( Atom* currAtom, PharmaFeature& niTrifluFeature );
           rg_FLAG                 isNiSulfinic( Atom* currAtom, PharmaFeature& niSulfinicFeature );
           rg_FLAG                 isNiSulfonic( Atom* currAtom, PharmaFeature& niSulfonicFeature );
           rg_FLAG                 isNiSulfuric( Atom* currAtom, PharmaFeature& niSulfuricFeature );
           rg_FLAG                 isNiPhosphinic( Atom* currAtom, PharmaFeature& niPhosphinicFeature );
           rg_FLAG                 isNiPhosphonic( Atom* currAtom, PharmaFeature& niPhosphonicFeature );
           rg_FLAG                 isNiPhosphoric( Atom* currAtom, PharmaFeature& niPhosphoricFeature );

       void                     defineHBAPharmaFeatures();
           rg_FLAG                 isHbaNlonpair( Atom* currAtom, PharmaFeature& hbaNlonpairFeature );
           rg_FLAG                 isHbaOlonpair( Atom* currAtom, PharmaFeature& hbaOlonpairFeature );
           rg_FLAG                 isHbaSlonpair( Atom* currAtom, PharmaFeature& hbaSlonpairFeature );

       void                     defineHBDPharmaFeatures();
           rg_FLAG                 isHbdHydroxyl( Atom* currAtom, PharmaFeature& hbdHydroxylFeature );
           rg_FLAG                 isHbdThiol( Atom* currAtom, PharmaFeature& hbdThiolFeature );
           rg_FLAG                 isHbdAcetylene( Atom* currAtom, PharmaFeature& hbdAcetyleneFeature );
           rg_FLAG                 isHbdNh( Atom* currAtom, PharmaFeature& hbdNhFeature );

       void                     defineHYPharmaFeatures();
           rg_FLAG                 isHyCRing( Atom* currAtom, list<PharmaFeature>& listHyCRingFeature );
           rg_FLAG                 isHyHalogen( Atom* currAtom, PharmaFeature& hyHalogenFeature );
           rg_FLAG                 isHyCGroup( Atom* currAtom, PharmaFeature& hyCGroupFeature );
               rg_FLAG                 isHyCRingChemType( Atom* currAtom );
               rg_FLAG                 isBondedWithCarbonAndHydrogen( Atom* currAtom );
               rg_FLAG                 isTargetAtomExistInAtomList( Atom* targetAtom, list<Atom*>* atomList );
               
       //void                    refineListOfPharmaFeatures();

               rg_FLAG                  isAllBondedAtomsNonNegative( Atom* currAtom );
               rg_FLAG                  isAllBondedAtomsNonPositive( Atom* currAtom );
               rg_FLAG                  isAllBondsSingle( Atom* currAtom );
               void                     countEachBondTypeFromBondListInAtom( Atom* currAtom, rg_INT& numOfSingleBond, rg_INT& numOfDoubleBond, rg_INT& numOfTripleBond, rg_INT& numOfAromaticBond );
               void                     countEachAtomTypeFromBondListInAtom( Atom* currAtom, rg_INT& numOfCarbon, rg_INT& numOfHydrogen, rg_INT& numOfNitrogen, rg_INT& numOfOxygen, rg_INT& numOfPhosphorus, rg_INT& numOfSulfur);
               rg_INT                   getNumOfBondedAtomTypeFromCurrAtom( Atom* currAtom, const AtomCode& targetAtomCode );
               void                     getEachpAtomFromBondListInAtom( Atom* currAtom, Atom** arrCarbon, Atom** arrHydrogen, Atom** arrNitrogen, Atom** arrOxygen, Atom** arrPhosphorus, Atom** arrSulfur );
               void                     getBondedAtomFromCurrAtom( Atom* currAtom, const AtomCode& targetAtomCode, const BondType& targetBondType, Atom** arrBondedAtom );
               void                     getBondedAtomFromCurrAtom( Atom* currAtom, const AtomCode& targetAtomCode, Atom** arrTargetAtom );
               void                     getBondedAtomFromCurrAtomExceptPrevAtom( Atom* currAtom, Atom* prevAtom, Atom** arrBondedAtom );
               void                     getBondedAtomFromCurrAtomExceptPrevAtom( Atom* currAtom, Atom* prevAtom, const AtomCode& targetAtomCode, Atom** arrTargetAtom );
               rg_FLAG                  isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( Atom* currAtom, Atom* prevAtom );
               void                     addPharmaFeature( rg_INT& featureID, PharmaFeature* aPharmaFeature );
               void                     addPharmaFeature( rg_INT& featureID, list<PharmaFeature>* aListPharmaFeature );
				rg_FLAG					 isOnRing( Atom* targetAtom, const rg_INT& numOfAtomsInRing, list<Atom**>& listOfAtomArrayInRing );
					void					getpBondsFromBondListInAtomExceptPrevBond( Atom* currAtom, ChemicalBond* prevBond, ChemicalBond** arrBonds );
                   rg_FLAG                 isValidBondForRing( const rg_INT& atomIDOfRing, ChemicalBond** arrBondsOfRingCandidate, ChemicalBond* targetBond );
                   rg_FLAG                 isIdenticalRing( const rg_INT& numOfAtomsInRing, ChemicalBond** arrBondsOfRingA, ChemicalBond** arrBondsOfRingB );

               bitset<MAX_CHEM_ATOM_NUM>   getBitSetKeyOfPharmaFeatureForMap( PharmaFeature* aPharmaFeature );



    ///// FOLLOWING 5 FUNCTIONS ARE TO BE MOVED TO PDB READER.
    void                     computeAndSetCenterOfMass();
    void                     computeAndSetMinEnclosingSphere();
    
    void                     evaluateHydrogenDonorAndAcceptorAtoms();
    void                     evaluateHydrogenDonorAtoms();
    void                     evaluateHydrogenAcceptorAtoms();
    


	//JKKIM ADDED
	void makePeoridicStructure();



    //  OPERATOR OVERLOADING
    Molecule& operator =( const Molecule& aMolecule );

    void syncPtrsWithSourceMolecule( Molecule* sourceMolecule );
        void initializeMapsForSync( AtomMap& mapOfAtom, ResidueMap& mapOfResidue, ChainMap& mapOfChain, ChemBondMap& mapOfChemBond, PharmaFeatureMap& mapOfPharmaFeature);
        void syncPtrsInAtom( rg_dList<Atom>* sourceAtoms, ResidueMap* mapOfResidue, ChemBondMap* mapOfChemBond, PharmaFeatureMap* mapOfPharmaFeature );
        void syncPtrsInResidue( rg_dList<Residue>* sourceResidues, AtomMap* mapOfAtom, ChainMap* mapOfChain );
        void syncPtrsInChain( rg_dList<Chain>* sourceChains, ResidueMap* mapOfResidue );
            void syncPtrsInSecondaryStructure( Chain* sourceChain, ResidueMap* mapOfResidue, Chain* targetChain );
        void syncPtrsInChemBond( rg_dList<ChemicalBond>* sourceBonds, AtomMap* mapOfAtom );
        void syncPtrsInPharmaFeatures( list<PharmaFeature>* sourcePharmaFeatures, AtomMap* mapOfAtom );
    
};



} // namespace GeometryTier
} // namespace V


#endif

