#ifndef _MOLECULE_IO_FUNCTIONS_H
#define _MOLECULE_IO_FUNCTIONS_H

#include "rg_Const.h"
#include "rg_Molecule.h"
#include "ConstForMoleculeIOFunctions.h"


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <list>
using namespace std;



namespace V {
namespace GeometryTier {



typedef pair<int, ResidueMap> ResidueMapWithChainID;    // first: chainID, second: ResidueMap
typedef map<string, int>      AtomSymbolMap;
typedef map<string, int>      ResidueSymbolMap;


namespace MoleculeIOFunctions
{
    
    // PDB IO Functions
    rg_FLAG readPDBFile( const char* pathFile, rg_dList<Molecule>& molecules );
    rg_FLAG readPDBFile( const rg_INT& targetModelID, const char* pathFile, Molecule& aMolecule );
        string  parseDepositionDate(const list<string>& recordsOfPDBFile);

        void    addRecordsOfPDBFIleToList( ifstream* fin, list<string>& recordsOfPDBFile );
        rg_FLAG setAtomRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidue, ResidueMap& mapOfNonStdResidue, ChainMap& mapOfChain, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule );
        
        //// FOLLOWING setAtomRecordsToMolecule FUNCTION IS REPLACED WITH setHeteroAtomRecordsToMolecule...
        //
        rg_FLAG setAtomRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidue, rg_dList<ResidueMapWithChainID>& mapsOfNonStdResidue, ChainMap& mapOfChain, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule );
            void    setResidueCodeToTargetResidueExceptNeucleicAcid( ResidueSymbolMap& mapOfResidueSymbol, const string& threeCodeNameOfRes, Residue* targetResidue );
            rg_FLAG setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( AtomSymbolMap& mapOfAtomSymbol, const string& atomName, Atom* targetAtom );
            RemoteIndicator  convertStringIntoRemoteIndicator( const string& strRemoteIndicator );
            BranchDesignator convertStringIntoBranchDesignator( const string& strBranchDesignator );

        rg_FLAG setHeteroAtomRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidueForHeteroRec, ChainMap& mapOfChainForHeteroRec, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule );

        rg_FLAG setConectRecordsToMolecule( list<string*>* connectRecLinesOfPDBFile, AtomMap& mapOfAtom, Molecule* aMolecule );
            rg_FLAG setConectRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, Molecule& aMolecule );
            rg_FLAG setConectRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, Molecule* aMolecule );
                void    setChemicalBondToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, Molecule* aMolecule );
                //JKKIM modified return void -> ChemicalBond*
                ChemicalBond* setChemicalBondToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, Molecule& aMolecule );
                //void    setChemicalBondToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, Molecule& aMolecule );
                rg_FLAG isChemicalBondExistInChemicalBondListInAtom( ChemicalBond* aChemicalBond, Atom* atomForBond, ChemicalBond* existingChemicalBond );
        void    setChemicalBondsToMoleculeForStandardResidues( Molecule* aMolecule );
        void    setChemicalBondsToMoleculeForStandardResidues( Molecule& aMolecule );
            void    setChemicalBondsBetweenStdResiduesInChain( Chain* aChain, Molecule* aMolecule );
            void    setChemicalBondsBetweenStdResiduesInChain( Chain* aChain, Molecule& aMolecule );
                void    setChemicalBondBetweenAminoAcids( Residue* firstAminoAcid, Residue* secondAminoAcid, Molecule* aMolecule );
                void    setChemicalBondBetweenAminoAcids( Residue* firstAminoAcid, Residue* secondAminoAcid, Molecule& aMolecule );
                void    setChemicalBondBetweenNucleicAcids( Residue* firstNucleicAcid, Residue* secondNucleicAcid, Molecule* aMolecule );
                void    setChemicalBondBetweenNucleicAcids( Residue* firstNucleicAcid, Residue* secondNucleicAcid, Molecule& aMolecule );
            void    setChemicalBondsBetweenAtomsInStdResidue( Residue* aResidue, Molecule* aMolecule );
            void    setChemicalBondsBetweenAtomsInStdResidue( Residue* aResidue, Molecule& aMolecule );
        void setChainsForNonStdResiduesForPDBFile( ResidueMap* mapOfNonStdResidue, Molecule& aMolecule );
        void setChainsForNonStdResiduesForPDBFile( rg_dList<ResidueMapWithChainID>* mapsOfNonStdResidue, Molecule& aMolecule );

        void setHelixRecordsToMolecule( list<string*>* helixRecLinesOfPDBFile, rg_dList<ResidueMapWithChainID>* mapsOfResidue, Molecule* aMolecule );
            void setHelixRecordsToMolecule( const string& strRecLine, rg_dList<ResidueMapWithChainID>* mapsOfResidue, rg_dList<Residue*>* residuesSortedBySequence );
        void setSheetRecordsToMolecule( list<string*>* sheetRecLinesOfPDBFile, rg_dList<ResidueMapWithChainID>* mapsOfResidue, Molecule* aMolecule );
            void setSheetRecordsToMolecule( const string& strRecLine, rg_dList<ResidueMapWithChainID>* mapsOfResidue, rg_dList<Residue*>* residuesSortedBySequence );
        void setTurnRecordsToMolecule( list<string*>* TurnRecLinesOfPDBFile, rg_dList<ResidueMapWithChainID>* mapsOfResidue, Molecule* aMolecule );
            void setTurnRecordsToMolecule( const string& strRecLine, rg_dList<ResidueMapWithChainID>* mapsOfResidue, rg_dList<Residue*>* residuesSortedBySequence );


        Residue* findResidueFromMapsOfResidue( const rg_INT& residueSequence, const rg_INT& chainIDinDecimal, rg_dList<ResidueMapWithChainID>* mapsOfResidue );
        void     findResiduesFromListSortedBySequence( rg_dList<Residue*>* residuesSortedBySequence, Residue* startResidue, Residue* endResidue, rg_dList<Residue*>& targetResidues );
        Atom*    findAtomFromMapsOfResidue( const string& atomName, const string& chainID, const rg_INT& residueSeqence, rg_dList<ResidueMapWithChainID>* mapsOfResidue );

    rg_FLAG writeRotatedPDBFileWithGivenAngle(const char* pathFile, Molecule* aMolecule, const rg_REAL& angleToRotate, const rg_INT& axis);
    rg_FLAG writePDBFile( const char* pathFile, Molecule* aMolecule );
        void    writeHeader( Molecule* aMolecule, ofstream& fout );
        void    writeAtomField( Molecule* aMolecule, ofstream& fout );
        void    writeAtomFieldInPDBFormat(const list<Atom*>& atoms, ofstream& fout);
        void    writeConnectField( Molecule* aMolecule, ofstream& fout );
            rg_FLAG isGeneralAtomInAminoOrDnaOrRnaResidue( Atom* anAtom );
        void    writeHeteroAtomField( Molecule* aMolecule, ofstream& fout );

        void    writeAtomField( Molecule* aMolecule, ofstream& fout, const rg_REAL& angleToRotate, const rg_INT& axis );
        void    writeHeteroAtomField( Molecule* aMolecule, ofstream& fout, const rg_REAL& angleToRotate, const rg_INT& axis );

	rg_FLAG writePDBFileForSCP( const char* pathFile, Molecule* aMolecule );
	rg_FLAG writePDBFileForSCP(const char* pathFile, Molecule* aMolecule, Molecule* refereneceMolecule);
	void reset_atom_serial_numbers(Molecule* aMolecule, const rg_INT& startNumber);
		void writeAtomFieldInOrderOfSequence( Molecule* aMolecule, ofstream& fout );
		void writeAtomFieldInOrderOfSequence2(Molecule* aMolecule, ofstream& fout);

    // JKKIM ADD 2013-04-15 ///////////////////////////////
	rg_FLAG readNextMoleculeFromGivenInputStream( ifstream& fin, Molecule& aMolecule, rg_FLAG& isNextExist );
		rg_FLAG setMoleculeRTIToMolecule( ifstream& fin, rg_INT* arrNumOfRecords, Molecule& aMolecule, string& strRecLine );
		rg_FLAG setAtomRTIToMolecule( ifstream& fin, const rg_INT& numOfAtomsGiven, AtomMap& mapOfAtom, ResidueMap& mapOfResidue, ResidueMap& mapOfNonStdResidue, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule& aMolecule, string& strRecLine );
		rg_FLAG setBondRTIToMolecule( ifstream& fin, const rg_INT& numOfBondsGiven, AtomMap& mapOfAtom, Molecule& aMolecule, string& strRecLine );
	///////////////////////////////////////////////////////

    // Tripos-mol2 IO Functions
    rg_FLAG readTriposMol2File( const char* pathFile, rg_dList<Molecule>& aMolecule );
    rg_FLAG readTriposMol2File( const rg_INT& targetModelID, const char* pathFile, Molecule& aMolecule );
        rg_BOOL isNewTriposRTILineStart( const string& strRecLine );
        void addRecordsOfMol2FIleToList( ifstream* fin, list<string>& recordsOfMol2File );
        rg_FLAG setMoleculeRTIToMolecule( list<string>::iterator& i_recLines, rg_INT* arrNumOfRecords, Molecule& aMolecule );
//        CAN'T READ MOL2 FILE... NEED TO DEBUG....
        rg_FLAG setAtomRTIToMolecule( list<string>::iterator& i_recLines, list<string>* recordsOfMol2File, const rg_INT& numOfAtomsGiven, AtomMap& mapOfAtom, ResidueMap& mapOfResidue, ResidueMap& mapOfNonStdResidue, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule& aMolecule );
//      rg_FLAG setAtomRTIToMolecule( list<string>::iterator& i_recLines, list<string>* recordsOfMol2File, const rg_INT& numOfAtomsGiven, AtomMap& mapOfAtom, ResidueMap& mapOfResidue, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule& aMolecule );
            void     setResidueCodeToTargetResidue( ResidueSymbolMap& mapOfResidueSymbol, const string& threeCodeNameOfRes, Residue* targetResidue );
            void     extractAtomRecordsFromAtomRTILine( const string&  strRecLine, rg_INT* arrIntRecs, rg_REAL* arrRealRecs, string* arrStrRecs );
            Residue* estimateResidueRecords( const rg_INT& substID, const string& substName, ResidueMap& mapOfResidue, ResidueSymbolMap& mapOfResidueSymbol, Molecule& aMolecule, rg_BOOL& isNewResidueCreated );
            Residue* getResiduePtrFromMap( ResidueMap* mapOfResidue,  const rg_INT& residueID );
            string   getResidueNameFromTriposSubstName( const string& substName );
            rg_FLAG  setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForMol2File( AtomSymbolMap& mapOfAtomSymbol, const string& atomName, const string& atomType, Atom* targetAtom );
            rg_FLAG  setChemicalPropertiesToNewAtomFromTriposAtomName( const string& atomName, const string& atomType, Atom* targetAtom );
            AmberAtomTypes getAmberAtomTypeFromSYBYL( const string& atomType );
        rg_FLAG setBondRTIToMolecule( list<string>::iterator& i_recLines, list<string>* recordsOfMol2File, const rg_INT& numOfBondsGiven, AtomMap& mapOfAtom, Molecule& aMolecule );
            void     extractBondRecordsFromBondRTILine( const string&  strRecLine, rg_INT& bondID, rg_INT& originAtomID, rg_INT& targetAtomID, string& bondType, string& statusBit );
            void setChemicalBondTypeForRTIBond(ChemicalBond* aChemicalBond, const string& bondType);
        void setAmberAtomTypeForHAtom(Molecule& aMolecule);
        void setChainsForNonStdResiduesForMol2File( ResidueMap* mapOfNonStdResidue, Molecule& aMolecule );
            
            /*AmberAtomTypes setAtomsTypeInAmber()( );*/

    void    filterEmptyChain( Molecule& aMolecule );
    void    evaluateAndSetChainCode( Molecule& aMolecule );
    void    evaluateAndSetResidueCodesForNeucleicAcid(  ResidueSymbolMap* mapOfResidueSymbolForDNA, ResidueSymbolMap* mapOfResidueSymbolForRNA, Molecule& aMolecule );
    // void    checkChainsInMoleculeForStandardResidue( Molecule& aMolecule );

    rg_FLAG isStringValidForAtomType( const string& aString );
    rg_FLAG isStringValidForRemoteIndicator( const string& aString );
    rg_FLAG isStringValidForBranchDesignator( const string& aString );
    rg_FLAG isStringValidForExtraBranchDesignator( const string& aString );


    rg_FLAG writerTriposMol2File( const char* pathFile, Molecule* aMolecule );
        void    writeTriposMoleculeFile( Molecule* aMolecule, ofstream& fout );
        void    writeTriposAtomField( Molecule* aMolecule, ofstream& fout );
        void    writeTriposBondField( Molecule* aMolecule, ofstream& fout );
        

    // CT-mol IO Functions
    rg_FLAG readCTMolFile( const rg_INT& targetMoleculeID, const char* pathFile, Molecule& aMolecule );
        void     readRecordsInHeaderBlock( list<string>::iterator& i_recLines, Molecule& aMolecule );
        rg_FLAG  readRecordsInCTAB( list<string>::iterator& i_recLines, AtomSymbolMap& mapOfAtomSymbol, AtomMap& mapOfAtom, Molecule& aMolecule );
        //AtomCode convertAtomSymbolToAtomCode( const string& atomSymbol );
        BondType convertIntegerToBondType( const rg_INT& intBondType );
        //Atom*    getpAtomFromMap( const rg_INT& atomID );
        void     setChemicalBondForMolFileToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, const rg_INT& serialFromInputFile, const BondType& typeOfBond, Molecule& aMolecule );
    
        void setAmberAtomTypeForAtomsInMolecule( Molecule& aMolecule );
        AmberAtomTypes getAmberAtomTypeForCarbon( Atom* carbon );
        AmberAtomTypes getAmberAtomTypeForNitrogen( Atom* nitrogen );
        AmberAtomTypes getAmberAtomTypeForOxygen( Atom* oxygen );
        AmberAtomTypes getAmberAtomTypeForHydrogen( Atom* hydrogen );

    void addRecordLinesOfFIleToList( ifstream* fin, list<string>& recordLinesOfFile );

    rg_FLAG writerCTMolFile( const char* pathFile, Molecule* aMolecule );


    rg_FLAG readChemOfficeChargeFile( const char* pathFile, Molecule& aMolecule );
        void    extractChargeRecordsFromRecLine( const string& strRecLine, string& atomType, rg_REAL& chargeValue, rg_INT& atomSerial );
        Atom*   getpAtomFromMap( AtomMap* mapOfAtom, const rg_INT& atomID );

    void mergeTwoMolecules( Molecule* molecule_A, Molecule* molecule_B, Molecule* mergedMolecule );


    
    void    initiallizeSymbolMapsOfAtomAndResidueForPDBFile( AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbolForProtein, ResidueSymbolMap& mapOfResidueSymbolForDNA, ResidueSymbolMap& mapOfResidueSymbolForRNA );
    void    initiallizeSymbolMapsOfAtomAndResidueForPDBFile( AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol ); // ONLY WORKS FOR PROTEINS. NOT DNA/RNA.
    void    initiallizeSymbolMapsOfAtomAndResidueForMol2File( AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol );
    void    initiallizeSymbolMapOfAtomForMolFile( AtomSymbolMap& mapOfAtomSymbol );

    rg_FLAG isValidMoleculeFileExtension( const char* pathFile, const MoleculeFileType& aMoleculeFileType );
    rg_BOOL checkAndReportErrorCodeForRecordType( const string& aRecType, const rg_FLAG& isRecLineOK, rg_FLAG& isFileRead );



    // .cfg via RMC simulation using the software RMCA
    rg_FLAG readRMCCFGFile( const string& pathFile, Molecule& aMolecule );
        void    readHeaderBlockOfRMCCFGFile( list<string>&              recordsOfMolFile, 
                                             list<string>::iterator&    i_record,
                                             vector<string>&            atomName,
                                             vector<float>&             atomCapacityInMolecule,
                                             int&                       totalNumAtoms,
                                             vector<int>&               numAtomsOfEachType,
                                             vector<rg_Point3D>&        definingVecotr);
        void    parseAtomNameAndCapacityInMolecule( const string&              record, 
                                             vector<string>&            atomName,
                                             vector<float>&             atomCapacityInMolecule);

    //  .dat for quasi-crystal data from Wenge Yang
    rg_FLAG readQuasiCrystalDataFile(const string& pathFile, Molecule& aMolecule);



    // PDB IO Functions
    rg_FLAG readPDBQFile(const string& pathFile, rg_dList<Molecule>& molecules);
        rg_FLAG setAtomRecordsToMoleculeForPDBQ(const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidue, ResidueMap& mapOfNonStdResidue, ChainMap& mapOfChain, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule);

    rg_FLAG writePDBQFile(const string& pathFile, const Molecule& aMolecule, const list<string>& frontRemarks);
    void    writePDBQFile(ofstream& fout, const Molecule& aMolecule, const list<string>& frontRemarks);
        void    writeHeaderOfPDBQFile(ofstream& fout, const Molecule& aMolecule, const list<string>& frontRemarks);
        void    writeAtomsOfPDBQFile(ofstream& fout, const Molecule& aMolecule);

};



} // namespace GeometryTier
} // namespace V


#endif

