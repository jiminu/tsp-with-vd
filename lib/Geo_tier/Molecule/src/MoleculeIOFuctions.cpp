#include "MoleculeIOFunctions.h"
#include "StringFunctions.h"
#include "SecondaryStructure.h"
#include "rg_TMatrix3D.h"
using namespace V::GeometryTier;


#include "FileOperations.h"
#include "FunctionsForMolecule.h"
#include <list>
using namespace std;



rg_FLAG MoleculeIOFunctions::readPDBFile( const char* pathFile, rg_dList<Molecule>& molecules )
{
    rg_FLAG isFileRead   = rg_FALSE;

    if( !isValidMoleculeFileExtension( pathFile, PDB_FILE ) ) {
        cerr << "Error: Invalid file extension !\n";
        return isFileRead;
    }
    
    
    ifstream fin( pathFile );
    
    if( fin.bad() ) {
        cerr << "Error: Could not open file!\n";
        return isFileRead;
    }
    
    list<string>  recordLinesOfPDBFile;
    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    addRecordsOfPDBFIleToList( &fin, recordLinesOfPDBFile );
    rg_INT numOfRecLines = recordLinesOfPDBFile.size();
    fin.close();

    string strPathFile = pathFile;
    rg_INT fileSize    = FileOperations::getFileSize( strPathFile );
    string depDate     = parseDepositionDate( recordLinesOfPDBFile );



    AtomSymbolMap    mapOfAtomSymbol;
    //ResidueSymbolMap mapOfResidueSymbol;
    ResidueSymbolMap mapOfResidueSymbolForProtein;
    ResidueSymbolMap mapOfResidueSymbolForDNA;
    ResidueSymbolMap mapOfResidueSymbolForRNA;
    
    //initiallizeSymbolMapsOfAtomAndResidueForPDBFile( mapOfAtomSymbol, mapOfResidueSymbol );
    initiallizeSymbolMapsOfAtomAndResidueForPDBFile( mapOfAtomSymbol,          mapOfResidueSymbolForProtein, 
                                                     mapOfResidueSymbolForDNA, mapOfResidueSymbolForRNA     );
    
    
    list<string*> connectRecLinesOfPDBFile;
    
    // Secondary Structure Section 
    list<string*> helixRecLinesOfPDBFile;
    list<string*> sheetRecLinesOfPDBFile;
    list<string*> turnRecLinesOfPDBFile;

    
    AtomMap          mapOfAtom;
    //ResidueMap       mapOfResidue;
    ResidueMap       mapOfNonStdResidue;    // ?
    
    rg_dList<ResidueMapWithChainID> mapsOfResidue;
    rg_dList<ResidueMapWithChainID> mapsOfNonStdResidue;
    
    ChainMap         mapOfChain;
    ChainMap         mapOfChainForHeteroRec;    // first: 0 for HOH, and 1~ for others.
    ChemBondMap      mapOfChemBond;
    
    rg_INT    countModelRec = 1;
    Molecule* aMolecule     = molecules.addTail( Molecule() );
    aMolecule->setModelSerialFromInputFile( countModelRec );
    aMolecule->setTimeStamp( depDate );
    aMolecule->setFileSize(  fileSize );
    
    list<string>::iterator i_recLines = recordLinesOfPDBFile.begin();

    while( i_recLines != recordLinesOfPDBFile.end() ) {
        string* currRecLine = &(*i_recLines);

        string recType = StringFunctions::subString( *currRecLine, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );


        if( recType == PDB_RECORD_TYPE_ATOM ) {
            //rg_FLAG isRecLineOK = setAtomRecordsToMolecule( *currRecLine, mapOfAtom, mapsOfResidue, mapOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbol, aMolecule );
            rg_FLAG isRecLineOK = setAtomRecordsToMolecule( *currRecLine, mapOfAtom, mapsOfResidue, mapOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbolForProtein, aMolecule );
            if( !checkAndReportErrorCodeForRecordType( recType, isRecLineOK, isFileRead ) ) {
                break;
            }
        }
        else if( recType == PDB_RECORD_TYPE_HETATM ) {
            //rg_FLAG isRecLineOK = setAtomRecordsToMolecule( *currRecLine, mapOfAtom, mapsOfResidue, mapOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbol, aMolecule );
            //rg_FLAG isRecLineOK = setAtomRecordsToMolecule( *currRecLine, mapOfAtom, mapsOfResidue, mapOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbolForProtein, aMolecule );
            
            //// setHeteroAtomRecordsToMolecule : TO BE COMPLETED
            rg_FLAG isRecLineOK = setAtomRecordsToMolecule( *currRecLine, mapOfAtom, mapsOfResidue, mapsOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbolForProtein, aMolecule );
            //rg_FLAG isRecLineOK = setHeteroAtomRecordsToMolecule( *currRecLine, mapOfAtom, mapsOfResidueForHeteroRec, mapOfChainForHeteroRec, mapOfAtomSymbol, mapOfResidueSymbolForProtein, aMolecule );
            if( !checkAndReportErrorCodeForRecordType( recType, isRecLineOK, isFileRead ) ) {
                break;
            }
        }
        else if( recType == PDB_RECORD_TYPE_CONECT ) {
            connectRecLinesOfPDBFile.push_back( &(*currRecLine) );
        }
        else if( recType == PDB_RECORD_TYPE_HELIX ) {
            helixRecLinesOfPDBFile.push_back( currRecLine );
        }
        else if( recType == PDB_RECORD_TYPE_SHEET ) {
            sheetRecLinesOfPDBFile.push_back( currRecLine );
        }
        else if( recType == PDB_RECORD_TYPE_TURN ) {
            turnRecLinesOfPDBFile.push_back( currRecLine );
        }
        else if( recType == PDB_RECORD_TYPE_MODEL ) {
            rg_INT modelSerial = atoi( StringFunctions::subString( *currRecLine, PDB_MODEL_RECORD_SERIAL_ST_POS, PDB_MODEL_RECORD_SERIAL_LENGTH ).c_str() );
            if( modelSerial != 0 ) {
                aMolecule->setModelSerialFromInputFile( modelSerial );
            }
        }    
        else if( recType == PDB_RECORD_TYPE_ENDMDL ) {
            countModelRec++;
            aMolecule = molecules.addTail( Molecule() );
            aMolecule->setModelSerialFromInputFile( countModelRec );
            aMolecule->setTimeStamp( depDate );
            aMolecule->setFileSize(  fileSize );

            mapOfAtom.clear();
            //mapOfResidue.clear();             //// WRONG SOURCE CODE !!! NOT WORKING WITH SEVERAL MODELS...
            //mapOfNonStdResidue.clear();       //// WRONG SOURCE CODE !!! NOT WORKING WITH SEVERAL MODELS...
            mapsOfResidue.removeAll();
            mapOfChain.clear();      
        }
		else if(recType == PDB_RECORD_TYPE_EXPDTA) {
			string method   = currRecLine->substr(10, currRecLine->length());
			method = StringFunctions::strTrimRight(method);
			aMolecule->setMethod(method);
		}
		else if(recType == PDB_RECORD_TYPE_REMARK) {
			string PDB_RECORD_RESOLUTION = "REMARK   2 RESOLUTION";

			string rec = StringFunctions::subString( *currRecLine, 0, 21 );

			if(rec == PDB_RECORD_RESOLUTION) {
				string resStr   = currRecLine->substr(22, currRecLine->length());		
				resStr = StringFunctions::strTrimLeft(resStr);
				
				const char* seps = " \t\n";
				char*	    token = NULL;

				char tempRes[200];
				strcpy(tempRes, resStr.c_str());
				token = strtok(tempRes, seps);   

				aMolecule->setResolution(atof(token));
			}			
		}
		else if(recType == PDB_RECORD_TYPE_CRYST) {
			aMolecule->setCryst(true);
			
			string crystStr   = currRecLine->substr(9, currRecLine->length());		
			crystStr = StringFunctions::strTrimLeft(crystStr);

			const char* seps = " \t\n";
			char*	    token = NULL;

			char tempCryst[200];
			strcpy(tempCryst, crystStr.c_str());
			token = strtok(tempCryst, seps);   

			for(int i = 0; i < 6; i++) {
				aMolecule->setCrystInfo(i, atof(token));
				token = strtok(NULL, seps);   
			}
		}
        
        i_recLines++;
    }

    //// WRONG SOURCE CODE !!! NOT WORKING WITH SEVERAL MODELS...
    //

    //// POST PROCESSES: KILL EMPTY MOLECULE MODEL, SET CONNECT RECORDS, AND ETC...
    //

    if( isFileRead ) {

        aMolecule = rg_NULL;
        molecules.reset4Loop();
        while ( molecules.setNext4Loop() ) {
            aMolecule = molecules.getpEntity();
        
            // KILL EMPTY MOLECULE MODEL
            if( aMolecule->getAtoms()->getSize() == 0 ) {
                molecules.killCurrent();
                continue;
            }
        

            // SET CONNECT RECORDS           
            list<string*>::iterator j_recLines = connectRecLinesOfPDBFile.begin();
            while( j_recLines != connectRecLinesOfPDBFile.end() ) {
                string* currRecLine = *j_recLines;

                rg_FLAG isRecLineOK = setConectRecordsToMolecule( *currRecLine, mapOfAtom, aMolecule );
                if( !checkAndReportErrorCodeForRecordType( PDB_RECORD_TYPE_CONECT, isRecLineOK, isFileRead ) ) {
                    break;
                }
                j_recLines++;        
            }
            setChemicalBondsToMoleculeForStandardResidues( *aMolecule );


            // SET CHAINS FOR NON STD RESIDUES            
            //setChainsForNonStdResiduesForPDBFile( &mapsOfNonStdResidue, *aMolecule );
            
            evaluateAndSetChainCode( *aMolecule );
            evaluateAndSetResidueCodesForNeucleicAcid( &mapOfResidueSymbolForDNA, &mapOfResidueSymbolForRNA, *aMolecule );


            // SET SECONDARY STRUCTURES
            // mapsOfResidue : HETATM is not considered for Chain.
            //setHelixRecordsToMolecule( &helixRecLinesOfPDBFile, &mapsOfResidue, aMolecule );
            //setSheetRecordsToMolecule( &sheetRecLinesOfPDBFile, &mapsOfResidue, aMolecule );
            //setTurnRecordsToMolecule( &turnRecLinesOfPDBFile, &mapsOfResidue, aMolecule );

            // ETC...
            string moleculeFileName = StringFunctions::getFileNameWithoutPath( string(pathFile) );
            aMolecule->setMoleculeFileName( moleculeFileName );

            aMolecule->evaluateHydrogenDonorAndAcceptorAtoms();
            aMolecule->computeAndSetCenterOfMass();
            aMolecule->computeAndSetMinEnclosingSphere();


			
			
        }


		
    }
    //
    ////

    return isFileRead;
}



rg_FLAG MoleculeIOFunctions::readPDBFile( const rg_INT& targetModelID, const char* pathFile, Molecule& aMolecule )
{
    rg_dList<Molecule> molecules;

    rg_FLAG isFileRead = readPDBFile( pathFile, molecules );

    molecules.reset4Loop();
    while ( molecules.setNext4Loop() ) {
        if( molecules.getpEntity()->getModelSerialFromInputFile() == targetModelID ) {
            aMolecule = molecules.getEntity();
            //aMolecule.syncPtrsWithSourceMolecule( molecules.getpEntity() );
            break;
        }
    }
    
    return isFileRead;
}



string  MoleculeIOFunctions::parseDepositionDate(const list<string>& recordsOfPDBFile)
{
    string depDate;
    list<string>::const_iterator i_recLines = recordsOfPDBFile.begin();

    while( i_recLines != recordsOfPDBFile.end() ) {
        string currRecLine = (*i_recLines);
        string recType     = StringFunctions::subString( currRecLine, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );

        if( recType == PDB_RECORD_TYPE_HEADER ) {
            depDate = StringFunctions::subString( currRecLine, PDB_HEADER_RECORD_DATE_ST_POS, PDB_HEADER_RECORD_DATE_LENGTH );
            break;
        }
		i_recLines++;
    }

    return depDate;
}



void MoleculeIOFunctions::addRecordsOfPDBFIleToList( ifstream* fin, list<string>& recordsOfPDBFile )
{
    string strRecLine  = "";
    
    while ( getline( *fin, strRecLine ) ) {
        recordsOfPDBFile.push_back( strRecLine );
    }
}



rg_FLAG MoleculeIOFunctions::setAtomRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidue, ResidueMap& mapOfNonStdResidue, ChainMap& mapOfChain, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule )
{

    rg_INT atomSerial = atoi( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_SERIAL_ST_POS, PDB_ATOM_RECORD_SERIAL_LENGTH ).c_str() );
    string atomName   = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_NAME_ST_POS, PDB_ATOM_RECORD_NAME_LENGTH );
    string altLoc     = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ALTLOC_ST_POS, PDB_ATOM_RECORD_ALTLOC_LENGTH );   // NOT USING
    string resName    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_RESNAME_ST_POS, PDB_ATOM_RECORD_RESNAME_LENGTH );
    string chainID    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_CHAINID_ST_POS, PDB_ATOM_RECORD_CHAINID_LENGTH );
    rg_INT resSeq     = atoi( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_RESSEQ_ST_POS, PDB_ATOM_RECORD_RESSEQ_LENGTH ).c_str() );
    string iCode      = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ICODE_ST_POS, PDB_ATOM_RECORD_ICODE_LENGTH );     // NOT USING
    float  x_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH ).c_str() );
    float  y_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH ).c_str() );
    float  z_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH ).c_str() );
    //rg_REAL  x_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH ).c_str() );
    //rg_REAL  y_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH ).c_str() );
    //rg_REAL  z_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH ).c_str() );
    float  occupancy  = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_OCCUPANCY_ST_POS, PDB_ATOM_RECORD_OCCUPANCY_LENGTH  ).c_str() );
    float  tempFactor = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_TEMPFACTOR_ST_POS, PDB_ATOM_RECORD_TEMPFACTOR_LENGTH ).c_str() );
    string segID      = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_SEGID_ST_POS, PDB_ATOM_RECORD_SEGID_LENGTH );     // NOT USING
    string element    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ELEMENT_ST_POS, PDB_ATOM_RECORD_ELEMENT_LENGTH ); // NOT USING
    string charge     = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_CHARGE_ST_POS, PDB_ATOM_RECORD_CHARGE_LENGTH );   // NOT USING

    
    if( altLoc == " " || altLoc == "A" ) {

        // Estimate Chain
        Chain* currChain = rg_NULL;
        rg_INT numOfChains = mapOfChain.size();
        rg_INT intChainID  = (int)chainID[0];

        

        if ( intChainID == 66 && resSeq == 121 )
        {
            int aaa =-0;
        }




        ResidueMapWithChainID* pResidueMapWithChainID = rg_NULL;

        ChainMap::iterator chainMap_i = mapOfChain.find( intChainID );
        
        if ( chainMap_i == mapOfChain.end() ) {
            currChain = aMolecule->addChain( Chain(numOfChains, aMolecule) );
            currChain->setChainIDFromInputFileInDecimal( intChainID );
            mapOfChain.insert( ChainMap::value_type(intChainID, currChain) );

            ResidueMapWithChainID tempResidueMapWithChainID;
            tempResidueMapWithChainID.first = numOfChains;
            pResidueMapWithChainID = mapsOfResidue.addTail( tempResidueMapWithChainID );
        }
        else {
            currChain = (*chainMap_i).second;

            mapsOfResidue.reset4Loop();
            while ( mapsOfResidue.setNext4Loop() ) {
                pResidueMapWithChainID = mapsOfResidue.getpEntity();

                if( pResidueMapWithChainID->first == currChain->getID() )
                    break;
            }
        }

        ResidueMap* aMapOfResidue = &(pResidueMapWithChainID->second);
        


        // Estimate Residue
        Residue* currResidue  = rg_NULL;
        rg_INT   numOfResidue = 0;
        mapsOfResidue.reset4Loop();
        while ( mapsOfResidue.setNext4Loop() ) {
            numOfResidue += mapsOfResidue.getpEntity()->second.size();
        }
        
        ResidueMap::iterator residueMap_i = aMapOfResidue->find( resSeq );
        
        if ( residueMap_i == aMapOfResidue->end() ) {
            currResidue = aMolecule->addResidue( Residue(numOfResidue) );
            currChain->addResidue( currResidue );
            currResidue->setChain( currChain );
            currResidue->setSequenceNumber( resSeq );
            currResidue->setResidueName( resName );
            setResidueCodeToTargetResidueExceptNeucleicAcid( mapOfResidueSymbol, resName, currResidue );  // NOT WORKING WITH DNA/RNA.
            aMapOfResidue->insert( ResidueMap::value_type(resSeq, currResidue) );
            
            if( !currResidue->isStandardResidue() ) {
                mapOfNonStdResidue.insert( ResidueMap::value_type( resSeq, currResidue ) );
            }
        }
        else {
            currResidue = (*residueMap_i).second;
        }

        //// Checking Oxygen in OXT row for duplication.
        //
        rg_BOOL isOXTDupplicated = rg_FALSE;
        
        if ( atomName == " OXT" ) {
            rg_dList<Atom*>* atomsInResidue = currResidue->getAtoms();
            atomsInResidue->reset4Loop();
            while ( atomsInResidue->setNext4Loop() ) {
                Atom* currAtom = atomsInResidue->getEntity();

                if ( currAtom->getAtomNameFromInputFile() == " O  " ) {
                    isOXTDupplicated = rg_TRUE;
                    break;
                }
            }
        }

        if ( isOXTDupplicated == rg_TRUE )
            return rg_TRUE;
        //
        ////


        // Create Atom : (1)Serial from inputfile, (2)AtomCode, (3)RemotenessIndicator, (4)BranchDesignagtor(include (5)ext.)
        //               (6)AtomTypeInAmber, (7)charge and (8)isAtomOnBackBone are considered.
        rg_INT numOfAtom = mapOfAtom.size();
        
        Atom* currAtom = aMolecule->addAtom( Atom(numOfAtom) );
        currResidue->addAtom( currAtom );
        currAtom->setResidue( currResidue );
        currAtom->setSerialFromInputFile( atomSerial );     //  (1)
        currAtom->setAtomNameFromInputFile( atomName );
        rg_FLAG isAtomNameOK = setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( mapOfAtomSymbol, atomName, currAtom ); // (2)~(8)

        // Hydrogen atom code verification
        if(element == string(" H"))
            currAtom->setAtomCode( H_ATOM );
        
        if( isAtomNameOK == rg_FALSE )
            return rg_FALSE;
        
        // set coordinates, radius, occupancy, and tempFactor.
        currAtom->setAtomBall( Sphere( x_coord, y_coord, z_coord, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );
        currAtom->getpChemicalProperties()->setOccupancy( occupancy );
        currAtom->getpChemicalProperties()->setTempFactor( tempFactor );
        
        mapOfAtom.insert( AtomMap::value_type(atomSerial, currAtom) );
    }
    
    return rg_TRUE;
}



rg_FLAG MoleculeIOFunctions::setAtomRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidue, rg_dList<ResidueMapWithChainID>& mapsOfNonStdResidue, ChainMap& mapOfChain, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule )
{
    rg_INT atomSerial = atoi( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_SERIAL_ST_POS, PDB_ATOM_RECORD_SERIAL_LENGTH ).c_str() );
    string atomName   = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_NAME_ST_POS, PDB_ATOM_RECORD_NAME_LENGTH );
    string altLoc     = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ALTLOC_ST_POS, PDB_ATOM_RECORD_ALTLOC_LENGTH );   // NOT USING
    string resName    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_RESNAME_ST_POS, PDB_ATOM_RECORD_RESNAME_LENGTH );
    string chainID    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_CHAINID_ST_POS, PDB_ATOM_RECORD_CHAINID_LENGTH );
    rg_INT resSeq     = atoi( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_RESSEQ_ST_POS, PDB_ATOM_RECORD_RESSEQ_LENGTH ).c_str() );
    string iCode      = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ICODE_ST_POS, PDB_ATOM_RECORD_ICODE_LENGTH );     // NOT USING
    rg_REAL  x_coord  = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH ).c_str() );
    rg_REAL  y_coord  = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH ).c_str() );
    rg_REAL  z_coord  = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH ).c_str() );
    //float  x_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH ).c_str() );
    //float  y_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH ).c_str() );
    //float  z_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH ).c_str() );
    float  occupancy  = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_OCCUPANCY_ST_POS, PDB_ATOM_RECORD_OCCUPANCY_LENGTH  ).c_str() );
    float  tempFactor = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_TEMPFACTOR_ST_POS, PDB_ATOM_RECORD_TEMPFACTOR_LENGTH ).c_str() );
    string segID      = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_SEGID_ST_POS, PDB_ATOM_RECORD_SEGID_LENGTH );     // NOT USING
    string element    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ELEMENT_ST_POS, PDB_ATOM_RECORD_ELEMENT_LENGTH ); // NOT USING
    string charge     = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_CHARGE_ST_POS, PDB_ATOM_RECORD_CHARGE_LENGTH );   // NOT USING

    
    if( altLoc == " " || altLoc == "A" ) {

        // Estimate Chain
        Chain* currChain = rg_NULL;
        rg_INT numOfChains = mapOfChain.size();
        rg_INT intChainID  = (int)chainID[0];

        
        ResidueMapWithChainID* pResidueMapWithChainID       = rg_NULL;
        ResidueMapWithChainID* pNonStdResidueMapWithChainID = rg_NULL;

        ChainMap::iterator chainMap_i = mapOfChain.find( intChainID );
        
        if ( chainMap_i == mapOfChain.end() ) {
            currChain = aMolecule->addChain( Chain(numOfChains, aMolecule) );
            currChain->setChainIDFromInputFileInDecimal( intChainID );
            mapOfChain.insert( ChainMap::value_type(intChainID, currChain) );

            ResidueMapWithChainID tempResidueMapWithChainID;
            tempResidueMapWithChainID.first = numOfChains;
            pResidueMapWithChainID = mapsOfResidue.addTail( tempResidueMapWithChainID );

            //// HETATM records...
            //
            string recType = StringFunctions::subString( strRecLine, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );
            if ( recType == PDB_RECORD_TYPE_HETATM ) {
                ResidueMapWithChainID tempNonStdResidueMapWithChainID;
                tempNonStdResidueMapWithChainID.first = numOfChains;
                pNonStdResidueMapWithChainID = mapsOfNonStdResidue.addTail( tempNonStdResidueMapWithChainID );
            }
            //
            ////
        }
        else {
            currChain = (*chainMap_i).second;

            mapsOfResidue.reset4Loop();
            while ( mapsOfResidue.setNext4Loop() ) {
                pResidueMapWithChainID = mapsOfResidue.getpEntity();

                if( pResidueMapWithChainID->first == currChain->getID() )
                    break;
            }

            //// HETATM records...
            //
            string recType = StringFunctions::subString( strRecLine, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );
            if ( recType == PDB_RECORD_TYPE_HETATM ) {

                if ( mapsOfNonStdResidue.getSize() == 0 ) {
                    ResidueMapWithChainID tempNonStdResidueMapWithChainID;
                    tempNonStdResidueMapWithChainID.first = numOfChains;
                    pNonStdResidueMapWithChainID = mapsOfNonStdResidue.addTail( tempNonStdResidueMapWithChainID );
                }
                else {
                    ResidueMapWithChainID tempNonStdResidueMapWithChainID;
                    tempNonStdResidueMapWithChainID.first = numOfChains;
                    pNonStdResidueMapWithChainID = mapsOfNonStdResidue.addTail( tempNonStdResidueMapWithChainID );
                    mapsOfNonStdResidue.reset4Loop();
                    while ( mapsOfNonStdResidue.setNext4Loop() ) {
                        pNonStdResidueMapWithChainID = mapsOfNonStdResidue.getpEntity();

                        if( pNonStdResidueMapWithChainID->first == currChain->getID() )
                            break;
                    }
                }
            }
            //
            ////
        }

        ResidueMap* aMapOfResidue = &(pResidueMapWithChainID->second);
        


        // Estimate Residue
        Residue* currResidue  = rg_NULL;
        rg_INT   numOfResidue = 0;
        mapsOfResidue.reset4Loop();
        while ( mapsOfResidue.setNext4Loop() ) {
            numOfResidue += mapsOfResidue.getpEntity()->second.size();
        }
        
        ResidueMap::iterator residueMap_i = aMapOfResidue->find( resSeq );
        
        if ( residueMap_i == aMapOfResidue->end() ) {
            currResidue = aMolecule->addResidue( Residue(numOfResidue) );
            currChain->addResidue( currResidue );
            currResidue->setChain( currChain );
            currResidue->setSequenceNumber( resSeq );
            currResidue->setResidueName( resName );
            setResidueCodeToTargetResidueExceptNeucleicAcid( mapOfResidueSymbol, resName, currResidue );  // NOT WORKING WITH DNA/RNA.
            aMapOfResidue->insert( ResidueMap::value_type(resSeq, currResidue) );
            
            if( !currResidue->isStandardResidue() ) {
                ResidueMap* aMapOfNonStdResidue = &(pNonStdResidueMapWithChainID->second);
                aMapOfNonStdResidue->insert( ResidueMap::value_type( resSeq, currResidue ) );
            }
        }
        else {
            currResidue = (*residueMap_i).second;
        }

        // 현재는 마지막 residue에 oxygen이 존재하지 않는 경우에만 OXT를 읽도록 되어 있음.
        // 존재할 경우에는 읽지 않음.
        //// Checking Oxygen in OXT row for duplication.
        //
        rg_BOOL isOXTDupplicated = rg_FALSE;
        
        if ( atomName == " OXT" ) {
            rg_dList<Atom*>* atomsInResidue = currResidue->getAtoms();
            atomsInResidue->reset4Loop();
            while ( atomsInResidue->setNext4Loop() ) {
                Atom* currAtom = atomsInResidue->getEntity();

                if ( currAtom->getAtomNameFromInputFile() == " O  " ) {
                    isOXTDupplicated = rg_TRUE;
                    break;
                }
            }
        }

        if ( isOXTDupplicated == rg_TRUE )
            return rg_TRUE;
        //
        ////


        // Create Atom : (1)Serial from inputfile, (2)AtomCode, (3)RemotenessIndicator, (4)BranchDesignagtor(include (5)ext.)
        //               (6)AtomTypeInAmber, (7)charge and (8)isAtomOnBackBone are considered.
        rg_INT numOfAtom = mapOfAtom.size();
        
        Atom* currAtom = aMolecule->addAtom( Atom(numOfAtom) );
        currResidue->addAtom( currAtom );
        currAtom->setResidue( currResidue );
        currAtom->setSerialFromInputFile( atomSerial );     //  (1)
        currAtom->setAtomNameFromInputFile( atomName );
        rg_FLAG isAtomNameOK = setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( mapOfAtomSymbol, atomName, currAtom ); // (2)~(8)
        
        // Hydrogen atom code verification
        if(element == string(" H"))
            currAtom->setAtomCode( H_ATOM );

        if( isAtomNameOK == rg_FALSE )
            return rg_FALSE;
        
        // set coordinates, radius, occupancy, and tempFactor.
        currAtom->setAtomBall( Sphere( x_coord, y_coord, z_coord, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );
        currAtom->getpChemicalProperties()->setOccupancy( occupancy );
        currAtom->getpChemicalProperties()->setTempFactor( tempFactor );
        
        mapOfAtom.insert( AtomMap::value_type(atomSerial, currAtom) );
    }
    
    return rg_TRUE;
}



void MoleculeIOFunctions::setResidueCodeToTargetResidueExceptNeucleicAcid( ResidueSymbolMap& mapOfResidueSymbol, const string& threeCodeNameOfRes, Residue* targetResidue )
{   
    ResidueCode convertedResidueCode = UNK_RESIDUE;

    if( threeCodeNameOfRes.size() == 3 ) {
        ResidueSymbolMap::iterator ResidueSymbolMap_i = mapOfResidueSymbol.find( threeCodeNameOfRes );
        
        if ( ResidueSymbolMap_i == mapOfResidueSymbol.end() ) {
            convertedResidueCode = UNK_RESIDUE;
        }
        else {
            convertedResidueCode = (ResidueCode)((*ResidueSymbolMap_i).second);
        }
    }
    
    targetResidue->setResidueCode( convertedResidueCode );
}



void MoleculeIOFunctions::setResidueCodeToTargetResidue( ResidueSymbolMap& mapOfResidueSymbol, const string& threeCodeNameOfRes, Residue* targetResidue )
{
    ResidueCode convertedResidueCode = UNK_RESIDUE;

    if( threeCodeNameOfRes.size() == 3 ) {
        ResidueSymbolMap::iterator ResidueSymbolMap_i = mapOfResidueSymbol.find( threeCodeNameOfRes );
        
        if ( ResidueSymbolMap_i == mapOfResidueSymbol.end() ) {
            convertedResidueCode = UNK_RESIDUE;
        }
        else {
            convertedResidueCode = (ResidueCode)((*ResidueSymbolMap_i).second);
        }
    }
    
    targetResidue->setResidueCode( convertedResidueCode );
}



rg_FLAG MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( AtomSymbolMap& mapOfAtomSymbol, const string& atomName, Atom* targetAtom )
{
    rg_FLAG isAtomNameOK = rg_TRUE;

////   SOME PDB FILES(ex: 1g9v) DO NOT FOLLOW THE CONDITIONS ON BELOW....
//
//     // Check if remoteness indicator is in numeric letter : remoteIndicator must be in alphabet letter !
//     if( StringFunctions::isStringLetterInNumber( atomName.substr( 2, 1 ) ) == true ) {
//         isAtomNameOK = rg_FALSE;
//     }
//     
//     // Check if brangeDesignator is in Alphabet : remoteIndicator must be in Numeric letter !
//     if( StringFunctions::isStringLetterInAlphabet( atomName.substr( 3, 1 ) ) == true ) {
//         isAtomNameOK = rg_FALSE;
//     }

    if( isAtomNameOK == rg_TRUE ) {

        // SET ATOM CODE
		string chemicalSymbol = atomName.substr( 0, 2 );
		AtomCode convertedAtomCode = UNK_ATOM;
		if(! doesStringContainHydrogenAtomName(atomName, targetAtom))
		{
			AtomSymbolMap::iterator AtomSymbolMap_i = mapOfAtomSymbol.find( chemicalSymbol );
			
			if ( AtomSymbolMap_i != mapOfAtomSymbol.end() ) {
				convertedAtomCode = (AtomCode)((*AtomSymbolMap_i).second);
			}
		}
		else
		{
			convertedAtomCode = H_ATOM;
		}

        targetAtom->setAtomCode( convertedAtomCode );
        


        ////////////////////////////////////////////////////////////////////////////////
        //
        string atomType         = StringFunctions::strTrim(ATOM_FEATURES[convertedAtomCode].symbol);
        rg_INT lengthOfAtomName = StringFunctions::strTrim(atomName).length();

        RemoteIndicator  remoteIndicator = UNK_REMOTE;
        BranchDesignator branchDesignator = UNK_BRANCH;
        BranchDesignator extraBranchDesignator = UNK_BRANCH;


        if (lengthOfAtomName > 0 && lengthOfAtomName < 5) {

            if (targetAtom->getResidue()->isStandardResidue()) {

                //string strRemoteAndBranch = StringFunctions::subString(atomName, atomType.length(), atomName.length());
                string::size_type pos = atomName.find(atomType);
                string::size_type len = atomType.length();
                string strRemoteAndBranch = StringFunctions::subString(atomName, pos + len, atomName.length());


                rg_INT lengthOfStrRemoteAndBranch = strRemoteAndBranch.length();

                switch (lengthOfStrRemoteAndBranch) {
                case 0: {
                    remoteIndicator = UNK_REMOTE;
                    branchDesignator = UNK_BRANCH;
                    if (pos == 1) {
                        extraBranchDesignator = convertStringIntoBranchDesignator(atomName.substr(0, 1));
                    }
                    else {
                        extraBranchDesignator = UNK_BRANCH;
                    }
                    break;
                }
                case 1: {
                    string strRemoteIndicator = strRemoteAndBranch.substr(0, 1);

                    if (isStringValidForRemoteIndicator(strRemoteIndicator)) {           // REMOTENESS IN CHAR
                        remoteIndicator = convertStringIntoRemoteIndicator(strRemoteIndicator);
                        if (pos == 1) {
                            extraBranchDesignator = convertStringIntoBranchDesignator(atomName.substr(0, 1));
                        }
                        else {
                            extraBranchDesignator = UNK_BRANCH;
                        }
                    }
                    break;
                }
                case 2: {
                    string strRemoteIndicator = strRemoteAndBranch.substr(0, 1);
                    string strBranchDesignator = strRemoteAndBranch.substr(1, 1);

                    if (isStringValidForRemoteIndicator(strRemoteIndicator) &&          // REMOTENESS IN CHAR
                        isStringValidForBranchDesignator(strBranchDesignator)) {        // BRANCH IN NUM
                        remoteIndicator = convertStringIntoRemoteIndicator(strRemoteIndicator);
                        branchDesignator = convertStringIntoBranchDesignator(strBranchDesignator);

                        if (pos == 1) {
                            extraBranchDesignator = convertStringIntoBranchDesignator(atomName.substr(0, 1));
                        }
                        else {
                            extraBranchDesignator = UNK_BRANCH;
                        }
                    }
                    break;
                }
                case 3: {
                    string strRemoteIndicator = strRemoteAndBranch.substr(0, 1);
                    string strBranchDesignator = strRemoteAndBranch.substr(1, 1);
                    string strExtraBranchDesignator = strRemoteAndBranch.substr(2, 1);

                    if (isStringValidForRemoteIndicator(strRemoteIndicator) &&        // REMOTENESS IN CHAR
                        isStringValidForBranchDesignator(strBranchDesignator) &&        // BRANCH IN NUM
                        isStringValidForBranchDesignator(strExtraBranchDesignator)) {   // BRANCH IN NUM

                        remoteIndicator = convertStringIntoRemoteIndicator(strRemoteIndicator);
                        branchDesignator = convertStringIntoBranchDesignator(strBranchDesignator);
                        extraBranchDesignator = convertStringIntoBranchDesignator(strExtraBranchDesignator);
                    }
                    break;
                }
                default: {
                    remoteIndicator = UNK_REMOTE;
                    branchDesignator = UNK_BRANCH;
                    extraBranchDesignator = UNK_BRANCH;
                    break;

                }
                }
            }
        }

        //
        ////////////////////////////////////////////////////////////////////////////////


        targetAtom->getpChemicalProperties()->setRemoteIndicator(remoteIndicator);
        targetAtom->getpChemicalProperties()->setBrangeDesignator(branchDesignator);
        targetAtom->getpChemicalProperties()->setExtraBrangeDesignator(extraBranchDesignator);


        /*
        // set ChemicalProperties
        BranchDesignator extraBranchDesignator = UNK_BRANCH;
        RemoteIndicator  remoteIndicator       = convertStringIntoRemoteIndicator( StringFunctions::subString( atomName, 2, 1 ) );
        BranchDesignator branchDesignator      = convertStringIntoBranchDesignator( StringFunctions::subString( atomName, 3, 1 ) );

        // Hydrogen have extraBrangeDesignator in number
        string firstLetterInAtomName = atomName.substr( 0, 1 );
        if( StringFunctions::isStringLetterInNumber( firstLetterInAtomName ) == true ) {
            chemicalSymbol.replace(0, 1, " ");
            extraBrangeDesignator = convertStringIntoBranchDesignator( firstLetterInAtomName );
            targetAtom->getpChemicalProperties()->setExtraBrangeDesignator( extraBrangeDesignator );
        }

        targetAtom->getpChemicalProperties()->setRemoteIndicator( remoteIndicator );
        targetAtom->getpChemicalProperties()->setBrangeDesignator( brangeDesignator );
        */


        // Check ResidueCode of targetAtom for AtomTypeInAmber : AtomTypeInAmber works only for Amino acid !
        ResidueCode resCodeForTargetAtom = targetAtom->getResidue()->getResidueCode();
        if ( (rg_INT)resCodeForTargetAtom > 30 )
            resCodeForTargetAtom = UNK_RESIDUE;





        // Set AtomTypeInAmber
        targetAtom->getpChemicalProperties()->setAtomTypeInAmber( AMBER_ATOM_TYPE[resCodeForTargetAtom]
                                                                                 [ATOM_FEATURES[convertedAtomCode].IDOfAtomInStandardResidue ]
                                                                                 [remoteIndicator]
                                                                                 [branchDesignator]
                                                                                 [extraBranchDesignator] );

        targetAtom->getpChemicalProperties()->setChargeInAmber( ELECTROSTATIC_CHARGE_COEFFS_NON_BONDED_ATOM_PAIR[resCodeForTargetAtom]
                                                                                                                [ATOM_FEATURES[convertedAtomCode].IDOfAtomInStandardResidue ]
                                                                                                                [remoteIndicator]
                                                                                                                [branchDesignator]
                                                                                                                [extraBranchDesignator] );

        // Check isOnBackBone
        rg_FLAG isAtomOnBackBone = rg_FALSE;

        // if residue is Amino acid
        if( (rg_INT)targetAtom->getResidue()->getResidueCode() > 0 && (rg_INT)targetAtom->getResidue()->getResidueCode() < 31 ) {
            if( atomName == " N  " || atomName == " CA " || atomName == " C  " )
                isAtomOnBackBone = rg_TRUE;
        }
        // else if residue is DNA or RNA
        else if( (rg_INT)targetAtom->getResidue()->getResidueCode() > 30 && (rg_INT)targetAtom->getResidue()->getResidueCode() < 51 ) {
            if( atomName == " P  " || atomName == " O5*" || atomName == " C5*" || 
                atomName == " C4*" || atomName == " C3*" || atomName == " O3*"  )
                isAtomOnBackBone = rg_TRUE;
        }
        
        targetAtom->getpChemicalProperties()->setIsOnBackBone( isAtomOnBackBone );

    }

    return isAtomNameOK;
}



RemoteIndicator MoleculeIOFunctions::convertStringIntoRemoteIndicator( const string& strRemoteIndicator )
{
    RemoteIndicator aRemoteIndicator = UNK_REMOTE;

    switch ( strRemoteIndicator[0] )  {
        case 'A':
            aRemoteIndicator = ALPHA_REMOTE;
            break;
        case 'B':
            aRemoteIndicator = BETA_REMOTE;
            break;
        case 'G':
            aRemoteIndicator = GAMMA_REMOTE;
            break;
        case 'D':
            aRemoteIndicator = DELTA_REMOTE;
            break;
        case 'E':
            aRemoteIndicator = EPSILON_REMOTE;
            break;
        case 'Z':
            aRemoteIndicator = ZETA_REMOTE;
            break;
        case 'H':
            aRemoteIndicator = ETA_REMOTE;
            break;
        case 'X':
        case 'T':
            aRemoteIndicator = TERMINATE_REMOTE;
            break;
        default:
            aRemoteIndicator = UNK_REMOTE;
            break;
    }
    
    return aRemoteIndicator;
}



BranchDesignator MoleculeIOFunctions::convertStringIntoBranchDesignator( const string& strBranchDesignator )
{
    BranchDesignator aBranchDesignator = UNK_BRANCH;

    switch ( strBranchDesignator[0] )  {
        case '1':
            aBranchDesignator = FIRST_BRANCH;
            break;
        case '2':
            aBranchDesignator = SECOND_BRANCH;
            break;
        case '3':
            aBranchDesignator = THIRD_BRANCH;
            break;
        case 'T':
            aBranchDesignator = FOURTH_BRANCH;
            break;
        default:
            aBranchDesignator = UNK_BRANCH;
            break;
    }
    
    return aBranchDesignator;
}



rg_FLAG MoleculeIOFunctions::setHeteroAtomRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidueForHeteroRec, ChainMap& mapOfChainForHeteroRec, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule )
{
    rg_INT atomSerial = atoi( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_SERIAL_ST_POS, PDB_ATOM_RECORD_SERIAL_LENGTH ).c_str() );
    string atomName   = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_NAME_ST_POS, PDB_ATOM_RECORD_NAME_LENGTH );
    string altLoc     = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ALTLOC_ST_POS, PDB_ATOM_RECORD_ALTLOC_LENGTH );   // NOT USING
    string resName    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_RESNAME_ST_POS, PDB_ATOM_RECORD_RESNAME_LENGTH );
    string chainID    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_CHAINID_ST_POS, PDB_ATOM_RECORD_CHAINID_LENGTH );
    rg_INT resSeq     = atoi( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_RESSEQ_ST_POS, PDB_ATOM_RECORD_RESSEQ_LENGTH ).c_str() );
    string iCode      = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ICODE_ST_POS, PDB_ATOM_RECORD_ICODE_LENGTH );     // NOT USING
    float  x_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH ).c_str() );
    float  y_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH ).c_str() );
    float  z_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH ).c_str() );
    float  occupancy  = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_OCCUPANCY_ST_POS, PDB_ATOM_RECORD_OCCUPANCY_LENGTH  ).c_str() );
    float  tempFactor = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_TEMPFACTOR_ST_POS, PDB_ATOM_RECORD_TEMPFACTOR_LENGTH ).c_str() );
    string segID      = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_SEGID_ST_POS, PDB_ATOM_RECORD_SEGID_LENGTH );     // NOT USING
    string element    = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_ELEMENT_ST_POS, PDB_ATOM_RECORD_ELEMENT_LENGTH ); // NOT USING
    string charge     = StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_CHARGE_ST_POS, PDB_ATOM_RECORD_CHARGE_LENGTH );   // NOT USING

    
    if( altLoc == " " || altLoc == "A" ) {
        
        // Estimate Chain
        Chain* currChain  = rg_NULL;
        rg_INT intChainID = (int)chainID[0];
        rg_INT numOfChains = aMolecule->getChains()->getSize();

        if( resName.compare( "HOH" ) == 0  ) {

            ChainMap::iterator chainMap_i = mapOfChainForHeteroRec.find( 0 );           // 0 for HOH chain
        
            if ( chainMap_i == mapOfChainForHeteroRec.end() ) {
                currChain = aMolecule->addChain( Chain( numOfChains, aMolecule ) );
                currChain->setChainIDFromInputFileInDecimal( intChainID );               
                mapOfChainForHeteroRec.insert( ChainMap::value_type( 0, currChain ) );  // 0 for HOH chain
            }
            else {
                currChain = (*chainMap_i).second;
            }
        }
        else {

            ResidueMapWithChainID* pResidueMapWithChainID = rg_NULL;

            mapsOfResidueForHeteroRec.reset4Loop();
            while ( mapsOfResidueForHeteroRec.setNext4Loop() ) {
                pResidueMapWithChainID = mapsOfResidueForHeteroRec.getpEntity();

                if( pResidueMapWithChainID->first == intChainID )
                    break;
            }

            ResidueMap* aMapOfResidue = &(pResidueMapWithChainID->second);

//             aMapOfResidue->find(  )
// 
// 
//             ResidueMapWithChainID tempResidueMapWithChainID;
//             tempResidueMapWithChainID.first = numOfChains;
//             pResidueMapWithChainID = mapsOfResidue.addTail( tempResidueMapWithChainID );
// 
// 
// 
//             ChainMap::iterator chainMap_i = mapOfChainForHeteroRec.find( intChainID );  // 0 for HOH chain
//         
//             if ( chainMap_i == mapOfChainForHeteroRec.end() ) {
//                 currChain = aMolecule->addChain( Chain( numOfChains, aMolecule ) );
//                 currChain->setChainIDFromInputFileInDecimal( intChainID );
//                 mapOfChainForHeteroRec.insert( ChainMap::value_type( numOfChains, currChain ) );  // 0 for HOH chain
//             }
//             else {
//                 currChain = (*chainMap_i).second;
//             }


        }
        
//         // Estimate Residue
//         Residue* currResidue  = rg_NULL;
//         rg_INT   numOfResidue = 0;
//         mapsOfResidue.reset4Loop();
//         while ( mapsOfResidue.setNext4Loop() ) {
//             numOfResidue += mapsOfResidue.getpEntity()->second.size();
//         }
//         
//         ResidueMap::iterator residueMap_i = aMapOfResidue->find( resSeq );
//         
//         if ( residueMap_i == aMapOfResidue->end() ) {
//             currResidue = aMolecule->addResidue( Residue(numOfResidue) );
//             currChain->addResidue( currResidue );
//             currResidue->setChain( currChain );
//             currResidue->setSequenceNumber( resSeq );
//             currResidue->setResidueName( resName );
//             setResidueCodeToTargetResidueExceptNeucleicAcid( mapOfResidueSymbol, resName, currResidue );  // NOT WORKING WITH DNA/RNA.
//             aMapOfResidue->insert( ResidueMap::value_type(resSeq, currResidue) );
//             
//             if( !currResidue->isStandardResidue() ) {
//                 mapOfNonStdResidue.insert( ResidueMap::value_type( resSeq, currResidue ) );
//             }
//         }
//         else {
//             currResidue = (*residueMap_i).second;
//         }
// 
// 
//         // Create Atom : (1)Serial from inputfile, (2)AtomCode, (3)RemotenessIndicator, (4)BranchDesignagtor(include (5)ext.)
//         //               (6)AtomTypeInAmber, (7)charge and (8)isAtomOnBackBone are considered.
//         rg_INT numOfAtom = mapOfAtom.size();
//         
//         Atom* currAtom = aMolecule->addAtom( Atom(numOfAtom) );
//         currResidue->addAtom( currAtom );
//         currAtom->setResidue( currResidue );
//         currAtom->setSerialFromInputFile( atomSerial );     //  (1)
//         currAtom->setAtomNameFromInputFile( atomName );
//         rg_FLAG isAtomNameOK = setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( mapOfAtomSymbol, atomName, currAtom ); // (2)~(8)
//         
//         if( isAtomNameOK == rg_FALSE )
//             return rg_FALSE;
//         
//         // set coordinates, radius, occupancy, and tempFactor.
//         currAtom->setAtomBall( Sphere( x_coord, y_coord, z_coord, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );
//         currAtom->getpChemicalProperties()->setOccupancy( occupancy );
//         currAtom->getpChemicalProperties()->setTempFactor( tempFactor );
//         
//         mapOfAtom.insert( AtomMap::value_type(atomSerial, currAtom) );
    }
    
    return rg_TRUE;
}



rg_FLAG MoleculeIOFunctions::setConectRecordsToMolecule( list<string*>* connectRecLinesOfPDBFile, AtomMap& mapOfAtom, Molecule* aMolecule )
{
    rg_FLAG isRecLinesOK = rg_TRUE;
    rg_FLAG isFileRead   = rg_TRUE;

    list<string*>::iterator i_recLines = connectRecLinesOfPDBFile->begin();
    while( i_recLines != connectRecLinesOfPDBFile->end() ) {
        string* currRecLine = *i_recLines;

        rg_FLAG isRecLineOK = setConectRecordsToMolecule( *currRecLine, mapOfAtom, aMolecule );
        if( !checkAndReportErrorCodeForRecordType( PDB_RECORD_TYPE_CONECT, isRecLineOK, isFileRead ) ) {
            break;
        }
        i_recLines++;        
    }

    return isFileRead;
}



rg_FLAG MoleculeIOFunctions::setConectRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, Molecule* aMolecule )
{
    rg_FLAG isRecLineOK = rg_TRUE;
    
    string strAtomSerialsForBond[5];
    
    strAtomSerialsForBond[0] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_ATM_SERIAL_ST_POS,    PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[1] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_A_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[2] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_B_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[3] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_C_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[4] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_D_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    
    
    rg_INT atomSerialsForBond[5] = { -1, -1, -1, -1, -1 };
    
    rg_INT i_bondedAtom = 0;
    for ( i_bondedAtom=0; i_bondedAtom<5; i_bondedAtom++ ) {
        if( strAtomSerialsForBond[i_bondedAtom] != "" ) {
            atomSerialsForBond[i_bondedAtom] = atoi( strAtomSerialsForBond[i_bondedAtom].c_str() );
        }
    }
    
    
    Atom* atomsForBond[5] = { rg_NULL, rg_NULL, rg_NULL, rg_NULL, rg_NULL };
    
    AtomMap::iterator AtomMap_i = mapOfAtom.find( atomSerialsForBond[0] );
    
    if ( AtomMap_i != mapOfAtom.end() ) {
        atomsForBond[0] = (*AtomMap_i).second;
    }
    
    
    if ( atomsForBond[0] == rg_NULL ) {
        // CONNECTION BETWEEN ATOMS OF "ALTLOC" IS NOT DEFINED.
        // isRecLineOK = rg_FALSE;
    }
    else {
        for ( i_bondedAtom=1; i_bondedAtom<5; i_bondedAtom++ ) {
            
            if( atomSerialsForBond[i_bondedAtom] == -1 )
                continue;
            
            AtomMap_i = mapOfAtom.find( atomSerialsForBond[i_bondedAtom] );
            
            if ( AtomMap_i != mapOfAtom.end() ) {
                atomsForBond[i_bondedAtom] = (*AtomMap_i).second;
            }
            else {
                continue;
                // CONNECTION BETWEEN ATOMS OF "ALTLOC" IS NOT DEFINED.
                // isRecLineOK = rg_FALSE;
                // break;
            }
            
            setChemicalBondToMoleculeWithoutDuplication( atomsForBond[0], atomsForBond[i_bondedAtom], aMolecule );
        }            
    }
    
    return isRecLineOK;
}



rg_FLAG MoleculeIOFunctions::setConectRecordsToMolecule( const string& strRecLine, AtomMap& mapOfAtom, Molecule& aMolecule )
{
    rg_FLAG isRecLineOK = rg_TRUE;

    string strAtomSerialsForBond[5];
    
    strAtomSerialsForBond[0] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_ATM_SERIAL_ST_POS,    PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[1] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_A_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[2] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_B_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[3] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_C_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );
    strAtomSerialsForBond[4] = StringFunctions::subString( strRecLine, PDB_CONNECT_RECORD_BATM_D_SERIAL_ST_POS, PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH );


    rg_INT atomSerialsForBond[5] = { -1, -1, -1, -1, -1 };

    rg_INT i_bondedAtom = 0;
    for ( i_bondedAtom=0; i_bondedAtom<5; i_bondedAtom++ ) {
        if( strAtomSerialsForBond[i_bondedAtom] != "" ) {
            atomSerialsForBond[i_bondedAtom] = atoi( strAtomSerialsForBond[i_bondedAtom].c_str() );
        }
    }


    Atom* atomsForBond[5] = { rg_NULL, rg_NULL, rg_NULL, rg_NULL, rg_NULL };

    AtomMap::iterator AtomMap_i = mapOfAtom.find( atomSerialsForBond[0] );
    
    if ( AtomMap_i != mapOfAtom.end() ) {
        atomsForBond[0] = (*AtomMap_i).second;
    }
    
    
    if ( atomsForBond[0] == rg_NULL ) {
        // CONNECTION BETWEEN ATOMS OF "ALTLOC" IS NOT DEFINED.
        // isRecLineOK = rg_FALSE;
    }
    else {
        for ( i_bondedAtom=1; i_bondedAtom<5; i_bondedAtom++ ) {

            if( atomSerialsForBond[i_bondedAtom] == -1 )
                continue;

            AtomMap_i = mapOfAtom.find( atomSerialsForBond[i_bondedAtom] );
            
            if ( AtomMap_i != mapOfAtom.end() ) {
                atomsForBond[i_bondedAtom] = (*AtomMap_i).second;
            }
            else {
                continue;
                // CONNECTION BETWEEN ATOMS OF "ALTLOC" IS NOT DEFINED.
                // isRecLineOK = rg_FALSE;
                // break;
            }

            setChemicalBondToMoleculeWithoutDuplication( atomsForBond[0], atomsForBond[i_bondedAtom], aMolecule );
        }            
    }

    return isRecLineOK;
}



void MoleculeIOFunctions::setChemicalBondToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, Molecule* aMolecule )
{
    ChemicalBond newChemicalBond( aMolecule->getChemicalBonds()->getSize(), firstAtom, secondAtom );
    
    ChemicalBond* existingChemBondA = rg_NULL;
    ChemicalBond* existingChemBondB = rg_NULL;
    
    
    rg_FLAG isNewChemBondExistInChemBondListInAtomA = isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, firstAtom, existingChemBondA );
    rg_FLAG isNewChemBondExistInChemBondListInAtomB = isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, secondAtom, existingChemBondB );
    
    
    //// If "existingChemBondA" is not "rg_NULL" then "existingChemBondB" MUST! not be "rg_NULL".
    //// "existingChemBondA" == "existingChemBondB" --> ALWAYS!!!
    //// But equal comparison between "existingChemBondA" and "existingChemBondB" is performed for bug test.
    
    if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
        ChemicalBond* pNewChemicalBond = aMolecule->addChemicalBond( newChemicalBond );
        firstAtom->addChemicalBond( pNewChemicalBond );
        secondAtom->addChemicalBond( pNewChemicalBond );
    }
    else if( isNewChemBondExistInChemBondListInAtomA == rg_TRUE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
        secondAtom->addChemicalBond( existingChemBondA );
    }
    else if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_TRUE ) {
        firstAtom->addChemicalBond( existingChemBondB );
    }
}


// JKKIM modified for mol2 bond info.
//void MoleculeIOFunctions::setChemicalBondToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, Molecule& aMolecule )
ChemicalBond* MoleculeIOFunctions::setChemicalBondToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, Molecule& aMolecule )
{
    ChemicalBond newChemicalBond( aMolecule.getChemicalBonds()->getSize(), firstAtom, secondAtom );
    
    ChemicalBond* existingChemBondA = rg_NULL;
    ChemicalBond* existingChemBondB = rg_NULL;

    ChemicalBond* pNewChemicalBond = rg_NULL;

    rg_FLAG isNewChemBondExistInChemBondListInAtomA = isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, firstAtom, existingChemBondA );
    rg_FLAG isNewChemBondExistInChemBondListInAtomB = isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, secondAtom, existingChemBondB );
 

    //// If "existingChemBondA" is not "rg_NULL" then "existingChemBondB" MUST! not be "rg_NULL".
    //// "existingChemBondA" == "existingChemBondB" --> ALWAYS!!!
    //// But equal comparison between "existingChemBondA" and "existingChemBondB" is performed for bug test.

    if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
        pNewChemicalBond = aMolecule.addChemicalBond( newChemicalBond );
        firstAtom->addChemicalBond( pNewChemicalBond );
        secondAtom->addChemicalBond( pNewChemicalBond );
    }
    else if( isNewChemBondExistInChemBondListInAtomA == rg_TRUE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
        secondAtom->addChemicalBond( existingChemBondA );
    }
    else if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_TRUE ) {
        firstAtom->addChemicalBond( existingChemBondB );
    }

    return pNewChemicalBond;
}



void MoleculeIOFunctions::setChemicalBondForMolFileToMoleculeWithoutDuplication( Atom* firstAtom, Atom* secondAtom, const rg_INT& serialFromInputFile, const BondType& typeOfBond, Molecule& aMolecule )
{
    ChemicalBond newChemicalBond( aMolecule.getChemicalBonds()->getSize(), firstAtom, secondAtom );
    newChemicalBond.setSerialFromInputFile( serialFromInputFile );
    newChemicalBond.setTypeOfBond( typeOfBond );
        

    ChemicalBond* existingChemBondA = rg_NULL;
    ChemicalBond* existingChemBondB = rg_NULL;
    
    
    rg_FLAG isNewChemBondExistInChemBondListInAtomA = isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, firstAtom, existingChemBondA );
    rg_FLAG isNewChemBondExistInChemBondListInAtomB = isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, secondAtom, existingChemBondB );
    
    
    //// If "existingChemBondA" is not "rg_NULL" then "existingChemBondB" MUST! not be "rg_NULL".
    //// "existingChemBondA" == "existingChemBondB" --> ALWAYS!!!
    //// But equal comparison between "existingChemBondA" and "existingChemBondB" is performed for bug test.
    
    if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
        ChemicalBond* pNewChemicalBond = aMolecule.addChemicalBond( newChemicalBond );
        firstAtom->addChemicalBond( pNewChemicalBond );
        secondAtom->addChemicalBond( pNewChemicalBond );
    }
    else if( isNewChemBondExistInChemBondListInAtomA == rg_TRUE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
        secondAtom->addChemicalBond( existingChemBondA );
    }
    else if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_TRUE ) {
        firstAtom->addChemicalBond( existingChemBondB );
    }
}



rg_FLAG MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( ChemicalBond* aChemicalBond, Atom* atomForBond, ChemicalBond* existingChemicalBond )
{
    existingChemicalBond           = rg_NULL;
    rg_FLAG isNewChemicalBondExist = rg_FALSE;
    
    rg_dList<ChemicalBond*>* bondListInAtom = atomForBond->getListChemicalBond();

    bondListInAtom->reset4Loop();
    while( bondListInAtom->setNext4Loop() ) {
        ChemicalBond* currBondInBondList = bondListInAtom->getEntity();
    
        if( *currBondInBondList == *aChemicalBond ) {
            isNewChemicalBondExist = rg_TRUE;
            existingChemicalBond   = currBondInBondList;
            break;
        }
    }
    
    return isNewChemicalBondExist;
}



void MoleculeIOFunctions::setChemicalBondsToMoleculeForStandardResidues( Molecule* aMolecule )
{
    rg_dList<Chain>* chainsInMolecule = aMolecule->getChains();
    
    chainsInMolecule->reset4Loop();
    
    while( chainsInMolecule->setNext4Loop() ) {
        
        Chain* currChain = chainsInMolecule->getpEntity();
        
        // Bonds between residues in Chain
        setChemicalBondsBetweenStdResiduesInChain( currChain, aMolecule );
        
        
        // Bonds between atoms in residue
        rg_dList<Residue*>* residuesInChain = currChain->getResidues();
        
        residuesInChain->reset4Loop();
        while( residuesInChain->setNext4Loop() ) {
            Residue* currResidue = residuesInChain->getEntity();
            
            setChemicalBondsBetweenAtomsInStdResidue( currResidue, aMolecule );
        }
    }
}



void MoleculeIOFunctions::setChemicalBondsToMoleculeForStandardResidues( Molecule& aMolecule )
{
    rg_dList<Chain>* chainsInMolecule = aMolecule.getChains();

    chainsInMolecule->reset4Loop();

    while( chainsInMolecule->setNext4Loop() ) {
        
        Chain* currChain = chainsInMolecule->getpEntity();
        
        // Bonds between residues in Chain
        setChemicalBondsBetweenStdResiduesInChain( currChain, aMolecule );

        
        // Bonds between atoms in residue
        rg_dList<Residue*>* residuesInChain = currChain->getResidues();
        
        residuesInChain->reset4Loop();
        while( residuesInChain->setNext4Loop() ) {
            Residue* currResidue = residuesInChain->getEntity();

            setChemicalBondsBetweenAtomsInStdResidue( currResidue, aMolecule );
        }
    }
}



void MoleculeIOFunctions::setChemicalBondsBetweenStdResiduesInChain( Chain* aChain, Molecule* aMolecule )
{
    rg_INT minResidueSeq = INT_MAX;
    rg_INT maxResidueSeq = 0;
    
    rg_dList<Residue*>* residuesInChain = aChain->getResidues();
    
    Residue* currResidue = rg_NULL;
    residuesInChain->reset4Loop();
    while( residuesInChain->setNext4Loop() ) {
        currResidue = residuesInChain->getEntity();
        
        rg_INT residueSequenceNumber = currResidue->getSequenceNumber();
        
        if ( residueSequenceNumber < minResidueSeq )
            minResidueSeq = residueSequenceNumber;
        
        if ( residueSequenceNumber > maxResidueSeq )
            maxResidueSeq = residueSequenceNumber;        
    }
    
    rg_INT numOfResSeq = maxResidueSeq - minResidueSeq + 1;
    Residue** residueArray = new Residue*[ numOfResSeq ];
    rg_INT i_residueSeq = 0;
    
    for ( i_residueSeq=0; i_residueSeq<numOfResSeq; i_residueSeq++ )
        residueArray[i_residueSeq] = rg_NULL;
    
    
    residuesInChain->reset4Loop();
    while( residuesInChain->setNext4Loop() ) {
        currResidue = residuesInChain->getEntity();
        
        residueArray[ currResidue->getSequenceNumber()-minResidueSeq ] = currResidue;
    }    
    
    for ( i_residueSeq=0; i_residueSeq<numOfResSeq-1; i_residueSeq++ )  {
        if ( residueArray[i_residueSeq] == rg_NULL || residueArray[i_residueSeq+1] == rg_NULL )
            continue;
        
        if ( residueArray[i_residueSeq]->isAminoResidue()   == rg_TRUE && 
            residueArray[i_residueSeq+1]->isAminoResidue() == rg_TRUE  ) {
            
            setChemicalBondBetweenAminoAcids( residueArray[i_residueSeq], residueArray[i_residueSeq+1], aMolecule );
            continue;
        }
        
        if ( residueArray[i_residueSeq]->isDNAResidue()   == rg_TRUE && 
            residueArray[i_residueSeq+1]->isDNAResidue() == rg_TRUE  ) {
            
            setChemicalBondBetweenNucleicAcids( residueArray[i_residueSeq], residueArray[i_residueSeq+1], aMolecule );
            continue;
        }
        
        if ( residueArray[i_residueSeq]->isRNAResidue()   == rg_TRUE && 
            residueArray[i_residueSeq+1]->isRNAResidue() == rg_TRUE  ) {
            
            setChemicalBondBetweenNucleicAcids( residueArray[i_residueSeq], residueArray[i_residueSeq+1], aMolecule );
        }
    }
    
    delete [] residueArray ;
}



void MoleculeIOFunctions::setChemicalBondsBetweenStdResiduesInChain( Chain* aChain, Molecule& aMolecule )
{
    rg_INT minResidueSeq = INT_MAX;
    rg_INT maxResidueSeq = 0;
    
    rg_dList<Residue*>* residuesInChain = aChain->getResidues();
    
    Residue* currResidue = rg_NULL;
    residuesInChain->reset4Loop();
    while( residuesInChain->setNext4Loop() ) {
        currResidue = residuesInChain->getEntity();

        rg_INT residueSequenceNumber = currResidue->getSequenceNumber();

        if ( residueSequenceNumber < minResidueSeq )
            minResidueSeq = residueSequenceNumber;
        
        if ( residueSequenceNumber > maxResidueSeq )
            maxResidueSeq = residueSequenceNumber;        
    }

    rg_INT numOfResSeq = maxResidueSeq - minResidueSeq + 1;
    Residue** residueArray = new Residue*[ numOfResSeq ];
    rg_INT i_residueSeq = 0;

    for ( i_residueSeq=0; i_residueSeq<numOfResSeq; i_residueSeq++ )
        residueArray[i_residueSeq] = rg_NULL;


    residuesInChain->reset4Loop();
    while( residuesInChain->setNext4Loop() ) {
        currResidue = residuesInChain->getEntity();
    
        residueArray[ currResidue->getSequenceNumber()-minResidueSeq ] = currResidue;
    }    

    for ( i_residueSeq=0; i_residueSeq<numOfResSeq-1; i_residueSeq++ )  {
        if ( residueArray[i_residueSeq] == rg_NULL || residueArray[i_residueSeq+1] == rg_NULL )
            continue;

        if ( residueArray[i_residueSeq]->isAminoResidue()   == rg_TRUE && 
             residueArray[i_residueSeq+1]->isAminoResidue() == rg_TRUE  ) {

            setChemicalBondBetweenAminoAcids( residueArray[i_residueSeq], residueArray[i_residueSeq+1], aMolecule );
            continue;
        }

        if ( residueArray[i_residueSeq]->isDNAResidue()   == rg_TRUE && 
             residueArray[i_residueSeq+1]->isDNAResidue() == rg_TRUE  ) {
            
            setChemicalBondBetweenNucleicAcids( residueArray[i_residueSeq], residueArray[i_residueSeq+1], aMolecule );
            continue;
        }

        if ( residueArray[i_residueSeq]->isRNAResidue()   == rg_TRUE && 
            residueArray[i_residueSeq+1]->isRNAResidue() == rg_TRUE  ) {
            
            setChemicalBondBetweenNucleicAcids( residueArray[i_residueSeq], residueArray[i_residueSeq+1], aMolecule );
        }
    }

    delete [] residueArray ;
}



void MoleculeIOFunctions::setChemicalBondBetweenAminoAcids( Residue* firstAminoAcid, Residue* secondAminoAcid, Molecule* aMolecule )
{
    Atom* firstAtom  = rg_NULL;
    Atom* secondAtom = rg_NULL;
    
    
    rg_dList<Atom*>* atomList = firstAminoAcid->getAtoms();
    
    Atom* tempAtom = rg_NULL;
    
    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " C  ") {
            firstAtom = tempAtom;
            break;
        }
    }
    
    atomList = secondAminoAcid->getAtoms();
    
    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " N  ") {
            secondAtom = tempAtom;
            break;
        }
    }
    
    if ( firstAtom != rg_NULL && secondAtom != rg_NULL )
        setChemicalBondToMoleculeWithoutDuplication( firstAtom, secondAtom, aMolecule );
}



void MoleculeIOFunctions::setChemicalBondBetweenAminoAcids( Residue* firstAminoAcid, Residue* secondAminoAcid, Molecule& aMolecule )
{
    Atom* firstAtom  = rg_NULL;
    Atom* secondAtom = rg_NULL;
    

    rg_dList<Atom*>* atomList = firstAminoAcid->getAtoms();

    Atom* tempAtom = rg_NULL;

    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();

        if( tempAtom->getAtomNameFromInputFile() == " C  ") {
            firstAtom = tempAtom;
            break;
        }
    }

    atomList = secondAminoAcid->getAtoms();

    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " N  ") {
            secondAtom = tempAtom;
            break;
        }
    }
        
    if ( firstAtom != rg_NULL && secondAtom != rg_NULL )
        setChemicalBondToMoleculeWithoutDuplication( firstAtom, secondAtom, aMolecule );
}



void MoleculeIOFunctions::setChemicalBondBetweenNucleicAcids( Residue* firstNucleicAcid, Residue* secondNucleicAcid, Molecule* aMolecule )
{
    Atom* firstAtom  = rg_NULL;
    Atom* secondAtom = rg_NULL;
    
    
    rg_dList<Atom*>* atomList = firstNucleicAcid->getAtoms();
    
    Atom* tempAtom = rg_NULL;
    
    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " O3*") {
            firstAtom = tempAtom;
            break;
        }
    }
    
    atomList = secondNucleicAcid->getAtoms();
    
    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " P  ") {
            secondAtom = tempAtom;
            break;
        }
    }
    
    if ( firstAtom != rg_NULL && secondAtom != rg_NULL )
        setChemicalBondToMoleculeWithoutDuplication( firstAtom, secondAtom, aMolecule );
}



void MoleculeIOFunctions::setChemicalBondBetweenNucleicAcids( Residue* firstNucleicAcid, Residue* secondNucleicAcid, Molecule& aMolecule )
{
    Atom* firstAtom  = rg_NULL;
    Atom* secondAtom = rg_NULL;
    
    
    rg_dList<Atom*>* atomList = firstNucleicAcid->getAtoms();
    
    Atom* tempAtom = rg_NULL;
    
    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " O3*") {
            firstAtom = tempAtom;
            break;
        }
    }
    
    atomList = secondNucleicAcid->getAtoms();
    
    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " P  ") {
            secondAtom = tempAtom;
            break;
        }
    }
    
    if ( firstAtom != rg_NULL && secondAtom != rg_NULL )
        setChemicalBondToMoleculeWithoutDuplication( firstAtom, secondAtom, aMolecule );
}



void MoleculeIOFunctions::setChemicalBondsBetweenAtomsInStdResidue( Residue* aResidue, Molecule* aMolecule )
{
    typedef map<string, Atom*> AtomNameMap;

    if( aResidue->isStandardResidue() == rg_TRUE ) {

        //// Initialize map of atom names ( key: (string)nameOfAtom, value: Atom* )
        //
        AtomNameMap mapOfAtomName;
        mapOfAtomName.clear();
        
        rg_dList<Atom*>* atomsInResidue = aResidue->getAtoms();
        
        atomsInResidue->reset4Loop();
        while( atomsInResidue->setNext4Loop() ) {
            Atom* currAtom = atomsInResidue->getEntity();
            mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
        }
        //
        ////


        rg_INT i_residueCode = (rg_INT)aResidue->getResidueCode();

        for( rg_INT i_bond=0; i_bond<RESIDUE_FEATURES[i_residueCode].numOfBonds; i_bond++ )  {
            
            string atomName[2];
            Atom*  atomPtr[2];

            atomName[0] = RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][0];
            atomName[1] = RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][1];
            
            atomPtr[0]  = rg_NULL;
            atomPtr[1]  = rg_NULL;

            AtomNameMap::iterator atomNameMap_i;
            
            for( rg_INT i_atomPair=0; i_atomPair<2; i_atomPair++ ) {

                atomNameMap_i = mapOfAtomName.find( atomName[i_atomPair] );

                if ( atomNameMap_i != mapOfAtomName.end() ) {
                    atomPtr[i_atomPair] = (*atomNameMap_i).second;
                }
            }

            if( atomPtr[0] == rg_NULL || atomPtr[1] == rg_NULL )
                continue;

            setChemicalBondToMoleculeWithoutDuplication( atomPtr[0], atomPtr[1], aMolecule );
        }
    }
}



void MoleculeIOFunctions::setChemicalBondsBetweenAtomsInStdResidue( Residue* aResidue, Molecule& aMolecule )
{
    typedef map<string, Atom*> AtomNameMap;
    
    if( aResidue->isStandardResidue() == rg_TRUE ) {
        
        //// Initialize map of atom names ( key: (string)nameOfAtom, value: Atom* )
        //
        AtomNameMap mapOfAtomName;
        mapOfAtomName.clear();
        
        rg_dList<Atom*>* atomsInResidue = aResidue->getAtoms();
        
        atomsInResidue->reset4Loop();
        while( atomsInResidue->setNext4Loop() ) {
            Atom* currAtom = atomsInResidue->getEntity();
            mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
        }
        //
        ////
        
        
        rg_INT i_residueCode = (rg_INT)aResidue->getResidueCode();
        
        for( rg_INT i_bond=0; i_bond<RESIDUE_FEATURES[i_residueCode].numOfBonds; i_bond++ )  {
            
            string atomName[2];
            Atom*  atomPtr[2];
            
            atomName[0] = RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][0];
            atomName[1] = RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][1];
            
            atomPtr[0]  = rg_NULL;
            atomPtr[1]  = rg_NULL;
            
            AtomNameMap::iterator atomNameMap_i;
            
            for( rg_INT i_atomPair=0; i_atomPair<2; i_atomPair++ ) {
                
                atomNameMap_i = mapOfAtomName.find( atomName[i_atomPair] );
                
                if ( atomNameMap_i != mapOfAtomName.end() ) {
                    atomPtr[i_atomPair] = (*atomNameMap_i).second;
                }
            }
            
            if( atomPtr[0] == rg_NULL || atomPtr[1] == rg_NULL )
                continue;
            
            setChemicalBondToMoleculeWithoutDuplication( atomPtr[0], atomPtr[1], aMolecule );
        }
    }
}



rg_FLAG MoleculeIOFunctions::readTriposMol2File( const char* pathFile, rg_dList<Molecule>& molecules )
{
    rg_BOOL isModelRead = rg_TRUE;
    Molecule* aMolecule;
    rg_INT modelID = 1;

    while(isModelRead) {
       aMolecule = molecules.pushBack(Molecule());

       isModelRead = readTriposMol2File(modelID++, pathFile, *aMolecule);        
    }    

    molecules.killTail();

    return rg_TRUE;
}


rg_FLAG MoleculeIOFunctions::readTriposMol2File( const rg_INT& targetModelID, const char* pathFile, Molecule& aMolecule )
{
    rg_INT modelID = 0;
    rg_FLAG isFileRead   = rg_FALSE;

    if( !isValidMoleculeFileExtension( pathFile, MOL2_FILE ) ) {
        cerr << "Error: Invalid file extension !\n";
        return isFileRead;
    }
    
    
    ifstream fin( pathFile );
    
    if( fin.bad() ) {
        cerr << "Error: Could not open file!\n";
        return isFileRead;
    }

    
//     AtomMap     * mapOfAtom     = new AtomMap;
//     ResidueMap  * mapOfResidue  = new ResidueMap;
//     ChainMap    * mapOfChain    = new ChainMap;
//     ChemBondMap * mapOfChemBond = new ChemBondMap;
//     
//     AtomSymbolMap    *mapOfAtomSymbol    = new AtomSymbolMap;
//     ResidueSymbolMap *mapOfResidueSymbol = new ResidueSymbolMap;
// 
// // CAN'T READ MOL2 FILE... NEED TO DEBUG....
//     ResidueMap  * mapOfNonStdResidue = new ResidueMap;

    //JKKIM Modified 2011-07-11 pointer->object
    AtomMap     mapOfAtom;
    ResidueMap  mapOfResidue;
    ChainMap    mapOfChain;
    ChemBondMap mapOfChemBond;
    
    AtomSymbolMap    mapOfAtomSymbol;
    ResidueSymbolMap mapOfResidueSymbol;

    ResidueMap  mapOfNonStdResidue;

    
    initiallizeSymbolMapsOfAtomAndResidueForMol2File( mapOfAtomSymbol, mapOfResidueSymbol );

    rg_INT arrNumOfRecords[5] = { 0, 0, 0, 0, 0 };

    string currRecType = "";
    string strRecLine  = "";

    if( !getline( fin, strRecLine ) )
        return isFileRead;

    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE

    list<string> recordsOfMol2File;
    addRecordsOfMol2FIleToList( &fin, recordsOfMol2File );
    fin.close();
    

    rg_INT countLines = 0;
    rg_INT numOfRecLines = recordsOfMol2File.size();
    rg_FLAG isRecLinesOK = rg_TRUE;

    list<string>::iterator i_recLines = recordsOfMol2File.begin();

    while( i_recLines != recordsOfMol2File.end() ) {
        
        countLines++;
        string* currRecLine = &(*i_recLines);

        string recType = StringFunctions::strTrim( *currRecLine );


        if ( recType == MOL2_RECORD_TYPE_MOLECULE ) {
            modelID++;

            if(modelID > targetModelID) {
                break;
            }

            if(modelID ==targetModelID) {
                isRecLinesOK = setMoleculeRTIToMolecule( i_recLines, arrNumOfRecords, aMolecule );
            }
            
        }        
        else if ( recType == MOL2_RECORD_TYPE_ATOM && arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID] > 0 ) {
            if(modelID == targetModelID) {
                isRecLinesOK = setAtomRTIToMolecule( i_recLines, &recordsOfMol2File, arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID], mapOfAtom, mapOfResidue, mapOfNonStdResidue, mapOfAtomSymbol, mapOfResidueSymbol, aMolecule );
            }
        }
        else if ( recType == MOL2_RECORD_TYPE_BOND && arrNumOfRecords[MOL2_NUM_OF_BONDS_ID] >0 ) {
            if(modelID == targetModelID) {
                isRecLinesOK = setBondRTIToMolecule( i_recLines, &recordsOfMol2File, arrNumOfRecords[MOL2_NUM_OF_BONDS_ID], mapOfAtom, aMolecule );
            }
        }
        else if ( recType == MOL2_RECORD_TYPE_SUBSTRUCTURE ) {
            int aaa = 0;
            
        }

        if( isRecLinesOK == rg_FALSE ) {
            isFileRead = rg_FALSE;
            break;
        }

        i_recLines++;
    }

    setAmberAtomTypeForHAtom(aMolecule);

//    int numOfNonStdResidue = mapOfNonStdResidue.size();
    
    //setChainsForNonStdResiduesForMol2File( &mapOfNonStdResidue, aMolecule );
    evaluateAndSetChainCode( aMolecule );

    // SET FILE NAME TO MOLECULE
    string fileName = StringFunctions::getFileNameWithoutPath( string(pathFile) ) ;
    aMolecule.setMoleculeFileName( fileName );
     

    if( isRecLinesOK ) {
        aMolecule.evaluateHydrogenDonorAndAcceptorAtoms();
        aMolecule.computeAndSetCenterOfMass();
        aMolecule.computeAndSetMinEnclosingSphere();
        isFileRead = rg_TRUE;
    }

    if(aMolecule.getAtoms()->getSize() == 0) {
        isFileRead = rg_FALSE;
    }

//     delete mapOfAtom;
//     delete mapOfResidue;
//     delete mapOfChain;
//     delete mapOfChemBond;
//     delete mapOfNonStdResidue;
// 
//     delete mapOfAtomSymbol;   
//     delete mapOfResidueSymbol;

    return isFileRead;
}



void MoleculeIOFunctions::addRecordsOfMol2FIleToList( ifstream* fin, list<string>& recordsOfMol2File )
{
    string strRecLine  = "";
    
    while ( getline( *fin, strRecLine ) ) {
        
        string trimedRecLine = StringFunctions::strTrim( strRecLine );
       
        if( trimedRecLine == "" )
            continue;

        recordsOfMol2File.push_back( trimedRecLine );
    }
}



rg_BOOL MoleculeIOFunctions::isNewTriposRTILineStart( const string& strRecLine )
{         
    if( strRecLine.substr( MOL2_RECORD_RTI_ST_POS, MOL2_RECORD_RTI_LENGTH ) == "@<TRIPOS>" )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG MoleculeIOFunctions::setMoleculeRTIToMolecule( list<string>::iterator& i_recLines, rg_INT* arrNumOfRecords, Molecule& aMolecule )
{
    rg_FLAG isMoleculeRTIOK = rg_TRUE;

    i_recLines++;
    string strRecLine  = *i_recLines;

    string  nameOfMolecule       = "";
    rg_BOOL isNameOfMoleculeRead = rg_FALSE;
    rg_BOOL isNumOfRecordsRead   = rg_FALSE;


    while( !isNewTriposRTILineStart( strRecLine ) ) {

        string firstLetter = strRecLine.substr(0, 1);
        if( firstLetter == " " || firstLetter == "" ) {
            
            i_recLines++;
            strRecLine = *i_recLines;
            continue;
        }


        // MOLECULE NAME
        if( isNameOfMoleculeRead == rg_FALSE /* && StringFunctions::isStringLetterInAlphabet( firstLetter ) */ ) {
            nameOfMolecule = StringFunctions::strTrim( strRecLine );
            
            aMolecule.setMoleculeName(nameOfMolecule);

            i_recLines++;
            strRecLine = StringFunctions::strTrim( *i_recLines );
            
            firstLetter = strRecLine.substr(0, 1);
            isNameOfMoleculeRead = rg_TRUE;
        }

        // NUM OF ATOMS/BONDS/SUBST/FEAT/SETS
        if( isNumOfRecordsRead == rg_FALSE && StringFunctions::isStringLetterInNumber( firstLetter ) ) {
            rg_INT i_numOfRec = 0;
            while( strRecLine.length() != 0 ) {

                rg_INT posOfBlank = strRecLine.find( " ", 0 );

                if(posOfBlank == -1) {
                    posOfBlank = strRecLine.find( "\t", 0 );
                }
                
                string numOfRecs = strRecLine.substr( 0, posOfBlank );
                
                arrNumOfRecords[i_numOfRec] = atoi( numOfRecs.c_str() );
                i_numOfRec++;

                strRecLine.erase( 0, numOfRecs.length() );
                strRecLine = StringFunctions::strTrimLeft( strRecLine );
            }

            isNumOfRecordsRead = rg_TRUE;
        }

        // MOLECULE TYPE : TO BE CONSIDERED...
        // CHARGE TYPE   : TO BE CONSIDERED...
            
        i_recLines++;
        strRecLine = *i_recLines;
        
    }

    i_recLines--;

    if ( isNameOfMoleculeRead == rg_FALSE && isNumOfRecordsRead == rg_FALSE )
        isMoleculeRTIOK = rg_FALSE;

    return isMoleculeRTIOK;
}



// CAN'T READ MOL2 FILE... NEED TO DEBUG....
rg_FLAG MoleculeIOFunctions::setAtomRTIToMolecule( 
    list<string>::iterator& i_recLines, 
    list<string>* recordsOfMol2File, 
    const rg_INT& numOfAtomsGiven, 
    AtomMap& mapOfAtom, 
    ResidueMap& mapOfResidue, 
    ResidueMap& mapOfNonStdResidue, 
    AtomSymbolMap& mapOfAtomSymbol, 
    ResidueSymbolMap& mapOfResidueSymbol, 
    Molecule& aMolecule )
{
    rg_FLAG isAtomRTIOK = rg_TRUE;

    i_recLines++;
    string strRecLine  = StringFunctions::strTrim( *i_recLines );
    
    rg_INT numOfChains = aMolecule.getChains()->getSize();
    Chain* newChain    = aMolecule.addChain( Chain(numOfChains, &aMolecule) );
    rg_BOOL isNewChainStart = rg_FALSE;


    rg_INT i_atom = 0;
    
    while( !isNewTriposRTILineStart( strRecLine ) ) {

        i_atom++;

        rg_INT  arrIntRecs[MOL2_ATOM_INT_REC_SIZE]  = { -1, -1 };
        rg_REAL arrRealRecs[MOL2_ATOM_REAL_REC_SIZE] = { 0.0, 0.0, 0.0 };
        string  arrStrRecs[MOL2_ATOM_STR_REC_SIZE]  = { "", "", "", "" };
        
        extractAtomRecordsFromAtomRTILine( strRecLine, arrIntRecs, arrRealRecs, arrStrRecs );

        if( arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID] == 0 && arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] == 0 ) {
            break;
        }
        

        // Estimate Chain : if name of prev atom was "OXT" then new Chain must be started.
        if( isNewChainStart == rg_TRUE ) {
            numOfChains = aMolecule.getChains()->getSize();
            newChain    = aMolecule.addChain( Chain(numOfChains, &aMolecule) );
            isNewChainStart = rg_FALSE;
        }

        
        // Estimate Residue : if residue is not in map, then create a new residue.
        Residue* currResidue = getResiduePtrFromMap( &mapOfResidue,  arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] );

        if ( currResidue == rg_NULL ) {
            rg_INT newResidueID = mapOfResidue.size();
            currResidue = aMolecule.addResidue( Residue(newResidueID) );
            currResidue->setSequenceNumber( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] );
            currResidue->setResidueName( arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME] );
            
            
            
            int aaa = MOL2_ATOM_INT_REC_SUBST_ID;

            string residueName = getResidueNameFromTriposSubstName( arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME] );
            setResidueCodeToTargetResidue( mapOfResidueSymbol, residueName, currResidue );    // NOT WORKING WITH DNA/RNA.
            mapOfResidue.insert( ResidueMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID], currResidue ) );


            ////////////////////////////////////////////////////////////////////////
            //
            //  2018.12.24 by Y.Cho
            newChain->addResidue( currResidue );
            currResidue->setChain( newChain );

            if (!currResidue->isStandardResidue()) {
                mapOfNonStdResidue.insert( ResidueMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID], currResidue ) );
            }

            //if( currResidue->isStandardResidue() ) {
            //    newChain->addResidue( currResidue );
            //    currResidue->setChain( newChain );
            //}
            //else {
            //    mapOfNonStdResidue.insert( ResidueMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID], currResidue ) );
            //}

            //
            ////////////////////////////////////////////////////////////////////////
        }

        // Estimate Atom
        rg_INT newAtomID = mapOfAtom.size();
        
        Atom* newAtom = aMolecule.addAtom( Atom( newAtomID ) );


        currResidue->addAtom( newAtom );
        newAtom->setResidue( currResidue );
        newAtom->setSerialFromInputFile( arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID] );
        newAtom->setAtomNameFromInputFile( arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME] );
        rg_FLAG isAtomNameOK = setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForMol2File( mapOfAtomSymbol, arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME], arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE], newAtom );
        newAtom->setAtomNameFromInputFile( newAtom->getAtomNameInPDBFormat());


        //if( arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME] == "OXT") {
        //    isNewChainStart = rg_TRUE;
        //}

        if ( !isAtomNameOK ) {
            isAtomRTIOK = rg_FALSE;
            break;
        }

        newAtom->setAtomBall( Sphere( arrRealRecs[MOL2_ATOM_REAL_REC_X_COORD], arrRealRecs[MOL2_ATOM_REAL_REC_Y_COORD], arrRealRecs[MOL2_ATOM_REAL_REC_Z_COORD], ATOM_FEATURES[newAtom->getAtomCode()].radius) );
        newAtom->getpChemicalProperties()->setCharge( arrRealRecs[MOL2_ATOM_REAL_REC_CHARGE] );      
        
        newAtom->getpChemicalProperties()->setSYBYLAtomType(arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE]);

        AmberAtomTypes amberAtomType = getAmberAtomTypeFromSYBYL(arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE]);

        newAtom->getpChemicalProperties()->setAtomTypeInAmber(amberAtomType);
        
        mapOfAtom.insert( AtomMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID], newAtom ) );


        i_recLines++;

        if ( i_recLines != recordsOfMol2File->end() ) {
            strRecLine = StringFunctions::strTrim( *i_recLines );
        }
        else {
            break;
        }

        if( i_atom == numOfAtomsGiven )
            break;
    }

    // Filter chains with no residue and reset ID.
    filterEmptyChain( aMolecule );

    if( i_atom != numOfAtomsGiven )
        isAtomRTIOK = rg_FALSE;

    i_recLines--;

    return isAtomRTIOK;
}



void MoleculeIOFunctions::extractAtomRecordsFromAtomRTILine( const string&  strRecLine, rg_INT* arrIntRecs, rg_REAL* arrRealRecs, string* arrStrRecs )
{

    istringstream recordFin(strRecLine);

    recordFin >> arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID];
    recordFin >> arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME];
    recordFin >> arrRealRecs[MOL2_ATOM_REAL_REC_X_COORD];
    recordFin >> arrRealRecs[MOL2_ATOM_REAL_REC_Y_COORD];
    recordFin >> arrRealRecs[MOL2_ATOM_REAL_REC_Z_COORD];

    recordFin >> arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE];
    recordFin >> arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID];
    recordFin >> arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME];
    recordFin >> arrRealRecs[MOL2_ATOM_REAL_REC_CHARGE];
    recordFin >> arrStrRecs[MOL2_ATOM_STR_REC_STATUS_BIT];

    int stop = 1;

/*
    rg_INT stPosOfSingleRec = 0;
    rg_INT edPosOfSingleRec = 0;

    // ATOM ID
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID] = atoi( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
//    cout << "ATOM ID : " << arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID] << endl;   
    

    // ATOM NAME
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME] = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
//    cout << "ATOM NAME : " << arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME] << endl;
    

    // X-COORD OF ATOM CENTER
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrRealRecs[MOL2_ATOM_REAL_REC_X_COORD] = atof( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
//    cout << "X-COORD OF ATOM CENTER : " << arrRealRecs[MOL2_ATOM_REAL_REC_X_COORD] << endl;
    
    
    // Y-COORD OF ATOM CENTER
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrRealRecs[MOL2_ATOM_REAL_REC_Y_COORD] = atof( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
//    cout << "Y-COORD OF ATOM CENTER : " << arrRealRecs[MOL2_ATOM_REAL_REC_Y_COORD] << endl;
    
    
    // Z-COORD OF ATOM CENTER
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrRealRecs[MOL2_ATOM_REAL_REC_Z_COORD] = atof( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
//    cout << "Z-COORD OF ATOM CENTER : " << arrRealRecs[MOL2_ATOM_REAL_REC_Z_COORD] << endl;
    
    
    // ATOM TYPE
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE] = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
//    cout << "ATOM TYPE : " << arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE] << endl;
    
    
    // SUBSTRUCTURE ID (RESIDUE ID)
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] = atoi( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
//    cout << "SUBSTRUCTURE ID (RESIDUE ID) : " << arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] << endl;
        
    
    // SUBSTRUCTURE NAME (RESIDUE NAME)
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME] = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
//    cout << "SUBSTRUCTURE NAME (RESIDUE NAME) : " << arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME] << endl;
    
    
    // CHARGE
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
       
    arrRealRecs[MOL2_ATOM_REAL_REC_CHARGE] = atof( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
//    cout << "CHARGE : " << arrRealRecs[MOL2_ATOM_REAL_REC_CHARGE] << endl;
    
    
//    // STATUS BIT
//    stPosOfSingleRec = edPosOfSingleRec;
//    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
//    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
//    
//    arrStrRecs[MOL2_ATOM_STR_REC_STATUS_BIT] = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
////    cout << "STATUS BIT : " << arrStrRecs[MOL2_ATOM_STR_REC_STATUS_BIT] << endl;
    */
}



Residue* MoleculeIOFunctions::getResiduePtrFromMap( ResidueMap* mapOfResidue,  const rg_INT& residueID )
{
    ResidueMap::iterator residueMap_i = mapOfResidue->find( residueID );

    if ( residueMap_i == mapOfResidue->end() )
        return rg_NULL;
    else
        return (*residueMap_i).second;
}



string MoleculeIOFunctions::getResidueNameFromTriposSubstName( const string& substName )
{
    rg_INT posOfEndOfResName = StringFunctions::getPositionStartWithNumber( substName );

    return substName.substr( 0, posOfEndOfResName );
}



rg_FLAG  MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForMol2File( AtomSymbolMap& mapOfAtomSymbol, const string& atomName, const string& atomType, Atom* targetAtom )
{
    rg_FLAG isAtomNameOK = rg_TRUE;
    
    ////   SOME PDB FILES(ex: 1g9v) DO NOT FOLLOW THE CONDITIONS ON BELOW....
    //
    //     // Check if remoteness indicator is in numeric letter : remoteIndicator must be in alphabet letter !
    //     if( StringFunctions::isStringLetterInNumber( atomName.substr( 2, 1 ) ) == true ) {
    //         isAtomNameOK = rg_FALSE;
    //     }
    //     
    //     // Check if brangeDesignator is in Alphabet : remoteIndicator must be in Numeric letter !
    //     if( StringFunctions::isStringLetterInAlphabet( atomName.substr( 3, 1 ) ) == true ) {
    //         isAtomNameOK = rg_FALSE;
    //     }

    rg_INT posOfAtomTypeEnd = atomType.find(".");
    string typeOfAtom;
    if( posOfAtomTypeEnd == string::npos ) {
        typeOfAtom = StringFunctions::strToUpper( atomType );
    }
    else {
        typeOfAtom = StringFunctions::strToUpper( atomType.substr(0, posOfAtomTypeEnd) );
    }

    // Dummy but carbon or oxygen
    if ( typeOfAtom == "DU" ) {
        typeOfAtom = atomName.substr(0, 1);
    }
    // Lone pair... : ex) Mercury
    if ( typeOfAtom == "LP") {
        typeOfAtom = atomName;
    }


    //if( targetAtom->getResidue()->isStandardResidue() ) {

    //    string typeOfAtomFromAtomName = StringFunctions::strToUpper( atomName.substr( 0, typeOfAtom.length() ) );
    //    
    //    if( typeOfAtomFromAtomName != typeOfAtom )
    //        isAtomNameOK = rg_FALSE;
    //}
    //else {
    //}



    if( isAtomNameOK == rg_TRUE ) {
        
        // SET ATOM CODE
        AtomSymbolMap::iterator AtomSymbolMap_i = mapOfAtomSymbol.find( typeOfAtom );
        
        AtomCode convertedAtomCode = UNK_ATOM;
        if ( AtomSymbolMap_i != mapOfAtomSymbol.end() ) {
            convertedAtomCode = (AtomCode)((*AtomSymbolMap_i).second);
        }
        
        targetAtom->setAtomCode( convertedAtomCode );
        
         // set ChemicalProperties
         isAtomNameOK = setChemicalPropertiesToNewAtomFromTriposAtomName( atomName, typeOfAtom, targetAtom );
    }
    
    return isAtomNameOK;
}



rg_FLAG MoleculeIOFunctions::setChemicalPropertiesToNewAtomFromTriposAtomName( const string& atomName, const string& atomType, Atom* targetAtom )
{
    rg_FLAG isRecordOK = rg_TRUE;
    
    rg_INT lengthOfAtomName = atomName.length();

    if ( lengthOfAtomName > 0 && lengthOfAtomName < 5 ) {

        RemoteIndicator  remoteIndicator       = UNK_REMOTE;
        BranchDesignator branchDesignator      = UNK_BRANCH;
        BranchDesignator extraBranchDesignator = UNK_BRANCH;

        if ( targetAtom->getResidue()->isStandardResidue() ) {

            //string strRemoteAndBranch = StringFunctions::subString(atomName, atomType.length(), atomName.length());
            string::size_type pos = atomName.find(atomType);
            string::size_type len = atomType.length();
            string strRemoteAndBranch = StringFunctions::subString(atomName, pos+len, atomName.length());


            rg_INT lengthOfStrRemoteAndBranch = strRemoteAndBranch.length();
            
            switch( lengthOfStrRemoteAndBranch ) {
            case 0 : {
                remoteIndicator       = UNK_REMOTE;
                branchDesignator      = UNK_BRANCH;
                if (pos == 1) {
                    extraBranchDesignator = convertStringIntoBranchDesignator(atomName.substr(0, 1));
                }
                else {
                    extraBranchDesignator = UNK_BRANCH;
                }
                break;
            }
            case 1 : { 
                string strRemoteIndicator = strRemoteAndBranch.substr(0, 1);
                
                if( isStringValidForRemoteIndicator( strRemoteIndicator ) ) {           // REMOTENESS IN CHAR
                    remoteIndicator = convertStringIntoRemoteIndicator( strRemoteIndicator );
                    if (pos == 1) {
                        extraBranchDesignator = convertStringIntoBranchDesignator(atomName.substr(0, 1));
                    }
                    else {
                        extraBranchDesignator = UNK_BRANCH;
                    }
                }
                else {
                    isRecordOK = rg_FALSE;
                }
                break;
            }            
            case 2 : {
                string strRemoteIndicator  = strRemoteAndBranch.substr(0, 1);
                string strBranchDesignator = strRemoteAndBranch.substr(1, 1);
                
                if( isStringValidForRemoteIndicator( strRemoteIndicator )   &&          // REMOTENESS IN CHAR
                    isStringValidForBranchDesignator( strBranchDesignator )  ) {        // BRANCH IN NUM
                    remoteIndicator = convertStringIntoRemoteIndicator( strRemoteIndicator );
                    branchDesignator = convertStringIntoBranchDesignator( strBranchDesignator );

                    if (pos == 1) {
                        extraBranchDesignator = convertStringIntoBranchDesignator( atomName.substr(0, 1) );
                    }
                    else {
                        extraBranchDesignator = UNK_BRANCH;
                    }
                    // 
                    //if( targetAtom->getAtomCode() == H_ATOM ) {
                    //    extraBranchDesignator = convertStringIntoBranchDesignator( strBranchDesignator );
                    //}
                    //else {
                    //    branchDesignator = convertStringIntoBranchDesignator( strBranchDesignator );
                    //}
                }
                else {
                    isRecordOK = rg_FALSE;
                }
                break;
                     }
            case 3 : {
                string strRemoteIndicator       = strRemoteAndBranch.substr(0, 1);
                string strBranchDesignator      = strRemoteAndBranch.substr(1, 1);
                string strExtraBranchDesignator = strRemoteAndBranch.substr(2, 1);
                
                if( isStringValidForRemoteIndicator( strRemoteIndicator )     &&        // REMOTENESS IN CHAR
                    isStringValidForBranchDesignator( strBranchDesignator )   &&        // BRANCH IN NUM
                    isStringValidForBranchDesignator( strExtraBranchDesignator )  ) {   // BRANCH IN NUM
                    
                    remoteIndicator       = convertStringIntoRemoteIndicator( strRemoteIndicator );
                    branchDesignator      = convertStringIntoBranchDesignator( strBranchDesignator );
                    extraBranchDesignator = convertStringIntoBranchDesignator( strExtraBranchDesignator );
                }
                break;
                     }
            default: {
                remoteIndicator       = UNK_REMOTE;
                branchDesignator      = UNK_BRANCH;
                extraBranchDesignator = UNK_BRANCH;
                break;
                
                }
            }
        }

        targetAtom->getpChemicalProperties()->setRemoteIndicator( remoteIndicator );
        targetAtom->getpChemicalProperties()->setBrangeDesignator( branchDesignator );
        targetAtom->getpChemicalProperties()->setExtraBrangeDesignator( extraBranchDesignator );


        ResidueCode resCodeForTargetAtom = targetAtom->getResidue()->getResidueCode();
        
        if ( !targetAtom->getResidue()->isAminoResidue() )
            resCodeForTargetAtom = UNK_RESIDUE;
        
        // Set AtomTypeInAmber
        targetAtom->getpChemicalProperties()->setAtomTypeInAmber( AMBER_ATOM_TYPE[resCodeForTargetAtom]
                                                                                 [ATOM_FEATURES[targetAtom->getAtomCode()].IDOfAtomInStandardResidue ]
                                                                                 [remoteIndicator]
                                                                                 [branchDesignator]
                                                                                 [extraBranchDesignator] );
        
        targetAtom->getpChemicalProperties()->setChargeInAmber( ELECTROSTATIC_CHARGE_COEFFS_NON_BONDED_ATOM_PAIR[resCodeForTargetAtom]
                                                                                                                [ATOM_FEATURES[targetAtom->getAtomCode()].IDOfAtomInStandardResidue ]
                                                                                                                [remoteIndicator]
                                                                                                                [branchDesignator]
                                                                                                                [extraBranchDesignator] );
        
        // Check isOnBackBone
        rg_FLAG isAtomOnBackBone = rg_FALSE;
        
        // if residue is Amino acid
        if( targetAtom->getResidue()->isAminoResidue()) {
            if( atomName == "N" || atomName == "CA" || atomName == "C" )
                isAtomOnBackBone = rg_TRUE;
        }
        // else if residue is DNA or RNA
        else if( targetAtom->getResidue()->isNeucleicAcidResidue() ) {
            if( atomName == "P" || atomName == "O5*" || atomName == "C5*" || 
                atomName == "C4*" || atomName == "C3*" || atomName == "O3*"  )
                isAtomOnBackBone = rg_TRUE;
        }
        
        targetAtom->getpChemicalProperties()->setIsOnBackBone( isAtomOnBackBone );
    }

    else {
        return isRecordOK = rg_FALSE;
    }

    return isRecordOK;
}



rg_FLAG MoleculeIOFunctions::setBondRTIToMolecule( list<string>::iterator& i_recLines, list<string>* recordsOfMol2File, const rg_INT& numOfBondsGiven, AtomMap& mapOfAtom, Molecule& aMolecule )
{
    rg_FLAG isBondRTIOK = rg_TRUE;
    
    i_recLines++;
    string strRecLine  = StringFunctions::strTrim( *i_recLines );

    rg_INT i_bond = 0;
    
    while( !isNewTriposRTILineStart( strRecLine ) ) {
        
        i_bond++;
    
        rg_INT bondID       = 0;
        rg_INT originAtomID = 0;
        rg_INT targetAtomID = 0;
        string bondType     = "";
        string statusBit    = "";

        Atom* originAtom = rg_NULL;
        Atom* targetAtom = rg_NULL;

        extractBondRecordsFromBondRTILine( strRecLine, bondID, originAtomID, targetAtomID, bondType, statusBit );

        ChemicalBond* pNewChemicalBond = rg_NULL;

        AtomMap::iterator AtomMap_i = mapOfAtom.find( originAtomID );
        AtomMap::iterator AtomMap_j = mapOfAtom.find( targetAtomID );

        
        if ( AtomMap_i != mapOfAtom.end() && AtomMap_j != mapOfAtom.end() ) {
            originAtom = (*AtomMap_i).second;
            targetAtom = (*AtomMap_j).second;

            
                
            pNewChemicalBond = setChemicalBondToMoleculeWithoutDuplication( originAtom, targetAtom, aMolecule );

            if(pNewChemicalBond != rg_NULL) {
                setChemicalBondTypeForRTIBond(pNewChemicalBond, bondType);
            }
        }
        else {
            // Just ignore when the originAtom or targetAtom does not exist in Atom field.
            // isBondRTIOK = rg_FALSE;
            //break;
        }

        i_recLines++;
        
        if ( i_recLines != recordsOfMol2File->end() ) {
            strRecLine = StringFunctions::strTrim( *i_recLines );
        }
        else {
            break;
        }

        if( i_bond == numOfBondsGiven )
            break;

    }

    if( i_bond != numOfBondsGiven )
        isBondRTIOK = rg_FALSE;
    
    i_recLines--;

    return isBondRTIOK;
}



void MoleculeIOFunctions::extractBondRecordsFromBondRTILine( const string&  strRecLine, rg_INT& bondID, rg_INT& originAtomID, rg_INT& targetAtomID, string& bondType, string& statusBit )
{
    rg_INT stPosOfSingleRec = 0;
    rg_INT edPosOfSingleRec = 0;
    
    // BOND ID
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    bondID = atoi( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
    
    
    // ORIGINAL ATOM ID
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    originAtomID = atoi( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
    
    
    // TARGET ATOM ID
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    targetAtomID = atoi( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );


    // TYPE OF BOND
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    bondType = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
    

    //// STATUS BIT
    //stPosOfSingleRec = edPosOfSingleRec;
    //StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    //edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    //
    //statusBit = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
}



AmberAtomTypes MoleculeIOFunctions::getAmberAtomTypeFromSYBYL( const string& atomType )
{
    rg_INT posOfAtomTypeEnd = atomType.find(".");
    string typeOfAtom[2] = {"",""};
    if( posOfAtomTypeEnd == -1 ) {
        typeOfAtom[0] = StringFunctions::strToUpper( atomType );
    }
    else {
        typeOfAtom[0] = StringFunctions::strToUpper( atomType.substr(0, posOfAtomTypeEnd) );
        typeOfAtom[1] = StringFunctions::strToUpper( atomType.substr(posOfAtomTypeEnd + 1, atomType.length()) );
    }

    AmberAtomTypes amberAtomType = UNK_ATM;

    if(typeOfAtom[0] == "C") {    
        if(typeOfAtom[1] == "AR") {                
            amberAtomType = CA_ATM;
        }
        else if(typeOfAtom[1] == "2") {                    
            amberAtomType = CM_ATM;
        }
        else if(typeOfAtom[1] == "3") {
            amberAtomType = CT_ATM;
        }
        else {
            amberAtomType = C_ATM;
        }
    }
    else if(typeOfAtom[0] == "F") {
        amberAtomType = F_ATM;
    }
    else if(typeOfAtom[0] == "H") {
        if(typeOfAtom[1] == "T3P") {
            amberAtomType = HW_ATM;
        }
        else {
            amberAtomType = H_ATM;
        }
    }
    else if(typeOfAtom[0] == "LP") {
        amberAtomType = IP_ATM;
    }
    else if(typeOfAtom[0] == "K") {
        amberAtomType = K_ATM;
    }
    else if(typeOfAtom[0] == "LI") {
        amberAtomType = Li_ATM;
    }
    else if(typeOfAtom[0] == "N") {
        if(typeOfAtom[1] == "3" || typeOfAtom[1] == "4") {                
            amberAtomType = N3_ATM;
        }
        else {
            amberAtomType = N_ATM;
        }
    }
    else if(typeOfAtom[0] == "O") {
        if(typeOfAtom[1] == "2" || typeOfAtom[1] == "CO2") {                
            amberAtomType = O_ATM;
        }
        else if(typeOfAtom[1] == "CO2") {                    
            amberAtomType = O2_ATM;
        }
        else if(typeOfAtom[1] == "3") {
            amberAtomType = OH_ATM;
        }
        else if(typeOfAtom[1] == "T3P") {
            amberAtomType = OW_ATM;
        }
        else {
            amberAtomType = OS_ATM;
        }
    }
    else if(typeOfAtom[0] == "P") {
        amberAtomType = P_ATM;
    }
    else if(typeOfAtom[0] == "RB") {
        amberAtomType = Rb_ATM;
    }
    else if(typeOfAtom[0] == "S") {
        amberAtomType = S_ATM;
    }
    else {
        amberAtomType = UNK_ATM;
    }
    
    
    return amberAtomType;
}


void MoleculeIOFunctions::setChemicalBondTypeForRTIBond(ChemicalBond* aChemicalBond, const string& bondType)
{
    if(bondType == "1") {
        aChemicalBond->setTypeOfBond(SINGLE_BOND);
    }
    else if(bondType == "2") {
        aChemicalBond->setTypeOfBond(DOUBLE_BOND);
    }
    else if(bondType == "3") {
        aChemicalBond->setTypeOfBond(TRIPLE_BOND);
    }
    else if(bondType == "am") {
        aChemicalBond->setTypeOfBond(AMIDE_BOND);
    }
    else if(bondType == "ar") {
        aChemicalBond->setTypeOfBond(AROMATIC_BOND);
    }
    else if(bondType == "du") {
        aChemicalBond->setTypeOfBond(DUMMY_BOND);
    }
    else if(bondType == "un") {
        aChemicalBond->setTypeOfBond(UNK_BOND);
    }
}


void MoleculeIOFunctions::setAmberAtomTypeForHAtom(Molecule& aMolecule)
{
    rg_dList<Atom>* atoms = aMolecule.getAtoms();

    atoms->reset4Loop();
    while(atoms->setNext4Loop()) {
        Atom* currAtom = atoms->getpEntity();


        if(currAtom->getChemicalProperties().getAtomTypeInAmber() != H_ATM) {
            continue;
        }
        
        rg_dList<ChemicalBond*>* chemicalBonds = currAtom->getListChemicalBond();

        if(chemicalBonds->getSize() == 0) {
            continue;
        }

        ChemicalBond* firstBond = chemicalBonds->getFirstEntity();

        Atom* oppositeAtom = rg_NULL;

        if(firstBond->getFirstAtom() == currAtom) {
            oppositeAtom = firstBond->getSecondAtom();
        }
        else if(firstBond->getSecondAtom() == currAtom) {
            oppositeAtom = firstBond->getFirstAtom();
        }

        AmberAtomTypes amberAtomTypeForOppositeAtom = oppositeAtom->getChemicalProperties().getAtomTypeInAmber();

        if(amberAtomTypeForOppositeAtom == CA_ATM) {
            currAtom->getpChemicalProperties()->setAtomTypeInAmber(HA_ATM);
        }
        else if(amberAtomTypeForOppositeAtom == C_ATM || amberAtomTypeForOppositeAtom == CM_ATM || amberAtomTypeForOppositeAtom == CT_ATM ) {
            currAtom->getpChemicalProperties()->setAtomTypeInAmber(HC_ATM);
        }        
        else if(amberAtomTypeForOppositeAtom == OH_ATM) {
            currAtom->getpChemicalProperties()->setAtomTypeInAmber(HO_ATM);
        }
        else if(amberAtomTypeForOppositeAtom == S_ATM) {
            currAtom->getpChemicalProperties()->setAtomTypeInAmber(HS_ATM);
        }
    }
}

void MoleculeIOFunctions::setChainsForNonStdResiduesForPDBFile( ResidueMap* mapOfNonStdResidue, Molecule& aMolecule )
{

    Chain* newChainForHOH = aMolecule.addChain( Chain(aMolecule.getChains()->getSize(), &aMolecule) );

    ResidueMap::iterator residueMap_i = mapOfNonStdResidue->begin();

    while ( residueMap_i != mapOfNonStdResidue->end() ) {
        Residue* currNonStdResidue = (*residueMap_i).second;
        rg_dList<Residue*>* residuesInOldChain = currNonStdResidue->getChain()->getResidues();

        

        if ( currNonStdResidue->getSequenceNumber() == 121 ) {
                int aaa = 0;
        }


        // Delete NonStdResidue ptr in Chain.
        residuesInOldChain->reset4Loop();
        while ( residuesInOldChain->setNext4Loop() ) {
            Residue* aResidueInOldChain = residuesInOldChain->getEntity();

            if ( aResidueInOldChain->getResidueCode() == HOH_RESIDUE &&
                aResidueInOldChain->getChain()->getChainIDFromInputFileInDecimal() == 66 ) {
                int aaa = 0;
            }

            if( currNonStdResidue == aResidueInOldChain ) {
                residuesInOldChain->killCurrent();
                break;
            }
        }
    
        // Filter HOH residue and set Chain ptr to Residue
        Chain* newChain = rg_NULL;
        if( currNonStdResidue->getResidueCode() == HOH_RESIDUE ) {
            newChain = newChainForHOH;
        }
        else {
            newChain = aMolecule.addChain( Chain(aMolecule.getChains()->getSize(), &aMolecule) );
        }
        
        newChain->addResidue( currNonStdResidue );
        currNonStdResidue->setChain( newChain );
        
        residueMap_i++;
    }


    // Filter chains with no residue and reset ID.
    filterEmptyChain( aMolecule );
}



void MoleculeIOFunctions::setChainsForNonStdResiduesForPDBFile( rg_dList<ResidueMapWithChainID>* mapsOfNonStdResidue, Molecule& aMolecule )
{
    Chain* newChainForHOH = aMolecule.addChain( Chain(aMolecule.getChains()->getSize(), &aMolecule) );
    
    mapsOfNonStdResidue->reset4Loop();
    while ( mapsOfNonStdResidue->setNext4Loop() ) {
        
        ResidueMap* mapOfNonStdResidue = &(mapsOfNonStdResidue->getpEntity()->second);

        ResidueMap::iterator residueMap_i = mapOfNonStdResidue->begin();

        while ( residueMap_i != mapOfNonStdResidue->end() ) {
            Residue* currNonStdResidue = (*residueMap_i).second;
            rg_dList<Residue*>* residuesInOldChain = currNonStdResidue->getChain()->getResidues();

            // Delete NonStdResidue ptr in Chain.
            residuesInOldChain->reset4Loop();
            while ( residuesInOldChain->setNext4Loop() ) {
                Residue* aResidueInOldChain = residuesInOldChain->getEntity();

                if( currNonStdResidue == aResidueInOldChain ) {
                    residuesInOldChain->killCurrent();
                    break;
                }
            }
    
            // Filter HOH residue and set Chain ptr to Residue
            Chain* newChain = rg_NULL;
            if( currNonStdResidue->getResidueCode() == HOH_RESIDUE ) {
                newChain = newChainForHOH;
            }
            else {
                newChain = aMolecule.addChain( Chain(aMolecule.getChains()->getSize(), &aMolecule) );
            }
        
            newChain->addResidue( currNonStdResidue );
            currNonStdResidue->setChain( newChain );
        
            residueMap_i++;
        }
    }


    // Filter chains with no residue and reset ID.
    filterEmptyChain( aMolecule );
}



void MoleculeIOFunctions::setHelixRecordsToMolecule( list<string*>* helixRecLinesOfPDBFile, rg_dList<ResidueMapWithChainID>* mapsOfResidue, Molecule* aMolecule )
{
    rg_dList<Residue*> residuesSortedBySequence = rg_NULL;
    aMolecule->getResiduesSortedBySequenceNumber( residuesSortedBySequence );

    list<string*>::iterator i_recLines = helixRecLinesOfPDBFile->begin();
    while( i_recLines != helixRecLinesOfPDBFile->end() ) {
        string* currRecLine = *i_recLines;
        setHelixRecordsToMolecule( *currRecLine, mapsOfResidue, &residuesSortedBySequence );
        i_recLines++;        
    }
}



void MoleculeIOFunctions::setHelixRecordsToMolecule( const string& strRecLine, rg_dList<ResidueMapWithChainID>* mapsOfResidue, rg_dList<Residue*>* residuesSortedBySequence )
{
    rg_INT serial      = atoi( StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_SERIAL_ST_POS, PDB_HELIX_RECORD_SERIAL_LENGTH ).c_str() );
    string helixID     = StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_HELIX_ID_ST_POS, PDB_HELIX_RECORD_HELIX_ID_LENGTH );
    string initChainID = StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_INIT_CHAIN_ID_ST_POS, PDB_HELIX_RECORD_INIT_CHAIN_ID_LENGTH );
    rg_INT initResSeq  = atoi( StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_INIT_RES_SEQ_ST_POS, PDB_HELIX_RECORD_INIT_RES_SEQ_LENGTH ).c_str() );
    string endChainID  = StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_END_CHAIN_ID_ST_POS, PDB_HELIX_RECORD_HELIX_CLASS_LENGTH );
    rg_INT endResSeq   = atoi( StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_END_RES_SEQ_ST_POS, PDB_HELIX_RECORD_END_RES_SEQ_LENGTH ).c_str() );
    rg_INT helixClass  = atoi( StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_HELIX_CLASS_ST_POS, PDB_HELIX_RECORD_HELIX_CLASS_LENGTH ).c_str() );
    string comments    = StringFunctions::subString( strRecLine, PDB_HELIX_RECORD_COMMENTS_ST_POS, PDB_HELIX_RECORD_COMMENTS_LENGTH );


    rg_INT initChainIDInDecimal = (int)initChainID[0];
    rg_INT endChainIDInDecimal  = (int)endChainID[0];

    // Find residues with initResidue and endResidue from mapsOfResidue
    Residue* startResidue = findResidueFromMapsOfResidue( initResSeq, initChainIDInDecimal, mapsOfResidue );
    Residue* endResidue   = findResidueFromMapsOfResidue( endResSeq,  endChainIDInDecimal,  mapsOfResidue );

    if( startResidue == rg_NULL || endResidue == rg_NULL )
        return;

    rg_dList<Residue*> residuesOfHelix;
    findResiduesFromListSortedBySequence( residuesSortedBySequence, startResidue, endResidue, residuesOfHelix );

    // Add a new helix to secondary structure
    SecondaryStructure* secStructForInitResidue = startResidue->getChain()->getSecondaryStructure();
    secStructForInitResidue->addHelix( Helix( serial, helixID, helixClass, &residuesOfHelix, comments ) );
}



void MoleculeIOFunctions::setSheetRecordsToMolecule( list<string*>* sheetRecLinesOfPDBFile, rg_dList<ResidueMapWithChainID>* mapsOfResidue, Molecule* aMolecule )
{
    rg_dList<Residue*> residuesSortedBySequence = rg_NULL;
    aMolecule->getResiduesSortedBySequenceNumber( residuesSortedBySequence );

    list<string*>::iterator i_recLines = sheetRecLinesOfPDBFile->begin();
    while( i_recLines != sheetRecLinesOfPDBFile->end() ) {
        string* currRecLine = *i_recLines;
        setSheetRecordsToMolecule( *currRecLine, mapsOfResidue, &residuesSortedBySequence );
        i_recLines++;        
    }
}



void MoleculeIOFunctions::setSheetRecordsToMolecule( const string& strRecLine, rg_dList<ResidueMapWithChainID>* mapsOfResidue, rg_dList<Residue*>* residuesSortedBySequence )
{
    rg_INT strandSerial = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_STRAND_SERIAL_ST_POS, PDB_SHEET_RECORD_STRAND_SERIAL_LENGTH ).c_str() );
    string sheetID      = StringFunctions::strTrim( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_SHEET_ID_ST_POS, PDB_SHEET_RECORD_SHEET_ID_LENGTH ) );
    rg_INT numOfStrands = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_NUM_STRANDS_ST_POS, PDB_SHEET_RECORD_NUM_STRANDS_LENGTH ).c_str() );

    string initChainID  = StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_INIT_CHAIN_ID_ST_POS, PDB_SHEET_RECORD_INIT_CHAIN_ID_LENGTH );
    rg_INT initResSeq   = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_INIT_RES_SEQ_ST_POS, PDB_SHEET_RECORD_INIT_RES_SEQ_LENGTH ).c_str() );
    string endChainID   = StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_END_CHAIN_ID_ST_POS, PDB_SHEET_RECORD_END_CHAIN_ID_LENGTH );
    rg_INT endResSeq    = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_END_RES_SEQ_ST_POS, PDB_SHEET_RECORD_END_RES_SEQ_LENGTH ).c_str() );
    
    rg_INT senseParallel = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_PARALLEL_ST_POS, PDB_SHEET_RECORD_PARALLEL_LENGTH ).c_str() );
    
    string currHBAtomName = StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_CURR_ATOM_NAME_ST_POS, PDB_SHEET_RECORD_CURR_ATOM_NAME_LENGTH );
    string currHBChainID  = StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_CURR_CHAIN_ID_ST_POS, PDB_SHEET_RECORD_CURR_CHAIN_ID_LENGTH );
    rg_INT currHBResSeq   = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_CURR_RES_SEQ_ST_POS, PDB_SHEET_RECORD_CURR_RES_SEQ_LENGTH ).c_str() );
    
    string prevHBAtomName = StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_PREV_ATOM_NAME_ST_POS, PDB_SHEET_RECORD_PREV_ATOM_NAME_LENGTH );
    string prevHBChainID  = StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_PREV_CHAIN_ID_ST_POS, PDB_SHEET_RECORD_PREV_CHAIN_ID_LENGTH );
    rg_INT prevHBResSeq   = atoi( StringFunctions::subString( strRecLine, PDB_SHEET_RECORD_PREV_RES_SEQ_ST_POS, PDB_SHEET_RECORD_PREV_RES_SEQ_LENGTH ).c_str() );


    // Find residues with initResidue and endResidue from mapsOfResidue
    rg_INT initChainIDInDecimal = (int)initChainID[0];
    rg_INT endChainIDInDecimal  = (int)endChainID[0];

    Residue* startResidue = findResidueFromMapsOfResidue( initResSeq, initChainIDInDecimal, mapsOfResidue );
    Residue* endResidue   = findResidueFromMapsOfResidue( endResSeq,  endChainIDInDecimal,  mapsOfResidue );

    if( startResidue == rg_NULL || endResidue == rg_NULL )
        return;

    // Add a new sheet or get an existing sheet from secStructForInitResidue.
    SecondaryStructure* secStructForInitResidue = startResidue->getChain()->getSecondaryStructure();
    
    Sheet* targetSheet = secStructForInitResidue->getSheet( sheetID );
    if ( targetSheet == rg_NULL ) {
        targetSheet = secStructForInitResidue->addSheet( Sheet( sheetID, numOfStrands ) );
    }
    targetSheet->setSenseOfStrandBySerial( strandSerial, senseParallel );


    // Get targetStrand from targetSheet and set properties.
    Strand* targetStrand = targetSheet->getStrandBySerial( strandSerial );
    if ( targetStrand == rg_NULL )
        return;

    rg_dList<Residue*> residuesOfStrand;
    findResiduesFromListSortedBySequence( residuesSortedBySequence, startResidue, endResidue, residuesOfStrand );

    targetStrand->setSerial( strandSerial );
    targetStrand->setSheet( targetSheet );
    targetStrand->setResidues( &residuesOfStrand );


    // Set Hydrogen bond atoms to targetSheet.
    if ( strandSerial < 2 )
        return;

    Atom* currHBAtom = findAtomFromMapsOfResidue( currHBAtomName, currHBChainID, currHBResSeq, mapsOfResidue );
    Atom* prevHBAtom = findAtomFromMapsOfResidue( prevHBAtomName, prevHBChainID, prevHBResSeq, mapsOfResidue );
    
    if ( currHBAtom == rg_NULL || prevHBAtom == rg_NULL )
        return;

    HydrogenBond startHBondBetweenStrands;
    rg_BOOL isOrdinaryHBond = startHBondBetweenStrands.evaluateAndSetDonorAndAcceptorForProteinAtoms( currHBAtom, prevHBAtom );
    
    if ( isOrdinaryHBond == rg_FALSE ) 
        startHBondBetweenStrands.setDonorAndAcceptor( currHBAtom, prevHBAtom );

    targetSheet->setHBondsBetweenStrandsBySerial( strandSerial, startHBondBetweenStrands );
}



void MoleculeIOFunctions::setTurnRecordsToMolecule( list<string*>* turnRecLinesOfPDBFile, rg_dList<ResidueMapWithChainID>* mapsOfResidue, Molecule* aMolecule )
{
    rg_dList<Residue*> residuesSortedBySequence = rg_NULL;
    aMolecule->getResiduesSortedBySequenceNumber( residuesSortedBySequence );

    list<string*>::iterator i_recLines = turnRecLinesOfPDBFile->begin();
    while( i_recLines != turnRecLinesOfPDBFile->end() ) {
        string* currRecLine = *i_recLines;
        setTurnRecordsToMolecule( *currRecLine, mapsOfResidue, &residuesSortedBySequence );
        i_recLines++;
    }
}



void MoleculeIOFunctions::setTurnRecordsToMolecule( const string& strRecLine, rg_dList<ResidueMapWithChainID>* mapsOfResidue, rg_dList<Residue*>* residuesSortedBySequence )
{
    rg_INT serial      = atoi( StringFunctions::subString( strRecLine, PDB_TURN_RECORD_SERIAL_ST_POS, PDB_TURN_RECORD_SERIAL_LENGTH ).c_str() );
    string turnID     = StringFunctions::subString( strRecLine, PDB_TURN_RECORD_TURN_ID_ST_POS, PDB_TURN_RECORD_TURN_ID_LENGTH );
    string initChainID = StringFunctions::subString( strRecLine, PDB_TURN_RECORD_INIT_CHAIN_ID_ST_POS, PDB_TURN_RECORD_INIT_CHAIN_ID_LENGTH );
    rg_INT initResSeq  = atoi( StringFunctions::subString( strRecLine, PDB_TURN_RECORD_INIT_RES_SEQ_ST_POS, PDB_TURN_RECORD_INIT_RES_SEQ_LENGTH ).c_str() );
    string endChainID  = StringFunctions::subString( strRecLine, PDB_TURN_RECORD_END_CHAIN_ID_ST_POS, PDB_TURN_RECORD_END_CHAIN_ID_LENGTH );
    rg_INT endResSeq   = atoi( StringFunctions::subString( strRecLine, PDB_TURN_RECORD_END_RES_SEQ_ST_POS, PDB_TURN_RECORD_END_RES_SEQ_LENGTH ).c_str() );
    string comments    = StringFunctions::subString( strRecLine, PDB_TURN_RECORD_COMMENTS_ST_POS, PDB_TURN_RECORD_COMMENTS_LENGTH );


    rg_INT initChainIDInDecimal = (int)initChainID[0];
    rg_INT endChainIDInDecimal  = (int)endChainID[0];

    // Find residues with initResidue and endResidue from mapsOfResidue
    Residue* startResidue = findResidueFromMapsOfResidue( initResSeq, initChainIDInDecimal, mapsOfResidue );
    Residue* endResidue   = findResidueFromMapsOfResidue( endResSeq,  endChainIDInDecimal,  mapsOfResidue );

    if( startResidue == rg_NULL || endResidue == rg_NULL )
        return;

    rg_dList<Residue*> residuesOfTurn;
    findResiduesFromListSortedBySequence( residuesSortedBySequence, startResidue, endResidue, residuesOfTurn );

    // Add a new turn to secondary structure
    SecondaryStructure* secStructForInitResidue = startResidue->getChain()->getSecondaryStructure();
    secStructForInitResidue->addTurn( Turn( serial, turnID, &residuesOfTurn, comments ) );
}



Residue* MoleculeIOFunctions::findResidueFromMapsOfResidue( const rg_INT& residueSequence, const rg_INT& chainIDinDecimal, rg_dList<ResidueMapWithChainID>* mapsOfResidue )
{
    Residue* targetResidue = rg_NULL;
    ResidueMapWithChainID* pResidueMapWithChainID = rg_NULL;

    mapsOfResidue->reset4Loop();
    while ( mapsOfResidue->setNext4Loop() ) {
        pResidueMapWithChainID = mapsOfResidue->getpEntity();

        ResidueMap* aMapOfResidue = &(pResidueMapWithChainID->second);

//// aMapOfResidue : HETATM is not considered for Chain.
//
//         ResidueMap::iterator residueMap_i = aMapOfResidue->begin();
//        
//         Chain* currChain = (*residueMap_i).second->getChain();
//         rg_INT chainIDOFirstResidue = currChain->getChainIDFromInputFileInDecimal();       
// 
//         if( chainIDOFirstResidue == chainIDinDecimal ) {
//             ResidueMap::iterator residueMap_i = aMapOfResidue->find( residueSequence );
//             
//              if ( residueMap_i != aMapOfResidue->end() ) {
//                 targetResidue = (*residueMap_i).second;
//                 break;
//              }
//         }

        ResidueMap::iterator residueMap_i = aMapOfResidue->find( residueSequence );
    
        if ( residueMap_i != aMapOfResidue->end() ) {
            if ( (*residueMap_i).second->getChain()->getChainIDFromInputFileInDecimal() == chainIDinDecimal ) {
                targetResidue = (*residueMap_i).second;
                break;
            }
        }
    }

    return targetResidue;
}



void MoleculeIOFunctions::findResiduesFromListSortedBySequence( rg_dList<Residue*>* residuesSortedBySequence, Residue* startResidue, Residue* endResidue, rg_dList<Residue*>& targetResidues )
{
    targetResidues.removeAll();

    residuesSortedBySequence->reset4Loop();
    while ( residuesSortedBySequence->setNext4Loop() ) {
        Residue* currResidue = residuesSortedBySequence->getEntity();
        if ( currResidue == startResidue ) {
            while ( currResidue != endResidue ) {
                targetResidues.addTail( currResidue );
                residuesSortedBySequence->setNext4Loop();
                currResidue = residuesSortedBySequence->getEntity();
            }
            targetResidues.addTail( currResidue );
            break;
        }
    }
}


Atom* MoleculeIOFunctions::findAtomFromMapsOfResidue( const string& atomName, const string& chainID, const rg_INT& residueSeqence, rg_dList<ResidueMapWithChainID>* mapsOfResidue )
{
    rg_INT chainIDInDecimal  = (int)chainID[0];

    Residue* targetResidue = findResidueFromMapsOfResidue( residueSeqence, chainIDInDecimal, mapsOfResidue );
 
    if ( targetResidue == rg_NULL ) {
        return rg_NULL;
    }
    else {
        return targetResidue->getAtom( atomName );
    }
}



void MoleculeIOFunctions::setChainsForNonStdResiduesForMol2File( ResidueMap* mapOfNonStdResidue, Molecule& aMolecule )
{
    ResidueMap::iterator residueMap_i = mapOfNonStdResidue->begin();

    while ( residueMap_i != mapOfNonStdResidue->end() ) {
        Residue* currNonStdResidue = (*residueMap_i).second;

        rg_INT numOfChains = aMolecule.getChains()->getSize();
        Chain* newChain    = aMolecule.addChain( Chain(numOfChains, &aMolecule) );

        newChain->addResidue( currNonStdResidue );
        currNonStdResidue->setChain( newChain );

        residueMap_i++;
    }
}



void MoleculeIOFunctions::filterEmptyChain( Molecule& aMolecule )
{
    rg_dList<Chain>* chains = aMolecule.getChains();

    chains->reset4Loop();
    while ( chains->setNext4Loop() ) {
        if( chains->getpEntity()->getResidues()->getSize() == 0 )
            chains->killCurrent();
    }
}



void MoleculeIOFunctions::evaluateAndSetChainCode( Molecule& aMolecule )
{
    rg_dList<Chain>* chainList = aMolecule.getChains();

    chainList->reset4Loop();
    while ( chainList->setNext4Loop() ) {
        Chain* currChain = chainList->getpEntity();
        currChain->evaluateAndSetChainCode();
    }
}



void MoleculeIOFunctions::evaluateAndSetResidueCodesForNeucleicAcid( ResidueSymbolMap* mapOfResidueSymbolForDNA, ResidueSymbolMap* mapOfResidueSymbolForRNA, Molecule& aMolecule )
{
    rg_dList<Chain>* chainList = aMolecule.getChains();

    chainList->reset4Loop();
    while ( chainList->setNext4Loop() ) {
        Chain* currChain = chainList->getpEntity();
        if ( currChain->getChainCode() == DNA_CHAIN ) {

            rg_dList<Residue*>* residues = currChain->getResidues();

            residues->reset4Loop();
            while ( residues->setNext4Loop() ) {
                Residue* currResidue = residues->getEntity();
                
                ResidueSymbolMap::iterator ResidueSymbolMap_i = mapOfResidueSymbolForDNA->find( currResidue->getResidueName() );
        
                if ( ResidueSymbolMap_i != mapOfResidueSymbolForDNA->end() )
                    currResidue->setResidueCode( (ResidueCode)((*ResidueSymbolMap_i).second) );
            }

        }
        else if ( currChain->getChainCode() == RNA_CHAIN ) {
            
            rg_dList<Residue*>* residues = currChain->getResidues();

            residues->reset4Loop();
            while ( residues->setNext4Loop() ) {
                Residue* currResidue = residues->getEntity();
                
                ResidueSymbolMap::iterator ResidueSymbolMap_i = mapOfResidueSymbolForRNA->find( currResidue->getResidueName() );
        
                if ( ResidueSymbolMap_i != mapOfResidueSymbolForRNA->end() )
                    currResidue->setResidueCode( (ResidueCode)((*ResidueSymbolMap_i).second) );
            }
        }
    }
}



// void MoleculeIOFunctions::checkChainsInMoleculeForStandardResidue( Molecule& aMolecule )
// {
//     rg_dList<Chain>* chainList = aMolecule.getChains();
// 
//     chainList->reset4Loop();
//     while ( chainList->setNext4Loop() ) {
//         Chain* currChain = chainList->getpEntity();
//         
//         rg_dList<Residue*>* residueListInChain = currChain->getResidues();
//         
//         rg_INT numOfAllResidue         = residueListInChain->getSize();
//         rg_INT numOfStandardResidue    = 0;
//         rg_INT numOfNonStandardResidue = 0;
//         rg_INT numOfHOH                = 0;
// 
//         residueListInChain->reset4Loop();
//         while ( residueListInChain->setNext4Loop() ) {
//             Residue* currResidue = residueListInChain->getEntity();
//             
//             if ( currResidue->isStandardResidue() )
//                 numOfStandardResidue++;
//             else
//                 numOfNonStandardResidue++;
// 
//             if ( currResidue->getResidueCode() == HOH_RESIDUE )
//                 numOfHOH++;
//         }
// 
//         if ( numOfStandardResidue == 0 ) {
//             if ( numOfNonStandardResidue > 0 ) {
//                 currChain->setChainCode( NON_STD_RES_CHAIN );
//             }
//             if ( numOfAllResidue == numOfHOH ) {
//                 currChain->setChainCode( HOH_CHAIN );
//             }
//         }
//         else {
//             currChain->setChainCode( STD_RES_CHAIN );
//         }       
//     }
// }



rg_FLAG MoleculeIOFunctions::isStringValidForAtomType( const string& aString )
{
    rg_FLAG isValid = rg_FALSE;
    
    if( (int)aString[0] >= 'A' && (int)aString[0] <= 'Z' ) {
        isValid = rg_TRUE;
    }
    
    return isValid;
}



rg_FLAG MoleculeIOFunctions::isStringValidForRemoteIndicator( const string& aString )
{
    rg_FLAG isValid = rg_FALSE;
    
    if( (int)aString[0] >= 'A' && (int)aString[0] <= 'Z' ) {
        isValid = rg_TRUE;
    }
    
    return isValid;
}



rg_FLAG MoleculeIOFunctions::isStringValidForBranchDesignator( const string& aString )
{
    rg_FLAG isValid = rg_FALSE;
    
    if( (int)aString[0] == '1' || (int)aString[0] == '2' || (int)aString[0] == '3' || (int)aString[0] == 'T' ) {
        isValid = rg_TRUE;
    }
    
    return isValid;
}



rg_FLAG MoleculeIOFunctions::isStringValidForExtraBranchDesignator( const string& aString )
{
    rg_FLAG isValid = rg_FALSE;
    
    if( (int)aString[0] == '1' || (int)aString[0] == '2' || (int)aString[0] == '3' ) {
        isValid = rg_TRUE;
    }
    
    return isValid;
}



// AmberAtomTypes MoleculeIOFunctions::setAtomsTypeInAmber();
// {
//         // Check ResidueCode of targetAtom for AtomTypeInAmber : AtomTypeInAmber works only for Amino acid !
//         ResidueCode resCodeForTargetAtom = targetAtom->getResidue()->getResidueCode();
//         
//         if ( (rg_INT)resCodeForTargetAtom > 22 )
//             resCodeForTargetAtom = UNK_RESIDUE;
//         
//         // Set AtomTypeInAmber
//         targetAtom->getpChemicalProperties()->setAtomTypeInAmber( AMBER_ATOM_TYPE[resCodeForTargetAtom]
//             [ATOM_FEATURES[convertedAtomCode].IDOfAtomInStandardResidue ]
//             [remoteIndicator]
//             [brangeDesignator]
//             [extraBrangeDesignator] );
//         
//         targetAtom->getpChemicalProperties()->setCharge( ELECTROSTATIC_CHARGE_COEFFS_NON_BONDED_ATOM_PAIR[resCodeForTargetAtom]
//             [ATOM_FEATURES[convertedAtomCode].IDOfAtomInStandardResidue ]
//             [remoteIndicator]
//             [brangeDesignator]
//             [extraBrangeDesignator] );
//         
//         // Check isOnBackBone
//         rg_FLAG isAtomOnBackBone = rg_FALSE;
//         
//         // if residue is Amino acid
//         if( (rg_INT)resCodeForTargetAtom > 0 && (rg_INT)resCodeForTargetAtom < 23 ) {
//             if( atomName == " N  " || atomName == " CA " || atomName == " C  " )
//                 isAtomOnBackBone = rg_TRUE;
//         }
//         // else if residue is DNA or RNA
//         else if( (rg_INT)resCodeForTargetAtom > 22 && (rg_INT)resCodeForTargetAtom < 43 ) {
//             if( atomName == " P  " || atomName == " O5*" || atomName == " C5*" || 
//                 atomName == " C4*" || atomName == " C3*" || atomName == " O3*"  )
//                 isAtomOnBackBone = rg_TRUE;
//         }
//         
//         targetAtom->getpChemicalProperties()->setIsOnBackBone( isAtomOnBackBone );


//     const rg_INT UNK_POS = -1;
//     
//     for ( int i=0; i<m_numAtoms; i++ )  {
//         rg_INT atomicNumber = m_atoms[i].getAtomicNumber();
//         string atomType     = m_atoms[i].getAtomType();
//         rg_INT posOfDot     = atomType.find(".");
//         rg_INT lengthOfType = atomType.length();
// 
//         if ( atomicNumber ==  6 ) {         // C
//             if ( posOfDot != UNK_POS ) {
//                 string bondType = atomType.substr( posOfDot+1, lengthOfType-posOfDot );
//                 if ( bondType == "1" )
//                     m_atoms[i].setAtomTypeInAmber( CT_ATM );
//                 else if ( bondType == "2" )
//                     m_atoms[i].setAtomTypeInAmber( CM_ATM );
//                 else if ( bondType == "ar" )
//                     m_atoms[i].setAtomTypeInAmber( CA_ATM );
//                 else 
//                     m_atoms[i].setAtomTypeInAmber(  C_ATM );
//             }
//             else // if ( posOfDot == UNK_POS )
//                 m_atoms[i].setAtomTypeInAmber(  C_ATM );
//         }
//         else if ( atomicNumber == 9   ) {   // F
//             m_atoms[i].setAtomTypeInAmber(  F_ATM );
//         }      
//         else if ( atomicNumber == 1   ) {   // H
//             if ( posOfDot == UNK_POS ) {
//                 for ( int j=0; j<m_numBonds; j++ ) {
//                     TriposAtom* origianlAtom = m_bonds[j].getOriginAtom();
//                     TriposAtom* targetAtom   = m_bonds[j].getTargetAtom();
//                     string      bondedAtomType;
//                     
//                     if ( m_atoms[i].getAtomID() == origianlAtom->getAtomID() )
//                         bondedAtomType = targetAtom->getAtomType();
//                     else if ( m_atoms[i].getAtomID() == targetAtom->getAtomID() )
//                         bondedAtomType = origianlAtom->getAtomType();
//                     else 
//                         continue;
// 
//                     rg_INT posOfDotInBondedAtom  = bondedAtomType.find(".");
//                     string bondedAtomSymbol;
// 
//                     if ( posOfDotInBondedAtom != UNK_POS )
//                         bondedAtomSymbol = bondedAtomType.substr( 0, posOfDotInBondedAtom-1 );
//                     else
//                         bondedAtomSymbol = bondedAtomType;
// 
//                     if ( bondedAtomSymbol == "N" ) {
//                         m_atoms[i].setAtomTypeInAmber( H_ATM );
//                         break;
//                     }
//                     else if ( bondedAtomSymbol == "O" ) {
//                         m_atoms[i].setAtomTypeInAmber( HO_ATM );
//                         break;
//                     }
//                     else if ( bondedAtomSymbol == "S" ) {
//                         m_atoms[i].setAtomTypeInAmber( HS_ATM );
//                         break;
//                     }
//                     else if ( bondedAtomSymbol == "C" && bondedAtomType == "C.ar" ) {
//                         m_atoms[i].setAtomTypeInAmber( HA_ATM );
//                         break;
//                     }
//                     else { 
//                         m_atoms[i].setAtomTypeInAmber( HC_ATM );
//                         break;
//                     }
//                 }
//             }
//             else {
//                 string bondType = atomType.substr( posOfDot+1, lengthOfType-posOfDot );
//                 if ( bondType == "t3p" )
//                     m_atoms[i].setAtomTypeInAmber( HW_ATM );
//                 else
//                     m_atoms[i].setAtomTypeInAmber( HC_ATM );
//             }
//         }
//         else if ( atomicNumber == 104 ) {   // LP
//             m_atoms[i].setAtomTypeInAmber( IP_ATM );
//         }
//         else if ( atomicNumber == 19  ) {   // k
//             m_atoms[i].setAtomTypeInAmber(  K_ATM );
//         }
//         else if ( atomicNumber == 3   ) {   // LI
//             m_atoms[i].setAtomTypeInAmber( Li_ATM );
//         }
//         else if ( atomicNumber == 7   ) {   // N
//             if ( posOfDot != UNK_POS ) {
//                 string bondType = atomType.substr( posOfDot+1, lengthOfType-posOfDot );
//                 if ( bondType == "4" )
//                     m_atoms[i].setAtomTypeInAmber( N3_ATM );
//                 else
//                     m_atoms[i].setAtomTypeInAmber(  N_ATM );
//             }
//             else // if ( posOfDot == UNK_POS )
//                 m_atoms[i].setAtomTypeInAmber(  N_ATM );
//         }
//         else if ( atomicNumber == 8   ) {   // O
//             if ( posOfDot != UNK_POS ) {
//                 string bondType = atomType.substr( posOfDot+1, lengthOfType-posOfDot );
//                 if ( bondType == "2" )
//                     m_atoms[i].setAtomTypeInAmber(  O_ATM );
//                 else if ( bondType == "3" )
//                     m_atoms[i].setAtomTypeInAmber( OS_ATM );
//                 else if ( bondType == "co2" )
//                     m_atoms[i].setAtomTypeInAmber( O2_ATM );
//                 else if ( bondType == "t3p" )
//                     m_atoms[i].setAtomTypeInAmber( OW_ATM );
//                 else 
//                     m_atoms[i].setAtomTypeInAmber( OH_ATM );
//             }
//             else // if ( posOfDot == UNK_POS )
//                 m_atoms[i].setAtomTypeInAmber( OH_ATM );
//         }
//         else if ( atomicNumber == 15  ) {   // P
//             m_atoms[i].setAtomTypeInAmber(  P_ATM );
//         }
//         else if ( atomicNumber == 37  ) {   // RB
//             m_atoms[i].setAtomTypeInAmber( Rb_ATM );
//         }
//         else if ( atomicNumber == 16  ) {   // S
//             m_atoms[i].setAtomTypeInAmber( S_ATM );
//         }
//     }
// }



rg_FLAG MoleculeIOFunctions::writePDBFile( const char* pathFile, Molecule* aMolecule )
{
    rg_FLAG isFileWritten   = rg_FALSE;
    
    ofstream  fout( pathFile );

    writeHeader( aMolecule, fout );
    writeAtomField( aMolecule, fout );
    //writeConnectField( aMolecule, fout );
    writeHeteroAtomField( aMolecule, fout );
    
    fout.close();


    return isFileWritten;
}



void    MoleculeIOFunctions::writeHeader( Molecule* aMolecule, ofstream& fout )
{
    
}



void    MoleculeIOFunctions::writeAtomField(Molecule* aMolecule, ofstream& fout)
{
    list<Atom*> atoms;

    rg_dList<Atom>* atomsInMolecule = aMolecule->getAtoms();
    atomsInMolecule->reset4Loop();
    while (atomsInMolecule->setNext4Loop()) {
        Atom* currAtom = atomsInMolecule->getpEntity();

        atoms.push_back(currAtom);
    }

    writeAtomFieldInPDBFormat(atoms, fout);
}



void    MoleculeIOFunctions::writeAtomFieldInPDBFormat(const list<Atom*>& atoms, ofstream& fout)
{
    fout.setf(ios::fixed, ios::floatfield);     // TO 3 PLACES OF DECIMALS 
    fout.precision(3);                          //

    for ( list<Atom*>::const_iterator i_atom=atoms.begin(); i_atom!=atoms.end(); ++i_atom) {
        Atom* currAtom = (*i_atom);

		// for SCP solver
		//if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
		//	continue;
        
        string recLine;
        string blankSpace  = "";
        string singleSpace = " ";

        string serial;
        string name;
        string altLoc;
        string resName;
        string chainID;
        string resSeq;
        string xCoord;
        string yCoord;
        string zCoord;
        string occupancy;
        string tempFactor;
        string element;
        string charge;

        if (currAtom->getResidue()->isStandardResidue()) {
            recLine = PDB_RECORD_TYPE_ATOM;
        }
        else {
            recLine = PDB_RECORD_TYPE_HETATM;
        }


        // SERIAL FIELD
        serial     = StringFunctions::convertIntegerToString( currAtom->getSerialFromInputFile() );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_SERIAL_LENGTH, serial );
        recLine    = recLine + blankSpace + serial + singleSpace;
        

        // NAME FIELD
        if( !currAtom->getResidue()->isStandardResidue() ) {
            string atomSymbol = (string)V::GeometryTier::ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol;
            name = atomSymbol + "  ";

            //string atomSymbol = StringFunctions::strTrim((string)ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol);
            //string atomName   = StringFunctions::strTrim( currAtom->getAtomNameFromInputFile() );

            //string::size_type posOfSymbol = atomName.find(atomSymbol);
            //if (posOfSymbol == 0) {
            //    name += " ";
            //}
            //name += atomName;

            //int numSpaces = PDB_ATOM_RECORD_NAME_LENGTH - name.length();
            //for (int i = 0; i < numSpaces; ++i) {
            //    name += " ";
            //}
        }

        else {
            string extraBranch = currAtom->getpChemicalProperties()->getExtraBrangeDesignatorInString();
            string atomName    = StringFunctions::strTrim( (string)ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol );
            string remoteInd   = currAtom->getpChemicalProperties()->getRemoteIndicatorInString();
            string branchDes   = currAtom->getpChemicalProperties()->getBrangeDesignatorInString();

            if( atomName.length() == 2 )
                extraBranch = "";

            name       = extraBranch + atomName + remoteInd + branchDes;
        }

        recLine    = recLine + name;

        // ALTLOC FIELD --> CURRENTLY, ALTLOC IS NOT APPLIED.
        altLoc     = " ";
        recLine    = recLine + altLoc;
        
        
        // RESIDUE NAME
        if (currAtom->getResidue()->isStandardResidue()) {
            resName = (string)V::GeometryTier::RESIDUE_FEATURES[currAtom->getResidue()->getResidueCode()].threeCodeName;
            //resName = currResidue->getResidueName();
        }
        else {
            string residueName = currAtom->getResidue()->getResidueName();
            if (residueName.length() >= PDB_ATOM_RECORD_RESNAME_LENGTH) {
                resName = residueName.substr(0, 3);
            }
            else {
                resName = residueName;
                int numSpaces = PDB_ATOM_RECORD_RESNAME_LENGTH - residueName.length();
                for (int i = 0; i < numSpaces; ++i) {
                    resName += " ";
                }
            }
        }

        //string nameOfResidue = currAtom->getResidue()->getResidueName();
        //if (nameOfResidue.length() == PDB_ATOM_RECORD_RESNAME_LENGTH) {
        //    resName = nameOfResidue;
        //}
        //else {
        //    resName = (string)(RESIDUE_FEATURES[currAtom->getResidue()->getResidueCode()].threeCodeName);
        //}
        recLine = recLine + resName + singleSpace;


        // CHAIN ID
        rg_INT chainIDInInt  = currAtom->getResidue()->getChain()->getID() + 65;
        char   chainIDInChar =  chainIDInInt;
        chainID    = chainID + chainIDInChar;
        recLine    = recLine + chainID;


        // RESIDUE SEQUENCE
        resSeq     = StringFunctions::convertIntegerToString( currAtom->getResidue()->getSequenceNumber() );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_RESSEQ_LENGTH, resSeq );
        recLine    = recLine + blankSpace + resSeq;


        // ICODE --> CURRENTLY, ICODE IS NOT APPLIED.
        recLine    = recLine + "    ";


        // X COORDINATE
        xCoord     = StringFunctions::convertDoubleToString( currAtom->getpAtomBall()->getCenter().getX(), 3 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_X_LENGTH, xCoord );
        recLine    = recLine + blankSpace + xCoord;


        // Y COORDINATE
        yCoord     = StringFunctions::convertDoubleToString( currAtom->getpAtomBall()->getCenter().getY(), 3 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_Y_LENGTH, yCoord );
        recLine    = recLine + blankSpace + yCoord;

        
        // Y COORDINATE
        zCoord     = StringFunctions::convertDoubleToString( currAtom->getpAtomBall()->getCenter().getZ(), 3 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_Z_LENGTH, zCoord );
        recLine    = recLine + blankSpace + zCoord;


        // OCCUPANCY
        occupancy  = StringFunctions::convertDoubleToString( currAtom->getpChemicalProperties()->getOccupancy(), 2 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_OCCUPANCY_LENGTH, occupancy );
        recLine    = recLine + blankSpace + occupancy;


        // TEMPFACTOR
        tempFactor = StringFunctions::convertDoubleToString( currAtom->getpChemicalProperties()->getTempFactor(), 2 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_TEMPFACTOR_LENGTH, tempFactor );
        recLine    = recLine + blankSpace + tempFactor;

        // SEG_ID --> CURRENTLY, SEG_ID IS NOT APPLIED.
        recLine    = recLine + "          ";


        // ELEMENT
        element = ATOM_FEATURES[currAtom->getAtomCode()].symbol;
        recLine    = recLine + element;


        // CHARGE --> CURRENTLY, CHARGE IS NOT APPLIED.
        recLine    = recLine + "   ";


        fout << recLine << endl;
    }
}



void    MoleculeIOFunctions::writeConnectField( Molecule* aMolecule, ofstream& fout )
{
    rg_dList<Atom>* atoms = aMolecule->getAtoms();
    atoms->reset4Loop();

    Atom* currAtom = rg_NULL;
    while ( atoms->setNext4Loop() ) {
        currAtom = atoms->getpEntity();
        
        rg_INT numOfBondsForCurrAtom = currAtom->getListChemicalBond()->getSize();
       
        if( numOfBondsForCurrAtom == 0 )
            continue;

        rg_INT numOfBondsCounted = 0;        

        string recLine    = PDB_RECORD_TYPE_CONECT;
        string atomSerial = StringFunctions::convertIntegerToString( currAtom->getSerialFromInputFile() );
        string blankSpace = StringFunctions::getBlankStringForFixedField( PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH, atomSerial );

        recLine = recLine + blankSpace + atomSerial;

        rg_dList<ChemicalBond*>* bondListForCurAtom = currAtom->getListChemicalBond();

        bondListForCurAtom->reset4Loop();
        while ( bondListForCurAtom->setNext4Loop() ) {    

            ChemicalBond* currBond      = bondListForCurAtom->getEntity();
            Atom*         connectedAtom = currBond->getBondedAtom( currAtom );

            rg_FLAG isCurrAtomGeneral = isGeneralAtomInAminoOrDnaOrRnaResidue( currAtom );
            rg_FLAG isConnAtomGeneral = isGeneralAtomInAminoOrDnaOrRnaResidue( connectedAtom );

            if( isCurrAtomGeneral == rg_TRUE && isConnAtomGeneral == rg_TRUE ) {
                continue;                
            }
            
            else if( ( isCurrAtomGeneral == rg_TRUE  && isConnAtomGeneral == rg_FALSE ) ||
                     ( isCurrAtomGeneral == rg_FALSE && isConnAtomGeneral == rg_TRUE  )  ) {
                
                atomSerial = StringFunctions::convertIntegerToString( connectedAtom->getSerialFromInputFile() );
                blankSpace = StringFunctions::getBlankStringForFixedField( PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH, atomSerial );

                recLine    = recLine + blankSpace + atomSerial;

                numOfBondsCounted++;

                if( numOfBondsCounted == 4 ) {
                    fout << recLine << endl;
                    atomSerial = StringFunctions::convertIntegerToString( currAtom->getSerialFromInputFile() );
                    blankSpace = StringFunctions::getBlankStringForFixedField( PDB_CONNECT_RECORD_ATM_SERIAL_LENGTH, atomSerial );
                    recLine = PDB_RECORD_TYPE_CONECT + blankSpace + atomSerial;
                    numOfBondsCounted = 0;
                }

            }
        }

        if( numOfBondsCounted > 0 )
            fout << recLine << endl;
    }
}



rg_FLAG MoleculeIOFunctions::isGeneralAtomInAminoOrDnaOrRnaResidue( Atom* anAtom )
{
    rg_FLAG isGeneralAtom = rg_FALSE;
    if( anAtom->getResidue()->isAminoResidue() == rg_TRUE ) {
        if( anAtom->getAtomCode() == C_ATOM || anAtom->getAtomCode() == N_ATOM ||
            anAtom->getAtomCode() == O_ATOM || anAtom->getAtomCode() == S_ATOM  ) {
            isGeneralAtom =  rg_TRUE;
        }
    }
    else if ( anAtom->getResidue()->isDNAResidue() == rg_TRUE || anAtom->getResidue()->isRNAResidue() == rg_TRUE ) {
        if( anAtom->getAtomCode() == C_ATOM || anAtom->getAtomCode() == N_ATOM ||
            anAtom->getAtomCode() == O_ATOM || anAtom->getAtomCode() == S_ATOM ||  anAtom->getAtomCode() == P_ATOM ) {
            isGeneralAtom =  rg_TRUE;
        }
    }

    return isGeneralAtom;
}



void    MoleculeIOFunctions::writeHeteroAtomField( Molecule* aMolecule, ofstream& fout )
{
    
}



rg_FLAG MoleculeIOFunctions::writerTriposMol2File( const char* pathFile, Molecule* aMolecule )
{
    rg_FLAG isFileWritten   = rg_TRUE;
    

    ofstream fout(pathFile);

    fout.setf(ios::fixed, ios::floatfield);     // TO 4 PLACES OF DECIMALS 
    fout.precision(4);                          //

    writeTriposMoleculeFile(aMolecule, fout);
    writeTriposAtomField(aMolecule, fout);
    writeTriposBondField(aMolecule, fout);

    fout.close();

    return isFileWritten;    
}


void MoleculeIOFunctions::writeTriposMoleculeFile( Molecule* aMolecule, ofstream& fout )
{
    fout << "@<TRIPOS>MOLECULE" << endl;
    
    fout << aMolecule->getMoleculeName() << endl;

    fout << aMolecule->getAtoms()->getSize() << "\t\t" << aMolecule->getChemicalBonds()->getSize() << "\t\t0\t\t0\t\t0" << endl;

    fout << "SMALL" << endl;

    fout << "USER_CHARGES" << endl << endl;
}

void MoleculeIOFunctions::writeTriposAtomField( Molecule* aMolecule, ofstream& fout )
{
    fout << "@<TRIPOS>ATOM" << endl;

    rg_dList<Atom>* atoms = aMolecule->getAtoms();
    
    atoms->reset4Loop();    
    while(atoms->setNext4Loop()) {
        Atom* currAtom = atoms->getpEntity();                
                    
        fout << currAtom->getSerialFromInputFile() << "\t\t"
             << currAtom->getAtomNameFromInputFile() << "\t\t";

        if(currAtom->getAtomBall().getCenter().getX() > 0) {
            fout << " ";
        }

        fout << currAtom->getAtomBall().getCenter().getX() << "\t\t";


        if(currAtom->getAtomBall().getCenter().getY() > 0) {
            fout << " ";
        }

        fout << currAtom->getAtomBall().getCenter().getY() << "\t\t";


        if(currAtom->getAtomBall().getCenter().getZ() > 0) {
            fout << " ";
        }
            
        fout << currAtom->getAtomBall().getCenter().getZ() << "\t\t";

        
        fout << currAtom->getpChemicalProperties()->getSYBYLAtomType() << "\t\t"
             << currAtom->getResidue()->getSequenceNumber() << "\t\t"
             << currAtom->getResidue()->getResidueName() << "\t\t";


        if(currAtom->getpChemicalProperties()->getCharge() > 0) {
            fout << " ";
        }

        fout << currAtom->getpChemicalProperties()->getCharge() << endl;
    }

}

void MoleculeIOFunctions::writeTriposBondField( Molecule* aMolecule, ofstream& fout )
{
    fout << "@<TRIPOS>BOND" << endl;


    rg_INT bondIndex = 1;
    rg_dList<ChemicalBond>* bonds = aMolecule->getChemicalBonds();

    bonds->reset4Loop();
    while (bonds->setNext4Loop()) {
        ChemicalBond* currBond = bonds->getpEntity();

        Atom* atom[2] = { currBond->getFirstAtom(), currBond->getSecondAtom() };
        if (atom[0]->getSerialFromInputFile() > atom[1]->getSerialFromInputFile()) {
            atom[0] = currBond->getSecondAtom();
            atom[1] = currBond->getFirstAtom();
        }

        fout << bondIndex << "\t" << atom[0]->getSerialFromInputFile() << "\t" << atom[1]->getSerialFromInputFile() << "\t";

        switch (currBond->getTypeOfBond())
        {
        case SINGLE_BOND:
            fout << "1" << endl;
            break;
        case DOUBLE_BOND:
            fout << "2" << endl;
            break;
        case TRIPLE_BOND:
            fout << "3" << endl;
            break;
        case AMIDE_BOND:
            fout << "am" << endl;
            break;
        case AROMATIC_BOND:
            fout << "ar" << endl;
            break;
        default:
            break;
        }

        ++bondIndex;
    }


    /*
    rg_INT bondIndex = 1;

    fout << "@<TRIPOS>BOND" << endl;

    rg_dList<Atom>* atoms = aMolecule->getAtoms();

    atoms->reset4Loop();

    while(atoms->setNext4Loop()) {
        Atom* currAtom = atoms->getpEntity();        

        rg_dList<ChemicalBond*>* bonds = currAtom->getListChemicalBond();

        bonds->reset4Loop();
        while(bonds->setNext4Loop()) {
            ChemicalBond* currBond = bonds->getEntity();

            Atom* oppositeAtom = rg_NULL;

            if(currBond->getFirstAtom() == currAtom) {
                oppositeAtom = currBond->getSecondAtom();
            }
            else if(currBond->getSecondAtom() == currAtom) {
                oppositeAtom = currBond->getFirstAtom();
            }

            if( oppositeAtom != rg_NULL && 
                oppositeAtom->getSerialFromInputFile() > currAtom->getSerialFromInputFile() ) {
                
                fout << bondIndex++ << "\t" << currAtom->getSerialFromInputFile() << "\t" << oppositeAtom->getSerialFromInputFile() << "\t";

                switch (currBond->getTypeOfBond())
                {
                    case SINGLE_BOND:
                        fout << "1" << endl;
                	    break;
                    case DOUBLE_BOND:
                        fout << "2" << endl;
                	    break;
                    case TRIPLE_BOND:
                        fout << "3" << endl;
                	    break;
                    case AMIDE_BOND:
                        fout << "am" << endl;
                	    break;
                    case AROMATIC_BOND:
                        fout << "ar" << endl;
                	    break;
                    default :
                        break;
                }
            }
        }
    }
    */
}


rg_FLAG MoleculeIOFunctions::readCTMolFile( const rg_INT& targetMoleculeID, const char* pathFile, Molecule& aMolecule )
{
    rg_FLAG isFileRead   = rg_FALSE;

    if( !isValidMoleculeFileExtension( pathFile, MOL_FILE ) ) {
        cerr << "Error: Invalid file extension !\n";
        return isFileRead;
    }
    
    ifstream fin( pathFile );
    
    if( fin.bad() ) {
        cerr << "Error: Could not open file!\n";
        return isFileRead;
    }

    string strRecLine  = "";
    if( !getline( fin, strRecLine ) )
        return isFileRead;
    
    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    
    list<string> recordsOfMolFile;
    addRecordLinesOfFIleToList( &fin, recordsOfMolFile );
    fin.close();

    AtomMap          mapOfAtom;   
    AtomSymbolMap    mapOfAtomSymbol;
    
    initiallizeSymbolMapOfAtomForMolFile( mapOfAtomSymbol );


    rg_INT countLines = 0;
    rg_INT numOfRecLines = recordsOfMolFile.size();
    
    list<string>::iterator i_recLines = recordsOfMolFile.begin();

    readRecordsInHeaderBlock( i_recLines , aMolecule );
    
    rg_FLAG isRecLinesOK = readRecordsInCTAB( i_recLines, mapOfAtomSymbol, mapOfAtom, aMolecule );

    if( isRecLinesOK ) {
     
        setAmberAtomTypeForAtomsInMolecule( aMolecule );       
        
        string moleculeFileName = StringFunctions::getFileNameWithoutPath( string(pathFile) );
        aMolecule.setMoleculeFileName( moleculeFileName );
        
        isFileRead = rg_TRUE;
    }
    
    return isFileRead;
}



void MoleculeIOFunctions::addRecordLinesOfFIleToList( ifstream* fin, list<string>& recordLinesOfFile )
{
    string strRecLine  = "";
    
    while ( getline( *fin, strRecLine ) ) {
        recordLinesOfFile.push_back( StringFunctions::strTrimRight( strRecLine ) );
    }
}



void MoleculeIOFunctions::readRecordsInHeaderBlock( list<string>::iterator& i_recLines, Molecule& aMolecule )
{
    rg_INT countLine = 0;
    while ( countLine < MOL_NUM_OF_LINES_IN_HEADER_BLOCK ) {
        countLine++;

        aMolecule.addHeaderRecords( StringFunctions::strTrimLeft( string(*i_recLines) ) );
        i_recLines++;
    }
}



rg_FLAG MoleculeIOFunctions::readRecordsInCTAB( list<string>::iterator& i_recLines, AtomSymbolMap& mapOfAtomSymbol, AtomMap& mapOfAtom, Molecule& aMolecule )
{   
    rg_FLAG isRecLineOK = rg_TRUE;

    Chain*   newChain   = aMolecule.addChain( Chain(0, &aMolecule) );
    Residue* newResidue = aMolecule.addResidue( Residue(0) );
    newChain->addResidue( newResidue );
    newResidue->setChain( newChain );
    newResidue->setResidueName( string("UNK") );
    
    // READ COUNTS LINE
    
    string strCountsLine = *i_recLines;
    
    rg_INT numOfAtoms = atoi( StringFunctions::subString( strCountsLine, MOL_COUNTS_LINE_ST_POS_NUM_ATOMS, MOL_COUNTS_LINE_LENGTH_NUM_ATOMS ).c_str() );
    rg_INT numOfBonds = atoi( StringFunctions::subString( strCountsLine, MOL_COUNTS_LINE_ST_POS_NUM_BONDS, MOL_COUNTS_LINE_LENGTH_NUM_BONDS ).c_str() );
    
    // READ ATOM BLOCK
    string strAtomLine = "";
    rg_INT i_lineNum =0;
    for( i_lineNum=1; i_lineNum<=numOfAtoms; i_lineNum++ ) {
        i_recLines++;
        strAtomLine = *i_recLines;
        
        rg_REAL  xCoord = atof( StringFunctions::subString( strAtomLine, MOL_ATOM_BLOCK_ST_POS_X_COORD, MOL_ATOM_BLOCK_LENGTH_X_COORD ).c_str() );
        rg_REAL  yCoord = atof( StringFunctions::subString( strAtomLine, MOL_ATOM_BLOCK_ST_POS_Y_COORD, MOL_ATOM_BLOCK_LENGTH_Y_COORD ).c_str() );
        rg_REAL  zCoord = atof( StringFunctions::subString( strAtomLine, MOL_ATOM_BLOCK_ST_POS_Z_COORD, MOL_ATOM_BLOCK_LENGTH_Z_COORD ).c_str() );
        string   atomSymbol = StringFunctions::strTrim( StringFunctions::subString( strAtomLine, MOL_ATOM_BLOCK_ST_POS_ATOM_SYMBOL, MOL_ATOM_BLOCK_LENGTH_ATOM_SYMBOL ) );

       
        AtomSymbolMap::iterator AtomSymbolMap_i = mapOfAtomSymbol.find( atomSymbol );
        
        AtomCode atomCode = UNK_ATOM;
        if ( AtomSymbolMap_i != mapOfAtomSymbol.end() ) {
            atomCode = (AtomCode)((*AtomSymbolMap_i).second);
        }        
        
        Atom* newAtom = aMolecule.addAtom( Atom(i_lineNum-1) );
        newAtom->setSerialFromInputFile( i_lineNum );
        newAtom->setAtomCode( atomCode );
        newAtom->setAtomBall( Sphere( rg_Point3D(xCoord, yCoord, zCoord), ATOM_FEATURES[atomCode].radius) );
        
        newAtom->setResidue( newResidue );
        newResidue->addAtom( newAtom );
        mapOfAtom.insert( AtomMap::value_type(i_lineNum, newAtom) );
    }
    
    // READ BOND BLOCK
    string strBondLine = "";
    for( i_lineNum=1; i_lineNum<=numOfBonds; i_lineNum++ ) {
        i_recLines++;
        strBondLine = *i_recLines;
        
        rg_INT firstAtomID  = atoi( StringFunctions::subString( strBondLine, MOL_BOND_BLOCK_ST_POS_FIRST_ATOM,  MOL_BOND_BLOCK_LENGTH_FIRST_ATOM ).c_str() );
        rg_INT secondAtomID = atoi( StringFunctions::subString( strBondLine, MOL_BOND_BLOCK_ST_POS_SECOND_ATOM, MOL_BOND_BLOCK_LENGTH_SECOND_ATOM ).c_str() );
        rg_INT intTypeBond  = atoi( StringFunctions::subString( strBondLine, MOL_BOND_BLOCK_ST_POS_BOND_TYPE,   MOL_BOND_BLOCK_LENGTH_BOND_TYPE ).c_str() );
        
        BondType typeOfBond = convertIntegerToBondType( intTypeBond );
        
        Atom* firstAtomForChemicalBond  = rg_NULL;
        Atom* secondAtomForChemicalBond = rg_NULL;

        AtomMap::iterator AtomMap_i = mapOfAtom.find( firstAtomID );
        AtomMap::iterator AtomMap_j = mapOfAtom.find( secondAtomID );
        
        if ( AtomMap_i != mapOfAtom.end() && AtomMap_j != mapOfAtom.end() ) {
            firstAtomForChemicalBond  = (*AtomMap_i).second;
            secondAtomForChemicalBond = (*AtomMap_j).second;

            setChemicalBondForMolFileToMoleculeWithoutDuplication( firstAtomForChemicalBond, secondAtomForChemicalBond, i_lineNum, typeOfBond, aMolecule );
        }
        else {
            isRecLineOK = rg_FALSE;
            break;
        }
    }
    
    if( isRecLineOK ) {
        aMolecule.evaluateHydrogenDonorAndAcceptorAtoms();
        aMolecule.computeAndSetCenterOfMass();
        aMolecule.computeAndSetMinEnclosingSphere();
    }

    return isRecLineOK;
}



BondType MoleculeIOFunctions::convertIntegerToBondType( const rg_INT& intBondType )
{
    BondType typeOfBond = UNK_BOND;
    
    switch ( intBondType )  {
    case 1:
        typeOfBond = SINGLE_BOND;
        break;
    case 2:
        typeOfBond = DOUBLE_BOND;
        break;
    case 3:
        typeOfBond = TRIPLE_BOND;
        break;
    case 4:
        typeOfBond = AROMATIC_BOND;
        break;
    case 5:
        typeOfBond = SINGLE_OR_DOUBLE_BOND;
        break;
    case 6:
        typeOfBond = SINGLE_OR_AROMATIC_BOND;
        break;
    case 7:
        typeOfBond = DOUBLE_OR_AROMATIC_BOND;
        break;
    default:
        typeOfBond = UNK_BOND;
        break;
    }
    return typeOfBond;
}



void MoleculeIOFunctions::setAmberAtomTypeForAtomsInMolecule( Molecule& aMolecule )
{
    rg_dList<Atom>* atoms = aMolecule.getAtoms();

    atoms->reset4Loop();
    while( atoms->setNext4Loop() ) {
        Atom* currAtom = atoms->getpEntity();
        AmberAtomTypes aAmberAtomType = UNK_ATM;
        
        switch ( currAtom->getAtomCode() )  {
            case C_ATOM:
                aAmberAtomType = getAmberAtomTypeForCarbon( currAtom );
                break;
            case N_ATOM:
                aAmberAtomType = getAmberAtomTypeForNitrogen( currAtom );
                break;
            case O_ATOM:
                aAmberAtomType = getAmberAtomTypeForOxygen( currAtom );
                break;
            case H_ATOM:
                aAmberAtomType = getAmberAtomTypeForHydrogen( currAtom );
                break;
            case F_ATOM:
                aAmberAtomType = F_ATM;
                break;
            case K_ATOM:
                aAmberAtomType = K_ATM;
                break;
            case LI_ATOM:
                aAmberAtomType = Li_ATM;
                break;
            case P_ATOM:
                aAmberAtomType = P_ATM;
                break;
            case RB_ATOM:
                aAmberAtomType = Rb_ATM;
                break;
            case S_ATOM:
                aAmberAtomType = S_ATM;
                break;
            default:
                aAmberAtomType = UNK_ATM;
                break;
        }
        currAtom->getpChemicalProperties()->setAtomTypeInAmber( aAmberAtomType );
    }
}



AmberAtomTypes MoleculeIOFunctions::getAmberAtomTypeForCarbon( Atom* carbon )
{
    rg_dList<ChemicalBond*>* bondList = carbon->getListChemicalBond();

    rg_INT numOfAllBonds      = bondList->getSize();
    rg_INT numOfSingleBonds   = 0;
    rg_INT numOfAromaticBonds = 0;
    rg_INT numOfDoubleBonds   = 0;
    rg_INT numOfOtherBonds    = 0;

    bondList->reset4Loop();
    while( bondList->setNext4Loop() ) {
        ChemicalBond* currBond = bondList->getEntity();
        switch ( currBond->getTypeOfBond() )  {
            case SINGLE_BOND:
                numOfSingleBonds++;
                break;
            case DOUBLE_BOND:
                numOfDoubleBonds++;
                break;
            case AROMATIC_BOND:
                numOfAromaticBonds++;
                break;
            default:
                numOfOtherBonds++;
                break;
        }
    }

    AmberAtomTypes aAmberAtomType = UNK_ATM;

    if( numOfAllBonds == numOfSingleBonds ) {
        aAmberAtomType = CT_ATM;
    }
    else {
        if( numOfAromaticBonds > 0 ) {
            aAmberAtomType = CA_ATM;
        }
        else {
            if( numOfDoubleBonds > 0 ) {
                aAmberAtomType = CM_ATM;
            }
            else {
                aAmberAtomType = C_ATM;
            }
        }
    }

    return aAmberAtomType;
}



AmberAtomTypes MoleculeIOFunctions::getAmberAtomTypeForNitrogen( Atom* nitrogen )
{
    rg_dList<ChemicalBond*>* bondList = nitrogen->getListChemicalBond();

    rg_INT numOfAllBonds      = bondList->getSize();

    AmberAtomTypes aAmberAtomType = UNK_ATM;

    if( numOfAllBonds == 4 ) {
        aAmberAtomType = N3_ATM;
    }
    else {
        aAmberAtomType = N_ATM;
    }

    return aAmberAtomType;
}



AmberAtomTypes MoleculeIOFunctions::getAmberAtomTypeForOxygen( Atom* oxygen )
{
    rg_dList<ChemicalBond*>* bondList = oxygen->getListChemicalBond();

    rg_INT numOfAllBonds      = bondList->getSize();

    AmberAtomTypes aAmberAtomType = UNK_ATM;

    if( numOfAllBonds == 1 ) {
        aAmberAtomType = O2_ATM;
    }
    else{
        rg_INT numOfSurfurConnected = 0;
        rg_INT numOfDoubleBonds     = 0;

        bondList->reset4Loop();
        while( bondList->setNext4Loop() ) {
            ChemicalBond* currBond = bondList->getEntity();
            
            if( currBond->getFirstAtom()->getAtomCode()  == S_ATOM || 
                currBond->getSecondAtom()->getAtomCode() == S_ATOM ) {
                numOfSurfurConnected++;
            }
            if( currBond->getTypeOfBond() == DOUBLE_BOND ) {
                numOfDoubleBonds++;
            }
        }

        if( numOfSurfurConnected > 0 ) {
            aAmberAtomType = OS_ATM;
        }
        else {
            if( numOfDoubleBonds > 0 ) {
                aAmberAtomType = O_ATM;
            }
            else {
                aAmberAtomType = OH_ATM;
            }
        }
    }

    return aAmberAtomType;
}



AmberAtomTypes MoleculeIOFunctions::getAmberAtomTypeForHydrogen( Atom* hydrogen )
{
    rg_dList<ChemicalBond*>* bondList = hydrogen->getListChemicalBond();

    rg_INT numOfAllBonds = bondList->getSize();

    AmberAtomTypes aAmberAtomType = UNK_ATM;

    rg_INT numOfNitrogenConnected = 0;
    rg_INT numOfOxygenConnected   = 0;
    rg_INT numOfSurfurConnected   = 0;
    rg_INT numOfCACarbonConnected   = 0;

    bondList->reset4Loop();
    while( bondList->setNext4Loop() ) {
        ChemicalBond* currBond = bondList->getEntity();   
       
        Atom* connectedAtom = rg_NULL;

        if( currBond->getFirstAtom()->getAtomCode() == H_ATOM ) {
            connectedAtom = currBond->getSecondAtom();
        }
        else {
            connectedAtom = currBond->getFirstAtom();
        }

        AtomCode atomCodeOfConnectedAtom = connectedAtom->getAtomCode();

        if( atomCodeOfConnectedAtom == N_ATOM ) {
            numOfNitrogenConnected++;
        }
        else if( atomCodeOfConnectedAtom == O_ATOM ) {
            numOfOxygenConnected++;
        }
        else if( atomCodeOfConnectedAtom == S_ATOM ) {
            numOfSurfurConnected++;
        }
        else if( atomCodeOfConnectedAtom == C_ATOM ) {
            if( getAmberAtomTypeForCarbon( connectedAtom ) == CA_ATM ) {
                numOfCACarbonConnected++;
            }
        }
        else {}
    }

    if( numOfNitrogenConnected > 0 ) {
        aAmberAtomType = H_ATM;
    }
    else {
        if( numOfOxygenConnected > 0 ) {
            aAmberAtomType = HO_ATM;
        }
        else {
            if( numOfSurfurConnected > 0 ) {
                aAmberAtomType = HS_ATM;
            }
            else {
                if( numOfCACarbonConnected > 0 ) {
                    aAmberAtomType = HA_ATM;
                }
                else {
                    aAmberAtomType = HC_ATM;
                }
            }
        }
    }

    return aAmberAtomType;
}



rg_FLAG MoleculeIOFunctions::writerCTMolFile( const char* pathFile, Molecule* aMolecule )
{
    rg_FLAG isFileWritten   = rg_FALSE;
    


    return isFileWritten;
    
}



void MoleculeIOFunctions::mergeTwoMolecules( Molecule* molecule_A, Molecule* molecule_B, Molecule* mergedMolecule )
{
    Molecule* twoMolecule[2];
    twoMolecule[0] = molecule_A;
    twoMolecule[1] = molecule_B;
    

    Chain*    newChains[2];
    newChains[0] = mergedMolecule->addChain( Chain(0, mergedMolecule) );
    newChains[1] = mergedMolecule->addChain( Chain(1, mergedMolecule) );

    newChains[0]->setChainIDFromInputFileInDecimal( 65 );
    newChains[1]->setChainIDFromInputFileInDecimal( 66 );

    newChains[0]->setChainCode( twoMolecule[0]->getChains()->getFirstpEntity()->getChainCode() );
    newChains[1]->setChainCode( twoMolecule[1]->getChains()->getFirstpEntity()->getChainCode() );



    rg_INT numOfResidue = 0;
    rg_INT numOfAtom    = 0;


    for( rg_INT i_chain=0; i_chain<2; i_chain++ ) {

        rg_dList<Residue>* residues = twoMolecule[i_chain]->getResidues();

        Residue* currResidue = rg_NULL;
        Residue* newResidue  = rg_NULL;

        residues->reset4Loop();
        while( residues->setNext4Loop() ) {
            currResidue = residues->getpEntity();
        
            newResidue = mergedMolecule->addResidue( Residue(numOfResidue) );
            newChains[i_chain]->addResidue( newResidue );
            newResidue->setChain( newChains[i_chain] );
            newResidue->setSequenceNumber( currResidue->getSequenceNumber() );
            newResidue->setResidueCode( currResidue->getResidueCode() );
            newResidue->setResidueName( currResidue->getResidueName() );

            rg_dList<Atom*>* atoms = currResidue->getAtoms();

            Atom* currAtom = rg_NULL;
            Atom* newAtom  = rg_NULL;

            atoms->reset4Loop();
            while( atoms->setNext4Loop() ) {
                currAtom = atoms->getEntity();
                
                newAtom = mergedMolecule->addAtom( Atom(numOfAtom) );
                newResidue->addAtom( newAtom );
                newAtom->setResidue( newResidue );
                newAtom->setSerialFromInputFile( currAtom->getSerialFromInputFile() );
                newAtom->setAtomNameFromInputFile( currAtom->getAtomNameFromInputFile() );
                newAtom->setAtomBall( currAtom->getAtomBall() );
                newAtom->setAtomCode( currAtom->getAtomCode() );
                newAtom->setChemicalProperties( currAtom->getChemicalProperties() );
                numOfAtom++;
            }
            numOfResidue++;
        }
    }
}



void MoleculeIOFunctions::initiallizeSymbolMapsOfAtomAndResidueForPDBFile( AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbolForProtein, ResidueSymbolMap& mapOfResidueSymbolForDNA, ResidueSymbolMap& mapOfResidueSymbolForRNA )
{
    mapOfAtomSymbol.clear();
    mapOfResidueSymbolForProtein.clear();
    mapOfResidueSymbolForDNA.clear();
    mapOfResidueSymbolForRNA.clear();
    
    for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        mapOfAtomSymbol.insert( AtomSymbolMap::value_type(atomSymbol, i_atomCode) );
    }
    
    rg_INT i_residueCode = 0;
    for ( i_residueCode = 0; i_residueCode<NUM_OF_AMINOACIDS+1; i_residueCode++ ) {
        string residueSymbol( RESIDUE_FEATURES[i_residueCode].threeCodeName );
        mapOfResidueSymbolForProtein.insert( ResidueSymbolMap::value_type(residueSymbol, i_residueCode) );
    }
    mapOfResidueSymbolForProtein.insert( ResidueSymbolMap::value_type( string( RESIDUE_FEATURES[HOH_RESIDUE].threeCodeName ), HOH_RESIDUE) );

    for ( i_residueCode = A_DNA_RESIDUE; i_residueCode<A_DNA_RESIDUE+NUM_OF_DNAS; i_residueCode++ ) {
        string residueSymbol( RESIDUE_FEATURES[i_residueCode].threeCodeName );
        mapOfResidueSymbolForDNA.insert( ResidueSymbolMap::value_type(residueSymbol, i_residueCode) );
    }

    for ( i_residueCode = A_RNA_RESIDUE; i_residueCode<A_RNA_RESIDUE+NUM_OF_RNAS; i_residueCode++ ) {
        string residueSymbol( RESIDUE_FEATURES[i_residueCode].threeCodeName );
        mapOfResidueSymbolForRNA.insert( ResidueSymbolMap::value_type(residueSymbol, i_residueCode) );
    }
}



void MoleculeIOFunctions::initiallizeSymbolMapsOfAtomAndResidueForPDBFile( AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol )
{
    mapOfAtomSymbol.clear();
    mapOfResidueSymbol.clear();
    
    for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        mapOfAtomSymbol.insert( AtomSymbolMap::value_type(atomSymbol, i_atomCode) );
    }
    
    for ( rg_INT i_residueCode = 0; i_residueCode<NUM_OF_RESIDUES; i_residueCode++ ) {
        string residueSymbol( RESIDUE_FEATURES[i_residueCode].threeCodeName );
        mapOfResidueSymbol.insert( ResidueSymbolMap::value_type(residueSymbol, i_residueCode) );
    }
}



void MoleculeIOFunctions::initiallizeSymbolMapsOfAtomAndResidueForMol2File( AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol )
{
    mapOfAtomSymbol.clear();
    mapOfResidueSymbol.clear();
    
    for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        mapOfAtomSymbol.insert( AtomSymbolMap::value_type( StringFunctions::strTrim( atomSymbol ), i_atomCode ) );
    }
    
    for ( rg_INT i_residueCode = 0; i_residueCode<NUM_OF_RESIDUES; i_residueCode++ ) {
        string residueSymbol( RESIDUE_FEATURES[i_residueCode].threeCodeName );
        mapOfResidueSymbol.insert( ResidueSymbolMap::value_type(residueSymbol, i_residueCode) );
    }   
}



void MoleculeIOFunctions::initiallizeSymbolMapOfAtomForMolFile( AtomSymbolMap& mapOfAtomSymbol )
{
    mapOfAtomSymbol.clear();
    
    for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        mapOfAtomSymbol.insert( AtomSymbolMap::value_type( StringFunctions::strTrim( atomSymbol ), i_atomCode ) );
    }
}



rg_FLAG MoleculeIOFunctions::isValidMoleculeFileExtension( const char* pathFile, const MoleculeFileType& aMoleculeFileType )
{
    rg_FLAG isValidExt = rg_FALSE;

    string nameOfExt = StringFunctions::getExtensionFromFileName( string(pathFile) );
    nameOfExt        = StringFunctions::strToLower( nameOfExt );

    switch( aMoleculeFileType ) {
        case PDB_FILE :
            if( nameOfExt == PDB_FILE_EXTENSION_1 || nameOfExt == PDB_FILE_EXTENSION_2 )
                isValidExt = rg_TRUE;
            break;
        case MOL2_FILE :
            if( nameOfExt == MOL2_FILE_EXTENSION )
                isValidExt = rg_TRUE;
            break;
        case MOL_FILE :
            if( nameOfExt == MOL_FILE_EXTENSION )
                isValidExt = rg_TRUE;
            break;
        case RMC_CFG_FILE:
            if (nameOfExt == RMC_CFG_FILE_EXTENSION)
                isValidExt = rg_TRUE;
            break;
        case PDBQ_FILE:
            if (nameOfExt == PDBQ_FILE_EXTENSION)
                isValidExt = rg_TRUE;
            break;
        default:
            isValidExt = rg_FALSE;
            break;
    }

    return isValidExt;
}



rg_BOOL MoleculeIOFunctions::checkAndReportErrorCodeForRecordType( const string& aRecType, const rg_FLAG& isRecLineOK, rg_FLAG& isFileRead )
{
    rg_BOOL isRecLineAcceptable = rg_TRUE;
    isFileRead = rg_TRUE;

    if( isRecLineOK == rg_FALSE ) {
        cerr << "Error: reading " << aRecType << "\n";
        isFileRead          = rg_FALSE;
        isRecLineAcceptable = rg_FALSE;
    }

    return isRecLineAcceptable;
}



rg_FLAG MoleculeIOFunctions::readChemOfficeChargeFile( const char* pathFile, Molecule& aMolecule )
{
    rg_FLAG isFileRead = rg_FALSE;


    if( aMolecule.getAtoms()->getSize() == 0 ) {
        cerr << "Error: Molecule is empty!\n";
        return isFileRead;
    }

    ifstream fin( pathFile );
    
    if( fin.bad() ) {
        cerr << "Error: Could not open file!\n";
        return isFileRead;
    }
    
    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    

    list<string> recordsOfChargeFile;
    addRecordLinesOfFIleToList( &fin, recordsOfChargeFile );
    fin.close();


    AtomMap mapOfAtom;   
    rg_dList<Atom>* listAtom = aMolecule.getAtoms();

    listAtom->reset4Loop();
    while( listAtom->setNext4Loop() ) {
        Atom* currAtom = listAtom->getpEntity();
        mapOfAtom.insert( AtomMap::value_type(currAtom->getSerialFromInputFile(), currAtom) );
    }

       
    string  atomType;
    rg_REAL chargeValue;
    rg_INT  atomSerial;

    list<string>::iterator i_recLines = recordsOfChargeFile.begin();
    
    rg_FLAG isRecLineOK = rg_TRUE;

    while( i_recLines != recordsOfChargeFile.end() ) {
        string currRecLine = *i_recLines;

        extractChargeRecordsFromRecLine( currRecLine, atomType, chargeValue, atomSerial );
    
        Atom* targetAtom = getpAtomFromMap( &mapOfAtom, atomSerial );

        if( targetAtom != rg_NULL ) {
            targetAtom->getpChemicalProperties()->setCharge( chargeValue );
        }
        else {
            isRecLineOK = rg_FALSE;
            break;
        }

        i_recLines++;
    }

    if( isRecLineOK == rg_TRUE )
        isFileRead = rg_TRUE;
    
    return isFileRead;
}



void MoleculeIOFunctions::extractChargeRecordsFromRecLine( const string& strRecLine, string& atomType, rg_REAL& chargeValue, rg_INT& atomSerial )
{
    rg_INT stPosOfSingleRec = 0;
    rg_INT edPosOfSingleRec = 0;
    
    // ATOM TYPE
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    atomType = StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec );
    
    
    // CHARGE VALUE
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    chargeValue = atof( StringFunctions::subString( strRecLine, stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str() );
    
    
    // ATOM SERIAL
    stPosOfSingleRec = edPosOfSingleRec;
    StringFunctions::updateStartPositionWithoutWhiteSpace( strRecLine, stPosOfSingleRec );
    edPosOfSingleRec = strRecLine.find_first_of( STR_SEPS, stPosOfSingleRec );
    
    string strAtomSerial = StringFunctions::subString( strRecLine,  stPosOfSingleRec, edPosOfSingleRec-stPosOfSingleRec ).c_str();

    rg_INT lengthOfAtomSerial = strAtomSerial.length();
    rg_INT stPosOfSerial = strAtomSerial.find("(");
    rg_INT edPosOfSerial = strAtomSerial.find(")");
    
    atomSerial = atoi( StringFunctions::subString( strRecLine, stPosOfSerial+1, lengthOfAtomSerial-edPosOfSerial+1 ).c_str());
}



Atom* MoleculeIOFunctions::getpAtomFromMap( AtomMap* mapOfAtom, const rg_INT& atomID )
{
    AtomMap::iterator map_i = mapOfAtom->find( atomID );
    
    if ( map_i == mapOfAtom->end() )
        return rg_NULL;
    else
        return (*map_i).second;
}




// rg_BOOL MoleculeIOFunctions::readMol2File(const char* pathFile, rg_dList<Molecule>& molecules)
// {
//     rg_FLAG isFileRead   = rg_FALSE;
// 
//     if( !isValidMoleculeFileExtension( pathFile, MOL2_FILE ) ) {
//         cerr << "Error: Invalid file extension !\n";
//         return isFileRead;
//     }
//     
//     
//     ifstream fin( pathFile );
//     
//     if( fin.bad() ) {
//         cerr << "Error: Could not open file!\n";
//         return isFileRead;
//     }
// 
//     string currRecType = "";
//     string strRecLine  = "";
// 
//     if( !getline( fin, strRecLine ) )
//         return isFileRead;
// 
//     fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
// 
//     list<string> recordsOfMol2File;
//     addRecordsOfMol2FIleToList( &fin, recordsOfMol2File );
//     fin.close();
// 
//     rg_INT countLines = 0;
//     rg_INT numOfRecLines = recordsOfMol2File.size();
//     rg_FLAG isRecLinesOK = rg_TRUE;
// 
//     list<string>::iterator i_recLines = recordsOfMol2File.begin();
// 
//     while( i_recLines != recordsOfMol2File.end() ) {
//         
//         countLines++;
//         string* currRecLine = &(*i_recLines);
// 
//         string recType = StringFunctions::strTrim( *currRecLine );
// 
// 
//         if ( recType == MOL2_RECORD_TYPE_MOLECULE ) {
//             isRecLinesOK = setMoleculeRTIToMolecule( i_recLines, arrNumOfRecords, aMolecule );
//         }
//         else if ( recType == MOL2_RECORD_TYPE_ATOM && arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID] > 0 ) {
// //            CAN'T READ MOL2 FILE... NEED TO DEBUG....
//             isRecLinesOK = setAtomRTIToMolecule( i_recLines, &recordsOfMol2File, arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID], *mapOfAtom, *mapOfResidue, *mapOfNonStdResidue, *mapOfAtomSymbol, *mapOfResidueSymbol, aMolecule );
// //            isRecLinesOK = setAtomRTIToMolecule( i_recLines, &recordsOfMol2File, arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID], *mapOfAtom, *mapOfResidue, *mapOfAtomSymbol, *mapOfResidueSymbol, aMolecule );
//         }
//         else if ( recType == MOL2_RECORD_TYPE_BOND && arrNumOfRecords[MOL2_NUM_OF_BONDS_ID] >0 ) {
//             isRecLinesOK = setBondRTIToMolecule( i_recLines, &recordsOfMol2File, arrNumOfRecords[MOL2_NUM_OF_BONDS_ID], *mapOfAtom, aMolecule );
//         }
//         else if ( recType == MOL2_RECORD_TYPE_SUBSTRUCTURE ) {
//             int aaa = 0;
//             
//         }
// 
//         if( isRecLinesOK == rg_FALSE ) {
//             isFileRead = rg_FALSE;
//             break;
//         }
// 
//         i_recLines++;
//     }
// }


rg_FLAG MoleculeIOFunctions::writeRotatedPDBFileWithGivenAngle( const char* pathFile, Molecule* aMolecule, const rg_REAL& angleToRotate, const rg_INT& axis )
{
    rg_FLAG isFileWritten   = rg_FALSE;
    
    ofstream  fout( pathFile );
    
    writeHeader( aMolecule, fout );
    writeAtomField( aMolecule, fout, angleToRotate, axis );
    writeConnectField( aMolecule, fout );
    writeHeteroAtomField( aMolecule, fout, angleToRotate, axis );
    
    fout.close();
    
    
    return isFileWritten;
}


void MoleculeIOFunctions::writeAtomField( Molecule* aMolecule, ofstream& fout, const rg_REAL& angleToRotate, const rg_INT& axis )
{
    fout.setf(ios::fixed, ios::floatfield);     // TO 3 PLACES OF DECIMALS 
    fout.precision(3);                          //

    rg_TMatrix3D rotationMatrix;
    switch(axis) {
        case 0:
            rotationMatrix.rotateX(angleToRotate);
            break;
        case 1:
            rotationMatrix.rotateY(angleToRotate);
            break;
        case 2:
            rotationMatrix.rotateZ(angleToRotate);
            break;
        default:            
            break;
    }
    
    
    rg_Point3D translationVector = aMolecule->getCenterOfMass();

    rg_dList<Atom>* atoms = aMolecule->getAtoms();
    atoms->reset4Loop();

    Atom* currAtom = rg_NULL;
    while ( atoms->setNext4Loop() ) {
        currAtom = atoms->getpEntity();
        
        string recLine     = PDB_RECORD_TYPE_ATOM;
        string blankSpace  = "";
        string singleSpace = " ";

        string serial;
        string name;
        string altLoc;
        string resName;
        string chainID;
        string resSeq;
        string xCoord;
        string yCoord;
        string zCoord;
        string occupancy;
        string tempFactor;
        string element;
        string charge;

        
        // SERIAL FIELD
        serial     = StringFunctions::convertIntegerToString( currAtom->getSerialFromInputFile() );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_SERIAL_LENGTH, serial );
        recLine    = recLine + blankSpace + serial + singleSpace;
        

        // NAME FIELD
        if( !currAtom->getResidue()->isStandardResidue() ) {
            name = currAtom->getAtomNameFromInputFile();
        }

        else {
            string extraBranch = currAtom->getpChemicalProperties()->getExtraBrangeDesignatorInString();
            string atomName    = StringFunctions::strTrim( (string)ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol );
            string remoteInd   = currAtom->getpChemicalProperties()->getRemoteIndicatorInString();
            string branchDes   = currAtom->getpChemicalProperties()->getBrangeDesignatorInString();

            if( atomName.length() == 2 )
                extraBranch = "";

            name       = extraBranch + atomName + remoteInd + branchDes;
        }

        recLine    = recLine + name;

        // ALTLOC FIELD --> CURRENTLY, ALTLOC IS NOT APPLIED.
        altLoc     = " ";
        recLine    = recLine + altLoc;
        
        
        // RESIDUE NAME
        resName = currAtom->getResidue()->getResidueName();
        recLine = recLine + resName + singleSpace;


        // CHAIN ID
        rg_INT chainIDInInt  = currAtom->getResidue()->getChain()->getID() + 65;
        char   chainIDInChar =  chainIDInInt;
        chainID    = chainID + chainIDInChar;
        recLine    = recLine + chainID;


        // RESIDUE SEQUENCE
        resSeq     = StringFunctions::convertIntegerToString( currAtom->getResidue()->getSequenceNumber() );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_RESSEQ_LENGTH, resSeq );
        recLine    = recLine + blankSpace + resSeq;


        // ICODE --> CURRENTLY, ICODE IS NOT APPLIED.
        recLine    = recLine + "    ";

        rg_Point3D coordinateToRotate(currAtom->getpAtomBall()->getCenter().getX(),
                                      currAtom->getpAtomBall()->getCenter().getY(),
                                      currAtom->getpAtomBall()->getCenter().getZ());

        
        coordinateToRotate -= translationVector;
        coordinateToRotate = rotationMatrix * coordinateToRotate;
        coordinateToRotate += translationVector;



        // X COORDINATE
        xCoord     = StringFunctions::convertDoubleToString( coordinateToRotate.getX(), 3 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_X_LENGTH, xCoord );
        recLine    = recLine + blankSpace + xCoord;


        // Y COORDINATE
        yCoord     = StringFunctions::convertDoubleToString( coordinateToRotate.getY(), 3 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_Y_LENGTH, yCoord );
        recLine    = recLine + blankSpace + yCoord;

        
        // Y COORDINATE
        zCoord     = StringFunctions::convertDoubleToString( coordinateToRotate.getZ(), 3 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_Z_LENGTH, zCoord );
        recLine    = recLine + blankSpace + zCoord;


        // OCCUPANCY
        occupancy  = StringFunctions::convertDoubleToString( currAtom->getpChemicalProperties()->getOccupancy(), 2 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_OCCUPANCY_LENGTH, occupancy );
        recLine    = recLine + blankSpace + occupancy;


        // TEMPFACTOR
        tempFactor = StringFunctions::convertDoubleToString( currAtom->getpChemicalProperties()->getTempFactor(), 2 );
        blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_TEMPFACTOR_LENGTH, tempFactor );
        recLine    = recLine + blankSpace + tempFactor;

        // SEG_ID --> CURRENTLY, SEG_ID IS NOT APPLIED.
        recLine    = recLine + "          ";


        // ELEMENT
        element = ATOM_FEATURES[currAtom->getAtomCode()].symbol;
        recLine    = recLine + element;


        // CHARGE --> CURRENTLY, CHARGE IS NOT APPLIED.
        recLine    = recLine + "   ";


        fout << recLine << endl;
    }
}

void MoleculeIOFunctions::writeHeteroAtomField( Molecule* aMolecule, ofstream& fout, const rg_REAL& angleToRotate, const rg_INT& axis )
{
    
}



rg_FLAG MoleculeIOFunctions::writePDBFileForSCP( const char* pathFile, Molecule* aMolecule )
{
	rg_FLAG isFileWritten   = rg_FALSE;

	ofstream  fout( pathFile );

	writeHeader( aMolecule, fout );
	//writeAtomFieldInOrderOfSequence( aMolecule, fout );
	writeAtomFieldInOrderOfSequence2(aMolecule, fout);
	//writeConnectField( aMolecule, fout ); // Java Viewer를 이용할 경우 이 함수를 comment out해야 함.
	writeHeteroAtomField( aMolecule, fout );

	fout.close();
	
	return isFileWritten;

}

rg_FLAG MoleculeIOFunctions::writePDBFileForSCP(const char * pathFile, Molecule * aMolecule, Molecule * refereneceMolecule)
{
	rg_FLAG isFileWritten = rg_FALSE;

	rg_INT startSerialNumber = refereneceMolecule->getAtoms()->getHead()->getpEntity()->getSerialFromInputFile();
	reset_atom_serial_numbers(aMolecule, startSerialNumber);

	ofstream  fout(pathFile);

	writeHeader(aMolecule, fout);
	//writeAtomFieldInOrderOfSequence( aMolecule, fout );
	writeAtomFieldInOrderOfSequence2(aMolecule, fout);
	//writeConnectField( aMolecule, fout ); // Java Viewer를 이용할 경우 이 함수를 comment out해야 함.
	writeHeteroAtomField(aMolecule, fout);

	fout.close();

	return isFileWritten;
}

void MoleculeIOFunctions::reset_atom_serial_numbers(Molecule * aMolecule, const rg_INT& startNumber)
{
	rg_INT serialNumber = startNumber;
	rg_dList<Residue>* residues = aMolecule->getResidues();
	residues->reset4Loop();
	while (residues->setNext4Loop())
	{
		Residue* currResidue = residues->getpEntity();
		rg_dList<Atom*>* atoms = currResidue->getAtoms();
		atoms->reset4Loop();
		while (atoms->setNext4Loop())
		{
			Atom* currAtom = atoms->getEntity();
			currAtom->setSerialFromInputFile(serialNumber++);
		}
	}
}

void MoleculeIOFunctions::writeAtomFieldInOrderOfSequence( Molecule* aMolecule, ofstream& fout )
{
	fout.setf(ios::fixed, ios::floatfield);     // TO 3 PLACES OF DECIMALS 
	fout.precision(3);                          //

	rg_dList<Residue>* residues = aMolecule->getResidues();

	residues->reset4Loop();
	while(residues->setNext4Loop())
	{
		Residue* currResidue = residues->getpEntity();
		
		//if(currResidue->getResidueCode() == HOH_RESIDUE)
		//	continue;
        if(! currResidue->isAminoResidue())
            continue;

		rg_dList<Atom*>* atoms = currResidue->getAtoms();

		atoms->reset4Loop();		
		while ( atoms->setNext4Loop() ) 
		{
			Atom* currAtom = rg_NULL;
			currAtom = atoms->getEntity();

			// for SCP solver
			//if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
			//	continue;

			// test
			//ResidueCode code = currAtom->getResidue()->getResidueCode();			
			//if(code != PRO_AMINO_RESIDUE )
			//	continue;

			string recLine     = PDB_RECORD_TYPE_ATOM;
			string blankSpace  = "";
			string singleSpace = " ";

			string serial;
			string name;
			string altLoc;
			string resName;
			string chainID;
			string resSeq;
			string xCoord;
			string yCoord;
			string zCoord;
			string occupancy;
			string tempFactor;
			string element;
			string charge;


			// SERIAL FIELD
			serial     = StringFunctions::convertIntegerToString( currAtom->getSerialFromInputFile() );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_SERIAL_LENGTH, serial );
			recLine    = recLine + blankSpace + serial + singleSpace;


			// NAME FIELD
			if( !currAtom->getResidue()->isStandardResidue() ) {
				name = currAtom->getAtomNameFromInputFile();
			}

			else {
				string extraBranch = currAtom->getpChemicalProperties()->getExtraBrangeDesignatorInString();
				string atomName    = StringFunctions::strTrim( (string)ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol );
				string remoteInd   = currAtom->getpChemicalProperties()->getRemoteIndicatorInString();
				string branchDes   = currAtom->getpChemicalProperties()->getBrangeDesignatorInString();

				if( atomName.length() == 2 )
					extraBranch = "";

				name       = extraBranch + atomName + remoteInd + branchDes;
			}

			recLine    = recLine + name;

			// ALTLOC FIELD --> CURRENTLY, ALTLOC IS NOT APPLIED.
			altLoc     = " ";
			recLine    = recLine + altLoc;


			// RESIDUE NAME
			resName = currAtom->getResidue()->getResidueName();
			recLine = recLine + resName + singleSpace;


			// CHAIN ID
			// May 27, 2016 by Joonghyun
			//rg_INT chainIDInInt  = currAtom->getResidue()->getChain()->getID() + 65;
			//char   chainIDInChar =  chainIDInInt;
			//chainID    = chainID + chainIDInChar;
//#ifdef _DEBUG
			chainID = currAtom->getResidue()->getChain()->getChainIDFromInputFileInString();
//#endif
			recLine    = recLine + chainID;

			// RESIDUE SEQUENCE
			resSeq     = StringFunctions::convertIntegerToString( currAtom->getResidue()->getSequenceNumber() );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_RESSEQ_LENGTH, resSeq );
			recLine    = recLine + blankSpace + resSeq;


			// ICODE --> CURRENTLY, ICODE IS NOT APPLIED.
			recLine    = recLine + "    ";


			// X COORDINATE
			xCoord     = StringFunctions::convertDoubleToString( currAtom->getpAtomBall()->getCenter().getX(), 3 );
			//xCoord     = StringFunctions::convertDoubleToString( floorf( currAtom->getpAtomBall()->getCenter().getX() * 1000) / 1000, 3 );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_X_LENGTH, xCoord );
			recLine    = recLine + blankSpace + xCoord;


			// Y COORDINATE
			yCoord     = StringFunctions::convertDoubleToString( currAtom->getpAtomBall()->getCenter().getY(), 3 );
			//yCoord     = StringFunctions::convertDoubleToString( floorf( currAtom->getpAtomBall()->getCenter().getY() * 1000) / 1000, 3 );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_Y_LENGTH, yCoord );
			recLine    = recLine + blankSpace + yCoord;


			// Y COORDINATE
			zCoord     = StringFunctions::convertDoubleToString( currAtom->getpAtomBall()->getCenter().getZ(), 3 );
			//zCoord     = StringFunctions::convertDoubleToString( floorf(currAtom->getpAtomBall()->getCenter().getZ() * 1000) / 1000, 3 );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_Z_LENGTH, zCoord );
			recLine    = recLine + blankSpace + zCoord;


			// OCCUPANCY
			occupancy  = StringFunctions::convertDoubleToString( currAtom->getpChemicalProperties()->getOccupancy(), 2 );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_OCCUPANCY_LENGTH, occupancy );
			recLine    = recLine + blankSpace + occupancy;


			// TEMPFACTOR
			tempFactor = StringFunctions::convertDoubleToString( currAtom->getpChemicalProperties()->getTempFactor(), 2 );
			blankSpace = StringFunctions::getBlankStringForFixedField( PDB_ATOM_RECORD_TEMPFACTOR_LENGTH, tempFactor );
			recLine    = recLine + blankSpace + tempFactor;

			// SEG_ID --> CURRENTLY, SEG_ID IS NOT APPLIED.
			recLine    = recLine + "          ";


			// ELEMENT
			element = ATOM_FEATURES[currAtom->getAtomCode()].symbol;
			recLine    = recLine + element;


			// CHARGE --> CURRENTLY, CHARGE IS NOT APPLIED.
			recLine    = recLine + "   ";


			fout << recLine << endl;
		}
	}
}



void MoleculeIOFunctions::writeAtomFieldInOrderOfSequence2(Molecule * aMolecule, ofstream & fout)
{
    fout.setf(ios::fixed, ios::floatfield);     // TO 3 PLACES OF DECIMALS 
    fout.precision(3);                          //

    rg_dList<Residue>* residues = aMolecule->getResidues();

    residues->reset4Loop();
    while (residues->setNext4Loop())
    {
        Residue* currResidue = residues->getpEntity();

        //if(currResidue->getResidueCode() == HOH_RESIDUE)
        //	continue;
        if (!currResidue->isAminoResidue())
        	continue;

        rg_dList<Atom*>* atoms = currResidue->getAtoms();

        atoms->reset4Loop();
        while (atoms->setNext4Loop())
        {
            Atom* currAtom = rg_NULL;
            currAtom = atoms->getEntity();

            // for SCP solver
            //if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
            //	continue;

            // test
            //ResidueCode code = currAtom->getResidue()->getResidueCode();			
            //if(code != PRO_AMINO_RESIDUE )
            //	continue;

            string recLine = PDB_RECORD_TYPE_ATOM;
            string blankSpace = "";
            string singleSpace = " ";

            string serial;
            string name;
            string altLoc;
            string resName;
            string chainID;
            string resSeq;
            string xCoord;
            string yCoord;
            string zCoord;
            string occupancy;
            string tempFactor;
            string element;
            string charge;


            // SERIAL FIELD
            serial = StringFunctions::convertIntegerToString(currAtom->getSerialFromInputFile());
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_SERIAL_LENGTH, serial);
            recLine = recLine + blankSpace + serial + singleSpace;


            // NAME FIELD
            if (!currAtom->getResidue()->isStandardResidue()) {
                name = currAtom->getAtomNameFromInputFile();
            }

            else {
                string extraBranch = currAtom->getpChemicalProperties()->getExtraBrangeDesignatorInString();
                string atomName = StringFunctions::strTrim((string)ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol);
                string remoteInd = currAtom->getpChemicalProperties()->getRemoteIndicatorInString();
                string branchDes = currAtom->getpChemicalProperties()->getBrangeDesignatorInString();

                if (atomName.length() == 2)
                    extraBranch = "";

                name = extraBranch + atomName + remoteInd + branchDes;
            }

            recLine = recLine + name;

            // ALTLOC FIELD --> CURRENTLY, ALTLOC IS NOT APPLIED.
            altLoc = " ";
            recLine = recLine + altLoc;


            // RESIDUE NAME
			if (currAtom->getResidue()->isStandardResidue()) {
				resName = (string)RESIDUE_FEATURES[currAtom->getResidue()->getResidueCode()].threeCodeName;
				//resName = currResidue->getResidueName();
			}
			else {
				string residueName = currAtom->getResidue()->getResidueName();
				if (residueName.length() >= PDB_ATOM_RECORD_RESNAME_LENGTH) {
					resName = residueName.substr(0, 3);
				}
				else {
					resName = residueName;
					int numSpaces = PDB_ATOM_RECORD_RESNAME_LENGTH - residueName.length();
					for (int i = 0; i < numSpaces; ++i) {
						resName += " ";
					}
				}
			}
            recLine = recLine + resName + singleSpace;


            // CHAIN ID
            // May 27, 2016 by Joonghyun
            //rg_INT chainIDInInt  = currAtom->getResidue()->getChain()->getID() + 65;
            //char   chainIDInChar =  chainIDInInt;
            //chainID    = chainID + chainIDInChar;
            //#ifdef _DEBUG
            chainID = currAtom->getResidue()->getChain()->getChainIDFromInputFileInString();
            //#endif
            recLine = recLine + chainID;

            // RESIDUE SEQUENCE
            resSeq = StringFunctions::convertIntegerToString(currAtom->getResidue()->getSequenceNumber());
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_RESSEQ_LENGTH, resSeq);
            recLine = recLine + blankSpace + resSeq;


            // ICODE --> CURRENTLY, ICODE IS NOT APPLIED.
            recLine = recLine + "    ";


            // X COORDINATE
            xCoord = StringFunctions::convertDoubleToString(currAtom->getpAtomBall()->getCenter().getX(), 3);
            //xCoord     = StringFunctions::convertDoubleToString( floorf( currAtom->getpAtomBall()->getCenter().getX() * 1000) / 1000, 3 );
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_X_LENGTH, xCoord);
            recLine = recLine + blankSpace + xCoord;


            // Y COORDINATE
            yCoord = StringFunctions::convertDoubleToString(currAtom->getpAtomBall()->getCenter().getY(), 3);
            //yCoord     = StringFunctions::convertDoubleToString( floorf( currAtom->getpAtomBall()->getCenter().getY() * 1000) / 1000, 3 );
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_Y_LENGTH, yCoord);
            recLine = recLine + blankSpace + yCoord;


            // Y COORDINATE
            zCoord = StringFunctions::convertDoubleToString(currAtom->getpAtomBall()->getCenter().getZ(), 3);
            //zCoord     = StringFunctions::convertDoubleToString( floorf(currAtom->getpAtomBall()->getCenter().getZ() * 1000) / 1000, 3 );
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_Z_LENGTH, zCoord);
            recLine = recLine + blankSpace + zCoord;


            // OCCUPANCY
            occupancy = StringFunctions::convertDoubleToString(currAtom->getpChemicalProperties()->getOccupancy(), 2);
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_OCCUPANCY_LENGTH, occupancy);
            recLine = recLine + blankSpace + occupancy;


            // TEMPFACTOR
            tempFactor = StringFunctions::convertDoubleToString(currAtom->getpChemicalProperties()->getTempFactor(), 2);
            blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_TEMPFACTOR_LENGTH, tempFactor);
            recLine = recLine + blankSpace + tempFactor;

            // SEG_ID --> CURRENTLY, SEG_ID IS NOT APPLIED.
            recLine = recLine + "          ";


            // ELEMENT
            element = ATOM_FEATURES[currAtom->getAtomCode()].symbol;
            recLine = recLine + element;


            // CHARGE --> CURRENTLY, CHARGE IS NOT APPLIED.
            recLine = recLine + "   ";


            fout << recLine << endl;
        }
    }
}


rg_FLAG MoleculeIOFunctions::readNextMoleculeFromGivenInputStream( ifstream& fin, Molecule& aMolecule, rg_FLAG& isNextExist )
{	
	isNextExist = rg_FALSE;
	string strRecLine  = "";

	rg_INT arrNumOfRecords[5] = { 0, 0, 0, 0, 0 };
	setMoleculeRTIToMolecule( fin, arrNumOfRecords, aMolecule, strRecLine);


	rg_FLAG isRecLinesOK   = rg_TRUE;

	AtomMap     * mapOfAtom     = new AtomMap;
	ResidueMap  * mapOfResidue  = new ResidueMap;
	ChainMap    * mapOfChain    = new ChainMap;
	ChemBondMap * mapOfChemBond = new ChemBondMap;

	AtomSymbolMap    *mapOfAtomSymbol    = new AtomSymbolMap;
	ResidueSymbolMap *mapOfResidueSymbol = new ResidueSymbolMap;

	// CAN'T READ MOL2 FILE... NEED TO DEBUG....
	ResidueMap *mapOfNonStdResidue = new ResidueMap;
	
	initiallizeSymbolMapsOfAtomAndResidueForMol2File( *mapOfAtomSymbol, *mapOfResidueSymbol );

	string currRecType = "";

	do {
		string recType = StringFunctions::strTrim( strRecLine );

		if ( recType == MOL2_RECORD_TYPE_ATOM && arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID] > 0 ) {
			//            CAN'T READ MOL2 FILE... NEED TO DEBUG....
			isRecLinesOK = setAtomRTIToMolecule( fin, arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID], *mapOfAtom, *mapOfResidue, *mapOfNonStdResidue, *mapOfAtomSymbol, *mapOfResidueSymbol, aMolecule, strRecLine );
			recType = StringFunctions::strTrim( strRecLine );
			//            isRecLinesOK = setAtomRTIToMolecule( i_recLines, &recordsOfMol2File, arrNumOfRecords[MOL2_NUM_OF_ATOMS_ID], *mapOfAtom, *mapOfResidue, *mapOfAtomSymbol, *mapOfResidueSymbol, aMolecule );
		}

		if(isRecLinesOK == rg_FALSE) {
			break;
		}

		if ( recType == MOL2_RECORD_TYPE_BOND && arrNumOfRecords[MOL2_NUM_OF_BONDS_ID] >0 ) {
			isRecLinesOK = setBondRTIToMolecule( fin, arrNumOfRecords[MOL2_NUM_OF_BONDS_ID], *mapOfAtom, aMolecule, strRecLine );
			recType = StringFunctions::strTrim( strRecLine );
		}

		if( isRecLinesOK == rg_FALSE ) {			
			break;
		}

		if ( recType == MOL2_RECORD_TYPE_MOLECULE ) {	
			isNextExist = rg_TRUE;
			break;
		}

	} while ( getline( fin, strRecLine ) );
	
	
		

	//    int numOfNonStdResidue = mapOfNonStdResidue.size();

	// CAN'T READ MOL2 FILE... NEED TO DEBUG....
	//setChainsForNonStdResiduesForMol2File( &mapOfNonStdResidue, aMolecule );
	setChainsForNonStdResiduesForMol2File( mapOfNonStdResidue, aMolecule );
//	checkChainsInMoleculeForStandardResidue( aMolecule );
    evaluateAndSetChainCode( aMolecule );

	// SET FILE NAME TO MOLECULE
//	string fileName = StringFunctions::getFileNameWithoutPath( string(pathFile) ) ;
//	aMolecule.setMoleculeFileName( fileName );


	if( isRecLinesOK ) {
		aMolecule.evaluateHydrogenDonorAndAcceptorAtoms();
		aMolecule.computeAndSetCenterOfMass();
		aMolecule.computeAndSetMinEnclosingSphere();	
	}

	delete mapOfAtom;
	delete mapOfResidue;
	delete mapOfChain;
	delete mapOfChemBond;
	delete mapOfNonStdResidue;

	delete mapOfAtomSymbol;   
	delete mapOfResidueSymbol;

	return isRecLinesOK;

}


rg_FLAG MoleculeIOFunctions::setMoleculeRTIToMolecule( ifstream& fin, rg_INT* arrNumOfRecords, Molecule& aMolecule, string& strRecLine )
{
	rg_FLAG isMoleculeRTIOK = rg_TRUE;

	string  nameOfMolecule       = "";
	rg_BOOL isNameOfMoleculeRead = rg_FALSE;
	rg_BOOL isNumOfRecordsRead   = rg_FALSE;


	while( getline(fin, strRecLine) ) {

		string firstLetter = strRecLine.substr(0, 1);
		if( firstLetter == " " || firstLetter == "" ) {
			continue;
		}

		if(isNewTriposRTILineStart( strRecLine )) {
			break;
		}


		// MOLECULE NAME
		if( isNameOfMoleculeRead == rg_FALSE /* && StringFunctions::isStringLetterInAlphabet( firstLetter ) */ ) {
			nameOfMolecule = StringFunctions::strTrim( strRecLine );

			aMolecule.setMoleculeFileName(nameOfMolecule);

			getline(fin, strRecLine);
			strRecLine = StringFunctions::strTrim( strRecLine );

			firstLetter = strRecLine.substr(0, 1);
			isNameOfMoleculeRead = rg_TRUE;
		}

		// NUM OF ATOMS/BONDS/SUBST/FEAT/SETS
		if( isNumOfRecordsRead == rg_FALSE && StringFunctions::isStringLetterInNumber( firstLetter ) ) {
			rg_INT i_numOfRec = 0;
			while( strRecLine.length() != 0 ) {

				rg_INT posOfBlank = strRecLine.find( " ", 0 );

				string numOfRecs = strRecLine.substr( 0, posOfBlank );

				arrNumOfRecords[i_numOfRec] = atoi( numOfRecs.c_str() );
				i_numOfRec++;

				strRecLine.erase( 0, numOfRecs.length() );
				strRecLine = StringFunctions::strTrimLeft( strRecLine );
			}

			isNumOfRecordsRead = rg_TRUE;
		}

		// MOLECULE TYPE : TO BE CONSIDERED...
		// CHARGE TYPE   : TO BE CONSIDERED...
	}

	if ( isNameOfMoleculeRead == rg_FALSE && isNumOfRecordsRead == rg_FALSE )
		isMoleculeRTIOK = rg_FALSE;

	return isMoleculeRTIOK;
}

rg_FLAG MoleculeIOFunctions::setAtomRTIToMolecule( ifstream& fin, const rg_INT& numOfAtomsGiven, AtomMap& mapOfAtom, ResidueMap& mapOfResidue, ResidueMap& mapOfNonStdResidue, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule& aMolecule, string& strRecLine )
{
	rg_FLAG isAtomRTIOK = rg_TRUE;

	getline(fin, strRecLine);

	strRecLine  = StringFunctions::strTrim( strRecLine );

	rg_INT numOfChains = aMolecule.getChains()->getSize();
	Chain* newChain    = aMolecule.addChain( Chain(numOfChains, &aMolecule) );
	rg_BOOL isNewChainStart = rg_FALSE;


	rg_INT i_atom = 0;

	do {
		if(isNewTriposRTILineStart( strRecLine )) {
			break;
		}

		i_atom++;

		rg_INT  arrIntRecs[MOL2_ATOM_INT_REC_SIZE]  = { -1, -1 };
		rg_REAL arrRealRecs[MOL2_ATOM_REAL_REC_SIZE] = { 0.0, 0.0, 0.0 };
		string  arrStrRecs[MOL2_ATOM_STR_REC_SIZE]  = { "", "", "", "" };

		extractAtomRecordsFromAtomRTILine( strRecLine, arrIntRecs, arrRealRecs, arrStrRecs );

		if( arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID] == 0 && arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] == 0 ) {
			break;
		}

		// Estimate Chain : if name of prev atom was "OXT" then new Chain must be started.
		if( isNewChainStart == rg_TRUE ) {
			numOfChains = aMolecule.getChains()->getSize();
			newChain    = aMolecule.addChain( Chain(numOfChains, &aMolecule) );
			isNewChainStart = rg_FALSE;
		}


		// Estimate Residue : if residue is not in map, then create a new residue.
		Residue* currResidue = getResiduePtrFromMap( &mapOfResidue,  arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] );

		if ( currResidue == rg_NULL ) {
			rg_INT newResidueID = mapOfResidue.size();
			currResidue = aMolecule.addResidue( Residue(newResidueID) );
			currResidue->setSequenceNumber( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID] );
			currResidue->setResidueName( arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME] );

			int aaa = MOL2_ATOM_INT_REC_SUBST_ID;

			string residueName = getResidueNameFromTriposSubstName( arrStrRecs[MOL2_ATOM_STR_REC_SUBST_NAME] );
			setResidueCodeToTargetResidue( mapOfResidueSymbol, residueName, currResidue );
			mapOfResidue.insert( ResidueMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID], currResidue ) );

			if( currResidue->isStandardResidue() ) {
				newChain->addResidue( currResidue );
				currResidue->setChain( newChain );
			}
			else {
				mapOfNonStdResidue.insert( ResidueMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_SUBST_ID], currResidue ) );
			}
		}

		// Estimate Atom
		rg_INT newAtomID = mapOfAtom.size();

		Atom* newAtom = aMolecule.addAtom( Atom( newAtomID ) );


		currResidue->addAtom( newAtom );
		newAtom->setResidue( currResidue );
		newAtom->setSerialFromInputFile( arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID] );
		newAtom->setAtomNameFromInputFile( arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME] );
		rg_FLAG isAtomNameOK = setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForMol2File( mapOfAtomSymbol, arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME], arrStrRecs[MOL2_ATOM_STR_REC_ATOM_TYPE], newAtom );

		if( arrStrRecs[MOL2_ATOM_STR_REC_ATOM_NAME] == "OXT") {
			isNewChainStart = rg_TRUE;
		}

		if ( !isAtomNameOK ) {
			isAtomRTIOK = rg_FALSE;
			break;
		}

		newAtom->setAtomBall( Sphere( arrRealRecs[MOL2_ATOM_REAL_REC_X_COORD], arrRealRecs[MOL2_ATOM_REAL_REC_Y_COORD], arrRealRecs[MOL2_ATOM_REAL_REC_Z_COORD], ATOM_FEATURES[newAtom->getAtomCode()].radius) );
		newAtom->getpChemicalProperties()->setCharge( arrRealRecs[MOL2_ATOM_REAL_REC_CHARGE] );

		mapOfAtom.insert( AtomMap::value_type( arrIntRecs[MOL2_ATOM_INT_REC_ATOM_ID], newAtom ) );


		if( i_atom == numOfAtomsGiven )
			break;

	} while( getline(fin, strRecLine) );

	
	// Filter chains with no residue and reset ID.
	filterEmptyChain( aMolecule );

	if( i_atom != numOfAtomsGiven )
		isAtomRTIOK = rg_FALSE;

	return isAtomRTIOK;
}


rg_FLAG MoleculeIOFunctions::setBondRTIToMolecule( ifstream& fin, const rg_INT& numOfBondsGiven, AtomMap& mapOfAtom, Molecule& aMolecule, string& strRecLine )
{
	rg_FLAG isBondRTIOK = rg_TRUE;

	
	getline(fin, strRecLine);

	strRecLine  = StringFunctions::strTrim( strRecLine );

	rg_INT i_bond = 0;

	do {
		if(isNewTriposRTILineStart( strRecLine )) {
			break;
		}

		i_bond++;

		rg_INT bondID       = 0;
		rg_INT originAtomID = 0;
		rg_INT targetAtomID = 0;
		string bondType     = "";
		string statusBit    = "";

		Atom* originAtom = rg_NULL;
		Atom* targetAtom = rg_NULL;

		extractBondRecordsFromBondRTILine( strRecLine, bondID, originAtomID, targetAtomID, bondType, statusBit );


		AtomMap::iterator AtomMap_i = mapOfAtom.find( originAtomID );
		AtomMap::iterator AtomMap_j = mapOfAtom.find( targetAtomID );

		if ( AtomMap_i != mapOfAtom.end() && AtomMap_j != mapOfAtom.end() ) {
			originAtom = (*AtomMap_i).second;
			targetAtom = (*AtomMap_j).second;

			setChemicalBondToMoleculeWithoutDuplication( originAtom, targetAtom, aMolecule );
		}
		else {
			// Just ignore when the originAtom or targetAtom does not exist in Atom field.
			// isBondRTIOK = rg_FALSE;
			//break;
		}


		if( i_bond == numOfBondsGiven )
			break;

	}
	while( getline(fin, strRecLine) );
	

	if( i_bond != numOfBondsGiven )
		isBondRTIOK = rg_FALSE;

	return isBondRTIOK;
}



///////////////////////////////////////////////////////////////////////////////
//    
// .cfg via RMC simulation using the software RMCA
//
rg_FLAG MoleculeIOFunctions::readRMCCFGFile( const string& pathFile, Molecule& aMolecule )
{
    rg_FLAG isFileRead   = rg_FALSE;

    if( !isValidMoleculeFileExtension( pathFile.c_str(), RMC_CFG_FILE ) ) {
        cerr << "Error: Invalid file extension !\n";
        return isFileRead;
    }
    
    ifstream fin( pathFile.c_str() );    
    if( fin.bad() ) {
        cerr << "Error: Could not open file!\n";
        return isFileRead;
    }

    string strRecLine  = "";
    if( !getline( fin, strRecLine ) ) {
        return isFileRead;
    }
    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    
    list<string> recordsOfMolFile;
    addRecordLinesOfFIleToList( &fin, recordsOfMolFile );
    fin.close();




    map<string, AtomCode> atomSymbolMap;   
    for ( int i_atomCode = UNK_ATOM; i_atomCode<NUM_OF_ATOMS; ++i_atomCode ) {
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        atomSymbolMap.insert( make_pair( StringFunctions::strTrim( atomSymbol ), (AtomCode)i_atomCode ) );
    }

    list<string>::iterator  i_record;
    vector<string>          atomName;
    vector<float>           atomCapacityInMolecule;
    int                     totalNumAtoms;
    vector<int>             numAtomsOfEachType;
    vector<rg_Point3D>      definingVecotr;
    readHeaderBlockOfRMCCFGFile( recordsOfMolFile, i_record, atomName, atomCapacityInMolecule, totalNumAtoms, numAtomsOfEachType, definingVecotr);


    rg_Point3D vector[3] = { definingVecotr[0], definingVecotr[1], definingVecotr[2] };


    while (true) {
        if (i_record->empty()) {
            ++i_record;
        }
        else {
            break;
        }
    }


    int i_type = 0;
    int count  = 1;
    for ( i_type=0; i_type<atomName.size(); ++i_type ) {
        AtomCode atomCode;
        map<string, AtomCode>::iterator i_code = atomSymbolMap.find( atomName[i_type] );
        if ( i_code == atomSymbolMap.end() ) {
            atomCode = UNK_ATOM;
        }
        else {
            atomCode = i_code->second;
        }

        Residue*   currResidue = aMolecule.addResidue( Residue(i_type+1, UNK_RESIDUE) );
        Chain*     currChain   = aMolecule.addChain( Chain(i_type+1, UNK_CHAIN) );
        currChain->addResidue( currResidue );
        currResidue->setChain( currChain );

        int numAtoms = numAtomsOfEachType[i_type];
        int i_atom   = 0;
        while ( i_atom<numAtoms && i_record!=recordsOfMolFile.end() ) {
            string currRecord = *i_record;
            if (currRecord.empty()) {
                continue;
            }

            istringstream recordFin(currRecord);
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            recordFin >> x >> y >> z;

            //string strX = StringFunctions::pop_numeric_char_at_front( currRecord );
            //string strY = StringFunctions::pop_numeric_char_at_front( currRecord );
            //string strZ = StringFunctions::pop_numeric_char_at_front( currRecord );

            //float  x = atof( strX.c_str() );
            //float  y = atof( strY.c_str() );
            //float  z = atof( strZ.c_str() );

            rg_Point3D position( x*vector[0].getX(), y*vector[1].getY(), z*vector[2].getZ() );
            double     radius = ATOM_FEATURES[atomCode].radius;

            Atom* currAtom = aMolecule.addAtom( Atom( count ) );
            currAtom->setAtomBall( Sphere( position, radius ) );
            currAtom->setAtomCode( atomCode );
            currAtom->setSerialFromInputFile( count );
            currAtom->setResidue( currResidue );

            currResidue->addAtom( currAtom );

            ++count;
            ++i_atom;
            ++i_record;
        }
    }

    return rg_TRUE;
}
        




void    MoleculeIOFunctions::readHeaderBlockOfRMCCFGFile( 
                                             list<string>&              recordsOfMolFile, 
                                             list<string>::iterator&    i_record,
                                             vector<string>&            atomName,
                                             vector<float>&             atomCapacityInMolecule,
                                             int&                       totalNumAtoms,
                                             vector<int>&               numAtomsOfEachType,
                                             vector<rg_Point3D>&        definingVecotr)
{
    const string KEYWORD_OF_LINE_FOR_ATOM_TYPE("S(Q)'s from MC-Simulation");
    const string KEYWORD_OF_LINE_FOR_TOTAL_ATOM_NUM("molecules (of all types)");
    const string KEYWORD_OF_LINE_FOR_ATOM_NUM_OF_TYPE("molecules of type");
    const string KEYWORD_OF_LINE_FOR_BOX("Defining vectors are:");

    string::size_type   indexKeyword;
    bool                continueHeaderParsing = true;

    i_record = recordsOfMolFile.begin();
    for (; continueHeaderParsing && i_record != recordsOfMolFile.end(); ++i_record) {
        string currRecord = *i_record;

        indexKeyword = currRecord.find(KEYWORD_OF_LINE_FOR_ATOM_TYPE);
        if (indexKeyword != string::npos) {
            parseAtomNameAndCapacityInMolecule(*i_record, atomName, atomCapacityInMolecule);
            continue;
        }

        indexKeyword = currRecord.find(KEYWORD_OF_LINE_FOR_TOTAL_ATOM_NUM);
        if (indexKeyword != string::npos) {
            totalNumAtoms = atoi(StringFunctions::pop_numeric_char_at_front(currRecord).c_str());
            continue;
        }

        indexKeyword = currRecord.find(KEYWORD_OF_LINE_FOR_BOX);
        if (indexKeyword != string::npos) {
            for (int i = 0; i<3; ++i) {
                ++i_record;
                string strDefiningVector = *i_record;
                string coordOfVector[3];
                coordOfVector[0] = StringFunctions::pop_numeric_char_at_front(strDefiningVector);
                coordOfVector[1] = StringFunctions::pop_numeric_char_at_front(strDefiningVector);
                coordOfVector[2] = StringFunctions::pop_numeric_char_at_front(strDefiningVector);

                definingVecotr.push_back(rg_Point3D(atof(coordOfVector[0].c_str()), atof(coordOfVector[1].c_str()), atof(coordOfVector[2].c_str())));
            }
            continue;
        }

        indexKeyword = currRecord.find(KEYWORD_OF_LINE_FOR_ATOM_NUM_OF_TYPE);
        if (indexKeyword != string::npos) {
            numAtomsOfEachType.push_back(atoi(StringFunctions::pop_numeric_char_at_front(currRecord).c_str()));
            if (numAtomsOfEachType.size() == atomName.size()) {
                ++i_record;
                ++i_record;
                ++i_record;
                continueHeaderParsing = false;
            }
            continue;
        }
    }


    /*
    int i=0;

    i_record = recordsOfMolFile.begin();
    ++i_record;

    parseAtomNameAndCapacityInMolecule( *i_record, atomName, atomCapacityInMolecule );

    for ( i=0; i<6; ++i ) {
        ++i_record;
    }

    string strTotalNumAtoms = *i_record;
    totalNumAtoms = atoi( StringFunctions::pop_numeric_char_at_front( strTotalNumAtoms ).c_str() );

    string strNumTypesOfMolecule = *(++i_record);
    int numTypesOfMolecule = atoi( StringFunctions::pop_numeric_char_at_front( strNumTypesOfMolecule ).c_str() );

    for ( i=0; i<6; ++i ) {
        ++i_record;
    }

    for ( i=0; i<3; ++i, ++i_record ) {
        string strDefiningVector = *i_record;
        string coordOfVector[3];
        coordOfVector[0] = StringFunctions::pop_numeric_char_at_front( strDefiningVector );
        coordOfVector[1] = StringFunctions::pop_numeric_char_at_front( strDefiningVector );
        coordOfVector[2] = StringFunctions::pop_numeric_char_at_front( strDefiningVector );

        definingVecotr.push_back( rg_Point3D( atof( coordOfVector[0].c_str() ), atof( coordOfVector[1].c_str() ), atof( coordOfVector[2].c_str() ) ) );
    }

    ++i_record;
    for (int i_type=0; i_type<numTypesOfMolecule; ++i_type ) {
        string strNumMolType = *i_record;
        numAtomsOfEachType.push_back( atoi( StringFunctions::pop_numeric_char_at_front( strNumMolType ).c_str() ) );
        for ( i=0; i<4; ++i ) {
            ++i_record;
        }
    }
    */
}


        
void    MoleculeIOFunctions::parseAtomNameAndCapacityInMolecule( const string&              record, 
                                             vector<string>&            atomName,
                                             vector<float>&             atomCapacityInMolecule)
{
    string limit = ", ";
    string currRecord = record;
    StringFunctions::pop_white_space_at_front(currRecord);
    string::size_type pos = currRecord.length();
    for (int i = 0; i < limit.size(); ++i) {
        string::size_type currPos = currRecord.find_first_of( limit[i] );
        if (pos > currPos) {
            pos = currPos;
        }
    }
    string atomCode = currRecord.substr(0, pos);

    while ( !atomCode.empty() ) {
        string strAtomName = StringFunctions::pop_alphabetic_char_at_front( atomCode );
        strAtomName = StringFunctions::strToUpper( strAtomName );
        atomName.push_back( strAtomName );


        string strCapacity = StringFunctions::pop_numeric_char_at_front( atomCode );
        if (strCapacity.empty()) {
            atomCapacityInMolecule.push_back(1);
        }
        else {
            atomCapacityInMolecule.push_back(atof(strCapacity.c_str()));
        }
    }
}



#ifdef WIN32
rg_FLAG MoleculeIOFunctions::readQuasiCrystalDataFile(const string& pathFile, Molecule& aMolecule)
{
    map<string, AtomCode> atomSymbolMap;
    for (int i_atomCode = UNK_ATOM; i_atomCode<NUM_OF_ATOMS; ++i_atomCode) {
        string atomSymbol(ATOM_FEATURES[i_atomCode].symbol);
        atomSymbolMap.insert(make_pair(StringFunctions::strTrim(atomSymbol), (AtomCode)i_atomCode));
    }

    map<int, float> atomRadiusMap;
    atomRadiusMap.insert(make_pair(13, 1.25f));
    atomRadiusMap.insert(make_pair(25, 1.40f));
    atomRadiusMap.insert(make_pair(46, 1.40f));


    const int   BUFFER_SIZE = 300;
    char	    lineData[BUFFER_SIZE];
    const char*	seps = " \t\n";
    char*	    token[10] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
    char*       context = NULL;

    ifstream fin;
    fin.open(pathFile.c_str());


    list<string> all_records;
    while (!fin.eof()) {
        fin.getline(lineData, BUFFER_SIZE);
        all_records.push_back(lineData);
    }
    all_records.pop_back();
    all_records.pop_back();
    all_records.pop_back();


    Residue*   currResidue = aMolecule.addResidue(Residue(1, UNK_RESIDUE));
    Chain*     currChain = aMolecule.addChain(Chain(1, UNK_CHAIN));
    currChain->addResidue(currResidue);
    currResidue->setChain(currChain);


    int count = 1;
    for (list<string>::iterator i_record = all_records.begin(); i_record != all_records.end(); ++i_record, ++count) {
        char* record = const_cast<char*>(i_record->c_str());
        token[0] = strtok_s(record, seps, &context);
        token[1] = strtok_s(NULL, seps, &context);
        token[2] = strtok_s(NULL, seps, &context);
        token[3] = strtok_s(NULL, seps, &context);
        token[4] = strtok_s(NULL, seps, &context);
        token[5] = strtok_s(NULL, seps, &context);
        token[6] = strtok_s(NULL, seps, &context);
        token[7] = strtok_s(NULL, seps, &context);
        token[8] = strtok_s(NULL, seps, &context);
        token[9] = strtok_s(NULL, seps, &context);


        if (token[0] == NULL || token[1] == NULL || token[2] == NULL || token[3] == NULL || token[4] == NULL ||
            token[5] == NULL || token[6] == NULL || token[7] == NULL || token[8] == NULL || token[9] == NULL) {
            continue;
        }


        string atomName = StringFunctions::strToUpper(token[0]);
        float  x = atof(token[1]);
        float  y = atof(token[2]);
        float  z = atof(token[3]);
        float  lattice1 = atof(token[4]);
        float  lattice2 = atof(token[5]);
        float  lattice3 = atof(token[6]);
        float  lattice4 = atof(token[7]);
        float  lattice5 = atof(token[8]);
        float  lattice6 = atof(token[9]);

        rg_Point3D position(x, y, z);
        float      radius = 0.0;
        AtomCode atomCode;
        map<string, AtomCode>::iterator i_code = atomSymbolMap.find(atomName);
        if (i_code == atomSymbolMap.end()) {
            atomCode = UNK_ATOM;
        }
        else {
            atomCode = i_code->second;
            map<int, float>::iterator i_radius = atomRadiusMap.find(atomCode);
            if (i_radius != atomRadiusMap.end()) {
                radius = i_radius->second;
            }
        }

        Atom* currAtom = aMolecule.addAtom(Atom(count));
        currAtom->setAtomBall(Sphere(position, radius));
        currAtom->setAtomCode(atomCode);
        currAtom->setAtomNameFromInputFile(atomName);
        currAtom->setSerialFromInputFile(count);
        currAtom->setResidue(currResidue);

        currResidue->addAtom(currAtom);


        //m_data.push_back(QuasiCrystalData(atomName, rg_Point3D(x, y, z), lattice1, lattice2, lattice3, lattice4, lattice5, lattice6));
    }


    if (count >= 4)
        return rg_TRUE;
    else
        return rg_FALSE;

}
#endif



// PDB IO Functions
rg_FLAG MoleculeIOFunctions::readPDBQFile(const string& pathFile, rg_dList<Molecule>& molecules)
{
    rg_FLAG isFileRead = rg_FALSE;

    if (!isValidMoleculeFileExtension(pathFile.c_str(), PDBQ_FILE)) {
        cerr << "Error: Invalid file extension !\n";
        return isFileRead;
    }


    ifstream fin(pathFile.c_str());

    if (fin.bad()) {
        cerr << "Error: Could not open file!\n";
        return isFileRead;
    }

    list<string>  recordLinesOfPDBFile;
    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    addRecordsOfPDBFIleToList(&fin, recordLinesOfPDBFile);
    rg_INT numOfRecLines = recordLinesOfPDBFile.size();
    fin.close();

    string strPathFile = pathFile;
    rg_INT fileSize = FileOperations::getFileSize(strPathFile);
    string depDate = parseDepositionDate(recordLinesOfPDBFile);



    AtomSymbolMap    mapOfAtomSymbol;
    //ResidueSymbolMap mapOfResidueSymbol;
    ResidueSymbolMap mapOfResidueSymbolForProtein;
    ResidueSymbolMap mapOfResidueSymbolForDNA;
    ResidueSymbolMap mapOfResidueSymbolForRNA;

    //initiallizeSymbolMapsOfAtomAndResidueForPDBFile( mapOfAtomSymbol, mapOfResidueSymbol );
    initiallizeSymbolMapsOfAtomAndResidueForPDBFile(mapOfAtomSymbol, mapOfResidueSymbolForProtein,
        mapOfResidueSymbolForDNA, mapOfResidueSymbolForRNA);


    list<string*> connectRecLinesOfPDBFile;

    // Secondary Structure Section 
    list<string*> helixRecLinesOfPDBFile;
    list<string*> sheetRecLinesOfPDBFile;
    list<string*> turnRecLinesOfPDBFile;


    AtomMap          mapOfAtom;
    //ResidueMap       mapOfResidue;
    ResidueMap       mapOfNonStdResidue;    // ?

    rg_dList<ResidueMapWithChainID> mapsOfResidue;
    rg_dList<ResidueMapWithChainID> mapsOfNonStdResidue;

    ChainMap         mapOfChain;
    ChainMap         mapOfChainForHeteroRec;    // first: 0 for HOH, and 1~ for others.
    ChemBondMap      mapOfChemBond;

    rg_INT    countModelRec = 1;
    Molecule* aMolecule = molecules.addTail(Molecule());
    aMolecule->setModelSerialFromInputFile(countModelRec);
    aMolecule->setTimeStamp(depDate);
    aMolecule->setFileSize(fileSize);
    aMolecule->setMoleculeFileName(pathFile);

    list<string>::iterator i_recLines = recordLinesOfPDBFile.begin();

    while (i_recLines != recordLinesOfPDBFile.end()) {
        string* currRecLine = &(*i_recLines);

        string recType = StringFunctions::subString(*currRecLine, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH);


        if (recType == PDB_RECORD_TYPE_ATOM) {
            rg_FLAG isRecLineOK = setAtomRecordsToMoleculeForPDBQ(*currRecLine, mapOfAtom, mapsOfResidue, mapOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbolForProtein, aMolecule);
            if (!checkAndReportErrorCodeForRecordType(recType, isRecLineOK, isFileRead)) {
                break;
            }
        }
        else if (recType == PDB_RECORD_TYPE_HETATM) {
            //// setHeteroAtomRecordsToMolecule : TO BE COMPLETED
            rg_FLAG isRecLineOK = setAtomRecordsToMoleculeForPDBQ(*currRecLine, mapOfAtom, mapsOfResidue, mapOfNonStdResidue, mapOfChain, mapOfAtomSymbol, mapOfResidueSymbolForProtein, aMolecule);
            if (!checkAndReportErrorCodeForRecordType(recType, isRecLineOK, isFileRead)) {
                break;
            }
        }
        else if (recType == PDB_RECORD_TYPE_CONECT) {
            connectRecLinesOfPDBFile.push_back(&(*currRecLine));
        }
        else if (recType == PDB_RECORD_TYPE_HELIX) {
            helixRecLinesOfPDBFile.push_back(currRecLine);
        }
        else if (recType == PDB_RECORD_TYPE_SHEET) {
            sheetRecLinesOfPDBFile.push_back(currRecLine);
        }
        else if (recType == PDB_RECORD_TYPE_TURN) {
            turnRecLinesOfPDBFile.push_back(currRecLine);
        }
        else if (recType == PDB_RECORD_TYPE_MODEL) {
            rg_INT modelSerial = atoi(StringFunctions::subString(*currRecLine, PDB_MODEL_RECORD_SERIAL_ST_POS, PDB_MODEL_RECORD_SERIAL_LENGTH).c_str());
            if (modelSerial != 0) {
                aMolecule->setModelSerialFromInputFile(modelSerial);
            }
        }
        else if (recType == PDB_RECORD_TYPE_ENDMDL) {
            countModelRec++;
            aMolecule = molecules.addTail(Molecule());
            aMolecule->setModelSerialFromInputFile(countModelRec);
            aMolecule->setTimeStamp(depDate);
            aMolecule->setFileSize(fileSize);

            mapOfAtom.clear();
            //mapOfResidue.clear();             //// WRONG SOURCE CODE !!! NOT WORKING WITH SEVERAL MODELS...
            //mapOfNonStdResidue.clear();       //// WRONG SOURCE CODE !!! NOT WORKING WITH SEVERAL MODELS...
            mapsOfResidue.removeAll();
            mapOfChain.clear();
        }
        else if (recType == PDB_RECORD_TYPE_EXPDTA) {
            string method = currRecLine->substr(10, currRecLine->length());
            method = StringFunctions::strTrimRight(method);
            aMolecule->setMethod(method);
        }
        else if (recType == PDB_RECORD_TYPE_REMARK) {
            aMolecule->addHeaderRecords(*currRecLine);
        }
        else if (recType == PDB_RECORD_TYPE_CRYST) {
            aMolecule->setCryst(true);

            string crystStr = currRecLine->substr(9, currRecLine->length());
            crystStr = StringFunctions::strTrimLeft(crystStr);

            const char* seps = " \t\n";
            char*	    token = NULL;

            char tempCryst[200];
            strcpy(tempCryst, crystStr.c_str());
            token = strtok(tempCryst, seps);

            for (int i = 0; i < 6; i++) {
                aMolecule->setCrystInfo(i, atof(token));
                token = strtok(NULL, seps);
            }
        }

        i_recLines++;
    }

    //// WRONG SOURCE CODE !!! NOT WORKING WITH SEVERAL MODELS...
    //

    //// POST PROCESSES: KILL EMPTY MOLECULE MODEL, SET CONNECT RECORDS, AND ETC...
    //

    if (isFileRead) {

        aMolecule = rg_NULL;
        molecules.reset4Loop();
        while (molecules.setNext4Loop()) {
            aMolecule = molecules.getpEntity();

            // KILL EMPTY MOLECULE MODEL
            if (aMolecule->getAtoms()->getSize() == 0) {
                molecules.killCurrent();
                continue;
            }


            // SET CONNECT RECORDS           
            list<string*>::iterator j_recLines = connectRecLinesOfPDBFile.begin();
            while (j_recLines != connectRecLinesOfPDBFile.end()) {
                string* currRecLine = *j_recLines;

                rg_FLAG isRecLineOK = setConectRecordsToMolecule(*currRecLine, mapOfAtom, aMolecule);
                if (!checkAndReportErrorCodeForRecordType(PDB_RECORD_TYPE_CONECT, isRecLineOK, isFileRead)) {
                    break;
                }
                j_recLines++;
            }
            setChemicalBondsToMoleculeForStandardResidues(*aMolecule);


            // SET CHAINS FOR NON STD RESIDUES            
            //setChainsForNonStdResiduesForPDBFile( &mapsOfNonStdResidue, *aMolecule );

            evaluateAndSetChainCode(*aMolecule);
            evaluateAndSetResidueCodesForNeucleicAcid(&mapOfResidueSymbolForDNA, &mapOfResidueSymbolForRNA, *aMolecule);


            // SET SECONDARY STRUCTURES
            // mapsOfResidue : HETATM is not considered for Chain.
            //setHelixRecordsToMolecule( &helixRecLinesOfPDBFile, &mapsOfResidue, aMolecule );
            //setSheetRecordsToMolecule( &sheetRecLinesOfPDBFile, &mapsOfResidue, aMolecule );
            //setTurnRecordsToMolecule( &turnRecLinesOfPDBFile, &mapsOfResidue, aMolecule );

            // ETC...

            aMolecule->evaluateHydrogenDonorAndAcceptorAtoms();
            aMolecule->computeAndSetCenterOfMass();
            aMolecule->computeAndSetMinEnclosingSphere();




        }



    }
    //
    ////

    return isFileRead;
}




rg_FLAG MoleculeIOFunctions::setAtomRecordsToMoleculeForPDBQ(const string& strRecLine, AtomMap& mapOfAtom, rg_dList<ResidueMapWithChainID>& mapsOfResidue, ResidueMap& mapOfNonStdResidue, ChainMap& mapOfChain, AtomSymbolMap& mapOfAtomSymbol, ResidueSymbolMap& mapOfResidueSymbol, Molecule* aMolecule)
{

    rg_INT atomSerial = atoi(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_SERIAL_ST_POS, PDB_ATOM_RECORD_SERIAL_LENGTH).c_str());
    string atomName = StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_NAME_ST_POS, PDB_ATOM_RECORD_NAME_LENGTH);
    string altLoc = StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_ALTLOC_ST_POS, PDB_ATOM_RECORD_ALTLOC_LENGTH);   // NOT USING
    string resName = StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_RESNAME_ST_POS, PDB_ATOM_RECORD_RESNAME_LENGTH);
    string chainID = StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_CHAINID_ST_POS, PDB_ATOM_RECORD_CHAINID_LENGTH);
    rg_INT resSeq = atoi(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_RESSEQ_ST_POS, PDB_ATOM_RECORD_RESSEQ_LENGTH).c_str());
    string iCode = StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_ICODE_ST_POS, PDB_ATOM_RECORD_ICODE_LENGTH);     // NOT USING
    float  x_coord = atof(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH).c_str());
    float  y_coord = atof(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH).c_str());
    float  z_coord = atof(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH).c_str());
    //rg_REAL  x_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_X_ST_POS, PDB_ATOM_RECORD_X_LENGTH ).c_str() );
    //rg_REAL  y_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Y_ST_POS, PDB_ATOM_RECORD_Y_LENGTH ).c_str() );
    //rg_REAL  z_coord    = atof( StringFunctions::subString( strRecLine, PDB_ATOM_RECORD_Z_ST_POS, PDB_ATOM_RECORD_Z_LENGTH ).c_str() );
    float  occupancy = atof(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_OCCUPANCY_ST_POS, PDB_ATOM_RECORD_OCCUPANCY_LENGTH).c_str());
    float  tempFactor = atof(StringFunctions::subString(strRecLine, PDB_ATOM_RECORD_TEMPFACTOR_ST_POS, PDB_ATOM_RECORD_TEMPFACTOR_LENGTH).c_str());
    float  partialCharge = atof(StringFunctions::subString(strRecLine, PDBQ_ATOM_RECORD_PARTIAL_CHARGE_ST_POS, PDBQ_ATOM_RECORD_PARTIAL_CHARGE_LENGTH).c_str());



    if (altLoc == " " || altLoc == "A") {

        // Estimate Chain
        Chain* currChain = rg_NULL;
        rg_INT numOfChains = mapOfChain.size();
        rg_INT intChainID = (int)chainID[0];


        ResidueMapWithChainID* pResidueMapWithChainID = rg_NULL;

        ChainMap::iterator chainMap_i = mapOfChain.find(intChainID);

        if (chainMap_i == mapOfChain.end()) {
            currChain = aMolecule->addChain(Chain(numOfChains, aMolecule));
            currChain->setChainIDFromInputFileInDecimal(intChainID);
            mapOfChain.insert(ChainMap::value_type(intChainID, currChain));

            ResidueMapWithChainID tempResidueMapWithChainID;
            tempResidueMapWithChainID.first = numOfChains;
            pResidueMapWithChainID = mapsOfResidue.addTail(tempResidueMapWithChainID);
        }
        else {
            currChain = (*chainMap_i).second;

            mapsOfResidue.reset4Loop();
            while (mapsOfResidue.setNext4Loop()) {
                pResidueMapWithChainID = mapsOfResidue.getpEntity();

                if (pResidueMapWithChainID->first == currChain->getID())
                    break;
            }
        }

        ResidueMap* aMapOfResidue = &(pResidueMapWithChainID->second);



        // Estimate Residue
        Residue* currResidue = rg_NULL;
        rg_INT   numOfResidue = 0;
        mapsOfResidue.reset4Loop();
        while (mapsOfResidue.setNext4Loop()) {
            numOfResidue += mapsOfResidue.getpEntity()->second.size();
        }

        ResidueMap::iterator residueMap_i = aMapOfResidue->find(resSeq);

        if (residueMap_i == aMapOfResidue->end()) {
            currResidue = aMolecule->addResidue(Residue(numOfResidue));
            currChain->addResidue(currResidue);
            currResidue->setChain(currChain);
            currResidue->setSequenceNumber(resSeq);
            currResidue->setResidueName(resName);
            setResidueCodeToTargetResidueExceptNeucleicAcid(mapOfResidueSymbol, resName, currResidue);  // NOT WORKING WITH DNA/RNA.
            aMapOfResidue->insert(ResidueMap::value_type(resSeq, currResidue));

            if (!currResidue->isStandardResidue()) {
                mapOfNonStdResidue.insert(ResidueMap::value_type(resSeq, currResidue));
            }
        }
        else {
            currResidue = (*residueMap_i).second;
        }

        //// Checking Oxygen in OXT row for duplication.
        //
        rg_BOOL isOXTDupplicated = rg_FALSE;

        if (atomName == " OXT") {
            rg_dList<Atom*>* atomsInResidue = currResidue->getAtoms();
            atomsInResidue->reset4Loop();
            while (atomsInResidue->setNext4Loop()) {
                Atom* currAtom = atomsInResidue->getEntity();

                if (currAtom->getAtomNameFromInputFile() == " O  ") {
                    isOXTDupplicated = rg_TRUE;
                    break;
                }
            }
        }

        if (isOXTDupplicated == rg_TRUE)
            return rg_TRUE;
        //
        ////


        // Create Atom : (1)Serial from inputfile, (2)AtomCode, (3)RemotenessIndicator, (4)BranchDesignagtor(include (5)ext.)
        //               (6)AtomTypeInAmber, (7)charge and (8)isAtomOnBackBone are considered.
        rg_INT numOfAtom = mapOfAtom.size();

        Atom* currAtom = aMolecule->addAtom(Atom(numOfAtom));
        currResidue->addAtom(currAtom);
        currAtom->setResidue(currResidue);
        currAtom->setSerialFromInputFile(atomSerial);     //  (1)
        currAtom->setAtomNameFromInputFile(atomName);
        rg_FLAG isAtomNameOK = setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile(mapOfAtomSymbol, atomName, currAtom); // (2)~(8)

                                                                                                                                        // Hydrogen atom code verification
        if (isAtomNameOK == rg_FALSE)
            return rg_FALSE;

        // set coordinates, radius, occupancy, and tempFactor.
        currAtom->setAtomBall(Sphere(x_coord, y_coord, z_coord, ATOM_FEATURES[currAtom->getAtomCode()].radius));
        currAtom->getpChemicalProperties()->setOccupancy(occupancy);
        currAtom->getpChemicalProperties()->setTempFactor(tempFactor);
        currAtom->getpChemicalProperties()->setCharge(partialCharge);

        mapOfAtom.insert(AtomMap::value_type(atomSerial, currAtom));
    }

    return rg_TRUE;
}



rg_FLAG MoleculeIOFunctions::writePDBQFile(const string& pathFile, const Molecule& aMolecule, const list<string>& frontRemarks)
{
    rg_FLAG isFileWritten = rg_FALSE;

    ofstream  fout(pathFile.c_str());

    writePDBQFile( fout, aMolecule, frontRemarks);

    fout.close();


    return isFileWritten;
}



void MoleculeIOFunctions::writePDBQFile(ofstream& fout, const Molecule& aMolecule, const list<string>& frontRemarks)
{
    string modelID = StringFunctions::convertIntegerToString(aMolecule.getModelSerialFromInputFile());
    string blankSpace = StringFunctions::getBlankStringForFixedField(PDB_MODEL_RECORD_SERIAL_LENGTH + 4, modelID);
    string modelRecord = PDB_RECORD_TYPE_MODEL + blankSpace + modelID;
    fout << modelRecord << endl;


    writeHeaderOfPDBQFile(fout, aMolecule, frontRemarks);

    fout << "ROOT" << endl;

    writeAtomsOfPDBQFile(fout, aMolecule);

    fout << "ENDROOT" << endl;

    fout << PDB_RECORD_TYPE_ENDMDL << endl;
}



void    MoleculeIOFunctions::writeHeaderOfPDBQFile(ofstream& fout, const Molecule& aMolecule, const list<string>& frontRemarks)
{
    for (list<string>::const_iterator i_str = frontRemarks.begin(); i_str != frontRemarks.end(); ++i_str) {
        fout << *i_str;
    }

    const rg_dList<string>&  headerRecords = aMolecule.getHeaderRecords();
    headerRecords.reset4Loop();
    while (headerRecords.setNext4Loop()) {
        fout << headerRecords.getEntity() << endl;
    }
}



void    MoleculeIOFunctions::writeAtomsOfPDBQFile(ofstream& fout, const Molecule& aMolecule)
{
    fout.setf(ios::fixed, ios::floatfield);     // TO 3 PLACES OF DECIMALS 
    fout.precision(3);                          //

    const rg_dList<Atom>& atoms = aMolecule.getAtomsInMolecule();

    Atom* currAtom = rg_NULL;
    atoms.reset4Loop();
    while (atoms.setNext4Loop()) {
        currAtom = atoms.getpEntity();

        // for SCP solver
        //if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
        //	continue;

        string recLine = PDB_RECORD_TYPE_ATOM;
        string blankSpace = "";
        string singleSpace = " ";

        string serial;
        string name;
        string altLoc;
        string resName;
        string chainID;
        string resSeq;
        string xCoord;
        string yCoord;
        string zCoord;
        string occupancy;
        string tempFactor;
        string partialCharge;


        // SERIAL FIELD
        serial = StringFunctions::convertIntegerToString(currAtom->getSerialFromInputFile());
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_SERIAL_LENGTH, serial);
        recLine = recLine + blankSpace + serial + singleSpace;


        // NAME FIELD
        if (!currAtom->getResidue()->isStandardResidue()) {
            name = currAtom->getAtomNameFromInputFile();
        }

        else {
            string extraBranch = currAtom->getpChemicalProperties()->getExtraBrangeDesignatorInString();
            string atomName = StringFunctions::strTrim((string)ATOM_FEATURES[(int)currAtom->getAtomCode()].symbol);
            string remoteInd = currAtom->getpChemicalProperties()->getRemoteIndicatorInString();
            string branchDes = currAtom->getpChemicalProperties()->getBrangeDesignatorInString();

            if (atomName.length() == 2)
                extraBranch = "";

            name = extraBranch + atomName + remoteInd + branchDes;
        }

        recLine = recLine + name;

        // ALTLOC FIELD --> CURRENTLY, ALTLOC IS NOT APPLIED.
        altLoc = " ";
        recLine = recLine + altLoc;


        // RESIDUE NAME
        resName = currAtom->getResidue()->getResidueName();
        recLine = recLine + resName + singleSpace;


        // CHAIN ID
        rg_INT chainIDInInt = currAtom->getResidue()->getChain()->getID() + 65;
        char   chainIDInChar = chainIDInInt;
        chainID = chainID + chainIDInChar;
        recLine = recLine + chainID;


        // RESIDUE SEQUENCE
        resSeq = StringFunctions::convertIntegerToString(currAtom->getResidue()->getSequenceNumber());
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_RESSEQ_LENGTH, resSeq);
        recLine = recLine + blankSpace + resSeq;


        // ICODE --> CURRENTLY, ICODE IS NOT APPLIED.
        recLine = recLine + "    ";


        // X COORDINATE
        xCoord = StringFunctions::convertDoubleToString(currAtom->getpAtomBall()->getCenter().getX(), 3);
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_X_LENGTH, xCoord);
        recLine = recLine + blankSpace + xCoord;


        // Y COORDINATE
        yCoord = StringFunctions::convertDoubleToString(currAtom->getpAtomBall()->getCenter().getY(), 3);
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_Y_LENGTH, yCoord);
        recLine = recLine + blankSpace + yCoord;


        // Y COORDINATE
        zCoord = StringFunctions::convertDoubleToString(currAtom->getpAtomBall()->getCenter().getZ(), 3);
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_Z_LENGTH, zCoord);
        recLine = recLine + blankSpace + zCoord;


        // OCCUPANCY
        occupancy = "";
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_OCCUPANCY_LENGTH, occupancy);
        recLine = recLine + blankSpace + occupancy;


        // TEMPFACTOR
        tempFactor = "";
        blankSpace = StringFunctions::getBlankStringForFixedField(PDB_ATOM_RECORD_TEMPFACTOR_LENGTH, tempFactor);
        recLine = recLine + blankSpace + tempFactor;

        // PARTIAL CHARGE.
        recLine = recLine + "    ";

        partialCharge = StringFunctions::convertDoubleToString(currAtom->getpChemicalProperties()->getCharge(), 3);
        blankSpace = StringFunctions::getBlankStringForFixedField(PDBQ_ATOM_RECORD_PARTIAL_CHARGE_LENGTH, partialCharge);
        recLine = recLine + blankSpace + partialCharge;


        fout << recLine << endl;
    }
}
