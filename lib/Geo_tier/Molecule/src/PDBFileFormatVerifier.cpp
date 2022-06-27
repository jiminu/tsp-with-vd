#include "PDBFileFormatVerifier.h"
#include "ConstForMoleculeIOFunctions.h"
#include "FunctionsForRotamerLibrary.h"
#include "FunctionsForMolecule.h"
#include "MoleculeIOFunctions.h"
#include "StringFunctions.h"
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

rg_BOOL PDBFileFormatVerifier::chceckWhetherPDBContainsMissingAtomsResiduesByLookingAtRemarkNumbers(const string& inputPDBFileNameWithPath,
                                                                                                    rg_BOOL     & bHaveMissingAtoms,
                                                                                                    rg_BOOL     & bHaveMissingResidues)
{
    // get records from PDB file
    list<string> records;
    rg_BOOL bValidFile = rg_FALSE;
    bValidFile = getRecordsOfPDBFile(inputPDBFileNameWithPath, records);

    rg_BOOL bHaveMissingAtomsResidues = rg_FALSE;
    if(! bValidFile)
    {
        bHaveMissingAtomsResidues = rg_TRUE;
        bHaveMissingAtoms = rg_TRUE;
        bHaveMissingResidues = rg_TRUE;
        return bHaveMissingAtomsResidues;
    }

    bHaveMissingAtoms = rg_FALSE;
    bHaveMissingAtoms = haveMissingAtoms(records);

    bHaveMissingResidues = rg_FALSE;
    bHaveMissingResidues = haveMissingResidues(records);

    if(bHaveMissingAtoms && bHaveMissingResidues)
        bHaveMissingAtomsResidues = rg_TRUE;

    return bHaveMissingAtomsResidues;
}

rg_BOOL PDBFileFormatVerifier::haveMissingAtoms(const list<string> recordsOfPDBFile)
{
    rg_BOOL bHaveMissingAtoms = rg_FALSE;

    list<string>::const_iterator itr = recordsOfPDBFile.begin();
    while( itr != recordsOfPDBFile.end() ) 
    {
        string currRecord = (*itr);
        string recType     = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );

        if( recType == PDB_RECORD_TYPE_REMARK ) 
        {
            string remarkNumStr = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_LENGTH + 1, PDB_RECORD_TYPE_LENGTH + 3 );
            rg_INT remarkNum    = atoi(remarkNumStr.c_str());

            if(remarkNum == 470)
            {
                bHaveMissingAtoms = rg_TRUE;
                break;
            }			
        }
        itr++;
    }

    return bHaveMissingAtoms;
}

rg_BOOL PDBFileFormatVerifier::haveMissingResidues(const list<string> recordsOfPDBFile)
{
    rg_BOOL bHaveMissingResidues = rg_FALSE;

    list<string>::const_iterator itr = recordsOfPDBFile.begin();
    while( itr != recordsOfPDBFile.end() ) 
    {
        string currRecord = (*itr);
        string recType     = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );

        if( recType == PDB_RECORD_TYPE_REMARK ) 
        {
            string remarkNumStr = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_LENGTH + 1, PDB_RECORD_TYPE_LENGTH + 3 );
            rg_INT remarkNum    = atoi(remarkNumStr.c_str());

            if(remarkNum == 465)
            {
                bHaveMissingResidues = rg_TRUE;
                break;
            }			
        }
        itr++;
    }

    return bHaveMissingResidues;
}

rg_BOOL PDBFileFormatVerifier::haveMissingAtoms(const string& inputPDBFileNameWithPath)
{
    // get records from PDB file
    list<string> records;
    rg_BOOL bValidFile = rg_FALSE;
    bValidFile = getRecordsOfPDBFile(inputPDBFileNameWithPath, records);

    rg_BOOL hHaveMissingAtoms    = rg_FALSE;
    if(! bValidFile)
    {
        hHaveMissingAtoms = rg_TRUE;
        return hHaveMissingAtoms;
    }

    list<string>::const_iterator itr = records.begin();
    while( itr != records.end() ) 
    {
        string currRecord = (*itr);
        string recType     = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );

        if( recType == PDB_RECORD_TYPE_REMARK ) 
        {
            string remarkNumStr = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_LENGTH + 1, PDB_RECORD_TYPE_LENGTH + 3 );
            rg_INT remarkNum    = atoi(remarkNumStr.c_str());

            if(remarkNum == 470)
            {
                hHaveMissingAtoms = rg_TRUE;
                break;
            }			
        }
        itr++;
    }

    return hHaveMissingAtoms;
}

rg_BOOL PDBFileFormatVerifier::haveMissingResidues(const string& inputPDBFileNameWithPath)
{
	// get records from PDB file
	list<string> records;
	rg_BOOL bValidFile = rg_FALSE;
	bValidFile = getRecordsOfPDBFile(inputPDBFileNameWithPath, records);

	// check if there exist the missing residues
	rg_BOOL bHaveMissingResidues = rg_FALSE;
	if(! bValidFile)
	{
		bHaveMissingResidues = rg_TRUE;
		return bHaveMissingResidues;
	}

	list<string>::const_iterator itr = records.begin();
	while( itr != records.end() ) 
	{
		string currRecord = (*itr);
		string recType     = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_ST_POS, PDB_RECORD_TYPE_LENGTH );

		if( recType == PDB_RECORD_TYPE_REMARK ) 
		{
			string remarkNumStr = StringFunctions::subString( currRecord, PDB_RECORD_TYPE_LENGTH + 1, PDB_RECORD_TYPE_LENGTH + 3 );
			rg_INT remarkNum    = atoi(remarkNumStr.c_str());

			if(remarkNum == 465)
			{
				bHaveMissingResidues = rg_TRUE;
				break;
			}			
		}
		itr++;
	}

	return bHaveMissingResidues;
}


rg_BOOL PDBFileFormatVerifier::haveHydrogenAtoms(const string& inputPDBFileNameWithPath)
{
    // read protein structure from a PDB file
    const rg_INT targetModelID = 1;
    Molecule protein;
    MoleculeIOFunctions::readPDBFile(targetModelID, 
                                     inputPDBFileNameWithPath.c_str(), 
                                     protein);

    rg_BOOL bHydrogenAtom = rg_FALSE;
    rg_dList<Atom>* atoms = protein.getAtoms();
    atoms->reset4Loop();
    while (atoms->setNext4Loop())
    {
        Atom* currAtom = atoms->getpEntity();
        if(currAtom->getAtomCode() == H_ATOM)
        {
            bHydrogenAtom = rg_TRUE;
            break;
        }
    }

    return bHydrogenAtom;
}


rg_BOOL PDBFileFormatVerifier::haveMissingResidues2(const string& inputPDBFileNameWithPath)
{
	// read protein structure from a PDB file
	const rg_INT targetModelID = 1;
	Molecule protein;
	MoleculeIOFunctions::readPDBFile(targetModelID, 
		                             inputPDBFileNameWithPath.c_str(), 
		                             protein);

	rg_BOOL bHaveMissingResidues = rg_FALSE;

	if(protein.haveMissingResidues())
		bHaveMissingResidues = rg_TRUE;

	return bHaveMissingResidues;
}

// Missing residues should be recognized by parsing the record "REMARK" in PDB file.
// However, currently we do not parse "REMARK".
// Hence, this function just check the sequence of residue numbers in order to find out the missing residues.
rg_INT  PDBFileFormatVerifier::findMissingResidues(const char* inputPDBFileNameWithPath, rg_dList<ChainIDWithSequenceNumbersOfMissingResidues>& chainIDsWithSeqNumbersOfMissingResidues)
{
	// read protein structure from a PDB file
	const rg_INT targetModelID = 1;
	Molecule protein;
	rg_FLAG bRead = rg_FALSE;
	bRead = MoleculeIOFunctions::readPDBFile(targetModelID, 
		                                     inputPDBFileNameWithPath, 
		                                     protein);	

	rg_INT numMissingResidues = 0;
	if(bRead)
	{
		protein.findMissingResidues(chainIDsWithSeqNumbersOfMissingResidues);

		chainIDsWithSeqNumbersOfMissingResidues.reset4Loop();
		while (chainIDsWithSeqNumbersOfMissingResidues.setNext4Loop())
		{
			ChainIDWithSequenceNumbersOfMissingResidues* currChainIDWithSequenceNumbersOfMissingResidues = rg_NULL;
			currChainIDWithSequenceNumbersOfMissingResidues = chainIDsWithSeqNumbersOfMissingResidues.getpEntity();
			rg_INT currNumMissingResidues = 0;
			currNumMissingResidues = currChainIDWithSequenceNumbersOfMissingResidues->getSeqNumbersOfMissingResidues()->getSize();
			numMissingResidues += currNumMissingResidues;
		}
	}
	
	return numMissingResidues;
}

rg_INT  PDBFileFormatVerifier::findMissingResidues(const Molecule& protein, rg_dList<ChainIDWithSequenceNumbersOfMissingResidues>& chainIDsWithSeqNumbersOfMissingResidues)
{
    rg_INT numMissingResidues = 0;
    protein.findMissingResidues(chainIDsWithSeqNumbersOfMissingResidues);

    chainIDsWithSeqNumbersOfMissingResidues.reset4Loop();
    while (chainIDsWithSeqNumbersOfMissingResidues.setNext4Loop())
    {
        ChainIDWithSequenceNumbersOfMissingResidues* currChainIDWithSequenceNumbersOfMissingResidues = rg_NULL;
        currChainIDWithSequenceNumbersOfMissingResidues = chainIDsWithSeqNumbersOfMissingResidues.getpEntity();
        rg_INT currNumMissingResidues = 0;
        currNumMissingResidues = currChainIDWithSequenceNumbersOfMissingResidues->getSeqNumbersOfMissingResidues()->getSize();
        numMissingResidues += currNumMissingResidues;
    }

    return numMissingResidues;
}


rg_BOOL PDBFileFormatVerifier::haveValidSidechainAtoms(Molecule& protein, rg_INT& seqNumOfFirstOccurredInvalidResidue)
{
	rg_BOOL bHaveValidSidechainAtoms = rg_TRUE;

	rg_dList<Residue*> orderedResidues;
	protein.getAminoResiduesOrderedByChainIDSequenceNumber( orderedResidues );

	orderedResidues.reset4Loop();
	while (orderedResidues.setNext4Loop())
	{
		Residue* currResidue = orderedResidues.getEntity();
		ResidueCode code = currResidue->getResidueCode();
		if(code == GLY_AMINO_RESIDUE)
			continue;

		RemoteIndicator startRemote = UNK_REMOTE;
		RemoteIndicator endRemote   = UNK_REMOTE;
		rg_BOOL bExist = rg_TRUE;
		bExist = FunctionsForRotamerLibrary
			::getStartNEndRemoteIndicatorsOfSidechainForResidue(code, startRemote, endRemote);


		rg_dList<Atom*> atomsOnSidechain;
		currResidue->getAtomsOnSideChain(&atomsOnSidechain);

		RemoteIndicator minRemote = TERMINATE_REMOTE;
		RemoteIndicator maxRemote = BETA_REMOTE;

		atomsOnSidechain.reset4Loop();
		while (atomsOnSidechain.setNext4Loop())
		{
			Atom* currAtom = atomsOnSidechain.getEntity();
			RemoteIndicator currRemote = UNK_REMOTE;
			currRemote = currAtom->getpChemicalProperties()->getRemoteIndicator();
			
			if(currRemote < minRemote)
				minRemote = currRemote;
			if(currRemote > maxRemote)
				maxRemote = currRemote;
		}

		if(minRemote != startRemote || maxRemote != endRemote)
		{
			bHaveValidSidechainAtoms = rg_FALSE;
			seqNumOfFirstOccurredInvalidResidue = getChainIDSeqNumOfResidue(currResidue);
			break;
		}
	}

	return bHaveValidSidechainAtoms;
}

rg_BOOL PDBFileFormatVerifier::compareTwoProteinsForSCP(Molecule& protein1, Molecule& protein2)
{
	protein1.moveOxygenInCarboxylGroupToBackbone();
	protein2.moveOxygenInCarboxylGroupToBackbone();

	// compare backbone atoms
	Atom** backboneAtoms1 = rg_NULL;
	rg_INT numBackboneAtoms1
		= protein1.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms1 );

	Atom** backboneAtoms2 = rg_NULL;
	rg_INT numBackboneAtoms2
		= protein2.getAtomsOnBackboneSortedByChainIDSequenceNumber( backboneAtoms2 );

	if(numBackboneAtoms1 != numBackboneAtoms2)
		return rg_FALSE;

	rg_INT i;
	for (i = 0;i < numBackboneAtoms1;i++)
	{
		rg_INT id1 = backboneAtoms1[ i ]->getID();
		rg_INT id2 = backboneAtoms2[ i ]->getID();
		AtomCode code1 = backboneAtoms1[ i ]->getAtomCode();
		AtomCode code2 = backboneAtoms2[ i ]->getAtomCode();
		Sphere ball1 = backboneAtoms1[ i ]->getAtomBall();
		Sphere ball2 = backboneAtoms2[ i ]->getAtomBall();

		rg_REAL dist = (ball1.getCenter()-ball2.getCenter()).magnitude();
		if(id1 != id2 || code1 != code2 || dist > 0.01 )
		{
			//return rg_FALSE;
			rg_BOOL bEqual = rg_FALSE;
		}
	}

	if(backboneAtoms1 != rg_NULL)
		delete[] backboneAtoms1;

	if(backboneAtoms2 != rg_NULL)
		delete[] backboneAtoms2;

	// compare side chains
	Residue** sortedResidues1 = rg_NULL;
	rg_INT numResidues1 = protein1.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues1);

	Residue** sortedResidues2 = rg_NULL;
	rg_INT numResidues2 = protein2.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues2);

	if(numResidues1 != numResidues2)
		return rg_FALSE;

	rg_INT j;
	for (j = 0;j < numResidues1;j++)
	{
		rg_INT id1 = sortedResidues1[ j ]->getID();
		rg_INT id2 = sortedResidues2[ j ]->getID();
		ResidueCode code1 = sortedResidues1[ j ]->getResidueCode();
		ResidueCode code2 = sortedResidues2[ j ]->getResidueCode();
		rg_INT seqNum1 = sortedResidues1[ j ]->getSequenceNumber();
		rg_INT seqNum2 = sortedResidues2[ j ]->getSequenceNumber();
		
		//rg_dList<Atom*>* atomSet1 = sortedResidues1[ j ]->getAtoms();
		//rg_dList<Atom*>* atomSet2 = sortedResidues2[ j ]->getAtoms();

		rg_dList<Atom*> atomSet1;
		rg_dList<Atom*> atomSet2;
		FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedResidues1[ j ], atomSet1);
		FunctionsForRotamerLibrary::getAtomsOnSidechain(sortedResidues2[ j ], atomSet2);
		
		rg_INT numAtoms1 = atomSet1.getSize();
		rg_INT numAtoms2 = atomSet2.getSize();

		if(id1 != id2 || code1 != code2 || seqNum1 != seqNum2 || numAtoms1 != numAtoms2)
			return rg_FALSE;

		atomSet1.reset4Loop();
		atomSet2.reset4Loop();
		while(atomSet1.setNext4Loop() && atomSet2.setNext4Loop())
		{
			Atom* currAtom1 = atomSet1.getEntity();
			Atom* currAtom2 = atomSet2.getEntity();

			rg_INT id1 = currAtom1->getID();
			rg_INT id2 = currAtom2->getID();
			AtomCode code1 = currAtom1->getAtomCode();
			AtomCode code2 = currAtom2->getAtomCode();
			Sphere ball1 = currAtom1->getAtomBall();
			Sphere ball2 = currAtom2->getAtomBall();

			rg_REAL dist = (ball1.getCenter()-ball2.getCenter()).magnitude();
			if(id1 != id2 || code1 != code2 || dist > 0.01 )
			{
				return rg_FALSE;
			}
		}
	}

	if(sortedResidues1 != rg_NULL)
		delete [] sortedResidues1;

	if(sortedResidues2 != rg_NULL)
		delete [] sortedResidues2;

	return rg_TRUE;
}

rg_BOOL PDBFileFormatVerifier::compareSequencesNResidueAtomsOfTwoProteinsForSCP(Molecule& protein1, 
	                                                                            Molecule& protein2,
																				rg_dList<Residue*>& residuesWithDiscrepancy1,
																				rg_dList<Residue*>& residuesWithDiscrepancy2)
{
	protein1.moveOxygenInCarboxylGroupToBackbone();
	protein2.moveOxygenInCarboxylGroupToBackbone();

	// compare sequence and residue atoms
	rg_dList<Residue*> sortedResidues1;
	protein1.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues1);

    //sortedResidues1.reset4Loop();
    //while(sortedResidues1.setNext4Loop())
    //{
    //    Residue* currResidue1 = sortedResidues1.getEntity();
    //    rg_dList<Atom*> atomsOfResidue1;
    //    currResidue1->getAtoms(atomsOfResidue1);
    //    atomsOfResidue1.reset4Loop();
    //    while (atomsOfResidue1.setNext4Loop())
    //    {
    //        Atom* currAtom = atomsOfResidue1.getEntity();
    //        if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
    //            atomsOfResidue1.killCurrent();
    //        if(currAtom->getAtomCode() == H_ATOM)
    //            atomsOfResidue1.killCurrent();
    //    }
    //}

	rg_dList<Residue*> sortedResidues2;
	protein2.getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues2);

    //sortedResidues2.reset4Loop();
    //while(sortedResidues2.setNext4Loop())
    //{
    //    Residue* currResidue2 = sortedResidues2.getEntity();
    //    rg_dList<Atom*> atomsOfResidue2;
    //    currResidue2->getAtoms(atomsOfResidue2);
    //    atomsOfResidue2.reset4Loop();
    //    while (atomsOfResidue2.setNext4Loop())
    //    {
    //        Atom* currAtom = atomsOfResidue2.getEntity();
    //        if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
    //            atomsOfResidue2.killCurrent();
    //        if(currAtom->getAtomCode() == H_ATOM)
    //            atomsOfResidue2.killCurrent();
    //    }
    //}


	sortedResidues1.reset4Loop();
	sortedResidues2.reset4Loop();
	while(sortedResidues1.setNext4Loop() && sortedResidues2.setNext4Loop())
	{
		Residue* currResidue1 = sortedResidues1.getEntity();
		Residue* currResidue2 = sortedResidues2.getEntity();

		ResidueCode code1 = currResidue1->getResidueCode();
		ResidueCode code2 = currResidue2->getResidueCode();

		if((rg_INT)code1 != (rg_INT)code2)
		{
			rg_INT seqNumberWithDiscrepancy1 = getChainIDSeqNumOfResidue( currResidue1 );
			rg_INT seqNumberWithDiscrepancy2 = getChainIDSeqNumOfResidue( currResidue2 );
			residuesWithDiscrepancy1.add(currResidue1);
			residuesWithDiscrepancy2.add(currResidue2);
		}
		else
		{
			rg_dList<Atom*> atomsOfResidue1;
			currResidue1->getAtoms(atomsOfResidue1);
            atomsOfResidue1.reset4Loop();
            while (atomsOfResidue1.setNext4Loop())
            {
                Atom* currAtom = atomsOfResidue1.getEntity();
                if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
                    atomsOfResidue1.killCurrent();
                if(currAtom->getAtomCode() == H_ATOM)
                    atomsOfResidue1.killCurrent();
            }

			rg_dList<Atom*> atomsOfResidue2;
			currResidue2->getAtoms(atomsOfResidue2);
            atomsOfResidue2.reset4Loop();
            while (atomsOfResidue2.setNext4Loop())
            {
                Atom* currAtom = atomsOfResidue2.getEntity();
                if(currAtom->getResidue()->getResidueCode() == HOH_RESIDUE)
                    atomsOfResidue2.killCurrent();
                if(currAtom->getAtomCode() == H_ATOM)
                    atomsOfResidue2.killCurrent();
            }

			rg_INT numAtoms1 = atomsOfResidue1.getSize();
			rg_INT numAtoms2 = atomsOfResidue2.getSize();

			if(numAtoms1 != numAtoms2)
			{
				rg_INT seqNumberWithDiscrepancy1 = getChainIDSeqNumOfResidue( currResidue1 );
				rg_INT seqNumberWithDiscrepancy2 = getChainIDSeqNumOfResidue( currResidue2 );
				residuesWithDiscrepancy1.add(currResidue1);
				residuesWithDiscrepancy2.add(currResidue2);
			}
			else
			{
				rg_dList<string> atomNamesOfResidue1;
				//getAtomNamesOfResidue(currResidue1, atomNamesOfResidue1);
                getAtomNames(atomsOfResidue1, atomNamesOfResidue1);
				rg_dList<string> atomNamesOfResidue2;
				//getAtomNamesOfResidue(currResidue2, atomNamesOfResidue2);
                getAtomNames(atomsOfResidue2, atomNamesOfResidue2);

				atomNamesOfResidue1.reset4Loop();
				while (atomNamesOfResidue1.setNext4Loop())
				{
					string currAtomName = atomNamesOfResidue1.getEntity();
					if( !atomNamesOfResidue2.isInList(currAtomName) )
					{
						residuesWithDiscrepancy1.add(currResidue1);
						residuesWithDiscrepancy2.add(currResidue2);
						break;
					}
				}

				atomNamesOfResidue2.reset4Loop();
				while (atomNamesOfResidue2.setNext4Loop())
				{
					string currAtomName = atomNamesOfResidue2.getEntity();
					if( !atomNamesOfResidue1.isInList(currAtomName) )
					{
						residuesWithDiscrepancy1.add(currResidue1);
						residuesWithDiscrepancy2.add(currResidue2);
						break;
					}
				}
			}
		}
	}

	if(residuesWithDiscrepancy1.getSize() == 0 && residuesWithDiscrepancy2.getSize() == 0)
		return rg_TRUE; // identical structure
	else
		return rg_FALSE;
}


rg_BOOL PDBFileFormatVerifier::getRecordsOfPDBFile(const string& inputPDBFileNameWithPath, list<string>& records)
{
	rg_BOOL bValidFile = rg_TRUE;

	ifstream fin( inputPDBFileNameWithPath.c_str() );

	if( fin.bad() ) 
	{
		cerr << "Error: Could not open file!\n";
		bValidFile = rg_FALSE;
		
		return bValidFile;
	}


	fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE

	// read each line of PDB file
	string strRecLine("");	
	while ( getline( fin, strRecLine ) ) 
	{
		records.push_back( strRecLine );
	}

	fin.close();

	return bValidFile;
}