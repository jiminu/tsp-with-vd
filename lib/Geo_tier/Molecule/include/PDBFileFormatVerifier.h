#ifndef _PDBFILEFORMATVERIFIER_H_
#define _PDBFILEFORMATVERIFIER_H_

#include "rg_Const.h"
#include "rg_Molecule.h"

using namespace V::GeometryTier;



class PDBFileFormatVerifier
{
public:
	// Temporary code
	// This should be replaced with code which handles "REMARK 465 and 470".
    static rg_BOOL chceckWhetherPDBContainsMissingAtomsResiduesByLookingAtRemarkNumbers(const string& inputPDBFileNameWithPath,
                                                                                        rg_BOOL     & bHaveMissingAtoms,
                                                                                        rg_BOOL     & bHaveMissingResidues);
        static rg_BOOL haveMissingAtoms(const list<string> recordsOfPDBFile);
        static rg_BOOL haveMissingResidues(const list<string> recordsOfPDBFile);
    static rg_BOOL haveMissingAtoms(const string& inputPDBFileNameWithPath);
	static rg_BOOL haveMissingResidues(const string& inputPDBFileNameWithPath);
    static rg_BOOL haveHydrogenAtoms(const string& inputPDBFileNameWithPath);
	static rg_BOOL haveMissingResidues2(const string& inputPDBFileNameWithPath);
	static rg_INT  findMissingResidues(const char* inputPDBFileNameWithPath, rg_dList<ChainIDWithSequenceNumbersOfMissingResidues>& chainIDsWithSeqNumbersOfMissingResidues);
    static rg_INT  findMissingResidues(const V::GeometryTier::Molecule& protein, rg_dList<ChainIDWithSequenceNumbersOfMissingResidues>& chainIDsWithSeqNumbersOfMissingResidues);
			
	static rg_BOOL haveValidSidechainAtoms(V::GeometryTier::Molecule& protein, rg_INT& seqNumOfFirstOccurredInvalidResidue);

	static rg_BOOL compareTwoProteinsForSCP(V::GeometryTier::Molecule& protein1, V::GeometryTier::Molecule& protein2);

	static rg_BOOL compareSequencesNResidueAtomsOfTwoProteinsForSCP(V::GeometryTier::Molecule& protein1,
                                                                    V::GeometryTier::Molecule& protein2,
																	rg_dList<V::GeometryTier::Residue*>& residuesWithDiscrepancy1,
																	rg_dList<V::GeometryTier::Residue*>& residuesWithDiscrepancy2);

private:
	static rg_BOOL getRecordsOfPDBFile(const string& inputPDBFileNameWithPath, list<string>& records);
};

#endif