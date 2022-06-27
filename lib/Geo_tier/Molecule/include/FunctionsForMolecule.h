#ifndef _FUNCTIONSFORMOLECULE_H
#define _FUNCTIONSFORMOLECULE_H

#include "rg_dList.h"
#include "rg_Const.h"
#include "rg_Molecule.h"
#include "BackboneDihedralAnglePair.h"
#include "SidechainDihedralAngleSet.h"
#include "ChainIDSeqNumRotLibIDPair.h"
//using namespace V::GeometryTier;

#include <vector>
using namespace std;

class RotamerSetOfResidue;



namespace V {
namespace GeometryTier {


class Atom;
class Residue;



rg_INT  compareResidueBySequenceNumber( const void* residueA, const void* residueB );

rg_INT  compareResidueByChainIDSequenceNumber( const void* residueA, const void* residueB );

rg_BOOL isThereBondBetween(Atom& atom1, Atom& atom2);

void    getAtomsOnBackboneCorrToResidues( rg_dList<Residue*>& targetListOfResidues, rg_dList<Atom*>& listOfAtomsOnBackbone );

void    getAtomsOnBackboneCorrTo18ResiduesExcludingGLYNALASortedBySequenceNumber(Molecule& molecule, rg_dList<Atom*>& listOfAtomsOnBackbone );

void    get18ResiduesExcludingGLYNALASortedBySequenceNumber(Molecule& molecule, rg_dList<Residue*>& targetListOfResidues );

void    get19ResiduesExcludingGLYSortedBySequenceNumber(Molecule& molecule, rg_dList<Residue*>& targetListOfResidues );

rg_INT  deleteAtomsInSideChainsCorrTo18ResiduesExcludingGLYNALA(Molecule& molecule, rg_BOOL replaceResidueCodeToUnknown, const RemoteIndicator& startRemoteIndicator);

rg_INT  recoverMissingAtomsIn18ResiduesExcludingGLYNALAWithoutCoordinates(Molecule& molecule, rg_dList<Residue*>& recoveredResidues);

rg_INT  recoverMissingAtomsIn18ResiduesExcludingGLYNALAWithoutCoordinates(Molecule& molecule);

rg_INT  getAtomsOfResidueWithGivenRemoteIndicatorNBranchDesignator(Residue* targetResidue, rg_dList<Atom*>& atomsOfResidueUnnecessaryAtomsRemoved, const RemoteIndicator& endRemote, const BranchDesignator& maxBranchDesignator);

rg_INT  getAtomsOnSidechainWithGivenRemoteIndicatorNBranchDesignator(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChainUnnecessaryAtomsRemoved, const RemoteIndicator& startRemote, const RemoteIndicator& endRemote, const BranchDesignator& maxBranchDesignator);

rg_INT  getAtomsOnSideChainInBreathFirstOrder(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChainInBreathFirstOrder, const RemoteIndicator& startRemote);

// currently not called from anywhere
rg_INT  getAtomsOnSideChainInBreathFirstOrder(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChainInBreathFirstOrder, const RemoteIndicator& startRemote, const RemoteIndicator& endRemote, const BranchDesignator& maxBranchDesignator);

rg_INT  getAtomsOnSideChainInBreathFirstOrder(rg_dList<Atom*>& atomsOnSideChain, rg_dList<Atom*>& atomsOnSideChainInBreathFirstOrder);

rg_INT  compareAtomsByRemoteIndicatorBranchDesignator(Atom* atom1, Atom* atom2);

rg_INT  getAtomsOnSideChainInBreathFirstOrder(Residue* targetResidue, rg_dList<Atom*>& atomsOnSideChain);

rg_BOOL getAtomsSortedByBranchDesignator(const RemoteIndicator& remote, Residue* targetResidue, rg_dList<Atom*>& atomsWithThisRemote);

Atom*   getNextAtomOnSideChainInBreathFirstOrder(Atom* targetAtom);

Atom*   getFirstAtomOnSideChainWithThisRemoteIndicatorInBreathFirstOrder(const RemoteIndicator& remote, Residue* targetResidue);

BranchDesignator getNextBranchDesignator(const BranchDesignator& target);
RemoteIndicator  getNextRemoteIndicator(const RemoteIndicator& target);

void    getAtomNamesOfResidue(const ResidueCode& code, rg_dList<string>& atomNames);

void    getAtomNamesOfResidueExceptHydrogenNOXT(const ResidueCode& code, rg_dList<string>& atomNames);

void    getAtomNamesOfSidechain(const ResidueCode& code, rg_dList<string>& atomNames);

void    getAtomNamesOfSidechainExceptHydrogenNOXT(const ResidueCode& code, rg_dList<string>& atomNames);

void    getAtomNamesOfResidue(Residue* residue, rg_dList<string>& atomNames);

void    getAtomNames(rg_dList<Atom*>& atoms, rg_dList<string>& atomNames);

rg_BOOL doesStringContainHydrogenAtomName(const string& atomName);
rg_BOOL doesStringContainHydrogenAtomName(const string& atomName, Atom* targetAtom);

rg_INT  getChainIDSeqNumOfResidue(Residue* residue);

rg_INT  getChainIDSeqNumOfResidue(const rg_INT& chainID, const rg_INT& seqNum);

void    getChainIDNSeqNumFromChainSeqID(const rg_INT& chainSeqID, rg_INT& chainID, rg_INT& seqNum);

void    getResidueCorrToChainIDSeqNum(Molecule* protein, const rg_INT& chainIDSeqNum, Residue*& targetResidue);

void    sortBackboneDihedralAnglePairs(rg_dList<BackboneDihedralAnglePair>& backboneDihedralAnglePairs);

bool    backboneDihedralAnglePairSortCriterion(const BackboneDihedralAnglePair& bdp1, const BackboneDihedralAnglePair& bdp2);

void    sortSidechainDihedralAngleSets(rg_dList<SidechainDihedralAngleSet>& sidechainDihedralAngleSets);

bool    sidechainDihedralAngleSetSortCriterion(const SidechainDihedralAngleSet& sds1, const SidechainDihedralAngleSet& sds2);

void    sortChainIDSeqNumRotLibIDPairSet(rg_dList<ChainIDSeqNumRotLibIDPair>& chainIDSeqNumRotLibIDPairSet);

bool    chainIDSeqNumRotLibIDPairSortCriterion(const ChainIDSeqNumRotLibIDPair& csr1, const ChainIDSeqNumRotLibIDPair& csr2);

rg_REAL computeVDWEnergyBetweenTwoAtoms(Atom* firstAtom, Atom* secondAtom);

rg_REAL computeVDWESHBEnergyOfProtein(Molecule& protein);

rg_REAL computeVDWESHBEnergyBetweenSidechainsNBackbone(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeVDWESHBEnergyBetweenSidechainPairs(Residue** residues, const rg_INT& numResidues);

rg_REAL computeVDWESHBEnergyOfComputedProteinStructureForSCPWithMissingSideChainAtoms(Molecule& targetProtein, Molecule& referenceProtein, rg_REAL& VDWEnergy, rg_REAL& electrostaticEnergy, rg_REAL& hydrogenBondingEnergy);

rg_REAL computeVDWESHBEnergyBetweenSidechainsNBackbone(Residue** sortedTargetResidues, 
                                                       const rg_INT& numTargetResidues, 
													   Residue** sortedReferenceResidues,
													   const rg_INT& numReferenceResidues,
                                                       Atom** backboneAtoms, 
                                                       const rg_INT& numBackboneAtoms,
                                                       rg_REAL& VDWEnergyBetweenSidechainsNBackbone,
                                                       rg_REAL& electrostaticEnergyBetweenSidechainsNBackbone,
                                                       rg_REAL& hydrogenBondingEnergyBetweenSidechainsNBackbone);

rg_REAL computeVDWESHBEnergyBetweenSidechainPairs(Residue** sortedTargetResidues, 
                                                  const rg_INT& numTargetResidues,
												  Residue** sortedReferenceResidues,
												  const rg_INT& numReferenceResidues,
                                                  rg_REAL& VDWEnergyBetweenSidechainPairs,
                                                  rg_REAL& electrostaticEnergyBetweenSidechainPairs,
                                                  rg_REAL& hydrogenBondingEnergyBetweenSidechainPairs);

rg_REAL computeVDWESHBEnergyOfProtein(Molecule& protein, rg_REAL& VDWEnergy, rg_REAL& electrostaticEnergy, rg_REAL& hydrogenBondingEnergy);

rg_REAL computeVDWESHBEnergyBetweenSidechainsNBackbone(Residue** residues, 
                                                       const rg_INT& numResidues, 
                                                       Atom** backboneAtoms, 
                                                       const rg_INT& numBackboneAtoms,
                                                       rg_REAL& VDWEnergyBetweenSidechainsNBackbone,
                                                       rg_REAL& electrostaticEnergyBetweenSidechainsNBackbone,
                                                       rg_REAL& hydrogenBondingEnergyBetweenSidechainsNBackbone);

rg_REAL computeVDWESHBEnergyBetweenSidechainPairs(Residue** residues, 
                                                  const rg_INT& numResidues,
                                                  rg_REAL& VDWEnergyBetweenSidechainPairs,
                                                  rg_REAL& electrostaticEnergyBetweenSidechainPairs,
                                                  rg_REAL& hydrogenBondingEnergyBetweenSidechainPairs);

rg_INT computeVDWESHBEnergyOfResiduesWithBackboneNOtherSidechains(Molecule& protein,
                                                                  rg_REAL*& engergyWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                  rg_REAL*& engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                  rg_REAL& energyOfProtein);

    rg_REAL computeVDWESHBEnergyBetweenSidechainsNBackbone(Residue** residues, 
                                                           const rg_INT& numResidues, 
                                                           Atom** backboneAtoms, 
                                                           const rg_INT& numBackboneAtoms,
                                                           rg_REAL*& engergyWithBackboneOfResidue_orderedByChainIDSeqNum);

    rg_REAL computeVDWESHBEnergyBetweenSidechainPairs(Residue** residues, 
                                                      const rg_INT& numResidues,
                                                      rg_REAL*& engergyWithOtherSidechainsOfResidue_orderedByChainIDSeqNum);

rg_INT computeXVolumeOfResidueWithBackboneNOtherSidechains(Molecule& protein,
                                                           rg_REAL*& xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                           rg_REAL*& xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                           rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                           rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                           rg_REAL& xvolumeOfProtein,
                                                           rg_INT&  numXOfProtein);
    
    rg_REAL computeXVolumeBetweenSidechainsNBackbone(Residue** residues, 
                                                     const rg_INT& numResidues, 
                                                     Atom** backboneAtoms, 
                                                     const rg_INT& numBackboneAtoms,
                                                     rg_REAL*& xvolumeWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                     rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                     rg_INT& numXBetweenSidechainsNBackbone);

    rg_REAL computeXVolumeBetweenSidechainPairs(Residue** residues, 
                                                const rg_INT& numResidues,
                                                rg_REAL*& xvolumeWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                rg_INT& numXBetweenSidechainPairs);


rg_INT computeXAreaOfResidueWithBackboneNOtherSidechains(Molecule& protein,
                                                         rg_REAL*& xAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_REAL*& xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                         rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                         rg_REAL& xAreaOfProtein,
                                                         rg_INT&  numXOfProtein);

    rg_REAL computeXAreaBetweenSidechainsNBackbone(Residue** residues, 
                                                   const rg_INT& numResidues, 
                                                   Atom** backboneAtoms, 
                                                   const rg_INT& numBackboneAtoms,
                                                   rg_REAL*& xAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                   rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                   rg_INT& numXBetweenSidechainsNBackbone);

    rg_REAL computeXAreaBetweenSidechainPairs(Residue** residues, 
                                              const rg_INT& numResidues,
                                              rg_REAL*& xAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                              rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                              rg_INT& numXBetweenSidechainPairs);


rg_INT computeXCircleAreaOfResidueWithBackboneNOtherSidechains(Molecule& protein,
                                                                rg_REAL*& xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL*& xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL*& xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL*& xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                                rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                                rg_REAL& sumXCircleAreaOfProtein,
                                                                rg_REAL& sumXCircleRadiusOfProtein,
                                                                rg_INT&  numXOfProtein);

    rg_REAL computeXCircleAreaBetweenSidechainsNBackbone(Residue** residues, 
                                                         const rg_INT& numResidues, 
                                                         Atom** backboneAtoms, 
                                                         const rg_INT& numBackboneAtoms,
                                                         rg_REAL*& xCircleAreaWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_REAL*& xCircleRadiusWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_INT*&  numXWithBackboneOfResidue_orderedByChainIDSeqNum,
                                                         rg_REAL& sumXCircleRadiusBetweenSidechainsNBackbone,
                                                         rg_INT& numXBetweenSidechainsNBackbone);

    rg_REAL computeXCircleAreaBetweenSidechainPairs(Residue** residues, 
                                                    const rg_INT& numResidues,
                                                    rg_REAL*& xCircleAreaWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                    rg_REAL*& xCircleRadiusWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                    rg_INT*&  numXWithOtherSidechainsOfResidue_orderedByChainIDSeqNum,
                                                    rg_REAL& sumXCircleRadiusBetweenSidechainPairs,
                                                    rg_INT& numXBetweenSidechainPairs);

rg_REAL computeVanDerWaalsEnergyOfProtein(Molecule& protein);

rg_REAL computeVanDerWaalsEnergyBetweenSidechainsNBackbone(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeVanDerWaalsEnergyBetweenSidechainPairs(Residue** residues, const rg_INT& numResidues);

rg_REAL computeElectrostaticEnergyOfProtein(Molecule& protein);

rg_REAL computeElectrostaticEnergyBetweenSidechainsNBackbone(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeElectrostaticEnergyBetweenSidechainPairs(Residue** residues, const rg_INT& numResidues);

rg_REAL computeHydrogenBondingEnergyOfProtein(Molecule& protein);

rg_REAL computeHydrogenBondingEnergyBetweenSidechainsNBackbone(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeHydrogenBondingEnergyBetweenSidechainPairs(Residue** residues, const rg_INT& numResidues);

rg_REAL computeVanDerWaalsEnergyOfProteinBySplineFitting(Molecule& protein);

rg_REAL computeVanDerWaalsEnergyBetweenSidechainsNBackboneBySplineFitting(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeVanDerWaalsEnergyBetweenSidechainNBackbone(Residue* residue, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeVanDerWaalsEnergyBetweenSidechainPairsBySplineFitting(Residue** residues, const rg_INT& numResidues);

rg_REAL computeVanDerWaalsEnergyBetweenTwoSidechains(Residue* residue1, Residue* residue2);

rg_REAL computeVanDerWaalsEnergyBetweenTwoResidues(Residue* residue1, Residue* residue2);

rg_REAL computeVanDerWaalsEnergyBetweenTwoAtomSets( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2 );

rg_REAL computeVanDerWaalsEnergyBetweenTwoAtomSetsBySplineFitting( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2 );

rg_REAL computeVanDerWaalsEnergyBetweenTwoAtoms( Atom* atom1, Atom* atom2 );

Vector3D evaluateGradientVectorForVanDerWaalsEnergyBetweenTwoAtomsWRTCoordinate(Atom* atom1, Atom* atom2);

Vector3D evaluateGradientVectorForNewScoreFunctionBetweenTwoAtomsWRTCoordinate(Atom* atom1, Atom* atom2);

Vector3D evaluateGradientVectorForPositionBetweenTwoAtomsBySplineFitting(Atom* atom1, Atom* atom2);

rg_REAL computeDerivativeOfVanDerWaalsEnergyBetweenTwoAtoms(Atom* atom1, Atom* atom2, const rg_REAL& distCenters);

// SCWRL3
// The following functions for computing energy are based on the article
// Yang Cao and Lin Song and Zhichao Miao and Yun Hu and Liqing Tian and and Taijiao Jiang: 
// Improved side-chain modeling by coupling clash-detection guided iterative search with rotamer relaxation,
// Bioinformatics, 27 (6), 785-790, 2011

rg_REAL  computeVanDerWaalsEnergyBySCWRL3(Molecule& protein);

rg_REAL  computeVanDerWaalsEnergyBetweenSidechainNBackboneBySCWRL3(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL  computeVanDerWaalsEnergyBetweenSidechainsBySCWRL3(Residue** residues, const rg_INT& numResidues);

rg_REAL  computeVanDerWaalsEnergyBetweenTwoAtomSetsBySCWRL3( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_INT& numXAtoms );

Vector3D evaluateGradientVectorForPositionBetweenTwoAtomsBySCWRL3(Atom* atom1, Atom* atom2);

rg_BOOL  detectClashBetweenRotamersByEvaluatingVanDerWaalsEnergy( Rotamer*  targetRotamer, 
	                                                              const rg_INT& targetResidueIndex, 
																  ManagerOfRotamers_Assigned_At_Residues* managerOfRotamers_Assigned_At_Residues, 
																  Backbone* backbone,
	                                                              rg_dList<rg_INDEX>& clashedResidueIndeices);

void    getRotamerSetsForResiduesSortedByChainIDSeqNumber(Rotamer* sortedAssignedRotamersOfProtein, const rg_INT& numRotamers, RotamerSetOfResidue*& sortedRotamerSetsOfResidues);

rg_REAL computeXVolumeOfProtein(Molecule& protein);

rg_REAL computeXVolumeBetweenSidechainsNBackbone(Residue** residues, const rg_INT& numResidues, Atom** backboneAtoms, const rg_INT& numBackboneAtoms);

rg_REAL computeXVolumeBetweenSidechainPairs(Residue** residues, const rg_INT& numResidues);

rg_BOOL haveIntersectionBetweenTwoResidues(Residue* residue1, Residue* residue2);

//void computeSumPairwiseXVolOfBallWithBallSet(const Sphere& ball1, rg_dList<Sphere>& ballSet2, rg_REAL& xVolume);

rg_REAL computeSumPairwiseXVolOfAtomWithAtomSet(Atom* atom1, rg_dList<Atom*>& atomSet2);

rg_REAL computeSumPairwiseXVolOfAtomWithAtomSet(Atom* atom1, rg_dList<Atom*>& atomSet2, rg_INT& numXAtoms);

rg_REAL computeSumPairwiseXVolBetweenTwoAtomSets( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_INT& numXAtoms );
rg_REAL computeSumLJPotentialEnergyBetweenTwoAtomSetsByMappingXVol2LJEnergy( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_INT& numXAtoms );

rg_REAL computeNewScoreFunction_old( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_INT& numXAtoms );
rg_REAL computeNewScoreFunction( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_INT& numXAtoms );

rg_REAL scale(const rg_REAL& valueBeforeScaling);

rg_REAL computeSumPairwiseXAreaBetweenTwoAtomSets( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_INT& numXAtoms );

rg_REAL computeSumPairwiseXCircleAreaBetweenTwoAtomSets( Atom** atomSet1, const rg_INT& numOfAtomSet1, Atom** atomSet2, const rg_INT& numOfAtomSet2, rg_REAL& sumXCircleRadii, rg_INT& numXAtoms );

rg_BOOL haveIntersectionBetweenAtomWithAtomSet(Atom* atom1, rg_dList<Atom*>& atomSet2);

void    computeSumPairwiseXVolOfAtomWithAtomSetWithoutChemBond(Atom* atom1, rg_dList<Atom*>& atomSet2, rg_REAL& xVolume);

rg_REAL computeRMSDBetweenTwoProteinsSharingCommonBackboneNHavingDifferentSidechains(Molecule& protein1, Molecule& protein2);

rg_BOOL computeRMSDsBetweenAminoResiduePairsOfTwoProteinsSharingFixedCommonBackbone(Molecule& protein1, 
                                                                                    Molecule& protein2, 
                                                                                    rg_REAL& rmsdOfProtein,
                                                                                    vector<rg_REAL>& RMSDsOrderedByChainIDSeqNumbers);

void writeDihedralAnglesOfAminoResidueSidechains(Molecule& molecule, 
                                                 const string& sidechainDihedralAngleFIle, 
                                                 const rg_BOOL& bAppend = rg_FALSE);
    void writeHeaderForDihedralAnglesOfAminoResidueSidechains(ofstream& fout);
    void writeBodyForDihedralAnglesOfAminoResidueSidechains(Molecule& molecule, ofstream& fout);

void writeDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(Molecule& referenceStructure, 
                                                                                          Molecule& targetStructure,
                                                                                          const string& sidechainDihedralAngleFIle,
                                                                                          const rg_BOOL& bAppend = rg_FALSE);
    void writeHeaderForDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(ofstream& fout);
    void writeBodyForDihedralAnglesOfAminoResidueSidechains_N_DihedralAngleDiffWithReferenceMolecule(Molecule& referenceStructure,
                                                                                                     Molecule& targetStructure,
                                                                                                     ofstream& fout);

void    analyzeSidechainDihedralAngles(const string& sidechainDihedralAngleFIle);

// CMKIM added on April 7th, 2011 for JHRYU
void    getChemicalBondsInResidueStartFromBetaCarbon( Residue* aResidue, rg_dList<ChemicalBond*>& targetListOfChemicalBonds );


// JKKIM added on 2011-07-12
void    getRotationalBonds(Molecule& aMolecule, rg_dList<ChemicalBond*>& rotationalBonds);
    void    removeBondOfHydrogenAtom(rg_dList<ChemicalBond*>& rotationalBonds);
    void    removeDoubleTripleAromaticAndAmideBond(rg_dList<ChemicalBond*>& rotationalBonds);
    void    removeThreeFoldSymmetryBond(rg_dList<ChemicalBond*>& rotationalBonds);
    void    removeBondOfAminoCarboxylAndHydroxylGroup(rg_dList<ChemicalBond*>& rotationalBonds);
    void    removeBondOfRing(rg_dList<ChemicalBond*>& rotationalBonds, rg_dList<ChemicalBond>* allChemicalBonds);
        rg_BOOL isChemicalBondOnRing(ChemicalBond* aChemicalBond, rg_dList<ChemicalBond>* allChemicalBonds);
            void getIncidentChemicalBondsAndKillIncidentBonds(Atom* atom, rg_dList<ChemicalBond*>& allChemicalBonds, rg_dList<ChemicalBond*>& incidentChemicalBonds );

void    getAtomSetConnectedWithFirstAtomOnGivenBond( Molecule* aMolecule, ChemicalBond* givenBond, rg_dList<Atom*>& atomSetConnectedWithFirstAtom );

void    rotateAtomsByGivenBond( Molecule* aMolecule, ChemicalBond* givenBond, const double& angleByRadian);



} // namespace GeometryTier
} // namespace V


#endif

