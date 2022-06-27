#pragma warning(disable:4390)

#include "rg_RelativeOp.h"
#include "rg_Molecule.h"
#include "ForceField.h"
#include "ConstForPharmacophore.h"
#include "MoleculeIOFunctions.h"
#include "FunctionsForMolecule.h"
#include "rg_TMatrix3D.h"
#include "FunctionsForRotamerLibrary.h"
using namespace V::GeometryTier;

#include <algorithm>
#include <queue>
using namespace std;



Molecule::Molecule()
{
    m_moleculeFileName = "";
    m_atomRadiusType   = VDW_BONDI_ATOM_RADIUS;

    m_modelSerialFromInputFile = 1;

    // JKKIM ADDED //
    m_moleculeName = "";

    // by Y.Cho at 2011-10-12
    m_timeStamp = "";
    m_fileSize  = 0;

	m_method = "NONE";
	m_resolution = -1.0;
}



Molecule::Molecule( const Molecule& aMolecule )
{
    duplicate(aMolecule);
 //   m_moleculeFileName         = aMolecule.m_moleculeFileName;
 //   m_atomRadiusType           = aMolecule.m_atomRadiusType;
 //   m_headerRecords            = aMolecule.m_headerRecords;
 //   m_atoms                    = aMolecule.m_atoms;
 //   m_chemicalBonds            = aMolecule.m_chemicalBonds;
 //   m_residues                 = aMolecule.m_residues;
 //   m_chains                   = aMolecule.m_chains;
 //   m_centerOfMass             = aMolecule.m_centerOfMass;
 //   m_pharmaFeatures           = aMolecule.m_pharmaFeatures;
 //   m_modelSerialFromInputFile = aMolecule.m_modelSerialFromInputFile;
 //   // JKKIM ADDED //
 //   m_moleculeName             = aMolecule.m_moleculeName;

 //   // by Y.Cho at 2011-10-12
 //   m_timeStamp                = aMolecule.m_timeStamp;
 //   m_fileSize                 = aMolecule.m_fileSize;

	//m_isCryst = false;

	//for(int i = 0; i < 6; i++) {
	//	m_cryst[i] = 0.0f;
	//}
}



Molecule::~Molecule()
{
}



void Molecule::duplicate(const Molecule& aMolecule)
{
    //  duplicate structure
    map<Atom*, Atom*>                   mapAtom2Atom;
    map<Residue*, Residue*>             mapResidue2Residue;
    map<Chain*, Chain*>                 mapChain2Chain;
    map<ChemicalBond*, ChemicalBond*>   mapBond2Bond;

    makeMapsOfEntitiesInMolecule(
        aMolecule, mapAtom2Atom, mapResidue2Residue, mapChain2Chain, mapBond2Bond);

    setEntitiesAndRelationshipAmongEntities(
        mapAtom2Atom, mapResidue2Residue, mapChain2Chain, mapBond2Bond);


    //  duplicate properties
    m_moleculeFileName = aMolecule.m_moleculeFileName;
    m_atomRadiusType   = aMolecule.m_atomRadiusType;
    m_headerRecords    = aMolecule.m_headerRecords;

    m_centerOfMass              = aMolecule.m_centerOfMass;
    m_minEnclosingSphere        = aMolecule.m_minEnclosingSphere;
    m_modelSerialFromInputFile  = aMolecule.m_modelSerialFromInputFile;

    // JKKIM ADDED //
    m_moleculeName = aMolecule.m_moleculeName;
    m_isCryst = aMolecule.m_isCryst;
    for (int i = 0; i < 6; i++) {
        m_cryst[i] = aMolecule.m_cryst[i];
    }

    // by Y.Cho at 2011-10-12
    m_timeStamp     = aMolecule.m_timeStamp;
    m_fileSize      = aMolecule.m_fileSize;
    m_method        = aMolecule.m_method;
    m_resolution    = aMolecule.m_resolution;


    for (list<PharmaFeature>::iterator i_feature = m_pharmaFeatures.begin(); i_feature != m_pharmaFeatures.end(); ++i_feature) {
        PharmaFeature* oldFeature = &(*i_feature);
        this->m_pharmaFeatures.push_back(PharmaFeature(oldFeature->getID()));

        list<Atom*> atomsInNewFeature;
        rg_dList<Atom*>* atomsInOldFeature = oldFeature->getAtoms();
        atomsInOldFeature->reset4Loop();
        while (atomsInOldFeature->setNext4Loop()) {
            Atom* oldAtom = atomsInOldFeature->getEntity();

            map<Atom*, Atom*>::iterator i_atom = mapAtom2Atom.find(oldAtom);
            if (i_atom != mapAtom2Atom.end()) {
                atomsInNewFeature.push_back(i_atom->second);
            }
        }

        this->m_pharmaFeatures.back().setTypesAndAtoms(oldFeature->getPharmaFeatureType(), oldFeature->getChemFuncGroupType(), &atomsInNewFeature);
        this->m_pharmaFeatures.back().setCentersOfPharmaFeature();
    }
}



void    Molecule::makeMapsOfEntitiesInMolecule(
            const Molecule& aMolecule,
            map<Atom*, Atom*>&           mapAtom2Atom,
            map<Residue*, Residue*>&     mapResidue2Residue,
            map<Chain*, Chain*>&         mapChain2Chain,
            map<ChemicalBond*, ChemicalBond*>&   mapBond2Bond)
{
    aMolecule.m_atoms.reset4Loop();
    while (aMolecule.m_atoms.setNext4Loop()) {
        Atom*   oldAtom = aMolecule.m_atoms.getpEntity();
        Atom*   newAtom = this->m_atoms.add( Atom( oldAtom->getID() ) );

        mapAtom2Atom.insert(make_pair(oldAtom, newAtom));
    }


    aMolecule.m_residues.reset4Loop();
    while (aMolecule.m_residues.setNext4Loop()) {
        Residue* oldResidue = aMolecule.m_residues.getpEntity();
        Residue* newResidue = this->m_residues.add(Residue(oldResidue->getID()));

        mapResidue2Residue.insert(make_pair(oldResidue, newResidue));
    }


    aMolecule.m_chains.reset4Loop();
    while (aMolecule.m_chains.setNext4Loop()) {
        Chain*  oldChain = aMolecule.m_chains.getpEntity();
        Chain*  newChain = this->m_chains.add(Chain(oldChain->getID()));

        mapChain2Chain.insert(make_pair(oldChain, newChain));
    }


    aMolecule.m_chemicalBonds.reset4Loop();
    while (aMolecule.m_chemicalBonds.setNext4Loop()) {
        ChemicalBond* oldBond = aMolecule.m_chemicalBonds.getpEntity();
        ChemicalBond* newBond = this->m_chemicalBonds.add(ChemicalBond(oldBond->getID()));

        mapBond2Bond.insert(make_pair(oldBond, newBond));
    }
}


void    Molecule::setEntitiesAndRelationshipAmongEntities(
            const map<Atom*, Atom*>&           mapAtom2Atom,
            const map<Residue*, Residue*>&     mapResidue2Residue,
            const map<Chain*, Chain*>&         mapChain2Chain,
            const map<ChemicalBond*, ChemicalBond*>&   mapBond2Bond)
{
    map<Atom*, Atom*>::const_iterator                   i_atom;
    map<Residue*, Residue*>::const_iterator             i_res;
    map<Chain*, Chain*>::const_iterator                 i_chain;
    map<ChemicalBond*, ChemicalBond*>::const_iterator   i_bond;


    for (i_atom = mapAtom2Atom.begin(); i_atom != mapAtom2Atom.end(); ++i_atom) {
        Atom* oldAtom = i_atom->first;
        Atom* newAtom = i_atom->second;

        newAtom->setID(                     oldAtom->getID() );
        newAtom->setAtomCode(               oldAtom->getAtomCode());
        newAtom->setAtomBall(               oldAtom->getAtomBall());
        newAtom->setChemicalProperties(     oldAtom->getChemicalProperties());
        newAtom->setSerialFromInputFile(    oldAtom->getSerialFromInputFile());
        newAtom->setAtomNameFromInputFile(  oldAtom->getAtomNameFromInputFile());

        i_res = mapResidue2Residue.find(oldAtom->getResidue());
        if (i_res != mapResidue2Residue.end()) {
            newAtom->setResidue(i_res->second);
        }

        rg_dList<ChemicalBond*>* bondsInOldAtom = oldAtom->getListChemicalBond();
        bondsInOldAtom->reset4Loop();
        while (bondsInOldAtom->setNext4Loop()) {
            ChemicalBond* oldBond = bondsInOldAtom->getEntity();
            i_bond = mapBond2Bond.find(oldBond);
            if (i_bond != mapBond2Bond.end()) {
                newAtom->addChemicalBond(i_bond->second);
            }
        }
    }


    for (i_res = mapResidue2Residue.begin(); i_res != mapResidue2Residue.end(); ++i_res) {
        Residue* oldResidue = i_res->first;
        Residue* newResidue = i_res->second;

        newResidue->setID(              oldResidue->getID());
        newResidue->setResidueCode(     oldResidue->getResidueCode());
        newResidue->setSequenceNumber(  oldResidue->getSequenceNumber());
        newResidue->setResidueName(     oldResidue->getResidueName());

        i_chain = mapChain2Chain.find(oldResidue->getChain());
        if (i_chain != mapChain2Chain.end()) {
            newResidue->setChain( i_chain->second );
        }

        rg_dList<Atom*>* atomsInOldResidue = oldResidue->getAtoms();
        atomsInOldResidue->reset4Loop();
        while (atomsInOldResidue->setNext4Loop()) {
            Atom* oldAtom = atomsInOldResidue->getEntity();

            i_atom = mapAtom2Atom.find(oldAtom);
            if (i_atom != mapAtom2Atom.end()) {
                newResidue->addAtom(i_atom->second);
            }
        }
    }


    for (i_chain = mapChain2Chain.begin(); i_chain != mapChain2Chain.end(); ++i_chain) {
        Chain* oldChain = i_chain->first;
        Chain* newChain = i_chain->second;

        newChain->setID(                            oldChain->getID());
        newChain->setChainIDFromInputFileInDecimal( oldChain->getChainIDFromInputFileInDecimal());
        newChain->setMolecule(                      this);
        newChain->setChainCode(                     oldChain->getChainCode());

        rg_dList<Residue*>* residuesInOldChain = oldChain->getResidues();
        residuesInOldChain->reset4Loop();
        while (residuesInOldChain->setNext4Loop()) {
            Residue* oldResidue = residuesInOldChain->getEntity();
            i_res = mapResidue2Residue.find(oldResidue);
            if (i_res != mapResidue2Residue.end()) {
                newChain->addResidue(i_res->second);
            }
        }
    }


    for (i_bond = mapBond2Bond.begin(); i_bond != mapBond2Bond.end(); ++i_bond) {
        ChemicalBond* oldBond = i_bond->first;
        ChemicalBond* newBond = i_bond->second;

        newBond->setID(                  oldBond->getID());
        newBond->setTypeOfBond(          oldBond->getTypeOfBond());
        newBond->setSerialFromInputFile( oldBond->getSerialFromInputFile());

        i_atom = mapAtom2Atom.find(oldBond->getFirstAtom());
        if (i_atom != mapAtom2Atom.end()) {
            newBond->setFirstAtom(i_atom->second);
        }

        i_atom = mapAtom2Atom.find(oldBond->getSecondAtom());
        if (i_atom != mapAtom2Atom.end()) {
            newBond->setSecondAtom(i_atom->second);
        }
    }
}



string Molecule::getMoleculeFileName() const
{
    return m_moleculeFileName;
}



AtomRadiusType Molecule::getAtomRadiusType() const
{
    return m_atomRadiusType;
}



string Molecule::getDescriptionOfAtomRadiusType() const
{
    return ATOM_RADIUS_DESCRIPTION[ (int)m_atomRadiusType ];
}



rg_dList<Atom>* Molecule::getAtoms()
{
    return &m_atoms;
}



void Molecule::getPtrAtoms(rg_dList<Atom*>& ptrAtoms) const
{    
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        ptrAtoms.add( m_atoms.getpEntity() );
    }
}



rg_dList<ChemicalBond>* Molecule::getChemicalBonds()
{
    return &m_chemicalBonds;
}



rg_dList<Residue>* Molecule::getResidues()
{
    return &m_residues;
}

rg_INT Molecule::getNonAminoResidues(Residue**& residueArr)
{
	rg_INT numNonAminoResidues = getNumNonAminoResidues();

	residueArr = new Residue*[numNonAminoResidues];

	rg_INDEX residueIndex = 0;
	m_residues.reset4Loop();
	while (m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		if (!currResidue->isAminoResidue())
			residueArr[residueIndex++] = currResidue;
	}

	return numNonAminoResidues;
}

rg_INT Molecule::getAminoResidues(Residue**& residueArr)
{
	rg_INT numResidues = getNumAminoResidues();

	residueArr = new Residue*[numResidues];

	rg_INDEX residueIndex = 0;
	m_residues.reset4Loop();
	while ( m_residues.setNext4Loop() ) 
	{
		Residue* currResidue = m_residues.getpEntity();
		if(currResidue->isAminoResidue())
			residueArr[residueIndex++] = currResidue;
	}

	return numResidues;
}

rg_INT Molecule::getNumNonAminoResidues()
{
	rg_INT numNonAminoResidues = 0;
	m_residues.reset4Loop();
	while (m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		if (!currResidue->isAminoResidue())
			numNonAminoResidues++;
	}
	return numNonAminoResidues;
}

rg_INT Molecule::getNumAminoResidues()
{
	rg_INT numAminoResidues = 0;
	m_residues.reset4Loop();
	while(m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		if(currResidue->isAminoResidue())
			numAminoResidues++;
	}
	return numAminoResidues;
}


rg_dList<Chain>* Molecule::getChains()
{
    return &m_chains;
}



Chain* Molecule::getChain( const rg_INT& chainIDInDecimal )
{
    Chain* targetChain = rg_NULL;

    m_chains.reset4Loop();
    while ( m_chains.setNext4Loop() ) {
        Chain* currChain = m_chains.getpEntity();
        if ( chainIDInDecimal == currChain->getChainIDFromInputFileInDecimal() ) {
            targetChain = currChain;
            break;
        }
    }
    return targetChain;
}



void Molecule::getResiduesSortedBySequenceNumber( rg_dList<Residue*>& targetListOfResidues )
{
    rg_INT numOfResidues = m_residues.getSize();

    Residue** sortedResidues = new Residue*[numOfResidues];

    rg_INT i_residue = 0;
    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        sortedResidues[i_residue] = m_residues.getpEntity();
        i_residue++;
    }

    qsort( (void *)sortedResidues, numOfResidues, sizeof(Residue*), compareResidueBySequenceNumber );

    targetListOfResidues.removeAll();
    
    for ( i_residue=0; i_residue<numOfResidues; i_residue++ ) {
        targetListOfResidues.addTail( sortedResidues[i_residue] );
    }
    
    delete [] sortedResidues;
}

void Molecule::getResiduesSortedBySequenceNumberForEachChain(rg_dList<Residue*>& targetListOfResidues)
{
	m_chains.reset4Loop();
	while (m_chains.setNext4Loop())
	{
		Chain* currChain = m_chains.getpEntity();
		rg_dList<Residue*> residues;
		currChain->getResiduesSortedBySequenceNumber(residues);
		targetListOfResidues.append(residues);				
	}
}

rg_INT Molecule::getResiduesSortedBySequenceNumberForEachChain(Residue**& sortedResidues)
{
	rg_dList<Residue*> sortedResidueList;
	getResiduesSortedBySequenceNumberForEachChain( sortedResidueList );

	rg_INT numResidues = sortedResidueList.getSize();
	sortedResidues = new Residue* [numResidues];

	rg_INDEX i = 0;

	sortedResidueList.reset4Loop();
	while (sortedResidueList.setNext4Loop())
	{
		sortedResidues[ i++ ] = sortedResidueList.getEntity();
	}
	return numResidues;
}

rg_INT Molecule::getAllResiduesOrderedByChainIDSequenceNumber(Residue **& sortedResidues)
{
	int numResidues = m_residues.getSize();
	Residue** residues_ordered = new Residue*[numResidues];

	int i = 0;
	m_residues.reset4Loop();
	while(m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		residues_ordered[i] = currResidue;
		i++;
	}

	if (numResidues > 0)
		qsort((void *)sortedResidues, numResidues, sizeof(Residue*), compareResidueByChainIDSequenceNumber);

	return numResidues;
}

rg_INT Molecule::getAminoResiduesOrderedByChainIDSequenceNumber( Residue**& sortedResidues )
{
	rg_INT numResidues = getAminoResidues( sortedResidues );	

    if(numResidues > 0)
	    qsort( (void *)sortedResidues, numResidues, sizeof(Residue*), compareResidueByChainIDSequenceNumber );

	return numResidues;
}

rg_INT Molecule::getNonAminoResiduesOrderedByChainIDSequenceNumber(Residue **& sortedResidues)
{
	rg_INT numResidues = getNonAminoResidues(sortedResidues);

	if (numResidues > 0)
		qsort((void *)sortedResidues, numResidues, sizeof(Residue*), compareResidueByChainIDSequenceNumber);

	return numResidues;
}

rg_INT Molecule::getAminoResiduesOrderedByChainIDSequenceNumber( rg_dList<Residue*>& residuesOrderedByChainIDSeqNum )
{
	rg_INT numResidues = 0;
	Residue** orderedResidues = rg_NULL;
	numResidues = getAminoResiduesOrderedByChainIDSequenceNumber(orderedResidues);

	rg_INDEX i;
	for (i = 0;i < numResidues;i++)
	{
		residuesOrderedByChainIDSeqNum.add( orderedResidues[ i ] );
	}

	if(orderedResidues != rg_NULL)
		delete [] orderedResidues;

	return numResidues;
}

void Molecule::getAtomsOnBackboneSortedBySequenceNumber( rg_dList<Atom*>& listOfAtomsOnBackbone )
{
    listOfAtomsOnBackbone.removeAll();

    rg_dList<Residue*> targetListOfResidues;
    getResiduesSortedBySequenceNumber( targetListOfResidues );

    targetListOfResidues.reset4Loop();
    while ( targetListOfResidues.setNext4Loop() ) {      

        rg_dList<Atom*> atomsOnBackboneInResidue;
        targetListOfResidues.getEntity()->getAtomsOnBackbone( &atomsOnBackboneInResidue );

        atomsOnBackboneInResidue.reset4Loop();
        while ( atomsOnBackboneInResidue.setNext4Loop() ) {
            listOfAtomsOnBackbone.addTail( atomsOnBackboneInResidue.getEntity() );
        }
    }
}

void Molecule::getAtomsOnBackboneSortedByChainIDSequenceNumber( rg_dList<Atom*>& listOfAtomsOnBackbone )
{
	Residue** sortedResidues = rg_NULL;
	rg_INT numResidues = getAminoResiduesOrderedByChainIDSequenceNumber(sortedResidues);

	listOfAtomsOnBackbone.removeAll();

	rg_INDEX i;
	for (i = 0;i < numResidues;i++)
	{
		rg_dList<Atom*> atomsOnBackboneInResidue;
		sortedResidues[ i ]->getAtomsOnBackbone(&atomsOnBackboneInResidue);

#ifdef _DEBUG
        if(atomsOnBackboneInResidue.getSize() != 4)
            int here = 4;
#endif
		listOfAtomsOnBackbone.append( atomsOnBackboneInResidue );
	}

	if(sortedResidues != rg_NULL)
		delete [] sortedResidues;
}


rg_INT Molecule::getAtomsOnBackboneSortedByChainIDSequenceNumber( Atom**& backboneAtoms )
{
	rg_dList<Atom*> atomsOnBackbone;
	getAtomsOnBackboneSortedByChainIDSequenceNumber(atomsOnBackbone);
	rg_INT numAtoms = atomsOnBackbone.getSize();

    if(numAtoms > 0)
	    backboneAtoms = new Atom*[ numAtoms ];
	rg_INDEX i = 0;

	atomsOnBackbone.reset4Loop();
	while (atomsOnBackbone.setNext4Loop())
	{
		Atom* currAtom = atomsOnBackbone.getEntity();
		backboneAtoms[ i++ ] = currAtom;
	}

	return numAtoms;
}


rg_INT Molecule::getSumAtomNumsOfAminoResidues() const
{
	rg_INT sumAtomNums = 0;
	m_residues.reset4Loop();
	while(m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		if(currResidue->isAminoResidue())
		{			
			rg_dList<Atom*> atomsOfResidue;
			FunctionsForRotamerLibrary::getAtomsOfResidue(currResidue, atomsOfResidue);
			rg_INT atomNums = 0;
			atomNums = atomsOfResidue.getSize();
			sumAtomNums += atomNums;
		}
	}
	return sumAtomNums;
}


rg_BOOL Molecule::haveMissingResidues() const
{
	rg_BOOL bHaveMissingResidues = rg_FALSE;

	m_chains.reset4Loop();
	while (m_chains.setNext4Loop())
	{
		Chain* currChain = rg_NULL;
		currChain = m_chains.getpEntity();
		if(currChain->haveMissingResidues())
		{
			bHaveMissingResidues = rg_TRUE;
			break;
		}
	}

	return bHaveMissingResidues;
}


rg_INT Molecule::findMissingResidues(rg_dList<ChainIDWithSequenceNumbersOfMissingResidues>& chainIDsWithSeqNumbersOfMissingResidues) const
{
	rg_BOOL bHaveMissingResidues = rg_FALSE;
	bHaveMissingResidues = haveMissingResidues();
	if(bHaveMissingResidues)
	{
		m_chains.reset4Loop();
		while (m_chains.setNext4Loop())
		{
			Chain* currChain = rg_NULL;
			currChain = m_chains.getpEntity();
			if(currChain->haveMissingResidues())
			{
				string chainIDFromInputPDBFile = currChain->getChainIDFromInputFileInString();
				rg_dList<rg_INT> seqNumbersOfMissingResidues;
				rg_INT numMissingResidues = 0;
				numMissingResidues = currChain->findMissingResidues(seqNumbersOfMissingResidues);
				if(numMissingResidues > 0)
				{
					ChainIDWithSequenceNumbersOfMissingResidues currChainIDWithSeqNumbersOfMissingResidues(chainIDFromInputPDBFile, seqNumbersOfMissingResidues);
					chainIDsWithSeqNumbersOfMissingResidues.add( currChainIDWithSeqNumbersOfMissingResidues );
				}
			}
		}

		return chainIDsWithSeqNumbersOfMissingResidues.getSize();
	}
	else
	{
		return 0;
	}
}


void Molecule::getListOfHAtomsFromHDonor( rg_dList<Atom*>* targetListOfHydrogenAtoms )  // TOBE CHECKED
{
    m_atoms.reset4Loop();
    
    while( m_atoms.setNext4Loop() ) {
        Atom* currAtom = m_atoms.getpEntity();
        if( currAtom->getAtomCode() == H_ATOM && currAtom->getListChemicalBond()->getSize() != 0 ) {
            Atom* firstBondedAtom = currAtom->getListChemicalBond()->getFirstEntity()->getFirstAtom();
            Atom* secondBondedAtom = currAtom->getListChemicalBond()->getFirstEntity()->getSecondAtom();


            if( ( firstBondedAtom->getAtomCode()  == O_ATOM || firstBondedAtom->getAtomCode()  == N_ATOM ) ||
                ( secondBondedAtom->getAtomCode() == O_ATOM || secondBondedAtom->getAtomCode() == N_ATOM )  )
                targetListOfHydrogenAtoms->addTail( currAtom );
        }
    }
}



void Molecule::getListOfHAcceptorAtoms( rg_dList<Atom*>* targetListOfHAcceptorAtoms )
{
    m_atoms.reset4Loop();
    
    while( m_atoms.setNext4Loop() ) {
        Atom* currAtom = m_atoms.getpEntity();
        if( currAtom->getChemicalProperties().getNumOfAcceptableHydrogen() >= 1 )
            targetListOfHAcceptorAtoms->addTail( currAtom );
    }
}



rg_Point3D Molecule::getCenterOfMass()
{
    return m_centerOfMass;
}



Sphere Molecule::getMinEnclosingSphere()
{
    return m_minEnclosingSphere;
}



rg_REAL Molecule::getMeanRadiusOfAtoms()
{
    rg_REAL sumOfRadius = 0.0;
    rg_INT  numOfAtoms  = m_atoms.getSize();

    m_atoms.reset4Loop();
    while( m_atoms.setNext4Loop() ) {
        Atom* currAtom = m_atoms.getpEntity();
        sumOfRadius += currAtom->getAtomBall().getRadius();
    }

    rg_REAL meanRadiusOfAtoms = sumOfRadius/numOfAtoms;

    return meanRadiusOfAtoms;
}



rg_REAL Molecule::getMaxRadiusOfAtoms()
{
    rg_REAL maxRadiusOfAtoms = rg_MAX_REAL*(-1);
    m_atoms.reset4Loop();
    while( m_atoms.setNext4Loop() ) {
        rg_REAL radiusOfCurrAtom = m_atoms.getpEntity()->getAtomBall().getRadius();
        if( rg_GT( radiusOfCurrAtom, maxRadiusOfAtoms ) )
            maxRadiusOfAtoms = radiusOfCurrAtom;
    }

    return maxRadiusOfAtoms;
}



rg_REAL Molecule::getMolecularWeight()
{
    rg_REAL molecularWeight = 0.0;

    m_atoms.reset4Loop();

    while( m_atoms.setNext4Loop() ) {
        molecularWeight +=  ATOM_FEATURES[ m_atoms.getpEntity()->getAtomCode() ].weight;
    }
    
    return molecularWeight;
}



//  by Youngsong Cho, 2012.03.14.
rg_INT Molecule::getNumStandardResidues() const
{
    rg_INT numStandardResidues = 0;

    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        Residue* currRes = m_residues.getpEntity();

        if ( currRes->isStandardResidue() ) {
            numStandardResidues++;
        }
    }

    return numStandardResidues;
}



rg_INT Molecule::getNumProteinChains() const
{
	rg_INT numProteinChains = 0;

	m_chains.reset4Loop();
	while (m_chains.setNext4Loop())
	{
		Chain* currChain = m_chains.getpEntity();
		if(currChain->getChainCode() == PROTEIN_CHAIN)
		{
			numProteinChains++;
		}
	}

	return numProteinChains;
}



rg_INT Molecule::getNumberOfChainsWithNonStandardResidues()
{
    rg_INT numberOfChainsWithNonStandardResidues = 0;

    m_chains.reset4Loop();
    while ( m_chains.setNext4Loop() ) {
        Chain* currChain = m_chains.getpEntity();
        if( currChain->getChainCode() == UNK_CHAIN ) {
            numberOfChainsWithNonStandardResidues++;
        }
    }

    return numberOfChainsWithNonStandardResidues;
}



rg_INT Molecule::getNumberOfChainsWithNonStandardResiduesAndHOH()
{
    rg_INT numberOfChainsWithNonStandardResiduesAndHOH = 0;
    
    m_chains.reset4Loop();
    while ( m_chains.setNext4Loop() ) {
        Chain* currChain = m_chains.getpEntity();
        if( currChain->getChainCode() == UNK_CHAIN ||
            currChain->getChainCode() == HOH_CHAIN          ) {
            numberOfChainsWithNonStandardResiduesAndHOH++;
        }
    }
    
    return numberOfChainsWithNonStandardResiduesAndHOH;
}



rg_INT Molecule::getModelSerialFromInputFile() const
{
    return m_modelSerialFromInputFile;
}


rg_INT  Molecule::getNumberOfAtomsExcludingHOH() const
{
    rg_INT numAtomsExcludingHOH = 0;

    m_atoms.reset4Loop();
    while (m_atoms.setNext4Loop()) {
        Atom* currAtom = m_atoms.getpEntity();
        Residue* currRes = currAtom->getResidue();

        if (currRes == rg_NULL) {
            continue;
        }
        if (currRes->getResidueCode() == HOH_RESIDUE) {
            continue;
        }

        ++numAtomsExcludingHOH;
    }

    return numAtomsExcludingHOH;
}



rg_INT  Molecule::getNumberOfAtomsInStandardResidues() const
{
    rg_INT numAtomsInStandardResidues = 0;

    m_atoms.reset4Loop();
    while (m_atoms.setNext4Loop()) {
        Atom* currAtom = m_atoms.getpEntity();
        Residue* currRes = currAtom->getResidue();
        if (currRes == rg_NULL) continue;

        if (currRes->isStandardResidue()) {
            ++numAtomsInStandardResidues;
        }
    }

    return numAtomsInStandardResidues;
}



rg_INT  Molecule::getNumberOfAtomsInNonStandardResidues() const
{
    rg_INT numAtomsInNonStandardResidues = 0;

    m_atoms.reset4Loop();
    while (m_atoms.setNext4Loop()) {
        Atom* currAtom = m_atoms.getpEntity();
        Residue* currRes = currAtom->getResidue();
        if (currRes == rg_NULL) continue;

        if (!currRes->isStandardResidue()) {
            ++numAtomsInNonStandardResidues;
        }
    }

    return numAtomsInNonStandardResidues;
}



rg_INT  Molecule::getNumberOfStandardResidues() const
{
    rg_INT numStandardResidues = 0;

    m_residues.reset4Loop();
    while (m_residues.setNext4Loop()) {
        Residue* currRes = m_residues.getpEntity();

        if (currRes->isStandardResidue()) {
            ++numStandardResidues;
        }
    }

    return numStandardResidues;
}



rg_INT  Molecule::getNumberOfNonStandardResidues() const
{
    rg_INT numNonStandardResidues = 0;

    m_residues.reset4Loop();
    while (m_residues.setNext4Loop()) {
        Residue* currRes = m_residues.getpEntity();

        if (currRes->getResidueCode() == HOH_RESIDUE) {
            continue;
        }
            
        if (!currRes->isStandardResidue()) {
            ++numNonStandardResidues;
        }
    }

    return numNonStandardResidues;
}



rg_INT  Molecule::getAtomsInStandardResidues(list<Atom*>& atomsInStandardResidues) const
{
    rg_INT numAtomsInStandardResidues = 0;

    m_atoms.reset4Loop();
    while (m_atoms.setNext4Loop()) {
        Atom* currAtom = m_atoms.getpEntity();
        Residue* currRes = currAtom->getResidue();
        if (currRes == rg_NULL) continue;

        if (currRes->isStandardResidue()) {
            atomsInStandardResidues.push_back(currAtom);
            ++numAtomsInStandardResidues;
        }
    }

    return numAtomsInStandardResidues;
}



rg_INT  Molecule::getAtomsInNonStandardResidues(list<Atom*>& atomsInNonStandardResidues) const
{
    rg_INT numAtomsInNonStandardResidues = 0;

    m_atoms.reset4Loop();
    while (m_atoms.setNext4Loop()) {
        Atom* currAtom = m_atoms.getpEntity();
        Residue* currRes = currAtom->getResidue();
        if (currRes == rg_NULL) continue;

        if (!currRes->isStandardResidue()) {
            atomsInNonStandardResidues.push_back(currAtom);
            ++numAtomsInNonStandardResidues;
        }
    }

    return numAtomsInNonStandardResidues;
}



rg_INT  Molecule::evaluateAtomFrequency(map<AtomCode, int>& atomFrequency) const
{
    m_atoms.reset4Loop();
    while (m_atoms.setNext4Loop()) {
        Atom* currAtom = m_atoms.getpEntity();

        if (currAtom->getResidue()->getResidueCode() == HOH_RESIDUE) {
            continue;
        }

        map<AtomCode, int>::iterator it_atom = atomFrequency.find(currAtom->getAtomCode());
        if (it_atom == atomFrequency.end()) {
            atomFrequency.insert(make_pair(currAtom->getAtomCode(), 1));
        }
        else {
            ++it_atom->second;
        }
    }

    return atomFrequency.size();
}



rg_INT  Molecule::evaluateAtomFrequency(map<string, int>& atomFrequency) const
{
    map<AtomCode, int> atomFrequencyByCode;
    evaluateAtomFrequency(atomFrequencyByCode);

    for (map<AtomCode, int>::iterator i_atom = atomFrequencyByCode.begin(); i_atom != atomFrequencyByCode.end(); ++i_atom) {
        string atomSymbol = ATOM_FEATURES[i_atom->first].symbol;
        atomFrequency.insert(make_pair(ATOM_FEATURES[i_atom->first].symbol, i_atom->second));
    }

    return atomFrequency.size();
}



// by Y.Cho at 2011-10-12
string Molecule::getTimeStamp() const
{
    return m_timeStamp;
}


// by Y.Cho at 2011-10-12
rg_INT Molecule::getFileSize() const
{
    return m_fileSize;
}



void Molecule::getConformation(ProteinConformation& proteinConformation)
{
    moveOxygenInCarboxylGroupToBackbone();
    // backbone conformation    
    Backbone backbone;
    getBackboneSortedByChainIDSeqNumber(backbone);
    BackboneConformation backboneConformation;
    backbone.getBackboneConformation(backboneConformation);

    // side-chain conformation
    Rotamer* rotamersOderedByChainIDSeqNumber = rg_NULL;
    rg_INT numRotamers = getRotamersCorrToSidechainsOfAllAminoResidueInstancesOrderedByChainIDSeqNumber(rotamersOderedByChainIDSeqNumber);
    SidechainConformation sidechainConformation(rotamersOderedByChainIDSeqNumber, numRotamers);

    proteinConformation.set(backboneConformation, sidechainConformation);

    if(rotamersOderedByChainIDSeqNumber != rg_NULL)
        delete [] rotamersOderedByChainIDSeqNumber;
}


void Molecule::getBackboneSortedByChainIDSeqNumber(Backbone& backbone)
{
	rg_dList<Atom*> atomsOnBackbone;
	getAtomsOnBackboneSortedByChainIDSequenceNumber(atomsOnBackbone);
	rg_INT numAtoms = atomsOnBackbone.getSize();

	Atom** backboneAtoms = new Atom*[ numAtoms ];
	rg_INDEX i = 0;

	atomsOnBackbone.reset4Loop();
	while (atomsOnBackbone.setNext4Loop())
	{
		Atom* currAtom = atomsOnBackbone.getEntity();
		backboneAtoms[ i++ ] = currAtom;
	}

	backbone.setBackboneAtoms(backboneAtoms, numAtoms);

	if(backboneAtoms != rg_NULL)
		delete[] backboneAtoms;
}

void Molecule::getRotamersCorrToSidechainsOfAllResidueInstancesOrderedByChainIDSeqNumber(ManagerOfRotamers_Assigned_At_Residues& managerOfRotamers_Assigned_At_Residues)
{
    //Residue** residuesOderedByChainIDSeqNumber = rg_NULL;
    //rg_INT numResidues = getAminoResiduesOrderedByChainIDSequenceNumber(residuesOderedByChainIDSeqNumber);

    //Rotamer* rotamersOderedByChainIDSeqNumber = rg_NULL;
    //rotamersOderedByChainIDSeqNumber = new Rotamer[ numResidues ];

    //rg_INDEX i;
    //for (i = 0;i < numResidues;i++)
    //{
    //    rotamersOderedByChainIDSeqNumber[ i ].setResidue( residuesOderedByChainIDSeqNumber[ i ] );
    //}

    Rotamer* rotamersOderedByChainIDSeqNumber = rg_NULL;
    rg_INT numResidues = 0;
    numResidues = getRotamersCorrToSidechainsOfAllAminoResidueInstancesOrderedByChainIDSeqNumber(rotamersOderedByChainIDSeqNumber);
    managerOfRotamers_Assigned_At_Residues.setOrderedRotamers(rotamersOderedByChainIDSeqNumber, numResidues);

    //if (residuesOderedByChainIDSeqNumber != rg_NULL)
    //{
    //    delete [] residuesOderedByChainIDSeqNumber;
    //}

    if(rotamersOderedByChainIDSeqNumber != rg_NULL)
    {
        delete [] rotamersOderedByChainIDSeqNumber;
    }
}

rg_INT Molecule::getVirtualRotamersCorrToSidechainsOfAllNonAminoResidueInstancesOrderedByChainIDSeqNumber(Rotamer*& rotamersOderedByChainIDSeqNumber)
{
	Residue** nonAminoResiduesOderedByChainIDSeqNumber = rg_NULL;
	rg_INT numResidues = getNonAminoResiduesOrderedByChainIDSequenceNumber(nonAminoResiduesOderedByChainIDSeqNumber);

	rg_INT numVirtualRotamers = numResidues;

	rotamersOderedByChainIDSeqNumber = rg_NULL;
	rotamersOderedByChainIDSeqNumber = new Rotamer[numVirtualRotamers];

	rg_INDEX i;
	for (i = 0; i < numVirtualRotamers; i++)
	{
		//rotamersOderedByChainIDSeqNumber[ i ].setResidue( residuesOderedByChainIDSeqNumber[ i ] );
		rotamersOderedByChainIDSeqNumber[i].setResidueWithDihedralAngles(nonAminoResiduesOderedByChainIDSeqNumber[i]);
	}

	if (nonAminoResiduesOderedByChainIDSeqNumber != rg_NULL)
	{
		delete[] nonAminoResiduesOderedByChainIDSeqNumber;
	}

	return numVirtualRotamers;
}

rg_INT Molecule::getRotamersCorrToSidechainsOfAllAminoResidueInstancesOrderedByChainIDSeqNumber(Rotamer*& rotamersOderedByChainIDSeqNumber)
{
	Residue** residuesOderedByChainIDSeqNumber = rg_NULL;
	rg_INT numResidues = getAminoResiduesOrderedByChainIDSequenceNumber(residuesOderedByChainIDSeqNumber);

	rg_INT numRotamers = numResidues;

	rotamersOderedByChainIDSeqNumber = rg_NULL;
	rotamersOderedByChainIDSeqNumber = new Rotamer[ numRotamers ];

	rg_INDEX i;
	for (i = 0;i < numRotamers;i++)
	{
		//rotamersOderedByChainIDSeqNumber[ i ].setResidue( residuesOderedByChainIDSeqNumber[ i ] );
		rotamersOderedByChainIDSeqNumber[ i ].setResidueWithDihedralAngles( residuesOderedByChainIDSeqNumber[ i ] );
	}

	if (residuesOderedByChainIDSeqNumber != rg_NULL)
	{
		delete [] residuesOderedByChainIDSeqNumber;
	}

	return numRotamers;
}

rg_INT Molecule::getRotamersCorrToSidechainsOfAllResidueInstancesOrderedByChainIDSeqNumber(Rotamer *& rotamersOderedByChainIDSeqNumber)
{
	Residue** residuesOderedByChainIDSeqNumber = rg_NULL;
	rg_INT numResidues = getAllResiduesOrderedByChainIDSequenceNumber(residuesOderedByChainIDSeqNumber);

	rg_INT numRotamers = numResidues;

	rotamersOderedByChainIDSeqNumber = rg_NULL;
	rotamersOderedByChainIDSeqNumber = new Rotamer[numRotamers];

	rg_INDEX i;
	for (i = 0; i < numRotamers; i++)
	{
		rotamersOderedByChainIDSeqNumber[i].setResidueWithDihedralAngles(residuesOderedByChainIDSeqNumber[i]);
	}

	if (residuesOderedByChainIDSeqNumber != rg_NULL)
	{
		delete[] residuesOderedByChainIDSeqNumber;
	}

	return numRotamers;
}

list<PharmaFeature>* Molecule::getListPharmaFeatures()
{
   return &m_pharmaFeatures;
}



///////////////////////////////////////////////////////////////////////////////
//
//  SET FUNCTION
//
void Molecule::setMoleculeFileName( const string& moleculeFileName )
{
    m_moleculeFileName = moleculeFileName;
}



void Molecule::setAtomRadiusType( const AtomRadiusType& atomRadiusType )
{
    m_atomRadiusType = atomRadiusType;
}



string* Molecule::addHeaderRecords( const string& recLine )
{
    return m_headerRecords.addTail( recLine );
}



Atom* Molecule::addAtom( const Atom& atom )
{
    return m_atoms.addTail( atom );
}



ChemicalBond* Molecule::addChemicalBond( const ChemicalBond& chemicalBond )
{
    return m_chemicalBonds.addTail( chemicalBond );
}



Residue* Molecule::addResidue( const Residue& residue )
{
    return m_residues.addTail( residue );
}



Chain* Molecule::addChain( const Chain& aChain )
{
    return m_chains.addTail( aChain );
}



void Molecule::setCenterOfMass( const rg_Point3D& centerOfMass )
{
    m_centerOfMass = centerOfMass;
}



void Molecule::setMinEnclosingSphere( const Sphere& minEnclosingSphere )
{
    m_minEnclosingSphere = minEnclosingSphere;
}



void Molecule::setModelSerialFromInputFile( const rg_INT& modelSerialFromInputFile )
{
    m_modelSerialFromInputFile = modelSerialFromInputFile;
}



// by Y.Cho at 2011-10-12
void Molecule::setTimeStamp( const string& timeStamp )
{
    m_timeStamp = timeStamp;
}



// by Y.Cho at 2011-10-12
void Molecule::setFileSize(  const rg_INT& filesize )
{
    m_fileSize = filesize;
}





void Molecule::moveOxygenInCarboxylGroupToBackbone()
{
    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        Residue* currResidue = m_residues.getpEntity();

        if ( currResidue->isAminoResidue() == rg_FALSE )
            continue;
        
        rg_dList<Atom*>* atomsOnResidue = currResidue->getAtoms();
        while ( atomsOnResidue->setNext4Loop() ) {
            Atom* currAtom = atomsOnResidue->getEntity();
        
            if( currAtom->getpChemicalProperties()->isOnBackBone() == rg_FALSE ) {
                if ( currAtom->getAtomNameFromInputFile() == " O  ") {
                    currAtom->getpChemicalProperties()->setIsOnBackBone( rg_TRUE );
                    break;
                }
            }
        }
    }
}



void Molecule::deleteAtomsInSideChains( rg_BOOL replaceResidueCodeToUnknown, const RemoteIndicator& startRemoteIndicator, const ResidueType& typeOfResidue )
{
    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        Residue* currResidue = m_residues.getpEntity();
        
        if ( currResidue->isResidueType( typeOfResidue ) )
            deleteAtomsInSideChain( replaceResidueCodeToUnknown, startRemoteIndicator, currResidue );
    }
}



void Molecule::deleteAtomsInSideChain( rg_BOOL replaceResidueCodeToUnknown, const RemoteIndicator& startRemoteIndicator, Residue* targetResidue )
{
    if ( !targetResidue->isStandardResidue() )
        return;
    
    rg_INT startRemote = (rg_INT) startRemoteIndicator;

    rg_dList<Atom*>* atomsInResidue = targetResidue->getAtoms();
    atomsInResidue->reset4Loop();
    while ( atomsInResidue->setNext4Loop() ) {
        Atom* currAtom = atomsInResidue->getEntity();
        rg_INT i_remote = (rg_INT) currAtom->getpChemicalProperties()->getRemoteIndicator();

        if ( !currAtom->getpChemicalProperties()->isOnBackBone() && i_remote >= startRemote )
            deleteAtom( replaceResidueCodeToUnknown, currAtom );
    }
}



void Molecule::deleteAtom( rg_BOOL replaceResidueCodeToUnknown, Atom* atomToDelete )
{
    // 1. Delete from AtomList in Residue
    Residue* residueOfAtom = atomToDelete->getResidue();
    rg_dList<Atom*>* atomsInResidue = residueOfAtom->getAtoms();
    atomsInResidue->kill( atomToDelete );

    // 2. Change the type of ResidueCode to "UNK_RESIDUE"   
    if( replaceResidueCodeToUnknown ) {
        residueOfAtom->setResidueCode( UNK_RESIDUE );
        residueOfAtom->setResidueName( "UNK" );
    }

    // 3. Change the chemical bond information
    rg_dList<ChemicalBond*>* bondsOfAtomToDelete = atomToDelete->getListChemicalBond();
    bondsOfAtomToDelete->reset4Loop();
    while ( bondsOfAtomToDelete->setNext4Loop() ) {
        ChemicalBond* currBond = bondsOfAtomToDelete->getEntity();
        Atom* bondedAtom = currBond->getBondedAtom( atomToDelete );
        
        // 3.1 Delete ChemicalBond from bondList in bondedAtom
        rg_dList<ChemicalBond*>* bondsOfBondedAtom = bondedAtom->getListChemicalBond();
        bondsOfBondedAtom->reset4Loop();
        while ( bondsOfBondedAtom->setNext4Loop() ) {
            ChemicalBond* currrBondOfBondedAtom = bondsOfBondedAtom->getEntity();
            
            if ( currrBondOfBondedAtom->getBondedAtom( bondedAtom ) == atomToDelete ) {
                bondsOfBondedAtom->killCurrent();
                break;
            }
        }
        
        // 3.2 Delete currBond from m_chemicalBonds
        m_chemicalBonds.kill( currBond );
    }

    // 4. Delete from m_atoms
    m_atoms.kill( atomToDelete );
    
    // 5. Reset ID of atoms : to synchronize with ParticleModel.
    rg_INT newAtomID = 0;
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        m_atoms.getpEntity()->setID( newAtomID );
        newAtomID++;
    }

    // 6. Reset ID of chemicalBonds
    rg_INT newBondID = 0;
    m_chemicalBonds.reset4Loop();
    while ( m_chemicalBonds.setNext4Loop() ) {
        m_chemicalBonds.getpEntity()->setID( newBondID );
        newBondID++;
    }
}



void Molecule::recoverMissingAtomsInResiduesWithoutCoordinates( const ResidueType& typeOfResidue )
{
    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        Residue* currResidue = m_residues.getpEntity();
        
        if ( currResidue->isResidueType( typeOfResidue ) )
            recoverMissingAtomsInResidueWithoutCoordinates( currResidue );
    }
}



void Molecule::recoverMissingAtomsInResidueWithoutCoordinates( Residue* targetResidue )
{
    if ( !targetResidue->isStandardResidue() ) {
        return;
    }

    typedef map<string, Atom*> AtomNameMap;

    //// Initialize map of atom names ( key: (string)nameOfAtom, value: Atom* )
    //
    AtomNameMap mapOfAtomName;
    mapOfAtomName.clear();
    
    rg_dList<Atom*>* atomsInResidue = targetResidue->getAtoms();
    
    atomsInResidue->reset4Loop();
    while( atomsInResidue->setNext4Loop() ) {
        Atom* currAtom = atomsInResidue->getEntity();
        mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
    }

    AtomSymbolMap mapOfAtomSymbol;
    mapOfAtomSymbol.clear();
    
    for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        mapOfAtomSymbol.insert( AtomSymbolMap::value_type(atomSymbol, i_atomCode) );
    }
    //
    ////


    ResidueCode aResidueCode = targetResidue->getResidueCode();

    rg_INT i_residueCode = (rg_INT)targetResidue->getResidueCode();

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
            else { // atom of "atomName[i_atomPair]" is missing : need to create Atom 
                if ( isHydrogenAtomName( atomName[i_atomPair] ) == rg_FALSE  ) {    // Hydrogen atoms will not recover'
                    
                    if ( strcmp(atomName[i_atomPair].c_str(), " OXT") == 0 )
                        continue;

                    rg_INT newID = m_atoms.getSize();
                    Atom*  currAtom  = addAtom( Atom(newID) );
                    targetResidue->addAtom( currAtom );
                    currAtom->setResidue( targetResidue );
                    currAtom->setSerialFromInputFile( getBiggestAtomSerialFromInputFile()+1 );
                    currAtom->setAtomNameFromInputFile( atomName[i_atomPair] );
                    MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( mapOfAtomSymbol, atomName[i_atomPair], currAtom );
                    currAtom->setAtomBall( Sphere( 9999.999, 9999.999, 9999.999, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );
                    mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
                    atomPtr[i_atomPair] = currAtom;
                }
            }
        }

        if( atomPtr[0] == rg_NULL || atomPtr[1] == rg_NULL )
            continue;
        

        //// Recover chemical bonds
        //
        
        ChemicalBond newChemicalBond( m_chemicalBonds.getSize(), atomPtr[0], atomPtr[1] );
    
        ChemicalBond* existingChemBondA = rg_NULL;
        ChemicalBond* existingChemBondB = rg_NULL;
    
    
        rg_FLAG isNewChemBondExistInChemBondListInAtomA = MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, atomPtr[0], existingChemBondA );
        rg_FLAG isNewChemBondExistInChemBondListInAtomB = MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, atomPtr[1], existingChemBondB );
    
    
        //// If "existingChemBondA" is not "rg_NULL" then "existingChemBondB" MUST! not be "rg_NULL".
        //// "existingChemBondA" == "existingChemBondB" --> ALWAYS!!!
        //// But equal comparison between "existingChemBondA" and "existingChemBondB" is performed for bug test.
    
        if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
            ChemicalBond* pNewChemicalBond = addChemicalBond( newChemicalBond );
            atomPtr[0]->addChemicalBond( pNewChemicalBond );
            atomPtr[1]->addChemicalBond( pNewChemicalBond );
        }
        else if( isNewChemBondExistInChemBondListInAtomA == rg_TRUE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
            atomPtr[1]->addChemicalBond( existingChemBondA );
        }
        else if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_TRUE ) {
            atomPtr[0]->addChemicalBond( existingChemBondB );
        }
        //
        ////
    }

}



rg_BOOL Molecule::isHydrogenAtomName( const string& atomName )
{
    size_t found_A  = atomName.find(" H");
    size_t found_B  = atomName.find("HH");
    size_t found_C  = atomName.find("HD");
    size_t found_D  = atomName.find("HE");

    if ( found_A != string::npos || found_B != string::npos ||
         found_C != string::npos || found_D != string::npos )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_INT Molecule::getBiggestAtomSerialFromInputFile()
{
    rg_INT biggestAtomSerial = 0;

    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        rg_INT currSerial = m_atoms.getpEntity()->getSerialFromInputFile();
        if( currSerial > biggestAtomSerial )
            biggestAtomSerial = currSerial;
    }
    return biggestAtomSerial;
}

void Molecule::recoverMissingAtomsInResidueWithGivenCoordinates( Residue* targetResidue, 
	                                                             const rg_Point3D& defaultAtomCenter, 
	                                                             AtomSymbolMap& mapOfAtomSymbols )
{
	if ( !targetResidue->isStandardResidue() ) 
	{
		return;
	}

	// get start and end index of bonded atom pairs in residue
	ResidueCode code = targetResidue->getResidueCode();
	rg_INDEX start_index = 0;	
	rg_INDEX end_index   = RESIDUE_FEATURES[code].numOfBonds - 1;

	// recover missing atoms
	recoverMissingAtomsInResidueWithGivenCoordinates(targetResidue,
		                                             defaultAtomCenter,
		                                             start_index,
		                                             end_index,
		                                             mapOfAtomSymbols);
	//// recover backbone atoms
	//recoverMissingAtomsOnBackboneOfResidueWithGivenCoordinates(targetResidue,
	//	                                                       defaultAtomCenter,
	//														   mapOfAtomSymbols);
	//// recover sidechain atoms
	//recoverMissingAtomsInSidechainOfResidueWithGivenCoordinates(targetResidue,
	//	                                                       defaultAtomCenter,
	//	                                                       mapOfAtomSymbols);
	//if ( !targetResidue->isStandardResidue() ) 
	//{
	//	return;
	//}

	//// Initialize map of atom names ( key: (string)nameOfAtom, value: Atom* )
	//typedef map<string, Atom*> AtomNameMap;
	//AtomNameMap mapOfAtomName;
	//mapOfAtomName.clear();

	//rg_dList<Atom*>* atomsInResidue = targetResidue->getAtoms();

	//atomsInResidue->reset4Loop();
	//while( atomsInResidue->setNext4Loop() ) 
	//{
	//	Atom* currAtom = atomsInResidue->getEntity();
	//	mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
	//}
	//////////////////////////////////////////////////////////////////////////

	//rg_INT i_residueCode = (rg_INT)targetResidue->getResidueCode();

	//for( rg_INT i_bond=0; i_bond<RESIDUE_FEATURES[i_residueCode].numOfBonds; i_bond++ )  
	//{
	//	string atomName[ 2 ] = { RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][0], 
	//		                     RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][1] };
	//	Atom*  atomPtr[ 2 ]  = { rg_NULL,    rg_NULL    };

	//	AtomNameMap::iterator atomNameMap_i;
	//	for( rg_INT i_atomPair=0; i_atomPair<2; i_atomPair++ ) 
	//	{
	//		atomNameMap_i = mapOfAtomName.find( atomName[i_atomPair] );

	//		if ( atomNameMap_i != mapOfAtomName.end() ) { 
	//			atomPtr[i_atomPair] = (*atomNameMap_i).second;
	//		}
	//		else 
	//		{ // atom of "atomName[i_atomPair]" is missing : need to create Atom 
	//			//if ( doesStringContainHydrogenAtomName( atomName[i_atomPair] )( atomName[i_atomPair] ) == rg_FALSE  ) 
	//			//{    // Hydrogen atoms will not recover'

	//				if ( strcmp(atomName[i_atomPair].c_str(), " OXT") == 0 )
	//					continue;

	//				rg_INT newID = m_atoms.getSize();
	//				Atom*  currAtom  = addAtom( Atom(newID) );
	//				targetResidue->addAtom( currAtom );
	//				currAtom->setResidue( targetResidue );
	//				currAtom->setSerialFromInputFile( getBiggestAtomSerialFromInputFile()+1 );
	//				currAtom->setAtomNameFromInputFile( atomName[i_atomPair] );
	//				MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( mapOfAtomSymbols, atomName[i_atomPair], currAtom );
	//				currAtom->setAtomBall( Sphere( defaultAtomCenter, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );
	//				mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
	//				atomPtr[i_atomPair] = currAtom;
	//			//}
	//		}
	//	}

	//	if( atomPtr[0] == rg_NULL || atomPtr[1] == rg_NULL )
	//		continue;

	//	// Recover chemical bonds
	//	ChemicalBond newChemicalBond( m_chemicalBonds.getSize(), atomPtr[0], atomPtr[1] );

	//	ChemicalBond* existingChemBondA = rg_NULL;
	//	ChemicalBond* existingChemBondB = rg_NULL;


	//	rg_FLAG isNewChemBondExistInChemBondListInAtomA = MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, atomPtr[0], existingChemBondA );
	//	rg_FLAG isNewChemBondExistInChemBondListInAtomB = MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, atomPtr[1], existingChemBondB );


	//	//// If "existingChemBondA" is not "rg_NULL" then "existingChemBondB" MUST! not be "rg_NULL".
	//	//// "existingChemBondA" == "existingChemBondB" --> ALWAYS!!!
	//	//// But equal comparison between "existingChemBondA" and "existingChemBondB" is performed for bug test.

	//	if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
	//		ChemicalBond* pNewChemicalBond = addChemicalBond( newChemicalBond );
	//		atomPtr[0]->addChemicalBond( pNewChemicalBond );
	//		atomPtr[1]->addChemicalBond( pNewChemicalBond );
	//	}
	//	else if( isNewChemBondExistInChemBondListInAtomA == rg_TRUE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
	//		atomPtr[1]->addChemicalBond( existingChemBondA );
	//	}
	//	else if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_TRUE ) {
	//		atomPtr[0]->addChemicalBond( existingChemBondB );
	//	}
	//	////////////////////////////
	//}
}

void Molecule::recoverMissingAtomsOnBackboneOfResidueWithGivenCoordinates( Residue* targetResidue, 
	                                                                       const rg_Point3D& defaultAtomCenter, 
	                                                                       AtomSymbolMap& mapOfAtomSymbols )
{
	if ( !targetResidue->isStandardResidue() ) 
	{
		return;
	}

	// get start and end index of bonded atom pairs in sidechain
	ResidueCode code = targetResidue->getResidueCode();
	rg_INDEX start_index = -1;
	rg_INDEX end_index   = -1;

	if( code == PRO_AMINO_RESIDUE || 
		code == PR0_AMINO_RESIDUE ||
		code == PRZ_AMINO_RESIDUE   )
	{
		start_index = 0;
		end_index   = 5;
	}
	else
	{
		start_index = 0;
		end_index   = 7;
	}

	// recover missing atoms
	recoverMissingAtomsInResidueWithGivenCoordinates(targetResidue,
		                                             defaultAtomCenter,
		                                             start_index,
		                                             end_index,
		                                             mapOfAtomSymbols);
}


void Molecule::recoverMissingAtomsInSidechainOfResidueWithGivenCoordinates( Residue* targetResidue, 
	                                                                        const rg_Point3D& defaultAtomCenter, 
	                                                                        AtomSymbolMap& mapOfAtomSymbols )
{
	if ( !targetResidue->isStandardResidue() ) 
	{
		return;
	}

	// get start and end index of bonded atom pairs in sidechain
	ResidueCode code = targetResidue->getResidueCode();
	rg_INDEX start_index = -1;	
	rg_INDEX end_index = RESIDUE_FEATURES[code].numOfBonds - 1;
	
	if( code == PRO_AMINO_RESIDUE || 
		code == PR0_AMINO_RESIDUE ||
		code == PRZ_AMINO_RESIDUE   )
	{
		start_index = 5;		
	}
	else
	{
		start_index = 7;
	}

	// recover missing atoms
	recoverMissingAtomsInResidueWithGivenCoordinates(targetResidue,
		                                             defaultAtomCenter,
													 start_index,
													 end_index,
													 mapOfAtomSymbols);
}

void Molecule::recoverMissingAtomsInResidueWithGivenCoordinates( Residue* targetResidue, 
	                                                             const rg_Point3D& defaultAtomCenter,
	                                                             const rg_INT& startIndexForPairOfBondedAtom,
	                                                             const rg_INT& endIndexForPairOfBondedAtom,
	                                                             AtomSymbolMap& mapOfAtomSymbols )
{
	// Initialize map of atom names ( key: (string)nameOfAtom, value: Atom* )
	typedef map<string, Atom*> AtomNameMap;
	AtomNameMap mapOfAtomName;
	mapOfAtomName.clear();

	rg_dList<Atom*>* atomsInResidue = targetResidue->getAtoms();

	atomsInResidue->reset4Loop();
	while( atomsInResidue->setNext4Loop() ) 
	{
		Atom* currAtom = atomsInResidue->getEntity();
		mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
	}
	////////////////////////////////////////////////////////////////////////

	rg_INT i_residueCode = (rg_INT)targetResidue->getResidueCode();

	for( rg_INT i_bond=startIndexForPairOfBondedAtom; i_bond<=endIndexForPairOfBondedAtom; i_bond++ )  
	{
		string atomName[ 2 ] = { RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][0], 
			RESIDUE_FEATURES[i_residueCode].pairOfBondedAtomName[i_bond][1] };
		Atom*  atomPtr[ 2 ]  = { rg_NULL,    rg_NULL    };

		AtomNameMap::iterator atomNameMap_i;
		for( rg_INT i_atomPair=0; i_atomPair<2; i_atomPair++ ) 
		{
			atomNameMap_i = mapOfAtomName.find( atomName[i_atomPair] );

			if ( atomNameMap_i != mapOfAtomName.end() ) { 
				atomPtr[i_atomPair] = (*atomNameMap_i).second;
			}
			else 
			{ // atom of "atomName[i_atomPair]" is missing : need to create Atom 
				//if ( doesStringContainHydrogenAtomName( atomName[i_atomPair] )( atomName[i_atomPair] ) == rg_FALSE  ) 
				//{    // Hydrogen atoms will not recover'

				if ( strcmp(atomName[i_atomPair].c_str(), " OXT") == 0 )
					continue;

				rg_INT newID = m_atoms.getSize();
				Atom*  currAtom  = addAtom( Atom(newID) );
				targetResidue->addAtom( currAtom );
				currAtom->setResidue( targetResidue );
				currAtom->setSerialFromInputFile( getBiggestAtomSerialFromInputFile()+1 );
				currAtom->setAtomNameFromInputFile( atomName[i_atomPair] );
				MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( mapOfAtomSymbols, atomName[i_atomPair], currAtom );
				currAtom->setAtomBall( Sphere( defaultAtomCenter, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );
				mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
				atomPtr[i_atomPair] = currAtom;
				//}
			}
		}

		if( atomPtr[0] == rg_NULL || atomPtr[1] == rg_NULL )
			continue;

		// Recover chemical bonds
		ChemicalBond newChemicalBond( m_chemicalBonds.getSize(), atomPtr[0], atomPtr[1] );

		ChemicalBond* existingChemBondA = rg_NULL;
		ChemicalBond* existingChemBondB = rg_NULL;


		rg_FLAG isNewChemBondExistInChemBondListInAtomA = MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, atomPtr[0], existingChemBondA );
		rg_FLAG isNewChemBondExistInChemBondListInAtomB = MoleculeIOFunctions::isChemicalBondExistInChemicalBondListInAtom( &newChemicalBond, atomPtr[1], existingChemBondB );


		//// If "existingChemBondA" is not "rg_NULL" then "existingChemBondB" MUST! not be "rg_NULL".
		//// "existingChemBondA" == "existingChemBondB" --> ALWAYS!!!
		//// But equal comparison between "existingChemBondA" and "existingChemBondB" is performed for bug test.

		if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
			ChemicalBond* pNewChemicalBond = addChemicalBond( newChemicalBond );
			atomPtr[0]->addChemicalBond( pNewChemicalBond );
			atomPtr[1]->addChemicalBond( pNewChemicalBond );
		}
		else if( isNewChemBondExistInChemBondListInAtomA == rg_TRUE && isNewChemBondExistInChemBondListInAtomB == rg_FALSE ) {
			atomPtr[1]->addChemicalBond( existingChemBondA );
		}
		else if( isNewChemBondExistInChemBondListInAtomA == rg_FALSE && isNewChemBondExistInChemBondListInAtomB == rg_TRUE ) {
			atomPtr[0]->addChemicalBond( existingChemBondB );
		}
		////////////////////////////
	}
}


void Molecule::resetIDs()
{
	resetSerialsFromInputFileOfAtoms();
}

void Molecule::resetSerialsFromInputFileOfAtoms()
{
	rg_INT newAtomSerial = 1;
	m_atoms.reset4Loop();
	while ( m_atoms.setNext4Loop() ) {
		m_atoms.getpEntity()->setSerialFromInputFile( newAtomSerial );
		newAtomSerial++;
	}
}


void Molecule::addOxygenOfCarboxylGroupWithCenter(Atom* carbonInCarboxylGroup, 
	                                              const rg_Point3D& oxygenCenter, 
												  Residue* targetResidue,
												  AtomSymbolMap& atomSymbolMap)
{
	// add new oxygen and set its properties
	rg_INT newID = m_atoms.getSize();
	Atom*  oxygenOfCarboxylGroup  = addAtom( Atom(newID) );
	targetResidue->addAtom( oxygenOfCarboxylGroup );
	oxygenOfCarboxylGroup->setResidue( targetResidue );
	oxygenOfCarboxylGroup->setSerialFromInputFile( getBiggestAtomSerialFromInputFile()+1 );
	string oxygenName(" O  ");
	oxygenOfCarboxylGroup->setAtomNameFromInputFile( oxygenName );
	MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile( atomSymbolMap, oxygenName, oxygenOfCarboxylGroup );
	oxygenOfCarboxylGroup->setAtomBall( Sphere( oxygenCenter, ATOM_FEATURES[oxygenOfCarboxylGroup->getAtomCode()].radius ) );

	// set chemical bonds
	MoleculeIOFunctions::setChemicalBondToMoleculeWithoutDuplication( carbonInCarboxylGroup, oxygenOfCarboxylGroup, *this);
}

Residue* Molecule::addNewRotamerWithFreezedBetaCarbon(Residue* prevResidue,
								                      Residue* nextResidue,
								                      const ResidueCode& residueCode,
								                      Chain*   chain,
								                      Residue* refResidue)
{
	// make new rotamer without specific atom coordinates
	Residue* targetResidue = 
		addNewResidueWithoutSpecificationOfAtomCoordinates(prevResidue, 
		                                                   nextResidue, 
														   residueCode, 
														   chain);
	// copy the coordinates of N, CA, C, O, CB
	rg_dList<Atom*> atomsOfRefResidue;
	atomsOfRefResidue.add(refResidue->getNitrogenInAminoGroupOfAminoResidue());
	atomsOfRefResidue.add(refResidue->getAlphaCarbonOfAminoResidue());
	atomsOfRefResidue.add(refResidue->getCarbonInCarboxylGroupOfAminoResidue());
	atomsOfRefResidue.add(refResidue->getOxygenInCarboxylGroupOfAminoResidue());
	atomsOfRefResidue.add(refResidue->getBetaCarbonInSideChainOfAminoResidue());

	rg_dList<Atom*>* freezedAtoms = targetResidue->getAtoms();
	freezedAtoms->reset4Loop();
	atomsOfRefResidue.reset4Loop();
	while (freezedAtoms->setNext4Loop() &&
		   atomsOfRefResidue.setNext4Loop())
	{
		Atom* currAtomOfTargetResidue = freezedAtoms->getEntity();
		Atom* currAtomOfRefResidue = atomsOfRefResidue.getEntity();
		currAtomOfTargetResidue->setAtomBall( currAtomOfRefResidue->getAtomBall() );
	}
	return targetResidue;
}

Residue* Molecule::addNewRotamerWithFreezedBetaCarbon(Residue* prevResidue,
								                      Residue* nextResidue,
								                      const ResidueCode& residueCode,
								                      Chain*   chain,
								                      rg_dList<Atom*>& freezedAtoms)
{
	// make residue and atoms
	Residue* targetResidue = 
		makeAminoResidueNConstituentAtoms(residueCode, chain);

	// copy the coordinates of N, CA, C, O, CB
	rg_INDEX index = 0;	
	rg_dList<Atom*>* atomsOfTargetResidue = targetResidue->getAtoms();
	freezedAtoms.reset4Loop();
	atomsOfTargetResidue->reset4Loop();
	while (freezedAtoms.setNext4Loop() && 
		   atomsOfTargetResidue->setNext4Loop() && 
		   index < 5)
	{
		Atom* currAtom = freezedAtoms.getEntity();
		Atom* currAtomInTarget = atomsOfTargetResidue->getEntity();
		currAtomInTarget->setAtomBall( currAtom->getAtomBall() );
		index++;
	}

	// disconnect chemical bond between the first and second residue
	disconnectBondBetween(prevResidue, nextResidue);

	// connect chemical bonds
	// (1) bond between the first and target residues
	MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(prevResidue, targetResidue, *this);
	// (2) bond between the second and target residues
	MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(targetResidue, nextResidue, *this);

	return targetResidue;
}


Residue* Molecule::addNewResidueWithoutSpecificationOfAtomCoordinates(Residue* prevResidue,
							                                          Residue* nextResidue,
								                                      const ResidueCode& residueCode,
								                                      Chain*   chain)
{
	if(prevResidue == rg_NULL && nextResidue == rg_NULL)
		return rg_NULL;

	// make residue and atoms
	Residue* targetResidue = 
		makeAminoResidueNConstituentAtoms(residueCode, chain);

	return addNewResidue(prevResidue, nextResidue, targetResidue);
}

Residue*  Molecule::addNewResidue(Residue* prevResidue, Residue* nextResidue,  Residue* newResidue)
{
	if(prevResidue == rg_NULL && nextResidue == rg_NULL)
		return rg_NULL;
	
	if(prevResidue != rg_NULL && nextResidue != rg_NULL)	
	{
		// disconnect chemical bond between previous and next residue
		disconnectBondBetween(prevResidue, nextResidue);
		// connect chemical bonds
		// (1) bond between previous and target residues
		MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(prevResidue, newResidue, *this);
		// (2) bond between next and target residues
		MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(newResidue, nextResidue, *this);		
	}
	// previous residue is the end of this chain
	else if(prevResidue != rg_NULL && nextResidue == rg_NULL)
	{
		// connect chemical bonds
		// (1) bond between previous and target residues
		MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(prevResidue, newResidue, *this);
	}
	// next residue is the start of this chain
	// else if(nextResidue != rg_NULL && prevResidue == rg_NULL)
	else
	{
		// connect chemical bonds
		// (1) bond between next and target residues
		MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(newResidue, nextResidue, *this);
	}
	return newResidue;
}

Residue* Molecule::makeAminoResidueNConstituentAtoms(const ResidueCode& residueCode,
		                                             Chain* chain)
{
    rg_INT minResidueSeq = INT_MAX;
    rg_INT maxResidueSeq = 0;    
    rg_dList<Residue*>* residuesInChain = chain->getResidues();
    
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
	
	// add residue
	rg_INT residueID = m_residues.getSize();
	Residue* targetResidue = m_residues.add(Residue());
	targetResidue->setID(residueID);
	targetResidue->setResidueCode(residueCode);
	targetResidue->setChain(chain);
	targetResidue->setSequenceNumber(maxResidueSeq + 1);
	targetResidue->setResidueName(string(RESIDUE_FEATURES[residueCode].threeCodeName));

	chain->addResidue(targetResidue);

	// make atoms of residue and add them into molecule
	rg_INT startAtomID = m_atoms.getSize();
	rg_INT startAtomSerial = m_atoms.getLastEntity().getSerialFromInputFile() + 1;

	makeAndAddAtomsOfResidue(targetResidue,
		                     startAtomID,
					         startAtomSerial);

	// make chemical bonds for atoms of residue
	MoleculeIOFunctions::setChemicalBondsBetweenAtomsInStdResidue(targetResidue, *this);

	//rg_dList<ChemicalBond>* chemBondsOfMolecule = molecule.getChemicalBonds();
	//rg_INT startChemBondID = chemBondsOfMolecule->getSize();

	// make chemical bonds of residue
// 	makeChemicalBondsOfResidue(targetResidue,
// 		                       atomsOfResidue,
// 							   startChemBondID,
// 							   chemicalBondsOfResidue);
	return targetResidue;
}

void Molecule::makeAndAddAtomsOfResidue(Residue* targetResidue, 
										const rg_INT& startAtomID,
										const rg_INT& startAtomSerial)
{
	// make map for atom names
	AtomSymbolMap mapOfAtomSymbol;
    for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) 
	{
        string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
        mapOfAtomSymbol.insert( AtomSymbolMap::value_type(atomSymbol, i_atomCode) );
    }	

	ResidueCode code = targetResidue->getResidueCode();
	// test
	if(code == ARG_AMINO_RESIDUE)
		int i = 111;

	// collect the types of atoms in target residue
	rg_INT maxNumBondsInResidue = 38;
	
	rg_INDEX start_index = -1;
	rg_INDEX end_index = -1;
	rg_INDEX i;
	for(i = 0;i < maxNumBondsInResidue;i++)
	{
		if(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ] == " CA ")
			start_index = i;
		if(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ] == " OXT")
			end_index = i - 1;
		if(start_index != -1 && end_index != -1)
			break;
	}
	
	//rg_dList<char*> atomNames;
	rg_dList<string> atomNames;
	for(i = start_index;i <= end_index;i++)
	{
		string firstAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]);
		string secondAtomName = string(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]);
		rg_INDEX firstIndex = firstAtomName.find("H");
		rg_INDEX secondIndex = secondAtomName.find("H");

		if(!atomNames.isInList(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]) &&
		   firstIndex != 0 && firstIndex != 1 )
			atomNames.add(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 0 ]);
		if(!atomNames.isInList(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]) &&
		   secondIndex != 0 && secondIndex != 1 )
			atomNames.add(RESIDUE_FEATURES[code].pairOfBondedAtomName[ i ][ 1 ]);
	}

	// make each atom one by one and add that into molecule
	rg_INT atomID = startAtomID;
	rg_INT atomSerial = startAtomSerial;
	atomNames.reset4Loop();
	while (atomNames.setNext4Loop())
	{		
		Atom* currAtom = m_atoms.add( Atom(atomID++) );		
		currAtom->setSerialFromInputFile(atomSerial++);
		//char* currAtomName = atomNames.getEntity();
		string currAtomName = atomNames.getEntity();
		currAtom->setAtomNameFromInputFile( string(currAtomName) );
		currAtom->setResidue(targetResidue);
		targetResidue->addAtom(currAtom);

		MoleculeIOFunctions::setAtomCodeAndChemicalPropertiesFromAtomNameToTargetAtomForPDBFile(mapOfAtomSymbol, currAtom->getAtomNameFromInputFile(), currAtom);
		currAtom->setAtomBall( Sphere( 9999.999, 9999.999, 9999.999, ATOM_FEATURES[currAtom->getAtomCode()].radius ) );		
	}        
}

void Molecule::makeChemicalBondsOfResidue(Residue* targetResidue, 
			                              rg_dList<Atom>& atomsOfResidue,
										  const rg_INT& startChemBondID,
										  rg_dList<ChemicalBond>& chemicalBondsOfResidue)
{
    typedef map<string, Atom*> AtomNameMap;
    AtomNameMap mapOfAtomName;

    atomsOfResidue.reset4Loop();
    while( atomsOfResidue.setNext4Loop() ) 
	{
        Atom* currAtom = atomsOfResidue.getpEntity();
        mapOfAtomName.insert( AtomNameMap::value_type( currAtom->getAtomNameFromInputFile(), currAtom ) );
    }	

	rg_INT bondID = startChemBondID;
    rg_INT i_residueCode = (rg_INT)targetResidue->getResidueCode();    
    for( rg_INT i_bond=0; i_bond<RESIDUE_FEATURES[i_residueCode].numOfBonds; i_bond++ )  
	{
        
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
        
		ChemicalBond newChemicalBond( bondID++, atomPtr[0], atomPtr[1] );
		ChemicalBond* pNewChemicalBond = chemicalBondsOfResidue.add(newChemicalBond);
		atomPtr[0]->addChemicalBond(pNewChemicalBond);
		atomPtr[1]->addChemicalBond(pNewChemicalBond);
    }
}

void Molecule::disconnectBondBetween(Residue* firstResidue, Residue* secondResidue)
{
	if(firstResidue == rg_NULL || secondResidue == rg_NULL)
		return;

    Atom* firstAtom  = rg_NULL;
    Atom* secondAtom = rg_NULL;    

    rg_dList<Atom*>* atomList = firstResidue->getAtoms();

    Atom* tempAtom = rg_NULL;

    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();

        if( tempAtom->getAtomNameFromInputFile() == " C  ") {
            firstAtom = tempAtom;
            break;
        }
    }

    atomList = secondResidue->getAtoms();

    atomList->reset4Loop();
    while( atomList->setNext4Loop() ) {
        tempAtom = atomList->getEntity();
        
        if( tempAtom->getAtomNameFromInputFile() == " N  ") {
            secondAtom = tempAtom;
            break;
        }
    }
        
    if ( firstAtom != rg_NULL && secondAtom != rg_NULL )
	{
		disconnectAndRemoveBondBetween(firstAtom, secondAtom);
	}
}

void Molecule::removeAminoResidueNConstituentAtoms(Residue* prevResidue,
		                                           Residue* nextResidue,
								                   Residue* targetResidue)
{
	// remove chemical bonds of target residue
	disconnectAndRemoveBondsOfResidue(targetResidue);

	// remove atoms of target residue
	removeAtomsOfResidue(targetResidue);

	// remove target residue from global residue list
	m_residues.reset4Loop();
	while (m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		if(currResidue == targetResidue)
		{
			m_residues.killCurrent();		
			break;
		}
	}

	// connect previous and next residues
	if(prevResidue != rg_NULL && nextResidue != rg_NULL)	
	{
		MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(prevResidue, nextResidue, *this);
	}
	// previous residue is the end of this chain
	else if(prevResidue != rg_NULL && nextResidue == rg_NULL)
	{
		// do nothing
	}
	// next residue is the start of this chain
	// else if(nextResidue != rg_NULL && prevResidue == rg_NULL)
	else
	{
		// do nothing
	}
}

void Molecule::removeAminoResidueNConstituentAtoms(Residue* targetResidue)
{
	// get previous and next residues
	Residue* prevResidue = targetResidue->getResidueConnectedByNitrogenInAminoGroupOfAminoResidue();
	Residue* nextResidue = targetResidue->getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue();	

	// remove chemical bonds of target residue
	disconnectAndRemoveBondsOfResidue(targetResidue);

	// remove atoms of target residue
	removeAtomsOfResidue(targetResidue);

	// remove target residue from global residue list
	m_residues.reset4Loop();
	while (m_residues.setNext4Loop())
	{
		Residue* currResidue = m_residues.getpEntity();
		if(currResidue == targetResidue)
		{
			m_residues.killCurrent();		
			break;
		}
	}

	// connect previous and next residues
	if(prevResidue != rg_NULL && nextResidue != rg_NULL)	
	{
		MoleculeIOFunctions::setChemicalBondBetweenAminoAcids(prevResidue, nextResidue, *this);
	}
	// previous residue is the end of this chain
	else if(prevResidue != rg_NULL && nextResidue == rg_NULL)
	{
		// do nothing
	}
	// next residue is the start of this chain
	// else if(nextResidue != rg_NULL && prevResidue == rg_NULL)
	else
	{
		// do nothing
	}
}

void Molecule::removeAtomsOfResidue(Residue* targetResidue)
{
	rg_dList<Atom*> atomsOfTargetResidue;
	targetResidue->getAtoms(atomsOfTargetResidue);
	atomsOfTargetResidue.reset4Loop();
	while (atomsOfTargetResidue.setNext4Loop())
	{
		Atom* currAtomOfTargetResidue = atomsOfTargetResidue.getEntity();
		m_atoms.reset4Loop();
		while (m_atoms.setNext4Loop())
		{
			Atom* currAtom = m_atoms.getpEntity();
			if(currAtom->getID() == currAtomOfTargetResidue->getID())
			{
				m_atoms.killCurrent();
				atomsOfTargetResidue.killCurrent();
				break;
			}
		}
	}
}

void Molecule::disconnectAndRemoveBondsOfResidue(Residue* targetResidue)
{
	rg_dList<Atom*> atomsOfTargetResidue;
	targetResidue->getAtoms(atomsOfTargetResidue);
	atomsOfTargetResidue.reset4Loop();
	while (atomsOfTargetResidue.setNext4Loop())
	{
		Atom* currAtom = atomsOfTargetResidue.getEntity();
		disconnectAndRemoveBondsOfAtom(currAtom);
	}
}

void Molecule::disconnectAndRemoveBondsOfAtom(Atom* targetAtom)
{
	rg_dList<ChemicalBond*>* bondListOfTargetAtom 
		= targetAtom->getListChemicalBond();

	// disconnect bonds with other atoms
	bondListOfTargetAtom->reset4Loop();
	while (bondListOfTargetAtom->setNext4Loop())
	{
		ChemicalBond* currBond = bondListOfTargetAtom->getEntity();

		Atom* bondedAtom = currBond->getBondedAtom(targetAtom);
		rg_dList<ChemicalBond*>* bondListOfBondedAtom 
			= bondedAtom->getListChemicalBond();

		bondListOfBondedAtom->reset4Loop();
		while (bondListOfBondedAtom->setNext4Loop())
		{
			ChemicalBond* currBondOfBondedAtom 
				= bondListOfBondedAtom->getEntity();
			if(currBond->getID() == currBondOfBondedAtom->getID())
			{
				bondListOfBondedAtom->killCurrent();
				break;
			}
		}
	}

	// remove chemical bonds from global bond list
	bondListOfTargetAtom->reset4Loop();
	while (bondListOfTargetAtom->setNext4Loop())
	{
		ChemicalBond* currBondOfTargetAtom 
			= bondListOfTargetAtom->getEntity();
		m_chemicalBonds.reset4Loop();
		while (m_chemicalBonds.setNext4Loop())
		{
			ChemicalBond* currBond = m_chemicalBonds.getpEntity();
			if(currBondOfTargetAtom->getID() == currBond->getID())
			{
				m_chemicalBonds.killCurrent();
				bondListOfTargetAtom->killCurrent();
				break;
			}
		}
	}
}

void Molecule::disconnectAndRemoveBondBetween(Atom* firstAtom, Atom* secondAtom)
{
	rg_dList<ChemicalBond*>* bondsOfFirstAtom = firstAtom->getListChemicalBond();
	bondsOfFirstAtom->reset4Loop();
	while (bondsOfFirstAtom->setNext4Loop())
	{
		ChemicalBond* currBond = bondsOfFirstAtom->getEntity();
		if(secondAtom == currBond->getBondedAtom(firstAtom))
		{
			bondsOfFirstAtom->killCurrent();
			break;
		}
	}
	rg_dList<ChemicalBond*>* bondsOfSecondAtom = secondAtom->getListChemicalBond();
	bondsOfSecondAtom->reset4Loop();
	while (bondsOfSecondAtom->setNext4Loop())
	{
		ChemicalBond* currBond = bondsOfSecondAtom->getEntity();
		if(firstAtom == currBond->getBondedAtom(secondAtom))
		{
			bondsOfSecondAtom->killCurrent();
			break;
		}
	}

	m_chemicalBonds.reset4Loop();
	while (m_chemicalBonds.setNext4Loop())
	{
		ChemicalBond* currBond = m_chemicalBonds.getpEntity();
		Atom* atom[ 2 ];
		currBond->getAtoms(atom[ 0 ], atom[ 1 ]);
		if( (firstAtom == atom[ 0 ] && secondAtom == atom[ 1 ]) ||
			(firstAtom == atom[ 1 ] && secondAtom == atom[ 0 ])    )
		{
			m_chemicalBonds.killCurrent();
			break;
		}
	}
}

Residue* Molecule::replaceAminoResidue(Residue* targetResidue, 
					                   Residue* substitue)
{
	// get previous and next residues
	Residue* prevResidue = rg_NULL;	
	Residue* nextResidue = rg_NULL;
	prevResidue = targetResidue->getResidueConnectedByNitrogenInAminoGroupOfAminoResidue();
	nextResidue = targetResidue->getResidueConnectedByCarbonInCarboxylGroupOfAminoResidue();	

	if(prevResidue == rg_NULL && nextResidue == rg_NULL)
		return rg_NULL;

	// remove target residue
	removeAminoResidueNConstituentAtoms(prevResidue,
		                                nextResidue,
		                                targetResidue);

	// add new residue
	Residue* replacedResidue = 
		addNewResidue(prevResidue,
		              nextResidue,
					  substitue);

	return targetResidue;
}

Residue* Molecule::replaceAminoResidue(Residue* prevResidue,
		                               Residue* nextResidue,
		                               Residue* targetResidue, 
					                   Residue* substitue)
{
	if(prevResidue == rg_NULL && nextResidue == rg_NULL)
		return rg_NULL;

	// remove target residue
	removeAminoResidueNConstituentAtoms(targetResidue);

	// add new residue
	Residue* replacedResidue = 
		addNewResidue(prevResidue,
		              nextResidue,
					  substitue);

	return targetResidue;
}

//void Molecule::replaceAminoResidue(Residue* targetResidue, const ResidueCode& newResidueCode, const rg_BOOL& bBetaCarbonComputed/* = rg_FALSE*/)
//copy N, Ca, C, O, and Cb(beta carbon) from OLD residue
void Molecule::replaceAminoResidue(Residue* targetResidue, const ResidueCode& newResidueCode)
{
	targetResidue->setResidueCode(newResidueCode);
	targetResidue->setResidueName(string(RESIDUE_FEATURES[newResidueCode].threeCodeName));
	Sphere betaCarbonBall = targetResidue->getBetaCarbonInSideChainOfAminoResidue()->getAtomBall();

	// remove atoms and chemical bonds on sidechain
	removeAtomsNTheirChemicalBondsOfSidechain(targetResidue);

	// make atoms and chemical bonds on sidechain
	if(newResidueCode == GLY_AMINO_RESIDUE)
		return;

	rg_Point3D defaultAtomCenter(9999.999, 9999.999, 9999.999);
	AtomSymbolMap mapOfAtomSymbols;
	mapOfAtomSymbols.clear();
	for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
		string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
		mapOfAtomSymbols.insert( AtomSymbolMap::value_type(atomSymbol, i_atomCode) );
	}

	recoverMissingAtomsInSidechainOfResidueWithGivenCoordinates(targetResidue,
		                                                        defaultAtomCenter,
		                                                        mapOfAtomSymbols);
	removeHydrongenInSidechain(targetResidue);
    // set beta carbon
    targetResidue->getBetaCarbonInSideChainOfAminoResidue()->getpAtomBall()->setSphere(betaCarbonBall.getCenter(), betaCarbonBall.getRadius());

	//if(! bBetaCarbonComputed)
	//    targetResidue->getBetaCarbonInSideChainOfAminoResidue()->getpAtomBall()->setSphere(betaCarbonBall.getCenter(), betaCarbonBall.getRadius());
 //   else
 //   {
 //       // construct candidates for coordinate of beta carbon and
 //       // choose the best one
 //   }
}

void Molecule::removeAtomsNTheirChemicalBondsOfSidechain(Residue* targetResidue)
{
	rg_dList<Atom*> atomsOfSidechain;
	targetResidue->getAtomsOnSideChain(&atomsOfSidechain);
	atomsOfSidechain.reset4Loop();
	while (atomsOfSidechain.setNext4Loop())
	{
		Atom* currAtom = atomsOfSidechain.getEntity();
		deleteAtom(rg_FALSE, currAtom);
	}
}

void Molecule::removeHydrongenInSidechain(Residue* targetResidue)
{
#ifdef _DEBUG
		if(targetResidue->getSequenceNumber() == 77)
			int here = 1;
#endif
	rg_dList<Atom*> atomsOfSidechain;
	targetResidue->getAtomsOnSidechain(atomsOfSidechain);
	atomsOfSidechain.reset4Loop();
	while(atomsOfSidechain.setNext4Loop())
	{
		Atom* currAtom = atomsOfSidechain.getEntity();
		if(currAtom->getAtomCode() == H_ATOM)
			deleteAtom(rg_FALSE, currAtom);
	}
}
//
//void Molecule::replaceAminoResidue(Residue* targetResidue, const ResidueCode& newResidueCode, const rg_BOOL& bBetaCarbonComputed/* = rg_FALSE*/)
//{
//	targetResidue->setResidueCode(newResidueCode);
//	targetResidue->setResidueName(string(RESIDUE_FEATURES[newResidueCode].threeCodeName));
//	Sphere betaCarbonBall = targetResidue->getBetaCarbonInSideChainOfAminoResidue()->getAtomBall();
//
//	// remove atoms and chemical bonds on sidechain
//	removeChemicalBondsOnSidechaOfResidue(targetResidue);
//	removeAtomsOnSidechainOfResidue(targetResidue);
//
//	// make atoms and chemical bonds on sidechain
//	if(newResidueCode == GLY_AMINO_RESIDUE)
//		return;
//
//	rg_Point3D defaultAtomCenter(9999.999, 9999.999, 9999.999);
//	AtomSymbolMap mapOfAtomSymbols;
//	mapOfAtomSymbols.clear();
//	for ( rg_INT i_atomCode = 0; i_atomCode<NUM_OF_ATOMS; i_atomCode++ ) {
//		string atomSymbol( ATOM_FEATURES[i_atomCode].symbol );
//		mapOfAtomSymbols.insert( AtomSymbolMap::value_type(atomSymbol, i_atomCode) );
//	}
//
//	recoverMissingAtomsInSidechainOfResidueWithGivenCoordinates(targetResidue,
//		defaultAtomCenter,
//		mapOfAtomSymbols);
//	// set beta carbon
//	targetResidue->getBetaCarbonInSideChainOfAminoResidue()->getpAtomBall()->setSphere(betaCarbonBall.getCenter(), betaCarbonBall.getRadius());
//}

void Molecule::removeChemicalBondsOnSidechaOfResidue(Residue* targetResidue)
{
	rg_dList<Atom*> atomsOfSidechain;
	targetResidue->getAtomsOnSideChain(&atomsOfSidechain);
	atomsOfSidechain.reset4Loop();
	while (atomsOfSidechain.setNext4Loop())
	{
		Atom* currAtom = atomsOfSidechain.getEntity();
		disconnectAndRemoveBondsOfAtom(currAtom);
	}
}

void Molecule::removeAtomsOnSidechainOfResidue(Residue* targetResidue)
{
	rg_dList<Atom*>* atomsOfTargetResidue = targetResidue->getAtoms();
	atomsOfTargetResidue->reset4Loop();
	while (atomsOfTargetResidue->setNext4Loop())
	{
		Atom* currAtomOfTargetResidue = atomsOfTargetResidue->getEntity();

		if( currAtomOfTargetResidue->getpChemicalProperties()->isOnBackBone() == rg_FALSE &&
			currAtomOfTargetResidue->getAtomNameFromInputFile() != " O  " ) 
		{
			m_atoms.reset4Loop();
			while (m_atoms.setNext4Loop())
			{
				Atom* currAtom = m_atoms.getpEntity();
				if(currAtom->getID() == currAtomOfTargetResidue->getID())
				{
					m_atoms.killCurrent();
					atomsOfTargetResidue->killCurrent();
					break;
				}
			}
		}
	}
}


Molecule& Molecule::operator=( const Molecule& aMolecule )
{
    if (this != &aMolecule) {
        duplicate(aMolecule);
    }
        
    return *this;

    //if( this == &aMolecule )
    //    return *this;

    //m_moleculeFileName         = aMolecule.m_moleculeFileName;
    //m_atomRadiusType           = aMolecule.m_atomRadiusType;
    //m_headerRecords            = aMolecule.m_headerRecords;

    //m_atoms                    = aMolecule.m_atoms;
    //m_residues                 = aMolecule.m_residues;
    //m_chains                   = aMolecule.m_chains;
    //
    //m_chemicalBonds            = aMolecule.m_chemicalBonds;
    //m_pharmaFeatures           = aMolecule.m_pharmaFeatures;
    //
    //m_centerOfMass             = aMolecule.m_centerOfMass;
    //m_minEnclosingSphere       = aMolecule.m_minEnclosingSphere;
    //m_modelSerialFromInputFile = aMolecule.m_modelSerialFromInputFile;

    //// JKKIM ADDED //
    //m_moleculeName             = aMolecule.m_moleculeName;

    //// by Y.Cho at 2011-10-12
    //m_timeStamp                = aMolecule.m_timeStamp;
    //m_fileSize                 = aMolecule.m_fileSize;

    //return *this;
}



void Molecule::syncPtrsWithSourceMolecule( Molecule* sourceMolecule )
{
    AtomMap          mapOfAtom;
    ResidueMap       mapOfResidue;
    ChainMap         mapOfChain;
    ChemBondMap      mapOfChemBond;
    PharmaFeatureMap mapOfPharmaFeature;
    
    initializeMapsForSync( mapOfAtom, mapOfResidue, mapOfChain, mapOfChemBond, mapOfPharmaFeature );
    
    syncPtrsInAtom( sourceMolecule->getAtoms(), &mapOfResidue, &mapOfChemBond, &mapOfPharmaFeature );
    syncPtrsInResidue( sourceMolecule->getResidues(), &mapOfAtom, &mapOfChain );
    syncPtrsInChain( sourceMolecule->getChains(), &mapOfResidue );
    syncPtrsInChemBond( sourceMolecule->getChemicalBonds(), &mapOfAtom );
    syncPtrsInPharmaFeatures( sourceMolecule->getListPharmaFeatures(), &mapOfAtom );
}



void Molecule::initializeMapsForSync( AtomMap& mapOfAtom, ResidueMap& mapOfResidue, ChainMap& mapOfChain, ChemBondMap& mapOfChemBond, PharmaFeatureMap& mapOfPharmaFeature)
{
    // INITIALIZE MAP OF ATOMS
    m_atoms.reset4Loop();
    while ( m_atoms.setNext4Loop() ) {
        Atom* currAtom = m_atoms.getpEntity();        
        mapOfAtom.insert( AtomMap::value_type( currAtom->getID(), currAtom ) );
    }

    // INITIALIZE MAP OF RESIDUES
    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        Residue* currResidue = m_residues.getpEntity();
        mapOfResidue.insert( ResidueMap::value_type( currResidue->getID(), currResidue ) );
    }

    // INITIALIZE MAP OF CHAINS
    m_chains.reset4Loop();
    while ( m_chains.setNext4Loop() ) {
        Chain* currChain = m_chains.getpEntity();
        mapOfChain.insert( ChainMap::value_type( currChain->getID(), currChain ) );
    }

    // INITIALIZE MAP OF CHEMICALBONDS
    m_chemicalBonds.reset4Loop();
    while ( m_chemicalBonds.setNext4Loop() ) {
        ChemicalBond* currBond = m_chemicalBonds.getpEntity();
        mapOfChemBond.insert( ChemBondMap::value_type( currBond->getID(), currBond ) );
    }

    // INITIALIZE MAP OF PHARMAFEATURES
    list<PharmaFeature>::iterator pharmafeature_i = m_pharmaFeatures.begin();
    
    while ( pharmafeature_i != m_pharmaFeatures.end() )    {
        PharmaFeature* currPharmaFeature = &(*pharmafeature_i);
        mapOfPharmaFeature.insert( PharmaFeatureMap::value_type( currPharmaFeature->getID(), currPharmaFeature ) );
        pharmafeature_i++;
    }
}



void Molecule::syncPtrsInAtom( rg_dList<Atom>* sourceAtoms, ResidueMap* mapOfResidue, ChemBondMap* mapOfChemBond, PharmaFeatureMap* mapOfPharmaFeature )
{
    m_atoms.reset4Loop();
    sourceAtoms->reset4Loop();
    while ( sourceAtoms->setNext4Loop() ) {
        m_atoms.setNext4Loop();

        Atom* currSourceAtom = sourceAtoms->getpEntity();
        Atom* currTargetAtom = m_atoms.getpEntity();

        // SYNC RESIDUE
        ResidueMap::iterator residueMap_i = mapOfResidue->find( currSourceAtom->getResidue()->getID() );
        Residue* residueOfCurrTargetAtom = (*residueMap_i).second;
        currTargetAtom->setResidue( residueOfCurrTargetAtom );


        // SYNC CHEMICAL BONDS
        rg_dList<ChemicalBond*>* chemBondsOfSourceAtom = currSourceAtom->getListChemicalBond();
        rg_dList<ChemicalBond*>* chemBondsOfTargetAtom = currTargetAtom->getListChemicalBond();
        chemBondsOfTargetAtom->removeAll();

        chemBondsOfSourceAtom->reset4Loop();
        while ( chemBondsOfSourceAtom->setNext4Loop() ) {
            ChemBondMap::iterator bondMap_i = mapOfChemBond->find( chemBondsOfSourceAtom->getEntity()->getID() );
            ChemicalBond* chemBondForTargetAtom = (*bondMap_i).second;
            chemBondsOfTargetAtom->addTail( chemBondForTargetAtom );
        }

        // SYNC PHARMA FEATURES
        rg_dList<PharmaFeature*>* pharmaOfSourceAtom = currSourceAtom->getpChemicalProperties()->getListOfPharmaFeatures();
        rg_dList<PharmaFeature*>* pharmaOfTargetAtom = currTargetAtom->getpChemicalProperties()->getListOfPharmaFeatures();
        pharmaOfTargetAtom->removeAll();

        pharmaOfSourceAtom->reset4Loop();
        while ( pharmaOfSourceAtom->setNext4Loop() ) {
            PharmaFeatureMap::iterator pharmaMap_i = mapOfPharmaFeature->find( pharmaOfSourceAtom->getEntity()->getID() );
            PharmaFeature* pharmaForTargetAtom = (*pharmaMap_i).second;
            pharmaOfTargetAtom->addTail( pharmaForTargetAtom );
        }
    }
}



void Molecule::syncPtrsInResidue( rg_dList<Residue>* sourceResidues, AtomMap* mapOfAtom, ChainMap* mapOfChain )
{
    m_residues.reset4Loop();
    sourceResidues->reset4Loop();
    while ( sourceResidues->setNext4Loop() ) {
        m_residues.setNext4Loop();

        Residue* currSourceResidue = sourceResidues->getpEntity();
        Residue* currTargetResidue = m_residues.getpEntity();

        // SYNC ATOMS
        rg_dList<Atom*>* atomsOfSourceResidue = currSourceResidue->getAtoms();
        rg_dList<Atom*>* atomsOfTargetResidue = currTargetResidue->getAtoms();
        atomsOfTargetResidue->removeAll();

        atomsOfSourceResidue->reset4Loop();
        while ( atomsOfSourceResidue->setNext4Loop() ) {
            AtomMap::iterator atomMap_i = mapOfAtom->find( atomsOfSourceResidue->getEntity()->getID() );
            Atom* atomForCurrTargetResidue = (*atomMap_i).second;
            atomsOfTargetResidue->addTail( atomForCurrTargetResidue );
        }

        // SYNC CHAIN
        ChainMap::iterator chainMap_i = mapOfChain->find( currSourceResidue->getChain()->getID() );
        Chain* chainOfCurrTargetResidue = (*chainMap_i).second;
        currTargetResidue->setChain( chainOfCurrTargetResidue );
    }
}



void Molecule::syncPtrsInChain( rg_dList<Chain>* sourceChains, ResidueMap* mapOfResidue )
{
    m_chains.reset4Loop();
    sourceChains->reset4Loop();
    while ( sourceChains->setNext4Loop() ) {
        m_chains.setNext4Loop();
        
        Chain* currSourceChain = sourceChains->getpEntity();
        Chain* currTargetChain = m_chains.getpEntity();
        
        // SYNC RESIDUES
        rg_dList<Residue*>* residuesOfSourceChain = currSourceChain->getResidues();
        rg_dList<Residue*>* residuesOfTargetChain = currTargetChain->getResidues();
        residuesOfTargetChain->removeAll();
        
        residuesOfSourceChain->reset4Loop();
        while ( residuesOfSourceChain->setNext4Loop() ) {
            ResidueMap::iterator residueMap_i = mapOfResidue->find( residuesOfSourceChain->getEntity()->getID() );
            Residue* residueForCurrTargetChain = (*residueMap_i).second;
            residuesOfTargetChain->addTail( residueForCurrTargetChain );
        }

        // SYNC SECONDARY STRUCTURES
        syncPtrsInSecondaryStructure( currSourceChain, mapOfResidue, currTargetChain );
        
        // SYNC MOLECULE
        currTargetChain->setMolecule( this );
    }
}



void Molecule::syncPtrsInSecondaryStructure( Chain* sourceChain, ResidueMap* mapOfResidue, Chain* targetChain )
{
    SecondaryStructure* sourceSecStruct = sourceChain->getSecondaryStructure();
    SecondaryStructure* targetSecStruct = targetChain->getSecondaryStructure();


    // SYNC HELIX
    rg_dList<Helix>* sourceHelices = sourceSecStruct->getHelices();
    rg_dList<Helix>* targetHelices = targetSecStruct->getHelices();
    
    sourceHelices->reset4Loop();
    targetHelices->reset4Loop();
    while ( sourceHelices->setNext4Loop() ) {
        targetHelices->setNext4Loop();

        Helix* currSourceHelix = sourceHelices->getpEntity();
        Helix* currTargetHelix = targetHelices->getpEntity();
        
        rg_dList<Residue*>* residuesInSourceHelix = currSourceHelix->getResidues();
        rg_dList<Residue*>* residuesInTargetHelix = currTargetHelix->getResidues();
        residuesInTargetHelix->removeAll();

        residuesInSourceHelix->reset4Loop();
        while ( residuesInSourceHelix->setNext4Loop() ) {
            Residue* currSourceResidue = residuesInSourceHelix->getEntity();
            
            ResidueMap::iterator residueMap_i = mapOfResidue->find( currSourceResidue->getID() );
            Residue* currTargetResidue = (*residueMap_i).second;
            residuesInTargetHelix->addTail( currTargetResidue );
        }
    }

    
    // SYNC SHEET   
    rg_dList<Sheet>* sourceSheets = sourceSecStruct->getSheets();
    rg_dList<Sheet>* targetSheets = targetSecStruct->getSheets();
    
    sourceSheets->reset4Loop();
    targetSheets->reset4Loop();
    while ( sourceSheets->setNext4Loop() ) {
        targetSheets->setNext4Loop();

        Sheet* currSourceSheet = sourceSheets->getpEntity();
        Sheet* currTargetSheet = targetSheets->getpEntity();
    
        Strand* strandsInSourceSheet = currSourceSheet->getStrands();
        Strand* strandsInTargetSheet = currTargetSheet->getStrands();

        rg_INT numOfStrands = currSourceSheet->getNumOfStrands();
        for ( rg_INT i_strand=0; i_strand<numOfStrands; i_strand++ ) {
            rg_dList<Residue*>* residuesInSourceStrand = strandsInSourceSheet[i_strand].getResidues();
            rg_dList<Residue*>* residuesInTargetStrand = strandsInTargetSheet[i_strand].getResidues();
            residuesInTargetStrand->removeAll();

            residuesInSourceStrand->reset4Loop();
            while ( residuesInSourceStrand->setNext4Loop() ) {
                Residue* currSourceResidue = residuesInSourceStrand->getEntity();
            
                ResidueMap::iterator residueMap_i = mapOfResidue->find( currSourceResidue->getID() );
                Residue* currTargetResidue = (*residueMap_i).second;
                residuesInTargetStrand->addTail( currTargetResidue );
            }
        }
    }


    // SYNC TURN
    rg_dList<Turn>* sourceTurns = sourceSecStruct->getTurns();
    rg_dList<Turn>* targetTurns = targetSecStruct->getTurns();
    
    sourceTurns->reset4Loop();
    targetTurns->reset4Loop();
    while ( sourceTurns->setNext4Loop() ) {
        targetTurns->setNext4Loop();

        Turn* currSourceTurn = sourceTurns->getpEntity();
        Turn* currTargetTurn = targetTurns->getpEntity();
        
        rg_dList<Residue*>* residuesInSourceTurn = currSourceTurn->getResidues();
        rg_dList<Residue*>* residuesInTargetTurn = currTargetTurn->getResidues();
        residuesInTargetTurn->removeAll();

        residuesInSourceTurn->reset4Loop();
        while ( residuesInSourceTurn->setNext4Loop() ) {
            Residue* currSourceResidue = residuesInSourceTurn->getEntity();
            
            ResidueMap::iterator residueMap_i = mapOfResidue->find( currSourceResidue->getID() );
            Residue* currTargetResidue = (*residueMap_i).second;
            residuesInTargetTurn->addTail( currTargetResidue );
        }
    }
}



void Molecule::syncPtrsInChemBond( rg_dList<ChemicalBond>* sourceBonds, AtomMap* mapOfAtom )
{
    m_chemicalBonds.reset4Loop();
    sourceBonds->reset4Loop();
    while ( sourceBonds->setNext4Loop() ) {
        m_chemicalBonds.setNext4Loop();

        ChemicalBond* currSourceBond = sourceBonds->getpEntity();
        ChemicalBond* currTargetBond = m_chemicalBonds.getpEntity();

        // SYNC ATOMS
        Atom* firstAtomFromSourceBond = currSourceBond->getFirstAtom();
        Atom* secondAtomFromSourceBond = currSourceBond->getSecondAtom();

        AtomMap::iterator atomMap_i = mapOfAtom->find( firstAtomFromSourceBond->getID() );
        Atom* atomForCurrTargetBond = (*atomMap_i).second;
        currTargetBond->setFirstAtom( atomForCurrTargetBond );

        atomMap_i = mapOfAtom->find( secondAtomFromSourceBond->getID() );
        atomForCurrTargetBond = (*atomMap_i).second;
        currTargetBond->setSecondAtom( atomForCurrTargetBond );
    }
}



void Molecule::syncPtrsInPharmaFeatures( list<PharmaFeature>* sourcePharmaFeatures, AtomMap* mapOfAtom )
{
    list<PharmaFeature>::iterator sourcePharmafeature_i = sourcePharmaFeatures->begin();
    list<PharmaFeature>::iterator targetPharmafeature_i = m_pharmaFeatures.begin();
    
    while ( sourcePharmafeature_i != sourcePharmaFeatures->end() )    {
        PharmaFeature* currSourcePharmaFeature = &(*sourcePharmafeature_i);
        PharmaFeature* currTargetPharmaFeature = &(*targetPharmafeature_i);
        
        rg_dList<Atom*>* atomsOfSourcePharma = currSourcePharmaFeature->getAtoms();
        rg_dList<Atom*>* atomsOfTargetPharma = currTargetPharmaFeature->getAtoms();
        atomsOfTargetPharma->removeAll();

        atomsOfSourcePharma->reset4Loop();
        while ( atomsOfSourcePharma->setNext4Loop() ) {
            AtomMap::iterator atomMap_i = mapOfAtom->find( atomsOfSourcePharma->getEntity()->getID() );
            Atom* atomForCurrTargetResidue = (*atomMap_i).second;
            atomsOfTargetPharma->addTail( atomForCurrTargetResidue );
        }

        sourcePharmafeature_i++;
        targetPharmafeature_i++;
    }

}



void Molecule::computeAndSetCenterOfMass()
{
    rg_Point3D sumOfAllAtomCenters;
    m_atoms.reset4Loop();
    while( m_atoms.setNext4Loop() ) {
        sumOfAllAtomCenters += m_atoms.getpEntity()->getAtomBall().getCenter();
    }
    m_centerOfMass = sumOfAllAtomCenters/m_atoms.getSize();
}



void Molecule::computeAndSetMinEnclosingSphere()
{
    EnclosingSphereOfSpheres enclosingSphereOfMolecule;

    list<Sphere> listOfAtomBalls;
    m_atoms.reset4Loop();
    while( m_atoms.setNext4Loop() ) {
        Sphere atomBall = m_atoms.getpEntity()->getAtomBall();
        listOfAtomBalls.push_back( atomBall );
    }
        
    enclosingSphereOfMolecule.setSpheres( listOfAtomBalls );
    enclosingSphereOfMolecule.computeEnclosingSphere();
    
    m_minEnclosingSphere = enclosingSphereOfMolecule;    
}



void Molecule::evaluateHydrogenDonorAndAcceptorAtoms()
{
    evaluateHydrogenDonorAtoms();
    evaluateHydrogenAcceptorAtoms();
}



void Molecule::evaluateHydrogenDonorAtoms()
{
// Donnor
// O, N, S

    m_atoms.reset4Loop();
    while( m_atoms.setNext4Loop() ) {
        Atom* currAtom = m_atoms.getpEntity();
        if( currAtom->getAtomCode() != N_ATOM || currAtom->getAtomCode() != O_ATOM || currAtom->getAtomCode() != S_ATOM ) {
            continue;
        }
        rg_INT numOfConnectedHydrogen = 0;
        rg_dList<ChemicalBond*>* listOfChemicalBond = currAtom->getListChemicalBond();
        listOfChemicalBond->reset4Loop();
        while( listOfChemicalBond->setNext4Loop() ) {
            ChemicalBond* currChemBond = listOfChemicalBond->getEntity();
            if( currChemBond->getFirstAtom()->getAtomCode()  == H_ATOM || 
                currChemBond->getSecondAtom()->getAtomCode() == H_ATOM    ) {
                ++numOfConnectedHydrogen;
            }
        }
        currAtom->getpChemicalProperties()->setNumOfConnectedHydrogen( numOfConnectedHydrogen );
    }
}



void Molecule::evaluateHydrogenAcceptorAtoms()
{
////////////////////////////// CONDITIONS TO BE A HYDROGEN ACCEPTOR ATOM //////////////////////////////
//                                                                                                   //
//  1. If a Residue is a AminoAcid
//          1) O_ATOM can accept 2 hydrogens
//          2) N_ATOM in HIS_AMINO_RESIDUE
//              i) No H_ATOM connected
//                  a) with DELTA_REMOTE and FIRST_BRANCH can accept 1 hydrogen
//                  b) with EPSILON_REMOTE and FIRST_BRANCH can accept 1 hydrogen
//          4) S_ATOM can accept 2 hydrogens
//
//  2. If a Residue is a NeucleicAcid
//          1) O_ATOM can accept 2 hydrogens
//          2) N_ATOM can accept 1 hydrogen
//          3) S_ATOM can accept 2 hydrogens
//
//  3. If a Residue is a LIGAND
//          1) O_ATOM can accept 2 hydrogens
//          2) N_ATOM with following conditions :
//            i) No H_ATOM connected
//                  a) one atom connected :
//                      * with SINGLE_BOND - can accept 1 hydrogen
//                      * with DOUBLE_BOND - can accept 2 hydrogens
//                      * with TRIPLE_BOND - can accept 1 hydrogen
//               b) two atoms connected :
//                      * one with SINGLE_BOND and another with DOUBL_BOND - can accept 1 hydrogen
//                      * with AROMATIC_BOND - can accept 1 hydrogen
//             ii) One H_ATOM connected : 
//                  a) all connected atoms are with SINGLE_BOND - can accept 1 hydrogen
//             iii) Two H_ATOM connected : 
//                  a) all connected atoms are with SINGLE_BOND - can accept 1 hydrogen
//          3) S_ATOM can accept 2 hydrogens
//          4) F_ATOM can accept 3 hydrogens
//          5) CL_ATOM can accept 3 hydrogens
//                                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////////

    rg_INT numOfAcceptableHydrogen = 0;

    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        
        Residue* currResidue = m_residues.getpEntity();
        rg_dList<Atom*>* atomsInCurrResidue = currResidue->getAtoms();
        atomsInCurrResidue->reset4Loop();
        
        if( currResidue->isAminoResidue() ) { 
            
            while ( atomsInCurrResidue->setNext4Loop() ) {
                numOfAcceptableHydrogen = 0;
                Atom*    currAtom = atomsInCurrResidue->getEntity();
                AtomCode atomCodeOfCurrAtom = currAtom->getAtomCode();
                ChemicalPropertiesOfAtom* chemPropOfCurrAtom = currAtom->getpChemicalProperties();

                if( atomCodeOfCurrAtom == O_ATOM || atomCodeOfCurrAtom == S_ATOM ) {
                    numOfAcceptableHydrogen = 2;
                }
                else if( atomCodeOfCurrAtom == N_ATOM && currResidue->getResidueCode() == HIS_AMINO_RESIDUE &&
                         chemPropOfCurrAtom->getBrangeDesignator() == FIRST_BRANCH ) {
                    if( chemPropOfCurrAtom->getNumOfConnectedHydrogen() == 0 ) {
                        if( chemPropOfCurrAtom->getRemoteIndicator()  == DELTA_REMOTE ||
                            chemPropOfCurrAtom->getRemoteIndicator()  == EPSILON_REMOTE )
                            numOfAcceptableHydrogen = 1;
                    }
                }
                chemPropOfCurrAtom->setNumOfAcceptableHydrogen( numOfAcceptableHydrogen );
            }
        }

        else if( currResidue->isNeucleicAcidResidue() ) {

            while ( atomsInCurrResidue->setNext4Loop() ) {
                numOfAcceptableHydrogen = 0;
                Atom*    currAtom = atomsInCurrResidue->getEntity();
                AtomCode atomCodeOfCurrAtom = currAtom->getAtomCode();

                if( atomCodeOfCurrAtom == O_ATOM || atomCodeOfCurrAtom == S_ATOM ) {
                    numOfAcceptableHydrogen = 2;
                }
                else if( atomCodeOfCurrAtom == N_ATOM ) {
                    numOfAcceptableHydrogen = 1;
                }
                currAtom->getpChemicalProperties()->setNumOfAcceptableHydrogen( numOfAcceptableHydrogen );                
            }
        }

        else if( !currResidue->isStandardResidue() ) {

            while ( atomsInCurrResidue->setNext4Loop() ) {
                numOfAcceptableHydrogen = 0;
                Atom*    currAtom = atomsInCurrResidue->getEntity();
                AtomCode atomCodeOfCurrAtom = currAtom->getAtomCode();
                ChemicalPropertiesOfAtom* chemPropOfCurrAtom = currAtom->getpChemicalProperties();

                if( atomCodeOfCurrAtom == O_ATOM || atomCodeOfCurrAtom == S_ATOM ) {
                    numOfAcceptableHydrogen = 2;
                }
                else if( atomCodeOfCurrAtom == N_ATOM ) {
                    rg_INT numOfConnectedHydrogen = chemPropOfCurrAtom->getNumOfConnectedHydrogen();
                    if( numOfConnectedHydrogen == 0 ) {
                        rg_dList<ChemicalBond*>* listOfChemicalBond = currAtom->getListChemicalBond();
                        rg_INT numOfChemicalBond = listOfChemicalBond->getSize();
                        if( numOfChemicalBond == 1 ) {
                            rg_INT typeOfBond = listOfChemicalBond->getFirstEntity()->getTypeOfBond();
                            if( typeOfBond == SINGLE_BOND )
                                numOfAcceptableHydrogen = 1;
                            else if( typeOfBond == DOUBLE_BOND )
                                numOfAcceptableHydrogen = 2;
                            else if( typeOfBond == TRIPLE_BOND )
                                numOfAcceptableHydrogen = 1;
                        }
                        else if( numOfChemicalBond == 2 ) {
                            rg_INT typeOfBondForFirstBond  = listOfChemicalBond->getFirstEntity()->getTypeOfBond();
                            rg_INT typeOfBondForSecondBond = listOfChemicalBond->getLastEntity()->getTypeOfBond();
                            if( (typeOfBondForFirstBond == SINGLE_BOND && typeOfBondForSecondBond == SINGLE_BOND) ||
                                (typeOfBondForFirstBond == SINGLE_BOND && typeOfBondForSecondBond == DOUBLE_BOND) ||
                                (typeOfBondForFirstBond == DOUBLE_BOND && typeOfBondForSecondBond == SINGLE_BOND) ||
                                (typeOfBondForFirstBond == AROMATIC_BOND && typeOfBondForSecondBond == AROMATIC_BOND) )
                                numOfAcceptableHydrogen = 1;
                        }
                        else if( numOfChemicalBond >= 3 ) {
                            rg_FLAG isChemBondSingleType = rg_TRUE;
                            listOfChemicalBond->reset4Loop();
                            while( listOfChemicalBond->setNext4Loop() ) {
                                ChemicalBond* currChemicalBond = listOfChemicalBond->getEntity();
                                if( currChemicalBond->getTypeOfBond() != SINGLE_BOND ) {
                                    isChemBondSingleType = rg_FALSE;
                                    break;
                                }
                            }
                            if( isChemBondSingleType == rg_TRUE )
                                numOfAcceptableHydrogen = 1;
                        }
                    }
                    else if( numOfConnectedHydrogen == 1 || numOfConnectedHydrogen == 2 ) {
                        rg_FLAG isChemBondSingleType = rg_TRUE;
                        rg_dList<ChemicalBond*>* listOfChemicalBond = currAtom->getListChemicalBond();
                        listOfChemicalBond->reset4Loop();
                        while( listOfChemicalBond->setNext4Loop() ) {
                            ChemicalBond* currChemicalBond = listOfChemicalBond->getEntity();
                            if( currChemicalBond->getTypeOfBond() != SINGLE_BOND ) {
                                isChemBondSingleType = rg_FALSE;
                                break;
                            }
                        }
                        if( isChemBondSingleType == rg_TRUE )
                            numOfAcceptableHydrogen = 1;
                    }
                }
                else if( atomCodeOfCurrAtom == F_ATOM || atomCodeOfCurrAtom == CL_ATOM ) {
                    numOfAcceptableHydrogen = 3;
                }
                chemPropOfCurrAtom->setNumOfAcceptableHydrogen( numOfAcceptableHydrogen );
            }
        }
    }
}



void Molecule::defineAndSetPharmaFeatures()
{
   if( m_chemicalBonds.getSize() == 0 )
       return;

   definePIPharmaFeatures();
   defineNIPharmaFeatures();
   defineHBAPharmaFeatures();
   defineHBDPharmaFeatures();
   defineHYPharmaFeatures();

   //refineListOfPharmaFeatures();
}



void Molecule::definePIPharmaFeatures()
{
   Atom*  currAtom        = rg_NULL;
   rg_INT pharmaFeatureID = m_pharmaFeatures.size();
   

   m_atoms.reset4Loop();
   while( m_atoms.setNext4Loop() ) {
		
       currAtom = m_atoms.getpEntity();
       
       PharmaFeature  aPiPositiveChargeFeature;
       if( isPiPositiveCharge( currAtom, aPiPositiveChargeFeature ) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aPiPositiveChargeFeature);
       }
		
       PharmaFeature  aPiAmineFeature;
       if( isPiAmine( currAtom, aPiAmineFeature ) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aPiAmineFeature );
       }
		
		PharmaFeature  aPiAmidineFeature;
       if( isPiAmidine( currAtom, aPiAmidineFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aPiAmidineFeature);
       }
		
		PharmaFeature  aPiGuanidineFeature;
       if( isPiGuanidine( currAtom, aPiGuanidineFeature ) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aPiGuanidineFeature );
       }
   }
}



rg_FLAG Molecule::isPiPositiveCharge( Atom* currAtom, PharmaFeature& piPositiveChargeFeature )
{
   rg_FLAG isPosChargeChemType = rg_FALSE;
   
   rg_REAL chargeOfCurrAtom     = currAtom->getChemicalProperties().getCharge();
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;

   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

   if( chargeOfCurrAtom >= 1.0 ) {
       if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
           isPosChargeChemType = rg_TRUE;
           
           piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
           return isPosChargeChemType;
       }
   }

   if( currAtom->getAtomCode() == N_ATOM ) {
       if( numOfBondForCurrAtom == 4 ) {
           if( numOfSingleBond == 4 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}
           }
       }
       
       else if( numOfBondForCurrAtom == 3 ) {
           if( numOfSingleBond == 2 && numOfDoubleBond == 1 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}
           }
           else if( numOfSingleBond == 1 && numOfAromaticBond == 2 ) {
					if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
						isPosChargeChemType = rg_TRUE;
                       piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
						return isPosChargeChemType;
					}
           }
       }

       else if( numOfBondForCurrAtom == 2 ) {
           if( numOfDoubleBond == 2 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}
           }
           else if( numOfSingleBond == 1 && numOfTripleBond == 1 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}
           }
       }
   }
   
   if( currAtom->getAtomCode() == O_ATOM ) {
       if( numOfBondForCurrAtom == 3 ) {
           if( numOfSingleBond == 1 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}           
           }
       }
       
       else if( numOfBondForCurrAtom == 2 ) {
           if( numOfSingleBond == 1 && numOfDoubleBond == 1 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}           
           }
           else if( numOfAromaticBond == 2 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}
           }
       }

       else if( numOfBondForCurrAtom == 1 ) {    
           if( numOfTripleBond == 1 ) {
				if( isAllBondedAtomsNonNegative( currAtom ) == rg_TRUE ) {
					isPosChargeChemType = rg_TRUE;
                   piPositiveChargeFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, POS_CHARGE_CHEM_TYPE, currAtom );
					return isPosChargeChemType;
				}
           }
       }
   }

   return isPosChargeChemType;
}



rg_FLAG Molecule::isPiAmine( Atom* currAtom, PharmaFeature& piAmineFeature  )
{
   rg_FLAG isAmineChemType = rg_FALSE;

   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;

   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;

   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );

   Atom** arrCarbon     = new Atom* [numOfCarbon];
   Atom** arrHydrogen   = new Atom* [numOfHydrogen];
   Atom** arrNitrogen   = new Atom* [numOfNitrogen];
   Atom** arrOxygen     = new Atom* [numOfOxygen];
   Atom** arrPhosphorus = new Atom* [numOfPhosphorus];
   Atom** arrSulfur     = new Atom* [numOfSulfur];

   getEachpAtomFromBondListInAtom( currAtom, arrCarbon, arrHydrogen, arrNitrogen, arrOxygen, arrPhosphorus, arrSulfur );


   if( currAtom->getAtomCode() == N_ATOM ) {

       if( numOfBondForCurrAtom == 3 && numOfSingleBond == 3 ) {

           if( ( numOfCarbon == 1 && numOfHydrogen == 2 ) || ( numOfCarbon == 2 && numOfHydrogen == 1 ) || ( numOfCarbon == 3 ) ) {
               
               rg_INT  numOfBondForBondedCarbon;
               rg_FLAG isAllSingleBond;
               rg_FLAG isAllBondedCarbonSatisfied = rg_TRUE;

               for( rg_INT i_bondedCarbon=0; i_bondedCarbon<numOfCarbon; i_bondedCarbon++ ) {
                   numOfBondForBondedCarbon = arrCarbon[i_bondedCarbon]->getListChemicalBond()->getSize();
                   isAllSingleBond          = isAllBondsSingle( arrCarbon[i_bondedCarbon] );
                   if( numOfBondForBondedCarbon == 4 && isAllSingleBond == rg_TRUE ) {
                       isAllBondedCarbonSatisfied = rg_TRUE;
                   }
                   else {
                       isAllBondedCarbonSatisfied = rg_FALSE;
                       break;
                   }
               }

               if( isAllBondedCarbonSatisfied == rg_TRUE ) {
                   isAmineChemType = rg_TRUE;
                   piAmineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, AMINE_CHEM_TYPE, currAtom );
               }
           }
       }
   }

   delete [] arrCarbon;
   delete [] arrHydrogen;
   delete [] arrNitrogen; 
   delete [] arrOxygen;
   delete [] arrPhosphorus;
   delete [] arrSulfur;

   return isAmineChemType;    
}



rg_FLAG Molecule::isPiAmidine( Atom* currAtom, PharmaFeature& piAmidineFeature )
{
   rg_FLAG isAmidineChemType = rg_FALSE;
	
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;

   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;

   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );


   if( currAtom->getAtomCode() == C_ATOM && 
       numOfBondForCurrAtom    == 3      && numOfNitrogen           == 2      &&
       numOfSingleBond         == 2      && numOfDoubleBond         == 1       ) {

       Atom** arrDoubleBondedNitrogens = new Atom* [numOfDoubleBond];
       arrDoubleBondedNitrogens[0] = rg_NULL;
       getBondedAtomFromCurrAtom( currAtom, N_ATOM, DOUBLE_BOND, arrDoubleBondedNitrogens );	// "currAtom" \BF\A1 "DOUBLE_BOND"\B5\C8 "N_ATOM"\C0\BB ã\BEƼ\AD "arrDoubleBondedNitrogens"\BF\A1 \C0\FA\C0\E5\C7϶\F3\B4\C2 \C7Լ\F6
       
       if( arrDoubleBondedNitrogens[0] == rg_NULL ) {
           isAmidineChemType = rg_FALSE;
           delete [] arrDoubleBondedNitrogens;
           return isAmidineChemType;
       }


       Atom** arrSingleBondedNitrogens = new Atom* [numOfSingleBond];
       getBondedAtomFromCurrAtom( currAtom, N_ATOM, SINGLE_BOND, arrSingleBondedNitrogens ); //

       rg_INT  numOfBondsForDoubleBondedNitrogen = arrDoubleBondedNitrogens[0]->getListChemicalBond()->getSize();

       rg_INT  numOfSingleBondForDoubleBondedNitrogen;
       rg_INT  numOfDoubleBondForDoubleBondedNitrogen;
       rg_INT  numOfTripleBondForDoubleBondedNitrogen;
       rg_INT  numOfAromaticBondForDoubleBondedNitrogen;

       countEachBondTypeFromBondListInAtom( arrDoubleBondedNitrogens[0], numOfSingleBondForDoubleBondedNitrogen, 
                                                                         numOfDoubleBondForDoubleBondedNitrogen, 
                                                                         numOfTripleBondForDoubleBondedNitrogen, 
                                                                         numOfAromaticBondForDoubleBondedNitrogen );

		rg_FLAG isDoubleBondedNitrogenSatisfied = rg_FALSE;
		
       if( numOfBondsForDoubleBondedNitrogen == 2 && numOfSingleBondForDoubleBondedNitrogen == 1 && numOfDoubleBondForDoubleBondedNitrogen ==1 ) {
           isDoubleBondedNitrogenSatisfied = rg_TRUE;
       }
       else {
           isDoubleBondedNitrogenSatisfied = rg_FALSE;
       }


       rg_INT  numOfBondsForSingleBondedNitrogen = arrSingleBondedNitrogens[0]->getListChemicalBond()->getSize();
		rg_FLAG isAllBondsSingleForSingleBondedNitrogen = isAllBondsSingle( arrSingleBondedNitrogens[0] );
		
		rg_FLAG isSingleBondedNitrogenSatisfied = rg_FALSE;
		
       if( numOfBondsForSingleBondedNitrogen == 3 && isAllBondsSingleForSingleBondedNitrogen == rg_TRUE ) {
           isSingleBondedNitrogenSatisfied = rg_TRUE;
       }
       else {
           isSingleBondedNitrogenSatisfied = rg_FALSE;
       }


       if( isDoubleBondedNitrogenSatisfied == rg_TRUE && isSingleBondedNitrogenSatisfied == rg_TRUE ) {

           Atom* prevAtom = currAtom;

           Atom** arrBondedAtomsOfSingleBondedNitrogen = new Atom* [2]; // EXCEPT PREV ATOM

           getBondedAtomFromCurrAtomExceptPrevAtom( arrSingleBondedNitrogens[0], prevAtom, arrBondedAtomsOfSingleBondedNitrogen );

           if( arrBondedAtomsOfSingleBondedNitrogen[0]->getAtomCode() == C_ATOM &&
               arrBondedAtomsOfSingleBondedNitrogen[1]->getAtomCode() == C_ATOM   ) {
               
               rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogen[0], arrSingleBondedNitrogens[0] );
               rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogen[1], arrSingleBondedNitrogens[0] );
				
               if( isValidCarbonA == rg_TRUE && isValidCarbonB == rg_TRUE ) {
                   
					piAmidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, AMIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0] );
					isAmidineChemType = rg_TRUE;

                   delete [] arrBondedAtomsOfSingleBondedNitrogen;
                   delete [] arrDoubleBondedNitrogens;
                   delete [] arrSingleBondedNitrogens;

                   return isAmidineChemType;
               }
           }

			if( arrBondedAtomsOfSingleBondedNitrogen[0]->getAtomCode() == H_ATOM || 
               arrBondedAtomsOfSingleBondedNitrogen[1]->getAtomCode() == H_ATOM  ) {
                   
					piAmidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, AMIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0] );
					isAmidineChemType = rg_TRUE;
           }

           delete [] arrBondedAtomsOfSingleBondedNitrogen;
       }

       delete [] arrDoubleBondedNitrogens;
       delete [] arrSingleBondedNitrogens;
   }

   return isAmidineChemType;
}		



rg_FLAG Molecule::isPiGuanidine( Atom* currAtom, PharmaFeature& piGuanidineFeature )
{
   rg_FLAG isGuanidineChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;

   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;

   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );


   if( currAtom->getAtomCode() == C_ATOM && 
       numOfBondForCurrAtom    == 3      && numOfNitrogen           == 3      &&
       numOfSingleBond         == 2      && numOfDoubleBond         == 1       ) {

       Atom** arrDoubleBondedNitrogens = new Atom* [numOfDoubleBond];
       Atom** arrSingleBondedNitrogens = new Atom* [numOfSingleBond];
       
       getBondedAtomFromCurrAtom( currAtom, N_ATOM, DOUBLE_BOND, arrDoubleBondedNitrogens );	// "currAtom" \BF\A1 "DOUBLE_BOND"\B5\C8 "N_ATOM"\C0\BB ã\BEƼ\AD "arrDoubleBondedNitrogens"\BF\A1 \C0\FA\C0\E5\C7϶\F3\B4\C2 \C7Լ\F6
       getBondedAtomFromCurrAtom( currAtom, N_ATOM, SINGLE_BOND, arrSingleBondedNitrogens ); //

       rg_INT  numOfBondsForDoubleBondedNitrogen = arrDoubleBondedNitrogens[0]->getListChemicalBond()->getSize();

       rg_INT  numOfSingleBondForDoubleBondedNitrogen;
       rg_INT  numOfDoubleBondForDoubleBondedNitrogen;
       rg_INT  numOfTripleBondForDoubleBondedNitrogen;
       rg_INT  numOfAromaticBondForDoubleBondedNitrogen;

       countEachBondTypeFromBondListInAtom( arrDoubleBondedNitrogens[0], numOfSingleBondForDoubleBondedNitrogen, 
                                                                         numOfDoubleBondForDoubleBondedNitrogen, 
                                                                         numOfTripleBondForDoubleBondedNitrogen, 
                                                                         numOfAromaticBondForDoubleBondedNitrogen );

		rg_FLAG isDoubleBondedNitrogenSatisfied = rg_FALSE;

       if( numOfBondsForDoubleBondedNitrogen == 2 && numOfSingleBondForDoubleBondedNitrogen == 1 && numOfDoubleBondForDoubleBondedNitrogen ==1 ) {
           isDoubleBondedNitrogenSatisfied = rg_TRUE;
       }
       else {
           isDoubleBondedNitrogenSatisfied = rg_FALSE;
       }


		rg_FLAG isSingleBondedNitrogenSatisfied = rg_FALSE;

       rg_INT i_bondedNitrogen=0;
       for( i_bondedNitrogen=0; i_bondedNitrogen<numOfSingleBond; i_bondedNitrogen++ ) {
           rg_INT  numOfBondsForSingleBondedNitrogen       = arrSingleBondedNitrogens[i_bondedNitrogen]->getListChemicalBond()->getSize();
           rg_FLAG isAllBondsSingleForSingleBondedNitrogen = isAllBondsSingle( arrSingleBondedNitrogens[i_bondedNitrogen] );

           if( numOfBondsForSingleBondedNitrogen == 3 && isAllBondsSingleForSingleBondedNitrogen == rg_TRUE ) {
               isSingleBondedNitrogenSatisfied = rg_TRUE;
           }
           else {
               isSingleBondedNitrogenSatisfied = rg_FALSE;
				break;          
           }
       }

       if( isDoubleBondedNitrogenSatisfied == rg_TRUE && isSingleBondedNitrogenSatisfied == rg_TRUE ) {

           Atom* prevAtom = currAtom;

           Atom** arrBondedAtomsOfSingleBondedNitrogenA = new Atom* [2]; // EXCEPT PREV ATOM
           Atom** arrBondedAtomsOfSingleBondedNitrogenB = new Atom* [2]; // EXCEPT PREV ATOM

           getBondedAtomFromCurrAtomExceptPrevAtom( arrSingleBondedNitrogens[0], prevAtom, arrBondedAtomsOfSingleBondedNitrogenA );
           getBondedAtomFromCurrAtomExceptPrevAtom( arrSingleBondedNitrogens[1], prevAtom, arrBondedAtomsOfSingleBondedNitrogenB );

           if( arrBondedAtomsOfSingleBondedNitrogenA[0]->getAtomCode() == C_ATOM &&
               arrBondedAtomsOfSingleBondedNitrogenA[1]->getAtomCode() == C_ATOM && 
               arrBondedAtomsOfSingleBondedNitrogenB[0]->getAtomCode() == C_ATOM && 
               arrBondedAtomsOfSingleBondedNitrogenB[1]->getAtomCode() == C_ATOM   ) {
               
               rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenA[0], arrSingleBondedNitrogens[0] );
               rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenA[1], arrSingleBondedNitrogens[0] );

               rg_FLAG isValidCarbonC = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenB[0], arrSingleBondedNitrogens[1] );
               rg_FLAG isValidCarbonD = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenB[1], arrSingleBondedNitrogens[1] );

               if( isValidCarbonA == rg_TRUE && isValidCarbonB == rg_TRUE &&
                   isValidCarbonC == rg_TRUE && isValidCarbonD == rg_TRUE  ) {
                   
					piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0] );
					isGuanidineChemType = rg_TRUE;

                   delete [] arrDoubleBondedNitrogens;
                   delete [] arrSingleBondedNitrogens;

                   delete [] arrBondedAtomsOfSingleBondedNitrogenA;
                   delete [] arrBondedAtomsOfSingleBondedNitrogenB;

                   return isGuanidineChemType;
               }

           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[0]->getAtomCode() == C_ATOM &&
               arrBondedAtomsOfSingleBondedNitrogenA[1]->getAtomCode() == C_ATOM && 
               arrBondedAtomsOfSingleBondedNitrogenB[0]->getAtomCode() == H_ATOM   ) {

               rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenA[0], arrSingleBondedNitrogens[0] );
               rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenA[1], arrSingleBondedNitrogens[0] );

               if( isValidCarbonA == rg_TRUE && isValidCarbonB == rg_TRUE ) {

					piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[1] );
					isGuanidineChemType = rg_TRUE;

                   delete [] arrDoubleBondedNitrogens;
                   delete [] arrSingleBondedNitrogens;
                   
                   delete [] arrBondedAtomsOfSingleBondedNitrogenA;
                   delete [] arrBondedAtomsOfSingleBondedNitrogenB;

                   return isGuanidineChemType;
               }
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[0]->getAtomCode() == C_ATOM &&
               arrBondedAtomsOfSingleBondedNitrogenA[1]->getAtomCode() == C_ATOM &&  
               arrBondedAtomsOfSingleBondedNitrogenB[1]->getAtomCode() == H_ATOM   ) {
               
               rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenA[0], arrSingleBondedNitrogens[0] );
               rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenA[1], arrSingleBondedNitrogens[0] );

               if( isValidCarbonA == rg_TRUE && isValidCarbonB == rg_TRUE ) {
                   
					piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[1] );
					isGuanidineChemType = rg_TRUE;

                   delete [] arrDoubleBondedNitrogens;
                   delete [] arrSingleBondedNitrogens;
                   
                   delete [] arrBondedAtomsOfSingleBondedNitrogenA;
                   delete [] arrBondedAtomsOfSingleBondedNitrogenB;

                   return isGuanidineChemType;
               }
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[0]->getAtomCode() == H_ATOM &&
               arrBondedAtomsOfSingleBondedNitrogenB[0]->getAtomCode() == C_ATOM && 
               arrBondedAtomsOfSingleBondedNitrogenB[1]->getAtomCode() == C_ATOM   ) {
               
               rg_FLAG isValidCarbonC = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenB[0], arrSingleBondedNitrogens[1] );
               rg_FLAG isValidCarbonD = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenB[1], arrSingleBondedNitrogens[1] );

               if( isValidCarbonC == rg_TRUE && isValidCarbonD == rg_TRUE ) {

					piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0] );
					isGuanidineChemType = rg_TRUE;

                   delete [] arrDoubleBondedNitrogens;
                   delete [] arrSingleBondedNitrogens;
                   
                   delete [] arrBondedAtomsOfSingleBondedNitrogenA;
                   delete [] arrBondedAtomsOfSingleBondedNitrogenB;

                   return isGuanidineChemType;
               }
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[1]->getAtomCode() == H_ATOM && 
               arrBondedAtomsOfSingleBondedNitrogenB[0]->getAtomCode() == C_ATOM && 
               arrBondedAtomsOfSingleBondedNitrogenB[1]->getAtomCode() == C_ATOM   ) {
               
               rg_FLAG isValidCarbonC = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenB[0], arrSingleBondedNitrogens[1] );
               rg_FLAG isValidCarbonD = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfSingleBondedNitrogenB[1], arrSingleBondedNitrogens[1] );

               if( isValidCarbonC == rg_TRUE && isValidCarbonD == rg_TRUE  ) {

					piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0] );
					isGuanidineChemType = rg_TRUE;

                   delete [] arrDoubleBondedNitrogens;
                   delete [] arrSingleBondedNitrogens;
                   
                   delete [] arrBondedAtomsOfSingleBondedNitrogenA;
                   delete [] arrBondedAtomsOfSingleBondedNitrogenB;

                   return isGuanidineChemType;
               }
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[0]->getAtomCode() == H_ATOM &&  
               arrBondedAtomsOfSingleBondedNitrogenB[0]->getAtomCode() == H_ATOM   ) {

				piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0], arrSingleBondedNitrogens[1] );
				isGuanidineChemType = rg_TRUE;

               delete [] arrDoubleBondedNitrogens;
               delete [] arrSingleBondedNitrogens;
               
               delete [] arrBondedAtomsOfSingleBondedNitrogenA;
               delete [] arrBondedAtomsOfSingleBondedNitrogenB;

               return isGuanidineChemType;
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[1]->getAtomCode() == H_ATOM &&  
               arrBondedAtomsOfSingleBondedNitrogenB[1]->getAtomCode() == H_ATOM   ) {

				piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0], arrSingleBondedNitrogens[1] );
				isGuanidineChemType = rg_TRUE;

               delete [] arrDoubleBondedNitrogens;
               delete [] arrSingleBondedNitrogens;
               
               delete [] arrBondedAtomsOfSingleBondedNitrogenA;
               delete [] arrBondedAtomsOfSingleBondedNitrogenB;

               return isGuanidineChemType;
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[0]->getAtomCode() == H_ATOM &&  
               arrBondedAtomsOfSingleBondedNitrogenB[1]->getAtomCode() == H_ATOM   ) {

               piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0], arrSingleBondedNitrogens[1] );
				isGuanidineChemType = rg_TRUE;

               delete [] arrDoubleBondedNitrogens;
               delete [] arrSingleBondedNitrogens;
               
               delete [] arrBondedAtomsOfSingleBondedNitrogenA;
               delete [] arrBondedAtomsOfSingleBondedNitrogenB;

				return isGuanidineChemType;
           }

           if( arrBondedAtomsOfSingleBondedNitrogenA[1]->getAtomCode() == H_ATOM &&  
               arrBondedAtomsOfSingleBondedNitrogenB[0]->getAtomCode() == H_ATOM   ) {

               piGuanidineFeature.setTypesAndAtoms( PI_FEATRURE_TYPE, GUANIDINE_CHEM_TYPE, arrDoubleBondedNitrogens[0], arrSingleBondedNitrogens[0], arrSingleBondedNitrogens[1] );
				isGuanidineChemType = rg_TRUE;

               delete [] arrDoubleBondedNitrogens;
               delete [] arrSingleBondedNitrogens;
               
               delete [] arrBondedAtomsOfSingleBondedNitrogenA;
               delete [] arrBondedAtomsOfSingleBondedNitrogenB;

               return isGuanidineChemType;
           }

           delete [] arrBondedAtomsOfSingleBondedNitrogenA;
           delete [] arrBondedAtomsOfSingleBondedNitrogenB;
       }
       
       delete [] arrDoubleBondedNitrogens;
       delete [] arrSingleBondedNitrogens;
   }

   return isGuanidineChemType;
}



void Molecule::defineNIPharmaFeatures()
{
   Atom*  currAtom        = rg_NULL;
   rg_INT pharmaFeatureID = m_pharmaFeatures.size();
	
   m_atoms.reset4Loop();
   while( m_atoms.setNext4Loop() ) {
		
       currAtom = m_atoms.getpEntity();
       
       PharmaFeature  aNiNegativeChargeFeature;
       if( isNiNegativeCharge( currAtom, aNiNegativeChargeFeature ) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiNegativeChargeFeature);
       }
		
       PharmaFeature  aNiCarboxylFeature;
       if( isNiCarboxyl( currAtom, aNiCarboxylFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiCarboxylFeature);
       }

       PharmaFeature  aNiTrifluFeature;
       if( isNiTriflu( currAtom, aNiTrifluFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiTrifluFeature);
       }
       
       PharmaFeature  aNiSulfinicFeature;
       if( isNiSulfinic( currAtom, aNiSulfinicFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiSulfinicFeature);
       }

       PharmaFeature  aNiSulfonicFeature;
       if( isNiSulfonic( currAtom, aNiSulfonicFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiSulfonicFeature);
       }

       PharmaFeature  aNiSulfuricFeature;
       if( isNiSulfuric( currAtom, aNiSulfuricFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiSulfuricFeature);
       }

       PharmaFeature  aNiPhosphinicFeature;
       if( isNiPhosphinic( currAtom, aNiPhosphinicFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiPhosphinicFeature);
       }

       PharmaFeature  aNiPhosphonicFeature;
       if( isNiPhosphonic( currAtom, aNiPhosphonicFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiPhosphonicFeature);
       }

       PharmaFeature  aNiPhosphoricFeature;
       if( isNiPhosphoric( currAtom, aNiPhosphoricFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aNiPhosphoricFeature);
       }
   }	
}



rg_FLAG Molecule::isNiNegativeCharge( Atom* currAtom, PharmaFeature& niNegativeChargeFeature )
{
   rg_FLAG isNegChargeChemType = rg_FALSE;
   
   rg_REAL chargeOfCurrAtom     = currAtom->getpChemicalProperties()->getCharge();
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   if( chargeOfCurrAtom <= -1.0 ) {
       if( isAllBondedAtomsNonPositive( currAtom ) == rg_TRUE ) {
           isNegChargeChemType = rg_TRUE;
           niNegativeChargeFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, NEG_CHARGE_CHEM_TYPE, currAtom );
           return isNegChargeChemType;
       }
   }
   
   if( currAtom->getAtomCode() == O_ATOM ) {
       if( numOfBondForCurrAtom == 1 && numOfSingleBond == 1 ) {
           if( isAllBondedAtomsNonPositive( currAtom ) == rg_TRUE ) {
               isNegChargeChemType = rg_TRUE;
               niNegativeChargeFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, NEG_CHARGE_CHEM_TYPE, currAtom );
           }
       }
   }
   
   return isNegChargeChemType;
}



rg_FLAG Molecule::isNiCarboxyl( Atom* currAtom, PharmaFeature& niCarboxylFeature )
{
   rg_FLAG isCarboxylChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;
   
   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );
   
   
   if( currAtom->getAtomCode() == C_ATOM && 
       numOfBondForCurrAtom    == 3      && numOfOxygen             == 2      &&
       numOfSingleBond         == 2      && numOfDoubleBond         == 1       ) {
       
       Atom** arrDoubleBondedOxygens = new Atom* [numOfDoubleBond];
       arrDoubleBondedOxygens[0] = rg_NULL;
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrDoubleBondedOxygens );	

       if( arrDoubleBondedOxygens[0] == rg_NULL ) {
           isCarboxylChemType = rg_FALSE;
           delete [] arrDoubleBondedOxygens;
           return isCarboxylChemType;
       }

       Atom** arrSingleBondedOxygens = new Atom* [numOfSingleBond];
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrSingleBondedOxygens ); 
       

       rg_INT  numOfBondsForDoubleBondedOxygen = arrDoubleBondedOxygens[0]->getListChemicalBond()->getSize();
       
       rg_INT  numOfSingleBondForDoubleBondedOxygen;
       rg_INT  numOfDoubleBondForDoubleBondedOxygen;
       rg_INT  numOfTripleBondForDoubleBondedOxygen;
       rg_INT  numOfAromaticBondForDoubleBondedOxygen;
       
       countEachBondTypeFromBondListInAtom( arrDoubleBondedOxygens[0], numOfSingleBondForDoubleBondedOxygen, 
                                                                       numOfDoubleBondForDoubleBondedOxygen, 
                                                                       numOfTripleBondForDoubleBondedOxygen, 
                                                                       numOfAromaticBondForDoubleBondedOxygen );
       
       rg_FLAG isDoubleBondedOxygenSatisfied = rg_FALSE;

       if( numOfBondsForDoubleBondedOxygen == 1 && numOfDoubleBondForDoubleBondedOxygen == 1 ) {   // BONDED WITH CURR C_ATOM
           isDoubleBondedOxygenSatisfied = rg_TRUE;
       }
       else {
           isDoubleBondedOxygenSatisfied = rg_FALSE;
       }
       
       rg_INT  numOfBondsForSingleBondedOxygen = arrSingleBondedOxygens[0]->getListChemicalBond()->getSize();
       rg_FLAG isAllBondsSingleForSingleBondedOxygen = isAllBondsSingle( arrSingleBondedOxygens[0] );
       
       rg_FLAG isSingleBondedOxygenSatisfied = rg_FALSE;
       
       if( numOfBondsForSingleBondedOxygen == 2 && isAllBondsSingleForSingleBondedOxygen == rg_TRUE ) {
           isSingleBondedOxygenSatisfied = rg_TRUE;
       }
       else {
           isSingleBondedOxygenSatisfied = rg_FALSE;
       }
       

       if( isDoubleBondedOxygenSatisfied == rg_TRUE && isSingleBondedOxygenSatisfied == rg_TRUE ) {
           
           Atom* prevAtom = currAtom;
           
           Atom** arrBondedAtomsOfSingleBondedOxygen = new Atom* [1]; // EXCEPT PREV ATOM
           
           getBondedAtomFromCurrAtomExceptPrevAtom( arrSingleBondedOxygens[0], prevAtom, arrBondedAtomsOfSingleBondedOxygen );
           
           if( arrBondedAtomsOfSingleBondedOxygen[0]->getAtomCode() == H_ATOM  ) {                
               niCarboxylFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, CARBOXYL_CHEM_TYPE, arrDoubleBondedOxygens[0], arrSingleBondedOxygens[0] );
               isCarboxylChemType = rg_TRUE;
           }

           delete [] arrBondedAtomsOfSingleBondedOxygen;
           
       }
       
       delete [] arrDoubleBondedOxygens;
       delete [] arrSingleBondedOxygens;
   }
   
   return isCarboxylChemType;
}



rg_FLAG Molecule::isNiTriflu( Atom* currAtom, PharmaFeature& niTrifluFeature )
{
   rg_FLAG isTrifluChemType = rg_FALSE;
	
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;

   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;

   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );


   if( currAtom->getAtomCode() == N_ATOM && numOfBondForCurrAtom    == 3  && numOfSingleBond    == 3  &&
       numOfCarbon             == 1      && numOfHydrogen           == 1  && numOfSulfur        == 1   ) {

       Atom** bondedCarbon = new Atom*;
       getBondedAtomFromCurrAtom( currAtom, C_ATOM, SINGLE_BOND, bondedCarbon );

       rg_INT  numOfBondsForBondedCarbon       = (*bondedCarbon)->getListChemicalBond()->getSize();
       rg_FLAG isAllBondsSingleForBondedCarbon = isAllBondsSingle( *bondedCarbon );
       
       if( numOfBondsForBondedCarbon == 4 && isAllBondsSingleForBondedCarbon == rg_TRUE ) {
           
           Atom*  prevAtom = currAtom;
           
           Atom** arrBondedAtomsOfBondedCarbons = new Atom* [3]; // EXCEPT PREV ATOM
           
           getBondedAtomFromCurrAtomExceptPrevAtom( *bondedCarbon, prevAtom, arrBondedAtomsOfBondedCarbons );
           
           if( arrBondedAtomsOfBondedCarbons[0]->getAtomCode() == F_ATOM &&
               arrBondedAtomsOfBondedCarbons[1]->getAtomCode() == F_ATOM &&
               arrBondedAtomsOfBondedCarbons[2]->getAtomCode() == F_ATOM    ) {

               Atom** bondedSulfur = new Atom*;        
               getBondedAtomFromCurrAtom( currAtom, S_ATOM, SINGLE_BOND, bondedSulfur );

               rg_INT  numOfBondsForBondedSulfur = (*bondedSulfur)->getListChemicalBond()->getSize();
       
               rg_INT  numOfSingleBondForBondedSulfur;
               rg_INT  numOfDoubleBondForBondedSulfur;
               rg_INT  numOfTripleBondForBondedSulfur;
               rg_INT  numOfAromaticBondForBondedSulfur;
       
               countEachBondTypeFromBondListInAtom( *bondedSulfur, numOfSingleBondForBondedSulfur, numOfDoubleBondForBondedSulfur, numOfTripleBondForBondedSulfur, numOfAromaticBondForBondedSulfur );
       
               if( numOfBondsForBondedSulfur == 4 && numOfSingleBondForBondedSulfur == 2 && numOfDoubleBondForBondedSulfur == 2 ) {
           
                   Atom** arrBondedAtomsOfBondedSulfur = new Atom* [3]; // EXCEPT PREV ATOM

                   getBondedAtomFromCurrAtomExceptPrevAtom( *bondedSulfur, prevAtom, arrBondedAtomsOfBondedSulfur );

                   rg_INT  numOfBondsForBondedAtomsOfBondedSulfurA = arrBondedAtomsOfBondedSulfur[0]->getListChemicalBond()->getSize();
                   
                   rg_INT  numOfSingleBondForBondedAtomsOfBondedSulfurA;
                   rg_INT  numOfDoubleBondForBondedAtomsOfBondedSulfurA;
                   rg_INT  numOfTripleBondForBondedAtomsOfBondedSulfurA;
                   rg_INT  numOfAromaticBondForBondedAtomsOfBondedSulfurA;
                   
                   countEachBondTypeFromBondListInAtom( arrBondedAtomsOfBondedSulfur[0], numOfSingleBondForBondedAtomsOfBondedSulfurA, numOfDoubleBondForBondedAtomsOfBondedSulfurA, numOfTripleBondForBondedAtomsOfBondedSulfurA, numOfAromaticBondForBondedAtomsOfBondedSulfurA );
                   
                   rg_INT  numOfBondsForBondedAtomsOfBondedSulfurB = arrBondedAtomsOfBondedSulfur[1]->getListChemicalBond()->getSize();
                   
                   rg_INT  numOfSingleBondForBondedAtomsOfBondedSulfurB;
                   rg_INT  numOfDoubleBondForBondedAtomsOfBondedSulfurB;
                   rg_INT  numOfTripleBondForBondedAtomsOfBondedSulfurB;
                   rg_INT  numOfAromaticBondForBondedAtomsOfBondedSulfurB;
                   
                   countEachBondTypeFromBondListInAtom( arrBondedAtomsOfBondedSulfur[1], numOfSingleBondForBondedAtomsOfBondedSulfurB, numOfDoubleBondForBondedAtomsOfBondedSulfurB, numOfTripleBondForBondedAtomsOfBondedSulfurB, numOfAromaticBondForBondedAtomsOfBondedSulfurB );

                   rg_INT  numOfBondsForBondedAtomsOfBondedSulfurC = arrBondedAtomsOfBondedSulfur[2]->getListChemicalBond()->getSize();
                   
                   rg_INT  numOfSingleBondForBondedAtomsOfBondedSulfurC;
                   rg_INT  numOfDoubleBondForBondedAtomsOfBondedSulfurC;
                   rg_INT  numOfTripleBondForBondedAtomsOfBondedSulfurC;
                   rg_INT  numOfAromaticBondForBondedAtomsOfBondedSulfurC;
                   
                   countEachBondTypeFromBondListInAtom( arrBondedAtomsOfBondedSulfur[2], numOfSingleBondForBondedAtomsOfBondedSulfurC, numOfDoubleBondForBondedAtomsOfBondedSulfurC, numOfTripleBondForBondedAtomsOfBondedSulfurC, numOfAromaticBondForBondedAtomsOfBondedSulfurC );

                   if( arrBondedAtomsOfBondedSulfur[0]->getAtomCode() == O_ATOM && arrBondedAtomsOfBondedSulfur[1]->getAtomCode() == O_ATOM && 
                       numOfBondsForBondedAtomsOfBondedSulfurA == 1             && numOfDoubleBondForBondedAtomsOfBondedSulfurA == 1        &&
                       numOfBondsForBondedAtomsOfBondedSulfurB == 1             && numOfDoubleBondForBondedAtomsOfBondedSulfurB == 1         ) {

                       niTrifluFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, TRIFLU_CHEM_TYPE, currAtom );
                       isTrifluChemType = rg_TRUE;

                       delete [] arrBondedAtomsOfBondedSulfur;
                       delete bondedSulfur;
                       delete [] arrBondedAtomsOfBondedCarbons;
                       delete bondedCarbon;
                                                  
                       return isTrifluChemType;
                   }

                   if( arrBondedAtomsOfBondedSulfur[0]->getAtomCode() == O_ATOM && arrBondedAtomsOfBondedSulfur[2]->getAtomCode() == O_ATOM &&
                       numOfBondsForBondedAtomsOfBondedSulfurA      == 1        && numOfDoubleBondForBondedAtomsOfBondedSulfurA == 1        &&
                       numOfBondsForBondedAtomsOfBondedSulfurC      == 1        && numOfDoubleBondForBondedAtomsOfBondedSulfurC == 1         ) {
                       
                       niTrifluFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, TRIFLU_CHEM_TYPE, currAtom );
                       isTrifluChemType = rg_TRUE;

                       delete [] arrBondedAtomsOfBondedSulfur;
                       delete bondedSulfur;
                       delete [] arrBondedAtomsOfBondedCarbons;
                       delete bondedCarbon;
                       
                       return isTrifluChemType;
                   }
                   
                   if( arrBondedAtomsOfBondedSulfur[1]->getAtomCode() == O_ATOM && arrBondedAtomsOfBondedSulfur[2]->getAtomCode() == O_ATOM &&
                       numOfBondsForBondedAtomsOfBondedSulfurB      == 1        && numOfDoubleBondForBondedAtomsOfBondedSulfurB == 1        &&
                       numOfBondsForBondedAtomsOfBondedSulfurC      == 1        && numOfDoubleBondForBondedAtomsOfBondedSulfurC == 1         ) {
                       
                       niTrifluFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, TRIFLU_CHEM_TYPE, currAtom );
                       isTrifluChemType = rg_TRUE;

                       delete [] arrBondedAtomsOfBondedSulfur;
                       delete bondedSulfur;
                       delete [] arrBondedAtomsOfBondedCarbons;
                       delete bondedCarbon;
                       
                       return isTrifluChemType;
                   }

                   delete [] arrBondedAtomsOfBondedSulfur;
               }
               
               delete bondedSulfur;
           }

           delete [] arrBondedAtomsOfBondedCarbons;
       }

       delete bondedCarbon;
   }

   return isTrifluChemType;
}



rg_FLAG Molecule::isNiSulfinic( Atom* currAtom, PharmaFeature& niSulfinicFeature )
{
   rg_FLAG isSulfinicChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;
   
   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );
   
   
   if( currAtom->getAtomCode() == S_ATOM && 
       numOfBondForCurrAtom    == 3      && numOfOxygen             == 2      &&
       numOfSingleBond         == 2      && numOfDoubleBond         == 1       ) {
       
       Atom** arrDoubleBondedOxygens = new Atom* [numOfDoubleBond];
       arrDoubleBondedOxygens[0] = rg_NULL;
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrDoubleBondedOxygens );

       if( arrDoubleBondedOxygens[0] == rg_NULL ) {
           isSulfinicChemType = rg_FALSE;
           delete [] arrDoubleBondedOxygens;
           return isSulfinicChemType;
       }


       Atom** arrSingleBondedOxygens = new Atom* [numOfSingleBond];
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrSingleBondedOxygens ); 
       
       rg_INT  numOfBondsForDoubleBondedOxygen = arrDoubleBondedOxygens[0]->getListChemicalBond()->getSize();
       
       rg_INT  numOfSingleBondForDoubleBondedOxygen;
       rg_INT  numOfDoubleBondForDoubleBondedOxygen;
       rg_INT  numOfTripleBondForDoubleBondedOxygen;
       rg_INT  numOfAromaticBondForDoubleBondedOxygen;
       
       countEachBondTypeFromBondListInAtom( arrDoubleBondedOxygens[0], numOfSingleBondForDoubleBondedOxygen, 
                                                                       numOfDoubleBondForDoubleBondedOxygen, 
                                                                       numOfTripleBondForDoubleBondedOxygen, 
                                                                       numOfAromaticBondForDoubleBondedOxygen );
       
       rg_FLAG isDoubleBondedOxygenSatisfied = rg_FALSE;
       
       if( numOfBondsForDoubleBondedOxygen == 1 && numOfDoubleBondForDoubleBondedOxygen == 1 ) {   // BONDED WITH CURR C_ATOM
           isDoubleBondedOxygenSatisfied = rg_TRUE;
       }
       else {
           isDoubleBondedOxygenSatisfied = rg_FALSE;
       }
       
       rg_INT  numOfBondsForSingleBondedOxygen = arrSingleBondedOxygens[0]->getListChemicalBond()->getSize();
       rg_FLAG isAllBondsSingleForSingleBondedOxygen = isAllBondsSingle( arrSingleBondedOxygens[0] );
       
       rg_FLAG isSingleBondedOxygenSatisfied = rg_FALSE;
       
       if( numOfBondsForSingleBondedOxygen == 2 && isAllBondsSingleForSingleBondedOxygen == rg_TRUE ) {
           isSingleBondedOxygenSatisfied = rg_TRUE;
       }
       else {
           isSingleBondedOxygenSatisfied = rg_FALSE;
       }
       
       
       if( isDoubleBondedOxygenSatisfied == rg_TRUE && isSingleBondedOxygenSatisfied == rg_TRUE ) {
           
           Atom* prevAtom = currAtom;
           
           Atom** arrBondedAtomsOfSingleBondedOxygen = new Atom* [1]; // EXCEPT PREV ATOM
           
           getBondedAtomFromCurrAtomExceptPrevAtom( arrSingleBondedOxygens[0], prevAtom, arrBondedAtomsOfSingleBondedOxygen );
           
           if( arrBondedAtomsOfSingleBondedOxygen[0]->getAtomCode() == H_ATOM  ) {                
               niSulfinicFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, SULFINIC_CHEM_TYPE, arrDoubleBondedOxygens[0], arrSingleBondedOxygens[0] );
               isSulfinicChemType = rg_TRUE;
           }
           
           delete [] arrBondedAtomsOfSingleBondedOxygen;
           
       }
       
       delete [] arrDoubleBondedOxygens;
       delete [] arrSingleBondedOxygens;
   }
   
   return isSulfinicChemType;
}



rg_FLAG Molecule::isNiSulfonic( Atom* currAtom, PharmaFeature& niSulfonicFeature ){

   rg_FLAG isSulfonicChemType = rg_FALSE;
	
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;

   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;

   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );


   if( currAtom->getAtomCode() == S_ATOM && 
       numOfBondForCurrAtom    == 4      && numOfOxygen     == 3      &&   numOfCarbon  == 1  &&
       numOfSingleBond         == 2      && numOfDoubleBond == 2       ) {

       Atom** arrCurrDoubleBondedOxygen = new Atom* [2];
       Atom** arrCurrSingleBondedOxygen = new Atom* [1];
       Atom** arrCurrSingleBondedCarbon = new Atom* [1];
       
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrCurrDoubleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrCurrSingleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, C_ATOM, SINGLE_BOND, arrCurrSingleBondedCarbon );

       rg_INT  numOfBondsForCurrDoubleBondedOxygenA = arrCurrDoubleBondedOxygen[0]->getListChemicalBond()->getSize();
       rg_INT  numOfBondsForCurrDoubleBondedOxygenB = arrCurrDoubleBondedOxygen[1]->getListChemicalBond()->getSize();

       if( numOfBondsForCurrDoubleBondedOxygenA == 1 && numOfBondsForCurrDoubleBondedOxygenB == 1) {

           rg_INT  numOfBondsForCurrSingleBondedOxygen = arrCurrSingleBondedOxygen[0]->getListChemicalBond()->getSize();

           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygen = isAllBondsSingle( arrCurrSingleBondedOxygen[0] );
           
           rg_INT  numOfCarbonForCurrSingleBondedOxygen;
           rg_INT  numOfHydrogenForCurrSingleBondedOxygen;
           rg_INT  numOfNitrogenForCurrSingleBondedOxygen;
           rg_INT  numOfOxygenForCurrSingleBondedOxygen;
           rg_INT  numOfPhosphorusForCurrSingleBondedOxygen;
           rg_INT  numOfSulfurForCurrSingleBondedOxygen;
           
           countEachAtomTypeFromBondListInAtom( arrCurrSingleBondedOxygen[0], numOfCarbonForCurrSingleBondedOxygen, 
                                                                              numOfHydrogenForCurrSingleBondedOxygen, 
                                                                              numOfNitrogenForCurrSingleBondedOxygen, 
                                                                              numOfOxygenForCurrSingleBondedOxygen, 
                                                                              numOfPhosphorusForCurrSingleBondedOxygen, 
                                                                              numOfSulfurForCurrSingleBondedOxygen );
           
           if( numOfBondsForCurrSingleBondedOxygen    == 2 && isAllBondsSingleForCurrSingleBondedOxygen == rg_TRUE &&
               numOfHydrogenForCurrSingleBondedOxygen == 1    ) {

               rg_FLAG isValidCurrSingleBondedCarbon = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrCurrSingleBondedCarbon[0], currAtom );

               if( isValidCurrSingleBondedCarbon == rg_TRUE) {

                   niSulfonicFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, SULFONIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrDoubleBondedOxygen[1], arrCurrSingleBondedOxygen[0] );
					isSulfonicChemType = rg_TRUE;

               }
           }
       }
   }

   return isSulfonicChemType;
}



rg_FLAG Molecule::isNiSulfuric( Atom* currAtom, PharmaFeature& niSulfuricFeature ){
   
   rg_FLAG isSulfuricChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;
   
   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );
   
   
   if( currAtom->getAtomCode() == S_ATOM && 
       numOfBondForCurrAtom    == 4      && numOfOxygen     == 4  && 
       numOfSingleBond         == 2      && numOfDoubleBond == 2     ) {
       
       Atom** arrCurrDoubleBondedOxygen = new Atom* [2];
       Atom** arrCurrSingleBondedOxygen = new Atom* [2];
       
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrCurrDoubleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrCurrSingleBondedOxygen );
       
       rg_INT  numOfBondsForCurrDoubleBondedOxygenA = arrCurrDoubleBondedOxygen[0]->getListChemicalBond()->getSize();
       rg_INT  numOfBondsForCurrDoubleBondedOxygenB = arrCurrDoubleBondedOxygen[1]->getListChemicalBond()->getSize();
       
       if( numOfBondsForCurrDoubleBondedOxygenA == 1 && numOfBondsForCurrDoubleBondedOxygenB == 1) {
           
           rg_INT  numOfBondsForCurrSingleBondedOxygenA = arrCurrSingleBondedOxygen[0]->getListChemicalBond()->getSize();
           rg_INT  numOfBondsForCurrSingleBondedOxygenB = arrCurrSingleBondedOxygen[1]->getListChemicalBond()->getSize();
           
           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenA = isAllBondsSingle( arrCurrSingleBondedOxygen[0] );
           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenB = isAllBondsSingle( arrCurrSingleBondedOxygen[1] );
          
           if( numOfBondsForCurrSingleBondedOxygenA       == 2       && numOfBondsForCurrSingleBondedOxygenB       == 2      &&
               isAllBondsSingleForCurrSingleBondedOxygenA == rg_TRUE && isAllBondsSingleForCurrSingleBondedOxygenB == rg_TRUE  ) {

               Atom* prevAtom = currAtom;
               
               Atom** arrBondedAtomsOfCurrSingleBondedOxygenA = new Atom* [1]; // EXCEPT PREV ATOM
               Atom** arrBondedAtomsOfCurrSingleBondedOxygenB = new Atom* [1];
               
               getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[0], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenA );
               getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[1], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenB );

               if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == C_ATOM &&
                   arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == H_ATOM    ) {

                   rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenA[0], 
                                                                                                   arrCurrSingleBondedOxygen[0] );

                   if( isValidCarbonA == rg_TRUE) {

                       niSulfuricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, SULFURIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrDoubleBondedOxygen[1], arrCurrSingleBondedOxygen[1] );
                       isSulfuricChemType = rg_TRUE;

                   }

               }

               else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == H_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == C_ATOM    ) {
                        
                        rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenB[0], 
                                                                                                        arrCurrSingleBondedOxygen[1] );
                   
                        if( isValidCarbonA == rg_TRUE) {
                       
                           niSulfuricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, SULFURIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrDoubleBondedOxygen[1], arrCurrSingleBondedOxygen[0] );
                           isSulfuricChemType = rg_TRUE;
                       
                       }
                   
               }
           }
       }
             
   }
   
   return isSulfuricChemType;
}



rg_FLAG Molecule::isNiPhosphinic( Atom* currAtom, PharmaFeature& niPhosphinicFeature ){
   
   rg_FLAG isPhosphinicChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;
   
   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );
   
   
   if( currAtom->getAtomCode() == P_ATOM && 
       numOfBondForCurrAtom    == 4      && numOfOxygen     == 2  &&  numOfCarbon  == 2  &&
       numOfSingleBond         == 3      && numOfDoubleBond == 1                            ) {
       
       Atom** arrCurrDoubleBondedOxygen = new Atom* [1];
       Atom** arrCurrSingleBondedOxygen = new Atom* [1];
       Atom** arrCurrSingleBondedCarbon = new Atom* [2];
       
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrCurrDoubleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrCurrSingleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, C_ATOM, SINGLE_BOND, arrCurrSingleBondedCarbon );
       
       rg_INT  numOfBondsForCurrDoubleBondedOxygen = arrCurrDoubleBondedOxygen[0]->getListChemicalBond()->getSize();
       
       if( numOfBondsForCurrDoubleBondedOxygen == 1 ) {
           
           rg_INT  numOfBondsForCurrSingleBondedOxygen = arrCurrSingleBondedOxygen[0]->getListChemicalBond()->getSize();
           
           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygen = isAllBondsSingle( arrCurrSingleBondedOxygen[0] );
           
           if( numOfBondsForCurrSingleBondedOxygen == 2  && isAllBondsSingleForCurrSingleBondedOxygen == rg_TRUE   ) {
               
               Atom* prevAtom = currAtom;
               
               Atom** arrBondedAtomsOfCurrSingleBondedOxygen = new Atom* [1]; // EXCEPT PREV ATOM
               
               getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[0], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygen );
               
               if( arrBondedAtomsOfCurrSingleBondedOxygen[0]->getAtomCode() == H_ATOM ) {
                   
                   rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrCurrSingleBondedCarbon[0], currAtom );
                   rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrCurrSingleBondedCarbon[1], currAtom );
                   
                   if( isValidCarbonA == rg_TRUE && isValidCarbonB == rg_TRUE ) {
                       
                       niPhosphinicFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHINIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[0] );
                       isPhosphinicChemType = rg_TRUE;
                       
                   }
                   
               }
               
           }
       }
       
   }
   
   return isPhosphinicChemType;
}



rg_FLAG Molecule::isNiPhosphonic( Atom* currAtom, PharmaFeature& niPhosphonicFeature ){
   
   rg_FLAG isPhosphonicChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;
   
   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );
   
   
   if( currAtom->getAtomCode() == P_ATOM && 
       numOfBondForCurrAtom    == 4      && numOfOxygen     == 3  &&  numOfCarbon  == 1 &&
       numOfSingleBond         == 3      && numOfDoubleBond == 1                           ) {
       
       Atom** arrCurrDoubleBondedOxygen = new Atom* [1];
       Atom** arrCurrSingleBondedOxygen = new Atom* [2];
       Atom** arrCurrSingleBondedCarbon = new Atom* [1];
       
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrCurrDoubleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrCurrSingleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, C_ATOM, SINGLE_BOND, arrCurrSingleBondedCarbon );
       
       rg_INT  numOfBondsForCurrDoubleBondedOxygen = arrCurrDoubleBondedOxygen[0]->getListChemicalBond()->getSize();
       
       if( numOfBondsForCurrDoubleBondedOxygen == 1 ) {

           rg_FLAG isValidCarbon = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrCurrSingleBondedCarbon[0], currAtom );

           if( isValidCarbon == rg_TRUE ){

               rg_INT  numOfBondsForCurrSingleBondedOxygenA = arrCurrSingleBondedOxygen[0]->getListChemicalBond()->getSize();
               rg_INT  numOfBondsForCurrSingleBondedOxygenB = arrCurrSingleBondedOxygen[1]->getListChemicalBond()->getSize();
               
               rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenA = isAllBondsSingle( arrCurrSingleBondedOxygen[0] );
               rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenB = isAllBondsSingle( arrCurrSingleBondedOxygen[1] );

               if( numOfBondsForCurrSingleBondedOxygenA       == 2       && numOfBondsForCurrSingleBondedOxygenB       == 2      &&
                   isAllBondsSingleForCurrSingleBondedOxygenA == rg_TRUE && isAllBondsSingleForCurrSingleBondedOxygenB == rg_TRUE  ) {

                   Atom* prevAtom = currAtom;
                   
                   Atom** arrBondedAtomsOfCurrSingleBondedOxygenA = new Atom* [1]; // EXCEPT PREV ATOM
                   Atom** arrBondedAtomsOfCurrSingleBondedOxygenB = new Atom* [1];
                   
                   getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[0], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenA );
                   getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[1], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenB );

                   if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == H_ATOM &&
                       arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == H_ATOM    ) {

                       niPhosphonicFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHONIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[0], arrCurrSingleBondedOxygen[1] );
                                         isPhosphonicChemType = rg_TRUE;

                                         return isPhosphonicChemType;
                   }

                   else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == C_ATOM &&
                            arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == H_ATOM    ) {

                            rg_FLAG isValidCarbon = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenA[0], 
                                                                                                           arrCurrSingleBondedOxygen[0] );

                            if( isValidCarbon == rg_TRUE ) {

                                niPhosphonicFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHONIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[1] );
                                isPhosphonicChemType = rg_TRUE;
                           
                                return isPhosphonicChemType;
                            }                        
                       
                   }

                   else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == H_ATOM &&
                            arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == C_ATOM    ) {
                       
                            rg_FLAG isValidCarbon = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenB[0], 
                                                                                                           arrCurrSingleBondedOxygen[0] );
                       
                            if( isValidCarbon == rg_TRUE ) {
                           
                                niPhosphonicFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHONIC_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[0] );
                                isPhosphonicChemType = rg_TRUE;
                           
                                return isPhosphonicChemType;
                            }                          
                   }
               }
           }
       }
   }
   
   return isPhosphonicChemType;
}



rg_FLAG Molecule::isNiPhosphoric( Atom* currAtom, PharmaFeature& niPhosphoricFeature ){
   
   rg_FLAG isPhosphoricChemType = rg_FALSE;
   
   rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
   
   rg_INT  numOfSingleBond;
   rg_INT  numOfDoubleBond;
   rg_INT  numOfTripleBond;
   rg_INT  numOfAromaticBond;
   
   countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
   
   rg_INT  numOfCarbon;
   rg_INT  numOfHydrogen;
   rg_INT  numOfNitrogen;
   rg_INT  numOfOxygen;
   rg_INT  numOfPhosphorus;
   rg_INT  numOfSulfur;
   
   countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );
   
   
   if( currAtom->getAtomCode() == P_ATOM && 
       numOfBondForCurrAtom    == 4      && numOfOxygen     == 4  && 
       numOfSingleBond         == 3      && numOfDoubleBond == 1     ) {
       
       Atom** arrCurrDoubleBondedOxygen = new Atom* [1];
       Atom** arrCurrSingleBondedOxygen = new Atom* [3];
       
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, DOUBLE_BOND, arrCurrDoubleBondedOxygen );
       getBondedAtomFromCurrAtom( currAtom, O_ATOM, SINGLE_BOND, arrCurrSingleBondedOxygen );
       
       rg_INT  numOfBondsForCurrDoubleBondedOxygen = arrCurrDoubleBondedOxygen[0]->getListChemicalBond()->getSize();
       
       if( numOfBondsForCurrDoubleBondedOxygen == 1 ) {
           
           rg_INT  numOfBondsForCurrSingleBondedOxygenA = arrCurrSingleBondedOxygen[0]->getListChemicalBond()->getSize();
           rg_INT  numOfBondsForCurrSingleBondedOxygenB = arrCurrSingleBondedOxygen[1]->getListChemicalBond()->getSize();
           rg_INT  numOfBondsForCurrSingleBondedOxygenC = arrCurrSingleBondedOxygen[2]->getListChemicalBond()->getSize();
           
           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenA = isAllBondsSingle( arrCurrSingleBondedOxygen[0] );
           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenB = isAllBondsSingle( arrCurrSingleBondedOxygen[1] );
           rg_FLAG isAllBondsSingleForCurrSingleBondedOxygenC = isAllBondsSingle( arrCurrSingleBondedOxygen[2] );
           
           if( numOfBondsForCurrSingleBondedOxygenA       == 2       &&  isAllBondsSingleForCurrSingleBondedOxygenA == rg_TRUE && 
               numOfBondsForCurrSingleBondedOxygenB       == 2       &&  isAllBondsSingleForCurrSingleBondedOxygenB == rg_TRUE &&
               numOfBondsForCurrSingleBondedOxygenC       == 2       &&  isAllBondsSingleForCurrSingleBondedOxygenC == rg_TRUE   ) {
               
               Atom* prevAtom = currAtom;
               
               Atom** arrBondedAtomsOfCurrSingleBondedOxygenA = new Atom* [1]; // EXCEPT PREV ATOM
               Atom** arrBondedAtomsOfCurrSingleBondedOxygenB = new Atom* [1];
               Atom** arrBondedAtomsOfCurrSingleBondedOxygenC = new Atom* [1];
               
               getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[0], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenA );
               getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[1], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenB );
               getBondedAtomFromCurrAtomExceptPrevAtom( arrCurrSingleBondedOxygen[2], prevAtom, arrBondedAtomsOfCurrSingleBondedOxygenC );
               
               if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == H_ATOM &&
                   arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == H_ATOM &&
                   arrBondedAtomsOfCurrSingleBondedOxygenC[0]->getAtomCode() == C_ATOM    ) {
                   
                   rg_FLAG isValidCarbonC = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenC[0], 
                                                                                                  arrCurrSingleBondedOxygen[2] );
                   
                   if( isValidCarbonC == rg_TRUE ) {
                       
                       niPhosphoricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHORIC_ESTER_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[0], arrCurrSingleBondedOxygen[1] );
                       isPhosphoricChemType = rg_TRUE;

                       return isPhosphoricChemType;
                       
                   }
                   
               }
               
               else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == H_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == C_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenC[0]->getAtomCode() == H_ATOM    ) {
                   
                        rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenB[0], 
                                                                                                        arrCurrSingleBondedOxygen[1] );
                   
                        if( isValidCarbonB == rg_TRUE ) {
                       
                            niPhosphoricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHORIC_ESTER_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[0], arrCurrSingleBondedOxygen[2] );
                            isPhosphoricChemType = rg_TRUE;
                       
                            return isPhosphoricChemType;
                       
                        }
                   
               }

               else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == C_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == H_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenC[0]->getAtomCode() == H_ATOM    ) {
                   
                        rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenA[0], 
                                                                                                        arrCurrSingleBondedOxygen[0] );
                   
                        if( isValidCarbonA == rg_TRUE ) {
                       
                            niPhosphoricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHORIC_ESTER_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[1], arrCurrSingleBondedOxygen[2] );
                            isPhosphoricChemType = rg_TRUE;
                       
                            return isPhosphoricChemType;
                        }
                   
               }

               else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == H_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == C_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenC[0]->getAtomCode() == C_ATOM    ) {
                   
                        rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenB[0], 
                                                                                                        arrCurrSingleBondedOxygen[1] );
                        rg_FLAG isValidCarbonC = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenC[0], 
                                                                                                        arrCurrSingleBondedOxygen[2] );
                   
                        if( isValidCarbonB == rg_TRUE && isValidCarbonC == rg_TRUE ) {
                       
                            niPhosphoricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHORIC_ESTER_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[0] );
                            isPhosphoricChemType = rg_TRUE;
                       
                            return isPhosphoricChemType;
                       
                        }
                   
               }

               else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == C_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == H_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenC[0]->getAtomCode() == C_ATOM    ) {
                   
                        rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenA[0], 
                                                                                                        arrCurrSingleBondedOxygen[0] );
                        rg_FLAG isValidCarbonC = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenC[0], 
                                                                                                        arrCurrSingleBondedOxygen[2] );
                   
                        if( isValidCarbonA == rg_TRUE && isValidCarbonC == rg_TRUE ) {
                       
                            niPhosphoricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHORIC_ESTER_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[1] );
                            isPhosphoricChemType = rg_TRUE;
                       
                            return isPhosphoricChemType;
                        }
               }

               else if( arrBondedAtomsOfCurrSingleBondedOxygenA[0]->getAtomCode() == C_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenB[0]->getAtomCode() == C_ATOM &&
                        arrBondedAtomsOfCurrSingleBondedOxygenC[0]->getAtomCode() == H_ATOM    ) {
                   
                        rg_FLAG isValidCarbonA = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenA[0], 
                                                                                                        arrCurrSingleBondedOxygen[0] );
                        rg_FLAG isValidCarbonB = isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( arrBondedAtomsOfCurrSingleBondedOxygenB[0], 
                                                                                                        arrCurrSingleBondedOxygen[1] );
                   
                        if( isValidCarbonA == rg_TRUE && isValidCarbonB == rg_TRUE ) {
                       
                            niPhosphoricFeature.setTypesAndAtoms( NI_FEATRURE_TYPE, PHOSPHORIC_ESTER_CHEM_TYPE, arrCurrDoubleBondedOxygen[0], arrCurrSingleBondedOxygen[2] );
                            isPhosphoricChemType = rg_TRUE;
                       
                            return isPhosphoricChemType;
                        }
               }

           }
       }
       
   }
   
   return isPhosphoricChemType;
}


void Molecule::defineHBAPharmaFeatures()
{
   Atom*  currAtom        = rg_NULL;
   rg_INT pharmaFeatureID = m_pharmaFeatures.size();
   
   m_atoms.reset4Loop();
   while( m_atoms.setNext4Loop() ) {
       
       currAtom = m_atoms.getpEntity();
       
       PharmaFeature  aHbaNlonpairFeature;
       if( isHbaNlonpair( currAtom, aHbaNlonpairFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbaNlonpairFeature);
       }

       PharmaFeature  aHbaOlonpairFeature;
       if( isHbaOlonpair( currAtom, aHbaOlonpairFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbaOlonpairFeature);
       }

       PharmaFeature  aHbaSlonpairFeature;
       if( isHbaSlonpair( currAtom, aHbaSlonpairFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbaSlonpairFeature);
       }

   }	
   
}



rg_FLAG Molecule::isHbaNlonpair( Atom* currAtom, PharmaFeature& hbaNlonpairFeature ) {

   rg_FLAG isNlonpairChemType = rg_FALSE;

   if( currAtom->getAtomCode() == N_ATOM ) {   

       rg_FLAG isCurrAtomAmine = rg_FALSE;

       rg_dList<PharmaFeature*>* pharmaFeaturesForCurrAtom = currAtom->getpChemicalProperties()->getListOfPharmaFeatures();

       pharmaFeaturesForCurrAtom->reset4Loop();
       while ( pharmaFeaturesForCurrAtom->setNext4Loop() ) {
           PharmaFeature* currPharmaFeature = pharmaFeaturesForCurrAtom->getEntity();

           if( currPharmaFeature->isChemFuncGroupType( AMINE_CHEM_TYPE ) == rg_TRUE ) {
               isCurrAtomAmine = rg_TRUE;
               break;
           }
       }

       if( isCurrAtomAmine == rg_FALSE ) {

           rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

           rg_INT  numOfSingleBond;
           rg_INT  numOfDoubleBond;
           rg_INT  numOfTripleBond;
           rg_INT  numOfAromaticBond;
           
           countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

           if( numOfBondForCurrAtom == 3 && numOfSingleBond == 3 ){

               hbaNlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, N_LONEPAIR_CHEM_TYPE, currAtom );
               isNlonpairChemType = rg_TRUE;
               
               return isNlonpairChemType;

           }

           else if( numOfBondForCurrAtom == 2 && numOfSingleBond == 1 && numOfDoubleBond == 1 ) {
               
                    hbaNlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, N_LONEPAIR_CHEM_TYPE, currAtom );
                    isNlonpairChemType = rg_TRUE;

                    return isNlonpairChemType;
           }

           else if( numOfBondForCurrAtom == 1 && numOfTripleBond == 1 ) {
               
                    hbaNlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, N_LONEPAIR_CHEM_TYPE, currAtom );
                    isNlonpairChemType = rg_TRUE;
               
                    return isNlonpairChemType;
           }

           else if( numOfBondForCurrAtom == 2 && numOfAromaticBond == 2 ) {
               
                    hbaNlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, N_LONEPAIR_CHEM_TYPE, currAtom );
                    isNlonpairChemType = rg_TRUE;
               
                    return isNlonpairChemType;
           }
       }    

   }

   return isNlonpairChemType;
}



rg_FLAG Molecule::isHbaOlonpair( Atom* currAtom, PharmaFeature& hbaOlonpairFeature ) {

   rg_FLAG isOlonpairChemType = rg_FALSE;

   if( currAtom->getAtomCode() == O_ATOM ) {
           
       rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
       
       rg_INT  numOfSingleBond;
       rg_INT  numOfDoubleBond;
       rg_INT  numOfTripleBond;
       rg_INT  numOfAromaticBond;
       
       countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
       
       if( numOfBondForCurrAtom == 2 && numOfSingleBond == 2 ){
           
           hbaOlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, O_LONEPAIR_CHEM_TYPE, currAtom );
           isOlonpairChemType = rg_TRUE;
           
           return isOlonpairChemType;  
       }

       else if( numOfBondForCurrAtom == 1 && numOfDoubleBond == 1 ){
           
                hbaOlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, O_LONEPAIR_CHEM_TYPE, currAtom );
                isOlonpairChemType = rg_TRUE;
           
                return isOlonpairChemType;  
       }  
   }
   return isOlonpairChemType;
}



rg_FLAG Molecule::isHbaSlonpair( Atom* currAtom, PharmaFeature& hbaSlonpairFeature ) {
   
   rg_FLAG isSlonpairChemType = rg_FALSE;
   
   if( currAtom->getAtomCode() == S_ATOM ) {
       
       rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

       if( numOfBondForCurrAtom < 3 ) {

           rg_INT  numOfSingleBond;
           rg_INT  numOfDoubleBond;
           rg_INT  numOfTripleBond;
           rg_INT  numOfAromaticBond;
           
           countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );

           if( numOfBondForCurrAtom == 2 && numOfSingleBond == 2 ) {
           
               hbaSlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, S_LONEPAIR_CHEM_TYPE, currAtom );
               isSlonpairChemType = rg_TRUE;
           
               return isSlonpairChemType;
           }

           else if( numOfBondForCurrAtom == 1 && numOfDoubleBond == 1 ) {
               
                    hbaSlonpairFeature.setTypesAndAtoms( HBA_FEATRURE_TYPE, S_LONEPAIR_CHEM_TYPE, currAtom );
                    isSlonpairChemType = rg_TRUE;
               
                    return isSlonpairChemType;
           }
       }
   }
   return isSlonpairChemType;
}



void Molecule::defineHBDPharmaFeatures()
{
   Atom*  currAtom        = rg_NULL;
   rg_INT pharmaFeatureID = m_pharmaFeatures.size();
   
   m_atoms.reset4Loop();
   while( m_atoms.setNext4Loop() ) {
       
       currAtom = m_atoms.getpEntity();

       PharmaFeature  aHbdHydroxylFeature;
       if( isHbdHydroxyl( currAtom, aHbdHydroxylFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbdHydroxylFeature);
       }      

       PharmaFeature  aHbdThiolFeature;
       if( isHbdThiol( currAtom, aHbdThiolFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbdThiolFeature);
       }

       PharmaFeature  aHbdAcetyleneFeature;
       if( isHbdAcetylene( currAtom, aHbdAcetyleneFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbdAcetyleneFeature);
       }

       PharmaFeature  aHbdNhFeature;
       if( isHbdNh( currAtom, aHbdNhFeature) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHbdNhFeature);
       }
       
   }
   
}



rg_FLAG Molecule::isHbdHydroxyl( Atom* currAtom, PharmaFeature& hbdHydroxylFeature ) {

   rg_FLAG isHydroxylChemType = rg_FALSE;

   rg_dList<PharmaFeature*>* listPharmaFeatureForCurrAtom = currAtom->getpChemicalProperties()->getListOfPharmaFeatures();
   
   rg_FLAG isNIFeature = rg_FALSE;
   
   listPharmaFeatureForCurrAtom->reset4Loop();
   while ( listPharmaFeatureForCurrAtom->setNext4Loop() ) {
       PharmaFeature* currPharmaFeature = listPharmaFeatureForCurrAtom->getEntity();
       if( currPharmaFeature->isChemFuncGroupType( NI_FEATRURE_TYPE ) == rg_TRUE ) {
           isNIFeature = rg_TRUE;
           break;
       }
   }
   
   if( isNIFeature == rg_FALSE ) {

       if( currAtom->getAtomCode() == O_ATOM ) {

           rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
           rg_FLAG isAllSingleBond      = isAllBondsSingle( currAtom );

           if( numOfBondForCurrAtom == 2 && isAllSingleBond == rg_TRUE ) {

               rg_INT  numOfCarbon;
               rg_INT  numOfHydrogen;
               rg_INT  numOfNitrogen;
               rg_INT  numOfOxygen;
               rg_INT  numOfPhosphorus;
               rg_INT  numOfSulfur;
               
               countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );

               if( numOfHydrogen == 1 ) {  
                   
                   hbdHydroxylFeature.setTypesAndAtoms( HBD_FEATRURE_TYPE, HYDROXYL_CHEM_TYPE, currAtom );
                   isHydroxylChemType = rg_TRUE;
                   
                   return isHydroxylChemType;
               }
        }
       }
   }

   else {
           isHydroxylChemType = rg_FALSE;
           return isHydroxylChemType;
   }

   return isHydroxylChemType;
}



rg_FLAG Molecule::isHbdThiol( Atom* currAtom, PharmaFeature& hbdThiolFeature ) {
   
   rg_FLAG isThiolChemType = rg_FALSE;
   
   if( currAtom->getAtomCode() == S_ATOM ) {
       
       rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
       rg_FLAG isAllSingleBond  = isAllBondsSingle( currAtom );
       
       if( numOfBondForCurrAtom == 2 && isAllSingleBond == rg_TRUE ) {
           
           hbdThiolFeature.setTypesAndAtoms( HBD_FEATRURE_TYPE, THIOL_CHEM_TYPE, currAtom );
           isThiolChemType = rg_TRUE;
           
           return isThiolChemType;
       }
   }
   return isThiolChemType;
}



rg_FLAG Molecule::isHbdAcetylene( Atom* currAtom, PharmaFeature& hbdAcetyleneFeature ) {
   
   rg_FLAG isAcetyleneChemType = rg_FALSE;
   
   if( currAtom->getAtomCode() == C_ATOM ) {
       
       rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();

       rg_INT  numOfSingleBond;
       rg_INT  numOfDoubleBond;
       rg_INT  numOfTripleBond;
       rg_INT  numOfAromaticBond;
       
       countEachBondTypeFromBondListInAtom( currAtom, numOfSingleBond, numOfDoubleBond, numOfTripleBond, numOfAromaticBond );
       
       if( numOfBondForCurrAtom == 2 && numOfSingleBond == 1 && numOfTripleBond == 1 ) {

           rg_INT  numOfCarbon;
           rg_INT  numOfHydrogen;
           rg_INT  numOfNitrogen;
           rg_INT  numOfOxygen;
           rg_INT  numOfPhosphorus;
           rg_INT  numOfSulfur;
           
           countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbon, numOfHydrogen, numOfNitrogen, numOfOxygen, numOfPhosphorus, numOfSulfur );

           if( numOfCarbon == 1 && numOfHydrogen == 1 ) {  
           
               hbdAcetyleneFeature.setTypesAndAtoms( HBD_FEATRURE_TYPE, ACETYLENE_CHEM_TYPE, currAtom );
               isAcetyleneChemType = rg_TRUE;
           
               return isAcetyleneChemType;
           }
       }
   }
   return isAcetyleneChemType;
}



rg_FLAG Molecule::isHbdNh( Atom* currAtom, PharmaFeature& hbdNhFeature ) {
   
   rg_FLAG isNhChemType = rg_FALSE;
   
   rg_dList<PharmaFeature*>* listPharmaFeatureForCurrAtom = currAtom->getpChemicalProperties()->getListOfPharmaFeatures();
   
   rg_FLAG isNIFeature = rg_FALSE;
   
   listPharmaFeatureForCurrAtom->reset4Loop();
   while ( listPharmaFeatureForCurrAtom->setNext4Loop() ) {
       PharmaFeature* currPharmaFeature = listPharmaFeatureForCurrAtom->getEntity();
       if( currPharmaFeature->isChemFuncGroupType( NI_FEATRURE_TYPE ) == rg_TRUE ) {
           isNIFeature = rg_TRUE;
           break;
       }
   }
   
   if( isNIFeature == rg_FALSE ) {
       
       if( currAtom->getAtomCode() == N_ATOM ) {

           Atom* bondedAtom = currAtom->getListChemicalBond()->getFirstEntity()->getBondedAtom( currAtom );

           rg_INT  numOfCarbonForBondedAtom;
           rg_INT  numOfHydrogenForBondedAtom;
           rg_INT  numOfNitrogenForBondedAtom;
           rg_INT  numOfOxygenForBondedAtom;
           rg_INT  numOfPhosphorusForBondedAtom;
           rg_INT  numOfSulfurForBondedAtom;
           
           countEachAtomTypeFromBondListInAtom( currAtom, numOfCarbonForBondedAtom, numOfHydrogenForBondedAtom, numOfNitrogenForBondedAtom, 
                                                numOfOxygenForBondedAtom, numOfPhosphorusForBondedAtom, numOfSulfurForBondedAtom );

           if( numOfHydrogenForBondedAtom >= 1 ) {
               hbdNhFeature.setTypesAndAtoms( HBD_FEATRURE_TYPE, NH_CHEM_TYPE, currAtom );
               isNhChemType = rg_TRUE;
               
               return isNhChemType;
           }
       }
   }
   
   else { 

           isNhChemType = rg_FALSE;
           return isNhChemType;
   }
   
   return isNhChemType;
}



void Molecule::defineHYPharmaFeatures()
{
   Atom*  currAtom        = rg_NULL;
   rg_INT pharmaFeatureID = m_pharmaFeatures.size();
   
   m_atoms.reset4Loop();
   while( m_atoms.setNext4Loop() ) {
       
       currAtom = m_atoms.getpEntity();

        list<PharmaFeature> listHyCRingFeature;
        if( isHyCRing( currAtom, listHyCRingFeature ) == rg_TRUE ) {
            addPharmaFeature( pharmaFeatureID, &listHyCRingFeature);
        }

        PharmaFeature  aHyHalogenFeature;
        if( isHyHalogen( currAtom, aHyHalogenFeature ) == rg_TRUE ) {
            addPharmaFeature( pharmaFeatureID, &aHyHalogenFeature);
        }
   }

   m_atoms.reset4Loop();
   while( m_atoms.setNext4Loop() ) {
       currAtom = m_atoms.getpEntity();

       PharmaFeature  aHyCGroupFeature;
       if( isHyCGroup( currAtom, aHyCGroupFeature ) == rg_TRUE ) {
           addPharmaFeature( pharmaFeatureID, &aHyCGroupFeature);
       }
   }

   
}



rg_FLAG Molecule::isHyCRing( Atom* currAtom, list<PharmaFeature>& listHyCRingFeature )
{
   rg_FLAG isCRingChemType = rg_FALSE;

   list<Atom**> listOfAtomArrayInPentaRing;
   list<Atom**> listOfAtomArrayInHexaRing;

   rg_FLAG isOnPentaRing = isOnRing( currAtom, 5, listOfAtomArrayInPentaRing );
   rg_FLAG isOnHexaRing  = isOnRing( currAtom, 6, listOfAtomArrayInHexaRing );

   if( isOnPentaRing == rg_TRUE ) {

       list<Atom**>::iterator i_listOfAtomArrayInPentaRing = listOfAtomArrayInPentaRing.begin();

       while ( i_listOfAtomArrayInPentaRing != listOfAtomArrayInPentaRing.end() ) {
   
           Atom** arrAtomsOnPentaRing = (Atom**)(*i_listOfAtomArrayInPentaRing);

           // Check all atoms are carbon
           rg_FLAG isAllAtomCarbon = rg_TRUE;
           for( rg_INT i_atomInPentaRing=0; i_atomInPentaRing<5; i_atomInPentaRing++ ) {
               if( arrAtomsOnPentaRing[i_atomInPentaRing]->getAtomCode() != C_ATOM ) {
                   isAllAtomCarbon = rg_FALSE;
                   break;
               }
           }
           
           if( isAllAtomCarbon == rg_TRUE ) {
               PharmaFeature hyCRingFeature;
               hyCRingFeature.setTypesAndAtoms( HY_FEATRURE_TYPE, C_RING_CHEM_TYPE, arrAtomsOnPentaRing[0], arrAtomsOnPentaRing[1], arrAtomsOnPentaRing[2], arrAtomsOnPentaRing[3], arrAtomsOnPentaRing[4] );
               listHyCRingFeature.push_back( hyCRingFeature );
               isCRingChemType = rg_TRUE;
           }
           
           i_listOfAtomArrayInPentaRing = listOfAtomArrayInPentaRing.erase( i_listOfAtomArrayInPentaRing );
           delete [] arrAtomsOnPentaRing;
       }
   }

   if( isOnHexaRing == rg_TRUE ) {

       list<Atom**>::iterator i_listOfAtomArrayInHexaRing = listOfAtomArrayInHexaRing.begin();

       while ( i_listOfAtomArrayInHexaRing != listOfAtomArrayInHexaRing.end() ) {
   
           Atom** arrAtomsOnHexaRing = (Atom**)(*i_listOfAtomArrayInHexaRing);

           // Check all atoms are carbon
           rg_FLAG isAllAtomCarbon = rg_TRUE;
           for( rg_INT i_atomInHexaRing=0; i_atomInHexaRing<6; i_atomInHexaRing++ ) {
               if( arrAtomsOnHexaRing[i_atomInHexaRing]->getAtomCode() != C_ATOM ) {
                   isAllAtomCarbon = rg_FALSE;
                   break;
               }
           }
           
           if( isAllAtomCarbon == rg_TRUE ) {
               PharmaFeature hyCRingFeature;
               hyCRingFeature.setTypesAndAtoms( HY_FEATRURE_TYPE, C_RING_CHEM_TYPE, arrAtomsOnHexaRing[0], arrAtomsOnHexaRing[1], arrAtomsOnHexaRing[2], arrAtomsOnHexaRing[3], arrAtomsOnHexaRing[4], arrAtomsOnHexaRing[5] );
               listHyCRingFeature.push_back( hyCRingFeature );
               isCRingChemType = rg_TRUE;
           }

           i_listOfAtomArrayInHexaRing = listOfAtomArrayInHexaRing.erase( i_listOfAtomArrayInHexaRing );
           delete [] arrAtomsOnHexaRing;
       }
   }

   return isCRingChemType;
}



rg_FLAG Molecule::isHyHalogen( Atom* currAtom, PharmaFeature& hyHalogenFeature ) {
   
   rg_FLAG isHalogenChemType = rg_FALSE;
  
   if( currAtom->getAtomCode() ==  F_ATOM || currAtom->getAtomCode() == CL_ATOM ||
       currAtom->getAtomCode() == BR_ATOM || currAtom->getAtomCode() ==  I_ATOM  ) {

       rg_INT  numOfBondForCurrAtom = currAtom->getListChemicalBond()->getSize();
       rg_FLAG isAllSingleBond      = isAllBondsSingle( currAtom );
       
       if( numOfBondForCurrAtom == 1 && isAllSingleBond == rg_TRUE ) {

           Atom* bondedAtom = currAtom->getListChemicalBond()->getFirstEntity()->getBondedAtom( currAtom );

           rg_INT  numOfSingleBondForBondedAtom;
           rg_INT  numOfDoubleBondForBondedAtom;
           rg_INT  numOfTripleBondForBondedAtom;
           rg_INT  numOfAromaticBondForBondedAtom;
       
           countEachBondTypeFromBondListInAtom( bondedAtom, numOfSingleBondForBondedAtom, 
                                                            numOfDoubleBondForBondedAtom, 
                                                            numOfTripleBondForBondedAtom, 
                                                            numOfAromaticBondForBondedAtom );

           if( numOfAromaticBondForBondedAtom > 0 ) {
               isHalogenChemType = rg_FALSE;
           }

           else { 
               hyHalogenFeature.setTypesAndAtoms( HY_FEATRURE_TYPE, HALOGEN_CHEM_TYPE, currAtom );
               isHalogenChemType = rg_TRUE;
           }
       }

       else {
           isHalogenChemType = rg_FALSE;
       }
   }
       
   return isHalogenChemType;
}



rg_FLAG Molecule::isHyCGroup( Atom* currAtom, PharmaFeature& hyCGroupFeature )
{
   if( currAtom->getAtomCode() != C_ATOM ) {
       return rg_FALSE;
   }

   rg_FLAG isCGroupChemType = rg_FALSE;

   typedef pair<Atom*, Atom*> AtomPair;    // first: prevAtom, second: currAtom

   list<Atom*> atomsOfFeatureCandidate;

   typedef list<AtomPair> AtomQueue;
   AtomQueue queueOfAtom;


   // Check for seed atom(currAtom) first.
   queueOfAtom.push_back( AtomPair( (Atom*)NULL, currAtom ) );
   
   AtomPair popedAtom = queueOfAtom.front();
   queueOfAtom.pop_front();

   rg_FLAG isCRingChemType   = isHyCRingChemType( popedAtom.second );
   rg_FLAG isBondedWithCAndH = isBondedWithCarbonAndHydrogen( popedAtom.second );

   if( isCRingChemType == rg_FALSE && isBondedWithCAndH == rg_TRUE ) {
       atomsOfFeatureCandidate.push_back( popedAtom.second );
       
       rg_INT numOfBondedCarbon = getNumOfBondedAtomTypeFromCurrAtom( popedAtom.second, C_ATOM );      

       if( numOfBondedCarbon == 0 )
           return rg_FALSE;

       Atom** arrBondedCarbon   = new Atom* [numOfBondedCarbon];
       getBondedAtomFromCurrAtom( popedAtom.second, C_ATOM, arrBondedCarbon );
   
       for( rg_INT i_bondedCarbon=0; i_bondedCarbon<numOfBondedCarbon; i_bondedCarbon++ ) {
           queueOfAtom.push_back( AtomPair( popedAtom.second, arrBondedCarbon[i_bondedCarbon] ) );
       }

       delete [] arrBondedCarbon;
   }

   // Check for rest of bonded atoms of currAtom
   while( queueOfAtom.size() != 0 ) {

       rg_INT sizeOfQueue = queueOfAtom.size();
       popedAtom = queueOfAtom.front();
       queueOfAtom.pop_front();
       sizeOfQueue = queueOfAtom.size();

       isCRingChemType   = isHyCRingChemType( popedAtom.second );
       isBondedWithCAndH = isBondedWithCarbonAndHydrogen( popedAtom.second );

       if( isCRingChemType == rg_FALSE && isBondedWithCAndH == rg_TRUE ) {
           
           list<Atom*>::iterator i_atomOfFeatureCandidate = find( atomsOfFeatureCandidate.begin(), atomsOfFeatureCandidate.end(), popedAtom.second );

           if( i_atomOfFeatureCandidate == atomsOfFeatureCandidate.end() ) {
               atomsOfFeatureCandidate.push_back( popedAtom.second );
           }
           
       
           rg_INT numOfBondedCarbon = getNumOfBondedAtomTypeFromCurrAtom( popedAtom.second, C_ATOM );

           if( numOfBondedCarbon == 1 )
               continue;

           Atom** arrBondedCarbon   = new Atom* [numOfBondedCarbon-1];
           getBondedAtomFromCurrAtomExceptPrevAtom( popedAtom.second, popedAtom.first, C_ATOM, arrBondedCarbon );
   
           for( rg_INT i_bondedCarbon=0; i_bondedCarbon<numOfBondedCarbon-1; i_bondedCarbon++ ) {

               i_atomOfFeatureCandidate = find( atomsOfFeatureCandidate.begin(), atomsOfFeatureCandidate.end(), arrBondedCarbon[i_bondedCarbon] );

               if( i_atomOfFeatureCandidate == atomsOfFeatureCandidate.end() ) {
                   queueOfAtom.push_back( AtomPair( popedAtom.second, arrBondedCarbon[i_bondedCarbon] ) );                    
               }
           }

           delete [] arrBondedCarbon;
       }
   }

   if( atomsOfFeatureCandidate.size() > 1 ) {
       isCGroupChemType = rg_TRUE;
       hyCGroupFeature.setTypesAndAtoms( HY_FEATRURE_TYPE, C_GROUP_CHEM_TYPE, &atomsOfFeatureCandidate );
   }

   return isCGroupChemType;
}



rg_FLAG Molecule::isHyCRingChemType( Atom* currAtom ) 
{
   rg_dList<PharmaFeature*>* pharmaFeatureOfCurrAtom = currAtom->getpChemicalProperties()->getListOfPharmaFeatures();

   rg_FLAG isCRingChemType = rg_FALSE;
   pharmaFeatureOfCurrAtom->reset4Loop();
   while( pharmaFeatureOfCurrAtom->setNext4Loop() ) {
       PharmaFeature* currPharmaFeature = pharmaFeatureOfCurrAtom->getEntity();

       if( currPharmaFeature->isChemFuncGroupType( C_RING_CHEM_TYPE ) == rg_TRUE ) {
           isCRingChemType = rg_TRUE;
           break;
       }
   }

   return isCRingChemType;
}



rg_FLAG Molecule::isBondedWithCarbonAndHydrogen( Atom* currAtom )
{
   rg_FLAG isCAndH = rg_TRUE;
   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  
   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       AtomCode bondedAtomCode = currBondList->getEntity()->getBondedAtom( currAtom )->getAtomCode();

       if( bondedAtomCode == C_ATOM  || bondedAtomCode == H_ATOM ) {
           isCAndH = rg_TRUE;
       }
       else {
           isCAndH = rg_FALSE;
           break;
       }
   }

   return isCAndH;
}



rg_FLAG Molecule::isAllBondedAtomsNonNegative( Atom* currAtom )
{
   rg_FLAG isNonNegative = rg_FALSE;
   
   ChemicalBond* currBond = rg_NULL;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {
       currBond = currBondList->getEntity();

       Atom* bondedAtom = currBond->getBondedAtom( currAtom );

       if( bondedAtom->getChemicalProperties().getCharge() > -1.0) {
           isNonNegative = rg_TRUE;
       }
       else {
           isNonNegative = rg_FALSE;
           break;
       }
   }

   return isNonNegative;
}



rg_FLAG Molecule::isAllBondedAtomsNonPositive( Atom* currAtom )
{
   rg_FLAG isNonPositive = rg_FALSE;
   
   ChemicalBond* currBond = rg_NULL;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {
       currBond = currBondList->getEntity();

       Atom* bondedAtom = currBond->getBondedAtom( currAtom );

       if( bondedAtom->getChemicalProperties().getCharge() < 1.0 ) {
           isNonPositive = rg_TRUE;
       }
       else {
           isNonPositive = rg_FALSE;
           break;
       }
   }

   return isNonPositive;
}



rg_FLAG Molecule::isAllBondsSingle( Atom* currAtom )
{
   rg_FLAG isSingleBond = rg_FALSE;
   
   ChemicalBond* currBond = rg_NULL;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {
       currBond = currBondList->getEntity();

       if( currBond->getTypeOfBond() == SINGLE_BOND || currBond->getTypeOfBond() == AMIDE_BOND ) {
           isSingleBond = rg_TRUE;
       }
       else {
           isSingleBond = rg_FALSE;
           break;
       }
   }

   return isSingleBond;
}



void Molecule::countEachBondTypeFromBondListInAtom( Atom* currAtom, rg_INT& numOfSingleBond, rg_INT& numOfDoubleBond, rg_INT& numOfTripleBond, rg_INT& numOfAromaticBond )
{
   ChemicalBond* currBond = rg_NULL;

   numOfSingleBond   = 0;
   numOfDoubleBond   = 0;
   numOfTripleBond   = 0;
   numOfAromaticBond = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {
       currBond = currBondList->getEntity();
       if( currBond->getTypeOfBond() == SINGLE_BOND || currBond->getTypeOfBond() == AMIDE_BOND ) {
           numOfSingleBond++;
       }
       if( currBond->getTypeOfBond() == DOUBLE_BOND ) {
           numOfDoubleBond++;
       }
       if( currBond->getTypeOfBond() == TRIPLE_BOND ) {
           numOfTripleBond++;
       }
       if( currBond->getTypeOfBond() == AROMATIC_BOND           || 
           currBond->getTypeOfBond() == SINGLE_OR_AROMATIC_BOND || 
           currBond->getTypeOfBond() == DOUBLE_OR_AROMATIC_BOND   ) {
           numOfAromaticBond++;
       }
   }
}



void Molecule::countEachAtomTypeFromBondListInAtom( Atom* currAtom, rg_INT& numOfCarbon, rg_INT& numOfHydrogen, rg_INT& numOfNitrogen, rg_INT& numOfOxygen, rg_INT& numOfPhosphorus, rg_INT& numOfSulfur )
{
   numOfCarbon     = 0;
   numOfHydrogen   = 0;
   numOfNitrogen   = 0;
   numOfOxygen     = 0;
   numOfPhosphorus = 0;
   numOfSulfur     = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       AtomCode bondedAtomCode = currBondList->getEntity()->getBondedAtom( currAtom )->getAtomCode();

       if( bondedAtomCode == C_ATOM ) {
           numOfCarbon++;
       }
       else if( bondedAtomCode == H_ATOM ) {
           numOfHydrogen++;
       }
       else if( bondedAtomCode == N_ATOM ) {
           numOfNitrogen++;
       }
       else if( bondedAtomCode == O_ATOM ) {
           numOfOxygen++;
       }
       else if( bondedAtomCode == P_ATOM ) {
           numOfPhosphorus++;
       }
       else if( bondedAtomCode == S_ATOM ) {
           numOfSulfur++;
       }
   }
}



void Molecule::getEachpAtomFromBondListInAtom( Atom* currAtom, Atom** arrCarbon, Atom** arrHydrogen, Atom** arrNitrogen, Atom** arrOxygen, Atom** arrPhosphorus, Atom** arrSulfur )
{
   rg_INT i_Carbon     = 0;
   rg_INT i_Hydrogen   = 0;
   rg_INT i_Nitrogen   = 0;
   rg_INT i_Oxygen     = 0;
   rg_INT i_Phosphorus = 0;
   rg_INT i_Sulfur     = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       Atom*    pBondedAtom  = currBondList->getEntity()->getBondedAtom( currAtom );
       AtomCode bondedAtomCode = pBondedAtom->getAtomCode();

       if( bondedAtomCode == C_ATOM ) {
           arrCarbon[i_Carbon] = pBondedAtom;
           i_Carbon++;
       }
       else if( bondedAtomCode == H_ATOM ) {
           arrHydrogen[i_Hydrogen] = pBondedAtom;
           i_Hydrogen++;
       }
       else if( bondedAtomCode == N_ATOM ) {
           arrNitrogen[i_Nitrogen] = pBondedAtom;
           i_Nitrogen++;
       }
       else if( bondedAtomCode == O_ATOM ) {
           arrOxygen[i_Oxygen] = pBondedAtom;
           i_Oxygen++;
       }
       else if( bondedAtomCode == P_ATOM ) {
           arrPhosphorus[i_Phosphorus] = pBondedAtom;
           i_Phosphorus++;
       }
       else if( bondedAtomCode == S_ATOM ) {
           arrSulfur[i_Sulfur] = pBondedAtom;
           i_Sulfur++;
       }
   }    
}



rg_INT Molecule::getNumOfBondedAtomTypeFromCurrAtom( Atom* currAtom, const AtomCode& targetAtomCode )
{
   rg_INT numOfTargetAtom = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {
       AtomCode bondedAtomCode = currBondList->getEntity()->getBondedAtom( currAtom )->getAtomCode();
  
       if( bondedAtomCode == targetAtomCode ) {
           numOfTargetAtom++;
       }
   }

   return numOfTargetAtom;
}



void Molecule::getBondedAtomFromCurrAtom( Atom* currAtom, const AtomCode& targetAtomCode, const BondType& targetBondType, Atom** arrBondedAtom )
{
   rg_INT i_bondedAtom = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       ChemicalBond* pCurrBond      = currBondList->getEntity();
       Atom*         pBondedAtom    = pCurrBond->getBondedAtom( currAtom );
       AtomCode      bondedAtomCode = pBondedAtom->getAtomCode();

       if( bondedAtomCode == targetAtomCode && pCurrBond->getTypeOfBond() == targetBondType ) {
           arrBondedAtom[i_bondedAtom] = pBondedAtom;
           i_bondedAtom++;
       }
   }    
}



void Molecule::getBondedAtomFromCurrAtom( Atom* currAtom, const AtomCode& targetAtomCode, Atom** arrTargetAtom )
{
   rg_INT i_targetAtom = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       Atom* pBondedAtom = currBondList->getEntity()->getBondedAtom( currAtom );
  
       if( pBondedAtom->getAtomCode() == targetAtomCode ) {
           arrTargetAtom[i_targetAtom] = pBondedAtom;
           i_targetAtom++;
       }
   }
}



void Molecule::getBondedAtomFromCurrAtomExceptPrevAtom( Atom* currAtom, Atom* prevAtom, Atom** arrTargetAtom )
{
   rg_INT i_targetAtom = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       Atom* pBondedAtom    = currBondList->getEntity()->getBondedAtom( currAtom );
       
       if ( pBondedAtom == prevAtom )
           continue;
   
       arrTargetAtom[i_targetAtom] = pBondedAtom;
       i_targetAtom++;
   }
}



void Molecule::getBondedAtomFromCurrAtomExceptPrevAtom( Atom* currAtom, Atom* prevAtom, const AtomCode& targetAtomCode, Atom** arrTargetAtom )
{
   rg_INT i_targetAtom = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       Atom* pBondedAtom = currBondList->getEntity()->getBondedAtom( currAtom );
  
       if ( pBondedAtom == prevAtom )
           continue;

       if( pBondedAtom->getAtomCode() == targetAtomCode ) {
           arrTargetAtom[i_targetAtom] = pBondedAtom;
           i_targetAtom++;
       }
   }    
}



rg_FLAG Molecule::isBondedAtomsExceptPrevAtomConsistOfCarbonAndHydrogen( Atom* currAtom, Atom* prevAtom )
{
   rg_FLAG isValid = rg_TRUE;
   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {
       ChemicalBond* pCurrBond      = currBondList->getEntity();
       Atom*         pBondedAtom    = pCurrBond->getBondedAtom( currAtom );
       
       if( pBondedAtom == prevAtom )
           continue;

       if( pBondedAtom->getAtomCode() == C_ATOM || pBondedAtom->getAtomCode() == H_ATOM) {
           isValid = rg_TRUE;
       }
       else {
           isValid = rg_FALSE;
           break;
       }
   }

   return isValid;    
}



void Molecule::addPharmaFeature( rg_INT& featureID, PharmaFeature* aPharmaFeature )
{
   PharmaFeature currPharmaFeature = *aPharmaFeature;

   currPharmaFeature.setID( featureID );
   featureID++;
   
   currPharmaFeature.setCentersOfPharmaFeature();
   m_pharmaFeatures.push_back( currPharmaFeature );

   list<PharmaFeature>::iterator i_currPharmaFeature = m_pharmaFeatures.end();
   i_currPharmaFeature--;

   rg_dList<Atom*>* pAtoms = aPharmaFeature->getAtoms();

   pAtoms->reset4Loop();
   while( pAtoms->setNext4Loop() ) {
       Atom* currAtomInFeature = pAtoms->getEntity();
       currAtomInFeature->getpChemicalProperties()->addPharmaFeature( &(*i_currPharmaFeature) );
   }
}



void Molecule::addPharmaFeature( rg_INT& featureID, list<PharmaFeature>* aListPharmaFeature )
{
   list<PharmaFeature>::iterator i_listPharmaFeature = aListPharmaFeature->begin();

   while ( i_listPharmaFeature != aListPharmaFeature->end() ) {
       
       addPharmaFeature( featureID, &(*i_listPharmaFeature) );
   
       i_listPharmaFeature++;
   }
}



rg_FLAG Molecule::isOnRing( Atom* targetAtom, const rg_INT& numOfAtomsInRing, list<Atom**>& listOfAtomArrayInRing )
{
   rg_INT i_atomIDInRing = 0;

   list<ChemicalBond**> listOfRingCandidates;
   	
   // targetAtom \BF\A1 \B4\EB\C7Ͽ\A9 \BF\AC\B0\E1\B5\C8 bond\B5\E9 \C0\FA\C0\E5.
   rg_dList<ChemicalBond*>* firstBonds = targetAtom->getListChemicalBond();
   firstBonds->reset4Loop();
	while( firstBonds->setNext4Loop() ) {
       ChemicalBond** ringCandidate = new ChemicalBond* [numOfAtomsInRing];

       for( i_atomIDInRing=0; i_atomIDInRing<numOfAtomsInRing; i_atomIDInRing++ ) {
           ringCandidate[i_atomIDInRing] = rg_NULL;
       }

       ringCandidate[0] = firstBonds->getEntity();
       listOfRingCandidates.push_back( ringCandidate );

       Atom** firstAtomInRingCand = new Atom* [numOfAtomsInRing];
       firstAtomInRingCand[0] = targetAtom;
       listOfAtomArrayInRing.push_back( firstAtomInRingCand );
	}
   
   for( i_atomIDInRing=1; i_atomIDInRing<numOfAtomsInRing; i_atomIDInRing++ ) {

       list<ChemicalBond**>::iterator i_listOfRingCandidates       = listOfRingCandidates.begin();
       list<Atom**>::iterator         i_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.begin();

       while ( i_listOfRingCandidates != listOfRingCandidates.end() ) {
          
           ChemicalBond** arrBondsOfRingCandidate = (ChemicalBond**)(*i_listOfRingCandidates);
           Atom**         arrFirstAtom            = (Atom**)(*i_listOfFirstAtomsInRingCand);
           
           Atom*  secondAtomInBond        = arrBondsOfRingCandidate[i_atomIDInRing-1]->getBondedAtom( arrFirstAtom[i_atomIDInRing-1] );          
           rg_INT numOfBondsForSecondAtom = secondAtomInBond->getListChemicalBond()->getSize();

           // Leaf \C0̸\E9, \C7ش\E7 ring candidate\B8\A6 \B8\AE\BD\BAƮ\BF\A1\BC\AD \BB\E8\C1\A6...
           if( numOfBondsForSecondAtom == 1 ) {
               
                i_listOfRingCandidates       = listOfRingCandidates.erase( i_listOfRingCandidates );
                i_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.erase( i_listOfFirstAtomsInRingCand );
                delete [] arrBondsOfRingCandidate;
                delete [] arrFirstAtom;

               continue;
           }

           ChemicalBond** bondsOfSecAtomExceptCurrBonds = new ChemicalBond* [numOfBondsForSecondAtom-1];

           getpBondsFromBondListInAtomExceptPrevBond( secondAtomInBond, arrBondsOfRingCandidate[i_atomIDInRing-1], bondsOfSecAtomExceptCurrBonds );
          
           // bondsOfSecAtomExceptCurrBonds \BF\A1 \C0\FA\C0\E5\B5\C8 2\B9\F8° bond \BA\CE\C5\CD \B0\E6\B7\CE \B0\BB\BD\C5 \BD\C3\C0\DB ( ù\B9\F8°\B4\C2 \C0ϴ\DC \B3\B2\B0ܳ\F6\BE\DF \C7Ѵ\D9. : \BF\F8\BA\BB \B0\E6\B7\CE \BA\B9\BB\E7\BF\EB )
           for( rg_INT i_bondsOfSecAtom=1; i_bondsOfSecAtom<numOfBondsForSecondAtom-1; i_bondsOfSecAtom++ ) {
               
               if ( isValidBondForRing( i_atomIDInRing, arrBondsOfRingCandidate, bondsOfSecAtomExceptCurrBonds[i_bondsOfSecAtom] ) == rg_TRUE ) {

                   // \BB\F5\B7ο\EE RingArray\B8\A6 \B8\B8\B5\E9\BE�?\BA\B9\BB\E7 \C8\C4 \B8\AE\BD\BAƮ\BF\A1 \C0\FA\C0\E5.
                   ChemicalBond** newRingCandidate       = new ChemicalBond* [numOfAtomsInRing];
                   Atom**         newFirstAtomInRingCand = new Atom* [numOfAtomsInRing];

                   for( rg_INT j_atomIdInring=0; j_atomIdInring<numOfAtomsInRing; j_atomIdInring++ ) {
                       newRingCandidate[j_atomIdInring]       = arrBondsOfRingCandidate[j_atomIdInring];
                       newFirstAtomInRingCand[j_atomIdInring] = arrFirstAtom[j_atomIdInring];
                   }

                   newRingCandidate[i_atomIDInRing]       = bondsOfSecAtomExceptCurrBonds[i_bondsOfSecAtom];
                   newFirstAtomInRingCand[i_atomIDInRing] = secondAtomInBond;
                   listOfRingCandidates.push_front( newRingCandidate );
                   listOfAtomArrayInRing.push_front( newFirstAtomInRingCand );
               }

           }

           // bondsOfSecAtomExceptCurrBonds \BF\A1 \C0\FA\C0\E5\B5\C8 1\B9\F8° bond \B0\E6\B7\CE \B0\BB\BD\C5 \BD\C3\C0\DB
           if ( isValidBondForRing( i_atomIDInRing, arrBondsOfRingCandidate, bondsOfSecAtomExceptCurrBonds[0] ) == rg_TRUE ) {
               arrBondsOfRingCandidate[i_atomIDInRing] = bondsOfSecAtomExceptCurrBonds[0]; // ù\B9\F8° bond\B4\C2 \C1\EF\BD\C3 \B0\E6\B7ο\A1 \C6\F7\C7Խ\C3\C4\D1\C1ش\D9.
               arrFirstAtom[i_atomIDInRing] = secondAtomInBond; // bond\C0\C7 firstAtom \B3־\EE\C1\DC.. 

               i_listOfRingCandidates++;
               i_listOfFirstAtomsInRingCand++;
           }
           else {
               i_listOfRingCandidates       = listOfRingCandidates.erase( i_listOfRingCandidates );
               i_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.erase( i_listOfFirstAtomsInRingCand );
               delete [] arrBondsOfRingCandidate;
               delete [] arrFirstAtom;

           }

           delete [] bondsOfSecAtomExceptCurrBonds;
       }

   }

   
   // Evaluate candidates of Ring
   list<ChemicalBond**>::iterator i_listOfRingCandidates       = listOfRingCandidates.begin();
   list<Atom**>::iterator         i_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.begin();

   while ( i_listOfRingCandidates != listOfRingCandidates.end() ) {
   
       ChemicalBond** arrBondsOfRingCandidate = (ChemicalBond**)(*i_listOfRingCandidates);
       Atom**         arrFirstAtom            = (Atom**)(*i_listOfFirstAtomsInRingCand);

       Atom* startAtom = arrFirstAtom[0];
       Atom* endAtom   = arrBondsOfRingCandidate[numOfAtomsInRing-1]->getBondedAtom( arrFirstAtom[numOfAtomsInRing-1] );

       if( startAtom != endAtom ) {
           i_listOfRingCandidates       = listOfRingCandidates.erase( i_listOfRingCandidates );
           i_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.erase( i_listOfFirstAtomsInRingCand );
           delete [] arrBondsOfRingCandidate;
           delete [] arrFirstAtom;
       }
       else {
           i_listOfRingCandidates++;
           i_listOfFirstAtomsInRingCand++;
       }    
   }

   // Delete duplicated Ring
   ChemicalBond** arrBondsOfRingCandidateA = rg_NULL;
   ChemicalBond** arrBondsOfRingCandidateB = rg_NULL;
   Atom**         arrFirstAtomB            = rg_NULL;

   for( i_listOfRingCandidates=listOfRingCandidates.begin(); i_listOfRingCandidates!=listOfRingCandidates.end(); i_listOfRingCandidates++) {
		arrBondsOfRingCandidateA = (ChemicalBond**)(*i_listOfRingCandidates);
       
       list<ChemicalBond**>::iterator j_listOfRingCandidates       = listOfRingCandidates.begin();
       list<Atom**>::iterator         j_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.begin();

       while( j_listOfRingCandidates != listOfRingCandidates.end() ) {

           if( j_listOfRingCandidates == i_listOfRingCandidates ) {
               j_listOfRingCandidates++;
               j_listOfFirstAtomsInRingCand++;
               continue;
           }

		    arrBondsOfRingCandidateB = (ChemicalBond**)(*j_listOfRingCandidates);
           arrFirstAtomB            = (Atom**)(*j_listOfFirstAtomsInRingCand);

           if( isIdenticalRing( numOfAtomsInRing, arrBondsOfRingCandidateA, arrBondsOfRingCandidateB ) == rg_TRUE ) {
               j_listOfRingCandidates       = listOfRingCandidates.erase( j_listOfRingCandidates );
               j_listOfFirstAtomsInRingCand = listOfAtomArrayInRing.erase( j_listOfFirstAtomsInRingCand );
               delete [] arrBondsOfRingCandidateB;
               delete [] arrFirstAtomB;
               break;
           }
           else {
               j_listOfRingCandidates++;
               j_listOfFirstAtomsInRingCand++;
           }
       }
   }


   // delete dynamic array stored in listOfRingCandidates
   i_listOfRingCandidates       = listOfRingCandidates.begin();
   while ( i_listOfRingCandidates != listOfRingCandidates.end() ) {
   
       ChemicalBond** arrBondsOfRingCandidate = (ChemicalBond**)(*i_listOfRingCandidates);
       delete [] arrBondsOfRingCandidate;
       i_listOfRingCandidates++;
   }


   if( listOfAtomArrayInRing.size() == 0 )
       return rg_FALSE;
   else
	    return rg_TRUE;
}



void Molecule::getpBondsFromBondListInAtomExceptPrevBond( Atom* currAtom, ChemicalBond* prevBond, ChemicalBond** arrBonds )
{
   rg_INT i_bond = 0;

   rg_dList<ChemicalBond*>* currBondList = currAtom->getListChemicalBond();  

   currBondList->reset4Loop();
   while( currBondList->setNext4Loop() ) {

       ChemicalBond* pCurrBond = currBondList->getEntity();
       
       if ( pCurrBond == prevBond )
           continue;
   
       arrBonds[i_bond] = pCurrBond;
       i_bond++;
   }
}



rg_FLAG Molecule::isValidBondForRing( const rg_INT& atomIDOfRing, ChemicalBond** arrBondsOfRingCandidate, ChemicalBond* targetBond )
{
   rg_FLAG isValidBond = rg_TRUE;

   for( rg_INT i_atomIDOfRing=0; i_atomIDOfRing<atomIDOfRing; i_atomIDOfRing++ ) {
       if( arrBondsOfRingCandidate[i_atomIDOfRing] == targetBond ) {
           isValidBond = rg_FALSE;
           break;
       }
   }
   
   return isValidBond;
}



rg_FLAG Molecule::isIdenticalRing( const rg_INT& numOfAtomsInRing, ChemicalBond** arrBondsOfRingA, ChemicalBond** arrBondsOfRingB )
{
   rg_FLAG isIdentical   = rg_TRUE;
   rg_INT  i_atomInRingA = numOfAtomsInRing-1;

   for( rg_INT i_atomInRingB=0; i_atomInRingB<numOfAtomsInRing; i_atomInRingB++ ) {
       if ( arrBondsOfRingA[i_atomInRingA] != arrBondsOfRingB[i_atomInRingB] ) {
           isIdentical = rg_FALSE;
           break;
       }
       i_atomInRingA--;
   }

   return isIdentical;   
}



bitset<MAX_CHEM_ATOM_NUM> Molecule::getBitSetKeyOfPharmaFeatureForMap( PharmaFeature* aPharmaFeature )
{
   bitset<MAX_CHEM_ATOM_NUM> bitSetKey;
   
   rg_dList<Atom*>* atomsInPharmaFeature = aPharmaFeature->getAtoms();

   atomsInPharmaFeature->reset4Loop();

   while( atomsInPharmaFeature->setNext4Loop() ) {
       bitSetKey.set( atomsInPharmaFeature->getEntity()->getID(), true );
   }

   return bitSetKey;
}

void Molecule::makePeoridicStructure()
{
	if(m_isCryst == false) {
		return;
	}

	Chain currChain = Chain(m_chains.getSize() + 1, this) ;
	
	Residue* currResidue = addResidue( Residue(m_residues.getSize() + 1) );
	currChain.addResidue( currResidue );
	
	currResidue->setSequenceNumber(m_residues.getSize() + 1);
	currResidue->setResidueName( "PRD" );
	
	int atomSerialNum = m_atoms.getSize() + 1;

	rg_dList<Atom> peoridicAtoms;	

	rg_Point3D latticeVectorX( m_cryst[0],  0.0,  0.0);
	rg_Point3D latticeVectorY( 0.0,  m_cryst[1],  0.0);
	rg_Point3D latticeVectorZ( 0.0,  0.0,  m_cryst[2]);

	rg_TMatrix3D tranformMatrix;

	tranformMatrix.makeIdentity();
	tranformMatrix.rotateZ((-(90-m_cryst[5]) * rg_PI) / 180.0);
	latticeVectorY = tranformMatrix * latticeVectorY;

	tranformMatrix.makeIdentity();
	tranformMatrix.rotateX((-(90.0-m_cryst[3]) * rg_PI) / 180.0);
	latticeVectorZ = tranformMatrix * latticeVectorZ;

	tranformMatrix.makeIdentity();
	tranformMatrix.rotateY((-(90.0-m_cryst[4]) * rg_PI) / 180.0);
	latticeVectorZ = tranformMatrix * latticeVectorZ;


	for(int i = 0; i < 3; i++){				
		for(int j = 0; j < 3; j++){
			for(int k = 0; k < 3; k++){
				if(i == 1 && j == 1 && k == 1) {
					continue;
				}

				m_atoms.reset4Loop();
				while(m_atoms.setNext4Loop()) {
					Atom currAtom = m_atoms.getEntity();
					currAtom.setID(atomSerialNum);
					currAtom.setSerialFromInputFile(atomSerialNum);
					atomSerialNum++;

					rg_Point3D delta((i - 1) * latticeVectorX + (j - 1) * latticeVectorY + (k - 1) * latticeVectorZ);
					rg_Point3D currPt = delta + currAtom.getAtomBall().getCenter();
					Sphere currBall(currPt, currAtom.getAtomBall().getRadius());

					currAtom.setAtomBall(currBall);						

					peoridicAtoms.add(currAtom);
				}				
			}
		}
	}

	peoridicAtoms.reset4Loop();
	while(peoridicAtoms.setNext4Loop()) {
		Atom currAtom = peoridicAtoms.getEntity();
		Atom* tempAtom = m_atoms.add(currAtom);

		tempAtom->setResidue(currResidue);
		currResidue->addAtom(tempAtom);
	}

	
	currResidue->setChain( addChain(currChain) );
}



//void Molecule::refineListOfPharmaFeatures()
//{
//   // 1. \C1ߺ\B9\B5\C8 PharmaFeature\B5\E9 \C0\E7\B0\C5 \C8\C4 \BA\B4\C7\D5
//   typedef map< string, PharmaFeature* > pharmaFeatureMap;
//	pharmaFeatureMap mapOfPharmaFeature;
//   pharmaFeatureMap::iterator i_mapOfPharmaFeature;
//
//   pair<pharmaFeatureMap::iterator, bool> mapElement;
//	
//   list<PharmaFeature>::iterator i_listOfPharmaFeatures = m_pharmaFeatures.begin();
//
//   while( i_listOfPharmaFeatures != m_pharmaFeatures.end() ) {
//       PharmaFeature* currPharmaFeature = (PharmaFeature*)(&(*i_listOfPharmaFeatures));
//
//       bitset<MAX_CHEM_ATOM_NUM> bitSetKey = getBitSetKeyOfPharmaFeatureForMap( currPharmaFeature );
//
//       mapElement = mapOfPharmaFeature.insert( pharmaFeatureMap::value_type( bitSetKey.to_string(), currPharmaFeature ) );
//
//       if( mapElement.second == false ) {
//
//           i_mapOfPharmaFeature = mapElement.first;
//           PharmaFeature* existPharmaFeature = (*i_mapOfPharmaFeature).second; // \B1\E2\C1\B8\BF\A1 \C0\FA\C0\E5\B5Ǿ\EE \C0ִ\C2 PharmaFeature\B8\A6 \B0\A1\C1\AE\BF\C8.
//
//           existPharmaFeature->mergePharmaFeatureType( currPharmaFeature->getPharmaFeatureType() );
//           existPharmaFeature->mergeChemFuncGroupType( currPharmaFeature->getChemFuncGroupType() );
//    
//           i_listOfPharmaFeatures = m_pharmaFeatures.erase( i_listOfPharmaFeatures );
//       }
//       else {
//           i_listOfPharmaFeatures++;
//       }
//   }
//
//   i_listOfPharmaFeatures = m_pharmaFeatures.begin();
//
//
//   // 2. \B1\E2\C1\B8\C0\C7 Atom\B5�?\BF\AC\B0\E1\B5Ǿ\EE \C0ִ\F8 PharmaFeature\B5\E9\C0\BB \BB\E8\C1\A6\C7Ѵ\D9.
//   m_atoms.reset4Loop();
//   while( m_atoms.setNext4Loop() ) {
//       m_atoms.getpEntity()->getpChemicalProperties()->getListOfPharmaFeatures()->removeAll();
//   }
//
//
//   // 3. PharmaFeatre\C0\C7 ID \BE\F7\B5\A5\C0\CCƮ \B9\D7 Atom\B5�?B0\D4 PharmaFeature \BC\C2\C6\C3
//   rg_INT IDOfUpdatedPharmaFeature = 0;
//   while( i_listOfPharmaFeatures != m_pharmaFeatures.end() ) {
//       PharmaFeature* currPharmaFeature = (PharmaFeature*)(&(*i_listOfPharmaFeatures));
//       currPharmaFeature->setID( IDOfUpdatedPharmaFeature );
//       IDOfUpdatedPharmaFeature++;
//
//       rg_dList<Atom*>* listOfAtoms =  currPharmaFeature->getAtoms();
//       listOfAtoms->reset4Loop();
//       while( listOfAtoms->setNext4Loop() ) {
//           Atom* currAtomInFeature = listOfAtoms->getEntity();
//           currAtomInFeature->getpChemicalProperties()->addPharmaFeature( currPharmaFeature );
//       }
//       i_listOfPharmaFeatures++;
//   }
//}
//
//
// void Molecule::setOOBBoxPoints( rg_Point3D* boxPoints )
// {
//     for(rg_INT i=0; i<8; i++) {
//         m_OOBBoxPoints[i] = boxPoints[i];
//     }
// }
// 
// 
// 
// rg_Point3D* Molecule::getOOBBoxPoints()
// {
//     return m_OOBBoxPoints;
// }

