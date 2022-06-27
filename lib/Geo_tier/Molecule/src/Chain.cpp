#include "Chain.h"
#include "rg_Molecule.h"
#include "Residue.h"
#include "FunctionsForMolecule.h"


V::GeometryTier::Chain::Chain()
: m_ID(0), m_chainCode(UNK_CHAIN), m_molecule(rg_NULL), m_chainIDFromInputFileInDecimal(32)
{
}



V::GeometryTier::Chain::Chain( const rg_INT& ID )
: m_ID(ID), m_chainCode(UNK_CHAIN), m_molecule(rg_NULL), m_chainIDFromInputFileInDecimal(32)
{   
}



V::GeometryTier::Chain::Chain( const rg_INT& ID, Molecule* aMolecule )
: m_ID(ID)
, m_chainCode(UNK_CHAIN)
, m_molecule(aMolecule)
, m_chainIDFromInputFileInDecimal(32)
{
}



V::GeometryTier::Chain::Chain( const rg_INT& ID, const ChainCode& aChainCode )
: m_ID(ID)
, m_chainCode(aChainCode)
, m_molecule(rg_NULL)
, m_chainIDFromInputFileInDecimal(32)
{
}



V::GeometryTier::Chain::Chain( const rg_INT& ID, const ChainCode& aChainCode, const rg_INT& chainIDFromInput )
: m_ID(ID)
, m_chainCode(aChainCode)
, m_molecule(rg_NULL)
, m_chainIDFromInputFileInDecimal(chainIDFromInput)
{   
}



V::GeometryTier::Chain::Chain( const V::GeometryTier::Chain& aChain )
{
    m_ID                              = aChain.m_ID;
    m_chainIDFromInputFileInDecimal   = aChain.m_chainIDFromInputFileInDecimal;
    
    m_molecule              = aChain.m_molecule;
    m_residues              = aChain.m_residues;

    m_chainCode             = aChain.m_chainCode;
    m_secondaryStructure    = aChain.m_secondaryStructure;
}



V::GeometryTier::Chain::~Chain()
{
}







string V::GeometryTier::Chain::getChainIDFromInputFileInString() const
{
    char   tempChar = m_chainIDFromInputFileInDecimal;
    string chainID( &tempChar, 1 );

    return chainID;
}



//  by Youngsong Cho, 2012.03.14.
rg_INT V::GeometryTier::Chain::getNumStandardResidues() const
{
    rg_INT numStandardResidues = 0;

    m_residues.reset4Loop();
    while ( m_residues.setNext4Loop() ) {
        Residue* currRes = m_residues.getEntity();

        if ( currRes->isStandardResidue() ) {
            numStandardResidues++;
        }
    }

    return numStandardResidues;
}



rg_BOOL V::GeometryTier::Chain::isLigandChain() const
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LIGAND: Druglikeness molecule according to Lipinski's rule of five                                 //
    //                                                                                                    //
    // * Summary of Lipinski's rule of five                                                               //
    // - Not more than 5 hydrogen bond donors (nitrogen or oxygen atoms with one or more hydrogen atoms)  //
    // - Not more than 10 hydrogen bond acceptors (nitrogen or oxygen atoms)                              //
    // - A molecular weight under 500 daltons                                                             //
    // - An octanol-water partition coefficient[3] log P less than 5                                      //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    return rg_FALSE;
}




rg_BOOL V::GeometryTier::Chain::haveMissingResidues()
{
	rg_dList<Residue*> residuesOrderedBySeqNumber;
	getAminoResiduesSortedBySequenceNumber(residuesOrderedBySeqNumber);

	rg_BOOL bHaveMissingResidues = rg_FALSE;

	if(residuesOrderedBySeqNumber.getSize() > 0)
	{
		rg_INT prevSeqNum = residuesOrderedBySeqNumber.getFirstEntity()->getSequenceNumber() - 1;
		residuesOrderedBySeqNumber.reset4Loop();
		while (residuesOrderedBySeqNumber.setNext4Loop())
		{
			Residue* currResidue = rg_NULL;
			currResidue = residuesOrderedBySeqNumber.getEntity();
			rg_INT currSeqNum = -1;
			currSeqNum = currResidue->getSequenceNumber();
			if(currSeqNum != prevSeqNum + 1)
			{
				bHaveMissingResidues = rg_TRUE;
				break;
			}
			else
			{
				prevSeqNum = currSeqNum;
			}
		}
	}
	
	return bHaveMissingResidues;
}


rg_INT V::GeometryTier::Chain::findMissingResidues(rg_dList<rg_INT>& seqNumbersOfMissingResidues)
{
	rg_dList<Residue*> residuesOrderedBySeqNumber;
	getAminoResiduesSortedBySequenceNumber(residuesOrderedBySeqNumber);

	rg_BOOL bHaveMissingResidues = rg_FALSE;

	rg_INT prevSeqNum = residuesOrderedBySeqNumber.getFirstEntity()->getSequenceNumber() - 1;
	residuesOrderedBySeqNumber.reset4Loop();
	while (residuesOrderedBySeqNumber.setNext4Loop())
	{
		Residue* currResidue = rg_NULL;
		currResidue = residuesOrderedBySeqNumber.getEntity();
		rg_INT currSeqNum = -1;
		currSeqNum = currResidue->getSequenceNumber();
		if(currSeqNum > prevSeqNum + 1)
		{
			rg_INT numMissingResidues = currSeqNum - (prevSeqNum + 1);
			rg_INDEX i;
			for (i = 0;i < numMissingResidues;i++)
			{
				seqNumbersOfMissingResidues.add( prevSeqNum + 1 + i );
			}
		}
		prevSeqNum = currSeqNum;
	}

	return seqNumbersOfMissingResidues.getSize();
}


void V::GeometryTier::Chain::getResiduesSortedBySequenceNumber( rg_dList<Residue*>& targetResidues )
{
    rg_INT    numOfResidues  = m_residues.getSize();
    Residue** sortedResidues = m_residues.getArray();

    qsort( (void *)sortedResidues, numOfResidues, sizeof(Residue*), compareResidueBySequenceNumber );

    for ( rg_INT i_residue=0; i_residue<numOfResidues; i_residue++ ) {
        targetResidues.addTail( sortedResidues[i_residue] );
    }
    
    delete [] sortedResidues;
}



void V::GeometryTier::Chain::getAminoResiduesSortedBySequenceNumber( rg_dList<Residue*>& aminoResiduesSortedBySeqNumbers )
{
	// collect amino residues
	rg_dList<Residue*> aminoResidues;
	m_residues.reset4Loop();
	while (m_residues.setNext4Loop())
	{
		Residue* currResidue = rg_NULL;
		currResidue = m_residues.getEntity();
		if(currResidue->isAminoResidue())
			aminoResidues.add( currResidue );
	}

	// get the collection in array
	rg_INT numResidues = aminoResidues.getSize();
	Residue** sortedAminoResidues = rg_NULL;
	sortedAminoResidues = aminoResidues.getArray();

	// sort
	qsort( (void *)sortedAminoResidues, numResidues, sizeof(Residue*), compareResidueBySequenceNumber );

	for ( rg_INT i=0; i<numResidues; i++ ) 
	{
		aminoResiduesSortedBySeqNumbers.addTail( sortedAminoResidues[i] );
	}

	if(sortedAminoResidues != rg_NULL)
		delete [] sortedAminoResidues;	
}


rg_INT V::GeometryTier::Chain::getAminoResiduesSortedBySequenceNumber( Residue**& aminoResiduesSortedBySeqNumbers )
{
    // collect amino residues
    rg_dList<Residue*> aminoResidues;
    m_residues.reset4Loop();
    while (m_residues.setNext4Loop())
    {
        Residue* currResidue = rg_NULL;
        currResidue = m_residues.getEntity();
        if(currResidue->isAminoResidue())
            aminoResidues.add( currResidue );
    }

    // get the collection in array
    rg_INT numResidues = aminoResidues.getSize();
    aminoResiduesSortedBySeqNumbers= rg_NULL;
    aminoResiduesSortedBySeqNumbers = aminoResidues.getArray();

    // sort
    qsort( (void *)aminoResiduesSortedBySeqNumbers, numResidues, sizeof(Residue*), compareResidueBySequenceNumber );

    return numResidues;
}


void V::GeometryTier::Chain::getResidueWithSSCsSortedBySequenceNumber( rg_dList<ResidueWithSSCode>& targetResidueWithSSCs )
{
    if ( m_chainCode != PROTEIN_CHAIN )
        return;

    rg_INT    numOfResidues  = m_residues.getSize();        
    Residue** sortedResidues = m_residues.getArray();

    qsort( (void *)sortedResidues, numOfResidues, sizeof(Residue*), compareResidueBySequenceNumber );

    for ( rg_INT i_residue=0; i_residue<numOfResidues; i_residue++ ) {
        
        ResidueWithSSCode* pResidueWithSSCode = rg_NULL;

        Helix* aHelix = m_secondaryStructure.getHelixOfResidue( sortedResidues[i_residue] );
        if( aHelix != rg_NULL ) {
            pResidueWithSSCode = targetResidueWithSSCs.addTail( ResidueWithSSCode( SSC_HELIX, aHelix, sortedResidues[i_residue] ) );
            continue;
        }

        Sheet* aSheet = m_secondaryStructure.getSheetOfResidue( sortedResidues[i_residue] );
        if( aSheet != rg_NULL ) {
            pResidueWithSSCode = targetResidueWithSSCs.addTail( ResidueWithSSCode( SSC_SHEET, aSheet, sortedResidues[i_residue] ) );
            continue;
        }

        Turn* aTurn = m_secondaryStructure.getTurnOfResidue( sortedResidues[i_residue] );
        if( aTurn != rg_NULL ) {
            pResidueWithSSCode = targetResidueWithSSCs.addTail( ResidueWithSSCode( SSC_TURN, aTurn, sortedResidues[i_residue] ) );
            continue;
        }

        if ( pResidueWithSSCode == rg_NULL ) {
            pResidueWithSSCode = targetResidueWithSSCs.addTail( ResidueWithSSCode( SSC_UNK, rg_NULL, sortedResidues[i_residue] ) );
        }
    }

    delete [] sortedResidues;
}



void V::GeometryTier::Chain::getAtomsOnBackboneSortedBySequenceNumber( rg_dList<Atom*>& targetAtoms )
{
    rg_dList<Residue*> targetResidues;
    getResiduesSortedBySequenceNumber( targetResidues );

    targetResidues.reset4Loop();
    while ( targetResidues.setNext4Loop() ) {      

        rg_dList<Atom*> atomsOnBackboneInResidue;
        targetResidues.getEntity()->getAtomsOnBackbone( &atomsOnBackboneInResidue );

        atomsOnBackboneInResidue.reset4Loop();
        while ( atomsOnBackboneInResidue.setNext4Loop() ) {
            targetAtoms.addTail( atomsOnBackboneInResidue.getEntity() );
        }
    }
}




void V::GeometryTier::Chain::evaluateAndSetChainCode()
{
    rg_INT numOfResidueInChain = m_residues.getSize();

    if( numOfResidueInChain == 0 )
        return;

    rg_INT firstResidueCode = m_residues.getFirstEntity()->getResidueCode();

    if( firstResidueCode >= ALA_AMINO_RESIDUE && firstResidueCode <= GLX_AMINO_RESIDUE ) {
        m_chainCode = PROTEIN_CHAIN;
    }
    else if( firstResidueCode == HOH_RESIDUE ) {
        m_chainCode = HOH_CHAIN;
    }
    else if( firstResidueCode == UNK_RESIDUE ) {
        
        string    thymidineOfDNA( RESIDUE_FEATURES[T_DNA_RESIDUE].threeCodeName );
        string modThymidineOfDNA( RESIDUE_FEATURES[TM_DNA_RESIDUE].threeCodeName );
        
        string    uridineOfRNA( RESIDUE_FEATURES[U_RNA_RESIDUE].threeCodeName );
        string modUridineOfRNA( RESIDUE_FEATURES[UM_RNA_RESIDUE].threeCodeName );

        
        m_residues.reset4Loop();
        while ( m_residues.setNext4Loop() ) {
            string threeCodeNameOfResidue = m_residues.getEntity()->getResidueName();
            
            if ( threeCodeNameOfResidue.compare( thymidineOfDNA )    == 0 || 
                 threeCodeNameOfResidue.compare( modThymidineOfDNA ) == 0  ) {
                
                m_chainCode = DNA_CHAIN;
                break;
            }
            else if ( threeCodeNameOfResidue.compare( uridineOfRNA )    == 0 || 
                      threeCodeNameOfResidue.compare( modUridineOfRNA ) == 0  ) {
                
                m_chainCode = RNA_CHAIN;
                break;
            }
        }
    }
    else {
        m_chainCode = UNK_CHAIN;
    }
}



V::GeometryTier::Residue* V::GeometryTier::Chain::addResidue(V::GeometryTier::Residue* aResidue )
{
    return *m_residues.addTail( aResidue );
}



V::GeometryTier::Chain& V::GeometryTier::Chain::operator=( const V::GeometryTier::Chain& aChain )
{
    if( this == &aChain )
        return *this;

    m_ID                              = aChain.m_ID;
    m_chainIDFromInputFileInDecimal   = aChain.m_chainIDFromInputFileInDecimal;
    
    m_molecule              = aChain.m_molecule;
    m_residues          = aChain.m_residues;

    m_chainCode             = aChain.m_chainCode;
    m_secondaryStructure    = aChain.m_secondaryStructure;
        

    return *this;
}
