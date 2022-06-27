#include "PocketPrimitive.h"
#include "ConstForBetaComplex.h"



//#include "GeometryForSphere.h"
/////////////////////////////////////////////////////////////////////////////////
//
// Constructor
PocketPrimitive::PocketPrimitive()
{
	m_sumOfAreaOfMBSFaces = 0.;

	m_chemicalPropertiesOfPocket.setIsThisComputed( rg_FALSE );
}




PocketPrimitive::PocketPrimitive(const PocketPrimitive& tempPocketPrimitive)
{
    m_ID = tempPocketPrimitive.m_ID;

	//m_faceOnOuterBetaShape = tempPocketPrimitive.m_faceOnOuterBetaShape;
	
	m_MBSFaces               = tempPocketPrimitive.m_MBSFaces;
	m_MBSVertices            = tempPocketPrimitive.m_MBSVertices;
    m_sumOfAreaOfMBSFaces    = tempPocketPrimitive.m_sumOfAreaOfMBSFaces;

	m_chemicalPropertiesOfPocket = tempPocketPrimitive.m_chemicalPropertiesOfPocket;
}




PocketPrimitive::~PocketPrimitive()
{
}




/////////////////////////////////////////////////////////////////////////////////
//
// Get Functions
// BetaFace* PocketPrimitive::getFacesOnOuterBetaShape()
// {
// 	return m_faceOnOuterBetaShape;
// }





rg_dList<MBSFace*>* PocketPrimitive::getMBSFaces()
{
	return &m_MBSFaces;
}




rg_dList<MBSVertex*>* PocketPrimitive::getMBSVertices()
{
    /*
    if ( m_MBSVertices.getSize() == 0 )  {
		fillMBSVertexList();
    }
    */

    fillMBSVertexList();
	return &m_MBSVertices;
}


void PocketPrimitive::fillMBSVertexList()
{
    m_MBSVertices.removeAll();

	m_MBSFaces.reset4Loop();
	while ( m_MBSFaces.setNext4Loop() )		{
		MBSFace* currFace = m_MBSFaces.getEntity();

		rg_dList< MBSVertex* > boundingVertices;
		currFace->searchBoundingVertices( boundingVertices );

		boundingVertices.reset4Loop();
		while ( boundingVertices.setNext4Loop() )	{
			MBSVertex* currBoundingVertex = boundingVertices.getEntity();

			m_MBSVertices.addWithoutSame( currBoundingVertex );
		}
	}
}





rg_INT  PocketPrimitive::getNumMBSFaces() const
{
    return m_MBSFaces.getSize();
}



rg_REAL PocketPrimitive::getMaximumDepth() const
{
    return m_maxDepth;
}


rg_REAL PocketPrimitive::getAverageDepth() const
{
    return m_aveDepth;
}


ChemicalPropertiesOfPocket*	PocketPrimitive::getpChemicalPropertiesOfPocket()
{
	if ( m_chemicalPropertiesOfPocket.isThisComputed() == rg_FALSE )
		computeChemicalProperties();

	return &m_chemicalPropertiesOfPocket;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Set Functions
// void PocketPrimitive::setFaceOnOuterBetaShape( BetaFace* tempFaceOfOuterBetaShape)
// {
// 	m_faceOnOuterBetaShape = tempFaceOfOuterBetaShape;
// }





void PocketPrimitive::addMBSFaceToThisPrimitive(MBSFace* newFace)
{
    m_MBSFaces.add( newFace );
}




void PocketPrimitive::addMBSVertexToThisPrimitive(MBSVertex* newVertex)
{
    m_MBSVertices.add( newVertex );
}





/////////////////////////////////////////////////////////////////////////////////
//
// Operator Overloadings
PocketPrimitive& PocketPrimitive::operator =(const PocketPrimitive& tempPocketPrimitive)
{
	if( this == &tempPocketPrimitive )
		return *this;

    m_ID = tempPocketPrimitive.m_ID;

    //m_faceOnOuterBetaShape = tempPocketPrimitive.m_faceOnOuterBetaShape;
	
	m_MBSFaces               = tempPocketPrimitive.m_MBSFaces;
	m_MBSVertices            = tempPocketPrimitive.m_MBSVertices;
    m_sumOfAreaOfMBSFaces    = tempPocketPrimitive.m_sumOfAreaOfMBSFaces;

	m_chemicalPropertiesOfPocket = tempPocketPrimitive.m_chemicalPropertiesOfPocket;

	return *this;
}




rg_FLAG PocketPrimitive::isItLengthZeroEdge( MBSEdge* edge )
{
    if ( edge->getStartVertex()->getOriginalBetaVertex() == edge->getEndVertex()->getOriginalBetaVertex() )
        return rg_TRUE;
    else
        return rg_FALSE;
}







rg_dList< void* >* PocketPrimitive::getBallPropertyList()
{
    /*
    if ( m_ballProperties.getSize() == 0 ) {
		fillBallPropertyList();    
    }
    */

    fillBallPropertyList();    
    return &m_ballProperties;
}


void PocketPrimitive::fillBallPropertyList()
{
    m_ballProperties.removeAll();

    /*
    if ( m_MBSVertices.getSize() == 0 ) {
		fillMBSVertexList();
    }
    */
    fillMBSVertexList();

	MBSVertex* currMBSVertex = rg_NULL;
    //Ball* currGeometry = rg_NULL;
    m_MBSVertices.reset4Loop();    
    while ( m_MBSVertices.setNext4Loop() )
    {
        currMBSVertex = m_MBSVertices.getEntity();

        BetaVertex* currBetaVertex = currMBSVertex->getOriginalBetaVertex();

        
		Ball* currGeometry = currBetaVertex->getBallProperty();  // FOR LINUX SYSTEM

		//#ifdef WIN32
		//	Particle* currGeometry = (Particle*) currBetaVertex->getBallProperty()->getProperty();  // FOR WINDOWS SYSTEM //BETA-MOL
  //      #else        
		//	Ball* currGeometry = currBetaVertex->getBallProperty();  // FOR LINUX SYSTEM
		//#endif

        m_ballProperties.addWithoutSame( currBetaVertex->getProperty() );        
    }

    // By Joonghyun on Jan. 21 2015
    // this function does not work because we use Atom in MGOS
	//sortPDBAtomsWRTAtomSerialNum();
}


void PocketPrimitive::sortPDBAtomsWRTAtomSerialNum()
{
	rg_INT n = m_ballProperties.getSize();
	bool swapped = true;

	do 
	{
		swapped = false;
		n--;

		m_ballProperties.reset4Loop();
		for ( int i=0; i < n; i++)
		{
			m_ballProperties.setNext4Loop();
			
			Atom* currAtom = (Atom*)m_ballProperties.getEntity();
			Atom* nextAtom = (Atom*)m_ballProperties.getNextEntity();

			if ( currAtom->getSerialFromInputFile() > nextAtom->getSerialFromInputFile() )
			{
				swapNode( m_ballProperties.getCurrentpNode(), m_ballProperties.getCurrentpNode()->getNext() );
				swapped = true;
				//swap하고 나서 currentNode가 제자리에 위치해 있는가 확인.
			}
		}	
	} while( swapped );
}

void PocketPrimitive::swapNode( rg_dNode< void* >* node1, rg_dNode< void* >* node2 )
{ //data가 pointer형이므로 value를 swap하는 것이 더 효율적이다.
	void* temp = node1->getEntity();
	node1->setEntity( node2->getEntity() );
	node2->setEntity( temp);
}

void PocketPrimitive::removeAll()
{
    m_MBSVertices.removeAll();
    m_MBSFaces.removeAll();
    m_ballProperties.removeAll();

    m_sumOfAreaOfMBSFaces = 0.;
	m_chemicalPropertiesOfPocket.setIsThisComputed( rg_FALSE );
}



rg_REAL PocketPrimitive::getAreaOfMBSFaces()
{
    if ( m_sumOfAreaOfMBSFaces != 0 )
        return m_sumOfAreaOfMBSFaces;

    m_MBSFaces.reset4Loop();
    while ( m_MBSFaces.setNext4Loop() )  {
        MBSFace* currMBSFace = m_MBSFaces.getEntity();

        m_sumOfAreaOfMBSFaces += computeAreaOfMBSFace( currMBSFace );
    }

    return m_sumOfAreaOfMBSFaces;
}



rg_REAL PocketPrimitive::computeAreaOfMBSFace( MBSFace* givenMBSFace )
{
    rg_dList< MBSVertex* > boundingVertices;
    givenMBSFace->searchBoundingVertices( boundingVertices );

    boundingVertices.reset4Loop();
    Ball* boundingSpheres[3];
    rg_INT i;
    for ( i=0; i<3; i++ )  {
        boundingVertices.setNext4Loop();

        boundingSpheres[i] = boundingVertices.getEntity()->getOriginalBetaVertex()->getBallProperty();
    }

    if ( boundingSpheres[0] == boundingSpheres[1] 
      || boundingSpheres[0] == boundingSpheres[2] 
      || boundingSpheres[1] == boundingSpheres[2] 
      )                                             //area = 0
      return 0;

    rg_Point3D vector1 = boundingSpheres[1]->getGeometry().getCenter()
                       - boundingSpheres[0]->getGeometry().getCenter();
    rg_Point3D vector2 = boundingSpheres[2]->getGeometry().getCenter()
                       - boundingSpheres[0]->getGeometry().getCenter();

    rg_Point3D crossProductOfVectors = vector1.crossProduct( vector2 );
    return crossProductOfVectors.magnitude() / 2.0;
}



rg_FLAG PocketPrimitive::isAreaOfMBSFaceZero( MBSFace* givenMBSFace )
{
    rg_FLAG isAreaZero = rg_FALSE;

    rg_dList< MBSVertex* > boundingVertices;
    givenMBSFace->searchBoundingVertices( boundingVertices );

    boundingVertices.reset4Loop();
    Ball* boundingSpheres[3];
    rg_INT i;
    for ( i=0; i<3; i++ )  {
        boundingVertices.setNext4Loop();

        boundingSpheres[i] = boundingVertices.getEntity()->getOriginalBetaVertex()->getBallProperty();
    }

    if ( boundingSpheres[0] == boundingSpheres[1] 
      || boundingSpheres[0] == boundingSpheres[2] 
      || boundingSpheres[1] == boundingSpheres[2] 
      )                                             //area = 0
      isAreaZero = rg_TRUE;

    return isAreaZero;
}


void PocketPrimitive::computeChemicalProperties()
{
    /*
	if ( m_ballProperties.getSize() == 0 )
		fillBallPropertyList();

	rg_INT numOfCarbonAtoms = 0;
	rg_INT numOfNitrogenAtoms = 0;
	rg_INT numOfOxygenAtoms = 0;	
	rg_INT numOfSulfurAtoms = 0;

	rg_INT numOfBackboneAtoms = 0;
	rg_INT numOfSideChainAtoms = 0;

	rg_INT numOfHydrophobicResidues = 0;
	rg_INT numOfPolarResidues = 0;
	rg_INT numOfChargedResidues = 0;

	rg_INT numOfAlphaHelices = 0;
	rg_INT numOfBetaSheets = 0;
	rg_INT numOfLoops = 0;

	rg_INT numOfHBondDonors = 0;
	rg_INT numOfHBondAcceptor = 0;

	

	Atom* currPDBAtom = rg_NULL;
	m_ballProperties.reset4Loop();
	while ( m_ballProperties.setNext4Loop() )
	{
		currPDBAtom = (Atom*)( m_ballProperties.getEntity() );

		rg_INT atomicNumber = currPDBAtom->getAtomicNumber();
		rg_INT residueTypeNumber = currPDBAtom->getResidue()->getID();
		rg_INT remotenessIndicator = currPDBAtom->getRemotenessIndicator(); 
					// 0:NONE, 1:ALPHA, 2:BETA, 3:GAMMA, 4:DELTA, 5:EPSILON, 6:ZETA, 7:ETA 



		switch( atomicNumber )
		{
		case 6: //carbon
			numOfCarbonAtoms++;
			break;

		case 7: //nitrogen
			numOfNitrogenAtoms++;
			break;

		case 8: //oxygen
			numOfOxygenAtoms++;
			break;
		case 16: //sulfur
			numOfSulfurAtoms++;
			break;

		default:
			break;
		}
					   	
		if ( currPDBAtom->isOnBackBone() == rg_TRUE )
			numOfBackboneAtoms++;
		else
			numOfSideChainAtoms++;

		switch( residueTypeNumber )
		{
		case 1: //Ala - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 3: //Cys - Polar
			numOfPolarResidues++;
			break;

		case 4: //Asp - Charged
			numOfChargedResidues++;
			break;

		case 5: //Glu - Charged
			numOfChargedResidues++;
			break;

		case 6: //Phe - Hydrophobic
			numOfHydrophobicResidues++;
			break;
		
		case 7: // Gly - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 8: //His - Polar
			numOfPolarResidues++;
			break;				

		case 9: // Ile - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 10: //Lys - Charged
			numOfChargedResidues++;
			break;

		case 11: //Leu - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 12: //Met - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 13: // Asn - Polar
			numOfPolarResidues++;
			break;

		case 14: //Pro - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 15: //Gln - Polar
			numOfPolarResidues++;
			break;

		case 16: //Arg - Charged
			numOfChargedResidues++;
			break;

		case 17: //Ser - Polar
			numOfPolarResidues++;
			break;

		case 18: //Thr - Polar
			numOfPolarResidues++;
			break;			

		case 19: //Val - Hydrophobic
			numOfHydrophobicResidues++;
			break;

		case 20: //Trp - Polar
			numOfPolarResidues++;
			break;

		case 21: //Tyr - Polar
			numOfPolarResidues++;
			break;

		default:
			break;	
		}

		//switch( currPDBAtom->getResidue()->getSecondaryStructure() ) //아직 계산이 안됨.

		//hydrogen bonding
			//H-donor
		if ( atomicNumber == 7) //nitrogen
		{
			if ( residueTypeNumber == 8 ) //HIS
			{
				if ( remotenessIndicator != 4 //Delta
				  && remotenessIndicator != 5) //Epsilon
				{
					numOfHBondDonors++;
				}
			}
			else
				numOfHBondDonors++;
		}
		else if ( atomicNumber == 16) //sulfur
		{
			if ( residueTypeNumber == 3) //Cys
				numOfHBondDonors++;
		}
		else if ( atomicNumber == 8) //oxygen
		{
			if ( remotenessIndicator == 3 ) //Gamma
			{
				if ( residueTypeNumber == 17 //Ser
				  || residueTypeNumber == 18 ) //Thr
					numOfHBondDonors++;
			}
			else if ( remotenessIndicator == 7) //ETA(H)
			{
				if ( residueTypeNumber == 21 ) //Tyr)
					numOfHBondDonors++;
			}
		}

		//hydrogen bonding
			//H-acceptor
		if ( atomicNumber == 8 ) //oxygen
			numOfHBondAcceptor++;
		else if ( atomicNumber == 16 ) //sulfur
		{
			if ( residueTypeNumber == 12) //Met
				numOfHBondAcceptor++;
		}		
	}

	
	m_chemicalPropertiesOfPocket.setIsThisComputed( rg_TRUE );
	
	m_chemicalPropertiesOfPocket.setNumOfTotalAtoms( m_ballProperties.getSize() );

	m_chemicalPropertiesOfPocket.setNumOfCarbonAtoms( numOfCarbonAtoms );
	m_chemicalPropertiesOfPocket.setNumOfNitrogenAtoms( numOfNitrogenAtoms );
	m_chemicalPropertiesOfPocket.setNumOfOxygenAtoms( numOfOxygenAtoms );
	m_chemicalPropertiesOfPocket.setNumOfSulfurAtoms( numOfSulfurAtoms );

	m_chemicalPropertiesOfPocket.setNumOfBackboneAtoms( numOfBackboneAtoms );
	m_chemicalPropertiesOfPocket.setNumOfSideChainAtoms( numOfSideChainAtoms );

	m_chemicalPropertiesOfPocket.setNumOfHydrophobicResidues( numOfHydrophobicResidues );
	m_chemicalPropertiesOfPocket.setNumOfPolarResidues( numOfPolarResidues );
	m_chemicalPropertiesOfPocket.setNumOfChargedResidues( numOfChargedResidues );

	m_chemicalPropertiesOfPocket.setNumOfAlphaHelices( numOfAlphaHelices );
	m_chemicalPropertiesOfPocket.setNumOfBetaSheets( numOfBetaSheets );
	m_chemicalPropertiesOfPocket.setNumOfLoops( numOfLoops );

	m_chemicalPropertiesOfPocket.setNumOfHBondDonors( numOfHBondDonors );
	m_chemicalPropertiesOfPocket.setNumOfHBondAcceptor( numOfHBondAcceptor );
*/
}

