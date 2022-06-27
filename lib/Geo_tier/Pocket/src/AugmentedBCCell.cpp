#include "AugmentedBCCell.h"

AugmentedBCCell::AugmentedBCCell()
: BetaCell()  
{
	m_boundingState = rg_UNKNOWN;	//EXTRANEOUS_SIMPLEX와 같으므로 주의.
	m_isArtificial  = rg_FALSE;
}

AugmentedBCCell::AugmentedBCCell(const rg_INT& id)
: BetaCell(id)
{
	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;
}

AugmentedBCCell::AugmentedBCCell( AugmentedBCFace* f0, AugmentedBCFace* f1, AugmentedBCFace* f2, AugmentedBCFace* f3,   
								 AugmentedBCVertex* v0, AugmentedBCVertex* v1, AugmentedBCVertex* v2, AugmentedBCVertex* v3 ) 						     	  
: BetaCell()
{
	BetaCell::setFace( 0, (BetaFace*)f0 );
	BetaCell::setFace( 1, (BetaFace*)f1 );
	BetaCell::setFace( 2, (BetaFace*)f2 );
	BetaCell::setFace( 3, (BetaFace*)f3 );

	BetaCell::setVertex( 0, (BetaVertex*)v0 );
	BetaCell::setVertex( 1, (BetaVertex*)v1 );
	BetaCell::setVertex( 2, (BetaVertex*)v2 );
	BetaCell::setVertex( 3, (BetaVertex*)v3 );

	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;
}


AugmentedBCCell::AugmentedBCCell( AugmentedBCFace* f0, AugmentedBCFace* f1, AugmentedBCFace* f2, AugmentedBCFace* f3 )
: BetaCell()
{
	BetaCell::setFace( 0, (BetaFace*)f0 );
	BetaCell::setFace( 1, (BetaFace*)f1 );
	BetaCell::setFace( 2, (BetaFace*)f2 );
	BetaCell::setFace( 3, (BetaFace*)f3 );

	BetaCell::setVertex( 0, findMateVertex(f0) );
	BetaCell::setVertex( 1, findMateVertex(f1) );
	BetaCell::setVertex( 2, findMateVertex(f2) );
	BetaCell::setVertex( 3, findMateVertex(f3) );

	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;
}
	

AugmentedBCCell::AugmentedBCCell(const AugmentedBCCell& augBCCell)
: BetaCell(augBCCell)
{
	m_boundingState = augBCCell.m_boundingState;
	m_isArtificial  = augBCCell.m_isArtificial;
}

AugmentedBCCell::~AugmentedBCCell()
{
}


rg_INT	AugmentedBCCell::getBoundingState() const //CHECK! overloading잘 되는지 check
{
	return m_boundingState;
}

rg_BOOL AugmentedBCCell::isArtificial() const
{
	return m_isArtificial;
}


void AugmentedBCCell::setBoundingState(const rg_INT& boundingState)
{
	m_boundingState = boundingState;
}


void AugmentedBCCell::setIsArtificial( const rg_BOOL& isArtificial)
{
	m_isArtificial = isArtificial;
}


AugmentedBCCell& AugmentedBCCell::operator =(const AugmentedBCCell& augBCCell)
{
	if ( this == &augBCCell )
		return *this;

	BetaCell::operator=(augBCCell);

	m_boundingState = augBCCell.m_boundingState;
	m_isArtificial  = augBCCell.m_isArtificial;
	
	return *this;
}




AugmentedBCFace* AugmentedBCCell::findNextIncidentFaceOfEdge(AugmentedBCEdge* givenEdge, AugmentedBCFace* firstIncidentFace)
{
	for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( BetaCell::getFace(i) == (BetaFace*)firstIncidentFace )
            continue;

        if ( BetaCell::getFace(i)->isThere((AugmentedBCEdge*)givenEdge) == rg_TRUE )  
            return (AugmentedBCFace*)(BetaCell::getFace(i));
    }

    return rg_NULL;
}


void AugmentedBCCell::findFacesIncidentToEdge(AugmentedBCEdge* givenEdge, AugmentedBCFace** incidentFaces)
{
	rg_INT index4Face = 0;

    for ( rg_INT i=0; i<4; i++ )  {
		AugmentedBCFace* currFace = (AugmentedBCFace*)BetaCell::getFace(i);
        if ( currFace->isThere((BetaEdge*)givenEdge) == rg_TRUE )  {
            incidentFaces[index4Face++] = currFace;
        }
    }
}


BetaVertex* AugmentedBCCell::findMateVertex( AugmentedBCFace* givenFace )
{
    if ( this->isThere( givenFace) == rg_FALSE )
        return rg_NULL;


    BetaVertex* mateVertex = rg_NULL;

    for ( rg_INT i=0; i < EIWDS_NUM_FACE_ON_CELL; i++) {
        AugmentedBCFace* currFace = (AugmentedBCFace*)BetaCell::getFace(i);
        if ( givenFace == currFace ) {
            continue;
        }

        rg_dList< BetaVertex* > boundingVertices;
        currFace->searchFiniteVerticesInIntraWorld( boundingVertices );

        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            BetaVertex* currVertex = boundingVertices.getEntity();

            if ( givenFace->isThere( currVertex ) == rg_FALSE ) {
                mateVertex = currVertex;
                break;
            }
        }

        if ( mateVertex != rg_NULL ) {
            break;
        }
    }

    return mateVertex;
}
