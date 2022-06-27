#include "AugmentedBCFace.h"

AugmentedBCFace::AugmentedBCFace()
: BetaFace()  
{
	m_boundingState = rg_UNKNOWN;	//EXTRANEOUS_SIMPLEX와 같으므로 주의.
	m_isArtificial  = rg_FALSE;

    m_indexOfExteriorRegionToLeftCell   = -1;
    m_indexOfExteriorRegionToRightCell  = -1;
}

AugmentedBCFace::AugmentedBCFace(const rg_INT& id)
: BetaFace(id)
{
	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;

    m_indexOfExteriorRegionToLeftCell   = -1;
    m_indexOfExteriorRegionToRightCell  = -1;
}

AugmentedBCFace::AugmentedBCFace(AugmentedBCEdge* e1, AugmentedBCEdge* e2, AugmentedBCEdge* e3 )
: BetaFace()
{
	BetaFace::setEdge(0, (BetaEdge*)e1);
	BetaFace::setEdge(1, (BetaEdge*)e2);
	BetaFace::setEdge(2, (BetaEdge*)e3);

	if ( BetaFace::getEdge(1)->isThere( BetaFace::getEdge(0)->getEndVertex() ) == rg_TRUE )
		BetaFace::setEdgeOrientation( 0, rg_TRUE );
    else
		BetaFace::setEdgeOrientation( 0, rg_FALSE );

    if ( BetaFace::getEdge(2)->isThere( BetaFace::getEdge(1)->getEndVertex() ) == rg_TRUE )
		BetaFace::setEdgeOrientation( 1, rg_TRUE );
    else
		BetaFace::setEdgeOrientation( 1, rg_FALSE );

    if ( BetaFace::getEdge(0)->isThere( BetaFace::getEdge(2)->getEndVertex() ) == rg_TRUE )
		BetaFace::setEdgeOrientation( 2, rg_TRUE );
    else
		BetaFace::setEdgeOrientation( 2, rg_FALSE );

	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;

    m_indexOfExteriorRegionToLeftCell   = -1;
    m_indexOfExteriorRegionToRightCell  = -1;
}

AugmentedBCFace::AugmentedBCFace(const AugmentedBCFace& augBCFace)
: BetaFace( augBCFace )
{
	m_boundingState = augBCFace.m_boundingState;
	m_isArtificial  = augBCFace.m_isArtificial;

    m_indexOfExteriorRegionToLeftCell   = augBCFace.m_indexOfExteriorRegionToLeftCell;
    m_indexOfExteriorRegionToRightCell  = augBCFace.m_indexOfExteriorRegionToRightCell;
}

AugmentedBCFace::~AugmentedBCFace()
{
}




rg_INT	AugmentedBCFace::getBoundingState() const //CHECK! overloading잘 되는지 check
{
	return m_boundingState;
}

rg_BOOL AugmentedBCFace::isArtificial() const
{
	return m_isArtificial;
}


rg_INT AugmentedBCFace::getIndexOfExteriorRegionToLeftCell() const
{
    return m_indexOfExteriorRegionToLeftCell;
}

rg_INT AugmentedBCFace::getIndexOfExteriorRegionToRightCell() const
{
    return m_indexOfExteriorRegionToRightCell;
}





void AugmentedBCFace::setBoundingState(const rg_INT& boundingState)
{
	m_boundingState = boundingState;
}


void AugmentedBCFace::setIsArtificial( const rg_BOOL& isArtificial)
{
	m_isArtificial = isArtificial;
}


void AugmentedBCFace::setIndexOfExteriorRegionToLeftCell( const rg_INT& indexOfExteriorRegionToLeftCell )
{
    m_indexOfExteriorRegionToLeftCell = indexOfExteriorRegionToLeftCell;
}

void AugmentedBCFace::setIndexOfExteriorRegionToRightCell( const rg_INT& indexOfExteriorRegionToRightCell )
{
    m_indexOfExteriorRegionToRightCell = indexOfExteriorRegionToRightCell;
}





AugmentedBCFace& AugmentedBCFace::operator =(const AugmentedBCFace& augBCFace)
{
	if ( this == &augBCFace )
		return *this;

	BetaFace::operator =(augBCFace);

	m_boundingState = augBCFace.m_boundingState;
	m_isArtificial  = augBCFace.m_isArtificial;

    m_indexOfExteriorRegionToLeftCell   = augBCFace.m_indexOfExteriorRegionToLeftCell;
    m_indexOfExteriorRegionToRightCell  = augBCFace.m_indexOfExteriorRegionToRightCell;

	
	return *this;
}




void AugmentedBCFace::reverse()
{
    //  1. reverse cells incident to this face
    AugmentedBCCell* tempCell = (AugmentedBCCell*) getRightCell();
    setRightCell( getLeftCell() );
    setLeftCell(  (BetaCell*) tempCell );



    //  2. reverse edges
    for (rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++)
    {
        if( BetaFace::getEdgeOrientation(i) == rg_TRUE ) {
            BetaFace::setEdgeOrientation(i, rg_FALSE);
        }
        else { //if( BetaFace::getEdgeOrientation(i) == rg_FALSE )
            BetaFace::setEdgeOrientation(i, rg_TRUE);
        }
    }

    //      switch m_edge[1], m_edge[2]
    //      switch m_orientation1[1], m_orientation[2]    
    BetaEdge* tempEdge        = BetaFace::getEdge(1);
    rg_BOOL   tempOrientation = BetaFace::getEdgeOrientation(1);

    BetaFace::setEdge(1, BetaFace::getEdge(2));
    BetaFace::setEdgeOrientation(1, BetaFace::getEdgeOrientation(2));
    
	BetaFace::setEdge(2, tempEdge);
    BetaFace::setEdgeOrientation(2, tempOrientation);


    //  3. reverse index of exterior region
    rg_INT tempIndexOfExteriorRegionForLeft = m_indexOfExteriorRegionToLeftCell;
    
    m_indexOfExteriorRegionToLeftCell       = m_indexOfExteriorRegionToRightCell;
    m_indexOfExteriorRegionToRightCell      = tempIndexOfExteriorRegionForLeft;
}



/*
void AugmentedBCFace::reverseEdgeOrientations()
{
    for (rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++)
    {
        if( BetaFace::getEdgeOrientation(i) == rg_TRUE )
            BetaFace::setEdgeOrientation(i, rg_FALSE);
        else if( BetaFace::getEdgeOrientation(i) == rg_FALSE )
            BetaFace::setEdgeOrientation(i, rg_TRUE);
    }

    //switch m_edge[1], m_edge[2]
    //switch m_orientation1[1], m_orientation[2]    
    BetaEdge* tempEdge        = BetaFace::getEdge(1);
    rg_BOOL   tempOrientation = BetaFace::getEdgeOrientation(1);

    BetaFace::setEdge(1, BetaFace::getEdge(2));
    BetaFace::setEdgeOrientation(1, BetaFace::getEdgeOrientation(2));
    
	BetaFace::setEdge(2, tempEdge);
    BetaFace::setEdgeOrientation(2, tempOrientation);	
}
*/


AugmentedBCEdge* AugmentedBCFace::findMateEdgeOfVertex(const AugmentedBCVertex* const vertex) const
{
	for (rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++)
    {
        if ( BetaFace::getEdge(i)->getStartVertex() != (BetaVertex*)vertex
          && BetaFace::getEdge(i)->getEndVertex()   != (BetaVertex*)vertex )
        {
            return (AugmentedBCEdge*)BetaFace::getEdge(i);
        }
    }

    return rg_NULL;
}



AugmentedBCEdge* AugmentedBCFace::findIncidentEdgeOfVertexNotGivenEdge(const AugmentedBCVertex* const givenVertex, 
														  const AugmentedBCEdge* const notThisEdge) const
{
	for (rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++)
    {
        if ( BetaFace::getEdge(i)->getStartVertex() == (BetaVertex*)givenVertex
          || BetaFace::getEdge(i)->getEndVertex()   == (BetaVertex*)givenVertex )
        {
            if ( BetaFace::getEdge(i) != (BetaEdge*)notThisEdge )
                return (AugmentedBCEdge*)BetaFace::getEdge(i);
        }
    }

    return rg_NULL;
}

