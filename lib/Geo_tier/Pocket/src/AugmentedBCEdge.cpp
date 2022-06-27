#include "AugmentedBCEdge.h"

#include "AugmentedBCFace.h"
#include "AugmentedBCCell.h"

AugmentedBCEdge::AugmentedBCEdge()
: BetaEdge()  
{
	m_boundingState = rg_UNKNOWN;	//EXTRANEOUS_SIMPLEX와 같으므로 주의.
	m_isArtificial  = rg_FALSE;

    m_indexOfExteriorRegion = -1;
}

AugmentedBCEdge::AugmentedBCEdge(const rg_INT& id)
: BetaEdge()
{
	BetaEdge::setID( id );

	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;

    m_indexOfExteriorRegion = -1;
}


AugmentedBCEdge::AugmentedBCEdge(AugmentedBCVertex* startVertex, AugmentedBCVertex* endVertex)
: BetaEdge((BetaVertex*)startVertex, (BetaVertex*)endVertex)
{
	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;

    m_indexOfExteriorRegion = -1;
}


AugmentedBCEdge::AugmentedBCEdge(const AugmentedBCEdge& augBCEdge)
: BetaEdge(augBCEdge)
{
	m_boundingState = augBCEdge.m_boundingState;
	m_isArtificial  = augBCEdge.m_isArtificial;

    m_indexOfExteriorRegion = augBCEdge.m_indexOfExteriorRegion;
}

AugmentedBCEdge::~AugmentedBCEdge()
{
}





rg_INT	AugmentedBCEdge::getBoundingState() const //CHECK! overloading잘 되는지 check
{
	return m_boundingState;
}

rg_BOOL AugmentedBCEdge::isArtificial() const
{
	return m_isArtificial;
}

rg_INT AugmentedBCEdge::getIndexOfExteriorRegion() const
{
    return m_indexOfExteriorRegion;    
}






void AugmentedBCEdge::setBoundingState(const rg_INT& boundingState)
{
	m_boundingState = boundingState;
}


void AugmentedBCEdge::setIsArtificial( const rg_BOOL& isArtificial)
{
	m_isArtificial = isArtificial;
}


void AugmentedBCEdge::setIndexOfExteriorRegion( const rg_INT& indexOfExteriorRegion )
{
    m_indexOfExteriorRegion = indexOfExteriorRegion;
}




AugmentedBCEdge& AugmentedBCEdge::operator =(const AugmentedBCEdge& augBCEdge)
{
	if ( this == &augBCEdge )
		return *this;

	BetaEdge::operator =(augBCEdge);

	m_boundingState = augBCEdge.m_boundingState;
	m_isArtificial  = augBCEdge.m_isArtificial;

    m_indexOfExteriorRegion = augBCEdge.m_indexOfExteriorRegion;

	
	return *this;
}



AugmentedBCVertex* AugmentedBCEdge::getBoundingVertexNotGivenVertex(const AugmentedBCVertex* const notThisVertex) const
{
	BetaVertex* returnVertex = rg_NULL;
    if ( BetaEdge::getStartVertex() == (BetaVertex*)notThisVertex )
        returnVertex = BetaEdge::getEndVertex();
    else
        returnVertex = BetaEdge::getStartVertex();

    return (AugmentedBCVertex*)returnVertex;
}





AugmentedBCCell* AugmentedBCEdge::getCCWNextCell(const AugmentedBCFace* const currFace) const
{
    if ( currFace->isThere( (BetaEdge*)this ) ) {
        BetaCell* ccwNextCell = rg_NULL;
        if ( currFace->getEdgeOrientation( (BetaEdge*)this ) ) {
            ccwNextCell = currFace->getRightCell();
        }
        else {
            ccwNextCell = currFace->getLeftCell();
        }

        return (AugmentedBCCell*)ccwNextCell;
    }
    else {
        return rg_NULL;
    }
}



AugmentedBCFace* AugmentedBCEdge::getCCWNextFace(const AugmentedBCFace* const currFace) const
{
    AugmentedBCCell* ccwNextCell = getCCWNextCell( currFace );

    if ( ccwNextCell != rg_NULL ) {
        AugmentedBCFace* incidentFace[2];

		rg_INT index4Face = 0;
		for ( rg_INT i=0; i<4; i++ )  {
			AugmentedBCFace* currBoundingFace = (AugmentedBCFace*)ccwNextCell->getFace(i);
			if ( currBoundingFace->isThere((BetaEdge*)this) == rg_TRUE )  {
				incidentFace[index4Face++] = currBoundingFace;
			}
		}

        //ccwNextCell->AugmentedBCCell::findFacesIncidentToEdge( this, incidentFace);
     
        AugmentedBCFace* ccwNextFace = rg_NULL;
        if ( currFace == incidentFace[0] ) {
            ccwNextFace = incidentFace[1];
        }
        else {
            ccwNextFace = incidentFace[0];
        }

		/*
		//added 2010-02-08, for elliptic edge in VD
		if ( ccwNextFace->getLeftCell() == rg_NULL
		  && ccwNextFace->getRightCell() == rg_NULL )
		  return rg_NULL;
		  */

        return ccwNextFace;
    }
    else {
        return rg_NULL;
    }
}



AugmentedBCCell* AugmentedBCEdge::getCWNextCell(const AugmentedBCFace* const currFace) const
{
    if ( currFace->isThere( (BetaEdge*)this ) ) {
        BetaCell* cwNextCell = rg_NULL;
        if ( currFace->getEdgeOrientation( (BetaEdge*)this ) ) {
            cwNextCell = currFace->getLeftCell();
        }
        else {
            cwNextCell = currFace->getRightCell();
        }

        return (AugmentedBCCell*)cwNextCell;
    }
    else {
        return rg_NULL;
    }
}



AugmentedBCFace* AugmentedBCEdge::getCWNextFace(const AugmentedBCFace* const currFace) const
{
    AugmentedBCCell* cwNextCell = getCWNextCell( currFace );

    if ( cwNextCell != rg_NULL ) {
        AugmentedBCFace* incidentFace[2];

		rg_INT index4Face = 0;
		for ( rg_INT i=0; i<4; i++ )  {
			AugmentedBCFace* currBoundingFace = (AugmentedBCFace*)cwNextCell->getFace(i);
			if ( currBoundingFace->isThere((BetaEdge*)this) == rg_TRUE )  {
				incidentFace[index4Face++] = currBoundingFace;
			}
		}

        //cwNextCell->AugmentedBCCell::findFacesIncidentToEdge( this, incidentFace);
     
        AugmentedBCFace* cwNextFace = rg_NULL;
        if ( currFace == incidentFace[0] ) {
            cwNextFace = incidentFace[1];
        }
        else {
            cwNextFace = incidentFace[0];
        }

		/*
		//added 2010-02-08, for elliptic edge in VD
		if ( cwNextFace->getLeftCell() == rg_NULL
		  && cwNextFace->getRightCell() == rg_NULL )
		  return rg_NULL;
		  */

        return cwNextFace;
    }
    else {
        return rg_NULL;
    }
}

