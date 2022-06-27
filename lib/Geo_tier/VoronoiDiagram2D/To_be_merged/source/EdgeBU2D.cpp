#include "EdgeBU2D.h"
#include "VertexBU2D.h"
#include "FaceBU2D.h"
using namespace BULL2D::GeometryTier;


EdgeBU2D::EdgeBU2D()
{
}



EdgeBU2D::EdgeBU2D(const rg_INT& ID)
: EdgeWEDS(ID) 
{
}



EdgeBU2D::EdgeBU2D(const EdgeBU2D& edge)
: EdgeWEDS(edge),
  m_betaSpan(edge.m_betaSpan)
{
}



EdgeBU2D::~EdgeBU2D()
{
}




BetaSpan EdgeBU2D::getBetaSpan() const
{
    return m_betaSpan;
}


rg_BOOL EdgeBU2D::isVirtual() const
{
    if ( getStartVertex()->isVirtual() || getEndVertex()->isVirtual() ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}

	
rg_BOOL EdgeBU2D::isOnSqueezedHull() const
{
    if(getLeftFace()->isVirtual() || getRightFace()->isVirtual() ) {
		return rg_TRUE;
    }
    else {
		return rg_FALSE;
    }
}



rg_BOOL EdgeBU2D::isAttached(const rg_Circle2D& minTangentCircle) const
{
    if ( isVirtual() ) {
        return rg_FALSE;
    }

    // check for whether two q-vertices of this q-edge 
    // define only a single q-edge in the whole quasi-triangulation or not.
    VertexBU2D* startVtx = getStartVertex();
    VertexBU2D* endVtx   = getEndVertex();
    rg_dList<VertexBU2D*> neighborVertexOfStartVertex;
    startVtx->getNeighborVertices( neighborVertexOfStartVertex );

    rg_INT numEdgesWithSameVertices = 0;
    neighborVertexOfStartVertex.reset4Loop();
    while ( neighborVertexOfStartVertex.setNext4Loop() ) {
        VertexBU2D* neighborVertex = neighborVertexOfStartVertex.getEntity();

        if ( endVtx == neighborVertex ) {
            numEdgesWithSameVertices++;
        }
    }


    //  When numEdgesWithSameVertices == 1, 
    //  two q-vertices of this q-edge 
    //  define only a single q-edge, this edge, in the whole quasi-triangulation.
    //  Therefore, the emptiness test of minimum tangent circle($\xi$) is needed only.
    if ( numEdgesWithSameVertices >= 1 ) {
        //  check the emptiness of minimum tangent circle(minTangentCircle).
        rg_dList<VertexBU2D*> neighborVtxs;
        getNeighborVertices( neighborVtxs );

        neighborVtxs.reset4Loop();
        while ( neighborVtxs.setNext4Loop() ) {
            VertexBU2D* currVtx = neighborVtxs.getEntity();

            if ( currVtx->isVirtual() ) {
                continue;
            }

            rg_Circle2D circle = currVtx->getCircle();
            if ( minTangentCircle.isIntersectWith( circle ) ) {
                return rg_TRUE;
            }
        }
    }

    //  When numEdgesWithSameVertices > 1, 
    //  two q-vertices of this q-edge define distinct multiple q-edges in the whole quasi-triangulation.
    //  In the primary structure (VD), this means two atoms define multiple v-faces.
    if ( numEdgesWithSameVertices > 1 ) {
        //  check signed area of incident faces.
        FaceBU2D* face[2] = {getLeftFace(), getRightFace()};
        for ( rg_INT i=0; i<2; i++ ) {
            if ( face[i]->isVirtual() ) {
                continue;
            }

            rg_REAL signedArea  = face[i]->computeSignedArea();
            if( signedArea < 0.0 ) {
                return rg_TRUE;
            }
        }
    }

    return rg_FALSE;

    //if ( isVirtual() ) {
    //    return rg_FALSE;
    //}

    //rg_BOOL isAttachedEdge = rg_FALSE;

    ////  intersection check.
    //rg_dList<VertexBU2D*> neighborVtxs;
    //getNeighborVertices( neighborVtxs );

    //neighborVtxs.reset4Loop();
    //while ( neighborVtxs.setNext4Loop() ) {
    //    VertexBU2D* currVtx = neighborVtxs.getEntity();

    //    if ( currVtx->isVirtual() ) {
    //        continue;
    //    }

    //    rg_Circle2D circle = currVtx->getCircle();
    //    if ( minTangentCircle.isIntersectWith( circle ) ) {
    //        isAttachedEdge = rg_TRUE;
    //        break;
    //    }
    //}

    //if ( isAttachedEdge ) {
    //    return isAttachedEdge;
    //}


    ////  signed area check of incident faces.
    //FaceBU2D* face[2] = {getLeftFace(), getRightFace()};
    //for ( rg_INT i=0; i<2; i++ ) {
    //    if ( face[i]->isVirtual() ) {
    //        continue;
    //    }

    //    rg_REAL signedArea  = face[i]->computeSignedArea();

    //    if( signedArea < 0.0 ) {
    //        isAttachedEdge = rg_TRUE;
    //        break;
    //    }
    //}

    //return isAttachedEdge;
}



void EdgeBU2D::computeBetaSpan()
{
    if ( isVirtual() ) {
        return;
    }

    rg_REAL minOfHigherSimplex = rg_REAL_INFINITY;
    rg_REAL maxOfHigherSimplex = -rg_REAL_INFINITY;

    FaceBU2D* face[2] = {getLeftFace(), getRightFace()};
    for ( rg_INT i=0; i<2; i++ ) {
        if ( face[i]->isVirtual() ) {
            continue;
        }

        rg_REAL radius = face[i]->getEmptyTangentCircle().getRadius();

        if( radius < minOfHigherSimplex ) {
			minOfHigherSimplex = radius;
        }
        if( radius > maxOfHigherSimplex ) {
			maxOfHigherSimplex = radius;
        }
    }

    rg_Circle2D circleOnStart    = getStartVertex()->getCircle();
	rg_Circle2D circleOnEnd      = getEndVertex()->getCircle();
	rg_Circle2D minTangentCircle = circleOnStart.computeMinTangentCircle(circleOnEnd);

	rg_REAL minOfCurrSimplex     = minTangentCircle.getRadius();


	if( isOnSqueezedHull() ) {
		if( isAttached(minTangentCircle) ) {
			m_betaSpan.setNumBetaInterVal( 2 );
			m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,  minOfHigherSimplex);
			m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex, rg_REAL_INFINITY);
        }
        else {
			m_betaSpan.setNumBetaInterVal( 3 );
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,  minOfCurrSimplex);
            m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,   minOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex, rg_REAL_INFINITY);			
        }
    }
    else {
		if( isAttached(minTangentCircle) ) {
			m_betaSpan.setNumBetaInterVal( 3 );
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,  minOfHigherSimplex);
            m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex, maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxOfHigherSimplex, rg_REAL_INFINITY);						
        }
        else {
			m_betaSpan.setNumBetaInterVal( 4 );
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,  minOfCurrSimplex);
            m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,   minOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex, maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfHigherSimplex, rg_REAL_INFINITY);			
        }
    }



}





EdgeBU2D& EdgeBU2D::operator =(const EdgeBU2D& edge)
{
    if ( this == &edge ) {
        return *this;
    }

    EdgeWEDS::operator =(edge);
    m_betaSpan = edge.m_betaSpan;

    return *this;
}



///////////////////////////////////////////////////////////////////////////
//  Topological operators
VertexBU2D* EdgeBU2D::getVertexOfLeftHand() const
{
    return (VertexBU2D*)(EdgeWEDS::getVertexOfLeftHand());
}



VertexBU2D* EdgeBU2D::getVertexOfRightHand() const
{
    return (VertexBU2D*)(EdgeWEDS::getVertexOfRightHand());
}



VertexBU2D* EdgeBU2D::getOppositeVertex( const VertexBU2D* const vertex ) const
{
    if( vertex == getStartVertex() )
    {
        return getEndVertex();
    }
    else if( vertex == getEndVertex() )
    {
        return getStartVertex();
    }
    else
    {
        return rg_NULL;
    }
}


    
rg_INT EdgeBU2D::getNeighborVertices(rg_dList<VertexBU2D*>& vertexList) const
{
    rg_dList<VertexWEDS*> vertexListInWEDS;
    EdgeWEDS::getNeighborVertices( vertexListInWEDS );
    
    vertexListInWEDS.reset4Loop();
    while ( vertexListInWEDS.setNext4Loop() ) {
        VertexWEDS* currVtx = vertexListInWEDS.getEntity();
        vertexList.add( (VertexBU2D*)currVtx );
    }

    return vertexList.getSize();
}



rg_INT EdgeBU2D::getIncidentEdges(rg_dList<EdgeBU2D*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    EdgeWEDS::getIncidentEdges( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgeListInWEDS.getEntity();
        edgeList.add( (EdgeBU2D*)currEdge );
    }

    return edgeList.getSize();
}



rg_INT EdgeBU2D::getIncidentFaces(rg_dList<FaceBU2D*>& faceList) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    EdgeWEDS::getIncidentFaces( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceWEDS* currFace = faceListInWEDS.getEntity();
        faceList.add( (FaceBU2D*)currFace );
    }

    return faceList.getSize();
}




rg_INT EdgeBU2D::getEdgesInStar(rg_dList<EdgeBU2D*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    EdgeWEDS::getEdgesInStar( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgeListInWEDS.getEntity();
        edgeList.add( (EdgeBU2D*)currEdge );
    }

    return edgeList.getSize();
}



rg_INT EdgeBU2D::getFacesInStar(rg_dList<FaceBU2D*>& faceList) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    EdgeWEDS::getFacesInStar( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceWEDS* currFace = faceListInWEDS.getEntity();
        faceList.add( (FaceBU2D*)currFace );
    }

    return faceList.getSize();
}



///////////////////////////////////////////////////////////////////////////
//  Beta-operator
rg_INT EdgeBU2D::searchEdgesInBetaComplex( const rg_REAL& beta, rg_dList<EdgeBU2D*>& betaEdgeList ) const
{
    EdgeBU2D* neighborEdge[4] = { getLeftHand(), getRightHand(), getLeftLeg(), getRightLeg() };

    for ( rg_INT i=0; i<4; i++ ) {
        if ( neighborEdge[i]->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaEdgeList.add( neighborEdge[i] );
        }
    }

    return betaEdgeList.getSize();
}



rg_INT EdgeBU2D::searchFacesInBetaComplex( const rg_REAL& beta, rg_dList<FaceBU2D*>& betaFaceList ) const
{
    FaceBU2D* incidentFace[2] = { getLeftFace(), getRightFace() };

    for ( rg_INT i=0; i<2; i++ ) {
        if ( incidentFace[i]->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaFaceList.add( incidentFace[i] );
        }
    }

    return betaFaceList.getSize();
}




rg_INT EdgeBU2D::searchEdgesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<EdgeBU2D*>& betaEdgeList ) const
{
    EdgeBU2D* neighborEdge[4] = { getLeftHand(), getRightHand(), getLeftLeg(), getRightLeg() };

    for ( rg_INT i=0; i<4; i++ ) {
        if ( neighborEdge[i]->getBoundingState(beta) == boundingState ) {
            betaEdgeList.add( neighborEdge[i] );
        }
    }

    return betaEdgeList.getSize();
}



rg_INT EdgeBU2D::searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<FaceBU2D*>& betaFaceList ) const
{
    FaceBU2D* incidentFace[2] = { getLeftFace(), getRightFace() };

    for ( rg_INT i=0; i<2; i++ ) {
        if ( incidentFace[i]->getBoundingState(beta) == boundingState ) {
            betaFaceList.add( incidentFace[i] );
        }
    }

    return betaFaceList.getSize();
}




rg_INT EdgeBU2D::searchEdgesInBetaShape(   const rg_REAL& beta, rg_dList<EdgeBU2D*>& betaEdgeList ) const
{
    EdgeBU2D* neighborEdge[4] = { getLeftHand(), getRightHand(), getLeftLeg(), getRightLeg() };

    for ( rg_INT i=0; i<4; i++ ) {
        rg_INT state = neighborEdge[i]->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaEdgeList.add( neighborEdge[i] );
        }
    }

    return betaEdgeList.getSize();
}



rg_BOOL EdgeBU2D::isMemberOf( const FaceBU2D* const face ) const
{
    if( face == getLeftFace() || face == getRightFace() )
    {
        return true;
    }
    else
    {
        return false;
    }
}



rg_BOOL EdgeBU2D::isAdjacentTo( const EdgeBU2D* const edge ) const
{
    if( edge == getLeftHand() || edge == getRightHand() || edge == getLeftLeg() || edge == getRightLeg() )
    {
        return true;
    }
    else
    {
        return false;
    }
}



rg_BOOL EdgeBU2D::isIncidentTo( const VertexBU2D* const vertex ) const
{
    if( vertex == getStartVertex() || vertex == getEndVertex() )
    {
        return true;
    }
    else
    {
        return false;
    }
}





