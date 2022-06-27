#include "EdgeBU2D.h"
#include "VertexBU2D.h"
#include "FaceBU2D.h"
using namespace V::GeometryTier;



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




rg_BetaSpan EdgeBU2D::getBetaSpan() const
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



void EdgeBU2D::computeBetaSpan_horizontal_or_vertical_edge()
{
    double r_minTangentCircle = (   1.0 
        - ((VertexBU2D*)m_startVertex)->getCircle().getRadius() 
        - ((VertexBU2D*)m_endVertex)->getCircle().getRadius() ) / 2.0;
    double minRadius_higherSimplex = 0.0;
    double maxRadius_higherSimplex = 0.0;

    FaceBU2D* leftFace  = (FaceBU2D*)m_leftFace;
    FaceBU2D* rightFace = (FaceBU2D*)m_rightFace;

    if ( isOnSqueezedHull() ) {
        if ( leftFace->isVirtual() ) {
            minRadius_higherSimplex = rightFace->getEmptyTangentCircle().getRadius();
        }
        else {
            minRadius_higherSimplex = leftFace->getEmptyTangentCircle().getRadius();
        }

        m_betaSpan.setNumBetaInterVal( 3 );
        m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,        r_minTangentCircle);
        m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   r_minTangentCircle,       minRadius_higherSimplex);
        m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minRadius_higherSimplex,  rg_REAL_INFINITY);
    }
    else {
        double r_leftFace   = leftFace->getEmptyTangentCircle().getRadius();
        double r_rightFace  = rightFace->getEmptyTangentCircle().getRadius();

        if ( r_leftFace < r_rightFace ) {
            minRadius_higherSimplex = r_leftFace;
            maxRadius_higherSimplex = r_rightFace;
        }
        else {
            minRadius_higherSimplex = r_rightFace;
            maxRadius_higherSimplex = r_leftFace;
        }

        m_betaSpan.setNumBetaInterVal( 4 );
        m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,        r_minTangentCircle);
        m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   r_minTangentCircle,       minRadius_higherSimplex);
        m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minRadius_higherSimplex,  maxRadius_higherSimplex);
        m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxRadius_higherSimplex,  rg_REAL_INFINITY);
    }
}



void EdgeBU2D::computeBetaSpan_diagonal_edge()
{
    const double SQUARE_ROOT_2             = 1.4142135623730950488016887242097;
    const double HALF_SQUARE_ROOT_2        = 0.70710678118654752440084436210485;

    bool    b_edgeCanBeAttached         = false;
    bool    b_edgeIsAttached            = false;

    double  r_startVertex               = ((VertexBU2D*)m_startVertex)->getCircle().getRadius();
    double  r_endVertex                 = ((VertexBU2D*)m_endVertex)->getCircle().getRadius();
    double  r_minTangentCircle          = ( SQUARE_ROOT_2 - r_startVertex - r_endVertex ) / 2.0;

    double  r_upperVertex               = ((VertexBU2D*)m_leftLeg->getEndVertex())->getCircle().getRadius();
    double  r_lowerVertex               = ((VertexBU2D*)m_rightHand->getStartVertex())->getCircle().getRadius();

    double  r_leftFace                  = ((FaceBU2D*)m_leftFace)->getEmptyTangentCircle().getRadius();
    double  r_rightFace                 = ((FaceBU2D*)m_rightFace)->getEmptyTangentCircle().getRadius();

    double  minRadius_higherSimplex, maxRadius_higherSimplex;
    if ( r_leftFace < r_rightFace ) {
        minRadius_higherSimplex = r_leftFace;
        maxRadius_higherSimplex = r_rightFace;
    }
    else {
        minRadius_higherSimplex = r_rightFace;
        maxRadius_higherSimplex = r_leftFace;
    }

    /*if ( r_minTangentCircle + r_upperVertex > HALF_SQUARE_ROOT_2 || r_minTangentCircle + r_lowerVertex > HALF_SQUARE_ROOT_2 ) {
        b_edgeCanBeAttached = true;
    }

    if ( b_edgeCanBeAttached ) {*/
        double  x_coord_minTangent          = ( r_minTangentCircle + r_startVertex ) * HALF_SQUARE_ROOT_2; // y_coord_minTangent == x_coord_minTangent
        double  x_coord_minTangent_square    = x_coord_minTangent * x_coord_minTangent;
        double  squaredDistance             = (2.0 * (x_coord_minTangent_square - x_coord_minTangent) + 1.0 );

        double  squaredSumRadii_upper       = r_minTangentCircle + r_upperVertex;
        squaredSumRadii_upper = squaredSumRadii_upper * squaredSumRadii_upper;
        double  squaredSumRadii_lower       = r_minTangentCircle + r_lowerVertex;
        squaredSumRadii_lower = squaredSumRadii_lower * squaredSumRadii_lower;

        if ( squaredDistance < squaredSumRadii_upper || squaredDistance < squaredSumRadii_lower ) {
            b_edgeIsAttached = true;
        }
    //}

    if ( b_edgeIsAttached ) {
        m_betaSpan.setNumBetaInterVal( 3 );
        m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,        minRadius_higherSimplex);
        m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minRadius_higherSimplex,  maxRadius_higherSimplex);
        m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxRadius_higherSimplex,  rg_REAL_INFINITY);
    }
    else {
        m_betaSpan.setNumBetaInterVal( 4 );
        m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY,        r_minTangentCircle);
        m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   r_minTangentCircle,       minRadius_higherSimplex);
        m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minRadius_higherSimplex,  maxRadius_higherSimplex);
        m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxRadius_higherSimplex,  rg_REAL_INFINITY);
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
    if ( vertex == m_startVertex ) {
        return (VertexBU2D*)m_endVertex;
    }
    else if ( vertex == m_endVertex ) {
        return (VertexBU2D*)m_startVertex;
    }
    else {
        return NULL;
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





