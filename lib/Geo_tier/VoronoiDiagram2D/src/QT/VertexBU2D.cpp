#include "VertexBU2D.h"
#include "EdgeBU2D.h"
#include "FaceBU2D.h"
using namespace V::GeometryTier;



VertexBU2D::VertexBU2D()
: m_disc(rg_NULL)
{
}



VertexBU2D::VertexBU2D(const rg_INT& ID, Disc* disc)
: VertexWEDS(ID), m_disc(disc)
{
}



VertexBU2D::VertexBU2D(const VertexBU2D& vertex)
: VertexWEDS(vertex), 
  m_betaSpan(vertex.m_betaSpan),
  m_disc(vertex.m_disc)
{
}



VertexBU2D::~VertexBU2D()
{
}




Disc*    VertexBU2D::getDisc()
{
    return m_disc;
}


rg_Point2D  VertexBU2D::getCoord() const
{
    return m_disc->getGeometry().getCenterPt();
}



rg_Circle2D VertexBU2D::getCircle() const
{
    return m_disc->getGeometry();
}



rg_BetaSpan VertexBU2D::getBetaSpan() const
{
    return m_betaSpan;
}




void     VertexBU2D::setDisc(Disc* disc)
{
    m_disc = disc;
}



rg_BOOL VertexBU2D::isOnSqueezedHull() const
{
    rg_dList<VertexBU2D*> neighborVertices;
	getNeighborVertices( neighborVertices );

    neighborVertices.reset4Loop();
    while ( neighborVertices.setNext4Loop() ) {
        VertexBU2D* currVtx = neighborVertices.getEntity();
        if ( currVtx->isVirtual() ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



void VertexBU2D::computeBetaSpan()
{
    if ( isVirtual() ) {
        return;
    }

    rg_dList<EdgeBU2D*> incidentEdges;
	getIncidentEdges( incidentEdges );

    rg_REAL minOfHigherSimplex = rg_REAL_INFINITY;
    rg_REAL maxOfHigherSimplex = -rg_REAL_INFINITY;

    incidentEdges.reset4Loop();
    while ( incidentEdges.setNext4Loop() ) {
        EdgeBU2D* currEdge = incidentEdges.getEntity();

        if ( currEdge->isVirtual() ) {
            continue;
        }

        Interval_REAL interval = currEdge->getBetaSpan().getBetaIntervalOfSingularRegularState();

        rg_REAL minValue = interval.getLowerValue();
        rg_REAL maxValue = interval.getUpperValue();
        if ( minValue < minOfHigherSimplex ) {
            minOfHigherSimplex = minValue;
        }

        if ( maxValue > maxOfHigherSimplex ) {
            maxOfHigherSimplex = maxValue;
        }                

    }


    if ( isOnSqueezedHull() ) {
        m_betaSpan.setNumBetaInterVal(2);
        m_betaSpan.setBetaInterval(0, SINGULAR_SIMPLEX, -rg_REAL_INFINITY,  minOfHigherSimplex);
        m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,  minOfHigherSimplex, rg_REAL_INFINITY);
    }
    else {
        m_betaSpan.setNumBetaInterVal(3);
        m_betaSpan.setBetaInterval(0, SINGULAR_SIMPLEX, -rg_REAL_INFINITY,  minOfHigherSimplex);
        m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,  minOfHigherSimplex, maxOfHigherSimplex);
        m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX, maxOfHigherSimplex, rg_REAL_INFINITY);
    }
}



VertexBU2D& VertexBU2D::operator =(const VertexBU2D& vertex)
{
    if ( this == &vertex ) {
        return *this;
    }

    VertexWEDS::operator =(vertex);
    m_betaSpan = vertex.m_betaSpan;
    m_disc     = vertex.m_disc;

    return *this;
}



///////////////////////////////////////////////////////////////////////////
//  Topological operators
rg_INT    VertexBU2D::getNeighborVertices(  rg_dList<VertexBU2D*>& vertexList) const
{
    rg_dList<VertexWEDS*> vertexListInWEDS;
    VertexWEDS::getNeighborVertices( vertexListInWEDS );
    
    vertexListInWEDS.reset4Loop();
    while ( vertexListInWEDS.setNext4Loop() ) {
        VertexWEDS* currVtx = vertexListInWEDS.getEntity();
        vertexList.add( (VertexBU2D*)currVtx );
    }

    return vertexList.getSize();
}



rg_INT    VertexBU2D::getIncidentEdges(     rg_dList<EdgeBU2D*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    VertexWEDS::getIncidentEdges( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgeListInWEDS.getEntity();
        edgeList.add( (EdgeBU2D*)currEdge );
    }

    return edgeList.getSize();
}



rg_INT    VertexBU2D::getIncidentFaces(     rg_dList<FaceBU2D*>& faceList) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    VertexWEDS::getIncidentFaces( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceWEDS* currFace = faceListInWEDS.getEntity();
        faceList.add( (FaceBU2D*)currFace );
    }

    return faceList.getSize();
}




rg_INT    VertexBU2D::getEdgesInStar(       rg_dList<EdgeBU2D*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    VertexWEDS::getEdgesInStar( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgeListInWEDS.getEntity();
        edgeList.add( (EdgeBU2D*)currEdge );
    }

    return edgeList.getSize();
}



rg_INT    VertexBU2D::getEdgesInShell(      rg_dList<EdgeBU2D*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    VertexWEDS::getEdgesInShell( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgeListInWEDS.getEntity();
        edgeList.add( (EdgeBU2D*)currEdge );
    }

    return edgeList.getSize();
}




EdgeBU2D* VertexBU2D::findConnectingEdge(   VertexBU2D* vertex) const
{
    EdgeWEDS* edge = VertexWEDS::findConnectingEdge( vertex );
    return (EdgeBU2D*)edge;

}



rg_INT    VertexBU2D::findSharingFace(      VertexBU2D* vertex, rg_dList<FaceBU2D*>& faceList) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    VertexWEDS::findSharingFace( vertex, faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceWEDS* currFace = faceListInWEDS.getEntity();
        faceList.add( (FaceBU2D*)currFace );
    }

    return faceList.getSize();
}



///////////////////////////////////////////////////////////////////////////
//  Beta-operator
rg_INT VertexBU2D::searchVerticesInBetaComplex(const rg_REAL& beta, rg_dList<VertexBU2D*>& betaVertexList ) const
{
    rg_dList<VertexWEDS*> vertexListInWEDS;
    VertexWEDS::getNeighborVertices( vertexListInWEDS );
    
    vertexListInWEDS.reset4Loop();
    while ( vertexListInWEDS.setNext4Loop() ) {
        VertexBU2D* currVtx = (VertexBU2D*)(vertexListInWEDS.getEntity());

        if ( currVtx->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaVertexList.add( currVtx );
        }
    }

    return betaVertexList.getSize();

}



rg_INT VertexBU2D::searchEdgesInBetaComplex(   const rg_REAL& beta, rg_dList<EdgeBU2D*>&   betaEdgeList   ) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    VertexWEDS::getIncidentEdges( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeBU2D* currEdge = (EdgeBU2D*)(edgeListInWEDS.getEntity());

        if ( currEdge->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaEdgeList.add( currEdge );
        }
    }

    return betaEdgeList.getSize();

}



rg_INT VertexBU2D::searchFacesInBetaComplex(   const rg_REAL& beta, rg_dList<FaceBU2D*>&   betaFaceList   ) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    VertexWEDS::getIncidentFaces( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceBU2D* currFace = (FaceBU2D*)(faceListInWEDS.getEntity());

        if ( currFace->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaFaceList.add( currFace );
        }
    }

    return betaFaceList.getSize();

}




rg_INT VertexBU2D::searchVerticesInBetaComplex(const rg_REAL& beta, const rg_INT& boundingState, rg_dList<VertexBU2D*>& betaVertexList ) const
{
    rg_dList<VertexWEDS*> vertexListInWEDS;
    VertexWEDS::getNeighborVertices( vertexListInWEDS );
    
    vertexListInWEDS.reset4Loop();
    while ( vertexListInWEDS.setNext4Loop() ) {
        VertexBU2D* currVtx = (VertexBU2D*)(vertexListInWEDS.getEntity());

        if ( currVtx->getBoundingState(beta) == boundingState ) {
            betaVertexList.add( currVtx );
        }
    }

    return betaVertexList.getSize();

}



rg_INT VertexBU2D::searchEdgesInBetaComplex(   const rg_REAL& beta, const rg_INT& boundingState, rg_dList<EdgeBU2D*>& betaEdgeList ) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    VertexWEDS::getIncidentEdges( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeBU2D* currEdge = (EdgeBU2D*)(edgeListInWEDS.getEntity());

        if ( currEdge->getBoundingState(beta) == boundingState ) {
            betaEdgeList.add( currEdge );
        }
    }

    return betaEdgeList.getSize();
}



rg_INT VertexBU2D::searchFacesInBetaComplex(   const rg_REAL& beta, const rg_INT& boundingState, rg_dList<FaceBU2D*>& betaFaceList ) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    VertexWEDS::getIncidentFaces( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceBU2D* currFace = (FaceBU2D*)(faceListInWEDS.getEntity());

        if ( currFace->getBoundingState(beta) == boundingState ) {
            betaFaceList.add( currFace );
        }
    }

    return betaFaceList.getSize();
}




rg_INT VertexBU2D::searchVerticesInBetaShape( const rg_REAL& beta, rg_dList<VertexBU2D*>& betaVertexList ) const
{
    rg_dList<VertexWEDS*> vertexListInWEDS;
    VertexWEDS::getNeighborVertices( vertexListInWEDS );
    
    vertexListInWEDS.reset4Loop();
    while ( vertexListInWEDS.setNext4Loop() ) {
        VertexBU2D* currVtx = (VertexBU2D*)(vertexListInWEDS.getEntity());

        rg_INT state = currVtx->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaVertexList.add( currVtx );
        }

    }

    return betaVertexList.getSize();
}



rg_INT VertexBU2D::searchEdgesInBetaShape(    const rg_REAL& beta, rg_dList<EdgeBU2D*>&   betaEdgeList   ) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    VertexWEDS::getIncidentEdges( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeBU2D* currEdge = (EdgeBU2D*)(edgeListInWEDS.getEntity());

        rg_INT state = currEdge->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaEdgeList.add( currEdge );
        }
    }

    return betaEdgeList.getSize();
}





