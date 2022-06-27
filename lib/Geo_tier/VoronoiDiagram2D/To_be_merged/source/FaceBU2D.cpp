#include "FaceBU2D.h"
#include "EdgeBU2D.h"
#include "VertexBU2D.h"
#include "rg_GeoFunc.h"

using namespace BULL2D::GeometryTier;



FaceBU2D::FaceBU2D()
{
}



FaceBU2D::FaceBU2D(const rg_INT& ID, const rg_Circle2D& emptyTangentCircle)
: FaceWEDS(ID), m_emptyTangentCircle(emptyTangentCircle)
{
}



FaceBU2D::FaceBU2D(const FaceBU2D& face)
: FaceWEDS(face), 
  m_betaSpan(face.m_betaSpan),
  m_emptyTangentCircle(face.m_emptyTangentCircle)
{
}



FaceBU2D::~FaceBU2D()
{
}




BetaSpan    FaceBU2D::getBetaSpan() const
{
    return m_betaSpan;
}



rg_Circle2D FaceBU2D::getEmptyTangentCircle() const
{
    return m_emptyTangentCircle;
}



rg_BOOL FaceBU2D::isVirtual() const
{
    rg_BOOL isVirtualFace = rg_FALSE;
    if ( getFirstEdge()->isVirtual() ) {
        isVirtualFace = rg_TRUE;
    }
    else {
        if ( this == getFirstEdge()->getLeftFace() ) {
            if ( getFirstEdge()->getLeftHand()->isVirtual() ) {
                isVirtualFace = rg_TRUE;
            }
        }
        else {
            if ( getFirstEdge()->getRightHand()->isVirtual() ) {
                isVirtualFace = rg_TRUE;
            }
        }
    }

    return isVirtualFace;
}



void        FaceBU2D::setEmptyTangentCircle(const rg_Circle2D& emptyTangentCircle)
{
    m_emptyTangentCircle = emptyTangentCircle;
}


	
void FaceBU2D::computeBetaSpan()
{
    if ( isVirtual() ) {
        return;
    }

    rg_REAL minOfCurrSimplex = m_emptyTangentCircle.getRadius();

    m_betaSpan.setNumBetaInterVal(2);
    m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, -rg_REAL_INFINITY, minOfCurrSimplex);
    m_betaSpan.setBetaInterval(1, INTERIOR_SIMPLEX,   minOfCurrSimplex,  rg_REAL_INFINITY);
}



FaceBU2D& FaceBU2D::operator =(const FaceBU2D& face)
{
    if ( this == &face ) {
        return *this;
    }

    FaceWEDS::operator =(face);
    m_betaSpan           = face.m_betaSpan;
    m_emptyTangentCircle = face.m_emptyTangentCircle;

    return *this;
}


    
rg_REAL FaceBU2D::computeSignedArea() const
{
    rg_REAL signedArea = 0.0;
    if ( !isVirtual() ) {
        rg_dList<VertexBU2D*> vertices;
        getBoundingVertices( vertices );

        rg_Point2D point[3] = { vertices.getFirstEntity()->getCoord(), 
                                vertices.getSecondEntity()->getCoord(), 
                                vertices.getLastEntity()->getCoord()  };

        signedArea  = rg_GeoFunc::computeSignedAreaOfTriangle(point[0], point[1], point[2]);
    }

    return signedArea;
}



///////////////////////////////////////////////////////////////////////////
//  Topological operators
rg_INT FaceBU2D::getBoundingVertices(rg_dList<VertexBU2D*>& boundingVertices) const
{
    rg_dList<VertexWEDS*> vertexListInWEDS;
    FaceWEDS::getBoundingVertices( vertexListInWEDS );
    
    vertexListInWEDS.reset4Loop();
    while ( vertexListInWEDS.setNext4Loop() ) {
        VertexWEDS* currVtx = vertexListInWEDS.getEntity();
        boundingVertices.add( (VertexBU2D*)currVtx );
    }

    return boundingVertices.getSize();
}



rg_INT FaceBU2D::getBoundingEdges(   rg_dList<EdgeBU2D*>&   boundingEdges) const
{
    rg_dList<EdgeWEDS*> edgeListInWEDS;
    FaceWEDS::getBoundingEdges( edgeListInWEDS );
    
    edgeListInWEDS.reset4Loop();
    while ( edgeListInWEDS.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgeListInWEDS.getEntity();
        boundingEdges.add( (EdgeBU2D*)currEdge );
    }

    return boundingEdges.getSize();
}



rg_INT FaceBU2D::getAdjacentFaces(   rg_dList<FaceBU2D*>&   faceList) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    FaceWEDS::getAdjacentFaces( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceWEDS* currFace = faceListInWEDS.getEntity();
        faceList.add( (FaceBU2D*)currFace );
    }

    return faceList.getSize();
}



///////////////////////////////////////////////////////////////////////////
//  Beta-operator
FaceBU2D* FaceBU2D::getOppositeFaceWithAnomalyInQuasiTriangulation() const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    FaceWEDS::getAdjacentFaces( faceListInWEDS );

    FaceWEDS* oppositeFaceWithAnomaly = rg_NULL;

    FaceWEDS** adjacentFace = faceListInWEDS.getArray();
    if ( adjacentFace[0] == adjacentFace[1] ) {
        oppositeFaceWithAnomaly = adjacentFace[0];
    }
    else if ( adjacentFace[0] == adjacentFace[2] ) {
        oppositeFaceWithAnomaly = adjacentFace[0];
    }
    else if ( adjacentFace[1] == adjacentFace[2] ) {
        oppositeFaceWithAnomaly = adjacentFace[1];
    }
    else {
        oppositeFaceWithAnomaly = rg_NULL;
    }

    delete [] adjacentFace;

    return (FaceBU2D*)oppositeFaceWithAnomaly;
}



FaceBU2D* FaceBU2D::getOppositeFaceWithAnomalyInBetaComplex(const rg_REAL& betaValue) const
{
    FaceBU2D* oppositeFaceInAnomaly = getOppositeFaceWithAnomalyInQuasiTriangulation();

    if ( oppositeFaceInAnomaly != rg_NULL ) {
        if (    this->getBoundingState(betaValue)                  != INTERIOR_SIMPLEX 
             || oppositeFaceInAnomaly->getBoundingState(betaValue) != INTERIOR_SIMPLEX ) {
            oppositeFaceInAnomaly = rg_NULL;
        }
    }

    return oppositeFaceInAnomaly;
}

rg_BOOL FaceBU2D::isAnomalyInQuasiTriangulation() const
{
    rg_BOOL isThis2MultiplicityAnomaly = rg_FALSE;

    FaceBU2D* oppositeFaceInAnomaly = getOppositeFaceWithAnomalyInQuasiTriangulation();
    if ( oppositeFaceInAnomaly != rg_NULL ) {
        isThis2MultiplicityAnomaly = rg_TRUE;
    }

    return isThis2MultiplicityAnomaly;
}



rg_BOOL FaceBU2D::isAnomalyInBetaComplex(const rg_REAL& betaValue) const
{
    rg_BOOL isThis2MultiplicityAnomaly = rg_FALSE;

    FaceBU2D* oppositeFaceInAnomaly = getOppositeFaceWithAnomalyInBetaComplex(betaValue);
    if ( oppositeFaceInAnomaly != rg_NULL ) {
        isThis2MultiplicityAnomaly = rg_TRUE;
    }

    return isThis2MultiplicityAnomaly;
}



rg_INT FaceBU2D::searchFacesInBetaComplex( const rg_REAL& beta, rg_dList<FaceBU2D*>& betaFaceList ) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    FaceWEDS::getAdjacentFaces( faceListInWEDS );
    
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



rg_INT FaceBU2D::searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<FaceBU2D*>& betaFaceList ) const
{
    rg_dList<FaceWEDS*> faceListInWEDS;
    FaceWEDS::getAdjacentFaces( faceListInWEDS );
    
    faceListInWEDS.reset4Loop();
    while ( faceListInWEDS.setNext4Loop() ) {
        FaceBU2D* currFace = (FaceBU2D*)(faceListInWEDS.getEntity());

        if ( currFace->getBoundingState(beta) == boundingState ) {
            betaFaceList.add( currFace );
        }
    }

    return betaFaceList.getSize();
}



rg_BOOL BULL2D::GeometryTier::FaceBU2D::isAdjacentTo( const FaceBU2D* const face ) const
{
    rg_BOOL isAdjacentToFace = false;

    rg_dList<FaceBU2D*> adjacentFaces;
    getAdjacentFaces( adjacentFaces );

    adjacentFaces.reset4Loop();
    while( adjacentFaces.setNext4Loop() ) {
        FaceBU2D* currFace = adjacentFaces.getEntity();

        if( currFace == face )
        {
            isAdjacentToFace = true;
            break;
        }
    }

    return isAdjacentToFace;
}



rg_BOOL BULL2D::GeometryTier::FaceBU2D::isIncidentTo( const EdgeBU2D* const edge ) const
{
    rg_BOOL isIncidentToEdge = false;

    rg_dList<EdgeBU2D*> boundingEdges;
    getBoundingEdges( boundingEdges );

    boundingEdges.reset4Loop();
    while( boundingEdges.setNext4Loop() ) {
        EdgeBU2D* currEdge = boundingEdges.getEntity();

        if( currEdge == edge )
        {
            isIncidentToEdge = true;
            break;
        }
    }

    return isIncidentToEdge;
}



rg_BOOL BULL2D::GeometryTier::FaceBU2D::isIncidentTo( const VertexBU2D* const vertex ) const
{
    rg_BOOL isIncidentToVertex = false;

    rg_dList<VertexBU2D*> boundingVertices;
    getBoundingVertices( boundingVertices );

    boundingVertices.reset4Loop();
    while( boundingVertices.setNext4Loop() ) {
        VertexBU2D* currVertex = boundingVertices.getEntity();

        if( currVertex == vertex )
        {
            isIncidentToVertex = true;
            break;
        }
    }

    return isIncidentToVertex;
}




