#include "VVertex2D.h"
#include "VEdge2D.h"
#include "VFace2D.h"
#include "Generator2D.h"
using namespace BULL2D::GeometryTier;



VVertex2D::VVertex2D()
{
    m_ID         = -1;
    m_firstVEdge = NULL;
    m_location   = rg_Point2D( 0.0, 0.0 );
    m_status     = WHITE_V;
    m_isFictitious = false;
}



VVertex2D::VVertex2D( const int& ID )
{
    m_ID         = ID;
    m_firstVEdge = NULL;
    m_location   = rg_Point2D( 0.0, 0.0 );
    m_status     = WHITE_V;
    m_isFictitious = false;
}



VVertex2D::VVertex2D( const int& ID, const rg_Point2D& location )
{
    m_ID         = ID;
    m_firstVEdge = NULL;
    m_location   = location;
    m_status     = WHITE_V;
    m_isFictitious = false;
}



VVertex2D::VVertex2D( const int& ID, const rg_Point2D& location, VEdge2D* const firstEdge )
{
    m_ID         = ID;
    m_firstVEdge = firstEdge;
    m_location   = location;
    m_status     = WHITE_V;
    m_isFictitious = false;
}



VVertex2D::VVertex2D( const VVertex2D& vertex )
{
    m_ID         = vertex.m_ID;
    m_firstVEdge = vertex.m_firstVEdge;
    m_location   = vertex.m_location;
    m_status     = vertex.m_status;
    m_isFictitious = vertex.m_isFictitious;
}



VVertex2D::~VVertex2D()
{
}



VVertex2D& VVertex2D::operator=( const VVertex2D& vertex )
{
    if( this == &vertex ) {
        return *this;
    }

    m_ID         = vertex.m_ID;
    m_firstVEdge = vertex.m_firstVEdge;
    m_location   = vertex.m_location;
    m_status     = vertex.m_status;
    m_isFictitious = vertex.m_isFictitious;

    return *this;
}




bool VVertex2D::operator==( const VVertex2D& vertex ) const
{
    if( this == & vertex ) {
        return true;
    }

    if( m_firstVEdge == vertex.m_firstVEdge &&
        m_ID == vertex.m_ID &&
        m_location == vertex.m_location &&
        m_status == vertex.m_status)
    {
        return true;
    }
    else {
        return false;
    }
}



void VVertex2D::getIncident3VEdges( list<VEdge2D*>& incidentVEdgesList ) const
{
    VEdge2D* currVEdge = m_firstVEdge;

    do {
        incidentVEdgesList.push_back( currVEdge );

        if( this == currVEdge->getStartVertex() ) {
            currVEdge = currVEdge->getLeftLeg();
        }
        else {
            currVEdge = currVEdge->getRightHand();
        }
    } while ( currVEdge != m_firstVEdge );
}



void VVertex2D::getIncident3VFaces( list<VFace2D*>& incidentVFacesList ) const
{
    VEdge2D* currVEdge = m_firstVEdge;

    do {
        if( this == currVEdge->getStartVertex() ) {
            incidentVFacesList.push_back( currVEdge->getLeftFace() );
            currVEdge = currVEdge->getLeftLeg();
        }
        else {
            incidentVFacesList.push_back( currVEdge->getRightFace() );
            currVEdge = currVEdge->getRightHand();
        }
    } while ( currVEdge != m_firstVEdge );
}



void VVertex2D::getAdjacent3VVertices( list<VVertex2D*>& adjacentVVerticesList ) const
{
    VEdge2D* currVEdge = m_firstVEdge;

    do {
        if( this == currVEdge->getStartVertex() ) {
            adjacentVVerticesList.push_back( currVEdge->getEndVertex() );
            currVEdge = currVEdge->getLeftLeg();
        }
        else {
            adjacentVVerticesList.push_back( currVEdge->getStartVertex() );
            currVEdge = currVEdge->getRightHand();
        }
    } while ( currVEdge != m_firstVEdge );
}



void VVertex2D::getDefining3Generators( list<Generator2D*>& definingGeneratorsList ) const
{
    list<VFace2D*> incidentVFacesList;
    getIncident3VFaces( incidentVFacesList );

    for( list<VFace2D*>::iterator it_face = incidentVFacesList.begin() ; it_face != incidentVFacesList.end() ; it_face++ ) {
        VFace2D* currVFace = *it_face;

        if( currVFace->getGenerator() != NULL ) {
            definingGeneratorsList.push_back( currVFace->getGenerator() );
        }
    }
}




double VVertex2D::computeMUValue( Generator2D* const newGenerator ) // voronoi diagram이 갖는 편이 나을듯
{

    list<Generator2D*> definingGenerators;
    getDefining3Generators( definingGenerators );

    if( definingGenerators.size() < 3 ) {
        return DBL_MAX;
    }
    
    Generator2D* definingGenerator = *definingGenerators.begin();
    double       rho               = definingGenerator->getDisk()->getCenterPt().distance(m_location) - definingGenerator->getDisk()->getRadius();
    double       delta             = newGenerator->getDisk()->getCenterPt().distance(m_location) - newGenerator->getDisk()->getRadius();
    double       rhoSquare         = rho * rho;

    double mu       = ( delta - rho ) / rhoSquare ;

    return mu;

}



bool VVertex2D::isInfinite() const
{
    list<VFace2D*> incidentFaces;
    getIncident3VFaces( incidentFaces );

    list<VFace2D*>::iterator it_incidentFace = incidentFaces.begin();
    for( ; it_incidentFace != incidentFaces.end() ; it_incidentFace++ ) {
        VFace2D* currFace = *it_incidentFace;
        if( currFace->isInfinite() ) {
            return true;
        }
    }

    return false;
}



double VVertex2D::computeRadiusOfTangentCircle() const
{
    if( isInfinite() ) {
        return DBL_MAX;
    }

    Generator2D*    definingGenerator = m_firstVEdge->getLeftFace()->getGenerator();
    double          distanceToCenterOfDefiningGenerator = definingGenerator->getDisk()->getCenterPt().distance( m_location );
    double          radiusOfDefiningGenerator = definingGenerator->getDisk()->getRadius();
    double          radiusOfTangentCircle = 0.0;

    radiusOfTangentCircle = distanceToCenterOfDefiningGenerator - radiusOfDefiningGenerator;

    return radiusOfTangentCircle;
}



bool VVertex2D::isMemberOf( const VFace2D* const face ) const
{
    bool isMemberOfFace = false;

    list<VFace2D*> incidentFaces;
    getIncident3VFaces( incidentFaces );

    list<VFace2D*>::iterator i_face;
    for( i_face = incidentFaces.begin() ; i_face != incidentFaces.end() ; i_face++ ) {
        VFace2D* currFace = *i_face;

        if( currFace == face )
        {
            isMemberOfFace = true;
            break;
        }
    }

    return isMemberOfFace;
}



bool VVertex2D::isMemberOf( const VEdge2D* const edge ) const
{
    bool isMemberOfEdge = false;

    list<VEdge2D*> incidentEdges;
    getIncident3VEdges( incidentEdges );

    list<VEdge2D*>::iterator i_edge;
    for( i_edge = incidentEdges.begin() ; i_edge != incidentEdges.end() ; i_edge++ ) {
        VEdge2D* currEdge = *i_edge;

        if( currEdge == edge )
        {
            isMemberOfEdge = true;
            break;
        }
    }

    return isMemberOfEdge;
}



bool VVertex2D::isAdjacentTo( const VVertex2D* const vertex ) const
{
    bool isAdjacentToVertex = false;

    list<VVertex2D*> adjacentVertices;
    getAdjacent3VVertices( adjacentVertices );

    list<VVertex2D*>::iterator i_vertex;
    for( i_vertex = adjacentVertices.begin() ; i_vertex != adjacentVertices.end() ; i_vertex++ ) {
        VVertex2D* currVertex = *i_vertex;

        if( currVertex == vertex )
        {
            isAdjacentToVertex = true;
            break;
        }
    }

    return isAdjacentToVertex;
}




