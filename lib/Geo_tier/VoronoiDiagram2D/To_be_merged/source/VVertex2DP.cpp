#include "VVertex2DP.h"
#include "VEdge2DP.h"
#include "VFace2DP.h"
#include "Generator2DP.h"

using namespace BULL2D::GeometryTier;


VVertex2DP::VVertex2DP()
{
    m_ID         = -1;
    m_firstVEdge = NULL;
    m_location   = rg_Point2D( 0.0, 0.0 );
    m_status     = WHITE_V_P;
}



VVertex2DP::VVertex2DP( const int& ID )
{
    m_ID         = ID;
    m_firstVEdge = NULL;
    m_location   = rg_Point2D( 0.0, 0.0 );
    m_status     = WHITE_V_P;
}



VVertex2DP::VVertex2DP( const int& ID, const rg_Point2D& location )
{
    m_ID         = ID;
    m_firstVEdge = NULL;
    m_location   = location;
    m_status     = WHITE_V_P;
}



VVertex2DP::VVertex2DP( const int& ID, const rg_Point2D& location, VEdge2DP* const firstEdge )
{
    m_ID         = ID;
    m_firstVEdge = firstEdge;
    m_location   = location;
    m_status     = WHITE_V_P;
}



VVertex2DP::VVertex2DP( const VVertex2DP& vertex )
{
    m_ID         = vertex.m_ID;
    m_firstVEdge = vertex.m_firstVEdge;
    m_location   = vertex.m_location;
    m_status     = vertex.m_status;
}



VVertex2DP::~VVertex2DP()
{
}



VVertex2DP& VVertex2DP::operator=( const VVertex2DP& vertex )
{
    if( this == &vertex ) {
        return *this;
    }

    m_ID         = vertex.m_ID;
    m_firstVEdge = vertex.m_firstVEdge;
    m_location   = vertex.m_location;
    m_status     = vertex.m_status;

    return *this;
}



void VVertex2DP::getIncident3VEdges( list<VEdge2DP*>& incidentVEdgesList ) const
{
    VEdge2DP* currVEdge = m_firstVEdge;

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



void VVertex2DP::getIncident3VFaces( list<VFace2DP*>& incidentVFacesList ) const
{
    VEdge2DP* currVEdge = m_firstVEdge;

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



void VVertex2DP::getAdjacent3VVertices( list<VVertex2DP*>& adjacentVVerticesList ) const
{
    VEdge2DP* currVEdge = m_firstVEdge;

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



void VVertex2DP::getDefining3Generators( list<Generator2DP*>& definingGeneratorsList ) const
{
    list<VFace2DP*> incidentVFacesList;
    getIncident3VFaces( incidentVFacesList );

    for( list<VFace2DP*>::iterator it_face = incidentVFacesList.begin() ; it_face != incidentVFacesList.end() ; it_face++ ) {
        VFace2DP* currVFace = *it_face;

        if( currVFace->getGenerator() != NULL ) {
            definingGeneratorsList.push_back( currVFace->getGenerator() );
        }
    }
}



double VVertex2DP::computeHValue( Generator2DP* const newGenerator )
{
    list<Generator2DP*> definingGenerators;
    getDefining3Generators( definingGenerators );

    if( definingGenerators.size() < 3 ) {
        return DBL_MAX;
    }
    else {
        list<Generator2DP*>::iterator it_generator = definingGenerators.begin();
        Generator2DP* firstGenerator     = *it_generator++;
        Generator2DP* secondGenerator    = *it_generator++;
        Generator2DP* thirdGenerator     = *it_generator;

        // calculate H value
        double x1 = 0.0;
        double y1 = 0.0;
        double z1 = 0.0;

        double x2 = 0.0;
        double y2 = 0.0;
        double z2 = 0.0;

        double x3 = 0.0;
        double y3 = 0.0;

        double x4 = 0.0;
        double y4 = 0.0;
        double z4 = 0.0;

        x1 = firstGenerator->getLocation()->getX();
        y1 = firstGenerator->getLocation()->getY();

        x2 = secondGenerator->getLocation()->getX();
        y2 = secondGenerator->getLocation()->getY();

        x3 = thirdGenerator->getLocation()->getX();
        y3 = thirdGenerator->getLocation()->getY();

        x4 = newGenerator->getLocation()->getX();
        y4 = newGenerator->getLocation()->getY();

        x1 = x1 - x3;
        y1 = y1 - y3;
        z1 = x1*x1 + y1*y1;

        x2 = x2 - x3;
        y2 = y2 - y3;
        z2 = x2*x2 + y2*y2;

        x4 = x4 - x3;
        y4 = y4 - y3;
        z4 = x4*x4 + y4*y4;

        double j2 = 0.0;
        double j3 = 0.0;
        double j4 = 0.0;

        j2 = y1*z2 - z1*y2;
        j3 = x1*z2 - z1*x2;
        j4 = x1*y2 - y1*x2;

        double hValue = j2*x4 - j3*y4 + j4*z4;

        return hValue;
    }
}



double VVertex2DP::computeRadiusOfTangentCircle() const
{
    if( isInfinite() ) {
        return DBL_MAX;
    }

    Generator2DP*   definingGenerator = m_firstVEdge->getLeftFace()->getGenerator();
    double          distanceToCenterOfDefiningGenerator = definingGenerator->getLocation()->distance( m_location );
    double          radiusOfTangentCircle = distanceToCenterOfDefiningGenerator;

    return radiusOfTangentCircle;
}



bool VVertex2DP::isInfinite() const
{
    list<VFace2DP*> incidentFaces;
    getIncident3VFaces( incidentFaces );

    list<VFace2DP*>::iterator i_incidentFace = incidentFaces.begin();
    for( ; i_incidentFace != incidentFaces.end() ; i_incidentFace++ ) {
        VFace2DP* currFace = *i_incidentFace;
        if( currFace->isInfinite() )
        {
            return true;
        }
    }

    return false;
}



bool VVertex2DP::isMemberOf( const VFace2DP* const face ) const
{
    bool    isMemberOfFace = false;

    list<VFace2DP*> incidentFaces;
    getIncident3VFaces( incidentFaces );

    list<VFace2DP*>::iterator i_face;
    for( i_face = incidentFaces.begin() ; i_face != incidentFaces.end() ; i_face++ ) {
        VFace2DP* currFace = *i_face;

        if( currFace == face )
        {
            isMemberOfFace = true;
            break;
        }
    }

    return isMemberOfFace;
}



bool VVertex2DP::isMemberOf( const VEdge2DP* const edge ) const
{
    bool    isMemberOfEdge = false;

    list<VEdge2DP*> incidentEdges;
    getIncident3VEdges( incidentEdges );

    list<VEdge2DP*>::iterator i_edge;
    for( i_edge = incidentEdges.begin() ; i_edge != incidentEdges.end() ; i_edge++ ) {
        VEdge2DP* currEdge = *i_edge;

        if( currEdge == edge )
        {
            isMemberOfEdge = true;
            break;
        }
    }

    return isMemberOfEdge;
}



bool VVertex2DP::isAdjacentTo( const VVertex2DP* const vertex ) const
{
    bool isAdjacentToVertex = false;

    list<VVertex2DP*> adjacentVertices;
    getAdjacent3VVertices( adjacentVertices );

    list<VVertex2DP*>::iterator i_vertex;
    for( i_vertex = adjacentVertices.begin() ; i_vertex != adjacentVertices.end() ; i_vertex++ ) {
        VVertex2DP* currVertex = *i_vertex;

        if( currVertex == vertex )
        {
            isAdjacentToVertex = true;
            break;
        }
    }

    return isAdjacentToVertex;
}



