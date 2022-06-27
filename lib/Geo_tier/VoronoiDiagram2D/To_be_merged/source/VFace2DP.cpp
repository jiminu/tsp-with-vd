#include "VFace2DP.h"
#include "VEdge2DP.h"
#include "VVertex2DP.h"
#include "Generator2DP.h"

using namespace BULL2D::GeometryTier;

VFace2DP::VFace2DP()
{
    m_ID            = -1;
    m_firstVEdge    = NULL;
    m_generator     = NULL;
}



VFace2DP::VFace2DP( const int& ID )
{
    m_ID            = ID;
    m_firstVEdge    = NULL;
    m_generator     = NULL;
}



VFace2DP::VFace2DP( const int& ID, VEdge2DP* const firstVEdge )
{
    m_ID            = ID;
    m_firstVEdge    = firstVEdge;
    m_generator     = NULL;
}



VFace2DP::VFace2DP( const int& ID, VEdge2DP* const firstVEdge, Generator2DP* const generator )
{
    m_ID            = ID;
    m_firstVEdge    = firstVEdge;
    m_generator     = generator;
}



VFace2DP::VFace2DP( const VFace2DP& face )
{
    m_ID            = face.m_ID;
    m_firstVEdge    = face.m_firstVEdge;
    m_generator     = face.m_generator;
}



VFace2DP::~VFace2DP()
{
}



VFace2DP& VFace2DP::operator=( const VFace2DP& face )
{
    if( this == &face ) {
        return *this;
    }

    m_ID            = face.m_ID;
    m_firstVEdge    = face.m_firstVEdge;
    m_generator     = face.m_generator;

    return *this;
}



void VFace2DP::getBoundaryVEdges( list<VEdge2DP*>& boundaryEdgesList ) const
{
    VEdge2DP* currVEdge = m_firstVEdge;

    do {
        boundaryEdgesList.push_back( currVEdge );

        if( this == currVEdge->getLeftFace() ) {
            currVEdge = currVEdge->getLeftHand();
        }
        else {
            currVEdge = currVEdge->getRightLeg();
        }
    } while( currVEdge != m_firstVEdge);
}



void VFace2DP::getBoundaryVVertices( list<VVertex2DP*>& boundaryVerticesList ) const
{
    VEdge2DP* currVEdge = m_firstVEdge;

    do {
        if( this == currVEdge->getLeftFace() ) {
            boundaryVerticesList.push_back( currVEdge->getStartVertex() );
            currVEdge = currVEdge->getLeftHand();
        }
        else {
            boundaryVerticesList.push_back( currVEdge->getEndVertex() );
            currVEdge = currVEdge->getRightLeg();
        }
    } while( currVEdge != m_firstVEdge);
}



void VFace2DP::getAdjacentVFaces( list<VFace2DP*>& adjacentFacesList ) const
{
    VEdge2DP* currVEdge = m_firstVEdge;

    do {
        if( this == currVEdge->getLeftFace() ) {
            adjacentFacesList.push_back( currVEdge->getRightFace() );
            currVEdge = currVEdge->getLeftHand();
        }
        else {
            adjacentFacesList.push_back( currVEdge->getLeftFace() );
            currVEdge = currVEdge->getRightLeg();
        }
    } while( currVEdge != m_firstVEdge);
}



bool VFace2DP::isAdjacentTo( const VFace2DP* const face ) const
{
    bool isAdjacentToFace = false;

    list<VFace2DP*> adjacentFaces;
    getAdjacentVFaces( adjacentFaces );

    list<VFace2DP*>::iterator i_face;
    for( i_face = adjacentFaces.begin() ; i_face != adjacentFaces.end() ; i_face++ ) {
        VFace2DP* currFace = *i_face;

        if( currFace == face )
        {
            isAdjacentToFace = true;
            break;
        }
    }

    return isAdjacentToFace;
}



bool VFace2DP::isIncidentTo( const VEdge2DP* const edge ) const
{
    bool isIncidentToEdge = false;

    list<VEdge2DP*> boundaryEdges;
    getBoundaryVEdges( boundaryEdges );

    list<VEdge2DP*>::iterator i_edge;
    for( i_edge = boundaryEdges.begin() ; i_edge != boundaryEdges.end() ; i_edge++ ) {
        VEdge2DP* currEdge = *i_edge;

        if( currEdge == edge )
        {
            isIncidentToEdge = true;
            break;
        }
    }

    return isIncidentToEdge;
}



bool VFace2DP::isIncidentTo( const VVertex2DP* const vertex ) const
{
    bool isIncidentToVertex = false;

    list<VVertex2DP*> boundaryVertices;
    getBoundaryVVertices( boundaryVertices );

    list<VVertex2DP*>::iterator i_vertex;
    for( i_vertex = boundaryVertices.begin() ; i_vertex != boundaryVertices.end() ; i_vertex++ ) {
        VVertex2DP* currVertex = *i_vertex;

        if( currVertex == vertex )
        {
            isIncidentToVertex = true;
            break;
        }
    }

    return isIncidentToVertex;
}



bool VFace2DP::isInfinite() const
{
    if( m_generator == NULL ) {
        return true;
    }
    else {
        return false;
    }
}



bool VFace2DP::isBounded() const
{
    if( isInfinite() ) {
        return false;
    }

    list<VFace2DP*> adjacentFaces;
    getAdjacentVFaces( adjacentFaces );

    list<VFace2DP*>::iterator it_adjFace = adjacentFaces.begin();
    for( ; it_adjFace != adjacentFaces.end() ; it_adjFace++ ) {
        VFace2DP* currAdjFace = *it_adjFace;
        if( currAdjFace->isInfinite() ) {
            return false;
        }
    }

    return true;
}



bool VFace2DP::isUnBounded() const
{
    if( isInfinite() ) {
        return false;
    }

    list<VFace2DP*> adjacentFaces;
    getAdjacentVFaces( adjacentFaces );

    list<VFace2DP*>::iterator it_adjFace = adjacentFaces.begin();
    for( ; it_adjFace != adjacentFaces.end() ; it_adjFace++ ) {
        VFace2DP* currAdjFace = *it_adjFace;
        if( currAdjFace->isInfinite() ) {
            return true;
        }
    }

    return false;
}


