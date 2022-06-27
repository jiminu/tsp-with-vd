#include "VFace2D.h"
#include "VEdge2D.h"
#include "VVertex2D.h"
#include "Generator2D.h"
using namespace BULL2D::GeometryTier;



VFace2D::VFace2D()
{
    m_ID            = -1;
    m_firstVEdge    = NULL;
    m_generator     = NULL;
}



VFace2D::VFace2D( const int& ID )
{
    m_ID            = ID;
    m_firstVEdge    = NULL;
    m_generator     = NULL;
}



VFace2D::VFace2D( const int& ID, VEdge2D* const firstVEdge )
{
    m_ID            = ID;
    m_firstVEdge    = firstVEdge;
    m_generator     = NULL;
}



VFace2D::VFace2D( const int& ID, VEdge2D* const firstVEdge, Generator2D* const generator )
{
    m_ID            = ID;
    m_firstVEdge    = firstVEdge;
    m_generator     = generator;
}



VFace2D::VFace2D( const VFace2D& face )
{
    m_ID            = face.m_ID;
    m_firstVEdge    = face.m_firstVEdge;
    m_generator     = face.m_generator;
}



VFace2D::~VFace2D()
{
}



VFace2D& VFace2D::operator=( const VFace2D& face )
{
    if( this == &face ) {
        return *this;
    }

    m_ID            = face.m_ID;
    m_firstVEdge    = face.m_firstVEdge;
    m_generator     = face.m_generator;

    return *this;
}




bool VFace2D::operator==( const VFace2D& face ) const
{
    if( this == &face ) {
        return true;
    }

    if( m_ID == face.m_ID &&
        m_firstVEdge == face.m_firstVEdge &&
        m_generator == face.m_generator )
    {
        return true;
    }
    else {
        return false;
    }
}



void VFace2D::getBoundaryVEdges( list<VEdge2D*>& boundaryEdgesList ) const
{
    VEdge2D* currVEdge = m_firstVEdge;

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



void VFace2D::getBoundaryVVertices( list<VVertex2D*>& boundaryVerticesList ) const
{
    VEdge2D* currVEdge = m_firstVEdge;

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



void VFace2D::getAdjacentVFaces( list<VFace2D*>& adjacentFacesList ) const
{
    VEdge2D* currVEdge = m_firstVEdge;

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



bool VFace2D::isInfinite() const
{
    if( m_generator == NULL ) {
        return true;
    }
    else {
        return false;
    }
}



bool VFace2D::isBounded() const
{
    if( isInfinite() ) {
        return false;
    }

    list<VFace2D*> adjacentFaces;
    getAdjacentVFaces( adjacentFaces );

    list<VFace2D*>::iterator it_adjFace = adjacentFaces.begin();
    for( ; it_adjFace != adjacentFaces.end() ; it_adjFace++ ) {
        VFace2D* currAdjFace = *it_adjFace;
        if( currAdjFace->isInfinite() ) {
            return false;
        }
    }

    return true;
}



bool VFace2D::isUnBounded() const
{
    if( isInfinite() ) {
        return false;
    }

    list<VFace2D*> adjacentFaces;
    getAdjacentVFaces( adjacentFaces );

    list<VFace2D*>::iterator it_adjFace = adjacentFaces.begin();
    for( ; it_adjFace != adjacentFaces.end() ; it_adjFace++ ) {
        VFace2D* currAdjFace = *it_adjFace;
        if( currAdjFace->isInfinite() ) {
            return true;
        }
    }

    return false;
}



bool VFace2D::isAdjacentTo( const VFace2D* const face ) const
{
    bool isAdjacentToFace = false;

    list<VFace2D*> adjacentFaces;
    getAdjacentVFaces( adjacentFaces );

    list<VFace2D*>::iterator i_face;
    for( i_face = adjacentFaces.begin() ; i_face != adjacentFaces.end() ; i_face++ ) {
        VFace2D* currFace = *i_face;

        if( currFace == face )
        {
            isAdjacentToFace = true;
            break;
        }
    }

    return isAdjacentToFace;
}



bool VFace2D::isIncidentTo( const VEdge2D* const edge ) const
{
    bool isIncidentToEdge = false;

    list<VEdge2D*> boundaryEdges;
    getBoundaryVEdges( boundaryEdges );

    list<VEdge2D*>::iterator i_edge;
    for( i_edge = boundaryEdges.begin() ; i_edge != boundaryEdges.end() ; i_edge++ ) {
        VEdge2D* currEdge = *i_edge;

        if( currEdge == edge )
        {
            isIncidentToEdge = true;
            break;
        }
    }

    return isIncidentToEdge;
}



bool VFace2D::isIncidentTo( const VVertex2D* const vertex ) const
{
    bool isIncidentToVertex = false;

    list<VVertex2D*> boundaryVertices;
    getBoundaryVVertices( boundaryVertices );

    list<VVertex2D*>::iterator i_vertex;
    for( i_vertex = boundaryVertices.begin() ; i_vertex != boundaryVertices.end() ; i_vertex++ ) {
        VVertex2D* currVertex = *i_vertex;

        if( currVertex == vertex )
        {
            isIncidentToVertex = true;
            break;
        }
    }

    return isIncidentToVertex;
}


