#include "VVertex2D.h"
#include "VEdge2D.h"
#include "VFace2D.h"
#include "Generator2D.h"
using namespace V::GeometryTier;


VVertex2D::VVertex2D()
{
    m_ID            = -1;
    m_firstVEdge    = NULL;
    m_circumcircle  = rg_Circle2D( 0.0, 0.0, 0.0 );
    m_status        = WHITE_V;
    m_isFictitious  = false;
    m_userData      = NULL;
}



VVertex2D::VVertex2D( const int& ID )
{
    m_ID         = ID;
    m_firstVEdge = NULL;
    m_circumcircle  = rg_Circle2D( 0.0, 0.0, 0.0 );
    m_status     = WHITE_V;
    m_isFictitious = false;
    m_userData      = NULL;
}



VVertex2D::VVertex2D( const int& ID, const rg_Point2D& location )
{
    m_ID            = ID;
    m_firstVEdge    = NULL;
    m_circumcircle  = rg_Circle2D( location.getX(), location.getY(), 0.0 );
    m_status        = WHITE_V;
    m_isFictitious  = false;
    m_userData      = NULL;
}



VVertex2D::VVertex2D( const int& ID, const rg_Circle2D& circumcircle )
{
    m_ID            = ID;
    m_firstVEdge    = NULL;
    m_circumcircle  = circumcircle;
    m_status        = WHITE_V;
    m_isFictitious  = false;
    m_userData      = NULL;
}



VVertex2D::VVertex2D( const int& ID, const rg_Point2D& location, VEdge2D* const firstEdge )
{
    m_ID            = ID;
    m_firstVEdge    = firstEdge;
    m_circumcircle  = rg_Circle2D( location.getX(), location.getY(), 0.0 );
    m_status        = WHITE_V;
    m_isFictitious  = false;
    m_userData      = NULL;
}



VVertex2D::VVertex2D( const int& ID, const rg_Circle2D& circumcircle, VEdge2D* const firstEdge )
{
    m_ID            = ID;
    m_firstVEdge    = firstEdge;
    m_circumcircle  = circumcircle;
    m_status        = WHITE_V;
    m_isFictitious  = false;
    m_userData      = NULL;
}



VVertex2D::VVertex2D( const VVertex2D& vertex )
{
    m_ID            = vertex.m_ID;
    m_firstVEdge    = vertex.m_firstVEdge;
    m_circumcircle  = vertex.m_circumcircle;
    m_status        = vertex.m_status;
    m_isFictitious  = vertex.m_isFictitious;
    m_userData      = vertex.m_userData;

    // CYSONG ADDED THE FOLLOWINGS. (SEP 16, 20)
    m_tangentCircle = vertex.m_tangentCircle;
}



VVertex2D::~VVertex2D()
{
}



VVertex2D& VVertex2D::operator=( const VVertex2D& vertex )
{
    if( this == &vertex ) {
        return *this;
    }

    m_ID            = vertex.m_ID;
    m_firstVEdge    = vertex.m_firstVEdge;
    m_circumcircle  = vertex.m_circumcircle;
    m_status        = vertex.m_status;
    m_isFictitious  = vertex.m_isFictitious;
    m_userData      = vertex.m_userData;
    return *this;
}




bool VVertex2D::operator==( const VVertex2D& vertex ) const
{
    if( this == & vertex ) {
        return true;
    }

    if( m_firstVEdge == vertex.m_firstVEdge &&
        m_ID == vertex.m_ID &&
        m_circumcircle == vertex.m_circumcircle &&
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

    return m_circumcircle.getRadius();
}



bool V::GeometryTier::VVertex2D::isMemberOf( const VFace2D* const face ) const
{
    bool ret = false;
    list<VFace2D*> incidentFaces;
    getIncident3VFaces( incidentFaces );
    for ( list<VFace2D*>::iterator i_face = incidentFaces.begin(); i_face != incidentFaces.end(); ++i_face ) {
        VFace2D* currFace = *i_face;
        if ( currFace == face ) {
            ret = true;
            break;
        }
    }

    return ret;
}



bool V::GeometryTier::VVertex2D::isMemberOf( const VEdge2D* const edge ) const
{
    bool ret = false;
    if ( edge->getStartVertex() == this || edge->getEndVertex() == this ) {
        ret = true;
    }

    return ret;
}



bool V::GeometryTier::VVertex2D::isAdjacentTo( const VVertex2D* const vertex ) const
{
    bool ret = false;
    list<VVertex2D*> adjacentVertices;
    getAdjacent3VVertices( adjacentVertices );
    for ( list<VVertex2D*>::iterator i_vtx = adjacentVertices.begin(); i_vtx != adjacentVertices.end(); ++i_vtx ) {
        VVertex2D* currVtx = *i_vtx;
        if ( currVtx == vertex ) {
            ret = true;
            break;
        }
    }

    return ret;
}
