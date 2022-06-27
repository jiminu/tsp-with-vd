#include "FaceWEDS.h"
#include "VertexWEDS.h"
#include "EdgeWEDS.h"


FaceWEDS::FaceWEDS()
: m_firstEdge(rg_NULL), m_visited(rg_FALSE)
{
}



FaceWEDS::FaceWEDS(const  rg_INT&	ID)
: TopologicalEntity(ID), m_firstEdge(rg_NULL), m_visited(rg_FALSE)
{
}



FaceWEDS::FaceWEDS(const  rg_INT&	ID, EdgeWEDS* firstEdge)
: TopologicalEntity(ID), m_firstEdge(firstEdge), m_visited(rg_FALSE)
{
}



FaceWEDS::FaceWEDS(const  FaceWEDS& face)
: TopologicalEntity(face), 
  m_firstEdge(face.m_firstEdge), m_visited(face.m_visited)
{
}



FaceWEDS::~FaceWEDS()
{
}




FaceWEDS& FaceWEDS::operator=(const FaceWEDS& face)
{
    if ( this == &face ) {
        return *this;
    }

    TopologicalEntity::operator =(face);
    m_firstEdge = face.m_firstEdge;
    m_visited   = face.m_visited;

    return *this;
}




//  Topological operators
rg_BOOL FaceWEDS::isThere(VertexWEDS* vertex) const
{
    EdgeWEDS* startEdge = m_firstEdge;
    EdgeWEDS* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            if ( vertex == currEdge->getEndVertex() ) {
                return rg_TRUE;
            }
            else {
                currEdge = currEdge->getLeftHand();
            }
        }
        else {
            if ( vertex == currEdge->getStartVertex() ) {
                return rg_TRUE;
            }
            else {            
                currEdge = currEdge->getRightLeg();
            }
        }
    } while ( startEdge != currEdge );   

    return rg_FALSE;
}



rg_BOOL FaceWEDS::isIncidentTo(EdgeWEDS* edge) const
{
    if ( this == edge->getLeftFace() || this == edge->getRightFace() ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_INT FaceWEDS::getBoundingVertices(rg_dList<VertexWEDS*>& boundingVertices) const
{
    EdgeWEDS* startEdge = m_firstEdge;
    EdgeWEDS* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            boundingVertices.add( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            boundingVertices.add( currEdge->getStartVertex() );
            currEdge = currEdge->getRightLeg();
        }
        else  {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return boundingVertices.getSize();
}



rg_INT FaceWEDS::getBoundingEdges(   rg_dList<EdgeWEDS*>&    boundingEdges) const
{
    EdgeWEDS* startEdge = m_firstEdge;
    EdgeWEDS* currEdge  = startEdge;

    do  {
        boundingEdges.add( currEdge );

        if ( this == currEdge->getLeftFace() ) {
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            currEdge = currEdge->getRightLeg();
        }
        else {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return boundingEdges.getSize();
}



rg_INT FaceWEDS::getAdjacentFaces(   rg_dList<FaceWEDS*>&   faceList) const
{
    EdgeWEDS* startEdge = m_firstEdge;
    EdgeWEDS* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            faceList.add( currEdge->getRightFace() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            faceList.add( currEdge->getLeftFace() );
            currEdge = currEdge->getRightLeg();
        }
        else  {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return faceList.getSize();
}



