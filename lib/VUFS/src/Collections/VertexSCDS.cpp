#include "VertexSCDS.h"
#include "FaceSCDS.h"



VertexSCDS::VertexSCDS()
: m_firstFace(rg_NULL), m_visited(rg_FALSE)
{
}



VertexSCDS::VertexSCDS(const rg_INT& ID)
: TopologicalEntity(ID),
  m_firstFace(rg_NULL), m_visited(rg_FALSE)
{
}



VertexSCDS::VertexSCDS(const VertexSCDS& vertex)
: TopologicalEntity(vertex),
  m_firstFace(vertex.m_firstFace), m_visited(vertex.m_visited)
{
}



VertexSCDS::~VertexSCDS()
{
}



void   VertexSCDS::setFirstFace(FaceSCDS* firstFace)
{
    m_firstFace = firstFace;
}




VertexSCDS& VertexSCDS::operator =(const VertexSCDS& vertex)
{
    if ( this == &vertex ) {
        return *this;
    }

    m_ID        = vertex.m_ID;
    m_firstFace = vertex.m_firstFace;
    m_visited   = vertex.m_visited;

    return *this;
}



void VertexSCDS::getIncidentFace(rg_dList<FaceSCDS*>& faceList)
{
    FaceSCDS* currFace = m_firstFace;

    do {
        faceList.add( currFace );

        rg_INT pos = currFace->getPosOfVertex(this);
        switch ( pos ) {
            case 0:
                currFace = currFace->getNeighbor(1);
                break;
            case 1:
                currFace = currFace->getNeighbor(2);
                break;
            case 2:
                currFace = currFace->getNeighbor(0);
               break;
            default:
                currFace = m_firstFace;
                break;
        }

    } while (currFace != m_firstFace);
}


