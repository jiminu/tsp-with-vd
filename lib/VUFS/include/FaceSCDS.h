#ifndef _FACESCDS_H
#define _FACESCDS_H


#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"

class VertexSCDS;



class FaceSCDS : public TopologicalEntity
{
protected:
    VertexSCDS* m_vertex[3];
    FaceSCDS*   m_neighbor[3];

    rg_BOOL     m_visited;

public:
    FaceSCDS();
    FaceSCDS(const rg_INT& ID);
    FaceSCDS(const FaceSCDS& face);
    ~FaceSCDS();

    VertexSCDS*     getVertex(const rg_INDEX& i) const;
    FaceSCDS*       getNeighbor(const rg_INDEX& i) const;
    inline VertexSCDS** getVertices() { return m_vertex; }
    inline FaceSCDS**   getNeighbors() { return m_neighbor; }

    inline rg_BOOL  isVisited() const                 { return m_visited; };
    inline void     isVisited(const rg_BOOL& visited) { m_visited = visited; };

    void            setVertex(  const rg_INDEX& i, VertexSCDS* vertex);
    void            setNeighbor(const rg_INDEX& i, FaceSCDS*   neighbor);
    void            setVertices(VertexSCDS* vertex1,   VertexSCDS* vertex2,   VertexSCDS* vertex3);
    void            setNeighbors(FaceSCDS*  neighbor1, FaceSCDS*   neighbor2, FaceSCDS*   neighbor3);

    FaceSCDS&       operator =(const FaceSCDS& face);


    ///////////////////////////////////////////////////////////////////////////

    FaceSCDS*       getMateFace(const rg_INT& pos) const;
    FaceSCDS*       getMateFace(VertexSCDS* vertex) const;

    rg_INT          getPosOfVertex(VertexSCDS* vertex) const;
    rg_INT          getPosOfNeighbor(FaceSCDS* neighbor) const;
    rg_INT          getThisPosInNeighbor(const rg_INT& neighborPos, FaceSCDS* neighbor) const;

    rg_BOOL         isThere(VertexSCDS* vertex) const;
    rg_BOOL         isThere(FaceSCDS*   neighbor) const;

    void            getBoundingVertices(rg_dList<VertexSCDS*> vertexList) const;
    void            getIncidentNeighbors(rg_dList<FaceSCDS*>  neighborList) const;
};

#endif

