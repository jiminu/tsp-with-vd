#ifndef _VERTEXSCDS_H
#define _VERTEXSCDS_H

#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"


class FaceSCDS;

class VertexSCDS : public TopologicalEntity
{
protected:
    FaceSCDS* m_firstFace;

    rg_BOOL   m_visited;

public:
    VertexSCDS();
    VertexSCDS(const rg_INT& ID);
    VertexSCDS(const VertexSCDS& vertex);
    ~VertexSCDS();

    inline FaceSCDS* getFirstFace() const              { return m_firstFace; };

    inline rg_BOOL   isVisited() const                 { return m_visited; };
    inline void      isVisited(const rg_BOOL& visited) { m_visited = visited; };

    void   setFirstFace(FaceSCDS* firstFace);

    VertexSCDS& operator =(const VertexSCDS& vertex);

    ///////////////////////////////////////////////////////////////////////////

    void getIncidentFace(rg_dList<FaceSCDS*>& faceList);

};


#endif

