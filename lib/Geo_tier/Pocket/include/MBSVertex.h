#ifndef _MBSVERTEX_H
#define _MBSVERTEX_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"
#include "Ball.h"

class MBSEdge;
class MBSFace;
class MBSShell;

namespace V{
    namespace GeometryTier{

class BetaVertex;

    }
}

using namespace V::GeometryTier;

class MBSVertex : public TopologicalEntity
{
private:
    MBSShell*   m_shell;

    MBSEdge*    m_firstEdge;

    BetaVertex* m_originalBetaVertex;
    rg_BOOL     m_isArtificial;

    rg_BOOL     m_visited;

public:
    MBSVertex();
    MBSVertex(const rg_INT& ID);
    MBSVertex(const MBSVertex& vertex);
    ~MBSVertex();

    MBSShell*       getShell() const;
    MBSEdge*        getFirstEdge() const;
    
    BetaVertex*     getOriginalBetaVertex() const;
    rg_BOOL         isArtificial() const;

    inline rg_BOOL  isVisited() const { return m_visited; }
    inline void     isVisited(const rg_BOOL& visited) { m_visited = visited; }


    void            setShell( MBSShell* shell );
    void            setFirstEdge(MBSEdge* firstEdge);
    
    void            setOriginalBetaVertex( BetaVertex* originalBetaVertex );
    void            isArtificial( const rg_BOOL& isArtificial );


    MBSVertex&      operator =(const MBSVertex& o_vertex);

    rg_FLAG         searchAdjacentVertices(rg_dList<MBSVertex*>& vertexList) const;
    rg_FLAG         searchIncidentEdges(   rg_dList<MBSEdge*>&   edgeList)   const;
    rg_FLAG         searchIncidentFaces(   rg_dList<MBSFace*>&   faceList)   const;

    
};

#endif
