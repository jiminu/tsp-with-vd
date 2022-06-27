#ifndef _MBSFACE_H
#define _MBSFACE_H

#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"

class MBSShell;
class MBSEdge;
class MBSVertex;

namespace V{
    namespace GeometryTier{
class BetaFace;
    }
}

using namespace V::GeometryTier;

class MBSFace : public TopologicalEntity
{
private:
    MBSEdge*            m_firstEdge;
    

    BetaFace*           m_betaFace;
    rg_BOOL             m_isArtificial;


    rg_BOOL             m_visited;

public:
    MBSFace();
    MBSFace(const rg_INT& ID);
    MBSFace(const MBSFace& face);
    ~MBSFace();

    MBSShell*           getShell() const;
    MBSEdge*            getFirstEdge() const;

    BetaFace*           getBetaFace() const;
    rg_BOOL             isArtificial() const;

    void                setFirstEdge( MBSEdge* firstEdge );

    void                setBetaFace( BetaFace* betaFace );
    void                isArtificial( const rg_BOOL& isArtificial );

    inline rg_BOOL      isVisited() const { return m_visited; }
    inline void         isVisited(const rg_BOOL& visited) { m_visited = visited; }
    
    MBSFace&             operator =(const MBSFace& o_face);

    
    void                searchBoundingVertices( rg_dList<MBSVertex*>& vertexList );
    void                searchBoundingEdges( rg_dList<MBSEdge*>& edgeList );
    rg_FLAG             searchIncidentFaces( rg_dList<MBSFace*>& faceList );

    rg_FLAG             isZeroVolumeFace();
    
};

#endif

