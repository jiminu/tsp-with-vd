#ifndef VFACE2DP_H
#define VFACE2DP_H


#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {


class VVertex2DP;
class VEdge2DP;
class Generator2DP;

class VFace2DP
{
private:
    int             m_ID;
    VEdge2DP*        m_firstVEdge;
    Generator2DP*    m_generator;

public:
    VFace2DP();
    VFace2DP( const int& ID );
    VFace2DP( const int& ID, VEdge2DP* const firstVEdge );
    VFace2DP( const int& ID, VEdge2DP* const firstVEdge, Generator2DP* const generator );
    VFace2DP( const VFace2DP& face );
    ~VFace2DP();

    inline int          getID() const { return m_ID; }
    inline VEdge2DP*     getFirstVEdge() const { return m_firstVEdge; }
    inline Generator2DP* getGenerator() const { return m_generator; }

    inline void         setID( const int& ID ) { m_ID = ID; }
    inline void         setFirstVEdge( VEdge2DP* const edge ) { m_firstVEdge = edge; }
    inline void         setGenerator( Generator2DP* const generator ) { m_generator = generator; }

    VFace2DP&            operator=( const VFace2DP& face );

    void                getBoundaryVEdges( list<VEdge2DP*>& boundaryEdgesList ) const;
    void                getBoundaryVVertices( list<VVertex2DP*>& boundaryVerticesList ) const;
    void                getAdjacentVFaces( list<VFace2DP*>& adjacentFacesList ) const;

    bool    isAdjacentTo( const VFace2DP* const face ) const;
    bool    isIncidentTo( const VEdge2DP* const edge ) const;
    bool    isIncidentTo( const VVertex2DP* const vertex ) const;
    bool    isInfinite() const;
    bool    isBounded() const;
    bool    isUnBounded() const;
};

} // namespace GeometryTier
} // BULL2D

#endif


