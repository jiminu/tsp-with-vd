#ifndef VFACE2D_H
#define VFACE2D_H

#include <list>
using namespace std;

namespace V {
namespace GeometryTier {

class VVertex2D;
class VEdge2D;
class Generator2D;

class VFace2D
{
private:
    int             m_ID;
    VEdge2D*        m_firstVEdge;
    Generator2D*    m_generator;

    void*           m_userData; // August 10, 2020 by CYSONG.


public:
    VFace2D();
    VFace2D( const int& ID );
    VFace2D( const int& ID, VEdge2D* const firstVEdge );
    VFace2D( const int& ID, VEdge2D* const firstVEdge, Generator2D* const generator );
    VFace2D( const VFace2D& face );
    ~VFace2D();

    inline int          getID() const { return m_ID; }
    inline VEdge2D*     getFirstVEdge() const { return m_firstVEdge; }
    inline Generator2D* getGenerator() const { return m_generator; }
    inline void*        getUserData()  const { return m_userData;}

    inline void         setID( const int& ID ) { m_ID = ID; }
    inline void         setFirstVEdge( VEdge2D* const edge ) { m_firstVEdge = edge; }
    inline void         setGenerator( Generator2D* const generator ) { m_generator = generator; }
    inline void         setUserData( void* userData ) { m_userData = userData; }

    VFace2D&            operator=( const VFace2D& face );
    bool                operator==( const VFace2D& face ) const;

    void                getBoundaryVEdges( list<VEdge2D*>& boundaryEdgesList ) const;
    void                getBoundaryVVertices( list<VVertex2D*>& boundaryVerticesList ) const;
    void                getAdjacentVFaces( list<VFace2D*>& adjacentFacesList ) const;
    void                getQuillVEdges( list<VEdge2D*>& quillEdge ) const;


    bool    isInfinite() const;
    bool    isBounded() const;
    bool    isUnBounded() const;

    bool    isAdjacentTo( const VFace2D* const face ) const;
    bool    isIncidentTo( const VEdge2D* const edge ) const;
    bool    isIncidentTo( const VVertex2D* const vertex ) const;
};


} // GeometryTier
} // V

#endif


