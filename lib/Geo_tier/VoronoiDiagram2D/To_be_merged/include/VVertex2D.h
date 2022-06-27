#ifndef VVERTEX2D_H
#define VVERTEX2D_H

#include "rg_Point2D.h"



#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class VEdge2D;
class VFace2D;
class Generator2D;

enum STATUS_OF_VVERTEX
{
    RED_V,        // vertex in the red ocean
    BLUE_V,       // vertex adjacent to red ocean
    PURPLE_V,     // vertex on the boundary of red ocean
    YELLOW_V,     // fictitious vertex (breakwater)
    WHITE_V       // default color
};

class VVertex2D
{
private:
    int                 m_ID;
    VEdge2D*            m_firstVEdge;
    rg_Point2D          m_location;
    STATUS_OF_VVERTEX   m_status;
    bool                m_isFictitious;

public:
    VVertex2D();
    VVertex2D( const int& ID );
    VVertex2D( const int& ID, const rg_Point2D& location );
    VVertex2D( const int& ID, const rg_Point2D& location, VEdge2D* const firstEdge );
    VVertex2D( const VVertex2D& vertex );
    ~VVertex2D(void);

    inline int               getID() const { return m_ID; }
    inline rg_Point2D        getLocation() const { return m_location; }
    inline VEdge2D*          getFirstVEdge() const { return m_firstVEdge; }
    inline STATUS_OF_VVERTEX getStatus() const { return m_status; }

    inline void              setID( const int& ID ) { m_ID = ID; }
    inline void              setLocation( const rg_Point2D& location ) { m_location = location; }
    inline void              setFirstVEdge( VEdge2D* const edge ) { m_firstVEdge = edge; }
    inline void              setStatus( const STATUS_OF_VVERTEX& status ) { m_status = status; }
    inline void              setFictitious() { m_isFictitious = true; }
    inline void              setRealVertex() { m_isFictitious = false; }

    VVertex2D&               operator=( const VVertex2D& vertex );
    bool                     operator==( const VVertex2D& vertex ) const;


    bool                     isFictitious() const { return m_isFictitious; }
    void                     getIncident3VEdges( list<VEdge2D*>& incidentVEdgesList ) const;
    void                     getIncident3VFaces( list<VFace2D*>& incidentVFacesList ) const;
    void                     getAdjacent3VVertices( list<VVertex2D*>& adjacentVVerticesList ) const;
    void                     getDefining3Generators( list<Generator2D*>& definingGeneratorsList ) const;
    double                   computeMUValue( Generator2D* const newGenerator );

    bool    isMemberOf( const VFace2D* const face ) const;
    bool    isMemberOf( const VEdge2D* const edge ) const;
    bool    isAdjacentTo( const VVertex2D* const vertex ) const;
    bool    isInfinite() const;
    double  computeRadiusOfTangentCircle() const;

};

}
}

#endif


