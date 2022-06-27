#ifndef VVERTEX2DP_H
#define VVERTEX2DP_H

#include "rg_Point2D.h"


#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class VEdge2DP;
class VFace2DP;
class Generator2DP;

enum STATUS_OF_VVERTEXP
{
    RED_V_P,        // vertex in the 
    BLUE_V_P,       // vertex adjacent 
    WHITE_V_P
};

class VVertex2DP
{
private:
    int                 m_ID;
    VEdge2DP*           m_firstVEdge;
    rg_Point2D          m_location;
    STATUS_OF_VVERTEXP  m_status;

public:
    VVertex2DP();
    VVertex2DP( const int& ID );
    VVertex2DP( const int& ID, const rg_Point2D& location );
    VVertex2DP( const int& ID, const rg_Point2D& location, VEdge2DP* const firstEdge );
    VVertex2DP( const VVertex2DP& vertex );
    ~VVertex2DP(void);

    inline int                  getID() const { return m_ID; }
    inline rg_Point2D           getLocation() const { return m_location; }
    inline VEdge2DP*            getFirstVEdge() const { return m_firstVEdge; }
    inline STATUS_OF_VVERTEXP   getStatus() const { return m_status; }

    inline void                 setID( const int& ID ) { m_ID = ID; }
    inline void                 setLocation( const rg_Point2D& location ) { m_location = location; }
    inline void                 setFirstVEdge( VEdge2DP* const edge ) { m_firstVEdge = edge; }
    inline void                 setStatus( const STATUS_OF_VVERTEXP& status ) { m_status = status; }

    VVertex2DP&                 operator=( const VVertex2DP& vertex );

    void                        getIncident3VEdges( list<VEdge2DP*>& incidentVEdgesList ) const;
    void                        getIncident3VFaces( list<VFace2DP*>& incidentVFacesList ) const;
    void                        getAdjacent3VVertices( list<VVertex2DP*>& adjacentVVerticesList ) const;
    void                        getDefining3Generators( list<Generator2DP*>& definingGeneratorsList ) const;
    double                      computeHValue( Generator2DP* const newGenerator );

    double                      computeRadiusOfTangentCircle() const;
    bool                        isInfinite() const;
    bool                        isMemberOf( const VFace2DP* const face ) const;
    bool                        isMemberOf( const VEdge2DP* const edge ) const;
    bool                        isAdjacentTo( const VVertex2DP* const vertex ) const;



};

} // namespace GeometryTier
} // namespace BULL2D


#endif


