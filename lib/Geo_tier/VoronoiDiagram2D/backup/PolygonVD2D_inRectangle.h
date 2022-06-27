#ifndef POLYGON_VD_IN_RECTANGLE_H
#define POLYGON_VD_IN_RECTANGLE_H

#include "PolygonVD2D.h"
#include "VoronoiDiagram2DC.h"
#include "Polygon2D.h"
#include "Generator2D.h"
#include "DiskGenerator2D.h"
#include "VertexGenerator2D.h"
#include "EdgeGenerator2D.h"
#include "Parabola2D.h"
#include "rg_RQBzCurve2D.h"
#include "rg_IntersectFunc.h"
#include <list>
#include <stack>
#include <queue>
#include <map>
#include <set>
#include <functional>

#include <cmath>

//#define DEBUG_VERTEX

namespace V {
namespace GeometryTier {


class PolygonVD2D_inRectangle : public PolygonVD2D
{
public:

private:
    Polygon2D m_rectContainer;
    list<Polygon2D> m_polygons;
    map<Generator2D*, Polygon2D*> m_mapGenerator2Polygon;


public:
    PolygonVD2D_inRectangle();
    ~PolygonVD2D_inRectangle();

    void constructVoronoiDiagram_rectContainer( const Polygon2D& polygon, const CHILDREN_DISK_GENERATION_METHOD& method = UNIFORM_OD );

    Polygon2D getRectContainer() const;
    void      setRectContainer( const Polygon2D& polygon );
    void constructVoronoiDiagramOfRectContainer();
    void insertOnePolygonIntoVoronoiDiagram( const Polygon2D& polygon );

private:
    
    void constructVoronoiDiagramOfRectContainer( const Polygon2D& inputPolygon );

    void insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer();

    // For construction of polygon Voronoi diagram in rectangular container
    DiskGenerator2D* insertOneChildDiskInPolygonVoronoiDiagram( const rg_Triplet<rg_Circle2D, int, void*>& diskIDUserData );

    //void identify_internal_and_on_and_external_Vvertices_and_Vedges();
    void identify_Vvertices_on_boundary( list< pair< VertexGenerator2D*, VVertex2D* > >& PVertex_N_onVVertexPair );
    void adjust_topology_of_interior_and_exterior_of_polygon_by_edge_flipping();
    void computeCoordOfNewVVertices_before_merging( list<VVertex2D*>& newVVertices );
    void refine_location_of_VVertices();
    void computeCoordOfNewVVertex( VVertex2D*& vertex );


    void replicatePolygonEntitiesAsVDGenerators( list<Generator2D*>& newGenerators );
    void generateChildrenDisksOfPolygonGenerators( const list<Generator2D*>& newGenerators );
    void insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer( const list<Generator2D*>& newGenerators );
    void mergeVFacesOfChildrenDisksOnEachPEdge( const list<Generator2D*>& newGenerators );
    void relocate_some_boundary_Vvertices_to_vertex_generator_centers( const list<Generator2D*>& newGenerators );

    void secure_VFace_topology_of_PVertices_by_flipping_the_mating_quill_VEdge_of_relocated_VVertex( const list<Generator2D*>& newGenerators );
    void connect_parent_generator_to_first_son_disk_generator( const list<Generator2D*>& newGenerators );
    void identify_internal_and_on_and_external_Vvertices_and_Vedges( const list<Generator2D*>& newGenerators );
    void identify_Vvertices_on_boundary( const list<Generator2D*>& newGenerators, list< pair< VertexGenerator2D*, VVertex2D* > >& PVertex_N_onVVertexPair );
    void refine_location_of_VVertices( const list<Generator2D*>& newGenerators );
    void adjust_topology_of_interior_and_exterior_of_polygon_by_edge_flipping( const list<Generator2D*>& newGenerators );
};


inline Polygon2D PolygonVD2D_inRectangle::getRectContainer() const { return m_rectContainer; }
inline void PolygonVD2D_inRectangle::setRectContainer( const Polygon2D& polygon ) { m_rectContainer = polygon; }

} // GeometryTier
} // V


#endif

