#ifndef POLYGON_VD_H
#define POLYGON_VD_H

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

//SH
//outside -> exterior ? function in function
//inside -> interiror ? function in function
// Priority of status: Infinite >> Exterior >> Interior
// Definition of status
// - Infinite: From infinite voroVertex
// - Exterior: Right side of polygon boundary (direction, loop?)
// - Interior: Left  "
//SH

//#define DEBUG_VERTEX

namespace V {
namespace GeometryTier {

class PolygonVD2D : public VoronoiDiagram2DC
{
public:
    enum CHILDREN_DISK_GENERATION_METHOD {FEWER_OD, UNIFORM_OD};
    enum VDEntity_Location_Status { OUTSIDE_POLYGON, INSIDE_POLYGON, ON_POLYGON_BOUNDARY, INFINITE_ENTITY, UNKNOWN_LOCATION };
    typedef map< VVertex2D*, VDEntity_Location_Status > VVertex_Location_Status_Map;
    typedef map< VEdge2D*,VDEntity_Location_Status >    VEdge_Location_Status_Map;

#ifdef CHECK_COMP_TIME
    double t_wave_propagation;
    double t_anomaly_check;
    double t_find_red_neighbor;
    double t_make_new_cell;
    double t_mark_in;
    double t_connect_topology;
    double t_merge_split;
    double t_comp_vertex;
    double t_others;

    int    num_DDD;
    int    num_DDL;
    int    num_DLV;
    int    num_DLL;
    double t_comp_V_DDD;
    double t_comp_V_DDL;
    double t_comp_V_DLV;
    double t_comp_V_DLL;

    double t_check_parents;
    double t_check_orientation;

#endif

protected:
    Polygon2D*                              m_rectContainer;
    list<Polygon2D>                         m_polygon;
    Polygon2D*                              m_ptr_last_inserted_polygon;

    map< Generator2D*, Polygon2D* >         m_mapGenerator2Polygon;
    map< Polygon2D*, list<Generator2D*> >   m_mapPolygon2Generator;

    VVertex_Location_Status_Map             m_VVertexToLocationStatus;
    VEdge_Location_Status_Map               m_VEdgeToLocationStatus;

    /*
    This map (m_map_DiskID2VFace) is used for disk packing in polygon.
    We supposed that all generators have unique ID.
    */
    map< int, VFace2D* >                    m_map_DiskID2VFace;
    list<rg_Circle2D>                       m_insertedDisks;

    bool                                    b_there_is_container;
    bool                                    b_VD_is_constructed;

public:
    PolygonVD2D();
    virtual ~PolygonVD2D();
#ifdef CHECK_COMP_TIME
    void init_();
    void print_();
#endif
    void constructVoronoiDiagram( const Polygon2D& polygon );
    void constructVoronoiDiagram( const list<Polygon2D>& polygons );
    void constructVoronoiDiagramOfRectangle( const rg_Point2D& minPt, const rg_Point2D& maxPt );
    void identify_internal_and_on_and_external_Vvertices_and_Vedges();
    void identify_internal_and_on_and_external_and_infinite_Vvertices_and_Vedges();

    void updateEdge( VEdge2D* target );

//private:
    void preprocess_for_construction( const Polygon2D& polygon, list<VertexGenerator2D*>& vertexGens, list<EdgeGenerator2D*>& edgeGens );
    void preprocess_for_construction( const list<Polygon2D>& polygons, list<VertexGenerator2D*>& vertexGens, list<EdgeGenerator2D*>& edgeGens );
    void replicatePolygonEntitiesAsVDGenerators( Polygon2D* const polygon, list<VertexGenerator2D*>& vertexGens, list<EdgeGenerator2D*>& edgeGens );
    void constructPhantomVoronoiDiagram();
    void createPhantomGenerators();
    void insertGeneratorsToPhantomVoronoiDiagram( const list<VertexGenerator2D*>& vertexGenerators, const list<EdgeGenerator2D*>& edgeGenerators );

    // The following function (collectVertex_N_EdgeGenerators) is not used in the latest version that use sets of VertexGenerator2D and EdgeGenerator2D as input arguments.
    void collectVertex_N_EdgeGenerators( Polygon2D* const polygon, list<VertexGenerator2D*>& vertexGenerators, list<EdgeGenerator2D*>& edgeGenerators );
    void insertVertexGeneratorSet( const list<VertexGenerator2D*>& vertexGenerators );
    void insertEdgeGeneratorSet( const list<EdgeGenerator2D*>& edgeGenerators );

    void updateVD_with_vertexGenerator( const VertexGenerator2D* const newVertexGen );
    void updateVD_with_edgeGenerator( const EdgeGenerator2D* const newEdgeGen );

    VVertex2D* findSeedRedVertex( Generator2D* const newGen, VFace2D* const anchorCell );
    double computeMUValue( VVertex2D* const vertex, Generator2D* const newGenerator );
    double computeMUValueFromEdgeGenerator( VVertex2D* const vertex, EdgeGenerator2D* const newEdgeGen );

    bool isAnomalizingEdge( VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    bool doAnomalyEdgeTest_for_DiskGenerator( VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    bool doAnomalyEdgeTest_for_EdgeGenerator( VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    bool doAnomalyEdgeTest_for_VertexGenerator( VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    bool areTwoPointsOnThisEdge( VEdge2D* edge, const rg_Point2D& pt1, const rg_Point2D& pt2 );


    void count_number_of_types_of_defining_generators( const VVertex2D* const vertex, int& numDisks, int& numLineSegs, int& numPoints );
    void computeCoordOfNewVertices( list<VVertex2D*>& newVertices );
    void computeCoordOf_A_newVertex( VVertex2D* const newVertex );
    int  computeTangentCircles_DDD( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_DDE( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_DDV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_DEE( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_DEV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_DVV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_EEE( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_EEV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_EVV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeTangentCircles_VVV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );

    int computeTangentCircles_of_two_disks_which_have_same_radii_and_a_line( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int computeTangentCircles_of_two_disks_and_a_line( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int computeTangentCircles_of_a_disk_and_two_lines( const rg_Circle2D& disk, const rg_Line2D& line1, const rg_Line2D& line2, vector<rg_Circle2D>& tangentCircles );
    int computeTangentCircles_of_two_disks_and_a_line_for_coordinate( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, vector<rg_Circle2D>& tangentCircles );
    int computeTangentCircles_of_two_disks_and_an_oriented_line( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    rg_Circle2D chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( const int& numTangentCircles, const rg_Circle2D& tangentCircle1, const rg_Circle2D& tangentCircle2, VVertex2D* const newVVertex );






    void classifyVertexGenerators();








    Polygon2D* getLastInsertedPolygon() { return m_ptr_last_inserted_polygon; };






public:
    // For inside-only packing
    // HAVE TO DO
    // In the function 'insertThisDiskInPolygonVoronoiDiagram', check if the input disk intersect with the polygon.
    // If there is no intersection, check the status of seed red vertex.
    // After we get new vertices, mark new vertices the same status with one of seed red vertex.
    DiskGenerator2D* insertThisDiskInPolygonVoronoiDiagram( const rg_Circle2D& disk );
    DiskGenerator2D* insertThisDiskInPolygonVoronoiDiagram( const rg_Circle2D& disk, list<VVertex2D*>& redVertices, list<VVertex2D*>& newVertices );
    DiskGenerator2D* insert_a_shrinked_disk_on_this_VVertex( const VVertex2D* const vertex, const double& shrink_ratio );
    DiskGenerator2D* insert_a_shrinked_disk_on_this_VVertex( const VVertex2D* const vertex, const double& shrink_ratio, list<VVertex2D*>& redVertices, list<VVertex2D*>& newVertices );
    DiskGenerator2D* insert_a_shrinked_disk_on_this_VVertex_outside_of_polygon( const VVertex2D* const vertex, const double& shrink_ratio, list<VVertex2D*>& redVertices, list<VVertex2D*>& newVertices );

    PolygonVD2D::VDEntity_Location_Status get_location_status_of_vvertex(const VVertex2D* const vertex);
    PolygonVD2D::VDEntity_Location_Status get_location_status_of_edge(const VEdge2D* const edge);

    void get_disk_generators_which_share_VFace_with_polygon_boundary(list<DiskGenerator2D*>& diskGenerators);
    void get_edge_generators(list<EdgeGenerator2D*>& edgeGenerators);
    void get_vertex_generators(list<VertexGenerator2D*>& vertexGenerators);

	void get_Voronoi_edges_inside_Polygon(list<const VEdge2D*>& VEdgesList);
    void get_Voronoi_edges_outside_Polygon( list<const VEdge2D*>& VEdgesList );
    
    rg_RQBzCurve2D get_geometry_of_edge(const VEdge2D* const edge);
    void           get_sampling_points_of_edge(const VEdge2D* const edge, const int& numSamplingPts, list<rg_Point2D>& samplingPts);

    DiskGenerator2D*    getGeneratorWhichHasThisID( const int& diskID );
    void                makeDiskID2VFaceMap();
    void                remove_red_entities();

protected:
        void set_type_of_PVertex(VertexGenerator2D* vertexGen, EdgeGenerator2D* prevEdgeGen, EdgeGenerator2D* nextEdgeGen);
        void find_two_VEdges_between_vertex_generator_N_incident_edge_generators(VertexGenerator2D *& vertexGen, EdgeGenerator2D *& prevEdgeGen, EdgeGenerator2D *& nextEdgeGen, VEdge2D*& edge_between_prevEdgeGen_N_vertexGen, VEdge2D*& edge_between_vertexGen_N_nextEdgeGen);
        void find_and_relocate_two_boundary_Vvertices_to_vertex_generator_center(VertexGenerator2D *& vertexGen, EdgeGenerator2D *& prevEdgeGen, EdgeGenerator2D *& nextEdgeGen);
    

    //void identify_internal_and_on_and_external_Vvertices_and_Vedges();
        void set_location_status_of_infinite_Vvertices_and_Vedges_with_outside(stack<VEdge2D*>& VEdgeStack);
        //void identify_Vvertices_on_boundary();
        void identify_Vvertices_on_boundary( list< pair< VertexGenerator2D*, VVertex2D* > >& PVertex_N_onVVertexPair );
        void identify_infinite_Vvertices_and_Vedges();
        void identify_outside_Vvertices_and_Vedges( const list< pair< VertexGenerator2D*, VVertex2D* > >& PVertex_N_onVVertexPair );
        void identify_and_set_location_status_of_both_outside_and_on_Vvertices_and_outside_Vedges(stack<VEdge2D*>& VEdgeStack);
        void set_location_status_of_remaining_Vvertices_and_Vedges_with_inside();
    void refine_location_of_interior_VVertices();
        PolygonVD2D::VDEntity_Location_Status get_location_status_of_Vvertex(const VVertex2D* const vertex);
        PolygonVD2D::VDEntity_Location_Status get_location_status_of_Vedge(const VEdge2D* const edge);
        rg_Circle2D refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle(const VVertex2D* const vertex);
        void count_number_of_occurred_generator_types_and_get_parent_generators_at_VVertex(const VVertex2D* const vertex, int& numberOfDisks, int& numberOfLineSegs, int& numberOfPoints, Generator2D* threeParentGenerators[]);
        void orthogonalize_VEdge_by_projecting_VVertex(VVertex2D* vertex, EdgeGenerator2D* edgeGen, const rg_Point2D& passingPt);
    void adjust_topology_of_interior_of_polygon_by_edge_flipping();
        void flatten_topology_of_VFace_boundary_of_reflex_vertex_gen_with_possible_flipping_out_of_VEdges(VertexGenerator2D* vertexGenerator);
        void flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges(EdgeGenerator2D* edgeGenerator);
        void flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges_exterior( EdgeGenerator2D* edgeGenerator );
        void find_interior_VEdge_chain_of_this_PEdge_CCW(EdgeGenerator2D* edgeGenerator, list<VEdge2D*>& VEdgeChainCCW);
        void find_exterior_VEdge_chain_of_this_PEdge_CCW( EdgeGenerator2D* edgeGenerator, list<VEdge2D*>& VEdgeChainCCW );
        void refine_VVertices_coordinates_of_flipped_VEdge(VEdge2D* edge);

    rg_RQBzCurve2D createRQBzCurveOfParabola(const rg_Circle2D& disk, const rg_Line2D& line);



    /////////////////////////////////////////////////////////////////
    //
    // For insertion
    //
    /////////////////////////////////////////////////////////////////

    DiskGenerator2D* createGeneratorOfInsertedDisk(const rg_Circle2D& disk);
    VFace2D* findAnchorCell(const rg_Point2D& pt);
    void makeNewCellAndConnectToCurrVD( DiskGenerator2D* newGenerator, list<VVertex2D*>& redVVertices, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& newVVertices, list<VEdge2D*>& newVEdges );
    VVertex2D* findSeedRedVertex( DiskGenerator2D* const newGenerator, Generator2D* const closestGenerator );
    VVertex2D* splitVEdgeAtFictitiousVVertex( VEdge2D* const currVEdge );
    bool isAnomalizingEdge_for_DiskGenerator(  VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    bool isAnomalizingEdge_parabola ( VEdge2D* const incidentEdge, const Parabola2D& parabola, Generator2D* const newGenerator, const double& radiusOfFocus );
    bool isAnomalizingEdge_hyperbola ( VEdge2D* const incidentEdge, const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Circle2D& disk3 );
    bool isAnomalizingEdge_line( const rg_Line2D& lineSeg, const rg_Circle2D& circumcircle_sp, const rg_Circle2D& circumcircle_ep, Generator2D* const newGenerator );

    void markLocationStatusAtNewEntities( list<VVertex2D*>& newVertices, list<VEdge2D*>& newEdges );
    void markLocationStatusAtNewEntities_outside( list<VVertex2D*>& newVertices, list<VEdge2D*>& newEdges );

    void computeCoordOfNewVVertices( list<VVertex2D*>& newVVertices );
    void computeCoordOfNewVVertices_outside( list<VVertex2D*>& newVVertices );
    void computeCoordOfNewVVertex_on_parabola( const Parabola2D& parabola, const rg_Circle2D& disk, const double& radiusOfFocus, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    int  computeCoordOfNewVVertex_on_line_of_two_polygon_edges( const rg_Line2D& geometry_of_edge, const rg_Circle2D& circumcircle_sp, const rg_Circle2D& circumcircle_ep, const rg_Circle2D& disk, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 );
    void computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( const rg_Point2D& SP, const rg_Point2D& dirVec, const rg_Circle2D& disk, rg_Circle2D& tangentCircle );
    
    bool findIntersectionPtWithinThisArrange( const rg_Line2D& line, const rg_Circle2D& disk_base, const rg_Circle2D& disk_another, const rg_Point2D& SP, const rg_Point2D& EP, rg_Circle2D& tangentCircle );

    bool this_circumcircle_has_right_orientation( const VVertex2D* const vertex, const rg_Circle2D& circumcircle );
    


    /////////////////////////////////////////////////////////////////
    //
    // For geometry
    //
    /////////////////////////////////////////////////////////////////

    rg_Point2D getTangentVector( const rg_Point2D& pt, const rg_Point2D& leftCenter, const rg_Point2D& rightCenter );
    rg_Point2D getPassingPt ( const rg_Circle2D& leftDisk, const rg_Circle2D& rightDisk );
};

inline  int     PolygonVD2D::computeTangentCircles_DDV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 ) { return computeTangentCircles_DDD( newVertex, tangentCircle1, tangentCircle2 ); }
inline  int     PolygonVD2D::computeTangentCircles_DVV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 ) { return computeTangentCircles_DDD( newVertex, tangentCircle1, tangentCircle2 ); }
inline  int     PolygonVD2D::computeTangentCircles_VVV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 ) { return computeTangentCircles_DDD( newVertex, tangentCircle1, tangentCircle2 ); }



} // GeometryTier
} // V

#endif

