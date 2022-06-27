#include "PolygonVD2D_inRectangle.h"
using namespace V::GeometryTier;

#include "rg_GeoFunc.h"
#include "Parabola2D.h"
#include "rg_TMatrix2D.h"
#include "rg_dList.h"

#if defined( DEBUG_VERTEX ) || defined ( CHECK_COMP_TIME )
#include <fstream>
#include <time.h>
using namespace std;
#endif


PolygonVD2D_inRectangle::PolygonVD2D_inRectangle()
{
    m_childDiskVD = rg_NULL;
}


PolygonVD2D_inRectangle::~PolygonVD2D_inRectangle()
{
    // commented by Joonghyun (already destructed by base class) 

    //for (list<Generator2D*>::iterator i_polygonGen = m_polygonGenerators.begin(); i_polygonGen != m_polygonGenerators.end(); ++i_polygonGen) {
    //    Generator2D* polygonGen = *i_polygonGen;
    //    delete polygonGen;
    //}

    //if ( m_childDiskVD != rg_NULL )
    //    delete m_childDiskVD;
}



void PolygonVD2D_inRectangle::constructVoronoiDiagram_rectContainer( const Polygon2D& polygon, const CHILDREN_DISK_GENERATION_METHOD& method )
{
    // construct VD of rectangle container
    constructVoronoiDiagramOfRectContainer( polygon );

    m_polygon = polygon;
    insertOnePolygonIntoVoronoiDiagram( polygon );
    /*
    m_polygon = polygon;
    list<Generator2D*> newPolygonGenerators;
    replicatePolygonEntitiesAsVDGenerators( newPolygonGenerators );

    generate_children_disks_of_polygon_generator(); //insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer

    insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer();

    //return;

    mergeVFacesOfChildrenDisksOnEachPEdge();
    
    relocate_some_boundary_Vvertices_to_vertex_generator_centers();
    
    // Children disk를 'as big as'로 만들었을 때, vertex generator가 anomaly가 될 수 있다.
    // vertex generator의 face는 topologically 삼각형 이상이어야 하기 때문에 다림질 하기 전에 삼각형을 만들어 준다.
    secure_VFace_topology_of_PVertices_by_flipping_the_mating_quill_VEdge_of_relocated_VVertex();


    // polygonVD에 있는 각 VertexGenerator와 EdgeGenerator 마다 childrenVD의 DiskGenerator 들 중 해당되는 Generator로 연결해준다.
    connect_parent_generator_to_first_son_disk_generator();

    identify_internal_and_on_and_external_Vvertices_and_Vedges();
    
    refine_location_of_VVertices();
    
    // VFace (즉, VEdge loop)를 다림질한다.
    adjust_topology_of_interior_and_exterior_of_polygon_by_edge_flipping();
    */
}



void PolygonVD2D_inRectangle::constructVoronoiDiagramOfRectContainer()
{
    list<rg_Point2D> exteriorBoundaryVertices;
    m_rectContainer.get_exterior_boundary_vertices( exteriorBoundaryVertices );
    rg_BoundingBox2D boundingBox;
    m_rectContainer.get_bounding_box( boundingBox );
    rg_Point2D minPt = boundingBox.getMinPt();
    rg_Point2D maxPt = boundingBox.getMaxPt();


    // create polygon generators for rectangle
    double minX = minPt.getX();
    double minY = minPt.getY();
    double maxX = maxPt.getX();
    double maxY = maxPt.getY();

    //rg_Point2D          boundaryPt[5] = { rg_Point2D( minX, maxY ), rg_Point2D( minX, minY ), rg_Point2D( maxX, minY ), rg_Point2D( maxX, maxY ), rg_Point2D( minX, maxY ) };
    rg_Point2D          boundaryPt[5] = { rg_Point2D( minX, maxY ), rg_Point2D( maxX, maxY ), rg_Point2D( maxX, minY ), rg_Point2D( minX, minY ), rg_Point2D( minX, maxY ) };
    VertexGenerator2D* vertexGenArray[5];
    EdgeGenerator2D* edgeGenArray[4];

    for ( int i = 0; i < 4; ++i ) {
        VertexGenerator2D* vertexGen = new VertexGenerator2D( boundaryPt[i], i * 2 );
        vertexGen->isThisFromContainer( true );
        m_polygonGenerators.push_back( vertexGen );
        vertexGenArray[i] = vertexGen;
        m_mapGenerator2Polygon[vertexGen] = &m_rectContainer;

        EdgeGenerator2D* edgeGen = new EdgeGenerator2D( boundaryPt[i], boundaryPt[i + 1], i * 2 + 1 );
        edgeGen->isThisFromContainer( true );
        m_polygonGenerators.push_back( edgeGen );
        edgeGenArray[i] = edgeGen;
        m_mapGenerator2Polygon[edgeGen] = &m_rectContainer;
    }
    vertexGenArray[4] = vertexGenArray[0];


    for ( int i = 0; i < 4; ++i ) {
        vertexGenArray[i]->set_next_edge_generator( edgeGenArray[i] );
        vertexGenArray[i + 1]->set_previous_edge_generator( edgeGenArray[i] );

        edgeGenArray[i]->set_start_vertex_generator( vertexGenArray[i] );
        edgeGenArray[i]->set_end_vertex_generator( vertexGenArray[i + 1] );
    }

    // create entities
    //  1. make elements of voronoi diagram of phantom generators

    VVertex2D* vertices[14] = {
        createVertex( 0 ),    // LT
        createVertex( 1 ),    // LB
        createVertex( 2 ),    // RB
        createVertex( 3 ),    // RT
        createVertex( 4 ),    // IN_1
        createVertex( 5 ),    // IN_2
        createVertex( 6 ),    // LTT
        createVertex( 7 ),    // LTL
        createVertex( 8 ),    // LBL
        createVertex( 9 ),    // LBB
        createVertex( 10 ),   // RBB
        createVertex( 11 ),   // RBR
        createVertex( 12 ),   // RTR
        createVertex( 13 ) };  // RTT

    VEdge2D* edges[21] = {
        createEdge( 0 ),        // INF_LT
        createEdge( 1 ),        // INF_L
        createEdge( 2 ),        // INF_LB
        createEdge( 3 ),        // INF_B
        createEdge( 4 ),        // INF_RB
        createEdge( 5 ),        // INF_R
        createEdge( 6 ),        // INF_RT
        createEdge( 7 ),        // INF_T
        createEdge( 8 ),        // UBND_TL
        createEdge( 9 ),        // UBND_LT
        createEdge( 10 ),       // UBND_LB
        createEdge( 11 ),       // UBND_BL
        createEdge( 12 ),       // UBND_BR
        createEdge( 13 ),       // UBND_RB
        createEdge( 14 ),       // UBND_RT
        createEdge( 15 ),       // UBND_TR
        createEdge( 16 ),       // IN_LT
        createEdge( 17 ),       // IN_LB
        createEdge( 18 ),       // IN_CENTER
        createEdge( 19 ),       // IN_RB
        createEdge( 20 ) };      // IN_RT

    VFace2D* faces[9] = {
        createFace( -1 ),       // INF
        createFace( 0 ),        // VTX_LT
        createFace( 1 ),        // EDGE_L
        createFace( 2 ),        // VTX_LB
        createFace( 3 ),        // EDGE_B
        createFace( 4 ),        // VTX_RB
        createFace( 5 ),        // EDGE_R
        createFace( 6 ),        // VTX_RT
        createFace( 7 ) };       // EDGE_T


    // IN_OUT marking.....
    m_VVertexToLocationStatus[vertices[0]] = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[1]] = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[2]] = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[3]] = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[4]] = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[5]] = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[6]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[7]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[8]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[9]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[10]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[11]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[12]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[13]] = INFINITE_ENTITY;

    m_VEdgeToLocationStatus[edges[0]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[1]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[2]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[3]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[4]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[5]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[6]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[7]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[8]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[9]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[10]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[11]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[12]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[13]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[14]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[15]] = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[16]] = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[17]] = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[18]] = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[19]] = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[20]] = OUTSIDE_POLYGON;

    // topology
    //vertexGenArray[0]->setOuterFace( faces[1] );
    //vertexGenArray[1]->setOuterFace( faces[3] );
    //vertexGenArray[2]->setOuterFace( faces[5] );
    //vertexGenArray[3]->setOuterFace( faces[7] );
    vertexGenArray[0]->setOuterFace( faces[1] );
    vertexGenArray[1]->setOuterFace( faces[7] );
    vertexGenArray[2]->setOuterFace( faces[5] );
    vertexGenArray[3]->setOuterFace( faces[3] );

    vertexGenArray[0]->setUserData( vertexGenArray[0] );
    vertexGenArray[1]->setUserData( vertexGenArray[1] );
    vertexGenArray[2]->setUserData( vertexGenArray[2] );
    vertexGenArray[3]->setUserData( vertexGenArray[3] );

    //edgeGenArray[0]->setOuterFace( faces[2] );
    //edgeGenArray[1]->setOuterFace( faces[4] );
    //edgeGenArray[2]->setOuterFace( faces[6] );
    //edgeGenArray[3]->setOuterFace( faces[8] );
    edgeGenArray[0]->setOuterFace( faces[8] );
    edgeGenArray[1]->setOuterFace( faces[6] );
    edgeGenArray[2]->setOuterFace( faces[4] );
    edgeGenArray[3]->setOuterFace( faces[2] );

    edgeGenArray[0]->setUserData( edgeGenArray[0] );
    edgeGenArray[1]->setUserData( edgeGenArray[1] );
    edgeGenArray[2]->setUserData( edgeGenArray[2] );
    edgeGenArray[3]->setUserData( edgeGenArray[3] );

    //faces[1]->setGenerator( (Generator2D*)vertexGenArray[0] );
    //faces[3]->setGenerator( (Generator2D*)vertexGenArray[1] );
    //faces[5]->setGenerator( (Generator2D*)vertexGenArray[2] );
    //faces[7]->setGenerator( (Generator2D*)vertexGenArray[3] );
    //faces[2]->setGenerator( (Generator2D*)edgeGenArray[0] );
    //faces[4]->setGenerator( (Generator2D*)edgeGenArray[1] );
    //faces[6]->setGenerator( (Generator2D*)edgeGenArray[2] );
    //faces[8]->setGenerator( (Generator2D*)edgeGenArray[3] );
    faces[1]->setGenerator( (Generator2D*)vertexGenArray[0] );
    faces[3]->setGenerator( (Generator2D*)vertexGenArray[3] );
    faces[5]->setGenerator( (Generator2D*)vertexGenArray[2] );
    faces[7]->setGenerator( (Generator2D*)vertexGenArray[1] );
    faces[2]->setGenerator( (Generator2D*)edgeGenArray[3] );
    faces[4]->setGenerator( (Generator2D*)edgeGenArray[2] );
    faces[6]->setGenerator( (Generator2D*)edgeGenArray[1] );
    faces[8]->setGenerator( (Generator2D*)edgeGenArray[0] );


    vertices[0]->setFirstVEdge( edges[16] );
    vertices[1]->setFirstVEdge( edges[17] );
    vertices[2]->setFirstVEdge( edges[19] );
    vertices[3]->setFirstVEdge( edges[20] );
    vertices[4]->setFirstVEdge( edges[18] );
    vertices[5]->setFirstVEdge( edges[18] );
    vertices[6]->setFirstVEdge( edges[0] );
    vertices[7]->setFirstVEdge( edges[1] );
    vertices[8]->setFirstVEdge( edges[2] );
    vertices[9]->setFirstVEdge( edges[3] );
    vertices[10]->setFirstVEdge( edges[4] );
    vertices[11]->setFirstVEdge( edges[5] );
    vertices[12]->setFirstVEdge( edges[6] );
    vertices[13]->setFirstVEdge( edges[7] );

    edges[0]->setTopology( vertices[6], vertices[7], faces[1], faces[0], edges[9], edges[1], edges[8], edges[7] );
    edges[1]->setTopology( vertices[7], vertices[8], faces[2], faces[0], edges[10], edges[2], edges[9], edges[0] );
    edges[2]->setTopology( vertices[8], vertices[9], faces[3], faces[0], edges[11], edges[3], edges[10], edges[1] );
    edges[3]->setTopology( vertices[9], vertices[10], faces[4], faces[0], edges[12], edges[4], edges[11], edges[2] );
    edges[4]->setTopology( vertices[10], vertices[11], faces[5], faces[0], edges[13], edges[5], edges[12], edges[3] );
    edges[5]->setTopology( vertices[11], vertices[12], faces[6], faces[0], edges[14], edges[6], edges[13], edges[4] );
    edges[6]->setTopology( vertices[12], vertices[13], faces[7], faces[0], edges[15], edges[7], edges[14], edges[5] );
    edges[7]->setTopology( vertices[13], vertices[6], faces[8], faces[0], edges[8], edges[0], edges[15], edges[6] );

    edges[8]->setTopology( vertices[0], vertices[6], faces[1], faces[8], edges[0], edges[7], edges[9], edges[16] );
    edges[9]->setTopology( vertices[0], vertices[7], faces[2], faces[1], edges[1], edges[0], edges[16], edges[8] );
    edges[10]->setTopology( vertices[1], vertices[8], faces[3], faces[2], edges[2], edges[1], edges[11], edges[17] );
    edges[11]->setTopology( vertices[1], vertices[9], faces[4], faces[3], edges[3], edges[2], edges[17], edges[10] );
    edges[12]->setTopology( vertices[2], vertices[10], faces[5], faces[4], edges[4], edges[3], edges[13], edges[19] );
    edges[13]->setTopology( vertices[2], vertices[11], faces[6], faces[5], edges[5], edges[4], edges[19], edges[12] );
    edges[14]->setTopology( vertices[3], vertices[12], faces[7], faces[6], edges[6], edges[5], edges[15], edges[20] );
    edges[15]->setTopology( vertices[3], vertices[13], faces[8], faces[7], edges[7], edges[6], edges[20], edges[14] );

    bool b_horizontal_rectangle = true;
    double diffX = maxX - minX;
    double diffY = maxY - minY;
    if ( diffX < diffY )
        b_horizontal_rectangle = false;

    if ( b_horizontal_rectangle ) {
        edges[16]->setTopology( vertices[0], vertices[4], faces[8], faces[2], edges[18], edges[17], edges[8], edges[9] );
        edges[17]->setTopology( vertices[1], vertices[4], faces[2], faces[4], edges[16], edges[18], edges[10], edges[11] );
        edges[18]->setTopology( vertices[4], vertices[5], faces[8], faces[4], edges[20], edges[19], edges[16], edges[17] );
        edges[19]->setTopology( vertices[5], vertices[2], faces[6], faces[4], edges[13], edges[12], edges[20], edges[18] );
        edges[20]->setTopology( vertices[5], vertices[3], faces[8], faces[6], edges[15], edges[14], edges[18], edges[19] );
    }
    else {
        edges[16]->setTopology( vertices[0], vertices[4], faces[8], faces[2], edges[20], edges[18], edges[8], edges[9] );
        edges[17]->setTopology( vertices[1], vertices[5], faces[2], faces[4], edges[18], edges[19], edges[10], edges[11] );
        edges[18]->setTopology( vertices[4], vertices[5], faces[6], faces[2], edges[19], edges[17], edges[20], edges[16] );
        edges[19]->setTopology( vertices[5], vertices[2], faces[6], faces[4], edges[13], edges[12], edges[18], edges[17] );
        edges[20]->setTopology( vertices[4], vertices[3], faces[8], faces[6], edges[15], edges[14], edges[16], edges[18] );
    }

    faces[0]->setFirstVEdge( edges[0] );
    faces[1]->setFirstVEdge( edges[0] );
    faces[2]->setFirstVEdge( edges[1] );
    faces[3]->setFirstVEdge( edges[2] );
    faces[4]->setFirstVEdge( edges[3] );
    faces[5]->setFirstVEdge( edges[4] );
    faces[6]->setFirstVEdge( edges[5] );
    faces[7]->setFirstVEdge( edges[6] );
    faces[8]->setFirstVEdge( edges[7] );


    // geometry

    vertices[0]->setCircumcircle( rg_Circle2D( minX, maxY, 0.0 ) );
    vertices[1]->setCircumcircle( rg_Circle2D( minX, minY, 0.0 ) );
    vertices[2]->setCircumcircle( rg_Circle2D( maxX, minY, 0.0 ) );
    vertices[3]->setCircumcircle( rg_Circle2D( maxX, maxY, 0.0 ) );

    double radius_circumcircle = 0.0;
    if ( b_horizontal_rectangle ) {
        radius_circumcircle = diffY / 2.0;
        vertices[4]->setCircumcircle( rg_Circle2D( minX + radius_circumcircle, minY + radius_circumcircle, radius_circumcircle ) );
        vertices[5]->setCircumcircle( rg_Circle2D( maxX - radius_circumcircle, minY + radius_circumcircle, radius_circumcircle ) );

    }
    else {
        radius_circumcircle = diffX / 2.0;
        vertices[4]->setCircumcircle( rg_Circle2D( minX + radius_circumcircle, maxY - radius_circumcircle, radius_circumcircle ) );
        vertices[5]->setCircumcircle( rg_Circle2D( minX + radius_circumcircle, minY + radius_circumcircle, radius_circumcircle ) );

    }

    vertices[6]->setCircumcircle( rg_Circle2D( minX, maxY + radius_circumcircle, 0.0 ) );
    vertices[7]->setCircumcircle( rg_Circle2D( minX - radius_circumcircle, maxY, 0.0 ) );
    vertices[8]->setCircumcircle( rg_Circle2D( minX - radius_circumcircle, minY, 0.0 ) );
    vertices[9]->setCircumcircle( rg_Circle2D( minX, minY - radius_circumcircle, 0.0 ) );
    vertices[10]->setCircumcircle( rg_Circle2D( maxX, minY - radius_circumcircle, 0.0 ) );
    vertices[11]->setCircumcircle( rg_Circle2D( maxX + radius_circumcircle, minY, 0.0 ) );
    vertices[12]->setCircumcircle( rg_Circle2D( maxX + radius_circumcircle, maxY, 0.0 ) );
    vertices[13]->setCircumcircle( rg_Circle2D( maxX, maxY + radius_circumcircle, 0.0 ) );
}



void PolygonVD2D_inRectangle::insertOnePolygonIntoVoronoiDiagram( const Polygon2D& polygon )
{
    m_polygons.push_back( polygon );

    list<Generator2D*> newGenerators;
    replicatePolygonEntitiesAsVDGenerators( newGenerators );

    generateChildrenDisksOfPolygonGenerators( newGenerators );

    insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer( newGenerators );

    mergeVFacesOfChildrenDisksOnEachPEdge( newGenerators );

    relocate_some_boundary_Vvertices_to_vertex_generator_centers( newGenerators );

    // Children disk를 'as big as'로 만들었을 때, vertex generator가 anomaly가 될 수 있다.
    // vertex generator의 face는 topologically 삼각형 이상이어야 하기 때문에 다림질 하기 전에 삼각형을 만들어 준다.
    secure_VFace_topology_of_PVertices_by_flipping_the_mating_quill_VEdge_of_relocated_VVertex( newGenerators );


    // polygonVD에 있는 각 VertexGenerator와 EdgeGenerator 마다 childrenVD의 DiskGenerator 들 중 해당되는 Generator로 연결해준다.
    connect_parent_generator_to_first_son_disk_generator( newGenerators );

    identify_internal_and_on_and_external_Vvertices_and_Vedges( newGenerators );

    //refine_location_of_VVertices( newGenerators );
    refine_location_of_VVertices();

    // VFace (즉, VEdge loop)를 다림질한다.
    adjust_topology_of_interior_and_exterior_of_polygon_by_edge_flipping( newGenerators );
}



void PolygonVD2D_inRectangle::constructVoronoiDiagramOfRectContainer( const Polygon2D& inputPolygon )
{
    rg_BoundingBox2D boundingBox;
    inputPolygon.get_bounding_box( boundingBox );

    double height = boundingBox.evaluateYLength();
    double width  = boundingBox.evaluateXLength();
    double offset_length = ( height > width ) ? width : height;
    //offset_length = 0.01;
    offset_length = offset_length * 0.5;

    rg_Point2D offset = rg_Point2D( offset_length, offset_length );
    rg_Point2D minPt = boundingBox.getMinPt() - offset;
    rg_Point2D maxPt = boundingBox.getMaxPt() + offset;


    // create polygon generators for rectangle
    double minX = minPt.getX();
    double minY = minPt.getY();
    double maxX = maxPt.getX();
    double maxY = maxPt.getY();

    //rg_Point2D          boundaryPt[5] = { rg_Point2D( minX, maxY ), rg_Point2D( minX, minY ), rg_Point2D( maxX, minY ), rg_Point2D( maxX, maxY ), rg_Point2D( minX, maxY ) };
    rg_Point2D          boundaryPt[5] = { rg_Point2D( minX, maxY ), rg_Point2D( maxX, maxY ), rg_Point2D( maxX, minY ), rg_Point2D( minX, minY ), rg_Point2D( minX, maxY ) };
    list<rg_Point2D> exteriorBoundaryPts;
    for ( int i = 0; i < 4; ++i ) {
        exteriorBoundaryPts.push_back( boundaryPt[i] );
    }
    m_rectContainer = Polygon2D( exteriorBoundaryPts );

    constructVoronoiDiagramOfRectContainer();



    /*
    VertexGenerator2D*  vertexGenArray[5];
    EdgeGenerator2D*    edgeGenArray[4];

    for ( int i = 0; i < 4; ++i ) {
        VertexGenerator2D* vertexGen = new VertexGenerator2D( boundaryPt[i], i * 2 );
        vertexGen->isThisFromContainer( true );
        m_polygonGenerators.push_back( vertexGen );
        vertexGenArray[i] = vertexGen;

        EdgeGenerator2D* edgeGen = new EdgeGenerator2D( boundaryPt[i], boundaryPt[i + 1], i * 2 + 1 );
        edgeGen->isThisFromContainer( true );
        m_polygonGenerators.push_back( edgeGen );
        edgeGenArray[i] = edgeGen;
    }
    vertexGenArray[4] = vertexGenArray[0];


    for ( int i = 0; i < 4; ++i ) {
        vertexGenArray[i]->set_next_edge_generator( edgeGenArray[i] );
        vertexGenArray[i + 1]->set_previous_edge_generator( edgeGenArray[i] );

        edgeGenArray[i]->set_start_vertex_generator( vertexGenArray[i] );
        edgeGenArray[i]->set_end_vertex_generator( vertexGenArray[i + 1] );
    }

    // create entities
    //  1. make elements of voronoi diagram of phantom generators

    VVertex2D* vertices[14] = {
        createVertex( 0 ),    // LT
        createVertex( 1 ),    // LB
        createVertex( 2 ),    // RB
        createVertex( 3 ),    // RT
        createVertex( 4 ),    // IN_1
        createVertex( 5 ),    // IN_2
        createVertex( 6 ),    // LTT
        createVertex( 7 ),    // LTL
        createVertex( 8 ),    // LBL
        createVertex( 9 ),    // LBB
        createVertex( 10 ),   // RBB
        createVertex( 11 ),   // RBR
        createVertex( 12 ),   // RTR
        createVertex( 13 ) };  // RTT

    VEdge2D* edges[21] = {
        createEdge( 0 ),        // INF_LT
        createEdge( 1 ),        // INF_L
        createEdge( 2 ),        // INF_LB
        createEdge( 3 ),        // INF_B
        createEdge( 4 ),        // INF_RB
        createEdge( 5 ),        // INF_R
        createEdge( 6 ),        // INF_RT
        createEdge( 7 ),        // INF_T
        createEdge( 8 ),        // UBND_TL
        createEdge( 9 ),        // UBND_LT
        createEdge( 10 ),       // UBND_LB
        createEdge( 11 ),       // UBND_BL
        createEdge( 12 ),       // UBND_BR
        createEdge( 13 ),       // UBND_RB
        createEdge( 14 ),       // UBND_RT
        createEdge( 15 ),       // UBND_TR
        createEdge( 16 ),       // IN_LT
        createEdge( 17 ),       // IN_LB
        createEdge( 18 ),       // IN_CENTER
        createEdge( 19 ),       // IN_RB
        createEdge( 20 ) };      // IN_RT

    VFace2D* faces[9] = {
        createFace( -1 ),       // INF
        createFace( 0 ),        // VTX_LT
        createFace( 1 ),        // EDGE_L
        createFace( 2 ),        // VTX_LB
        createFace( 3 ),        // EDGE_B
        createFace( 4 ),        // VTX_RB
        createFace( 5 ),        // EDGE_R
        createFace( 6 ),        // VTX_RT
        createFace( 7 ) };       // EDGE_T


    // IN_OUT marking.....
    m_VVertexToLocationStatus[vertices[0]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[1]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[2]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[3]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[4]]  = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[5]]  = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[6]]  = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[7]]  = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[8]]  = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[9]]  = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[10]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[11]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[12]] = INFINITE_ENTITY;
    m_VVertexToLocationStatus[vertices[13]] = INFINITE_ENTITY;

    m_VEdgeToLocationStatus[edges[0]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[1]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[2]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[3]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[4]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[5]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[6]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[7]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[8]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[9]]       = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[10]]      = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[11]]      = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[12]]      = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[13]]      = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[14]]      = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[15]]      = INFINITE_ENTITY;
    m_VEdgeToLocationStatus[edges[16]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[17]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[18]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[19]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[20]]      = OUTSIDE_POLYGON;

    // topology
    //vertexGenArray[0]->setOuterFace( faces[1] );
    //vertexGenArray[1]->setOuterFace( faces[3] );
    //vertexGenArray[2]->setOuterFace( faces[5] );
    //vertexGenArray[3]->setOuterFace( faces[7] );
    vertexGenArray[0]->setOuterFace( faces[1] );
    vertexGenArray[1]->setOuterFace( faces[7] );
    vertexGenArray[2]->setOuterFace( faces[5] );
    vertexGenArray[3]->setOuterFace( faces[3] );

    vertexGenArray[0]->setUserData( vertexGenArray[0] );
    vertexGenArray[1]->setUserData( vertexGenArray[1] );
    vertexGenArray[2]->setUserData( vertexGenArray[2] );
    vertexGenArray[3]->setUserData( vertexGenArray[3] );

    //edgeGenArray[0]->setOuterFace( faces[2] );
    //edgeGenArray[1]->setOuterFace( faces[4] );
    //edgeGenArray[2]->setOuterFace( faces[6] );
    //edgeGenArray[3]->setOuterFace( faces[8] );
    edgeGenArray[0]->setOuterFace( faces[8] );
    edgeGenArray[1]->setOuterFace( faces[6] );
    edgeGenArray[2]->setOuterFace( faces[4] );
    edgeGenArray[3]->setOuterFace( faces[2] );

    edgeGenArray[0]->setUserData( edgeGenArray[0] );
    edgeGenArray[1]->setUserData( edgeGenArray[1] );
    edgeGenArray[2]->setUserData( edgeGenArray[2] );
    edgeGenArray[3]->setUserData( edgeGenArray[3] );

    //faces[1]->setGenerator( (Generator2D*)vertexGenArray[0] );
    //faces[3]->setGenerator( (Generator2D*)vertexGenArray[1] );
    //faces[5]->setGenerator( (Generator2D*)vertexGenArray[2] );
    //faces[7]->setGenerator( (Generator2D*)vertexGenArray[3] );
    //faces[2]->setGenerator( (Generator2D*)edgeGenArray[0] );
    //faces[4]->setGenerator( (Generator2D*)edgeGenArray[1] );
    //faces[6]->setGenerator( (Generator2D*)edgeGenArray[2] );
    //faces[8]->setGenerator( (Generator2D*)edgeGenArray[3] );
    faces[1]->setGenerator( (Generator2D*)vertexGenArray[0] );
    faces[3]->setGenerator( (Generator2D*)vertexGenArray[3] );
    faces[5]->setGenerator( (Generator2D*)vertexGenArray[2] );
    faces[7]->setGenerator( (Generator2D*)vertexGenArray[1] );
    faces[2]->setGenerator( (Generator2D*)edgeGenArray[3] );
    faces[4]->setGenerator( (Generator2D*)edgeGenArray[2] );
    faces[6]->setGenerator( (Generator2D*)edgeGenArray[1] );
    faces[8]->setGenerator( (Generator2D*)edgeGenArray[0] );


    vertices[0]->setFirstVEdge( edges[16] );
    vertices[1]->setFirstVEdge( edges[17] );
    vertices[2]->setFirstVEdge( edges[19] );
    vertices[3]->setFirstVEdge( edges[20] );
    vertices[4]->setFirstVEdge( edges[18] );
    vertices[5]->setFirstVEdge( edges[18] );
    vertices[6]->setFirstVEdge( edges[0] );
    vertices[7]->setFirstVEdge( edges[1] );
    vertices[8]->setFirstVEdge( edges[2] );
    vertices[9]->setFirstVEdge( edges[3] );
    vertices[10]->setFirstVEdge( edges[4] );
    vertices[11]->setFirstVEdge( edges[5] );
    vertices[12]->setFirstVEdge( edges[6] );
    vertices[13]->setFirstVEdge( edges[7] );

    edges[0]->setTopology( vertices[6], vertices[7], faces[1], faces[0], edges[9], edges[1], edges[8], edges[7] );
    edges[1]->setTopology( vertices[7], vertices[8], faces[2], faces[0], edges[10], edges[2], edges[9], edges[0] );
    edges[2]->setTopology( vertices[8], vertices[9], faces[3], faces[0], edges[11], edges[3], edges[10], edges[1] );
    edges[3]->setTopology( vertices[9], vertices[10], faces[4], faces[0], edges[12], edges[4], edges[11], edges[2] );
    edges[4]->setTopology( vertices[10], vertices[11], faces[5], faces[0], edges[13], edges[5], edges[12], edges[3] );
    edges[5]->setTopology( vertices[11], vertices[12], faces[6], faces[0], edges[14], edges[6], edges[13], edges[4] );
    edges[6]->setTopology( vertices[12], vertices[13], faces[7], faces[0], edges[15], edges[7], edges[14], edges[5] );
    edges[7]->setTopology( vertices[13], vertices[6], faces[8], faces[0], edges[8], edges[0], edges[15], edges[6] );

    edges[8]->setTopology( vertices[0], vertices[6], faces[1], faces[8], edges[0], edges[7], edges[9], edges[16] );
    edges[9]->setTopology( vertices[0], vertices[7], faces[2], faces[1], edges[1], edges[0], edges[16], edges[8] );
    edges[10]->setTopology( vertices[1], vertices[8], faces[3], faces[2], edges[2], edges[1], edges[11], edges[17] );
    edges[11]->setTopology( vertices[1], vertices[9], faces[4], faces[3], edges[3], edges[2], edges[17], edges[10] );
    edges[12]->setTopology( vertices[2], vertices[10], faces[5], faces[4], edges[4], edges[3], edges[13], edges[19] );
    edges[13]->setTopology( vertices[2], vertices[11], faces[6], faces[5], edges[5], edges[4], edges[19], edges[12] );
    edges[14]->setTopology( vertices[3], vertices[12], faces[7], faces[6], edges[6], edges[5], edges[15], edges[20] );
    edges[15]->setTopology( vertices[3], vertices[13], faces[8], faces[7], edges[7], edges[6], edges[20], edges[14] );

    bool b_horizontal_rectangle = true;
    double diffX = maxX - minX;
    double diffY = maxY - minY;
    if ( diffX < diffY )
        b_horizontal_rectangle = false;

    if ( b_horizontal_rectangle ) {
        edges[16]->setTopology( vertices[0], vertices[4], faces[8], faces[2], edges[18], edges[17], edges[8], edges[9] );
        edges[17]->setTopology( vertices[1], vertices[4], faces[2], faces[4], edges[16], edges[18], edges[10], edges[11] );
        edges[18]->setTopology( vertices[4], vertices[5], faces[8], faces[4], edges[20], edges[19], edges[16], edges[17] );
        edges[19]->setTopology( vertices[5], vertices[2], faces[6], faces[4], edges[13], edges[12], edges[20], edges[18] );
        edges[20]->setTopology( vertices[5], vertices[3], faces[8], faces[6], edges[15], edges[14], edges[18], edges[19] );
    }
    else {
        edges[16]->setTopology( vertices[0], vertices[4], faces[8], faces[2], edges[20], edges[18], edges[8], edges[9] );
        edges[17]->setTopology( vertices[1], vertices[5], faces[2], faces[4], edges[18], edges[19], edges[10], edges[11] );
        edges[18]->setTopology( vertices[4], vertices[5], faces[6], faces[2], edges[19], edges[17], edges[20], edges[16] );
        edges[19]->setTopology( vertices[5], vertices[2], faces[6], faces[4], edges[13], edges[12], edges[18], edges[17] );
        edges[20]->setTopology( vertices[4], vertices[3], faces[8], faces[6], edges[15], edges[14], edges[16], edges[18] );
    }

    faces[0]->setFirstVEdge( edges[0] );
    faces[1]->setFirstVEdge( edges[0] );
    faces[2]->setFirstVEdge( edges[1] );
    faces[3]->setFirstVEdge( edges[2] );
    faces[4]->setFirstVEdge( edges[3] );
    faces[5]->setFirstVEdge( edges[4] );
    faces[6]->setFirstVEdge( edges[5] );
    faces[7]->setFirstVEdge( edges[6] );
    faces[8]->setFirstVEdge( edges[7] );


    // geometry

    vertices[0]->setCircumcircle( rg_Circle2D( minX, maxY, 0.0 ) );
    vertices[1]->setCircumcircle( rg_Circle2D( minX, minY, 0.0 ) );
    vertices[2]->setCircumcircle( rg_Circle2D( maxX, minY, 0.0 ) );
    vertices[3]->setCircumcircle( rg_Circle2D( maxX, maxY, 0.0 ) );

    double radius_circumcircle = 0.0;
    if ( b_horizontal_rectangle ) {
        radius_circumcircle = diffY / 2.0;
        vertices[4]->setCircumcircle( rg_Circle2D( minX + radius_circumcircle, minY + radius_circumcircle, radius_circumcircle ) );
        vertices[5]->setCircumcircle( rg_Circle2D( maxX - radius_circumcircle, minY + radius_circumcircle, radius_circumcircle ) );

    }
    else {
        radius_circumcircle = diffX / 2.0;
        vertices[4]->setCircumcircle( rg_Circle2D( minX + radius_circumcircle, maxY - radius_circumcircle, radius_circumcircle ) );
        vertices[5]->setCircumcircle( rg_Circle2D( minX + radius_circumcircle, minY + radius_circumcircle, radius_circumcircle ) );

    }

    vertices[6]->setCircumcircle( rg_Circle2D( minX, maxY + radius_circumcircle, 0.0 ) );
    vertices[7]->setCircumcircle( rg_Circle2D( minX - radius_circumcircle, maxY, 0.0 ) );
    vertices[8]->setCircumcircle( rg_Circle2D( minX - radius_circumcircle, minY, 0.0 ) );
    vertices[9]->setCircumcircle( rg_Circle2D( minX, minY - radius_circumcircle, 0.0 ) );
    vertices[10]->setCircumcircle( rg_Circle2D( maxX, minY - radius_circumcircle, 0.0 ) );
    vertices[11]->setCircumcircle( rg_Circle2D( maxX + radius_circumcircle, minY, 0.0 ) );
    vertices[12]->setCircumcircle( rg_Circle2D( maxX + radius_circumcircle, maxY, 0.0 ) );
    vertices[13]->setCircumcircle( rg_Circle2D( maxX, maxY + radius_circumcircle, 0.0 ) );
    */
}



void PolygonVD2D_inRectangle::insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer()
{
    list< rg_Triplet<rg_Circle2D, int, void*> > disk_UserID_UserData_TripletList;
    int childDiskID = m_polygonGenerators.back()->getID();


    for ( list<Generator2D*>::iterator i_gen = m_polygonGenerators.begin(); i_gen != m_polygonGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        Generator2D::Generator_Type type = currGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

            if ( vertexGen->isThisFromContainer() ) {
                break;
            }

            vertexGen->setID( childDiskID + 1 );
            rg_Circle2D childDisk = vertexGen->get_child_disk();
            disk_UserID_UserData_TripletList.push_back( rg_Triplet<rg_Circle2D, int, void*>( childDisk, ++childDiskID, vertexGen ) );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)currGen;

            if ( edgeGen->isThisFromContainer() ) {
                break;
            }

            edgeGen->setID( childDiskID + 1 );
            list<rg_Circle2D>& childrenDisks = edgeGen->get_children_disks();
            for ( list<rg_Circle2D>::iterator i_childDisk = childrenDisks.begin(); i_childDisk != childrenDisks.end(); ++i_childDisk ) {
                disk_UserID_UserData_TripletList.push_back( rg_Triplet<rg_Circle2D, int, void*>( *i_childDisk, ++childDiskID, edgeGen ) );
            }
        }
        break;

        default:
            break;
        }
    }

    int iii = 0;

    list< pair< DiskGenerator2D*, void* > > childrenGens_to_parentGen;
    for ( list< rg_Triplet<rg_Circle2D, int, void*> >::iterator i_diskTriplet = disk_UserID_UserData_TripletList.begin(); i_diskTriplet != disk_UserID_UserData_TripletList.end(); ++i_diskTriplet ) {
        rg_Triplet<rg_Circle2D, int, void*> currTriplet = *i_diskTriplet;
        DiskGenerator2D* currGen = insertOneChildDiskInPolygonVoronoiDiagram( currTriplet );
        childrenGens_to_parentGen.push_back( make_pair( currGen, currTriplet.m_third ) );

        ++iii;
        if ( iii == 20 ) {
            int stopHere1 = 1;
            //break;
        }
    }

    for ( list< pair< DiskGenerator2D*, void* > >::iterator i_genPair = childrenGens_to_parentGen.begin(); i_genPair != childrenGens_to_parentGen.end(); ++i_genPair ) {
        i_genPair->first->setUserData( i_genPair->second );
    }
}



DiskGenerator2D* PolygonVD2D_inRectangle::insertOneChildDiskInPolygonVoronoiDiagram( const rg_Triplet<rg_Circle2D, int, void*>& diskIDUserData )
{
    /*
    각 Polygon Vertex 마다 하나의 VVertex (on-vertex)가 놓인다.
    on-vertex는 findSeedRedVertex, wavePropagation_ver1 함수에서 절대로 red vertex가 되어서는 안된다.
    하지만, 수치적인 문제로 violate 될 수 있다.
    (현재의 코드는 red vertex가 되지 않도록 강제해놓았기 때문에 문제가 없음).
    (Rectangle에서는 이런일이 없으나, reflex vertex에서 생길 수 있는 현상임).
    */
    void* userData = diskIDUserData.m_third;
    DiskGenerator2D* newGenerator = new DiskGenerator2D( diskIDUserData.m_second, diskIDUserData.m_first, diskIDUserData.m_third );
    newGenerator->setUserData( newGenerator );
    m_generators.push_back( newGenerator );


    // This is same in the incremental process of the TOI_D2 algorithm.
    // From here (1)
    VFace2D* anchorCell = findAnchorCell( newGenerator->getDisk().getCenterPt() );
    VVertex2D* seedRedVertex = findSeedRedVertex( newGenerator, (Generator2D*)anchorCell->getGenerator() );

    list<VVertex2D*> redVVertices;
    list<VVertex2D*> blueVVertices;
    list<VVertex2D*> fictitiousVVertices;

    seedRedVertex->setStatus( RED_V );
    redVVertices.push_back( seedRedVertex );

    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );

    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeNewCellAndConnectToCurrVD( newGenerator, redVVertices, blueVVertices, newVVertices, newVEdges );
    // To here (1)


    // This is same in the incremental process of the TOI_D2 algorithm.
    // From here (2)
    connectCurrVDToNewVEdgeLoop( newVEdges );
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
    // To here (2)

    // This is unique in the VD of simple polygon
    // Geometry computation: 그 외의 모든 것은 topology operation임.
    // Edge geometry는 계산하지 않음.
    computeCoordOfNewVVertices_before_merging( newVVertices );


    // disk packing에서 사용되는 map
    // Voronoi diagram에서는 disk의 precision을 single로 유지한다.
    // disk packing문제에서는 double precision의 disk를 가지고 intersection 여부를 따져야 하기 때문에
    // VD 바깥쪽에 있는 double precision의 disk와 VD 내부의 disk를 mapping하고 있다.
    m_map_DiskID2VFace.insert( make_pair( newGenerator->getID(), newGenerator->getOuterFace() ) );

    //newGenerator->setUserData( userData );

    return newGenerator;
}


/*
void PolygonVD2D_inRectangle::identify_internal_and_on_and_external_Vvertices_and_Vedges()
{
    list< pair< VertexGenerator2D*, VVertex2D* > > PVertex_N_onVVertexPair;
    identify_Vvertices_on_boundary( PVertex_N_onVVertexPair );

    stack<VEdge2D*> outerEdges;
    stack<VEdge2D*> innerEdges;

    // collect initial stacks of outer edges and inner edges.
    for ( list< pair< VertexGenerator2D*, VVertex2D* > >::iterator i_pair = PVertex_N_onVVertexPair.begin(); i_pair != PVertex_N_onVVertexPair.end(); ++i_pair ) {
        VertexGenerator2D*  currPVertex     = i_pair->first;
        VVertex2D*          currOnVVertex   = i_pair->second;

        EdgeGenerator2D*    prevPEdge       = currPVertex->get_previous_edge_generator();
        EdgeGenerator2D*    nextPEdge       = currPVertex->get_next_edge_generator();

        VertexGenerator2D::Vertex_Type vertexType = currPVertex->vertex_type();

        // find THAT edge
        list<VEdge2D*> incidentEdges;
        currOnVVertex->getIncident3VEdges( incidentEdges );
        for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
            VEdge2D* currEdge = *i_edge;

            bool b_this_edge_is_defined_by_two_PEdges = true;
            if ( currEdge->getLeftFace()->getGenerator()->getUserData() == currPVertex || currEdge->getRightFace()->getGenerator()->getUserData() == currPVertex ) {
                b_this_edge_is_defined_by_two_PEdges = false;
            }
            else {
                b_this_edge_is_defined_by_two_PEdges = true;
            }


            if ( vertexType == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
                if ( b_this_edge_is_defined_by_two_PEdges ) {
                    outerEdges.push( currEdge );
                }
                else {
                    innerEdges.push( currEdge );
                }
            }
            else {
                if ( b_this_edge_is_defined_by_two_PEdges ) {
                    innerEdges.push( currEdge );
                }
                else {
                    outerEdges.push( currEdge );
                }
            }
        }
    }

    // outer edges
    while ( !outerEdges.empty() )
    {
        VEdge2D* outerEdge = outerEdges.top();
        outerEdges.pop();

        VVertex2D* startVertex  = outerEdge->getStartVertex();
        VVertex2D* endVertex    = outerEdge->getEndVertex();
        VVertex_Of_FirstSonVD_To_Location_Status_Map::iterator i_startVertexToLocationStatus    = m_VVertexToLocationStatus.find( startVertex );
        VVertex_Of_FirstSonVD_To_Location_Status_Map::iterator i_endVertexToLocationStatus      = m_VVertexToLocationStatus.find( endVertex );

        // If start vertex was visited before
        if ( i_startVertexToLocationStatus != m_VVertexToLocationStatus.end() )
        {
            if ( i_endVertexToLocationStatus != m_VVertexToLocationStatus.end() )
            {
                m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );
                continue;
            }
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( endVertex, OUTSIDE_POLYGON ) );
            m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            endVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }

        // If end vertex was visited before
        if ( i_endVertexToLocationStatus != m_VVertexToLocationStatus.end() )
        {
            if ( i_startVertexToLocationStatus != m_VVertexToLocationStatus.end() )
            {
                m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );
                continue;
            }
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( startVertex, OUTSIDE_POLYGON ) );
            m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            startVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }
    }

    set_location_status_of_remaining_Vvertices_and_Vedges_with_inside();
}
*/


void PolygonVD2D_inRectangle::identify_Vvertices_on_boundary( list< pair< VertexGenerator2D*, VVertex2D* > >& PVertex_N_onVVertexPair )
{
    for ( list<Generator2D*>::iterator i_gen = m_polygonGenerators.begin(); i_gen != m_polygonGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() != Generator2D::Generator_Type::VERTEX_G )
            continue;

        VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

        if ( vertexGen->isThisFromContainer() ) {
            continue;
        }

        EdgeGenerator2D* prevEdgeGen = vertexGen->get_previous_edge_generator();
        EdgeGenerator2D* nextEdgeGen = vertexGen->get_next_edge_generator();

        VEdge2D* VEdge_between_prevEdgeGen_N_vertexGen = rg_NULL;
        VEdge2D* VEdge_between_vertexGen_N_nextEdgeGen = rg_NULL;
        find_two_VEdges_between_vertex_generator_N_incident_edge_generators( vertexGen, prevEdgeGen, nextEdgeGen, VEdge_between_prevEdgeGen_N_vertexGen, VEdge_between_vertexGen_N_nextEdgeGen );

        VVertex2D* startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen = VEdge_between_prevEdgeGen_N_vertexGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen = VEdge_between_prevEdgeGen_N_vertexGen->getEndVertex();
        VVertex2D* startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getEndVertex();

        VVertex2D* commonVVertex = rg_NULL;
        if ( startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen || startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
        {
            commonVVertex = startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
        }
        else if ( endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen || endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
        {
            commonVVertex = endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
        }
        else
        {
        }

        if ( commonVVertex != rg_NULL )
        {
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( commonVVertex, ON_POLYGON_BOUNDARY ) );
            PVertex_N_onVVertexPair.push_back( make_pair( vertexGen, commonVVertex ) );
        }
    }
}



void PolygonVD2D_inRectangle::adjust_topology_of_interior_and_exterior_of_polygon_by_edge_flipping()
{
    for ( list<Generator2D*>::iterator i_generators = m_polygonGenerators.begin(); i_generators != m_polygonGenerators.end(); ++i_generators ) {
        Generator2D* currGenerator = *i_generators;
        Generator2D::Generator_Type type = currGenerator->getType();

        if ( type == Generator2D::Generator_Type::VERTEX_G )
        {
            VertexGenerator2D* currVertexGen = (VertexGenerator2D*)currGenerator;

            if ( currVertexGen->isThisFromContainer() ) {
                continue;
            }

            // This function can handle reflex and non-reflex vertex.
            flatten_topology_of_VFace_boundary_of_reflex_vertex_gen_with_possible_flipping_out_of_VEdges( currVertexGen );
        }
        else if ( type == Generator2D::Generator_Type::EDGE_G ) {
            EdgeGenerator2D* currEdgeGen = (EdgeGenerator2D*)currGenerator;

            if ( currEdgeGen->isThisFromContainer() ) {
                continue;
            }

            flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges( currEdgeGen );
            flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges_exterior( currEdgeGen );
        }
    }
}

void PolygonVD2D_inRectangle::computeCoordOfNewVVertices_before_merging( list<VVertex2D*>& newVVertices )
{
    // This funtion is for computing the coordinate of new V-vertices 
    // while we insert children disks of PVertices and PEdges of input polygon 
    // into the Voronoi diagram of the rectangular container.

    // The generators which define a new VVertex
    // 1. include at lease one disk
    // 2. do not include any PVertex
    // so there are only DDD, DDL, DLL cases.
    // (D: Disk, L: Line, P: Point)

    // All PEdge is of the rectangular container and its direction is reversed 
    // (because the status of Voronoi entities which are inside of the container are OUTSIDE_POLYGON).
#ifdef DEBUG_VERTEX
    ofstream fout_debug_vertex( "vertex_debugging.txt" );
#endif

    for ( list<VVertex2D*>::iterator i_vtx = newVVertices.begin(); i_vtx != newVVertices.end(); ++i_vtx ) {
        VVertex2D* vertex = *i_vtx;

#ifdef DEBUG_VERTEX
        fout_debug_vertex << vertex->getID() << endl;
        rg_Point2D foorprint_debug[2];
        rg_Circle2D insertedDisk;
#endif

#ifdef CHECK_COMP_TIME
        clock_t startTime, endTime;
        startTime = clock();
#endif

        int numDisks, numLines, numPoints;
        Generator2D* parentGens[3];

        numDisks = numLines = numPoints = 0;
        list<Generator2D*> diskGens;
        vertex->getDefining3Generators( diskGens );

        int index_gen = 0;
        for ( list<Generator2D*>::iterator i_gen = diskGens.begin(); i_gen != diskGens.end(); ++i_gen, ++index_gen ) {
            Generator2D* parentGen = (Generator2D*)( *i_gen )->getUserData();
            Generator2D::Generator_Type type = parentGen->getType();

            switch ( type )
            {
            case Generator2D::Generator_Type::DISK_G:
            {
                ++numDisks;
                parentGens[index_gen] = parentGen;
            }
            break;

            case Generator2D::Generator_Type::EDGE_G:
            {
                if ( ( (EdgeGenerator2D*)parentGen )->isThisFromContainer() ) {
                    ++numLines;
                    parentGens[index_gen] = parentGen;
                }
                else {
                    ++numDisks;
                    parentGens[index_gen] = *i_gen;
                }
            }
            break;

            case Generator2D::Generator_Type::VERTEX_G:
            {
                if ( ( (VertexGenerator2D*)parentGen )->isThisFromContainer() ) {
                    ++numPoints;
                    parentGens[index_gen] = parentGen;
                }
                else {
                    ++numDisks;
                    parentGens[index_gen] = *i_gen;
                }
            }
            break;

            default:
                break;
            }
        }

#ifdef CHECK_COMP_TIME
        endTime = clock();
        t_check_parents = t_check_parents + endTime - startTime;
#endif

        //rg_Circle2D refinedTangentCircle;
        rg_Circle2D circumcircle_of_vertex;
        rg_Circle2D tangentCircle[2];

        switch ( numDisks )
        {
        case 3: // numDisks == 3, numLines == 0, numPoints == 0
        {
#ifdef CHECK_COMP_TIME
            ++num_DDD;
            startTime = clock();
#endif
            rg_Circle2D disks[3];
            disks[0] = ( (DiskGenerator2D*)parentGens[0] )->getDisk();
            disks[1] = ( (DiskGenerator2D*)parentGens[1] )->getDisk();
            disks[2] = ( (DiskGenerator2D*)parentGens[2] )->getDisk();

            int numTangentCircles = rg_Circle2D::makeCircumcircle( disks[0], disks[1], disks[2], tangentCircle[0], tangentCircle[1] );
#ifdef CHECK_COMP_TIME
            endTime = clock();
            t_comp_V_DDD = t_comp_V_DDD + endTime - startTime;

            startTime = clock();
#endif
            if ( numTangentCircles == 1 ) {
                circumcircle_of_vertex = tangentCircle[0];
            }
            else {
                if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                    rg_Circle2D tempCircle = tangentCircle[0];
                    tangentCircle[0] = tangentCircle[1];
                    tangentCircle[1] = tempCircle;
                }

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
            }
#ifdef CHECK_COMP_TIME
            endTime = clock();
            t_check_orientation = t_check_orientation + endTime - startTime;
#endif
        }
        break;



        case 2:
        {
            // numDisks == 2, numLines == 1, numPoints == 0
            if ( numLines == 1 )
            {
#ifdef CHECK_COMP_TIME
                ++num_DDL;
                startTime = clock();
#endif
                rg_Circle2D disks[2];
                rg_Line2D lineSeg3;
                rg_INT diskIndex = 0;
                for ( rg_INT i = 0; i < 3; i++ )
                {
                    if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                        disks[diskIndex++] = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                    if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                        ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg3 );
                }

                int numTangentCircles = 0;
                rg_Circle2D tangentCircle[2];
                lineSeg3 = lineSeg3.get_reversed_line2D();
                numTangentCircles = computeCoordOfNewVVertex_of_two_disks_and_a_line( disks[0], disks[1], lineSeg3, tangentCircle[0], tangentCircle[1] );

#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_comp_V_DDL = t_comp_V_DDL + endTime - startTime;

                startTime = clock();
#endif
                if ( numTangentCircles == 1 ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                        rg_Circle2D tempCircle = tangentCircle[0];
                        tangentCircle[0] = tangentCircle[1];
                        tangentCircle[1] = tempCircle;
                    }

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        circumcircle_of_vertex = tangentCircle[1];
                    }
                }
#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_check_orientation = t_check_orientation + endTime - startTime;
#endif
            }

            // numDisks == 2, numLines == 0, numPoints == 1
            else
            {
#ifdef CHECK_COMP_TIME
                ++num_DDD;
                startTime = clock();
#endif
                rg_Circle2D disks[2];
                rg_Point2D point3;
                rg_INT diskIndex = 0;
                for ( rg_INT i = 0; i < 3; i++ )
                {
                    if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                        disks[diskIndex++] = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                    if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                        ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point3 );
                }

                int numTangentCircles = rg_Circle2D::makeCircumcircle( disks[0], disks[1], rg_Circle2D( point3, 0.0 ), tangentCircle[0], tangentCircle[1] );
#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_comp_V_DDD = t_comp_V_DDD + endTime - startTime;

                startTime = clock();
#endif
                if ( numTangentCircles == 1 ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                        rg_Circle2D tempCircle = tangentCircle[0];
                        tangentCircle[0] = tangentCircle[1];
                        tangentCircle[1] = tempCircle;
                    }

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        circumcircle_of_vertex = tangentCircle[1];
                    }
                }
#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_check_orientation = t_check_orientation + endTime - startTime;
#endif
            }
        }
        break;



        case 1:
        {
            switch ( numLines )
            {
            case 2: // numDisks == 1, numLines == 2, numPoints == 0
            {
#ifdef CHECK_COMP_TIME
                ++num_DLL;
                startTime = clock();
#endif
                rg_Circle2D disk1;
                rg_Line2D lineSeg[2];
                rg_INT lineSegIndex = 0;
                for ( rg_INT i = 0; i < 3; i++ )
                {
                    if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                        disk1 = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                    if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                        ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg[lineSegIndex++] );
                }

                lineSeg[0] = lineSeg[0].get_reversed_line2D();
                lineSeg[1] = lineSeg[1].get_reversed_line2D();

                VEdge2D* edge_of_two_polygon_edges = NULL;
                list<VEdge2D*> incidentEdges;
                vertex->getIncident3VEdges( incidentEdges );
                for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
                    VEdge2D* currEdge = *i_edge;

                    Generator2D* leftParent = (Generator2D*)currEdge->getLeftFace()->getGenerator()->getUserData();
                    Generator2D* rightParent = (Generator2D*)currEdge->getRightFace()->getGenerator()->getUserData();

                    if ( leftParent->getType() == Generator2D::Generator_Type::EDGE_G && rightParent->getType() == Generator2D::Generator_Type::EDGE_G ) {
                        edge_of_two_polygon_edges = currEdge;
                        break;
                    }
                }

                VVertex2D* vertex_which_has_coordinate = NULL;
                if ( edge_of_two_polygon_edges->getStartVertex() == vertex ) {
                    vertex_which_has_coordinate = edge_of_two_polygon_edges->getEndVertex();
                }
                else {
                    vertex_which_has_coordinate = edge_of_two_polygon_edges->getStartVertex();
                }

                rg_Point2D dirVec_line1 = lineSeg[0].getEP() - lineSeg[0].getSP();
                rg_Point2D dirVec_line2 = lineSeg[1].getSP() - lineSeg[1].getEP();
                dirVec_line1 = dirVec_line1.getUnitVector();
                dirVec_line2 = dirVec_line2.getUnitVector();

                rg_Point2D dirVec_edge;
                double cross_product = dirVec_line1 * dirVec_line2;
                if ( rg_ZERO( cross_product ) ) {
                    dirVec_edge = dirVec_line1;
                }
                //else {
                //    dirvec_edge = dirvec_line1 + dirvec_line2;
                //    dirvec_edge = dirvec_edge.getunitvector();
                //}
                else if ( rg_POS( cross_product ) ) {
                    dirVec_edge = dirVec_line1 + dirVec_line2;
                    dirVec_edge = dirVec_edge.getUnitVector();
                }
                else { // if ( rg_NE( cross_product ) ) {
                    dirVec_edge = -dirVec_line1 - dirVec_line2;
                    dirVec_edge = dirVec_edge.getUnitVector();
                }

                rg_Point2D SP = vertex_which_has_coordinate->getLocation();
                rg_Point2D EP = SP + dirVec_edge;
                rg_Line2D geometry_edge( SP, EP );

                rg_Point2D footprint_SP, footprint_EP;
                lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( SP, footprint_SP );
                lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( EP, footprint_EP );

                double radius_SP = SP.distance( footprint_SP );
                double radius_EP = EP.distance( footprint_EP );





                computeCoordOfNewVVertex_on_line_of_two_polygon_edges( geometry_edge, rg_Circle2D( SP, radius_SP ), rg_Circle2D( EP, radius_EP ), disk1, tangentCircle[0], tangentCircle[1] );
#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_comp_V_DLL = t_comp_V_DLL + endTime - startTime;

                startTime = clock();
#endif
                if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                    rg_Circle2D tempCircle = tangentCircle[0];
                    tangentCircle[0] = tangentCircle[1];
                    tangentCircle[1] = tempCircle;
                }

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_check_orientation = t_check_orientation + endTime - startTime;
#endif

#ifdef DEBUG_VERTEX
                lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( circumcircle_of_vertex.getCenterPt(), foorprint_debug[0] );
                lineSeg[1].compute_perpendicular_footprint_of_point_onto_entire_line( circumcircle_of_vertex.getCenterPt(), foorprint_debug[1] );
                insertedDisk = disk1;
#endif
            }
            break;

            case 1: // numDisks == 1, numLines == 1, numPoints == 1
            {
#ifdef CHECK_COMP_TIME
                startTime = clock();
#endif
                EdgeGenerator2D* edgeGen = rg_NULL;
                VertexGenerator2D* vertexGen = rg_NULL;

                rg_Circle2D disk1;
                rg_Line2D  lineSeg2;
                rg_Point2D point3;
                for ( rg_INT i = 0; i < 3; i++ )
                {
                    if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() ) {
                        disk1 = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                    }
                    if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                    {
                        edgeGen = (EdgeGenerator2D*)parentGens[i];
                        ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg2 );
                    }
                    if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                    {
                        vertexGen = (VertexGenerator2D*)parentGens[i];
                        ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point3 );
                    }
                }
#ifdef CHECK_COMP_TIME
                endTime = clock();
                double localCompTime = endTime - startTime;
#endif

                lineSeg2 = lineSeg2.get_reversed_line2D();
                if ( vertexGen->get_next_edge_generator() != edgeGen && vertexGen->get_previous_edge_generator() != edgeGen ) {
#ifdef CHECK_COMP_TIME
                    num_DDL++;
                    startTime = clock();
#endif
                    int numTangentCircles = 0;
                    rg_Circle2D tangentCircle[2];
                    numTangentCircles = computeCoordOfNewVVertex_of_two_disks_and_a_line( rg_Circle2D( point3, 0.0 ), disk1, lineSeg2, tangentCircle[0], tangentCircle[1] );
#ifdef CHECK_COMP_TIME
                    endTime = clock();
                    t_comp_V_DDL = t_comp_V_DDL + endTime - startTime + localCompTime;

                    startTime = clock();
#endif
                    if ( numTangentCircles == 1 ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                            rg_Circle2D tempCircle = tangentCircle[0];
                            tangentCircle[0] = tangentCircle[1];
                            tangentCircle[1] = tempCircle;
                        }

                        if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                            circumcircle_of_vertex = tangentCircle[0];
                        }
                        else {
                            circumcircle_of_vertex = tangentCircle[1];
                        }
                    }
#ifdef CHECK_COMP_TIME
                    endTime = clock();
                    t_check_orientation = t_check_orientation + endTime - startTime;
#endif

                    ////////////////////// CURVE - CURVE INTERSECTION ////////////////////////////
                    /*
                    rg_RQBzCurve2D curve1 = createRQBzCurveOfParabola(rg_Circle2D(point3, 0.0), lineSeg2);
                    rg_RQBzCurve2D curve2 = createRQBzCurveOfParabola(disk1, lineSeg2);

                    rg_dList <rg_Point2D> intersectList;
                    rg_IntersectFunc::intersectRQBzCurveVsRQBzCurve(curve1, curve2, intersectList);

                    rg_Point2D circumcircle_center;

                    if ( intersectList.getSize() == 0 ) {
                        cout << "CANNOT COMPUTE COORDINATE OF VERTEX " << vertex->getID() << endl;
                    }


                    if ( intersectList.getSize() == 1 ) {
                        circumcircle_center = intersectList.getFirstEntity();
                    }
                    else {
                        intersectList.reset4Loop();
                        while ( intersectList.setNext4Loop() ) {
                            rg_Point2D intersectPt = intersectList.getEntity();
                            double radius = intersectPt.distance(point3);
                            if ( this_circumcircle_has_right_orientation( vertex, rg_Circle2D( intersectPt, radius ) ) ) {
                                circumcircle_center = intersectPt;
                                break;
                            }
                        }
                    }
                    double radius = circumcircle_center.distance(point3);
                    circumcircle_of_vertex = rg_Circle2D( circumcircle_center, radius );
                    */


                    ////////////////////// PARABOLA IMPLICIT FORM ////////////////////////////////
                    /*
                    double      radius_of_focus = 0.0;
                    rg_Point2D  focus           = point3;
                    rg_Line2D   directrix       = lineSeg2;
                    rg_Point2D  normalVec       = directrix.getNormalVector().getUnitVector();

                    directrix.setSP(directrix.getSP() - normalVec * radius_of_focus);
                    directrix.setEP(directrix.getEP() - normalVec * radius_of_focus);
                    Parabola2D parabola(focus, directrix);

                    computeCoordOfNewVVertex_on_parabola(parabola, disk1, radius_of_focus, tangentCircle[0], tangentCircle[1]);

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        circumcircle_of_vertex = tangentCircle[1];
                    }
                    */


                    ///////////////////////// IMPLICIT FORM & NUMERICAL //////////////////////////
                    /*
                    //rg_Point2D targetCenter = circumcircle_of_vertex.getCenterPt();
                    rg_Point2D targetCenter;

                    {
                        rg_Point2D init_tangentPt[3];
                        lineSeg2.compute_footprint_of_point_onto_line_segment(point3, init_tangentPt[0]);
                        lineSeg2.compute_footprint_of_point_onto_line_segment(disk1.getCenterPt(), init_tangentPt[1]);
                        init_tangentPt[0] = (init_tangentPt[0] + init_tangentPt[1]) / 2.0;
                        init_tangentPt[1] = point3;
                        disk1.compute_perpendicular_footprint_of_point_onto_circle(init_tangentPt[0], init_tangentPt[2]);

                        targetCenter = (init_tangentPt[0] + init_tangentPt[1] + init_tangentPt[2]) / 3.0;
                    }

                    rg_Point2D currentCenter;
                    rg_Circle2D refinedTangentCircle;
                    rg_INT numIterations = 0;
                    do
                    {
                        currentCenter = targetCenter;

                        rg_Point2D footprint[3];
                        lineSeg2.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[0]);
                        footprint[1] = point3;
                        disk1.compute_perpendicular_footprint_of_point_onto_circle(currentCenter, footprint[2]);

                        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
                        targetCenter = refinedTangentCircle.getCenterPt();
                        ++numIterations;
                    } while (targetCenter.distance(currentCenter) >= 0.001 && numIterations <= 20);

                    double radius = refinedTangentCircle.getCenterPt().distance(point3);
                    refinedTangentCircle.setRadius(radius);


                    circumcircle_of_vertex = refinedTangentCircle;
                    */

#ifdef DEBUG_VERTEX
                    lineSeg2.compute_perpendicular_footprint_of_point_onto_entire_line( circumcircle_of_vertex.getCenterPt(), foorprint_debug[0] );
                    foorprint_debug[1] = point3;
                    insertedDisk = disk1;
#endif
                }
                else
                {
#ifdef CHECK_COMP_TIME
                    ++num_DLV;
                    startTime = clock();
#endif
                    rg_Point2D SP = point3;
                    rg_Point2D dirVec = lineSeg2.getNormalVector();
                    dirVec = dirVec.getUnitVector();

                    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP, dirVec, disk1, tangentCircle[0] );
                    circumcircle_of_vertex = tangentCircle[0];
#ifdef CHECK_COMP_TIME
                    endTime = clock();
                    t_comp_V_DLV = t_comp_V_DLV + endTime - startTime + localCompTime;
#endif

#ifdef DEBUG_VERTEX
                    foorprint_debug[0] = point3;
                    foorprint_debug[1] = point3;
                    insertedDisk = disk1;
#endif
                }
            }
            break;
            case 0: // numDisks == 1, numLines == 0, numPoints == 2
            {
#ifdef CHECK_COMP_TIME
                ++num_DDD;
                startTime = clock();
#endif
                rg_Circle2D disk1;
                rg_Point2D point[2];
                rg_INT pointIndex = 0;
                for ( rg_INT i = 0; i < 3; i++ )
                {
                    if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                        disk1 = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                    if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                        ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point[pointIndex++] );
                }

                int numTangentCircles = rg_Circle2D::makeCircumcircle( disk1, rg_Circle2D( point[0], 0.0 ), rg_Circle2D( point[1], 0.0 ), tangentCircle[0], tangentCircle[1] );

#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_comp_V_DDD = t_comp_V_DDD + endTime - startTime;

                startTime = clock();
#endif
                if ( numTangentCircles == 1 ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                        rg_Circle2D tempCircle = tangentCircle[0];
                        tangentCircle[0] = tangentCircle[1];
                        tangentCircle[1] = tempCircle;
                    }

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        circumcircle_of_vertex = tangentCircle[1];
                    }
                }
#ifdef CHECK_COMP_TIME
                endTime = clock();
                t_check_orientation = t_check_orientation + endTime - startTime;
#endif

#ifdef DEBUG_VERTEX
                foorprint_debug[0] = point[0];
                foorprint_debug[1] = point[1];
                insertedDisk = disk1;
#endif
            }
            break;
            }
        }
        break;
        }

#ifdef DEBUG_VERTEX
        fout_debug_vertex << foorprint_debug[0].getX() << "\t" << foorprint_debug[0].getY() << endl;
        fout_debug_vertex << foorprint_debug[1].getX() << "\t" << foorprint_debug[1].getY() << endl;
        fout_debug_vertex << insertedDisk.getCenterPt().getX() << "\t" << insertedDisk.getCenterPt().getY() << "\t" << insertedDisk.getRadius() << endl;
        fout_debug_vertex << circumcircle_of_vertex.getCenterPt().getX() << "\t" << circumcircle_of_vertex.getCenterPt().getY() << "\t" << circumcircle_of_vertex.getRadius() << endl << endl << endl;
#endif
        vertex->setCircumcircle( circumcircle_of_vertex );
        vertex->setStatus( WHITE_V );
    }
#ifdef DEBUG_VERTEX
    fout_debug_vertex.close();
#endif
}



void PolygonVD2D_inRectangle::refine_location_of_VVertices()
{
    for ( list<VVertex2D*>::iterator i_vtx = m_VVertices.begin(); i_vtx != m_VVertices.end(); ++i_vtx ) {
        VVertex2D* vertex = *i_vtx;

        if ( vertex->isInfinite() )
            continue;

        if ( get_location_status_of_Vvertex( vertex ) == ON_POLYGON_BOUNDARY )
            continue;

        computeCoordOfNewVVertex( vertex );
    }
}



void PolygonVD2D_inRectangle::computeCoordOfNewVVertex( VVertex2D*& vertex )
{
    // We suppose that VVertices are classified by their status (OUTSIDE, ON, INSIDE of polygon).
    VDEntity_Location_Status status = get_location_status_of_Vvertex( vertex );

    int numDisks, numLines, numPoints;
    Generator2D* parentGens[3];

    numDisks = numLines = numPoints = 0;
    list<Generator2D*> diskGens;
    vertex->getDefining3Generators( diskGens );

    int index_gen = 0;
    for ( list<Generator2D*>::iterator i_gen = diskGens.begin(); i_gen != diskGens.end(); ++i_gen, ++index_gen ) {
        Generator2D* parentGen = (Generator2D*)( *i_gen )->getUserData();
        Generator2D::Generator_Type type = parentGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::DISK_G:
        {
            ++numDisks;
            parentGens[index_gen] = parentGen;
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            ++numLines;
            parentGens[index_gen] = parentGen;
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            ++numPoints;
            parentGens[index_gen] = parentGen;
        }
        break;

        default:
            break;
        }
    }

    rg_Circle2D circumcircle_of_vertex;
    rg_Circle2D tangentCircle[2];

    switch ( numDisks )
    {
    case 3: // numDisks == 3, numLines == 0, numPoints == 0
    {
        rg_Circle2D disks[3];
        disks[0] = ( (DiskGenerator2D*)parentGens[0] )->getDisk();
        disks[1] = ( (DiskGenerator2D*)parentGens[1] )->getDisk();
        disks[2] = ( (DiskGenerator2D*)parentGens[2] )->getDisk();

        int numTangentCircles = rg_Circle2D::makeCircumcircle( disks[0], disks[1], disks[2], tangentCircle[0], tangentCircle[1] );

        if ( numTangentCircles == 1 ) {
            circumcircle_of_vertex = tangentCircle[0];
        }
        else {
            if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                rg_Circle2D tempCircle = tangentCircle[0];
                tangentCircle[0] = tangentCircle[1];
                tangentCircle[1] = tempCircle;
            }

            if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                circumcircle_of_vertex = tangentCircle[0];
            }
            else {
                circumcircle_of_vertex = tangentCircle[1];
            }
        }
    }
    break;



    case 2:
    {
        // numDisks == 2, numLines == 1, numPoints == 0
        if ( numLines == 1 )
        {
            rg_Circle2D disks[2];
            rg_Line2D lineSeg3;
            rg_INT diskIndex = 0;
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                    disks[diskIndex++] = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                    ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg3 );
            }

            int numTangentCircles = 0;
            rg_Circle2D tangentCircle[2];
            if ( status == VDEntity_Location_Status::OUTSIDE_POLYGON ) {
                lineSeg3 = lineSeg3.get_reversed_line2D();
            }
            numTangentCircles = computeCoordOfNewVVertex_of_two_disks_and_a_line( disks[0], disks[1], lineSeg3, tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 1 ) {
                circumcircle_of_vertex = tangentCircle[0];
            }
            else {
                if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                    rg_Circle2D tempCircle = tangentCircle[0];
                    tangentCircle[0] = tangentCircle[1];
                    tangentCircle[1] = tempCircle;
                }

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
            }
        }

        // numDisks == 2, numLines == 0, numPoints == 1
        else
        {
            rg_Circle2D disks[2];
            rg_Point2D point3;
            rg_INT diskIndex = 0;
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                    disks[diskIndex++] = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                    ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point3 );
            }

            int numTangentCircles = rg_Circle2D::makeCircumcircle( disks[0], disks[1], rg_Circle2D( point3, 0.0 ), tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 1 ) {
                circumcircle_of_vertex = tangentCircle[0];
            }
            else {
                if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                    rg_Circle2D tempCircle = tangentCircle[0];
                    tangentCircle[0] = tangentCircle[1];
                    tangentCircle[1] = tempCircle;
                }

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
            }
        }
    }
    break;



    case 1:
    {
        switch ( numLines )
        {
        case 2: // numDisks == 1, numLines == 2, numPoints == 0
        {
            rg_Circle2D disk1;
            rg_Line2D lineSeg[2];
            rg_INT lineSegIndex = 0;
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                    disk1 = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                    ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg[lineSegIndex++] );
            }

            if ( status == VDEntity_Location_Status::OUTSIDE_POLYGON ) {
                lineSeg[0] = lineSeg[0].get_reversed_line2D();
                lineSeg[1] = lineSeg[1].get_reversed_line2D();
            }

            VEdge2D* edge_of_two_polygon_edges = NULL;
            list<VEdge2D*> incidentEdges;
            vertex->getIncident3VEdges( incidentEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;

                Generator2D* leftParent = (Generator2D*)currEdge->getLeftFace()->getGenerator()->getUserData();
                Generator2D* rightParent = (Generator2D*)currEdge->getRightFace()->getGenerator()->getUserData();

                if ( leftParent->getType() == Generator2D::Generator_Type::EDGE_G && rightParent->getType() == Generator2D::Generator_Type::EDGE_G ) {
                    edge_of_two_polygon_edges = currEdge;
                    break;
                }
            }

            VVertex2D* vertex_which_has_coordinate = NULL;
            if ( edge_of_two_polygon_edges->getStartVertex() == vertex ) {
                vertex_which_has_coordinate = edge_of_two_polygon_edges->getEndVertex();
            }
            else {
                vertex_which_has_coordinate = edge_of_two_polygon_edges->getStartVertex();
            }

            rg_Point2D dirVec_line1 = lineSeg[0].getEP() - lineSeg[0].getSP();
            rg_Point2D dirVec_line2 = lineSeg[1].getSP() - lineSeg[1].getEP();
            dirVec_line1 = dirVec_line1.getUnitVector();
            dirVec_line2 = dirVec_line2.getUnitVector();

            rg_Point2D dirVec_edge;
            double cross_product = dirVec_line1 * dirVec_line2;
            if ( rg_ZERO( cross_product ) ) {
                dirVec_edge = dirVec_line1;
            }
            else if ( rg_POS( cross_product ) ) {
                dirVec_edge = dirVec_line1 + dirVec_line2;
                dirVec_edge = dirVec_edge.getUnitVector();
            }
            else { // if ( rg_NE( cross_product ) ) {
                dirVec_edge = -dirVec_line1 - dirVec_line2;
                dirVec_edge = dirVec_edge.getUnitVector();
            }

            rg_Point2D SP = vertex_which_has_coordinate->getLocation();
            rg_Point2D EP = SP + dirVec_edge;
            rg_Line2D geometry_edge( SP, EP );

            rg_Point2D footprint_SP, footprint_EP;
            lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( SP, footprint_SP );
            lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( EP, footprint_EP );

            double radius_SP = SP.distance( footprint_SP );
            double radius_EP = EP.distance( footprint_EP );

            computeCoordOfNewVVertex_on_line_of_two_polygon_edges( geometry_edge, rg_Circle2D( SP, radius_SP ), rg_Circle2D( EP, radius_EP ), disk1, tangentCircle[0], tangentCircle[1] );

            if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                rg_Circle2D tempCircle = tangentCircle[0];
                tangentCircle[0] = tangentCircle[1];
                tangentCircle[1] = tempCircle;
            }

            if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                circumcircle_of_vertex = tangentCircle[0];
            }
            else {
                circumcircle_of_vertex = tangentCircle[1];
            }
        }
        break;

        case 1: // numDisks == 1, numLines == 1, numPoints == 1
        {
            EdgeGenerator2D* edgeGen = rg_NULL;
            VertexGenerator2D* vertexGen = rg_NULL;

            rg_Circle2D disk1;
            rg_Line2D  lineSeg2;
            rg_Point2D point3;
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() ) {
                    disk1 = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                }
                if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                {
                    edgeGen = (EdgeGenerator2D*)parentGens[i];
                    ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg2 );
                }
                if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                {
                    vertexGen = (VertexGenerator2D*)parentGens[i];
                    ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point3 );
                }
            }

            if ( status == VDEntity_Location_Status::OUTSIDE_POLYGON ) {
                lineSeg2 = lineSeg2.get_reversed_line2D();
            }

            if ( vertexGen->get_next_edge_generator() != edgeGen && vertexGen->get_previous_edge_generator() != edgeGen ) {
                int numTangentCircles = 0;
                rg_Circle2D tangentCircle[2];
                numTangentCircles = computeCoordOfNewVVertex_of_two_disks_and_a_line( rg_Circle2D( point3, 0.0 ), disk1, lineSeg2, tangentCircle[0], tangentCircle[1] );

                if ( numTangentCircles == 1 ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                        rg_Circle2D tempCircle = tangentCircle[0];
                        tangentCircle[0] = tangentCircle[1];
                        tangentCircle[1] = tempCircle;
                    }

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        circumcircle_of_vertex = tangentCircle[1];
                    }
                }
            }
            else
            {
                rg_Point2D SP = point3;
                rg_Point2D dirVec = lineSeg2.getNormalVector();
                dirVec = dirVec.getUnitVector();

                computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP, dirVec, disk1, tangentCircle[0] );
                circumcircle_of_vertex = tangentCircle[0];
            }
        }
        break;
        case 0: // numDisks == 1, numLines == 0, numPoints == 2
        {
            rg_Circle2D disk1;
            rg_Point2D point[2];
            rg_INT pointIndex = 0;
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::DISK_G == parentGens[i]->getType() )
                    disk1 = ( (DiskGenerator2D*)parentGens[i] )->getDisk();
                if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                    ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point[pointIndex++] );
            }

            int numTangentCircles = rg_Circle2D::makeCircumcircle( disk1, rg_Circle2D( point[0], 0.0 ), rg_Circle2D( point[1], 0.0 ), tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 1 ) {
                circumcircle_of_vertex = tangentCircle[0];
            }
            else {
                if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                    rg_Circle2D tempCircle = tangentCircle[0];
                    tangentCircle[0] = tangentCircle[1];
                    tangentCircle[1] = tempCircle;
                }

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
            }
        }
        break;
        }
    }
    break;



    case 0:
    {
        switch ( numLines )
        {
        case 3: // numDisks == 0, numLines == 3, numPoints == 0
        {
            rg_Line2D lineSeg[3];
            ( (EdgeGenerator2D*)parentGens[0] )->get_geometry( lineSeg[0] );
            ( (EdgeGenerator2D*)parentGens[1] )->get_geometry( lineSeg[1] );
            ( (EdgeGenerator2D*)parentGens[2] )->get_geometry( lineSeg[2] );

            if ( status == VDEntity_Location_Status::OUTSIDE_POLYGON ) {
                lineSeg[0] = lineSeg[0].get_reversed_line2D();
                lineSeg[1] = lineSeg[1].get_reversed_line2D();
                lineSeg[2] = lineSeg[2].get_reversed_line2D();
            }
            circumcircle_of_vertex = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators( lineSeg[0], lineSeg[1], lineSeg[2], vertex->getCircumcircle() );
        }
        break;
        case 2: // numDisks == 0, numLines == 2, numPoints == 1
        {
            // THERE ARE TWO CASES:
            // CASE I: VERTEX IS OFF BOTH EDGES 
            //	CASE I-1: WE WANT SMALL CIRCLE
            //	CASE I-2: WE WANT LARGE CIRCLE
            // CASE II: VERTEX BOUNDS ONE OF TWO EDGES ON THE EDGE
            EdgeGenerator2D* edgeGen[2] = { rg_NULL, rg_NULL };
            VertexGenerator2D* vertexGen = rg_NULL;

            rg_Line2D  lineSeg[2];
            rg_INT lineSegindex = 0;
            rg_Point2D point3;
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                {
                    edgeGen[lineSegindex] = (EdgeGenerator2D*)parentGens[i];
                    ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg[lineSegindex++] );
                }
                if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                {
                    ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point3 );
                    vertexGen = (VertexGenerator2D*)parentGens[i];
                }
            }

            if ( status == VDEntity_Location_Status::OUTSIDE_POLYGON ) {
                lineSeg[0] = lineSeg[0].get_reversed_line2D();
                lineSeg[1] = lineSeg[1].get_reversed_line2D();
            }

            bool bReflexVertexBoundsEdge = true;
            rg_Line2D perpendicularBisector;
            if ( vertexGen->get_next_edge_generator() == edgeGen[0] || vertexGen->get_previous_edge_generator() == edgeGen[0] )
            {
                perpendicularBisector.setSP( point3 );
                perpendicularBisector.setEP( point3 + lineSeg[0].getNormalVector() );
            }
            else if ( vertexGen->get_next_edge_generator() == edgeGen[1] || vertexGen->get_previous_edge_generator() == edgeGen[1] )
            {
                perpendicularBisector.setSP( point3 );
                perpendicularBisector.setEP( point3 + lineSeg[1].getNormalVector() );
            }
            else
            {
                bReflexVertexBoundsEdge = false;
            }

            // this code should be rewritten as follows:
            // precise method: compute intersection between both "bisector line" of two lineSeg generators 
            // and "bisector parabola" of lineSeg and point generators

            // CASE II
            if ( bReflexVertexBoundsEdge )
            {
                rg_Line2D bisectorLine = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0].get_reversed_line2D(), lineSeg[1] );
                bool bTwoLinesAreParallel = false;
                rg_Point2D refinedCenter = bisectorLine.compute_intersection_with_line( perpendicularBisector, bTwoLinesAreParallel );
                rg_REAL    refinedRadius = refinedCenter.distance( point3 );
                circumcircle_of_vertex.setCircle( refinedCenter, refinedRadius );
            }
            // CASE I
            else
                circumcircle_of_vertex = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators( lineSeg[0], lineSeg[1], point3, vertex->getCircumcircle() );
        }
        break;
        case 1: // numDisks == 0, numLines == 1, numPoints == 2
        {
            // THERE ARE THREE CASES 
            // CASE I: ADD THE CASE THAT LINE SEPERATES TWO POINTS: INFEASIBLE CASE

            EdgeGenerator2D* edgeGen = rg_NULL;
            VertexGenerator2D* vertexGen[2] = { rg_NULL, rg_NULL };

            rg_Line2D  lineSeg;
            rg_INT pointIndex = 0;
            rg_Point2D point[2];
            for ( rg_INT i = 0; i < 3; i++ )
            {
                if ( Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType() )
                {
                    edgeGen = (EdgeGenerator2D*)parentGens[i];
                    ( (EdgeGenerator2D*)parentGens[i] )->get_geometry( lineSeg );
                }
                if ( Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType() )
                {
                    vertexGen[pointIndex] = (VertexGenerator2D*)parentGens[i];
                    ( (VertexGenerator2D*)parentGens[i] )->get_geometry( point[pointIndex++] );
                }
            }

            if ( status == VDEntity_Location_Status::OUTSIDE_POLYGON ) {
                lineSeg = lineSeg.get_reversed_line2D();
            }

            bool bfirstReflexVertexBoundsEdge = false;
            bool bSecondReflexVertexBoundsEdge = false;
            rg_Line2D perpendicularBisectorBetweenReflexVertexAndEdge;
            if ( vertexGen[0]->get_next_edge_generator() == edgeGen || vertexGen[0]->get_previous_edge_generator() == edgeGen )
            {
                perpendicularBisectorBetweenReflexVertexAndEdge.setSP( point[0] );
                perpendicularBisectorBetweenReflexVertexAndEdge.setEP( point[0] + lineSeg.getNormalVector() );
                bfirstReflexVertexBoundsEdge = true;
            }
            if ( vertexGen[1]->get_next_edge_generator() == edgeGen || vertexGen[1]->get_previous_edge_generator() == edgeGen )
            {
                perpendicularBisectorBetweenReflexVertexAndEdge.setSP( point[1] );
                perpendicularBisectorBetweenReflexVertexAndEdge.setEP( point[1] + lineSeg.getNormalVector() );
                bSecondReflexVertexBoundsEdge = true;
            }

            // CASE II: ONE VERTEX BOUNDS EDGE(LINESEGMENT) ON THE EDGE
            if ( bfirstReflexVertexBoundsEdge || bSecondReflexVertexBoundsEdge )
            {
                rg_Line2D perpendicularBisectorBetweenTwoReflexVertices = rg_GeoFunc::compute_bisector_line_between_two_points( point[0], point[1] );
                bool bTwoLinesAreParllel = false;
                rg_Point2D refinedCenter = perpendicularBisectorBetweenReflexVertexAndEdge.compute_intersection_with_line( perpendicularBisectorBetweenTwoReflexVertices, bTwoLinesAreParllel );

                rg_Point2D edgeDirVec = lineSeg.evaluateVector();
                rg_Point2D perpendicularBisectorTowardPolygonInside = refinedCenter - perpendicularBisectorBetweenReflexVertexAndEdge.getSP();
                if ( edgeDirVec.operator*( perpendicularBisectorTowardPolygonInside ) > 0.0 )
                {
                    rg_REAL    refinedRadius = refinedCenter.distance( point[0] );
                    circumcircle_of_vertex.setCircle( refinedCenter, refinedRadius );
                }
                else
                {
                    // MAYBE THIS CASE WILL NOT BE ENCOUNTERED ... REMOVED THIS BLOCK
                    //AfxMessageBox(_T("Orthogonalize VEdge by projection"));

                    if ( bfirstReflexVertexBoundsEdge )
                        orthogonalize_VEdge_by_projecting_VVertex( const_cast<VVertex2D*>( vertex ), edgeGen, point[0] );
                    else
                        orthogonalize_VEdge_by_projecting_VVertex( const_cast<VVertex2D*>( vertex ), edgeGen, point[1] );

                    circumcircle_of_vertex = vertex->getCircumcircle();
                }
            }
            // CASE III: TWO VERTICES ARE OFF THE EDGE(LINESEGMENT) IN THE SAME HALFSPACE
            else
            {
                //refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(lineSeg, point[0], point[1], vertex->getCircumcircle());
                int numTangentCircles = 0;
                rg_Circle2D tangentCircle[2];
                numTangentCircles = computeCoordOfNewVVertex_of_two_disks_and_a_line( rg_Circle2D( point[0], 0.0 ), rg_Circle2D( point[1], 0.0 ), lineSeg, tangentCircle[0], tangentCircle[1] );

                if ( numTangentCircles == 1 ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                        rg_Circle2D tempCircle = tangentCircle[0];
                        tangentCircle[0] = tangentCircle[1];
                        tangentCircle[1] = tempCircle;
                    }

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        circumcircle_of_vertex = tangentCircle[0];
                    }
                    else {
                        circumcircle_of_vertex = tangentCircle[1];
                    }
                }
            }

            // CASE IV: ONE VERTEX BOUNDS EDGE(LINESEGMENT) ON THE EDGE, THE OTHER VERTEX DOES NOT BOUND THE EDGE BUT ON THE LINE EXTENDING THE EDGE
            // THIS CASE CAN EXIST BUT INFEASIBLE SOLUTION
        }
        break;
        case 0: // numberOfEllipses == 0, numberOfLineSegs == 0, numberOfPoints == 3
        {
            rg_Point2D point[3];
            ( (VertexGenerator2D*)parentGens[0] )->get_geometry( point[0] );
            ( (VertexGenerator2D*)parentGens[1] )->get_geometry( point[1] );
            ( (VertexGenerator2D*)parentGens[2] )->get_geometry( point[2] );
            rg_Circle2D nullCircle;
            rg_Circle2D::makeCircumcircle( rg_Circle2D( point[0], 0.0 ), rg_Circle2D( point[1], 0.0 ), rg_Circle2D( point[2], 0.0 ), circumcircle_of_vertex, nullCircle );
        }
        break;
        }
    }
    break;
    default:
    break;
    }

    vertex->setCircumcircle( circumcircle_of_vertex );
}



void PolygonVD2D_inRectangle::replicatePolygonEntitiesAsVDGenerators( list<Generator2D*>& newGenerators )
{
    /*
    Polygon boundary entities(P_Vertex, P_Edge)들을 VD를 위한 generator로 변환하고 기본적인 topology를 연결해 준다.
    */
    list<rg_Point2D> boundaryVertices;
    m_polygons.back().get_boundary_vertices( boundaryVertices );

    list<rg_Point2D>::iterator i_boundaryVertices = boundaryVertices.begin();
    list<rg_Point2D>::iterator j_boundaryVertices = ++boundaryVertices.begin();

    vector<VertexGenerator2D*> vertexGenArray;
    vector<EdgeGenerator2D*>   edgeGenArray;

    vertexGenArray.resize( boundaryVertices.size() + 1 );
    edgeGenArray.resize( boundaryVertices.size() );

    int index_gen = 0;
    int memberGeneratorUserID = 0;
    while ( j_boundaryVertices != boundaryVertices.end() )
    {
        rg_Point2D startPoint = *i_boundaryVertices;
        rg_Point2D endPoint = *j_boundaryVertices;

        VertexGenerator2D* vertexGen = new VertexGenerator2D( startPoint, memberGeneratorUserID++ );
        m_polygonGenerators.push_back( vertexGen );
        vertexGenArray[index_gen] = vertexGen;
        m_mapGenerator2Polygon[vertexGen] = &m_polygons.back();
        newGenerators.push_back( vertexGen );

        EdgeGenerator2D* edgeGen = new EdgeGenerator2D( startPoint, endPoint, memberGeneratorUserID++ );
        m_polygonGenerators.push_back( edgeGen );
        edgeGenArray[index_gen] = edgeGen;
        m_mapGenerator2Polygon[edgeGen] = &m_polygons.back();
        newGenerators.push_back( edgeGen );

        ++i_boundaryVertices;
        ++j_boundaryVertices;
        ++index_gen;
    }

    rg_Point2D startPoint = boundaryVertices.back();
    rg_Point2D endPoint = boundaryVertices.front();
    VertexGenerator2D* vertexGen = new VertexGenerator2D( startPoint, memberGeneratorUserID++ );
    m_polygonGenerators.push_back( vertexGen );
    newGenerators.push_back( vertexGen );
    vertexGenArray[index_gen] = vertexGen;
    vertexGenArray[index_gen + 1] = vertexGenArray[0];

    EdgeGenerator2D* edgeGen = new EdgeGenerator2D( startPoint, endPoint, memberGeneratorUserID );
    m_polygonGenerators.push_back( edgeGen );
    newGenerators.push_back( edgeGen );
    edgeGenArray[index_gen] = edgeGen;


    // VertexGenerator와 EdgeGenerator 사이의 topology를 연결해 준다.
    int numBoundaryVertices = boundaryVertices.size();
    for ( int i = 0; i < numBoundaryVertices; ++i ) {
        vertexGenArray[i]->set_next_edge_generator( edgeGenArray[i] );
        vertexGenArray[i + 1]->set_previous_edge_generator( edgeGenArray[i] );

        edgeGenArray[i]->set_start_vertex_generator( vertexGenArray[i] );
        edgeGenArray[i]->set_end_vertex_generator( vertexGenArray[i + 1] );
    }
}



void PolygonVD2D_inRectangle::generateChildrenDisksOfPolygonGenerators( const list<Generator2D*>& newGenerators )
{
    double minEdgeLength = m_polygons.back().find_minimum_edge_length();

    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        Generator2D::Generator_Type type = currGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

            if ( !vertexGen->isThisFromContainer() ) {
                vertexGen->set_child_disk( rg_Circle2D( vertexGen->get_point(), minEdgeLength / 4.0 ) );
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)currGen;

            if ( !edgeGen->isThisFromContainer() ) {
                list<rg_Circle2D> childrenDisks;
                generate_children_disks_of_polygon_edge_generator( edgeGen, minEdgeLength, childrenDisks );
                edgeGen->set_children_disks( childrenDisks );
            }
        }
        break;

        default:
            break;
        }
    }
}



void PolygonVD2D_inRectangle::insertChildrenDisksIntoVoronoiDiagramOfRectangularContainer( const list<Generator2D*>& newGenerators )
{
    list< rg_Triplet<rg_Circle2D, int, void*> > disk_UserID_UserData_TripletList;
    int childDiskID = m_polygonGenerators.back()->getID();


    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        Generator2D::Generator_Type type = currGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

            if ( vertexGen->isThisFromContainer() ) {
                break;
            }

            vertexGen->setID( childDiskID + 1 );
            rg_Circle2D childDisk = vertexGen->get_child_disk();
            disk_UserID_UserData_TripletList.push_back( rg_Triplet<rg_Circle2D, int, void*>( childDisk, ++childDiskID, vertexGen ) );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)currGen;

            if ( edgeGen->isThisFromContainer() ) {
                break;
            }

            edgeGen->setID( childDiskID + 1 );
            list<rg_Circle2D>& childrenDisks = edgeGen->get_children_disks();
            for ( list<rg_Circle2D>::iterator i_childDisk = childrenDisks.begin(); i_childDisk != childrenDisks.end(); ++i_childDisk ) {
                disk_UserID_UserData_TripletList.push_back( rg_Triplet<rg_Circle2D, int, void*>( *i_childDisk, ++childDiskID, edgeGen ) );
            }
        }
        break;

        default:
            break;
        }
    }


    list< pair< DiskGenerator2D*, void* > > childrenGens_to_parentGen;
    for ( list< rg_Triplet<rg_Circle2D, int, void*> >::iterator i_diskTriplet = disk_UserID_UserData_TripletList.begin(); i_diskTriplet != disk_UserID_UserData_TripletList.end(); ++i_diskTriplet ) {
        rg_Triplet<rg_Circle2D, int, void*> currTriplet = *i_diskTriplet;
        DiskGenerator2D* currGen = insertOneChildDiskInPolygonVoronoiDiagram( currTriplet );
        childrenGens_to_parentGen.push_back( make_pair( currGen, currTriplet.m_third ) );
    }

    for ( list< pair< DiskGenerator2D*, void* > >::iterator i_genPair = childrenGens_to_parentGen.begin(); i_genPair != childrenGens_to_parentGen.end(); ++i_genPair ) {
        i_genPair->first->setUserData( i_genPair->second );
    }
}



void PolygonVD2D_inRectangle::mergeVFacesOfChildrenDisksOnEachPEdge( const list<Generator2D*>& newGenerators )
{
    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() != Generator2D::Generator_Type::EDGE_G ) {
            continue;
        }

        EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)currGen;

        if ( edgeGen->isThisFromContainer() ) {
            continue;
        }

        rg_INT firstSonID = edgeGen->getID();
        VFace2D* firstSonVFace = getGeneratorWhichHasThisID( firstSonID )->getOuterFace();

        int numChildrenDisks = edgeGen->get_children_disks().size();

        for ( int i = 1; i < numChildrenDisks; ++i ) {
            VFace2D* currentVFace = getGeneratorWhichHasThisID( firstSonID + i )->getOuterFace();
            mergeSecondVFaceToTheFirstIfTheyAreAdjacent( firstSonVFace, currentVFace );
        }
    }

    removeAllExtraneousVVerticesAndVEdges();
}



void PolygonVD2D_inRectangle::relocate_some_boundary_Vvertices_to_vertex_generator_centers( const list<Generator2D*>& newGenerators )
{
    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() != Generator2D::Generator_Type::VERTEX_G )
            continue;

        VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

        if ( vertexGen->isThisFromContainer() ) {
            continue;
        }

        EdgeGenerator2D* prevEdgeGen = vertexGen->get_previous_edge_generator();
        EdgeGenerator2D* nextEdgeGen = vertexGen->get_next_edge_generator();

        set_type_of_PVertex( vertexGen, prevEdgeGen, nextEdgeGen );

        VEdge2D* VEdge_between_prevEdgeGen_N_vertexGen = rg_NULL;
        VEdge2D* VEdge_between_vertexGen_N_nextEdgeGen = rg_NULL;
        find_two_VEdges_between_vertex_generator_N_incident_edge_generators( vertexGen, prevEdgeGen, nextEdgeGen, VEdge_between_prevEdgeGen_N_vertexGen, VEdge_between_vertexGen_N_nextEdgeGen );

        VVertex2D* startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen = VEdge_between_prevEdgeGen_N_vertexGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen = VEdge_between_prevEdgeGen_N_vertexGen->getEndVertex();
        VVertex2D* startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getEndVertex();


        VFace2D* childVFaceOfVertexGen = getGeneratorWhichHasThisID( vertexGen->getID() )->getOuterFace();
        VVertex2D* commonVVertexToBeRelocated = rg_NULL;
        if ( vertexGen->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE )
        {
            if ( VEdge_between_prevEdgeGen_N_vertexGen->getLeftFace() == childVFaceOfVertexGen )
            {
                if ( endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen
                    || endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
                    commonVVertexToBeRelocated = endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
            }
            else
            {
                if ( startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen
                    || startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
                    commonVVertexToBeRelocated = startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
            }
        }
        else
        {
            if ( VEdge_between_prevEdgeGen_N_vertexGen->getLeftFace() == childVFaceOfVertexGen )
            {
                if ( startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen
                    || startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
                    commonVVertexToBeRelocated = startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
            }
            else
            {
                if ( endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen
                    || endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
                    commonVVertexToBeRelocated = endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
            }
        }

        if ( commonVVertexToBeRelocated != rg_NULL )
        {
            commonVVertexToBeRelocated->setCircumcircle( rg_Circle2D( vertexGen->get_point(), 0.0 ) );
        }
        else
        {
            find_and_relocate_two_boundary_Vvertices_to_vertex_generator_center( vertexGen, prevEdgeGen, nextEdgeGen );
        }
    }
}



void PolygonVD2D_inRectangle::secure_VFace_topology_of_PVertices_by_flipping_the_mating_quill_VEdge_of_relocated_VVertex( const list<Generator2D*>& newGenerators )
{
    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() != Generator2D::Generator_Type::VERTEX_G )
            continue;

        VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

        if ( vertexGen->isThisFromContainer() ) {
            continue;
        }

        VFace2D* childVFaceOfVertexGen = getGeneratorWhichHasThisID( vertexGen->getID() )->getOuterFace();

        list<VEdge2D*>  boundingVEdgesOfVertexGen;
        childVFaceOfVertexGen->getBoundaryVEdges( boundingVEdgesOfVertexGen );

        if ( boundingVEdgesOfVertexGen.size() > 2 )
            continue;

        EdgeGenerator2D* prevPEdgeGen = vertexGen->get_previous_edge_generator();
        EdgeGenerator2D* nextPEdgeGen = vertexGen->get_next_edge_generator();

        VEdge2D* VEdge_between_prevPEdgeGen_N_vertexGen = rg_NULL;
        VEdge2D* VEdge_between_vertexGen_N_nextPEdgeGen = rg_NULL;
        find_two_VEdges_between_vertex_generator_N_incident_edge_generators( vertexGen, prevPEdgeGen, nextPEdgeGen, VEdge_between_prevPEdgeGen_N_vertexGen, VEdge_between_vertexGen_N_nextPEdgeGen );

        VEdge2D* targetVEdge = rg_NULL;
        if ( vertexGen->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE )
        {
            if ( VEdge_between_prevPEdgeGen_N_vertexGen->getLeftFace() == childVFaceOfVertexGen )
            {
                targetVEdge = VEdge_between_prevPEdgeGen_N_vertexGen->getRightLeg();
            }
            else
            {
                targetVEdge = VEdge_between_prevPEdgeGen_N_vertexGen->getLeftHand();
            }
        }
        else
        {
            if ( VEdge_between_prevPEdgeGen_N_vertexGen->getLeftFace() == childVFaceOfVertexGen )
            {
                targetVEdge = VEdge_between_prevPEdgeGen_N_vertexGen->getRightHand();
            }
            else
            {
                targetVEdge = VEdge_between_prevPEdgeGen_N_vertexGen->getLeftLeg();
            }
        }
        targetVEdge->flip();
    }
}



void PolygonVD2D_inRectangle::connect_parent_generator_to_first_son_disk_generator( const list<Generator2D*>& newGenerators )
{
    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        Generator2D::Generator_Type type = currGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;
            if ( vertexGen->isThisFromContainer() ) {
                break;
            }
            DiskGenerator2D* firstSonDiskGenerator = getGeneratorWhichHasThisID( vertexGen->getID() );
            vertexGen->set_first_son_disk_generator( firstSonDiskGenerator );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)currGen;
            if ( edgeGen->isThisFromContainer() ) {
                break;
            }
            DiskGenerator2D* firstSonDiskGenerator = getGeneratorWhichHasThisID( edgeGen->getID() );
            edgeGen->set_first_son_disk_generator( firstSonDiskGenerator );
        }
        break;

        default:
            break;
        }
    }
}



void PolygonVD2D_inRectangle::identify_internal_and_on_and_external_Vvertices_and_Vedges( const list<Generator2D*>& newGenerators )
{
    list< pair< VertexGenerator2D*, VVertex2D* > > PVertex_N_onVVertexPair;
    identify_Vvertices_on_boundary( newGenerators, PVertex_N_onVVertexPair );

    stack<VEdge2D*> outerEdges;
    stack<VEdge2D*> innerEdges;

    // collect initial stacks of outer edges and inner edges.
    for ( list< pair< VertexGenerator2D*, VVertex2D* > >::iterator i_pair = PVertex_N_onVVertexPair.begin(); i_pair != PVertex_N_onVVertexPair.end(); ++i_pair ) {
        VertexGenerator2D* currPVertex = i_pair->first;
        VVertex2D* currOnVVertex = i_pair->second;

        EdgeGenerator2D* prevPEdge = currPVertex->get_previous_edge_generator();
        EdgeGenerator2D* nextPEdge = currPVertex->get_next_edge_generator();

        VertexGenerator2D::Vertex_Type vertexType = currPVertex->vertex_type();

        // find THAT edge
        list<VEdge2D*> incidentEdges;
        currOnVVertex->getIncident3VEdges( incidentEdges );
        for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
            VEdge2D* currEdge = *i_edge;

            bool b_this_edge_is_defined_by_two_PEdges = true;
            if ( currEdge->getLeftFace()->getGenerator()->getUserData() == currPVertex || currEdge->getRightFace()->getGenerator()->getUserData() == currPVertex ) {
                b_this_edge_is_defined_by_two_PEdges = false;
            }
            else {
                b_this_edge_is_defined_by_two_PEdges = true;
            }


            if ( vertexType == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
                if ( b_this_edge_is_defined_by_two_PEdges ) {
                    outerEdges.push( currEdge );
                }
                else {
                    innerEdges.push( currEdge );
                }
            }
            else {
                if ( b_this_edge_is_defined_by_two_PEdges ) {
                    innerEdges.push( currEdge );
                }
                else {
                    outerEdges.push( currEdge );
                }
            }
        }
    }

    // outer edges
    while ( !outerEdges.empty() )
    {
        VEdge2D* outerEdge = outerEdges.top();
        outerEdges.pop();

        VVertex2D* startVertex = outerEdge->getStartVertex();
        VVertex2D* endVertex = outerEdge->getEndVertex();
        VVertex_Of_FirstSonVD_To_Location_Status_Map::iterator i_startVertexToLocationStatus = m_VVertexToLocationStatus.find( startVertex );
        VVertex_Of_FirstSonVD_To_Location_Status_Map::iterator i_endVertexToLocationStatus = m_VVertexToLocationStatus.find( endVertex );

        // If start vertex was visited before
        if ( i_startVertexToLocationStatus != m_VVertexToLocationStatus.end() )
        {
            if ( i_endVertexToLocationStatus != m_VVertexToLocationStatus.end() )
            {
                m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );
                continue;
            }
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( endVertex, OUTSIDE_POLYGON ) );
            m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            endVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }

        // If end vertex was visited before
        if ( i_endVertexToLocationStatus != m_VVertexToLocationStatus.end() )
        {
            if ( i_startVertexToLocationStatus != m_VVertexToLocationStatus.end() )
            {
                m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );
                continue;
            }
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( startVertex, OUTSIDE_POLYGON ) );
            m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            startVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }
    }

    set_location_status_of_remaining_Vvertices_and_Vedges_with_inside();
}



void PolygonVD2D_inRectangle::identify_Vvertices_on_boundary( const list<Generator2D*>& newGenerators, list<pair<VertexGenerator2D*, VVertex2D*>>& PVertex_N_onVVertexPair )
{
    for ( list<Generator2D*>::const_iterator i_gen = newGenerators.begin(); i_gen != newGenerators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() != Generator2D::Generator_Type::VERTEX_G )
            continue;

        VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

        if ( vertexGen->isThisFromContainer() ) {
            continue;
        }

        EdgeGenerator2D* prevEdgeGen = vertexGen->get_previous_edge_generator();
        EdgeGenerator2D* nextEdgeGen = vertexGen->get_next_edge_generator();

        VEdge2D* VEdge_between_prevEdgeGen_N_vertexGen = rg_NULL;
        VEdge2D* VEdge_between_vertexGen_N_nextEdgeGen = rg_NULL;
        find_two_VEdges_between_vertex_generator_N_incident_edge_generators( vertexGen, prevEdgeGen, nextEdgeGen, VEdge_between_prevEdgeGen_N_vertexGen, VEdge_between_vertexGen_N_nextEdgeGen );

        VVertex2D* startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen = VEdge_between_prevEdgeGen_N_vertexGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen = VEdge_between_prevEdgeGen_N_vertexGen->getEndVertex();
        VVertex2D* startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getEndVertex();

        VVertex2D* commonVVertex = rg_NULL;
        if ( startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen || startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
        {
            commonVVertex = startVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
        }
        else if ( endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen || endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen == endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen )
        {
            commonVVertex = endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen;
        }
        else
        {
        }

        if ( commonVVertex != rg_NULL )
        {
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( commonVVertex, ON_POLYGON_BOUNDARY ) );
            PVertex_N_onVVertexPair.push_back( make_pair( vertexGen, commonVVertex ) );
        }
    }
}



void PolygonVD2D_inRectangle::refine_location_of_VVertices( const list<Generator2D*>& newGenerators )
{
}



void PolygonVD2D_inRectangle::adjust_topology_of_interior_and_exterior_of_polygon_by_edge_flipping( const list<Generator2D*>& newGenerators )
{
    for ( list<Generator2D*>::const_iterator i_generators = newGenerators.begin(); i_generators != newGenerators.end(); ++i_generators ) {
        Generator2D* currGenerator = *i_generators;
        Generator2D::Generator_Type type = currGenerator->getType();

        if ( type == Generator2D::Generator_Type::VERTEX_G )
        {
            VertexGenerator2D* currVertexGen = (VertexGenerator2D*)currGenerator;

            //if ( currVertexGen->vertex_type() != VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE )
            //{
            //    continue;
            //}
            flatten_topology_of_VFace_boundary_of_reflex_vertex_gen_with_possible_flipping_out_of_VEdges( currVertexGen );
        }
        else if ( type == Generator2D::Generator_Type::EDGE_G ) {
            EdgeGenerator2D* currEdgeGen = (EdgeGenerator2D*)currGenerator;
            flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges( currEdgeGen );
        }
    }
}





