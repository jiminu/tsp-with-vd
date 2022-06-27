#include "PolygonVD2D.h"
using namespace V::GeometryTier;

#include "rg_GeoFunc.h"
#include "Parabola2D.h"
#include "rg_TMatrix2D.h"
#include "rg_dList.h"
#include "BucketForDisk.h"
#include "BucketForChildrenDiskGeneration.h"
#include <time.h>
#include <fstream>
#if defined( DEBUG_VERTEX ) || defined ( CHECK_COMP_TIME )
#include <fstream>
#include <time.h>
using namespace std;
#endif


PolygonVD2D::PolygonVD2D()
{
    m_rectContainer             = NULL;
    m_ptr_last_inserted_polygon = NULL;
    b_there_is_container        = false;
    b_VD_is_constructed         = false;
}



PolygonVD2D::~PolygonVD2D()
{
    if ( m_rectContainer != NULL ) {
        delete m_rectContainer;
    }
    m_rectContainer = NULL;
}



#ifdef CHECK_COMP_TIME
void PolygonVD2D::init_()
{
    t_wave_propagation = 0.0;
    t_anomaly_check = 0.0;
    t_find_red_neighbor = 0.0;
    t_make_new_cell = 0.0;
    t_mark_in = 0.0;
    t_connect_topology = 0.0;
    t_merge_split = 0.0;
    t_comp_vertex = 0.0;
    t_others = 0.0;

    num_DDD = 0;
    num_DDL = 0;
    num_DLV = 0;
    num_DLL = 0;
    t_comp_V_DDD = 0.0;
    t_comp_V_DDL = 0.0;
    t_comp_V_DLV = 0.0;
    t_comp_V_DLL = 0.0;
}

void PolygonVD2D::print_()
{
    ofstream fout( "comp_time_polyVD.txt" );
    fout << t_wave_propagation << "\t"
        << t_make_new_cell << "\t"
        << t_mark_in << "\t"
        << t_connect_topology << "\t"
        << t_merge_split << "\t"
        << t_comp_vertex << "\t"
        << t_others << endl;

    fout << num_DDD << "\t"
        << num_DDL << "\t"
        << num_DLV << "\t"
        << num_DLL << "\t"
        << t_comp_V_DDD << "\t"
        << t_comp_V_DDL << "\t"
        << t_comp_V_DLV << "\t"
        << t_comp_V_DLL << endl;

    fout.close();
}
#endif


void PolygonVD2D::constructVoronoiDiagram( const Polygon2D & polygon )
{
    list<VertexGenerator2D*>    vertexGens;
    list<EdgeGenerator2D*>      edgeGens;
    preprocess_for_construction( polygon, vertexGens, edgeGens );


    constructPhantomVoronoiDiagram();

    
    insertGeneratorsToPhantomVoronoiDiagram( vertexGens, edgeGens );


    removeAllExtraneousVVerticesAndVEdges();
}



void PolygonVD2D::constructVoronoiDiagram( const list<Polygon2D>& polygons )
{
    list<VertexGenerator2D*>    vertexGens;
    list<EdgeGenerator2D*>      edgeGens;
    preprocess_for_construction( polygons, vertexGens, edgeGens );


    constructPhantomVoronoiDiagram();


    insertGeneratorsToPhantomVoronoiDiagram( vertexGens, edgeGens );


    removeAllExtraneousVVerticesAndVEdges();
}



void PolygonVD2D::constructVoronoiDiagramOfRectangle( const rg_Point2D& minPt, const rg_Point2D& maxPt )
{
    b_there_is_container = false;

    Polygon2D polygon;
    polygon.create_polygon_as_rectangle(minPt, maxPt);
    m_polygon.push_back( polygon );
    Polygon2D* ptrRectangle = &m_polygon.back();


    // create polygon generators for rectangle
    double minX = minPt.getX();
    double minY = minPt.getY();
    double maxX = maxPt.getX();
    double maxY = maxPt.getY();

    rg_Point2D          boundaryPt[5] = {rg_Point2D(minX, maxY), rg_Point2D(minX, minY), rg_Point2D(maxX, minY), rg_Point2D(maxX, maxY), rg_Point2D(minX, maxY)};
    VertexGenerator2D*  vertexGenArray[5];
    EdgeGenerator2D*    edgeGenArray[4];

    list<Generator2D*> generatorsOfRectangle;
    for ( int i = 0; i < 4; ++i ) {
        VertexGenerator2D* vertexGen = new VertexGenerator2D(boundaryPt[i], i*2);
        m_generators.push_back(vertexGen);
        vertexGenArray[i] = vertexGen;
        m_mapGenerator2Polygon[vertexGen] = ptrRectangle;
        generatorsOfRectangle.push_back( vertexGen );

        EdgeGenerator2D* edgeGen = new EdgeGenerator2D(boundaryPt[i], boundaryPt[i+1], i*2+1);
        m_generators.push_back(edgeGen);
        edgeGenArray[i] = edgeGen;
        m_mapGenerator2Polygon[edgeGen] = ptrRectangle;
        generatorsOfRectangle.push_back( edgeGen );
    }
    vertexGenArray[4] = vertexGenArray[0];
    m_mapPolygon2Generator[ptrRectangle] = generatorsOfRectangle;


    for ( int i = 0; i < 4; ++i ) {
        vertexGenArray[i]->set_next_edge_generator(edgeGenArray[i]);
        vertexGenArray[i+1]->set_previous_edge_generator(edgeGenArray[i]);

        edgeGenArray[i]->set_start_vertex_generator(vertexGenArray[i]);
        edgeGenArray[i]->set_end_vertex_generator(vertexGenArray[i+1]);
    }

    // create entities
    //  1. make elements of voronoi diagram of phantom generators

    VVertex2D* vertices[14] = { 
        createVertex(0),    // LT
        createVertex(1),    // LB
        createVertex(2),    // RB
        createVertex(3),    // RT
        createVertex(4),    // IN_1
        createVertex(5),    // IN_2
        createVertex(6),    // LTT
        createVertex(7),    // LTL
        createVertex(8),    // LBL
        createVertex(9),    // LBB
        createVertex(10),   // RBB
        createVertex(11),   // RBR
        createVertex(12),   // RTR
        createVertex(13)};  // RTT

    VEdge2D* edges[21] = { 
        createEdge(0),        // INF_LT
        createEdge(1),        // INF_L
        createEdge(2),        // INF_LB
        createEdge(3),        // INF_B
        createEdge(4),        // INF_RB
        createEdge(5),        // INF_R
        createEdge(6),        // INF_RT
        createEdge(7),        // INF_T
        createEdge(8),        // UBND_TL
        createEdge(9),        // UBND_LT
        createEdge(10),       // UBND_LB
        createEdge(11),       // UBND_BL
        createEdge(12),       // UBND_BR
        createEdge(13),       // UBND_RB
        createEdge(14),       // UBND_RT
        createEdge(15),       // UBND_TR
        createEdge(16),       // IN_LT
        createEdge(17),       // IN_LB
        createEdge(18),       // IN_CENTER
        createEdge(19),       // IN_RB
        createEdge(20)};      // IN_RT

    VFace2D* faces[9] = {  
        createFace(-1),       // INF
        createFace(0),        // VTX_LT
        createFace(1),        // EDGE_L
        createFace(2),        // VTX_LB
        createFace(3),        // EDGE_B
        createFace(4),        // VTX_RB
        createFace(5),        // EDGE_R
        createFace(6),        // VTX_RT
        createFace(7)};       // EDGE_T


    // IN_OUT marking.....
    m_VVertexToLocationStatus[vertices[0]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[1]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[2]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[3]]  = ON_POLYGON_BOUNDARY;
    m_VVertexToLocationStatus[vertices[4]]  = INSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[5]]  = INSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[6]]  = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[7]]  = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[8]]  = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[9]]  = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[10]] = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[11]] = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[12]] = OUTSIDE_POLYGON;
    m_VVertexToLocationStatus[vertices[13]] = OUTSIDE_POLYGON;

    m_VEdgeToLocationStatus[edges[0]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[1]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[2]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[3]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[4]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[5]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[6]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[7]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[8]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[9]]       = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[10]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[11]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[12]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[13]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[14]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[15]]      = OUTSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[16]]      = INSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[17]]      = INSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[18]]      = INSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[19]]      = INSIDE_POLYGON;
    m_VEdgeToLocationStatus[edges[20]]      = INSIDE_POLYGON;

    // topology
    vertexGenArray[0]->setOuterFace(faces[1]);
    vertexGenArray[1]->setOuterFace(faces[3]);
    vertexGenArray[2]->setOuterFace(faces[5]);
    vertexGenArray[3]->setOuterFace(faces[7]);

    vertexGenArray[0]->setUserData(vertexGenArray[0]);
    vertexGenArray[1]->setUserData(vertexGenArray[1]);
    vertexGenArray[2]->setUserData(vertexGenArray[2]);
    vertexGenArray[3]->setUserData(vertexGenArray[3]);

    edgeGenArray[0]->setOuterFace(faces[2]);
    edgeGenArray[1]->setOuterFace(faces[4]);
    edgeGenArray[2]->setOuterFace(faces[6]);
    edgeGenArray[3]->setOuterFace(faces[8]);

    edgeGenArray[0]->setUserData(edgeGenArray[0]);
    edgeGenArray[1]->setUserData(edgeGenArray[1]);
    edgeGenArray[2]->setUserData(edgeGenArray[2]);
    edgeGenArray[3]->setUserData(edgeGenArray[3]);

    faces[1]->setGenerator((Generator2D*)vertexGenArray[0]);
    faces[3]->setGenerator((Generator2D*)vertexGenArray[1]);
    faces[5]->setGenerator((Generator2D*)vertexGenArray[2]);
    faces[7]->setGenerator((Generator2D*)vertexGenArray[3]);
    faces[2]->setGenerator((Generator2D*)edgeGenArray[0]);
    faces[4]->setGenerator((Generator2D*)edgeGenArray[1]);
    faces[6]->setGenerator((Generator2D*)edgeGenArray[2]);
    faces[8]->setGenerator((Generator2D*)edgeGenArray[3]);


    vertices[0]->setFirstVEdge(edges[16]);
    vertices[1]->setFirstVEdge(edges[17]);
    vertices[2]->setFirstVEdge(edges[19]);
    vertices[3]->setFirstVEdge(edges[20]);
    vertices[4]->setFirstVEdge(edges[18]);
    vertices[5]->setFirstVEdge(edges[18]);
    vertices[6]->setFirstVEdge(edges[0]);
    vertices[7]->setFirstVEdge(edges[1]);
    vertices[8]->setFirstVEdge(edges[2]);
    vertices[9]->setFirstVEdge(edges[3]);
    vertices[10]->setFirstVEdge(edges[4]);
    vertices[11]->setFirstVEdge(edges[5]);
    vertices[12]->setFirstVEdge(edges[6]);
    vertices[13]->setFirstVEdge(edges[7]);

    edges[0]->setTopology(vertices[6],   vertices[7],  faces[1], faces[0], edges[9], edges[1], edges[8], edges[7]);
    edges[1]->setTopology(vertices[7],   vertices[8],  faces[2], faces[0], edges[10], edges[2], edges[9], edges[0]);
    edges[2]->setTopology(vertices[8],   vertices[9],  faces[3], faces[0], edges[11], edges[3], edges[10], edges[1]);
    edges[3]->setTopology(vertices[9],   vertices[10], faces[4], faces[0], edges[12], edges[4], edges[11], edges[2]);
    edges[4]->setTopology(vertices[10],  vertices[11], faces[5], faces[0], edges[13], edges[5], edges[12], edges[3]);
    edges[5]->setTopology(vertices[11],  vertices[12], faces[6], faces[0], edges[14], edges[6], edges[13], edges[4]);
    edges[6]->setTopology(vertices[12],  vertices[13], faces[7], faces[0], edges[15], edges[7], edges[14], edges[5]);
    edges[7]->setTopology(vertices[13],  vertices[6],  faces[8], faces[0], edges[8], edges[0], edges[15], edges[6]);

    edges[8]->setTopology(vertices[0],   vertices[6],  faces[1], faces[8], edges[0], edges[7], edges[9],  edges[16]);
    edges[9]->setTopology(vertices[0],   vertices[7],  faces[2], faces[1], edges[1], edges[0], edges[16], edges[8]);
    edges[10]->setTopology(vertices[1], vertices[8], faces[3], faces[2],   edges[2], edges[1], edges[11], edges[17]);
    edges[11]->setTopology(vertices[1], vertices[9], faces[4], faces[3],   edges[3], edges[2], edges[17], edges[10]);
    edges[12]->setTopology(vertices[2], vertices[10], faces[5], faces[4],  edges[4], edges[3], edges[13], edges[19]);
    edges[13]->setTopology(vertices[2], vertices[11], faces[6], faces[5],  edges[5], edges[4], edges[19], edges[12]);
    edges[14]->setTopology(vertices[3], vertices[12], faces[7], faces[6],  edges[6], edges[5], edges[15], edges[20]);
    edges[15]->setTopology(vertices[3], vertices[13], faces[8], faces[7],  edges[7], edges[6], edges[20], edges[14]);

    bool b_horizontal_rectangle = true;
    double diffX = maxX - minX;
    double diffY = maxY - minY;
    if ( diffX < diffY )
        b_horizontal_rectangle = false;

    if ( b_horizontal_rectangle ) {
        edges[16]->setTopology(vertices[0], vertices[4], faces[8], faces[2], edges[18], edges[17], edges[8],  edges[9]);
        edges[17]->setTopology(vertices[1], vertices[4], faces[2], faces[4], edges[16], edges[18], edges[10], edges[11]);
        edges[18]->setTopology(vertices[4], vertices[5], faces[8], faces[4], edges[20], edges[19], edges[16], edges[17]);
        edges[19]->setTopology(vertices[5], vertices[2], faces[6], faces[4], edges[13], edges[12], edges[20], edges[18]);
        edges[20]->setTopology(vertices[5], vertices[3], faces[8], faces[6], edges[15], edges[14], edges[18], edges[19]);
    }
    else {
        edges[16]->setTopology(vertices[0], vertices[4], faces[8], faces[2], edges[20], edges[18], edges[8],  edges[9]);
        edges[17]->setTopology(vertices[1], vertices[5], faces[2], faces[4], edges[18], edges[19], edges[10], edges[11]);
        edges[18]->setTopology(vertices[4], vertices[5], faces[6], faces[2], edges[19], edges[17], edges[20], edges[16]);
        edges[19]->setTopology(vertices[5], vertices[2], faces[6], faces[4], edges[13], edges[12], edges[18], edges[17]);
        edges[20]->setTopology(vertices[4], vertices[3], faces[8], faces[6], edges[15], edges[14], edges[16], edges[18]);
    }    

    faces[0]->setFirstVEdge(edges[0]);
    faces[1]->setFirstVEdge(edges[0]);
    faces[2]->setFirstVEdge(edges[1]);
    faces[3]->setFirstVEdge(edges[2]);
    faces[4]->setFirstVEdge(edges[3]);
    faces[5]->setFirstVEdge(edges[4]);
    faces[6]->setFirstVEdge(edges[5]);
    faces[7]->setFirstVEdge(edges[6]);
    faces[8]->setFirstVEdge(edges[7]);


    // geometry

    vertices[0]->setCircumcircle( rg_Circle2D(minX, maxY, 0.0) );
    vertices[1]->setCircumcircle( rg_Circle2D(minX, minY, 0.0) );
    vertices[2]->setCircumcircle( rg_Circle2D(maxX, minY, 0.0) );
    vertices[3]->setCircumcircle( rg_Circle2D(maxX, maxY, 0.0) );

    double radius_circumcircle = 0.0;
    if ( b_horizontal_rectangle ) {
        radius_circumcircle = diffY / 2.0;
        vertices[4]->setCircumcircle( rg_Circle2D(minX+radius_circumcircle, minY+radius_circumcircle, radius_circumcircle) );
        vertices[5]->setCircumcircle( rg_Circle2D(maxX-radius_circumcircle, minY+radius_circumcircle, radius_circumcircle) );

    }
    else {
        radius_circumcircle = diffX / 2.0;
        vertices[4]->setCircumcircle( rg_Circle2D(minX+radius_circumcircle, maxY-radius_circumcircle, radius_circumcircle) );
        vertices[5]->setCircumcircle( rg_Circle2D(minX+radius_circumcircle, minY+radius_circumcircle, radius_circumcircle) );

    }    

    vertices[6]->setCircumcircle( rg_Circle2D(minX, maxY+radius_circumcircle, 0.0) );
    vertices[7]->setCircumcircle( rg_Circle2D(minX-radius_circumcircle, maxY, 0.0) );
    vertices[8]->setCircumcircle( rg_Circle2D(minX-radius_circumcircle, minY, 0.0) );
    vertices[9]->setCircumcircle( rg_Circle2D(minX, minY-radius_circumcircle, 0.0) );
    vertices[10]->setCircumcircle( rg_Circle2D(maxX, minY-radius_circumcircle, 0.0) );
    vertices[11]->setCircumcircle( rg_Circle2D(maxX+radius_circumcircle, minY, 0.0) );
    vertices[12]->setCircumcircle( rg_Circle2D(maxX+radius_circumcircle, maxY, 0.0) );
    vertices[13]->setCircumcircle( rg_Circle2D(maxX, maxY+radius_circumcircle, 0.0) );
}



void PolygonVD2D::preprocess_for_construction( const Polygon2D& polygon, list<VertexGenerator2D*>& vertexGens, list<EdgeGenerator2D*>& edgeGens )
{
    m_polygon.push_back( polygon );
    replicatePolygonEntitiesAsVDGenerators( &m_polygon.back(), vertexGens, edgeGens );
}



void V::GeometryTier::PolygonVD2D::preprocess_for_construction( const list<Polygon2D>& polygons, list<VertexGenerator2D*>& vertexGens, list<EdgeGenerator2D*>& edgeGens )
{
    for ( list<Polygon2D>::const_iterator i_polygon = polygons.begin(); i_polygon != polygons.end(); ++i_polygon ) {
        m_polygon.push_back( *i_polygon );
        
        list<VertexGenerator2D*> currVertexGenerators;
        list<EdgeGenerator2D*>   currEdgeGenerators;
        replicatePolygonEntitiesAsVDGenerators( &m_polygon.back(), currVertexGenerators, currEdgeGenerators );

        vertexGens.insert( vertexGens.end(), currVertexGenerators.begin(), currVertexGenerators.end() );
        edgeGens.insert( edgeGens.end(), currEdgeGenerators.begin(), currEdgeGenerators.end() );
    }
}



/*
This is the block comment for the following function.
To be extracted as a manual description.

This function makes VertexGenerator2D and EdgeGenerator2D objects 
from input polygon (represented as an ordered set of vertices).

INPUT: an ordered/oriented set of vertices.
OUTPUT: to pieces of member data set by this function: VertexGenerator2D and EdgeGenerator2D.
*/
void PolygonVD2D::replicatePolygonEntitiesAsVDGenerators( Polygon2D* const polygon, list<VertexGenerator2D*>& vertexGens, list<EdgeGenerator2D*>& edgeGens )
{    
    list< list< rg_Point2D > > boundaryVertices;
    polygon->get_boundary_vertices( boundaryVertices );


    list<Generator2D*> newGenerators;
    for ( list< list< rg_Point2D > >::iterator i_shell = boundaryVertices.begin(); i_shell != boundaryVertices.end(); ++i_shell ) {
        list< rg_Point2D > shellBoundaryVertices = *i_shell;

        vector<VertexGenerator2D*> vertexGenArray;
        vector<EdgeGenerator2D*>   edgeGenArray;
        vertexGenArray.resize( shellBoundaryVertices.size() + 1 );
        edgeGenArray.resize( shellBoundaryVertices.size() );

        list<rg_Point2D>::iterator i_boundaryVertices = shellBoundaryVertices.begin();
        list<rg_Point2D>::iterator j_boundaryVertices = ++shellBoundaryVertices.begin();

        int index_gen = 0;
        int memberGeneratorUserID = 0;
        if ( m_generators.size() != 0 ) {
            memberGeneratorUserID = m_generators.back()->getID() + 1;
        }

        while ( j_boundaryVertices != shellBoundaryVertices.end() )
        {
            rg_Point2D startPoint = *i_boundaryVertices;
            rg_Point2D endPoint = *j_boundaryVertices;

            VertexGenerator2D* vertexGen = new VertexGenerator2D( startPoint, memberGeneratorUserID++ );
            vertexGen->setDisk( rg_Circle2D( startPoint, 0.0 ) );
            m_generators.push_back( vertexGen );
            vertexGenArray[index_gen] = vertexGen;
            m_mapGenerator2Polygon[vertexGen] = polygon;
            newGenerators.push_back( vertexGen );
            vertexGens.push_back( vertexGen );

            EdgeGenerator2D* edgeGen = new EdgeGenerator2D( startPoint, endPoint, memberGeneratorUserID++ );
            m_generators.push_back( edgeGen );
            edgeGenArray[index_gen] = edgeGen;
            m_mapGenerator2Polygon[edgeGen] = polygon;
            newGenerators.push_back( edgeGen );
            edgeGens.push_back( edgeGen );

            ++i_boundaryVertices;
            ++j_boundaryVertices;
            ++index_gen;
        }

        rg_Point2D startPoint = shellBoundaryVertices.back();
        rg_Point2D endPoint = shellBoundaryVertices.front();

        VertexGenerator2D* vertexGen = new VertexGenerator2D( startPoint, memberGeneratorUserID++ );
        vertexGen->setDisk( rg_Circle2D( startPoint, 0.0 ) );
        m_generators.push_back( vertexGen );
        m_mapGenerator2Polygon[vertexGen] = polygon;
        newGenerators.push_back( vertexGen );
        vertexGenArray[index_gen] = vertexGen;
        vertexGenArray[index_gen + 1] = vertexGenArray[0];
        vertexGens.push_back( vertexGen );

        EdgeGenerator2D* edgeGen = new EdgeGenerator2D( startPoint, endPoint, memberGeneratorUserID );
        m_generators.push_back( edgeGen );
        m_mapGenerator2Polygon[edgeGen] = polygon;
        newGenerators.push_back( edgeGen );
        edgeGenArray[index_gen] = edgeGen;
        edgeGens.push_back( edgeGen );




        // VertexGenerator?? EdgeGenerator ?????? topology?? ?????? ???.
        int numBoundaryVertices = shellBoundaryVertices.size();
        for ( int i = 0; i < numBoundaryVertices; ++i ) {
            vertexGenArray[i]->set_next_edge_generator( edgeGenArray[i] );
            vertexGenArray[i + 1]->set_previous_edge_generator( edgeGenArray[i] );

            edgeGenArray[i]->set_start_vertex_generator( vertexGenArray[i] );
            edgeGenArray[i]->set_end_vertex_generator( vertexGenArray[i + 1] );
        }


    }
    
    m_mapPolygon2Generator[polygon] = newGenerators;
}



void PolygonVD2D::constructPhantomVoronoiDiagram()
{
    createPhantomGenerators();
    constructVoronoiDiagramForPhantomGenerators();
}



void PolygonVD2D::createPhantomGenerators()
{
    rg_BoundingBox2D boundingBox;
    for ( list<Polygon2D>::iterator i_polygon = m_polygon.begin(); i_polygon != m_polygon.end(); ++i_polygon ) {
        Polygon2D* currPolygon = &( *i_polygon );
        rg_BoundingBox2D currBoundingBox;
        currPolygon->get_bounding_box( currBoundingBox );

        boundingBox.contain( currBoundingBox );
    }

    rg_Point2D coordinate[3];
    computeCoordOfPhantomGenerators_beforeOffset( boundingBox, coordinate[0], coordinate[1], coordinate[2] );

    m_phantomCircles.push_back( rg_Circle2D( coordinate[0], 0.0 ) );
    m_phantomCircles.push_back( rg_Circle2D( coordinate[1], 0.0 ) );
    m_phantomCircles.push_back( rg_Circle2D( coordinate[2], 0.0 ) );

    VertexGenerator2D* phantom1 = new VertexGenerator2D( coordinate[0], -3 );
    VertexGenerator2D* phantom2 = new VertexGenerator2D( coordinate[1], -2 );
    VertexGenerator2D* phantom3 = new VertexGenerator2D( coordinate[2], -1 );

    phantom1->setDisk( rg_Circle2D( coordinate[0], 0.0 ) );
    phantom2->setDisk( rg_Circle2D( coordinate[1], 0.0 ) );
    phantom3->setDisk( rg_Circle2D( coordinate[2], 0.0 ) );

    m_phantomGenerators.push_back( phantom1 );
    m_phantomGenerators.push_back( phantom2 );
    m_phantomGenerators.push_back( phantom3 );
}



void PolygonVD2D::insertGeneratorsToPhantomVoronoiDiagram( const list<VertexGenerator2D*>& vertexGenerators, const list<EdgeGenerator2D*>& edgeGenerators )
{
    insertVertexGeneratorSet( vertexGenerators );

    insertEdgeGeneratorSet( edgeGenerators );
}



void PolygonVD2D::collectVertex_N_EdgeGenerators( Polygon2D* const polygon, list<VertexGenerator2D*>& vertexGenerators, list<EdgeGenerator2D*>& edgeGenerators )
{
    list<Generator2D*> gensOfInputPolygon = m_mapPolygon2Generator[polygon];
    for ( list<Generator2D*>::iterator i_gen = gensOfInputPolygon.begin(); i_gen != gensOfInputPolygon.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() == Generator2D::VERTEX_G ) {
            vertexGenerators.push_back( (VertexGenerator2D*)currGen );
        }
        else { // if ( currGen->getType() == Generator2D::EDGE_G ) {
            edgeGenerators.push_back( (EdgeGenerator2D*)currGen );
        }
    }
}



void PolygonVD2D::insertVertexGeneratorSet( const list<VertexGenerator2D*>& vertexGenerators )
{
    int a = 0;
    for ( list<VertexGenerator2D*>::const_iterator i_gen = vertexGenerators.begin(); i_gen != vertexGenerators.end(); ++i_gen ) {
        // ���� �ð��� ��� �̰Ŵ� pointVD�� �ƴϰ� diskVD�Դϴ�.
        // ���߿� �ð����� �����? �ʿ��ϸ� pointVD�� ��ġ����.
        VertexGenerator2D* vertexGen = *i_gen;
        updateVD_with_vertexGenerator( vertexGen );
    }
}



void PolygonVD2D::insertEdgeGeneratorSet( const list<EdgeGenerator2D*>& edgeGenerators )
{
    int a = 0;
    for ( list<EdgeGenerator2D*>::const_iterator i_gen = edgeGenerators.begin(); i_gen != edgeGenerators.end(); ++i_gen ) {
        EdgeGenerator2D* edgeGen = *i_gen;
        updateVD_with_edgeGenerator( edgeGen );
    }
}



void PolygonVD2D::updateVD_with_vertexGenerator( const VertexGenerator2D* const newVertexGen )
{
    updateVoronoiDiagramWithNewGenerator( const_cast<VertexGenerator2D*>(newVertexGen) );
}



void PolygonVD2D::updateVD_with_edgeGenerator( const EdgeGenerator2D* const newEdgeGen )
{
    // The VCell containing the start vertex of the edge is the anchor cell.
    VFace2D* anchorCell = newEdgeGen->get_start_vertex_generator()->getOuterFace();



    // There is always a seed red vertex.
    list<VVertex2D*> redVertices;
    VVertex2D* seedRedVertex = findSeedRedVertex( const_cast<EdgeGenerator2D*>(newEdgeGen), anchorCell );
    seedRedVertex->setStatus( RED_V );
    redVertices.push_back( seedRedVertex );



    // An anomaly edge can be detected in this step.
    // See the comment in mergeSplitVEdgesByFictitiousVVertex.
    list<VVertex2D*> fictitiousVertices;
    list<VVertex2D*> blueVertices;
    wavePropagation_ver1( const_cast<EdgeGenerator2D*>( newEdgeGen ), redVertices, blueVertices, fictitiousVertices );



    // Crossing edge has one red and one blue VVertices.
    list<VEdge2D*> intersectingEdges;
    findCrossingEdges( redVertices, intersectingEdges );



    // This operation must be done before makeVEdgeLoopForNewGeneratorAndConnectToCurrVD
    // even if the opposite seems natural.
    reconfigureByConvertingBlueVVerticesToWhite( blueVertices );



    list<VVertex2D*> newVertices;
    list<VEdge2D*>   newEdges;
    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( const_cast<EdgeGenerator2D*>( newEdgeGen ), intersectingEdges, newVertices, newEdges );



    connectCurrVDToNewVEdgeLoop( newEdges );


    // In WavePropagation step, an anomaly edge could be identified and was split.
    // The split edges should be merged for next iteration.
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVertices );



    computeCoordOfNewVVertices( newVertices );
}



VVertex2D* PolygonVD2D::findSeedRedVertex( Generator2D* const newGen, VFace2D* const anchorCell )
{
    list<VVertex2D*> boundingVertices;
    anchorCell->getBoundaryVVertices( boundingVertices );

    double minMUValue = DBL_MAX;
    VVertex2D* seedRedVertex = NULL;

    for ( list<VVertex2D*>::iterator i_vtx = boundingVertices.begin(); i_vtx != boundingVertices.end(); ++i_vtx ) {
        VVertex2D* currVertex = *i_vtx;

        if ( get_location_status_of_Vvertex( currVertex ) == ON_POLYGON_BOUNDARY ) {
            continue;
        }

        double currMUValue = computeMUValue( currVertex, newGen );
        if ( currMUValue < minMUValue ) {
            minMUValue = currMUValue;
            seedRedVertex = currVertex;
        }
    }

    // When VertexGenerator of new polygon is inserted, self anomaly can exist.
    if ( newGen->getType() != Generator2D::EDGE_G && minMUValue >= 0.0 ) {
        list<VEdge2D*> boundingEdges;
        anchorCell->getBoundaryVEdges( boundingEdges );

        for ( list<VEdge2D*>::iterator i_edge = boundingEdges.begin(); i_edge != boundingEdges.end(); ++i_edge ) {
            VEdge2D* currEdge = *i_edge;
            if ( isAnomalizingEdge( currEdge, newGen ) ) {
                VVertex2D* fictitiousVertex = splitVEdgeAtFictitiousVVertex( currEdge );
                seedRedVertex = fictitiousVertex;
                break;
            }
        }
    }

    return seedRedVertex;
}



double PolygonVD2D::computeMUValue( VVertex2D* const vertex, Generator2D* const newGenerator )
{
    double MUValue = 0.0;
    Generator2D::Generator_Type type = newGenerator->getType();
    switch ( type ) {
    case Generator2D::EDGE_G:
    {
        MUValue = computeMUValueFromEdgeGenerator( vertex, (EdgeGenerator2D*)newGenerator );
    }
    break;

    case Generator2D::VERTEX_G:
    case Generator2D::DISK_G:
    {
        VVertex_Location_Status_Map::iterator i_status = m_VVertexToLocationStatus.find( vertex );
        if ( i_status != m_VVertexToLocationStatus.end() ) {
            if ( i_status->second == VDEntity_Location_Status::ON_POLYGON_BOUNDARY )
                return DBL_MAX;
            else
                return VoronoiDiagram2DC::computeMUValue( vertex, newGenerator );
        }
        else {
            return VoronoiDiagram2DC::computeMUValue( vertex, newGenerator );
        }
    }
    break;

    default:
        break;
    }

    return MUValue;
}



double PolygonVD2D::computeMUValueFromEdgeGenerator( VVertex2D* const vertex, EdgeGenerator2D* const newEdgeGen )
{
    if ( vertex->isInfinite() ) {
        return DBL_MAX;
    }

    rg_Line2D lineSeg = newEdgeGen->get_geometry();
    rg_Point2D projectionPt;
    lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line( vertex->getLocation(), projectionPt );

    double MUValue = 0.0;
    if ( lineSeg.does_contain_in_line_segment( projectionPt ) ) {
        double rho = vertex->getCircumcircle().getRadius();
        double delta = lineSeg.getDistance( vertex->getLocation() );
        MUValue = delta - rho;
    }
    else {
        MUValue = DBL_MAX;
    }

    return MUValue;
}



bool PolygonVD2D::isAnomalizingEdge( VEdge2D* const incidentEdge, Generator2D* const newGenerator )
{
    if ( incidentEdge->getStartVertex()->isFictitious() || incidentEdge->getEndVertex()->isFictitious() ) {
        return false;
    }

    bool b_this_edge_is_an_anomalyEdge = false;

    STATUS_OF_VVERTEX color_startVertex = incidentEdge->getStartVertex()->getStatus();
    STATUS_OF_VVERTEX color_endVertex   = incidentEdge->getEndVertex()->getStatus();

    if ( color_startVertex == WHITE_V ) {
        double MU_startVertex = computeMUValue( incidentEdge->getStartVertex(), newGenerator );

        if ( MU_startVertex < 0.0 ) {
            color_startVertex = RED_V;
        }
        else {
            color_startVertex = BLUE_V;
        }
    }
    if ( color_endVertex == WHITE_V ) {
        double MU_endVertex = computeMUValue( incidentEdge->getEndVertex(), newGenerator );

        if ( MU_endVertex < 0.0 ) {
            color_endVertex = RED_V;
        }
        else {
            color_endVertex = BLUE_V;
        }
    }

    if ( color_startVertex != color_endVertex ) {
        return false;
    }

    if ( incidentEdge->getStartVertex()->isFictitious() || incidentEdge->getEndVertex()->isFictitious() ) {
        return false;
    }



    Generator2D::Generator_Type type = newGenerator->getType();
    switch ( type )
    {
    case Generator2D::DISK_G:
    {
        b_this_edge_is_an_anomalyEdge = doAnomalyEdgeTest_for_DiskGenerator( incidentEdge, newGenerator );
    }
    break;


    case Generator2D::EDGE_G:
    {
        b_this_edge_is_an_anomalyEdge = doAnomalyEdgeTest_for_EdgeGenerator( incidentEdge, newGenerator );
    }
    break;


    case Generator2D::VERTEX_G:
    {
        b_this_edge_is_an_anomalyEdge = doAnomalyEdgeTest_for_VertexGenerator( incidentEdge, newGenerator );
    }
    break;


    default:
        break;
    }

    return b_this_edge_is_an_anomalyEdge;
}



bool PolygonVD2D::doAnomalyEdgeTest_for_DiskGenerator( VEdge2D* const incidentEdge, Generator2D* const newGenerator )
{
    Generator2D* leftGen = (Generator2D*)incidentEdge->getLeftFace()->getGenerator();
    Generator2D* rightGen = (Generator2D*)incidentEdge->getRightFace()->getGenerator();

    Generator2D::Generator_Type leftType = leftGen->getType();
    Generator2D::Generator_Type rightType = rightGen->getType();

    bool this_edge_is_an_anomalyEdge = false;

    switch ( leftType ) {
    case Generator2D::Generator_Type::DISK_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            this_edge_is_an_anomalyEdge = VoronoiDiagram2DC::isAnomalizingEdge( incidentEdge, newGenerator );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;
            double radius_of_disk = diskGen->getDisk().getRadius();

            rg_Point2D focus = diskGen->getDisk().getCenterPt();
            rg_Line2D  directrix = edgeGen->get_geometry();
            rg_Point2D footprint_focus;
            directrix.compute_perpendicular_footprint_of_point_onto_entire_line( focus, footprint_focus );
            rg_Point2D vec_to_directrix = footprint_focus - focus;
            vec_to_directrix = vec_to_directrix.getUnitVector();
            directrix.setSP( directrix.getSP() + vec_to_directrix * radius_of_disk );
            directrix.setEP( directrix.getEP() + vec_to_directrix * radius_of_disk );
            Parabola2D parabola( focus, directrix );

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola( incidentEdge, parabola, newGenerator, radius_of_disk );
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola( incidentEdge, newGenerator->getDisk(), diskGen->getDisk(), rg_Circle2D( vertexGen->get_point(), 0.0 ) );
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::EDGE_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;
            double radius_of_disk = diskGen->getDisk().getRadius();

            rg_Point2D focus = diskGen->getDisk().getCenterPt();
            rg_Line2D  directrix = edgeGen->get_geometry();
            rg_Point2D footprint_focus;
            directrix.compute_perpendicular_footprint_of_point_onto_entire_line( focus, footprint_focus );
            rg_Point2D vec_to_directrix = footprint_focus - focus;
            vec_to_directrix = vec_to_directrix.getUnitVector();
            directrix.setSP( directrix.getSP() + vec_to_directrix * radius_of_disk );
            directrix.setEP( directrix.getEP() + vec_to_directrix * radius_of_disk );
            Parabola2D parabola( focus, directrix );

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola( incidentEdge, parabola, newGenerator, radius_of_disk );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            rg_Circle2D circumcircle_sp = incidentEdge->getStartVertex()->getCircumcircle();
            rg_Circle2D circumcircle_ep = incidentEdge->getEndVertex()->getCircumcircle();

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_line( rg_Line2D( circumcircle_sp.getCenterPt(), circumcircle_ep.getCenterPt() ), circumcircle_sp, circumcircle_ep, newGenerator );
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            if ( edgeGen->get_start_vertex_generator() == vertexGen || edgeGen->get_end_vertex_generator() == vertexGen )
                return false;

            rg_Point2D focus = vertexGen->get_point();
            rg_Line2D  directrix = edgeGen->get_geometry();
            Parabola2D parabola( focus, directrix );

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola( incidentEdge, parabola, newGenerator, 0.0 );
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::VERTEX_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola( incidentEdge, newGenerator->getDisk(), diskGen->getDisk(), rg_Circle2D( vertexGen->get_point(), 0.0 ) );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;

            if ( edgeGen->get_start_vertex_generator() == vertexGen || edgeGen->get_end_vertex_generator() == vertexGen )
                return false;

            rg_Point2D focus = vertexGen->get_point();
            rg_Line2D  directrix = edgeGen->get_geometry();
            Parabola2D parabola( focus, directrix );

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola( incidentEdge, parabola, newGenerator, 0.0 );
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen[2] = { (VertexGenerator2D*)leftGen, (VertexGenerator2D*)rightGen };

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola( incidentEdge, newGenerator->getDisk(), rg_Circle2D( vertexGen[0]->get_point(), 0.0 ), rg_Circle2D( vertexGen[1]->get_point(), 0.0 ) );
        }
        break;

        default:
            break;
        }
    }
    break;

    default:
        break;
    }

    return this_edge_is_an_anomalyEdge;
}



bool PolygonVD2D::doAnomalyEdgeTest_for_EdgeGenerator( VEdge2D* const incidentEdge, Generator2D* const newGenerator )
{
    EdgeGenerator2D*    newEdgeGenerator = (EdgeGenerator2D*)newGenerator;
    rg_Line2D           newLineSeg = newEdgeGenerator->get_geometry();

    Generator2D* leftGen =  (Generator2D*)incidentEdge->getLeftFace()->getGenerator();
    Generator2D* rightGen = (Generator2D*)incidentEdge->getRightFace()->getGenerator();

    Generator2D::Generator_Type leftType = leftGen->getType();
    Generator2D::Generator_Type rightType = rightGen->getType();

    bool this_edge_is_an_anomalyEdge = false;

    switch ( leftType ) {
    case Generator2D::Generator_Type::DISK_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            rg_Circle2D tangentCircle[2];
            int numTangentCircles = 0;
            if ( rg_EQ( leftGen->getDisk().getRadius(), rightGen->getDisk().getRadius() ) ) {
                numTangentCircles = computeTangentCircles_of_two_disks_which_have_same_radii_and_a_line( leftGen->getDisk(), rightGen->getDisk(), newLineSeg, tangentCircle[0], tangentCircle[1] );
            }
            else {
                numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( leftGen->getDisk(), rightGen->getDisk(), newLineSeg, tangentCircle[0], tangentCircle[1] );
            }

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            vector<rg_Circle2D> tangentCircles;
            int numTangentCircles = computeTangentCircles_of_a_disk_and_two_lines( leftGen->getDisk(), ( (EdgeGenerator2D*)rightGen )->get_geometry(), newLineSeg, tangentCircles );

            int numValidCircumcircles = 0;
            rg_Line2D lineSeg[2] = { ( (EdgeGenerator2D*)rightGen )->get_geometry(), newLineSeg };
            switch ( numTangentCircles ) {
            case 2:
            {
                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }

                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }
            }
            break;

            default:
                break;
            }

            if ( numValidCircumcircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircles[0].getCenterPt(), tangentCircles[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            DiskGenerator2D*    diskGen     = (DiskGenerator2D*)leftGen;
            VertexGenerator2D*  vertexGen   = (VertexGenerator2D*)rightGen;

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = 0;
            if ( vertexGen->get_previous_edge_generator() == newGenerator || vertexGen->get_next_edge_generator() == newGenerator ) {
                rg_Point2D  point = vertexGen->get_point();
                rg_Line2D   lineSeg = ((EdgeGenerator2D*)newGenerator)->get_geometry();
                rg_Circle2D disk = diskGen->getDisk();

                rg_Point2D  SP = point;
                rg_Point2D  dirVec = lineSeg.getNormalVector().getUnitVector();

                double distanceBetweenSP_N_disk = SP.distance( disk.getCenterPt() ) - disk.getRadius();
                if ( rg_ZERO( distanceBetweenSP_N_disk ) ) {
                    numTangentCircles = 1;
                }

                double signedDistance = lineSeg.signed_distance( disk.getCenterPt() );

                if ( abs( signedDistance ) < disk.getRadius() ) {
                    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP, dirVec, disk, tangentCircle[0] );
                    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP, -dirVec, disk, tangentCircle[1] );

                    return 2;
                }
                else {
                    numTangentCircles = 1;
                }
            }
            else {
                //rg_Circle2D tangentCircle[2];
                //numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( diskGen->getDisk(), rg_Circle2D( vertexGen->get_point(), 0.0 ), newLineSeg, tangentCircle[0], tangentCircle[1] );
                if ( rg_EQ( leftGen->getDisk().getRadius(), rightGen->getDisk().getRadius() ) ) {
                    numTangentCircles = computeTangentCircles_of_two_disks_which_have_same_radii_and_a_line( leftGen->getDisk(), rightGen->getDisk(), newLineSeg, tangentCircle[0], tangentCircle[1] );
                }
                else {
                    numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( leftGen->getDisk(), rightGen->getDisk(), newLineSeg, tangentCircle[0], tangentCircle[1] );
                }
            }

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::EDGE_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;

            vector<rg_Circle2D> tangentCircles;
            int numTangentCircles = computeTangentCircles_of_a_disk_and_two_lines( diskGen->getDisk(), ( (EdgeGenerator2D*)edgeGen )->get_geometry(), newLineSeg, tangentCircles );

            int numValidCircumcircles = 0;
            rg_Line2D lineSeg[2] = { ( (EdgeGenerator2D*)edgeGen )->get_geometry(), newLineSeg };
            switch ( numTangentCircles ) {
            case 2:
            {
                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }

                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }
            }
            break;

            default:
                break;
            }

            if ( numValidCircumcircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircles[0].getCenterPt(), tangentCircles[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            this_edge_is_an_anomalyEdge = false;
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            vector<rg_Circle2D> tangentCircles;
            int numTangentCircles = 0;

            bool vertexGenIsConnectedWithNewEdgeGenerator = false;
            bool vertexGenIsConnectedWithOldEdgeGenerator = false;
            if ( vertexGen->get_previous_edge_generator() == edgeGen || vertexGen->get_next_edge_generator() == edgeGen ) {
                vertexGenIsConnectedWithOldEdgeGenerator = true;
            }
            if ( vertexGen->get_previous_edge_generator() == newEdgeGenerator || vertexGen->get_next_edge_generator() == newEdgeGenerator ) {
                vertexGenIsConnectedWithNewEdgeGenerator = true;
            }


            rg_Line2D lineSeg[2] = { ( (EdgeGenerator2D*)edgeGen )->get_geometry(), newLineSeg };
            if ( vertexGenIsConnectedWithNewEdgeGenerator && vertexGenIsConnectedWithOldEdgeGenerator ) {
                numTangentCircles = 1;
            }
            else if (!vertexGenIsConnectedWithNewEdgeGenerator && !vertexGenIsConnectedWithOldEdgeGenerator ) {
                numTangentCircles = computeTangentCircles_of_a_disk_and_two_lines( rg_Circle2D( vertexGen->get_point(), 0.0 ), ( (EdgeGenerator2D*)edgeGen )->get_geometry(), newLineSeg, tangentCircles );
            }
            else {
                if ( lineSeg[0].is_parallel_to( lineSeg[1] ) ) {
                    numTangentCircles = 1;
                }
                else {
                    rg_Point2D  SP = vertexGen->get_point();
                    rg_Point2D  dirVec = vertexGenIsConnectedWithOldEdgeGenerator ? lineSeg[0].getNormalVector().getUnitVector() : lineSeg[1].getNormalVector().getUnitVector();
                    rg_Line2D   bisector_of_consecutive_vertex_and_edge( SP, SP + dirVec );

                    rg_Line2D bisector_1 = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1] );
                    rg_Line2D bisector_2 = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1].get_reversed_line2D() );

                    bool dumpFlag = false;
                    rg_Point2D intersectionPt[2];
                    intersectionPt[0] = bisector_1.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );
                    intersectionPt[1] = bisector_2.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );

                    tangentCircles.resize( 2 );
                    numTangentCircles = 2;

                    tangentCircles[0] = rg_Circle2D( intersectionPt[0], 0.0 );
                    tangentCircles[1] = rg_Circle2D( intersectionPt[1], 0.0 );
                }
            }

            int numValidCircumcircles = 0;
            switch ( numTangentCircles ) {
            case 2:
            {
                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }

                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }
            }
            break;

            default:
                break;
            }

            if ( numValidCircumcircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircles[0].getCenterPt(), tangentCircles[1].getCenterPt() );
            }
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::VERTEX_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;

            if ( vertexGen->get_previous_edge_generator() == newGenerator || vertexGen->get_next_edge_generator() == newGenerator ) {
                this_edge_is_an_anomalyEdge = false;
                break;
            }

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( diskGen->getDisk(), rg_Circle2D( vertexGen->get_point(), 0.0 ), newLineSeg, tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;

            if ( vertexGen->get_previous_edge_generator() == edgeGen || vertexGen->get_next_edge_generator() == edgeGen || vertexGen->get_previous_edge_generator() == newEdgeGenerator || vertexGen->get_next_edge_generator() == newEdgeGenerator ) {
                this_edge_is_an_anomalyEdge = false;
                break;
            }

            vector<rg_Circle2D> tangentCircles;
            int numTangentCircles = computeTangentCircles_of_a_disk_and_two_lines( rg_Circle2D( vertexGen->get_point(), 0.0 ), ( (EdgeGenerator2D*)edgeGen )->get_geometry(), newLineSeg, tangentCircles );

            int numValidCircumcircles = 0;
            rg_Line2D lineSeg[2] = { ( (EdgeGenerator2D*)edgeGen )->get_geometry(), newLineSeg };
            switch ( numTangentCircles ) {
            case 2:
            {
                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }

                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }
            }
            break;

            default:
                break;
            }

            if ( numValidCircumcircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircles[0].getCenterPt(), tangentCircles[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen[2] = { (VertexGenerator2D*)leftGen, (VertexGenerator2D*)rightGen };

            if ( vertexGen[0]->get_previous_edge_generator() == newEdgeGenerator || vertexGen[0]->get_next_edge_generator() == newEdgeGenerator || vertexGen[1]->get_previous_edge_generator() == newEdgeGenerator || vertexGen[1]->get_next_edge_generator() == newEdgeGenerator ) {
                this_edge_is_an_anomalyEdge = false;
                break;
            }

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( rg_Circle2D( vertexGen[0]->get_point(), 0.0 ), rg_Circle2D( vertexGen[1]->get_point(), 0.0 ), newLineSeg, tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        default:
            break;
        }
    }
    break;

    default:
        break;
    }

    return this_edge_is_an_anomalyEdge;
}



bool PolygonVD2D::doAnomalyEdgeTest_for_VertexGenerator( VEdge2D* const incidentEdge, Generator2D* const newGenerator )
{
    VertexGenerator2D* newVertexGen = (VertexGenerator2D*)newGenerator;

    Generator2D* leftGen = (Generator2D*)incidentEdge->getLeftFace()->getGenerator();
    Generator2D* rightGen = (Generator2D*)incidentEdge->getRightFace()->getGenerator();

    Generator2D::Generator_Type leftType = leftGen->getType();
    Generator2D::Generator_Type rightType = rightGen->getType();

    bool this_edge_is_an_anomalyEdge = false;

    switch ( leftType ) {
    case Generator2D::Generator_Type::DISK_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola( incidentEdge, leftGen->getDisk(), rightGen->getDisk(), rg_Circle2D( newVertexGen->get_point(), 0.0 ) );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;

            if ( newVertexGen->get_previous_edge_generator() == edgeGen || newVertexGen->get_next_edge_generator() == edgeGen ) {
                this_edge_is_an_anomalyEdge = false;
                break;
            }

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( diskGen->getDisk(), rg_Circle2D( newVertexGen->get_point(), 0.0 ), edgeGen->get_geometry(), tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola( incidentEdge, diskGen->getDisk(), rg_Circle2D( vertexGen->get_point(), 0.0 ), rg_Circle2D( newVertexGen->get_point(), 0.0 ) );
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::EDGE_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;

            if ( newVertexGen->get_previous_edge_generator() == edgeGen || newVertexGen->get_next_edge_generator() == edgeGen ) {
                this_edge_is_an_anomalyEdge = false;
                break;
            }

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( diskGen->getDisk(), rg_Circle2D( newVertexGen->get_point(), 0.0 ), edgeGen->get_geometry(), tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            if ( newVertexGen->get_previous_edge_generator() == leftGen || newVertexGen->get_next_edge_generator() == leftGen || newVertexGen->get_previous_edge_generator() == rightGen || newVertexGen->get_next_edge_generator() == rightGen ) {
                this_edge_is_an_anomalyEdge = false;
                break;
            }

            rg_Line2D leftLine  = ( (EdgeGenerator2D*)leftGen )->get_geometry();
            rg_Line2D rightLine = ( (EdgeGenerator2D*)rightGen )->get_geometry();
            vector<rg_Circle2D> tangentCircles;
            int numTangentCircles = computeTangentCircles_of_a_disk_and_two_lines( rg_Circle2D( newVertexGen->get_point(), 0.0 ), leftLine, rightLine, tangentCircles );

            int numValidCircumcircles = 0;
            rg_Line2D lineSeg[2] = { leftLine, rightLine };
            switch ( numTangentCircles ) {
            case 2:
            {
                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[0].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }

                if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() )
                    && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( tangentCircles[1].getCenterPt() ) )
                {
                    ++numValidCircumcircles;
                }
            }
            break;

            default:
                break;
            }

            if ( numValidCircumcircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircles[0].getCenterPt(), tangentCircles[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            if ( edgeGen->get_start_vertex_generator() == vertexGen || edgeGen->get_end_vertex_generator() == vertexGen || edgeGen->get_start_vertex_generator() == newVertexGen || edgeGen->get_end_vertex_generator() == newVertexGen ) {
                return false;
            }

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( rg_Circle2D( vertexGen->get_point(), 0.0 ), rg_Circle2D( newVertexGen->get_point(), 0.0 ), edgeGen->get_geometry(), tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::VERTEX_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola( incidentEdge, diskGen->getDisk(), rg_Circle2D( vertexGen->get_point(), 0.0 ), rg_Circle2D( newVertexGen->get_point(), 0.0 ) );
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;
            

            if ( edgeGen->get_start_vertex_generator() == vertexGen || edgeGen->get_end_vertex_generator() == vertexGen || edgeGen->get_start_vertex_generator() == newVertexGen || edgeGen->get_end_vertex_generator() == newVertexGen ) {
                return false;
            }

            rg_Circle2D tangentCircle[2];
            int numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( rg_Circle2D( vertexGen->get_point(), 0.0 ), rg_Circle2D( newVertexGen->get_point(), 0.0 ), edgeGen->get_geometry(), tangentCircle[0], tangentCircle[1] );

            if ( numTangentCircles == 2 ) {
                this_edge_is_an_anomalyEdge = areTwoPointsOnThisEdge( incidentEdge, tangentCircle[0].getCenterPt(), tangentCircle[1].getCenterPt() );
            }
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            // this_edge_is_an_anomalyEdge = false;
        }
        break;

        default:
            break;
        }
    }
    break;

    default:
        break;
    }

    return this_edge_is_an_anomalyEdge;
}



bool PolygonVD2D::areTwoPointsOnThisEdge( VEdge2D* edge, const rg_Point2D& pt1, const rg_Point2D& pt2 )
{
    Generator2D* leftGen = edge->getLeftFace()->getGenerator();
    Generator2D::Generator_Type type = leftGen->getType();

    rg_Point2D startPt = edge->getStartVertex()->getLocation();
    rg_Point2D endPt = edge->getEndVertex()->getLocation();

    rg_Point2D focus;

    switch ( type )
    {
    case Generator2D::DISK_G:
    {
        focus = leftGen->getDisk().getCenterPt();
    }
    break;

    case Generator2D::VERTEX_G:
    {
        focus = ( (VertexGenerator2D*)leftGen )->get_point();
    }
    break;

    case Generator2D::EDGE_G:
    {
        rg_Line2D lineSeg = ( (EdgeGenerator2D*)leftGen )->get_geometry();
        rg_Point2D startPt_projected, endPt_projected;
        lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line( startPt, startPt_projected );
        lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line( endPt, endPt_projected );

        rg_Point2D pt1_projected, pt2_projected;
        lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line( pt1, pt1_projected );
        lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line( pt2, pt2_projected );

        bool pt1_is_on_this_edge = false;
        bool pt2_is_on_this_edge = false;

        double param_sp = lineSeg.get_parameter_of_point_on_line_segment( startPt_projected );
        double param_ep = lineSeg.get_parameter_of_point_on_line_segment( endPt_projected );
        double param_pt1 = lineSeg.get_parameter_of_point_on_line_segment( pt1_projected );
        double param_pt2 = lineSeg.get_parameter_of_point_on_line_segment( pt2_projected );

        // This case deals with the case V-edge is defined by two connected VertexGenerator and EdgeGenerator.
        // We assume that the focus is the mid-point of EdgeGenerator and then we do the same process with other cases.
        if ( ( rg_EQ( param_sp, param_pt1 ) && rg_EQ( param_sp, param_pt2 ) ) || ( rg_EQ( param_ep, param_pt1 ) && rg_EQ( param_ep, param_pt2 ) ) ) {
            //focus = ( startPt_projected + endPt_projected ) / 2.0;
            focus = ( lineSeg.getSP() + lineSeg.getEP() ) / 2.0;
        }
        else {
            if ( param_sp > param_ep ) {
                double tempParam = param_sp;
                param_sp = param_ep;
                param_ep = tempParam;
            }
            if ( rg_GE( param_pt1, param_sp ) && rg_LE( param_pt1, param_ep ) && rg_GE( param_pt2, param_sp ) && rg_LE( param_pt2, param_ep ) ) {
                return true;
            }
            else {
                return false;
            }
        }
    }
    break;

    default:
        break;
    }

    rg_Point2D CCW_startVector  = ( startPt - focus ).getUnitVector();
    rg_Point2D CCW_endVector    = ( endPt - focus ).getUnitVector();
    rg_Point2D vec_to_pt1       = ( pt1 - focus ).getUnitVector();
    rg_Point2D vec_to_pt2       = ( pt2 - focus ).getUnitVector();

    double angleToCCW = angleFromVec1toVec2( CCW_startVector, CCW_endVector );
    double angleToPt1 = angleFromVec1toVec2( CCW_startVector, vec_to_pt1 );
    double angleToPt2 = angleFromVec1toVec2( CCW_startVector, vec_to_pt2 );

    if ( rg_EQ( CCW_startVector.getX(), vec_to_pt1.getX() ) && rg_EQ( CCW_startVector.getY(), vec_to_pt1.getY() ) ) {
        angleToPt1 = 0.0;
    }
    if ( rg_EQ( CCW_startVector.getX(), vec_to_pt2.getX() ) && rg_EQ( CCW_startVector.getY(), vec_to_pt2.getY() ) ) {
        angleToPt2 = 0.0;
    }

    if ( rg_LE( angleToPt1, angleToCCW ) && rg_LE( angleToPt2, angleToCCW ) ) {
        return true;
    }
    else {
        return false;
    }
}



void PolygonVD2D::count_number_of_types_of_defining_generators( const VVertex2D* const vertex, int& numDisks, int& numLineSegs, int& numPoints )
{
    numDisks = numLineSegs = numPoints = 0;
    list<Generator2D*> diskGens;
    vertex->getDefining3Generators( diskGens );

    int index_gen = 0;
    for ( list<Generator2D*>::iterator i_gen = diskGens.begin(); i_gen != diskGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        Generator2D::Generator_Type type = currGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::DISK_G:
        {
            ++numDisks;
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            ++numLineSegs;
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            ++numPoints;
        }
        break;

        default:
            break;
        }
    }
}



void PolygonVD2D::computeCoordOfNewVertices( list<VVertex2D*>& newVertices )
{
    for ( list<VVertex2D*>::iterator i_vtx = newVertices.begin(); i_vtx != newVertices.end(); ++i_vtx ) {
        VVertex2D* vertex = *i_vtx;
        computeCoordOf_A_newVertex( vertex );
    }
}



void PolygonVD2D::computeCoordOf_A_newVertex( VVertex2D* const newVertex )
{
#ifdef DEBUG_VERTEX
    ofstream fout_debug_vertex( "vertex_debugging.txt" );
#endif

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
    count_number_of_types_of_defining_generators( newVertex, numDisks, numLines, numPoints );

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
        int numTangentCircles = computeTangentCircles_DDD( newVertex, tangentCircle[0], tangentCircle[1] );
#ifdef CHECK_COMP_TIME
        endTime = clock();
        t_comp_V_DDD = t_comp_V_DDD + endTime - startTime;

        startTime = clock();
#endif
        circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

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
            int numTangentCircles = computeTangentCircles_DDE( newVertex, tangentCircle[0], tangentCircle[1] );

#ifdef CHECK_COMP_TIME
            endTime = clock();
            t_comp_V_DDL = t_comp_V_DDL + endTime - startTime;

            startTime = clock();
#endif
            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );
            
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
            int numTangentCircles = computeTangentCircles_DDV( newVertex, tangentCircle[0], tangentCircle[1] );
#ifdef CHECK_COMP_TIME
            endTime = clock();
            t_comp_V_DDD = t_comp_V_DDD + endTime - startTime;

            startTime = clock();
#endif
            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

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
            int numTangentCircles = computeTangentCircles_DEE( newVertex, tangentCircle[0], tangentCircle[1] );

#ifdef CHECK_COMP_TIME
            endTime = clock();
            t_comp_V_DLL = t_comp_V_DLL + endTime - startTime;

            startTime = clock();
#endif
            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

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
            int numTangentCircles = computeTangentCircles_DEV( newVertex, tangentCircle[0], tangentCircle[1] );

            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

        }
        break;

        case 0: // numDisks == 1, numLines == 0, numPoints == 2
        {
#ifdef CHECK_COMP_TIME
            ++num_DDD;
            startTime = clock();
#endif
            int numTangentCircles = computeTangentCircles_DVV( newVertex, tangentCircle[0], tangentCircle[1] );

#ifdef CHECK_COMP_TIME
            endTime = clock();
            t_comp_V_DDD = t_comp_V_DDD + endTime - startTime;

            startTime = clock();
#endif
            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

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



    case 0:
    {
        switch ( numLines )
        {
        case 3: // numDisks == 0, numLines == 3, numPoints == 0
        {
            int numTangentCircles = computeTangentCircles_EEE( newVertex, tangentCircle[0], tangentCircle[1] );

            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

        }
        break;

        case 2: // numDisks == 0, numLines == 2, numPoints == 1
        {
            int numTangentCircles = computeTangentCircles_EEV( newVertex, tangentCircle[0], tangentCircle[1] );

            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );

        }
        break;

        case 1: // numDisks == 0, numLines == 1, numPoints == 2
        {
            int numTangentCircles = computeTangentCircles_EVV( newVertex, tangentCircle[0], tangentCircle[1] );

            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );
        }
        break;

        case 0: // numDisks == 0, numLines == 0, numPoints == 3
        {
            int numTangentCircles = computeTangentCircles_VVV( newVertex, tangentCircle[0], tangentCircle[1] );

            circumcircle_of_vertex = chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( numTangentCircles, tangentCircle[0], tangentCircle[1], newVertex );
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
    if ( rg_ZERO( circumcircle_of_vertex.getCenterPt().getX() ) && rg_ZERO( circumcircle_of_vertex.getCenterPt().getY() ) ) {
        int stopHere = 1;
    }

    newVertex->setCircumcircle( circumcircle_of_vertex );
    newVertex->setStatus( WHITE_V );

#ifdef DEBUG_VERTEX
    fout_debug_vertex.close();
#endif
}



int PolygonVD2D::computeTangentCircles_DDD( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    rg_Circle2D disks[3];
    int index = 0;
    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen, ++index ) {
        Generator2D* currGen = *i_gen;
        if ( currGen->getType() == Generator2D::DISK_G ) {
            disks[index] = currGen->getDisk();
        }
        else { // currGen->getType() == Generator2D::VERTEX_G
            rg_Point2D centerPt = ( (VertexGenerator2D*)currGen )->get_point();
            disks[index] = rg_Circle2D( centerPt, 0.0 );
        }
    }

    return rg_Circle2D::makeCircumcircle( disks[0], disks[1], disks[2], tangentCircle1, tangentCircle2 );
}



int PolygonVD2D::computeTangentCircles_DDE( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    rg_Circle2D disks[2];
    rg_Line2D   lineSeg;
    rg_INT      diskIndex = 0;

    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        if ( Generator2D::DISK_G == currGen->getType() ) {
            disks[diskIndex++] = ( (DiskGenerator2D*)currGen )->getDisk();
        }
        else if ( Generator2D::VERTEX_G == currGen->getType() ) {
            disks[diskIndex++] = rg_Circle2D( ( (VertexGenerator2D*)currGen )->get_point(), 0.0 );
        }
        else { // ( Generator2D::EDGE_G == currGen->getType() ) {
            ( (EdgeGenerator2D*)currGen )->get_geometry( lineSeg );
        }
    }

    //rg_Circle2D candidateTangentCircle[2];
    vector<rg_Circle2D> candidateTangentCircle;
    int numCircumcircles = 0;
    if ( rg_EQ( disks[0].getRadius(), disks[1].getRadius() ) ) {
        rg_Circle2D candidateTangentCircle_array[2];
        numCircumcircles = computeTangentCircles_of_two_disks_which_have_same_radii_and_a_line( disks[0], disks[1], lineSeg, candidateTangentCircle_array[0], candidateTangentCircle_array[1] );

        for ( int i = 0; i < numCircumcircles; ++i ) {
            candidateTangentCircle.push_back( candidateTangentCircle_array[i] );
        }
    }
    else {
        //vector<rg_Circle2D> candidateTangentCircles_vector;
        //numCircumcircles = computeTangentCircles_of_two_disks_and_a_line( disks[0], disks[1], lineSeg, candidateTangentCircle[0], candidateTangentCircle[1] );
        numCircumcircles = computeTangentCircles_of_two_disks_and_a_line_for_coordinate( disks[0], disks[1], lineSeg, candidateTangentCircle );
        
        //if ( numCircumcircles > 0 ) {
        //    candidateTangentCircle[0] = candidateTangentCircles_vector[0];
        //}
        //if ( numCircumcircles > 1 ) {
        //    candidateTangentCircle[1] = candidateTangentCircles_vector[1];
        //}
    }


    list<rg_Circle2D> validTangentCircles;
    for ( int i = 0; i < numCircumcircles; ++i ) {
        if ( lineSeg.Is_perpendicular_footprint_of_point_on_line_segment( candidateTangentCircle[i].getCenterPt() ) )
        {
            validTangentCircles.push_back( candidateTangentCircle[i] );
        }
    }

    //switch ( numCircumcircles ) {
    //case 2:
    //{
    //    if ( lineSeg.Is_perpendicular_footprint_of_point_on_line_segment( candidateTangentCircle[0].getCenterPt() ) )
    //    {
    //        validTangentCircles.push_back( candidateTangentCircle[0] );
    //    }
    //    if ( lineSeg.Is_perpendicular_footprint_of_point_on_line_segment( candidateTangentCircle[1].getCenterPt() ) )
    //    {
    //        validTangentCircles.push_back( candidateTangentCircle[1] );
    //    }
    //}
    //break;
    //
    //case 1:
    //{
    //    if ( lineSeg.Is_perpendicular_footprint_of_point_on_line_segment( candidateTangentCircle[0].getCenterPt() ) )
    //    {
    //        validTangentCircles.push_back( candidateTangentCircle[0] );
    //    }
    //}
    //break;
    //
    //default:
    //    break;
    //}


    switch ( validTangentCircles.size() ) {
    case 2:
    {
        tangentCircle1 = validTangentCircles.front();
        tangentCircle2 = validTangentCircles.back();
    }
    break;

    case 1:
    {
        tangentCircle1 = validTangentCircles.front();
    }
    break;

    default:
        break;
    }

    return validTangentCircles.size();
}



int PolygonVD2D::computeTangentCircles_DEE( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    rg_Circle2D disk;
    rg_Line2D   lineSeg[2];
    int         lineSegIndex = 0;

    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        if ( Generator2D::DISK_G == currGen->getType() ) {
            disk = ( (DiskGenerator2D*)currGen )->getDisk();
        }
        else if ( Generator2D::VERTEX_G == currGen->getType() ) {
            disk = rg_Circle2D( ( (VertexGenerator2D*)currGen )->get_point(), 0.0 );
        }
        else { // ( Generator2D::EDGE_G == currGen->getType() ) {
            lineSeg[lineSegIndex++] = ( (EdgeGenerator2D*)currGen )->get_geometry();
        }
    }

    vector<rg_Circle2D> candidateTangentCircles;
    int numCircumcircles = computeTangentCircles_of_a_disk_and_two_lines( disk, lineSeg[0], lineSeg[1], candidateTangentCircles );

    list<rg_Circle2D> validTangentCircles;
    for ( int i = 0; i < numCircumcircles; ++i ) {
        if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( candidateTangentCircles[i].getCenterPt() )
            && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( candidateTangentCircles[i].getCenterPt() ) )
        {
            validTangentCircles.push_back( candidateTangentCircles[i] );
        }
    }

    switch ( validTangentCircles.size() ) {
    case 2:
    {
        tangentCircle1 = validTangentCircles.front();
        tangentCircle2 = validTangentCircles.back();
    }
    break;

    case 1:
    {
        tangentCircle1 = validTangentCircles.front();
    }
    break;

    case 0:
    {
        // SPECIAL CASE
        // We assume that case 0 is occured by numerical error so circumcircle which has closest footprint to both line segments.
        rg_Circle2D circumcircle_which_has_closest_footprint;
        double minDistanceBetweenLine_N_footprint = DBL_MAX;

        for ( int i = 0; i < numCircumcircles; ++i ) {
            rg_Point2D footprint_line0;
            lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( candidateTangentCircles[i].getCenterPt(), footprint_line0 );
            double dist_line0 = 0.0;
            if ( !lineSeg[0].does_contain_in_line_segment( footprint_line0 ) ) {
                double dist_SP_0 = lineSeg[0].getSP().distance( footprint_line0 );
                double dist_EP_0 = lineSeg[0].getEP().distance( footprint_line0 );
                dist_line0 = min( dist_SP_0, dist_EP_0 );
            }

            rg_Point2D footprint_line1;
            lineSeg[1].compute_perpendicular_footprint_of_point_onto_entire_line( candidateTangentCircles[i].getCenterPt(), footprint_line1 );
            double dist_line1 = 0.0;
            if ( !lineSeg[1].does_contain_in_line_segment( footprint_line1 ) ) {
                double dist_SP_1 = lineSeg[1].getSP().distance( footprint_line1 );
                double dist_EP_1 = lineSeg[1].getEP().distance( footprint_line1 );
                dist_line1 = min( dist_SP_1, dist_EP_1 );
            }

            double currDist = max( dist_line0, dist_line1 );
            if ( currDist < minDistanceBetweenLine_N_footprint ) {
                minDistanceBetweenLine_N_footprint = currDist;
                circumcircle_which_has_closest_footprint = candidateTangentCircles[i];
            }
        }

        //tangentCircle1 = circumcircle_which_has_closest_footprint;
        tangentCircle1 = candidateTangentCircles[0];
        return 1;
    }
    default:
        break;
    }

    return validTangentCircles.size();
}



int PolygonVD2D::computeTangentCircles_DEV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    DiskGenerator2D*    diskGen     = NULL;
    EdgeGenerator2D*    edgeGen     = NULL;
    VertexGenerator2D*  vertexGen   = NULL;

    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        if ( Generator2D::DISK_G == currGen->getType() ) {
            diskGen = (DiskGenerator2D*)currGen;
        }
        else if ( Generator2D::VERTEX_G == currGen->getType() ) {
            vertexGen = (VertexGenerator2D*)currGen;
        }
        else { // ( Generator2D::EDGE_G == currGen->getType() ) {
            edgeGen = (EdgeGenerator2D*)currGen;
        }
    }
#ifdef CHECK_COMP_TIME
    endTime = clock();
    double localCompTime = endTime - startTime;
#endif

    if ( vertexGen->get_next_edge_generator() != edgeGen && vertexGen->get_previous_edge_generator() != edgeGen ) {
        return computeTangentCircles_DDE( newVertex, tangentCircle1, tangentCircle2 );
    }

    rg_Point2D  point   = vertexGen->get_point();
    rg_Line2D   lineSeg = edgeGen->get_geometry();
    rg_Circle2D disk    = diskGen->getDisk();

    rg_Point2D  SP = point;
    rg_Point2D  dirVec = lineSeg.getNormalVector().getUnitVector();

    double distanceBetweenSP_N_disk = SP.distance( disk.getCenterPt() ) - disk.getRadius();
    if ( rg_ZERO( distanceBetweenSP_N_disk ) ) {
        tangentCircle1.setCenterPt( SP );
        tangentCircle1.setRadius( 0.0 );
        return 1;
    }

    double signedDistance = lineSeg.signed_distance( disk.getCenterPt() );

    if ( abs( signedDistance ) < disk.getRadius() ) {
        computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP,  dirVec, disk, tangentCircle1 );
        computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP, -dirVec, disk, tangentCircle2 );

        return 2;
    }
    else {
        if ( signedDistance < 0.0 ) {
            dirVec = -dirVec;
        }

        computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( SP, dirVec, disk, tangentCircle1 );

        return 1;
    }
}



int PolygonVD2D::computeTangentCircles_EEE( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    rg_Line2D lineSeg[3];
    rg_INT lineSegIndex = 0;

    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        lineSeg[lineSegIndex++] = ( (EdgeGenerator2D*)currGen )->get_geometry();
    }



    rg_Line2D baseLine  = lineSeg[0];
    rg_Line2D line1     = lineSeg[1];
    rg_Line2D line2     = lineSeg[2];

    // Check if there is a pair of parallel line segments.
    bool there_is_parallel_line_pair = false;
    if ( lineSeg[0].is_parallel_to( lineSeg[1] ) ) {
        baseLine = lineSeg[2];
        line1 = lineSeg[0];
        line2 = lineSeg[1];
    }
    else if ( lineSeg[1].is_parallel_to( lineSeg[2] ) ) {
        /*
        baseLine = lineSeg[0];
        line1 = lineSeg[1];
        line2 = lineSeg[2];
        */
    }
    else if ( lineSeg[2].is_parallel_to( lineSeg[0] ) ) {
        baseLine = lineSeg[1];
        line1 = lineSeg[0];
        line2 = lineSeg[2];
    }
    
    rg_Line2D bisector_1  = rg_GeoFunc::compute_bisector_line_between_two_line_segments( baseLine, line1 );
    rg_Line2D bisector_1_ = rg_GeoFunc::compute_bisector_line_between_two_line_segments( baseLine, line1.get_reversed_line2D() );
    rg_Line2D bisector_2  = rg_GeoFunc::compute_bisector_line_between_two_line_segments( baseLine, line2 );
    rg_Line2D bisector_2_ = rg_GeoFunc::compute_bisector_line_between_two_line_segments( baseLine, line2.get_reversed_line2D() );



    // Compute intersection pts ( 1 vs 2, 1 vs 2_, 1_ vs 2, 1_ vs 2_ )
    bool bisectors_are_parallel[4] = { false, false, false, false };
    rg_Point2D intersectionPt[4];

    intersectionPt[0] = bisector_1.compute_intersection_with_line( bisector_2, bisectors_are_parallel[0] );
    intersectionPt[1] = bisector_1.compute_intersection_with_line( bisector_2_, bisectors_are_parallel[1] );
    intersectionPt[2] = bisector_1_.compute_intersection_with_line( bisector_2, bisectors_are_parallel[2] );
    intersectionPt[3] = bisector_1_.compute_intersection_with_line( bisector_2_, bisectors_are_parallel[3] );




    list<rg_Point2D> validPts;
    // Check validity
    for ( int i = 0; i < 4; ++i ) {
        if ( bisectors_are_parallel[i] ) {
            continue;
        }

        if ( lineSeg[0].Is_perpendicular_footprint_of_point_on_line_segment( intersectionPt[i] )
            && lineSeg[1].Is_perpendicular_footprint_of_point_on_line_segment( intersectionPt[i] )
            && lineSeg[2].Is_perpendicular_footprint_of_point_on_line_segment( intersectionPt[i] ) ) 
        {
            validPts.push_back( intersectionPt[i] );
        }
    }

    switch ( validPts.size() )
    {
    case 2:
    {
        double clearance1 = lineSeg[0].getDistance ( validPts.front() );
        tangentCircle1.setCenterPt( validPts.front() );
        tangentCircle1.setRadius( clearance1 );

        double clearance2 = lineSeg[1].getDistance ( validPts.back() );
        tangentCircle2.setCenterPt( validPts.back() );
        tangentCircle2.setRadius( clearance2 );
    }
    break;

    case 1:
    {
        double clearance1 = lineSeg[0].getDistance ( validPts.front() );
        tangentCircle1.setCenterPt( validPts.front() );
        tangentCircle1.setRadius( clearance1 );
    }

    default:
        break;
    }



    //rg_Line2D bisector_of_line0_N_line1 = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1].get_reversed_line2D() );
    //rg_Line2D bisector_of_line0_N_line2 = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[2].get_reversed_line2D() );

    //bool b_two_bisector_are_parallel = false;
    //rg_Point2D intersectionPt = bisector_of_line0_N_line1.compute_intersection_with_line( bisector_of_line0_N_line2, b_two_bisector_are_parallel );
    //double clearance = lineSeg[0].getDistance( intersectionPt );

    //tangentCircle1.setCenterPt( intersectionPt );
    //tangentCircle1.setRadius( clearance );

    return validPts.size();
}



int PolygonVD2D::computeTangentCircles_EEV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    // Get generators by their type.
    EdgeGenerator2D* edgeGen[2] = { NULL, NULL };
    VertexGenerator2D*  vertexGen = NULL;
    int                 lineSegIndex = 0;

    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        if ( Generator2D::VERTEX_G == currGen->getType() ) {
            vertexGen = (VertexGenerator2D*)currGen;
        }
        else { // ( Generator2D::EDGE_G == currGen->getType() ) {
            edgeGen[lineSegIndex++] = (EdgeGenerator2D*)currGen;
        }
    }

    // Classify LLP case into three configurations according to their adjacency.
    bool b_vertex_is_incident_to_lineSeg_0 = false;
    bool b_vertex_is_incident_to_lineSeg_1 = false;
    if ( edgeGen[0]->get_start_vertex_generator() == vertexGen || edgeGen[0]->get_end_vertex_generator() == vertexGen ) {
        b_vertex_is_incident_to_lineSeg_0 = true;
    }
    if ( edgeGen[1]->get_start_vertex_generator() == vertexGen || edgeGen[1]->get_end_vertex_generator() == vertexGen ) {
        b_vertex_is_incident_to_lineSeg_1 = true;
    }


    rg_Point2D  point = vertexGen->get_point();
    rg_Line2D   lineSeg[2] = { edgeGen[0]->get_geometry(), edgeGen[1]->get_geometry().get_reversed_line2D() };

    // CASE I: All three are disconnected to each other.
    if ( !b_vertex_is_incident_to_lineSeg_0 && !b_vertex_is_incident_to_lineSeg_1 ) {
        return computeTangentCircles_DEE( newVertex, tangentCircle1, tangentCircle2 );
    }
    // CASE II: All three are connected.
    else if ( b_vertex_is_incident_to_lineSeg_0 && b_vertex_is_incident_to_lineSeg_1 ) {
        tangentCircle1.setCenterPt( vertexGen->get_point() );
        tangentCircle1.setRadius( 0.0 );
        return 1;
    }

    // CASE III: The VertexGenerator is connected with one EdgeGenerator 
    //           but is disconnected with the other.
    rg_Point2D  SP = point;
    rg_Point2D  dirVec = b_vertex_is_incident_to_lineSeg_0 ? lineSeg[0].getNormalVector().getUnitVector() : lineSeg[1].getNormalVector().getUnitVector();
    rg_Line2D   bisector_of_consecutive_vertex_and_edge( SP, SP + dirVec );
    
    rg_Line2D   non_adjacent_lineSeg = b_vertex_is_incident_to_lineSeg_0 ? lineSeg[1] : lineSeg[0];

    list<rg_Point2D> validPts;
    if ( lineSeg[0].is_parallel_to( lineSeg[1] ) ) {
        rg_Line2D   bisector_of_two_lineSegs = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1] );
        
        bool        dumpFlag = false;
        rg_Point2D  intersectionPt = bisector_of_two_lineSegs.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );

        
        rg_Point2D  projectionPt;
        non_adjacent_lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line( intersectionPt, projectionPt );

        if ( !non_adjacent_lineSeg.does_contain_in_line_segment( projectionPt ) ) {
            bisector_of_two_lineSegs = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1].get_reversed_line2D() );
            intersectionPt = bisector_of_two_lineSegs.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );
        }

        validPts.push_back( intersectionPt );
    }
    else {
        rg_Line2D bisector_1 = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1] );
        rg_Line2D bisector_2 = rg_GeoFunc::compute_bisector_line_between_two_line_segments( lineSeg[0], lineSeg[1].get_reversed_line2D() );

        bool dumpFlag = false;
        rg_Point2D intersectionPt[2];
        intersectionPt[0] = bisector_1.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );
        intersectionPt[1] = bisector_2.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );

        if ( non_adjacent_lineSeg.Is_perpendicular_footprint_of_point_on_line_segment( intersectionPt[0] ) ) {
            validPts.push_back( intersectionPt[0] );
        }
        if ( non_adjacent_lineSeg.Is_perpendicular_footprint_of_point_on_line_segment( intersectionPt[1] ) ) {
            validPts.push_back( intersectionPt[1] );
        }
    }


    switch ( validPts.size() )
    {
    case 2:
    {
        double clearance = lineSeg[0].getDistance( validPts.front() );
        tangentCircle1.setCenterPt( validPts.front() );
        tangentCircle1.setRadius( clearance );

        clearance = lineSeg[0].getDistance( validPts.back() );
        tangentCircle2.setCenterPt( validPts.back() );
        tangentCircle2.setRadius( clearance );
    }
    break;

    case 1:
    {
        double clearance = lineSeg[0].getDistance( validPts.front() );
        tangentCircle1.setCenterPt( validPts.front() );
        tangentCircle1.setRadius( clearance );
    }
    break;

    default:
        break;
    }

    return validPts.size();
}



int PolygonVD2D::computeTangentCircles_EVV( VVertex2D* const newVertex, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    list<Generator2D*> definingGens;
    newVertex->getDefining3Generators( definingGens );

    EdgeGenerator2D*    edgeGen = NULL;
    VertexGenerator2D* vertexGen[2] = { NULL, NULL };
    int                 vertexIndex = 0;

    for ( list<Generator2D*>::iterator i_gen = definingGens.begin(); i_gen != definingGens.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;
        if ( Generator2D::VERTEX_G == currGen->getType() ) {
            vertexGen[vertexIndex++] = (VertexGenerator2D*)currGen;
        }
        else { // ( Generator2D::EDGE_G == currGen->getType() ) {
            edgeGen = (EdgeGenerator2D*)currGen;
        }
    }

    bool b_edge_is_incident_to_vertex_0 = false;
    bool b_edge_is_incident_to_vertex_1 = false;
    if ( edgeGen->get_start_vertex_generator() == vertexGen[0] || edgeGen->get_end_vertex_generator() == vertexGen[0] ) {
        b_edge_is_incident_to_vertex_0 = true;
    }
    if ( edgeGen->get_start_vertex_generator() == vertexGen[1] || edgeGen->get_end_vertex_generator() == vertexGen[1] ) {
        b_edge_is_incident_to_vertex_1 = true;
    }



    if ( !b_edge_is_incident_to_vertex_0 && !b_edge_is_incident_to_vertex_1 ) {
        return computeTangentCircles_DDE( newVertex, tangentCircle1, tangentCircle2 );
    }
    else if ( b_edge_is_incident_to_vertex_0 && b_edge_is_incident_to_vertex_1 ) {
        // impossible
        // WARNING!!!
        return 0;
    }

    rg_Point2D  point[2] = { vertexGen[0]->get_point(), vertexGen[1]->get_point() };
    rg_Line2D   lineSeg = edgeGen->get_geometry();

    rg_Line2D   bisector_of_two_points = rg_GeoFunc::compute_bisector_line_between_two_points( point[0], point[1] );

    rg_Point2D  SP      = b_edge_is_incident_to_vertex_0 ? point[0] : point[1];
    rg_Point2D  dirVec  = lineSeg.getNormalVector().getUnitVector();
    rg_Line2D   bisector_of_consecutive_vertex_and_edge( SP, SP + dirVec );

    bool        dumpFlag = false;
    rg_Point2D  intersectionPt = bisector_of_two_points.compute_intersection_with_line( bisector_of_consecutive_vertex_and_edge, dumpFlag );
    double      clearance = lineSeg.getDistance( intersectionPt );

    tangentCircle1.setCenterPt( intersectionPt );
    tangentCircle1.setRadius( clearance );

    return 1;
}



int PolygonVD2D::computeTangentCircles_of_two_disks_which_have_same_radii_and_a_line( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    double      radius_of_disk = disk1.getRadius();

    double signedDistance_disk1_center = line.signed_distance( disk1.getCenterPt() );
    double signedDistance_disk2_center = line.signed_distance( disk2.getCenterPt() );

    bool b_disk1_is_intersected_with_line = rg_LE( abs( signedDistance_disk1_center ), radius_of_disk ) ? true : false;
    bool b_disk2_is_intersected_with_line = rg_LE( abs( signedDistance_disk2_center ), radius_of_disk ) ? true : false;

    if ( !b_disk1_is_intersected_with_line && !b_disk2_is_intersected_with_line && ( signedDistance_disk1_center * signedDistance_disk2_center <= 0.0 ) ) {
        return 0;
    }

    rg_Point2D  dirVec_bisector = ( disk2.getCenterPt() - disk1.getCenterPt() ).getPerpendicularVector().getUnitVector();
    rg_Point2D  SP              = ( disk1.getCenterPt() + disk2.getCenterPt() ) / 2.0;
    rg_Point2D  EP              = SP + dirVec_bisector;

    // if SP and EP are in different half-space of input line, we use opposite direction vector of bisector to get another EP in same half-space.
    double signedDistance_SP = line.signed_distance( SP );
    double signedDistance_EP = line.signed_distance( EP );
    if ( signedDistance_SP * signedDistance_EP < 0.0 ) {
        dirVec_bisector = -dirVec_bisector;
        EP = SP + dirVec_bisector;
    }
    // We assume input disks have zero-radii.
    // For this assumption, we add that radius to clearance.
    double      clearance_SP    = abs( signedDistance_SP ) + radius_of_disk;
    double      clearance_EP    = abs( line.signed_distance( EP ) ) + radius_of_disk;
    double      diff_clearance  = clearance_EP - clearance_SP;

    // We assume that we transform input data to make SP is on origin point and EP has ( 1.0, 0.0 ) coordinate.
    // Then, the center if disk1 has ( 0.0, y_1 ) coodinate, where y is distance between the center of disk1 and SP.
    // However, the clearnace of SP and EP are not affected by transformation.
    double      y_1             = disk1.getCenterPt().distance( SP );

    // A point on the bisector has ( t, 0 ) coordinate.
    // The problem becomes Finding T which has distance D and clearance C and D == C.
    // D = sqrt ( t*t + y_1*y_1 ), and C = clearance_SP + t*diff_clearance.
    // ( -half_b +/- sqrt( half_b^2 - a*c ) ) / a
    double      a               = diff_clearance * diff_clearance - 1.0;
    double      half_b          = diff_clearance * clearance_SP;
    double      c               = ( clearance_SP * clearance_SP ) - ( y_1 * y_1 );
    double      bb_ac           = half_b * half_b - a * c;

    if ( bb_ac < 0.0 ) {
        return 0;
    }

    if ( rg_ZERO( a ) ) {
        double t = -c / half_b / 2.0;
        tangentCircle1.setCenterPt( SP + t * dirVec_bisector );
        tangentCircle1.setRadius( clearance_SP + t * diff_clearance - radius_of_disk );

        return 1;
    }

    double      squareRoot      = sqrt( bb_ac );
    if ( rg_ZERO( squareRoot ) ) {
        double  t = -half_b / a;
        tangentCircle1.setCenterPt( SP + t * dirVec_bisector );
        tangentCircle1.setRadius( clearance_SP + t * diff_clearance - radius_of_disk );

        return 1;
    }

    double t_1 = ( -half_b + squareRoot ) / a;
    double t_2 = ( -half_b - squareRoot ) / a;
    tangentCircle1.setCenterPt( SP + t_1 * dirVec_bisector );
    tangentCircle1.setRadius( clearance_SP + t_1 * diff_clearance - radius_of_disk );
    tangentCircle2.setCenterPt( SP + t_2 * dirVec_bisector );
    tangentCircle2.setRadius( clearance_SP + t_2 * diff_clearance - radius_of_disk );

    return 2;
}



int PolygonVD2D::computeTangentCircles_of_two_disks_and_a_line( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    //terminate_for_debugging = false;
    // Previous numerial computation
    /*
    rg_Circle2D disk_base;
    rg_Circle2D another_disk;

    {
        double distanceBetweenDisk1Center_N_line = line.getDistance(disk1.getCenterPt());
        double distanceBetweenDisk2Center_N_line = line.getDistance(disk2.getCenterPt());

        if ( distanceBetweenDisk1Center_N_line + disk1.getRadius() < distanceBetweenDisk2Center_N_line + disk2.getRadius() ) {
            disk_base = disk1;
            another_disk = disk2;
        }
        else {
            disk_base = disk2;
            another_disk = disk1;
        }
    }

    rg_Point2D dirVec = line.getNormalVector().getUnitVector();
    rg_Point2D footprint_d;
    line.compute_perpendicular_footprint_of_point_onto_entire_line(disk_base.getCenterPt(), footprint_d);

    rg_Circle2D tangentCircle_base;
    rg_Circle2D tangentCircle_another;
    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( footprint_d, dirVec, disk_base,    tangentCircle_base);
    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( footprint_d, dirVec, another_disk, tangentCircle_another);

    //if ( rg_GE(tangentCircle_base.getRadius(), tangentCircle_another.getRadius()) ) {
    //    return 0;
    //}

    rg_Point2D  midPt       = (line.getSP() + line.getEP()) / 2.0;
    double      halfLength  = line.getLength() / 2.0;
    double      distanceBetween_midPT_N_footprint_d = footprint_d.distance(midPt);

    bool        b_footprint_is_before_line  = false;
    bool        b_footprint_is_on_line      = false;
    bool        b_footprint_is_after_line   = false;

    if ( distanceBetween_midPT_N_footprint_d < halfLength ) {
        b_footprint_is_on_line = true;
    }
    else {
        double diff, diff_f;
        if ( rg_EQ( line.getSP().getX(), line.getEP().getX() ) ) {
            diff    = line.getEP().getY()  - line.getSP().getY();
            diff_f  = footprint_d.getY() - line.getSP().getY();
        }
        else {
            diff    = line.getEP().getX()  - line.getSP().getX();
            diff_f  = footprint_d.getX() - line.getSP().getX();
        }

        if ( signbit(diff) == signbit(diff_f) )
            b_footprint_is_after_line;
        else
            b_footprint_is_before_line;
    }

    int numTangentCircles = 0;
    if ( b_footprint_is_on_line ) {
        bool there_is_intersection1 = findIntersectionPtWithinThisArrange(line, disk_base, another_disk, line.getSP(), footprint_d, tangentCircle1);
        bool there_is_intersection2 = findIntersectionPtWithinThisArrange(line, disk_base, another_disk, footprint_d, line.getEP(), tangentCircle2);

        if ( there_is_intersection1 )
            ++numTangentCircles;

        if ( there_is_intersection2 )
            ++numTangentCircles;

        if ( numTangentCircles == 1 && there_is_intersection2 ) {
            tangentCircle1 = tangentCircle2;
            tangentCircle2 = rg_Circle2D();
        }
    }
    else {
        bool there_is_intersection  = findIntersectionPtWithinThisArrange(line, disk_base, another_disk, line.getSP(), line.getEP(), tangentCircle1);

        if ( there_is_intersection )
            ++numTangentCircles;
    }

    return numTangentCircles;
    */

    // NOTE: If the signed distance between input line and the center of one disk is negative,
    //       input line will be reversed.
    //       Translate the line in order to reflect the shrinkage of the disk to a point.
    rg_Line2D lineWithCorrectDirection;
    double signed_distance_center_1 = line.signed_distance( disk1.getCenterPt() );
    double signed_distance_center_2 = line.signed_distance( disk2.getCenterPt() );

    bool b_disk1_is_intersected_with_line = rg_LE( abs( signed_distance_center_1 ), disk1.getRadius() ) ? true : false;
    bool b_disk2_is_intersected_with_line = rg_LE( abs( signed_distance_center_2 ), disk2.getRadius() ) ? true : false;

    if ( !b_disk1_is_intersected_with_line && !b_disk2_is_intersected_with_line && ( signed_distance_center_1 * signed_distance_center_2 <= 0.0 ) ) {
        return 0;
    }

    if ( abs( signed_distance_center_1 ) > disk1.getRadius() ) {
        if ( signed_distance_center_1 < 0.0 ) {
            lineWithCorrectDirection = line.get_reversed_line2D();
        }
        else {
            lineWithCorrectDirection = line;
        }
    }
    else if ( abs( signed_distance_center_2 ) > disk2.getRadius() ) {
        if ( signed_distance_center_2 < 0.0 ) {
            lineWithCorrectDirection = line.get_reversed_line2D();
        }
        else {
            lineWithCorrectDirection = line;
        }
    }
    else {
        // XXX # tangent circle could be 4.
        lineWithCorrectDirection = line.get_reversed_line2D();
        lineWithCorrectDirection = line;
    }

    // Version: Mathematically correct solution without using transformation matrix
    rg_Point2D  SP = lineWithCorrectDirection.getSP();
    rg_Point2D  EP = lineWithCorrectDirection.getEP();
    rg_Point2D  vec_dir = ( EP - SP ).getUnitVector();
    rg_Point2D  vec_norm = lineWithCorrectDirection.getNormalVector().getUnitVector();

    // Create the directrix lines which reflect the shrinkage of disks.
    rg_Line2D   moved_line1, moved_line2;
    moved_line1.setSP( SP - vec_norm * disk1.getRadius() );
    moved_line1.setEP( EP - vec_norm * disk1.getRadius() );
    moved_line2.setSP( SP - vec_norm * disk2.getRadius() );
    moved_line2.setEP( EP - vec_norm * disk2.getRadius() );


    // We define a parabola between Y = 4 * p * X^2 when the minimum of Y coinsides the origin.
    // p_1 and p_2 play the role of "p" in the formula.
    rg_Point2D  footprint1;
    rg_Point2D  footprint2;
    moved_line1.compute_perpendicular_footprint_of_point_onto_entire_line( disk1.getCenterPt(), footprint1 );
    moved_line2.compute_perpendicular_footprint_of_point_onto_entire_line( disk2.getCenterPt(), footprint2 );

    rg_Point2D  min_loc_bisector1 = ( footprint1 + disk1.getCenterPt() ) / 2.0; // min_loc_bisector1
    rg_Point2D  min_loc_bisector2 = ( footprint2 + disk2.getCenterPt() ) / 2.0;

    double      p_1 = disk1.getCenterPt().distance( min_loc_bisector1 );
    double      p_2 = disk2.getCenterPt().distance( min_loc_bisector2 );

    //////// HERE 1 ////////////////////////////
    //
    // HERE 1 cannot compute right V-vertex coordinate where one disk center is below the input line.
    /*
    //rg_Point2D  footPrint1;
    //rg_Point2D  footPrint2;
    //line.compute_perpendicular_footprint_of_point_onto_entire_line( disk1.getCenterPt(), footPrint1 );
    //line.compute_perpendicular_footprint_of_point_onto_entire_line( disk2.getCenterPt(), footPrint2 );

    //double      distance_between_line_N_disk1 = disk1.distance( footPrint1 );
    //double      distance_between_line_N_disk2 = disk2.distance( footPrint2 );

    //rg_Point2D  vertex1     = footPrint1 + vec_norm * distance_between_line_N_disk1 / 2.0;
    //rg_Point2D  vertex2     = footPrint2 + vec_norm * distance_between_line_N_disk2 / 2.0;

    //double      p_1         = distance_between_line_N_disk1 / 2.0 + disk1.getRadius();
    //double      p_2         = distance_between_line_N_disk2 / 2.0 + disk2.getRadius();
    */
    //
    //
    ///////// HERE 1 ////////////////////////////


    // Suppose that we transform line and two disks like vertex1 is origin and line is parallel to X-axis.
    // We compute intersection point with above conditions, and then we back transform the result(s).
    // We use determinant to get difference between vertex1 and vertex (transformed).
    // vertex1  = (x_1, y_1)
    // vertex2  = (x_2, y_2)
    // vec_dir  = (x_d, y_d)
    // vec_norm = (x_n, y_n)
    // x_1  +  x_d * diff_x_transformed  +  x_n * diff_y_transformed = x_2
    // y_1  +  y_d * diff_x_transformed  +  y_n * diff_y_transformed = y_2
    //
    //         | x_d   x_n |                          1            |  y_n  -x_n |
    //  A   =  |           |     A^-1   =  ----------------------- |            |    A * transformed_difference = difference   transformed_difference = A^-1 * difference
    //         | y_d   y_n |                x_d * y_n - x_n * y_d  | -y_d   x_d |
    double      x_d = vec_dir.getX();
    double      y_d = vec_dir.getY();
    double      x_n = vec_norm.getX();
    double      y_n = vec_norm.getY();

    double      diff_x_vertex = min_loc_bisector2.getX() - min_loc_bisector1.getX();
    double      diff_y_vertex = min_loc_bisector2.getY() - min_loc_bisector1.getY();
    double      x_line_start = lineWithCorrectDirection.getSP().getX() - min_loc_bisector1.getX();
    double      y_line_start = lineWithCorrectDirection.getSP().getY() - min_loc_bisector1.getY();
    double      x_line_end = lineWithCorrectDirection.getEP().getX() - min_loc_bisector1.getX();
    double      y_line_end = lineWithCorrectDirection.getEP().getY() - min_loc_bisector1.getY();

    double      ad_bc = vec_dir.getX() * vec_norm.getY() - vec_norm.getX() * vec_dir.getY();
    double      diff_x_vertex_transformed = ( y_n * diff_x_vertex - x_n * diff_y_vertex ) / ad_bc;
    double      diff_y_vertex_transformed = ( -y_d * diff_x_vertex + x_d * diff_y_vertex ) / ad_bc;
    double      x_line_start_transformed = ( y_n * x_line_start - x_n * y_line_start ) / ad_bc;
    double      x_line_end_transformed = ( y_n * x_line_end - x_n * y_line_end ) / ad_bc;

    rg_Point2D  intersectionPt[2];
    int         numIntersectionPt = 0;

    // The equation of intersection point is 
    // ( 1 - (p_2 / p_1) ) * x^2   -   2 * x_diff * x   +   (x_diff^2 + 4*p_2*y_diff)   =   0   
    // --> a * x^2  +  2 * b * x  +  c  =  0
    double      a_for_root = 1.0 - ( p_2 / p_1 );
    double      b_for_root = -diff_x_vertex_transformed;
    double      c_for_root = diff_x_vertex_transformed * diff_x_vertex_transformed + 4.0 * p_2 * diff_y_vertex_transformed;

    double      b_square_minus_ac = b_for_root * b_for_root - a_for_root * c_for_root;

    double      x_transformed[2] = { 0.0, 0.0 };
    double      y_transformed[2] = { 0.0, 0.0 };

    // We can get y_transformed by putting x_transformed to equation x^2 = 4 * p_1 * y;
    // OCT 5, 2020 discussion
    if ( rg_EQ( a_for_root, 0.0 ) ) {
        // This case includes CASE 1 and 2.
        numIntersectionPt = 1;
        x_transformed[0] = -c_for_root / 2.0 / b_for_root;
        y_transformed[0] = x_transformed[0] * x_transformed[0] / 4.0 / p_1;
    }
    else {
        if ( rg_EQ( b_square_minus_ac, 0.0 ) ) {
            // Two parabolas have a single, tangential contact.
            // The narrow parabola is contained in the wide one.
            numIntersectionPt = 1;
            x_transformed[0] = -b_for_root / a_for_root;
            y_transformed[0] = x_transformed[0] * x_transformed[0] / 4.0 / p_1;
        }
        else if ( rg_GT( b_square_minus_ac, 0.0 ) ) {
            numIntersectionPt = 2;
            double sqrt_b_square_minus_ac = sqrt( b_square_minus_ac );

            x_transformed[0] = ( -b_for_root + sqrt_b_square_minus_ac ) / a_for_root;
            x_transformed[1] = ( -b_for_root - sqrt_b_square_minus_ac ) / a_for_root;

            y_transformed[0] = x_transformed[0] * x_transformed[0] / 4.0 / p_1;
            y_transformed[1] = x_transformed[1] * x_transformed[1] / 4.0 / p_1;
        }
        else { // if ( rg_LT( b_square_minus_ac, 0.0 ) ) {
            // We believe this case occurs only when the bigger disk containes the smaller one.
            // By definition, this case is excluded from the problem scope.
            return 0; // There have to be at least one intersection point where they share one directrix
        }
    }


    // Check if projection of intersection point is on the line segments
    if ( numIntersectionPt == 2 ) {
        //if ( x_transformed[1] < x_line_start_transformed || x_transformed[1] > x_line_end_transformed ) {
        if ( rg_LT( x_transformed[1], x_line_start_transformed ), rg_GT( x_transformed[1], x_line_end_transformed) ) {
            --numIntersectionPt;
        }
    }


    //if ( x_transformed[0] < x_line_start_transformed || x_transformed[0] > x_line_end_transformed ) {
    if ( rg_LT( x_transformed[0], x_line_start_transformed ), rg_GT( x_transformed[0], x_line_end_transformed ) ) {
        --numIntersectionPt;
        x_transformed[0] = x_transformed[1];
        y_transformed[0] = y_transformed[1];
    }


    // back transform
    if ( numIntersectionPt > 0 ) {
        intersectionPt[0] = min_loc_bisector1 + x_transformed[0] * vec_dir + y_transformed[0] * vec_norm;
        tangentCircle1.setCenterPt( intersectionPt[0] );
        tangentCircle1.setRadius( y_transformed[0] + p_1 - disk1.getRadius() );
    }

    if ( numIntersectionPt == 2 ) {
        intersectionPt[1] = min_loc_bisector1 + x_transformed[1] * vec_dir + y_transformed[1] * vec_norm;
        tangentCircle2.setCenterPt( intersectionPt[1] );
        tangentCircle2.setRadius( y_transformed[1] + p_1 - disk1.getRadius() );
    }

    if ( numIntersectionPt == 0 )
    {
        int stopHere = 1;
    }

    return numIntersectionPt;

}



int PolygonVD2D::computeTangentCircles_of_a_disk_and_two_lines( const rg_Circle2D& disk, const rg_Line2D& line1, const rg_Line2D& line2, vector<rg_Circle2D>& tangentCircles )
{
    rg_Line2D bisector_between_two_line_segments[2] = {
        rg_GeoFunc::compute_bisector_line_between_two_line_segments( line1, line2 ),
        rg_GeoFunc::compute_bisector_line_between_two_line_segments( line1, line2.get_reversed_line2D() ) };

    rg_Point2D  SP[2] = { bisector_between_two_line_segments[0].getSP(), bisector_between_two_line_segments[1].getSP() };
    rg_Point2D  EP[2] = { bisector_between_two_line_segments[0].getEP(), bisector_between_two_line_segments[1].getEP() };

    double radius_SP[2] = { line1.signed_distance( SP[0] ), line1.signed_distance( SP[1] ) };
    double radius_EP[2] = { line1.signed_distance( EP[0] ), line1.signed_distance( EP[1] ) };

    double signed_distance_between_centerPt_N_line1 = line1.signed_distance( disk.getCenterPt() );
    double signed_distance_between_centerPt_N_line2 = line2.signed_distance( disk.getCenterPt() );

    bool twoLinesAreParallel = false;
    rg_Point2D intersectionPtOfTwoLines = line1.compute_intersection_with_line( line2, twoLinesAreParallel );

    bool disk_is_intersected_with_line1 = ( abs( signed_distance_between_centerPt_N_line1 ) < disk.getRadius() ) ? true : false;
    bool disk_is_intersected_with_line2 = ( abs( signed_distance_between_centerPt_N_line2 ) < disk.getRadius() ) ? true : false;

    if ( twoLinesAreParallel ) {
        rg_Circle2D circumcircle1, circumcircle2;
        if ( disk_is_intersected_with_line1 || disk_is_intersected_with_line2 ) {
            computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[0], rg_Circle2D( SP[0], abs( radius_SP[0] ) ), rg_Circle2D( EP[0], abs( radius_EP[0] ) ), disk, circumcircle1, circumcircle2 );
            tangentCircles.push_back( circumcircle1 );
            tangentCircles.push_back( circumcircle2 );
            return 2;
        }
        else {
            double signedDistanceBetweenSP_N_line2 = line2.signed_distance( SP[0] );
            if ( signed_distance_between_centerPt_N_line1 * radius_SP[0] > 0.0 && signed_distance_between_centerPt_N_line2 * signedDistanceBetweenSP_N_line2 > 0.0 ) {
                computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[0], rg_Circle2D( SP[0], abs( radius_SP[0] ) ), rg_Circle2D( EP[0], abs( radius_EP[0] ) ), disk, circumcircle1, circumcircle2 );
                tangentCircles.push_back( circumcircle1 );
                tangentCircles.push_back( circumcircle2 );
                return 2;
            }
            else {
                return 0;
            }
        }
    }
    else {
        bool disk_includes_intersection_pt = ( intersectionPtOfTwoLines.distance( disk.getCenterPt() ) < disk.getRadius() ) ? true : false;
        // intersection ?? ???? case ?��? ???

        if ( !disk_is_intersected_with_line1 && !disk_is_intersected_with_line2 ) {
            int index_of_bisector = 0;
            if ( signed_distance_between_centerPt_N_line1 * signed_distance_between_centerPt_N_line2 < 0.0 ) {
                index_of_bisector = 0;
            }
            else {
                index_of_bisector = 1;
            }

            rg_Circle2D circumcircle1, circumcircle2;
            if ( signed_distance_between_centerPt_N_line1 < 0.0 ) {
                computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[index_of_bisector], rg_Circle2D( SP[index_of_bisector], -radius_SP[index_of_bisector] ), rg_Circle2D( EP[index_of_bisector], -radius_EP[index_of_bisector] ), disk, circumcircle1, circumcircle2 );
            }
            else {
                computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[index_of_bisector], rg_Circle2D( SP[index_of_bisector], radius_SP[index_of_bisector] ), rg_Circle2D( EP[index_of_bisector], radius_EP[index_of_bisector] ), disk, circumcircle1, circumcircle2 );
            }
            tangentCircles.push_back( circumcircle1 );
            tangentCircles.push_back( circumcircle2 );

            return 2;
        }
        else {
            rg_Circle2D circumcircle1_bisector0_positive;
            rg_Circle2D circumcircle2_bisector0_positive;
            int numCircumcircles_b0_P = computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[0], rg_Circle2D( SP[0], radius_SP[0] ), rg_Circle2D( EP[0], radius_EP[0] ), disk, circumcircle1_bisector0_positive, circumcircle2_bisector0_positive );

            rg_Circle2D circumcircle1_bisector0_negative;
            rg_Circle2D circumcircle2_bisector0_negative;
            int numCircumcircles_b0_N = computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[0], rg_Circle2D( SP[0], -radius_SP[0] ), rg_Circle2D( EP[0], -radius_EP[0] ), disk, circumcircle1_bisector0_negative, circumcircle2_bisector0_negative );

            rg_Circle2D circumcircle1_bisector1_positive;
            rg_Circle2D circumcircle2_bisector1_positive;
            int numCircumcircles_b1_P = computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[1], rg_Circle2D( SP[1], radius_SP[1] ), rg_Circle2D( EP[1], radius_EP[1] ), disk, circumcircle1_bisector1_positive, circumcircle2_bisector1_positive );

            rg_Circle2D circumcircle1_bisector1_negative;
            rg_Circle2D circumcircle2_bisector1_negative;
            int numCircumcircles_b1_N = computeCoordOfNewVVertex_on_line_of_two_polygon_edges( bisector_between_two_line_segments[1], rg_Circle2D( SP[1], -radius_SP[1] ), rg_Circle2D( EP[1], -radius_EP[1] ), disk, circumcircle1_bisector1_negative, circumcircle2_bisector1_negative );

            if ( numCircumcircles_b0_P != 0 ) {
                if ( numCircumcircles_b0_P > 0 && circumcircle1_bisector0_positive.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle1_bisector0_positive );
                }
                if ( numCircumcircles_b0_P > 1 && circumcircle2_bisector0_positive.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle2_bisector0_positive );
                }
            }

            if ( numCircumcircles_b0_N != 0 ) {
                if ( numCircumcircles_b0_N > 0 && circumcircle1_bisector0_negative.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle1_bisector0_negative );
                }
                if ( numCircumcircles_b0_N > 1 && circumcircle2_bisector0_negative.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle2_bisector0_negative );
                }
            }

            if ( numCircumcircles_b1_P != 0 ) {
                if ( numCircumcircles_b1_P > 0 && circumcircle1_bisector1_positive.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle1_bisector1_positive );
                }
                if ( numCircumcircles_b1_P > 1 && circumcircle2_bisector1_positive.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle2_bisector1_positive );
                }
            }

            if ( numCircumcircles_b1_N != 0 ) {
                if ( numCircumcircles_b1_N > 0 && circumcircle1_bisector1_negative.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle1_bisector1_negative );
                }
                if ( numCircumcircles_b1_N > 1 && circumcircle2_bisector1_negative.getRadius() >= 0.0 ) {
                    tangentCircles.push_back( circumcircle2_bisector1_negative );
                }
            }

            return tangentCircles.size();
        }
    }
}



int V::GeometryTier::PolygonVD2D::computeTangentCircles_of_two_disks_and_a_line_for_coordinate( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, vector<rg_Circle2D>& tangentCircles )
{
    tangentCircles.clear();

    // NOTE: If the signed distance between input line and the center of one disk is negative,
    //       input line will be reversed.
    //       Translate the line in order to reflect the shrinkage of the disk to a point.
    rg_Line2D lineWithCorrectDirection;
    double signed_distance_center_1 = line.signed_distance( disk1.getCenterPt() );
    double signed_distance_center_2 = line.signed_distance( disk2.getCenterPt() );

    bool b_disk1_is_intersected_with_line = rg_LT( abs( signed_distance_center_1 ), disk1.getRadius() ) ? true : false;
    bool b_disk2_is_intersected_with_line = rg_LT( abs( signed_distance_center_2 ), disk2.getRadius() ) ? true : false;

    if ( !b_disk1_is_intersected_with_line && !b_disk2_is_intersected_with_line && ( signed_distance_center_1 * signed_distance_center_2 <= 0.0 ) ) {
        return 0;
    }

    //if ( abs( signed_distance_center_1 ) > disk1.getRadius() ) {
    if ( !b_disk1_is_intersected_with_line ) {
        if ( signed_distance_center_1 < 0.0 ) {
            lineWithCorrectDirection = line.get_reversed_line2D();
        }
        else {
            lineWithCorrectDirection = line;
        }

        rg_Circle2D currTangentCircles[2];
        int numTangentCircles = computeTangentCircles_of_two_disks_and_an_oriented_line( disk1, disk2, lineWithCorrectDirection, currTangentCircles[0], currTangentCircles[1] );

        for ( int i = 0; i < numTangentCircles; ++i ) {
            tangentCircles.push_back( currTangentCircles[i] );
        }
    }
    //else if ( abs( signed_distance_center_2 ) > disk2.getRadius() ) {
    else if ( !b_disk2_is_intersected_with_line ) {
        if ( signed_distance_center_2 < 0.0 ) {
            lineWithCorrectDirection = line.get_reversed_line2D();
        }
        else {
            lineWithCorrectDirection = line;
        }

        rg_Circle2D currTangentCircles[2];
        int numTangentCircles = computeTangentCircles_of_two_disks_and_an_oriented_line( disk1, disk2, lineWithCorrectDirection, currTangentCircles[0], currTangentCircles[1] );

        for ( int i = 0; i < numTangentCircles; ++i ) {
            tangentCircles.push_back( currTangentCircles[i] );
        }
    }
    else {
        rg_Circle2D tangentCircles_forward_directed_line[2];
        int numTangentCircles_forward_directed_line = computeTangentCircles_of_two_disks_and_an_oriented_line( disk1, disk2, line, tangentCircles_forward_directed_line[0], tangentCircles_forward_directed_line[1] );
        for ( int i = 0; i < numTangentCircles_forward_directed_line; ++i ) {
            tangentCircles.push_back( tangentCircles_forward_directed_line[i] );
        }

        rg_Circle2D tangentCircles_backward_directed_line[2];
        int numTangentCircles_backward_directed_line = computeTangentCircles_of_two_disks_and_an_oriented_line( disk1, disk2, line.get_reversed_line2D(), tangentCircles_backward_directed_line[0], tangentCircles_backward_directed_line[1] );
        for ( int i = 0; i < numTangentCircles_backward_directed_line; ++i ) {
            tangentCircles.push_back( tangentCircles_backward_directed_line[i] );
        }
    }

    return tangentCircles.size();
}



int V::GeometryTier::PolygonVD2D::computeTangentCircles_of_two_disks_and_an_oriented_line( const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D& line, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    // Version: Mathematically correct solution without using transformation matrix
    rg_Point2D  SP = line.getSP();
    rg_Point2D  EP = line.getEP();
    rg_Point2D  vec_dir = ( EP - SP ).getUnitVector();
    rg_Point2D  vec_norm = line.getNormalVector().getUnitVector();

    // Create the directrix lines which reflect the shrinkage of disks.
    rg_Line2D   moved_line1, moved_line2;
    moved_line1.setSP( SP - vec_norm * disk1.getRadius() );
    moved_line1.setEP( EP - vec_norm * disk1.getRadius() );
    moved_line2.setSP( SP - vec_norm * disk2.getRadius() );
    moved_line2.setEP( EP - vec_norm * disk2.getRadius() );


    // We define a parabola between Y = 4 * p * X^2 when the minimum of Y coinsides the origin.
    // p_1 and p_2 play the role of "p" in the formula.
    rg_Point2D  footprint1;
    rg_Point2D  footprint2;
    moved_line1.compute_perpendicular_footprint_of_point_onto_entire_line( disk1.getCenterPt(), footprint1 );
    moved_line2.compute_perpendicular_footprint_of_point_onto_entire_line( disk2.getCenterPt(), footprint2 );

    rg_Point2D  min_loc_bisector1 = ( footprint1 + disk1.getCenterPt() ) / 2.0; // min_loc_bisector1
    rg_Point2D  min_loc_bisector2 = ( footprint2 + disk2.getCenterPt() ) / 2.0;

    double      p_1 = disk1.getCenterPt().distance( min_loc_bisector1 );
    double      p_2 = disk2.getCenterPt().distance( min_loc_bisector2 );

    // Suppose that we transform line and two disks like vertex1 is origin and line is parallel to X-axis.
    // We compute intersection point with above conditions, and then we back transform the result(s).
    // We use determinant to get difference between vertex1 and vertex (transformed).
    // vertex1  = (x_1, y_1)
    // vertex2  = (x_2, y_2)
    // vec_dir  = (x_d, y_d)
    // vec_norm = (x_n, y_n)
    // x_1  +  x_d * diff_x_transformed  +  x_n * diff_y_transformed = x_2
    // y_1  +  y_d * diff_x_transformed  +  y_n * diff_y_transformed = y_2
    //
    //         | x_d   x_n |                          1            |  y_n  -x_n |
    //  A   =  |           |     A^-1   =  ----------------------- |            |    A * transformed_difference = difference   transformed_difference = A^-1 * difference
    //         | y_d   y_n |                x_d * y_n - x_n * y_d  | -y_d   x_d |
    double      x_d = vec_dir.getX();
    double      y_d = vec_dir.getY();
    double      x_n = vec_norm.getX();
    double      y_n = vec_norm.getY();

    double      diff_x_vertex = min_loc_bisector2.getX() - min_loc_bisector1.getX();
    double      diff_y_vertex = min_loc_bisector2.getY() - min_loc_bisector1.getY();
    double      x_line_start = SP.getX() - min_loc_bisector1.getX();
    double      y_line_start = SP.getY() - min_loc_bisector1.getY();
    double      x_line_end = EP.getX() - min_loc_bisector1.getX();
    double      y_line_end = EP.getY() - min_loc_bisector1.getY();

    double      ad_bc = vec_dir.getX() * vec_norm.getY() - vec_norm.getX() * vec_dir.getY();
    double      diff_x_vertex_transformed = ( y_n * diff_x_vertex - x_n * diff_y_vertex ) / ad_bc;
    double      diff_y_vertex_transformed = ( -y_d * diff_x_vertex + x_d * diff_y_vertex ) / ad_bc;
    double      x_line_start_transformed = ( y_n * x_line_start - x_n * y_line_start ) / ad_bc;
    double      x_line_end_transformed = ( y_n * x_line_end - x_n * y_line_end ) / ad_bc;

    rg_Point2D  intersectionPt[2];
    int         numIntersectionPt = 0;

    // The equation of intersection point is 
    // ( 1 - (p_2 / p_1) ) * x^2   -   2 * x_diff * x   +   (x_diff^2 + 4*p_2*y_diff)   =   0   
    // --> a * x^2  +  2 * b * x  +  c  =  0
    double      a_for_root = 1.0 - ( p_2 / p_1 );
    double      b_for_root = -diff_x_vertex_transformed;
    double      c_for_root = diff_x_vertex_transformed * diff_x_vertex_transformed + 4.0 * p_2 * diff_y_vertex_transformed;

    double      b_square_minus_ac = b_for_root * b_for_root - a_for_root * c_for_root;

    double      x_transformed[2] = { 0.0, 0.0 };
    double      y_transformed[2] = { 0.0, 0.0 };

    // We can get y_transformed by putting x_transformed to equation x^2 = 4 * p_1 * y;
    // OCT 5, 2020 discussion
    if ( rg_EQ( a_for_root, 0.0 ) ) {
        // This case includes CASE 1 and 2.
        numIntersectionPt = 1;
        x_transformed[0] = -c_for_root / 2.0 / b_for_root;
        y_transformed[0] = x_transformed[0] * x_transformed[0] / 4.0 / p_1;
    }
    else {
        if ( rg_EQ( b_square_minus_ac, 0.0 ) ) {
            // Two parabolas have a single, tangential contact.
            // The narrow parabola is contained in the wide one.
            numIntersectionPt = 1;
            x_transformed[0] = -b_for_root / a_for_root;
            y_transformed[0] = x_transformed[0] * x_transformed[0] / 4.0 / p_1;
        }
        else if ( rg_GT( b_square_minus_ac, 0.0 ) ) {
            numIntersectionPt = 2;
            double sqrt_b_square_minus_ac = sqrt( b_square_minus_ac );

            x_transformed[0] = ( -b_for_root + sqrt_b_square_minus_ac ) / a_for_root;
            x_transformed[1] = ( -b_for_root - sqrt_b_square_minus_ac ) / a_for_root;

            y_transformed[0] = x_transformed[0] * x_transformed[0] / 4.0 / p_1;
            y_transformed[1] = x_transformed[1] * x_transformed[1] / 4.0 / p_1;
        }
        else { // if ( rg_LT( b_square_minus_ac, 0.0 ) ) {
            // We believe this case occurs only when the bigger disk containes the smaller one.
            // By definition, this case is excluded from the problem scope.
            return 0; // There have to be at least one intersection point where they share one directrix
        }
    }


    // Check if projection of intersection point is on the line segments
    if ( numIntersectionPt == 2 ) {
        //if ( x_transformed[1] < x_line_start_transformed || x_transformed[1] > x_line_end_transformed ) {
        if ( rg_LT( x_transformed[1], x_line_start_transformed ), rg_GT( x_transformed[1], x_line_end_transformed ) ) {
            --numIntersectionPt;
        }
    }


    //if ( x_transformed[0] < x_line_start_transformed || x_transformed[0] > x_line_end_transformed ) {
    if ( rg_LT( x_transformed[0], x_line_start_transformed ), rg_GT( x_transformed[0], x_line_end_transformed ) ) {
        --numIntersectionPt;
        x_transformed[0] = x_transformed[1];
        y_transformed[0] = y_transformed[1];
    }


    // back transform
    if ( numIntersectionPt > 0 ) {
        intersectionPt[0] = min_loc_bisector1 + x_transformed[0] * vec_dir + y_transformed[0] * vec_norm;
        tangentCircle1.setCenterPt( intersectionPt[0] );
        tangentCircle1.setRadius( y_transformed[0] + p_1 - disk1.getRadius() );
    }

    if ( numIntersectionPt == 2 ) {
        intersectionPt[1] = min_loc_bisector1 + x_transformed[1] * vec_dir + y_transformed[1] * vec_norm;
        tangentCircle2.setCenterPt( intersectionPt[1] );
        tangentCircle2.setRadius( y_transformed[1] + p_1 - disk1.getRadius() );
    }

    if ( numIntersectionPt == 0 )
    {
        int stopHere = 1;
    }

    return numIntersectionPt;
}



rg_Circle2D V::GeometryTier::PolygonVD2D::chooseTheCorrectTangentCircle_whichHasCorrectTopologicalOrientation( const int& numTangentCircles, const rg_Circle2D& tangentCircle1, const rg_Circle2D& tangentCircle2, VVertex2D* const newVVertex )
{
    rg_Circle2D correctTangentCircle;

    if ( numTangentCircles == 1 ) {
        correctTangentCircle = tangentCircle1;
    }
    else {
        rg_Circle2D firstTangentCircle, secondTangentCircle;
        if ( tangentCircle1.getRadius() > tangentCircle2.getRadius() ) {
            firstTangentCircle = tangentCircle2;
            secondTangentCircle = tangentCircle1;
        }
        else {
            firstTangentCircle = tangentCircle1;
            secondTangentCircle = tangentCircle2;
        }

        VVertex2D* currVVertex = const_cast<VVertex2D*>( newVVertex );
        if ( this_circumcircle_has_right_orientation( currVVertex, firstTangentCircle ) ) {
            correctTangentCircle = firstTangentCircle;
        }
        else {
            correctTangentCircle = secondTangentCircle;
        }
    }

    return correctTangentCircle;
}



void PolygonVD2D::classifyVertexGenerators()
{
    for ( list<Generator2D*>::iterator i_gen = m_generators.begin(); i_gen != m_generators.end(); ++i_gen ) {
        Generator2D* currGen = *i_gen;

        if ( currGen->getType() != Generator2D::Generator_Type::VERTEX_G )
            continue;

        VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;

        if ( vertexGen->isThisFromContainer() ) {
            continue;
        }

        EdgeGenerator2D* prevEdgeGen = vertexGen->get_previous_edge_generator();
        EdgeGenerator2D* nextEdgeGen = vertexGen->get_next_edge_generator();

        if ( prevEdgeGen == NULL || nextEdgeGen == NULL ) {
            continue;
        }

        rg_Point2D currPt = vertexGen->get_point();
        rg_Point2D prevPt = prevEdgeGen->get_start_vertex_point();
        rg_Point2D nextPt = nextEdgeGen->get_end_vertex_point();

        if ( rg_GeoFunc::left_turn( prevPt, currPt, nextPt ) ) {
            vertexGen->set_vertex_type( VertexGenerator2D::Vertex_Type::NON_REFLEX_FROM_POLYGON_INSIDE );
        }
        else {
            vertexGen->set_vertex_type( VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE );
        }
    }
}



DiskGenerator2D * PolygonVD2D::insertThisDiskInPolygonVoronoiDiagram( const rg_Circle2D & disk )
{
    /*
    ?? Polygon Vertex ???? ????? VVertex (on-vertex)?? ???��?.
    on-vertex?? findSeedRedVertex, wavePropagation_ver1 ??????? ????? red vertex?? ????? ????.
    ??????, ??????? ?????? violate ?? ?? ???. 
    (?????? ???? red vertex?? ???? ????? ?????????? ?????? ?????? ????).
    (Rectangle?????? ??????? ??????, reflex vertex???? ???? ?? ??? ??????).
    */
    DiskGenerator2D* newGenerator  = createGeneratorOfInsertedDisk(disk);
    
    // This is same in the incremental process of the TOI_D2 algorithm.
    // From here (1)
    VFace2D*         anchorCell    = findAnchorCell(disk.getCenterPt());    
    VVertex2D*       seedRedVertex = findSeedRedVertex(newGenerator, (Generator2D*)anchorCell->getGenerator());

    list<VVertex2D*> redVVertices;
    list<VVertex2D*> blueVVertices;
    list<VVertex2D*> fictitiousVVertices;

    seedRedVertex->setStatus(RED_V);
    redVVertices.push_back(seedRedVertex);

    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );

    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeNewCellAndConnectToCurrVD( newGenerator, redVVertices, blueVVertices, newVVertices, newVEdges );
    // To here (1)


    // This is unique in the VD of simple polygon.
    // polygonVD?? VVertex?? VEdge?? Polygon?? ????? ??????? ???????.
    // Polygon ???��? disk?? insert??? ????????? new VVertices?? new VEdges?? ?????? Polygon ????? ?????????.
    // ??? ????????? newVVertices?? newVEdges?? polygon ????? ???? marking?? ???.
    markLocationStatusAtNewEntities( newVVertices, newVEdges );

    
    // This is same in the incremental process of the TOI_D2 algorithm.
    // From here (2)
    connectCurrVDToNewVEdgeLoop( newVEdges );
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
    // To here (2)

    // This is unique in the VD of simple polygon
    // Geometry computation: ?? ???? ??? ???? topology operation??.
    // Edge geometry?? ??????? ????.
    computeCoordOfNewVVertices( newVVertices );

    
    // This is same in the incremental process of the TOI_D2 algorithm.
    removeAllExtraneousVVerticesAndVEdges();

    
    // disk packing???? ????? map
    // Voronoi diagram?????? disk?? precision?? single?? ???????.
    // disk packing?????????? double precision?? disk?? ?????? intersection ???��? ?????? ??? ??????
    // VD ?????? ??? double precision?? disk?? VD ?????? disk?? mapping??? ???.
    m_map_DiskID2VFace.insert( make_pair( newGenerator->getID(), newGenerator->getOuterFace() ) );

    return newGenerator;
}



DiskGenerator2D * PolygonVD2D::insertThisDiskInPolygonVoronoiDiagram( const rg_Circle2D & disk, list<VVertex2D*>& redVertices, list<VVertex2D*>& newVertices )
{
    /*
    ?? Polygon Vertex ???? ????? VVertex (on-vertex)?? ???��?.
    on-vertex?? findSeedRedVertex, wavePropagation_ver1 ??????? ????? red vertex?? ????? ????.
    ??????, ??????? ?????? violate ?? ?? ???.
    (?????? ???? red vertex?? ???? ????? ?????????? ?????? ?????? ????).
    (Rectangle?????? ??????? ??????, reflex vertex???? ???? ?? ??? ??????).
    */
    DiskGenerator2D* newGenerator = createGeneratorOfInsertedDisk( disk );

    // This is same in the incremental process of the TOI_D2 algorithm.
    // From here (1)
    VFace2D*         anchorCell = findAnchorCell( disk.getCenterPt() );
    VVertex2D*       seedRedVertex = findSeedRedVertex( newGenerator, (Generator2D*)anchorCell->getGenerator() );

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


    // This is unique in the VD of simple polygon.
    // polygonVD?? VVertex?? VEdge?? Polygon?? ????? ??????? ???????.
    // Polygon ???��? disk?? insert??? ????????? new VVertices?? new VEdges?? ?????? Polygon ????? ?????????.
    // ??? ????????? newVVertices?? newVEdges?? polygon ????? ???? marking?? ???.
    markLocationStatusAtNewEntities( newVVertices, newVEdges );


    // This is same in the incremental process of the TOI_D2 algorithm.
    // From here (2)
    connectCurrVDToNewVEdgeLoop( newVEdges );
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
    // To here (2)

    // This is unique in the VD of simple polygon
    // Geometry computation: ?? ???? ??? ???? topology operation??.
    // Edge geometry?? ??????? ????.
    computeCoordOfNewVVertices( newVVertices );


    // This is same in the incremental process of the TOI_D2 algorithm.
    removeAllExtraneousVVerticesAndVEdges();


    // disk packing???? ????? map
    // Voronoi diagram?????? disk?? precision?? single?? ???????.
    // disk packing?????????? double precision?? disk?? ?????? intersection ???��? ?????? ??? ??????
    // VD ?????? ??? double precision?? disk?? VD ?????? disk?? mapping??? ???.
    m_map_DiskID2VFace.insert( make_pair( newGenerator->getID(), newGenerator->getOuterFace() ) );


    redVertices.insert( redVertices.end(), redVVertices.begin(), redVVertices.end() );
    newVertices.insert( newVertices.end(), newVVertices.begin(), newVVertices.end() );

    return newGenerator;
}



DiskGenerator2D* PolygonVD2D::insert_a_shrinked_disk_on_this_VVertex( const VVertex2D* const vertex, const double& shrink_ratio )
{
    VVertex2D*      seedRedVertex   = const_cast<VVertex2D*>(vertex);
    rg_Point2D      diskCenter      = seedRedVertex->getCircumcircle().getCenterPt();
    double          radius          = seedRedVertex->getCircumcircle().getRadius() * shrink_ratio;
    rg_Circle2D     disk(diskCenter, radius);

    DiskGenerator2D* newGenerator  = createGeneratorOfInsertedDisk(disk);
    
    list<VVertex2D*> redVVertices;
    list<VVertex2D*> blueVVertices;
    list<VVertex2D*> fictitiousVVertices;

    seedRedVertex->setStatus(RED_V);
    redVVertices.push_back(seedRedVertex);
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );

    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeNewCellAndConnectToCurrVD( newGenerator, redVVertices, blueVVertices, newVVertices, newVEdges );

    markLocationStatusAtNewEntities( newVVertices, newVEdges );

    connectCurrVDToNewVEdgeLoop( newVEdges );

    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );

    computeCoordOfNewVVertices( newVVertices );

    removeAllExtraneousVVerticesAndVEdges();

    m_map_DiskID2VFace.insert( make_pair( newGenerator->getID(), newGenerator->getOuterFace() ) );

    return newGenerator;
}



DiskGenerator2D * PolygonVD2D::insert_a_shrinked_disk_on_this_VVertex( const VVertex2D * const vertex, const double & shrink_ratio, list<VVertex2D*>& redVertices, list<VVertex2D*>& newVertices )
{
#ifdef CHECK_COMP_TIME
    clock_t startTime, endTime;
    startTime = clock();
#endif
    VVertex2D*      seedRedVertex = const_cast<VVertex2D*>( vertex );
    rg_Point2D      diskCenter = seedRedVertex->getCircumcircle().getCenterPt();
    double          radius = seedRedVertex->getCircumcircle().getRadius() * shrink_ratio;
    rg_Circle2D     disk( diskCenter, radius );

    DiskGenerator2D* newGenerator = createGeneratorOfInsertedDisk( disk );

    list<VVertex2D*> redVVertices;
    list<VVertex2D*> blueVVertices;
    list<VVertex2D*> fictitiousVVertices;
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_others = t_others + endTime - startTime;

    startTime = clock();
#endif
    seedRedVertex->setStatus( RED_V );
    redVVertices.push_back( seedRedVertex );
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_wave_propagation = t_wave_propagation + endTime - startTime;

    startTime = clock();
#endif
    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeNewCellAndConnectToCurrVD( newGenerator, redVVertices, blueVVertices, newVVertices, newVEdges );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_make_new_cell = t_make_new_cell + endTime - startTime;

    startTime = clock();
#endif
    markLocationStatusAtNewEntities( newVVertices, newVEdges );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_mark_in = t_mark_in + endTime - startTime;

    startTime = clock();
#endif
    connectCurrVDToNewVEdgeLoop( newVEdges );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_connect_topology = t_connect_topology + endTime - startTime;

    startTime = clock();
#endif
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_merge_split = t_merge_split + endTime - startTime;

    startTime = clock();
#endif
    computeCoordOfNewVVertices( newVVertices );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_comp_vertex = t_comp_vertex + endTime - startTime;

    //removeAllExtraneousVVerticesAndVEdges();

    startTime = clock();
#endif
    m_map_DiskID2VFace.insert( make_pair( newGenerator->getID(), newGenerator->getOuterFace() ) );


    redVertices.insert( redVertices.end(), redVVertices.begin(), redVVertices.end() );
    newVertices.insert( newVertices.end(), newVVertices.begin(), newVVertices.end() );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_others = t_others + endTime - startTime;
#endif

    return newGenerator;
}



DiskGenerator2D* PolygonVD2D::insert_a_shrinked_disk_on_this_VVertex_outside_of_polygon( const VVertex2D* const vertex, const double& shrink_ratio, list<VVertex2D*>& redVertices, list<VVertex2D*>& newVertices )
{
#ifdef CHECK_COMP_TIME
    clock_t startTime, endTime;
    startTime = clock();
#endif
    VVertex2D*      seedRedVertex   = const_cast<VVertex2D*>( vertex );
    rg_Point2D      diskCenter      = seedRedVertex->getCircumcircle().getCenterPt();
    double          radius          = seedRedVertex->getCircumcircle().getRadius() * shrink_ratio;
    rg_Circle2D     disk( diskCenter, radius );

    DiskGenerator2D* newGenerator   = createGeneratorOfInsertedDisk( disk );

    list<VVertex2D*> redVVertices;
    list<VVertex2D*> blueVVertices;
    list<VVertex2D*> fictitiousVVertices;
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_others = t_others + endTime - startTime;

    startTime = clock();
#endif
    seedRedVertex->setStatus( RED_V );
    redVVertices.push_back( seedRedVertex );
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_wave_propagation = t_wave_propagation + endTime - startTime;

    startTime = clock();
#endif
    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeNewCellAndConnectToCurrVD( newGenerator, redVVertices, blueVVertices, newVVertices, newVEdges );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_make_new_cell = t_make_new_cell + endTime - startTime;

    startTime = clock();
#endif
    markLocationStatusAtNewEntities_outside( newVVertices, newVEdges );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_mark_in = t_mark_in + endTime - startTime;

    startTime = clock();
#endif
    connectCurrVDToNewVEdgeLoop( newVEdges );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_connect_topology = t_connect_topology + endTime - startTime;

    startTime = clock();
#endif
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_merge_split = t_merge_split + endTime - startTime;

    startTime = clock();
#endif
    //computeCoordOfNewVVertices_outside( newVVertices );
    computeCoordOfNewVertices(newVVertices);
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_comp_vertex = t_comp_vertex + endTime - startTime;

    //removeAllExtraneousVVerticesAndVEdges();

    startTime = clock();
#endif
    m_map_DiskID2VFace.insert( make_pair( newGenerator->getID(), newGenerator->getOuterFace() ) );


    redVertices.insert( redVertices.end(), redVVertices.begin(), redVVertices.end() );
    newVertices.insert( newVertices.end(), newVVertices.begin(), newVVertices.end() );
#ifdef CHECK_COMP_TIME
    endTime = clock();
    t_others = t_others + endTime - startTime;
#endif

    return newGenerator;
}



PolygonVD2D::VDEntity_Location_Status PolygonVD2D::get_location_status_of_vvertex( const VVertex2D * const vertex )
{
    VVertex_Location_Status_Map::iterator i_vertex_status = m_VVertexToLocationStatus.find(const_cast<VVertex2D*>(vertex));
    if ( i_vertex_status != m_VVertexToLocationStatus.end() )
        return i_vertex_status->second;
    else
        return VDEntity_Location_Status::UNKNOWN_LOCATION;
}



PolygonVD2D::VDEntity_Location_Status PolygonVD2D::get_location_status_of_edge( const VEdge2D * const edge )
{
    VEdge_Location_Status_Map::iterator i_edge_status = m_VEdgeToLocationStatus.find(const_cast<VEdge2D*>(edge));
    if ( i_edge_status != m_VEdgeToLocationStatus.end() )
        return i_edge_status->second;
    else
        return VDEntity_Location_Status::UNKNOWN_LOCATION;
}



void PolygonVD2D::get_disk_generators_which_share_VFace_with_polygon_boundary(list<DiskGenerator2D*>& diskGenerators)
{
    list<EdgeGenerator2D*> edgeGenerators;
    get_edge_generators(edgeGenerators);

    list<EdgeGenerator2D*>::iterator i_edgeGen = edgeGenerators.begin();
    for (; i_edgeGen != edgeGenerators.end(); ++i_edgeGen)
    {
        VFace2D* VFaceOfCurrEdgeGen = (*i_edgeGen)->get_first_son_disk_generator()->getOuterFace();
        list<VEdge2D*> boundingVEdges;
        VFaceOfCurrEdgeGen->getBoundaryVEdges(boundingVEdges);

        list<VEdge2D*>::iterator i_VEdge = boundingVEdges.begin();
        for (; i_VEdge != boundingVEdges.end(); ++i_VEdge)
        {
            VEdge2D* currVEdge = (*i_VEdge);

            if (get_location_status_of_edge(currVEdge) != VDEntity_Location_Status::INSIDE_POLYGON)
                continue;

            VFace2D* oppositeVFaceOfCurrVEdge = rg_NULL;
            if (currVEdge->getRightFace() == VFaceOfCurrEdgeGen)
            {
                oppositeVFaceOfCurrVEdge = currVEdge->getLeftFace();
            }
            else if (currVEdge->getLeftFace() == VFaceOfCurrEdgeGen)
            {
                oppositeVFaceOfCurrVEdge = currVEdge->getRightFace();
            }
            else {}

            Generator2D::Generator_Type genType = ((Generator2D*)(oppositeVFaceOfCurrVEdge->getGenerator()->getUserData()))->getType();
            if (genType == Generator2D::Generator_Type::DISK_G)
            {
                DiskGenerator2D* currDiskGen = (DiskGenerator2D*)(oppositeVFaceOfCurrVEdge->getGenerator()->getUserData());
                diskGenerators.push_back(currDiskGen);
            }
        }
    }
    diskGenerators.sort();
    diskGenerators.unique();
}



void PolygonVD2D::get_edge_generators(list<EdgeGenerator2D*>& edgeGenerators)
{
    list<Generator2D*>::const_iterator i_generator = m_generators.begin();
    for (; i_generator != m_generators.end(); ++i_generator)
    {
        Generator2D::Generator_Type genType = ((Generator2D*)(*i_generator))->getType();
        if (genType == Generator2D::Generator_Type::EDGE_G)
            edgeGenerators.push_back((EdgeGenerator2D*)(*i_generator));
    }
}



void PolygonVD2D::get_vertex_generators(list<VertexGenerator2D*>& vertexGenerators) 
{
    list<Generator2D*>::const_iterator i_generator = m_generators.begin();
    for (; i_generator != m_generators.end(); ++i_generator)
    {
        Generator2D::Generator_Type genType = ((Generator2D*)(*i_generator))->getType();
        if (genType == Generator2D::Generator_Type::VERTEX_G)
            vertexGenerators.push_back((VertexGenerator2D*)(*i_generator));
    }
}




void PolygonVD2D::get_Voronoi_edges_inside_Polygon(list<const VEdge2D*>& VEdgesList)
{
	list<const VEdge2D*> edges;
	VoronoiDiagram2DC::getVoronoiEdges(edges);
	for (list<const VEdge2D*>::iterator i_edge = edges.begin(); i_edge != edges.end(); ++i_edge)
	{
		const VEdge2D* currEdge = *i_edge;

		// discard VEdge outside polygon
		if (get_location_status_of_edge(currEdge) == VDEntity_Location_Status::INSIDE_POLYGON)
        //if ( get_location_status_of_edge( currEdge ) != VDEntity_Location_Status::INFINITE_ENTITY )
			VEdgesList.push_back(currEdge);
	}
}



void PolygonVD2D::get_Voronoi_edges_outside_Polygon( list<const VEdge2D*>& VEdgesList )
{
    list<const VEdge2D*> edges;
    VoronoiDiagram2DC::getVoronoiEdges( edges );
    for ( list<const VEdge2D*>::iterator i_edge = edges.begin(); i_edge != edges.end(); ++i_edge )
    {
        const VEdge2D* currEdge = *i_edge;

        if ( get_location_status_of_edge( currEdge ) == VDEntity_Location_Status::OUTSIDE_POLYGON )
            VEdgesList.push_back( currEdge );
    }
}



rg_RQBzCurve2D PolygonVD2D::get_geometry_of_edge( const VEdge2D * const edge )
{
    Generator2D* leftGenerator  = edge->getLeftFace()->getGenerator();
    Generator2D* rightGenerator = edge->getRightFace()->getGenerator();

    Generator2D::Generator_Type leftType    = leftGenerator->getType();
    Generator2D::Generator_Type rightType   = rightGenerator->getType();

    rg_RQBzCurve2D geometry_of_edge;

    switch ( leftType ) {
    case Generator2D::Generator_Type::VERTEX_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            rg_Point2D sp = edge->getStartVertex()->getLocation();
            rg_Point2D ep = edge->getEndVertex()->getLocation();
            rg_Point2D cp = ( sp + ep ) / 2.0;

            geometry_of_edge.setCtrlPt(0, sp);
            geometry_of_edge.setCtrlPt(1, cp);
            geometry_of_edge.setCtrlPt(2, ep);
            geometry_of_edge.setWeight(0, 1.0);
            geometry_of_edge.setWeight(1, 1.0);
            geometry_of_edge.setWeight(2, 1.0);
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGenerator;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGenerator;

            if ( vertexGen->get_previous_edge_generator() == edgeGen || vertexGen->get_next_edge_generator() == edgeGen ) {
                rg_Point2D sp = edge->getStartVertex()->getLocation();
                rg_Point2D ep = edge->getEndVertex()->getLocation();
                rg_Point2D cp = ( sp + ep ) / 2.0;

                geometry_of_edge.setCtrlPt(0, sp);
                geometry_of_edge.setCtrlPt(1, cp);
                geometry_of_edge.setCtrlPt(2, ep);
                geometry_of_edge.setWeight(0, 1.0);
                geometry_of_edge.setWeight(1, 1.0);
                geometry_of_edge.setWeight(2, 1.0);
            }
            else {
                rg_Point2D focus     = vertexGen->get_point();
                rg_Line2D  directrix = edgeGen->get_geometry();
                Parabola2D parabola(focus, directrix);

                rg_Point2D sp        = edge->getStartVertex()->getLocation();
                rg_Point2D ep        = edge->getEndVertex()->getLocation();
                rg_Point2D tvs       = parabola.get_tangent_vector(sp);
                rg_Point2D tve       = parabola.get_tangent_vector(ep);
                rg_Point2D passingPt = parabola.get_passing_point(sp, ep);

                geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passingPt);

                if ( geometry_of_edge.getWeight(1) < 0 )
                    geometry_of_edge.setWeight(1, geometry_of_edge.getWeight(1) * -1.0);
            }
        }
        break;

        case Generator2D::Generator_Type::DISK_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGenerator;
            DiskGenerator2D*   diskGen   = (DiskGenerator2D*)rightGenerator;

            rg_Point2D  vertex = vertexGen->get_point();
            rg_Circle2D disk   = diskGen->getDisk();

            rg_Point2D sp, ep, tvs, tve, passPoint;

            if( rg_ZERO( disk.getRadius(), resNeg6 ) )
            {
                sp = edge->getStartVertex()->getLocation();
                ep = edge->getEndVertex()->getLocation();

                geometry_of_edge.setDegree(1);
                geometry_of_edge.setCtrlPt(0, sp);
                geometry_of_edge.setCtrlPt(0, ep);
            }
            else
            {
                sp = edge->getStartVertex()->getLocation();
                ep = edge->getEndVertex()->getLocation();
                tvs = getTangentVector( sp, vertex, disk.getCenterPt() );
                tve = getTangentVector( ep, vertex, disk.getCenterPt() );
                passPoint = getPassingPt(rg_Circle2D(vertex, 0.0), disk);

                geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passPoint);
                if( geometry_of_edge.getWeight(1) < 0 )
                    geometry_of_edge.setWeight( 1, geometry_of_edge.getWeight(1) * -1 );
            }
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::EDGE_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGenerator;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGenerator;

            if ( vertexGen->get_previous_edge_generator() == edgeGen || vertexGen->get_next_edge_generator() == edgeGen ) {
                rg_Point2D sp = edge->getStartVertex()->getLocation();
                rg_Point2D ep = edge->getEndVertex()->getLocation();
                rg_Point2D cp = ( sp + ep ) / 2.0;

                geometry_of_edge.setCtrlPt(0, sp);
                geometry_of_edge.setCtrlPt(1, cp);
                geometry_of_edge.setCtrlPt(2, ep);
                geometry_of_edge.setWeight(0, 1.0);
                geometry_of_edge.setWeight(1, 1.0);
                geometry_of_edge.setWeight(2, 1.0);
            }
            else {
                rg_Point2D focus     = vertexGen->get_point();
                rg_Line2D  directrix = edgeGen->get_geometry();
                Parabola2D parabola(focus, directrix);

                rg_Point2D sp        = edge->getStartVertex()->getLocation();
                rg_Point2D ep        = edge->getEndVertex()->getLocation();
                rg_Point2D tvs       = parabola.get_tangent_vector(sp);
                rg_Point2D tve       = parabola.get_tangent_vector(ep);
                rg_Point2D passingPt = parabola.get_passing_point(sp, ep);

                geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passingPt);

                if ( geometry_of_edge.getWeight(1) < 0 )
                    geometry_of_edge.setWeight(1, geometry_of_edge.getWeight(1) * -1.0);
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            rg_Point2D sp = edge->getStartVertex()->getLocation();
            rg_Point2D ep = edge->getEndVertex()->getLocation();
            rg_Point2D cp = ( sp + ep ) / 2.0;

            geometry_of_edge.setCtrlPt(0, sp);
            geometry_of_edge.setCtrlPt(1, cp);
            geometry_of_edge.setCtrlPt(2, ep);
            geometry_of_edge.setWeight(0, 1.0);
            geometry_of_edge.setWeight(1, 1.0);
            geometry_of_edge.setWeight(2, 1.0);
        }
        break;

        case Generator2D::Generator_Type::DISK_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGenerator;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGenerator;

            rg_Circle2D disk     = diskGen->getDisk();
            rg_Line2D   line     = edgeGen->get_geometry();

            rg_Point2D focus     = disk.getCenterPt();

            rg_Point2D sp = edge->getStartVertex()->getLocation();
            rg_Point2D ep = edge->getEndVertex()->getLocation();

            rg_Line2D directrix;
            {
                double distanceToSP = sp.distance( focus );
                double distanceToEP = ep.distance( focus );

                if ( !rg_EQ( distanceToSP, disk.getRadius() ) ) {
                    rg_Point2D footprint_of_sp;
                    line.compute_perpendicular_footprint_of_point_onto_entire_line( sp, footprint_of_sp );
                    rg_Point2D vec_sp_to_line = ( footprint_of_sp - sp ).getUnitVector();

                    if ( distanceToSP < disk.getRadius() ) {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() - vec_sp_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() - vec_sp_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                    else {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() + vec_sp_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() + vec_sp_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                }
                else if ( !rg_EQ( distanceToEP, disk.getRadius() ) ) {
                    rg_Point2D footprint_of_ep;
                    line.compute_perpendicular_footprint_of_point_onto_entire_line( ep, footprint_of_ep );
                    rg_Point2D vec_ep_to_line = ( footprint_of_ep - ep ).getUnitVector();

                    if ( distanceToEP < disk.getRadius() ) {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() - vec_ep_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() - vec_ep_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                    else {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() + vec_ep_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() + vec_ep_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                }
                else {
                    rg_Point2D normVector = line.getNormalVector().getUnitVector();
                    double signed_distance_to_focus = line.signed_distance( focus );

                    if ( signed_distance_to_focus < 0.0 ) {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() + normVector * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() + normVector * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                    else {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() - normVector * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() - normVector * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                }
            }



            Parabola2D parabola( focus, directrix );
            rg_Point2D tvs       = parabola.get_tangent_vector(sp);
            rg_Point2D tve       = parabola.get_tangent_vector(ep);
            rg_Point2D passingPt = parabola.get_passing_point(sp, ep);

            geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passingPt);

            if ( geometry_of_edge.getWeight(1) < 0 )
                geometry_of_edge.setWeight(1, geometry_of_edge.getWeight(1) * -1.0);
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::DISK_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGenerator;
            DiskGenerator2D*   diskGen   = (DiskGenerator2D*)leftGenerator;

            rg_Point2D  vertex = vertexGen->get_point();
            rg_Circle2D disk   = diskGen->getDisk();

            rg_Point2D sp, ep, tvs, tve, passPoint;

            if( rg_ZERO( disk.getRadius(), resNeg6 ) )
            {
                sp = edge->getStartVertex()->getLocation();
                ep = edge->getEndVertex()->getLocation();
                rg_Point2D cp = ( sp + ep ) / 2.0;

                geometry_of_edge.setCtrlPt(0, sp);
                geometry_of_edge.setCtrlPt(1, cp);
                geometry_of_edge.setCtrlPt(2, ep);
                geometry_of_edge.setWeight(0, 1.0);
                geometry_of_edge.setWeight(1, 1.0);
                geometry_of_edge.setWeight(2, 1.0);
            }
            else
            {
                sp = edge->getStartVertex()->getLocation();
                ep = edge->getEndVertex()->getLocation();
                tvs = getTangentVector( sp, vertex, disk.getCenterPt() );
                tve = getTangentVector( ep, vertex, disk.getCenterPt() );
                passPoint = getPassingPt(rg_Circle2D(vertex, 0.0), disk);

                geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passPoint);
                if( geometry_of_edge.getWeight(1) < 0 )
                    geometry_of_edge.setWeight( 1, geometry_of_edge.getWeight(1) * -1 );
            }
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGenerator;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGenerator;

            rg_Circle2D disk     = diskGen->getDisk();
            rg_Line2D   line     = edgeGen->get_geometry();

            rg_Point2D focus     = disk.getCenterPt();

            rg_Point2D sp = edge->getStartVertex()->getLocation();
            rg_Point2D ep = edge->getEndVertex()->getLocation();
            
            rg_Line2D directrix;
            {
                double distanceToSP = sp.distance( focus );
                double distanceToEP = ep.distance( focus );

                if ( !rg_EQ( distanceToSP, disk.getRadius() ) ) {
                    rg_Point2D footprint_of_sp;
                    line.compute_perpendicular_footprint_of_point_onto_entire_line( sp, footprint_of_sp );
                    rg_Point2D vec_sp_to_line = ( footprint_of_sp - sp ).getUnitVector();

                    if ( distanceToSP < disk.getRadius() ) {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() - vec_sp_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() - vec_sp_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                    else {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() + vec_sp_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() + vec_sp_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                }
                else if ( !rg_EQ( distanceToEP, disk.getRadius() ) ) {
                    rg_Point2D footprint_of_ep;
                    line.compute_perpendicular_footprint_of_point_onto_entire_line( ep, footprint_of_ep );
                    rg_Point2D vec_ep_to_line = ( footprint_of_ep - ep ).getUnitVector();

                    if ( distanceToEP < disk.getRadius() ) {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() - vec_ep_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() - vec_ep_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                    else {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() + vec_ep_to_line * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() + vec_ep_to_line * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                }
                else {
                    rg_Point2D normVector = line.getNormalVector().getUnitVector();
                    double signed_distance_to_focus = line.signed_distance( focus );

                    if ( signed_distance_to_focus < 0.0 ) {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() + normVector * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() + normVector * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                    else {
                        rg_Point2D sp_directrix = edgeGen->get_start_vertex_point() - normVector * disk.getRadius();
                        rg_Point2D ep_directrix = edgeGen->get_end_vertex_point() - normVector * disk.getRadius();
                        directrix = rg_Line2D( sp_directrix, ep_directrix );
                    }
                }
            }
            
            

            Parabola2D parabola(focus, directrix);
            rg_Point2D tvs       = parabola.get_tangent_vector(sp);
            rg_Point2D tve       = parabola.get_tangent_vector(ep);
            rg_Point2D passingPt = parabola.get_passing_point(sp, ep);

            geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passingPt);

            if ( geometry_of_edge.getWeight(1) < 0 )
                geometry_of_edge.setWeight(1, geometry_of_edge.getWeight(1) * -1.0);
        }
        break;

        case Generator2D::Generator_Type::DISK_G:
        {
            DiskGenerator2D* leftDiskGen    = (DiskGenerator2D*)rightGenerator;
            DiskGenerator2D* rightDiskGen   = (DiskGenerator2D*)leftGenerator;

            rg_Circle2D disk1 = leftDiskGen->getDisk();
            rg_Circle2D disk2 = rightDiskGen->getDisk();

            rg_Point2D sp, ep, tvs, tve, passPoint;

            if( rg_EQ( disk1.getRadius(), disk2.getRadius(), resNeg6 ) )
            {
                sp = edge->getStartVertex()->getLocation();
                ep = edge->getEndVertex()->getLocation();
                rg_Point2D cp = ( sp + ep ) / 2.0;

                geometry_of_edge.setCtrlPt(0, sp);
                geometry_of_edge.setCtrlPt(1, cp);
                geometry_of_edge.setCtrlPt(2, ep);
                geometry_of_edge.setWeight(0, 1.0);
                geometry_of_edge.setWeight(1, 1.0);
                geometry_of_edge.setWeight(2, 1.0);
            }
            else
            {
                sp = edge->getStartVertex()->getLocation();
                ep = edge->getEndVertex()->getLocation();
                tvs = getTangentVector( sp, disk1.getCenterPt(), disk2.getCenterPt() );
                tve = getTangentVector( ep, disk1.getCenterPt(), disk2.getCenterPt() );
                passPoint = getPassingPt(disk1, disk2);

                geometry_of_edge.makeRQBezier(sp, tvs, ep, tve, passPoint);
                if( geometry_of_edge.getWeight(1) < 0 )
                    geometry_of_edge.setWeight( 1, geometry_of_edge.getWeight(1) * -1 );


                rg_Point2D controlPt = geometry_of_edge.getCtrlPt(1);
                double x_cp = controlPt.getX();
                double y_cp = controlPt.getY();
                double minValue = -32.0000008867;
                double maxValue = -26.4277301698;
                if ( x_cp < minValue || x_cp > maxValue || y_cp < minValue || y_cp > maxValue ) {
                    rg_Point2D sp_2 = edge->getStartVertex()->getLocation();
                    rg_Point2D ep_2 = edge->getEndVertex()->getLocation();
                    rg_Point2D tvs_2 = getTangentVector( sp, disk1.getCenterPt(), disk2.getCenterPt() );
                    rg_Point2D tve_2 = getTangentVector( ep, disk1.getCenterPt(), disk2.getCenterPt() );
                    rg_Point2D passPoint_2 = getPassingPt(disk1, disk2);

                    rg_RQBzCurve2D geometry_of_edge2;
                    geometry_of_edge2.makeRQBezier(sp, tvs, ep, tve, passPoint);
                    if( geometry_of_edge2.getWeight(1) < 0 )
                        geometry_of_edge2.setWeight( 1, geometry_of_edge2.getWeight(1) * -1 );
                    int stop = 1;
                }
                
            }
        }
        break;

        default:
            break;
        }
    }

    default:
        break;
    }


    return geometry_of_edge;
}



void PolygonVD2D::get_sampling_points_of_edge( const VEdge2D * const edge, const int& numSamplingPts, list<rg_Point2D>& samplingPts )
{
    rg_RQBzCurve2D curve = get_geometry_of_edge(edge);

    int numPtsOnCurve = numSamplingPts;
    if ( numPtsOnCurve < 1 )
        numPtsOnCurve = 2;

    int divisor = numPtsOnCurve - 1;
    for ( int i=0; i<numPtsOnCurve; i++ ) {            					
        double parameter = ((double)i)/divisor;

        rg_Point2D pt = curve.evaluatePt(parameter);
        samplingPts.push_back(pt);				
    }
}



DiskGenerator2D * PolygonVD2D::getGeneratorWhichHasThisID( const int & diskID )
{
    map< int, VFace2D* >::iterator i_map_diskID2VFace;
    i_map_diskID2VFace = m_map_DiskID2VFace.find(diskID);

   Generator2D* generatorWhichHasInputID = NULL;

    bool thereIsDiskGeneratorWhichHasInputID = ( i_map_diskID2VFace != m_map_DiskID2VFace.end() );
    if ( thereIsDiskGeneratorWhichHasInputID ) {
        generatorWhichHasInputID = i_map_diskID2VFace->second->getGenerator();
    }

    return (DiskGenerator2D*)generatorWhichHasInputID;
}



void PolygonVD2D::makeDiskID2VFaceMap()
{
    m_map_DiskID2VFace.clear();
    
    for ( list<Generator2D*>::iterator i_gen = m_generators.begin(); i_gen != m_generators.end(); ++i_gen ) {
        Generator2D*    currGenerator   = *i_gen;
        int             currDiskID      = currGenerator->getID();
        VFace2D*        currVFace       = currGenerator->getOuterFace();

        m_map_DiskID2VFace.insert( make_pair(currDiskID, currVFace) );
    }
}



void PolygonVD2D::remove_red_entities()
{
    removeAllExtraneousVVerticesAndVEdges();
}



void PolygonVD2D::set_type_of_PVertex( VertexGenerator2D * vertexGen, EdgeGenerator2D * prevEdgeGen, EdgeGenerator2D * nextEdgeGen )
{
    rg_Line2D prevEdge;
    prevEdgeGen->get_geometry(prevEdge);
    rg_Line2D nextEdge;
    nextEdgeGen->get_geometry(nextEdge);

    rg_Point2D vecPrevEdge = prevEdge.evaluateVector();
    rg_Point2D vecNextEdge = nextEdge.evaluateVector();

    if (vecPrevEdge.operator*(vecNextEdge) > 0.0)
        vertexGen->set_vertex_type(VertexGenerator2D::Vertex_Type::NON_REFLEX_FROM_POLYGON_INSIDE);
    else
        vertexGen->set_vertex_type(VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE);
}



void PolygonVD2D::find_two_VEdges_between_vertex_generator_N_incident_edge_generators( VertexGenerator2D *& vertexGen, EdgeGenerator2D *& prevEdgeGen, EdgeGenerator2D *& nextEdgeGen, VEdge2D *& edge_between_prevEdgeGen_N_vertexGen, VEdge2D *& edge_between_vertexGen_N_nextEdgeGen )
{
    VFace2D* faceOfVertexGen = rg_NULL;
    VFace2D* faceOfPrevEdgeGen = rg_NULL;
    VFace2D* faceOfNextEdgeGen = rg_NULL;

    // get child VFace of vertex, previous edge, and next edge generator
    faceOfVertexGen   = vertexGen->getOuterFace();
    faceOfPrevEdgeGen = prevEdgeGen->getOuterFace();
    faceOfNextEdgeGen = nextEdgeGen->getOuterFace();


    // find two VEdges ( prevEdgeGen-vertexGen, vertexGen-nextEdgeGen )

    list<VEdge2D*> boundaryEdgesOfVertexGen;
    faceOfVertexGen->getBoundaryVEdges(boundaryEdgesOfVertexGen);

    list<VEdge2D*>::iterator i_edge;
    for (i_edge = boundaryEdgesOfVertexGen.begin(); i_edge != boundaryEdgesOfVertexGen.end(); ++i_edge) {
        VEdge2D* currBoundaryEdge = *i_edge;

        if (currBoundaryEdge->getLeftFace() == faceOfVertexGen) {
            if (currBoundaryEdge->getRightFace() == faceOfPrevEdgeGen) {
                edge_between_prevEdgeGen_N_vertexGen = currBoundaryEdge;
            }
            else if (currBoundaryEdge->getRightFace() == faceOfNextEdgeGen) {
                edge_between_vertexGen_N_nextEdgeGen = currBoundaryEdge;
            }
        }
        else {
            if (currBoundaryEdge->getLeftFace() == faceOfPrevEdgeGen) {
                edge_between_prevEdgeGen_N_vertexGen = currBoundaryEdge;
            }
            else if (currBoundaryEdge->getLeftFace() == faceOfNextEdgeGen) {
                edge_between_vertexGen_N_nextEdgeGen = currBoundaryEdge;
            }
        }
    }
}



void PolygonVD2D::find_and_relocate_two_boundary_Vvertices_to_vertex_generator_center( VertexGenerator2D *& vertexGen, EdgeGenerator2D *& prevEdgeGen, EdgeGenerator2D *& nextEdgeGen )
{
    VFace2D* childVFaceOfVertexGen = rg_NULL;
    VFace2D* childVFaceOfPrevEdgeGen = rg_NULL;
    VFace2D* childVFaceOfNextEdgeGen = rg_NULL;

    // get child VFace of vertex, previous edge, and next edge generator
    childVFaceOfVertexGen   = getGeneratorWhichHasThisID(vertexGen->getID())->getOuterFace();
    childVFaceOfPrevEdgeGen = getGeneratorWhichHasThisID(prevEdgeGen->getID())->getOuterFace();
    childVFaceOfNextEdgeGen = getGeneratorWhichHasThisID(nextEdgeGen->getID())->getOuterFace();


    // find two VEdges ( prevEdgeGen-vertexGen, vertexGen-nextEdgeGen )
    VEdge2D* edge_between_prevEdgeGen_N_vertexGen = rg_NULL;
    VEdge2D* edge_between_vertexGen_N_nextEdgeGen = rg_NULL;

    find_two_VEdges_between_vertex_generator_N_incident_edge_generators(vertexGen, prevEdgeGen, nextEdgeGen, edge_between_prevEdgeGen_N_vertexGen, edge_between_vertexGen_N_nextEdgeGen);


    // check whether nextEdgeGen turns right of left
    //////////////////////////////////////////////////////////////
    //                                                          //
    //   nextQuillEdge        trappedEdges     prevQuillEdge    //
    //          -----------*---*---*---*----*--------------     //
    //                    /                 ^                   //
    //   endBoundaryEdge /                  | startBoundaryEdge //
    //                  /                   |                   //
    //     endFace     v                    |     startFace     //
    //                                                          //
    //////////////////////////////////////////////////////////////

    // collect entities for refine Vvertex
    rg_Line2D   lineOfStartFace;
    rg_Line2D   lineOfEndFace;
    VEdge2D*    startBoundaryEdge   = rg_NULL;
    VEdge2D*    endBoundaryEdge     = rg_NULL;
    VFace2D*    startFace           = rg_NULL;
    VFace2D*    endFace             = rg_NULL;
    VEdge2D*    prevQuillEdge_CCW   = rg_NULL;
    VEdge2D*    nextQuillEdge_CCW   = rg_NULL;
    list<VEdge2D*> trappedEdgesInCCWOrder;

    if (vertexGen->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE) { // nextEdgeGen turns right
                                                                                                  //if ( anglePrevEdgeGenToNextEdgeGen < rg_PI ) { 
        startBoundaryEdge = edge_between_prevEdgeGen_N_vertexGen;
        endBoundaryEdge = edge_between_vertexGen_N_nextEdgeGen;
        startFace = childVFaceOfPrevEdgeGen;
        endFace = childVFaceOfNextEdgeGen;
        prevEdgeGen->get_geometry(lineOfStartFace);
        nextEdgeGen->get_geometry(lineOfEndFace);
    }
    else {
        startBoundaryEdge = edge_between_vertexGen_N_nextEdgeGen;
        endBoundaryEdge = edge_between_prevEdgeGen_N_vertexGen;
        startFace = childVFaceOfNextEdgeGen;
        endFace = childVFaceOfPrevEdgeGen;
        nextEdgeGen->get_geometry(lineOfStartFace);
        prevEdgeGen->get_geometry(lineOfEndFace);
    }


    VEdge2D* currTrappedEdge = rg_NULL;
    if (startBoundaryEdge->getLeftFace() == childVFaceOfVertexGen) {
        currTrappedEdge = startBoundaryEdge->getLeftHand();
        prevQuillEdge_CCW = startBoundaryEdge->getRightHand();
    }
    else {
        currTrappedEdge = startBoundaryEdge->getRightLeg();
        prevQuillEdge_CCW = startBoundaryEdge->getLeftLeg();
    }

    if (endBoundaryEdge->getLeftFace() == childVFaceOfVertexGen) {
        nextQuillEdge_CCW = endBoundaryEdge->getRightLeg();
    }
    else {
        nextQuillEdge_CCW = endBoundaryEdge->getLeftHand();
    }

    while (currTrappedEdge != endBoundaryEdge) {
        trappedEdgesInCCWOrder.push_back(currTrappedEdge);

        if (currTrappedEdge->getLeftFace() == childVFaceOfVertexGen) {
            currTrappedEdge = currTrappedEdge->getLeftHand();
        }
        else {
            currTrappedEdge = currTrappedEdge->getRightLeg();
        }
    }


    // find edge to be split
    VEdge2D* edgeToBeSplit = rg_NULL;
    list<VEdge2D*>::iterator i_trappedEdge;
    for (i_trappedEdge = trappedEdgesInCCWOrder.begin(); i_trappedEdge != trappedEdgesInCCWOrder.end(); ++i_trappedEdge) {
        VEdge2D*   currEdge = *i_trappedEdge;
        rg_Point2D startPt = currEdge->getStartVertex()->getLocation();
        rg_Point2D endPt = currEdge->getEndVertex()->getLocation();

        rg_Point2D footprintOnLineOfStartFace;
        rg_Point2D foorprintOnLineOfEndFace;

        // check whether start and edge points are closer to same line segment
        bool    startPtIsCloserToLineOfStartFace = false;
        bool    endPtIsCloserToLineOfStartFace = false;

        // check start point
        lineOfStartFace.compute_footprint_of_point_onto_line_segment(startPt, footprintOnLineOfStartFace);
        lineOfEndFace.compute_footprint_of_point_onto_line_segment(startPt, foorprintOnLineOfEndFace);

        double distanceToStartLine = startPt.distance(footprintOnLineOfStartFace);
        double distanceToEndLine = startPt.distance(foorprintOnLineOfEndFace);

        if (distanceToStartLine < distanceToEndLine) {
            startPtIsCloserToLineOfStartFace = true;
        }

        // check end point
        lineOfStartFace.compute_footprint_of_point_onto_line_segment(endPt, footprintOnLineOfStartFace);
        lineOfEndFace.compute_footprint_of_point_onto_line_segment(endPt, foorprintOnLineOfEndFace);

        distanceToStartLine = endPt.distance(footprintOnLineOfStartFace);
        distanceToEndLine = endPt.distance(foorprintOnLineOfEndFace);

        if (distanceToStartLine < distanceToEndLine) {
            endPtIsCloserToLineOfStartFace = true;
        }



        if (startPtIsCloserToLineOfStartFace != endPtIsCloserToLineOfStartFace) {
            edgeToBeSplit = currEdge;
            break;
        }

    }

    // flip trappedEdges in CCW order before edgeToBeSplit
    for (i_trappedEdge = trappedEdgesInCCWOrder.begin(); i_trappedEdge != trappedEdgesInCCWOrder.end(); ++i_trappedEdge) {
        VEdge2D* currEdge = *i_trappedEdge;

        if (currEdge == edgeToBeSplit) {
            break;
        }

        currEdge->flip();
    }

    // flip trappedEdges in CW order before edgeToBeSplit
    i_trappedEdge = trappedEdgesInCCWOrder.end();
    --i_trappedEdge;
    for (; i_trappedEdge != trappedEdgesInCCWOrder.begin(); --i_trappedEdge) {
        VEdge2D* currEdge = *i_trappedEdge;

        if (currEdge == edgeToBeSplit) {
            break;
        }

        currEdge->flip();
    }

    edgeToBeSplit->flip();
    if (edgeToBeSplit->getMateFace(edgeToBeSplit->getStartVertex()) == childVFaceOfVertexGen) {
        edgeToBeSplit->getEndVertex()->setCircumcircle(rg_Circle2D(vertexGen->get_point(), 0.0));
    }
    else {
        edgeToBeSplit->getStartVertex()->setCircumcircle(rg_Circle2D(vertexGen->get_point(), 0.0));
    }
}



void PolygonVD2D::identify_internal_and_on_and_external_Vvertices_and_Vedges()
{
    identify_internal_and_on_and_external_and_infinite_Vvertices_and_Vedges();
    /*
    classifyVertexGenerators();

    list< pair< VertexGenerator2D*, VVertex2D* > > PVertex_N_onVVertexPair;
    identify_Vvertices_on_boundary( PVertex_N_onVVertexPair );

    // infinite
    stack<VEdge2D*> outerEdges;
    stack<VEdge2D*> innerEdges;

    // collect initial stacks of outer edges and inner edges.
    for ( list< pair< VertexGenerator2D*, VVertex2D* > >::iterator i_pair = PVertex_N_onVVertexPair.begin(); i_pair != PVertex_N_onVVertexPair.end(); ++i_pair ) {
        VertexGenerator2D* currPVertex = i_pair->first;
        VVertex2D* currOnVVertex = i_pair->second;

        EdgeGenerator2D* prevPEdge = currPVertex->get_previous_edge_generator();
        EdgeGenerator2D* nextPEdge = currPVertex->get_next_edge_generator();

        if ( prevPEdge == NULL || nextPEdge == NULL ) {
            continue;
        }

        VertexGenerator2D::Vertex_Type vertexType = currPVertex->vertex_type();

        // find initial stack of edge
        list<VEdge2D*> incidentEdges;
        currOnVVertex->getIncident3VEdges( incidentEdges );
        for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
            VEdge2D* currEdge = *i_edge;

            bool b_this_edge_is_defined_by_two_PEdges = true;
            if ( currEdge->getLeftFace()->getGenerator() == currPVertex || currEdge->getRightFace()->getGenerator() == currPVertex ) {
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
        VVertex_Location_Status_Map::iterator i_startVertexToLocationStatus = m_VVertexToLocationStatus.find( startVertex );
        VVertex_Location_Status_Map::iterator i_endVertexToLocationStatus = m_VVertexToLocationStatus.find( endVertex );

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
    */
}



void V::GeometryTier::PolygonVD2D::identify_internal_and_on_and_external_and_infinite_Vvertices_and_Vedges()
{
    classifyVertexGenerators();

    list< pair< VertexGenerator2D*, VVertex2D* > > PVertex_N_onVVertexPair;
    identify_Vvertices_on_boundary( PVertex_N_onVVertexPair );

    identify_infinite_Vvertices_and_Vedges();

    identify_outside_Vvertices_and_Vedges( PVertex_N_onVVertexPair );
    
    set_location_status_of_remaining_Vvertices_and_Vedges_with_inside();
}



void V::GeometryTier::PolygonVD2D::updateEdge(VEdge2D* target)
{
    rg_RQBzCurve2D geometry = get_geometry_of_edge( target );
    target->setGeometry(geometry);
}



void PolygonVD2D::set_location_status_of_infinite_Vvertices_and_Vedges_with_outside( stack<VEdge2D*>& VEdgeStack )
{
    VFace2D* infiniteFace = NULL;

    for ( list<VFace2D*>::iterator i_face = m_VFaces.begin(); i_face != m_VFaces.end(); ++i_face ) {
        VFace2D* currFace = *i_face;
        if ( currFace->isInfinite() ) {
            infiniteFace = currFace;
            break;
        }
    }

    // mark 'outside' tag for infinite VVertex
    list<VVertex2D*> infiniteVVertices;
    infiniteFace->getBoundaryVVertices(infiniteVVertices);
    for ( list<VVertex2D*>::iterator i_infVtx = infiniteVVertices.begin(); i_infVtx != infiniteVVertices.end(); ++i_infVtx ) {
        VVertex2D* currInfVtx = *i_infVtx;
        m_VVertexToLocationStatus.insert(pair<VVertex2D*, VDEntity_Location_Status>(currInfVtx, OUTSIDE_POLYGON));
    }


    // mark 'outside' tag for infinite VEdges
    list<VEdge2D*> infiniteVEdges;
    infiniteFace->getBoundaryVEdges(infiniteVEdges);
    for ( list<VEdge2D*>::iterator i_infEdge = infiniteVEdges.begin(); i_infEdge != infiniteVEdges.end(); ++i_infEdge ) {
        VEdge2D* currInfEdge = *i_infEdge;
        m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(currInfEdge, OUTSIDE_POLYGON));
    }

    // get_geometry non-infinite VEdges incident to infinite VVertices	
    for ( list<VVertex2D*>::iterator i_infVtx = infiniteVVertices.begin(); i_infVtx != infiniteVVertices.end(); ++i_infVtx ) {
        VVertex2D* currInfVtx = *i_infVtx;

        list<VEdge2D*> incidentVEdges;
        currInfVtx->getIncident3VEdges(incidentVEdges);
        for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
            VEdge2D* currEdge = *i_edge;

            if ( !currEdge->isInfinite() ) {
                VEdgeStack.push(currEdge);
                break;
            }
        }
    }
}



void PolygonVD2D::identify_Vvertices_on_boundary( list<pair<VertexGenerator2D*, VVertex2D*>>& PVertex_N_onVVertexPair )
{
    for ( list<Generator2D*>::iterator i_gen = m_generators.begin(); i_gen != m_generators.end(); ++i_gen ) {
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
        VVertex2D* endVertex_of_VEdge_between_prevEdgeGen_N_vertexGen   = VEdge_between_prevEdgeGen_N_vertexGen->getEndVertex();
        VVertex2D* startVertex_of_VEdge_between_vertexGen_N_nextEdgeGen = VEdge_between_vertexGen_N_nextEdgeGen->getStartVertex();
        VVertex2D* endVertex_of_VEdge_between_vertexGen_N_nextEdgeGen   = VEdge_between_vertexGen_N_nextEdgeGen->getEndVertex();

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



void PolygonVD2D::identify_infinite_Vvertices_and_Vedges()
{
    VFace2D* infiniteFace = NULL;

    // There is one infinite VFace in the list of VFaces.
    for ( list<VFace2D*>::iterator i_face = m_VFaces.begin(); i_face != m_VFaces.end(); ++i_face ) {
        VFace2D* currFace = *i_face;
        if ( currFace->isInfinite() ) {
            infiniteFace = currFace;
            break;
        }
    }

    stack<VEdge2D*> infiniteEdges;

    // collect initial stacks of infinite edges.
    list<VEdge2D*> boundaryEdges;
    infiniteFace->getBoundaryVEdges( boundaryEdges );
    for ( list<VEdge2D*>::iterator i_edge = boundaryEdges.begin(); i_edge != boundaryEdges.end(); ++i_edge ) {
        VEdge2D* currEdge = *i_edge;
        infiniteEdges.push( currEdge );
    }


    // outer edges
    while ( !infiniteEdges.empty() )
    {
        VEdge2D* currInfEdge = infiniteEdges.top();
        infiniteEdges.pop();

        if ( m_VEdgeToLocationStatus.find( currInfEdge ) != m_VEdgeToLocationStatus.end() ) {
            continue;
        }

        m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( currInfEdge, INFINITE_ENTITY ) );

        VVertex2D* startVertex = currInfEdge->getStartVertex();
        VVertex2D* endVertex = currInfEdge->getEndVertex();
        VVertex_Location_Status_Map::iterator i_startVertexToLocationStatus = m_VVertexToLocationStatus.find( startVertex );
        VVertex_Location_Status_Map::iterator i_endVertexToLocationStatus = m_VVertexToLocationStatus.find( endVertex );

        if ( i_startVertexToLocationStatus == m_VVertexToLocationStatus.end() ) {
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( startVertex, INFINITE_ENTITY ) );

            list<VEdge2D*> incidentVEdges;
            startVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    infiniteEdges.push( currEdge );
            }
        }

        if ( i_endVertexToLocationStatus == m_VVertexToLocationStatus.end() ) {
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( endVertex, INFINITE_ENTITY ) );

            list<VEdge2D*> incidentVEdges;
            endVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    infiniteEdges.push( currEdge );
            }
        }
    }
}



void V::GeometryTier::PolygonVD2D::identify_outside_Vvertices_and_Vedges( const list< pair< VertexGenerator2D*, VVertex2D* > >& PVertex_N_onVVertexPair )
{
    stack<VEdge2D*> outerEdges;
    //stack<VEdge2D*> innerEdges;

    // collect initial stacks of outer edges.
    for ( list< pair< VertexGenerator2D*, VVertex2D* > >::const_iterator i_pair = PVertex_N_onVVertexPair.begin(); i_pair != PVertex_N_onVVertexPair.end(); ++i_pair ) {
        VertexGenerator2D* currPVertex = i_pair->first;
        VVertex2D* currOnVVertex = i_pair->second;

        EdgeGenerator2D* prevPEdge = currPVertex->get_previous_edge_generator();
        EdgeGenerator2D* nextPEdge = currPVertex->get_next_edge_generator();

        if ( prevPEdge == NULL || nextPEdge == NULL ) {
            continue;
        }

        VertexGenerator2D::Vertex_Type vertexType = currPVertex->vertex_type();

        // find initial stack of edge
        list<VEdge2D*> incidentEdges;
        currOnVVertex->getIncident3VEdges( incidentEdges );
        for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
            VEdge2D* currEdge = *i_edge;

            if ( m_VEdgeToLocationStatus.find( currEdge ) != m_VEdgeToLocationStatus.end() ) {
                continue;
            }

            bool b_this_edge_is_defined_by_two_PEdges = true;
            if ( currEdge->getLeftFace()->getGenerator() == currPVertex || currEdge->getRightFace()->getGenerator() == currPVertex ) {
                b_this_edge_is_defined_by_two_PEdges = false;
            }
            else {
                b_this_edge_is_defined_by_two_PEdges = true;
            }


            if ( vertexType == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
                if ( b_this_edge_is_defined_by_two_PEdges ) {
                    outerEdges.push( currEdge );
                }
                //else {
                //    innerEdges.push( currEdge );
                //}
            }
            else {
                if ( !b_this_edge_is_defined_by_two_PEdges ) {
                    outerEdges.push( currEdge );
                }
                //else {
                //    innerEdges.push( currEdge );
                //}
            }
        }
    }

    // outer edges
    while ( !outerEdges.empty() )
    {
        VEdge2D* outerEdge = outerEdges.top();
        outerEdges.pop();

        m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

        VVertex2D* startVertex = outerEdge->getStartVertex();
        VVertex2D* endVertex = outerEdge->getEndVertex();
        VVertex_Location_Status_Map::iterator i_startVertexToLocationStatus = m_VVertexToLocationStatus.find( startVertex );
        VVertex_Location_Status_Map::iterator i_endVertexToLocationStatus = m_VVertexToLocationStatus.find( endVertex );

        if ( i_startVertexToLocationStatus == m_VVertexToLocationStatus.end() ) {
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( startVertex, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            startVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }

        if ( i_endVertexToLocationStatus == m_VVertexToLocationStatus.end() ) {
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( endVertex, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            endVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }

        /*
        // If start vertex was visited before
        if ( i_startVertexToLocationStatus != m_VVertexToLocationStatus.end() )
        {
            if ( i_endVertexToLocationStatus != m_VVertexToLocationStatus.end() )
            {
                //m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );
                continue;
            }
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( endVertex, OUTSIDE_POLYGON ) );
            //m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

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
                //m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );
                continue;
            }
            m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( startVertex, OUTSIDE_POLYGON ) );
            //m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( outerEdge, OUTSIDE_POLYGON ) );

            list<VEdge2D*> incidentVEdges;
            startVertex->getIncident3VEdges( incidentVEdges );
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find( currEdge ) == m_VEdgeToLocationStatus.end() )
                    outerEdges.push( currEdge );
            }
        }
        */
    }
}



void PolygonVD2D::identify_and_set_location_status_of_both_outside_and_on_Vvertices_and_outside_Vedges(stack<VEdge2D*>& VEdgeStack)
{
    while (!VEdgeStack.empty())
    {
        VEdge2D* currentVEdge = VEdgeStack.top();
        VEdgeStack.pop();

        VVertex2D* startVertex = currentVEdge->getStartVertex();
        VVertex2D* endVertex   = currentVEdge->getEndVertex();
        VVertex_Location_Status_Map::iterator i_startVertexToLocationStatus = m_VVertexToLocationStatus.find(startVertex);
        VVertex_Location_Status_Map::iterator i_endVertexToLocationStatus   = m_VVertexToLocationStatus.find(endVertex);

        // If start vertex was visited before
        if (i_startVertexToLocationStatus != m_VVertexToLocationStatus.end())
        {
            if (i_endVertexToLocationStatus != m_VVertexToLocationStatus.end())
            {
                m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(currentVEdge, OUTSIDE_POLYGON));
                continue;
            }
            m_VVertexToLocationStatus.insert(pair<VVertex2D*, VDEntity_Location_Status>(endVertex, OUTSIDE_POLYGON));
            m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(currentVEdge, OUTSIDE_POLYGON));

            list<VEdge2D*> incidentVEdges;
            endVertex->getIncident3VEdges(incidentVEdges);
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find(currEdge) == m_VEdgeToLocationStatus.end() )
                    VEdgeStack.push(currEdge);
            }
        }

        // If end vertex was visited before
        if (i_endVertexToLocationStatus != m_VVertexToLocationStatus.end())
        {
            if (i_startVertexToLocationStatus != m_VVertexToLocationStatus.end())
            {
                m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(currentVEdge, OUTSIDE_POLYGON));
                continue;
            }
            m_VVertexToLocationStatus.insert(pair<VVertex2D*, VDEntity_Location_Status>(startVertex, OUTSIDE_POLYGON));
            m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(currentVEdge, OUTSIDE_POLYGON));

            list<VEdge2D*> incidentVEdges;
            startVertex->getIncident3VEdges(incidentVEdges);
            for ( list<VEdge2D*>::iterator i_edge = incidentVEdges.begin(); i_edge != incidentVEdges.end(); ++i_edge ) {
                VEdge2D* currEdge = *i_edge;
                if ( m_VEdgeToLocationStatus.find(currEdge) == m_VEdgeToLocationStatus.end() )
                    VEdgeStack.push(currEdge);
            }
        }
    }
}



void PolygonVD2D::set_location_status_of_remaining_Vvertices_and_Vedges_with_inside()
{
    // We should collect new VVertices and new VEdges to accelerate
    for ( list<VVertex2D*>::iterator i_vtx = m_VVertices.begin(); i_vtx != m_VVertices.end(); ++i_vtx ) {
        VVertex2D* vertex = *i_vtx;
        if (m_VVertexToLocationStatus.find(vertex) == m_VVertexToLocationStatus.end())
            m_VVertexToLocationStatus.insert(pair<VVertex2D*, VDEntity_Location_Status>(vertex, INSIDE_POLYGON));
    }

    for ( list<VEdge2D*>::iterator i_edge = m_VEdges.begin(); i_edge != m_VEdges.end(); ++i_edge ) {
        VEdge2D* edge = *i_edge;
        if (m_VEdgeToLocationStatus.find(edge) == m_VEdgeToLocationStatus.end())
            m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(edge, INSIDE_POLYGON));
    }
}



void PolygonVD2D::refine_location_of_interior_VVertices()
{
    for ( list<VVertex2D*>::iterator i_vtx = m_VVertices.begin(); i_vtx != m_VVertices.end(); ++i_vtx ) {
        VVertex2D* vertex = *i_vtx;

        if ( vertex->isInfinite() )
            continue;

        if ( get_location_status_of_Vvertex(vertex) == ON_POLYGON_BOUNDARY )
            continue;

        bool tangentCircleLocatedPolygonOutward = false;
        rg_Circle2D correctedEmptyTangnetCircle = refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle(vertex);
        vertex->setCircumcircle(correctedEmptyTangnetCircle);
    }
}



PolygonVD2D::VDEntity_Location_Status PolygonVD2D::get_location_status_of_Vvertex( const VVertex2D * const vertex )
{
    VVertex_Location_Status_Map::iterator i_VVertexToLocationStatus = m_VVertexToLocationStatus.find(const_cast<VVertex2D*>(vertex));
    if (i_VVertexToLocationStatus != m_VVertexToLocationStatus.end())
        return i_VVertexToLocationStatus->second;
    else
        return VDEntity_Location_Status::UNKNOWN_LOCATION;
}



PolygonVD2D::VDEntity_Location_Status PolygonVD2D::get_location_status_of_Vedge( const VEdge2D * const edge )
{
    VEdge_Location_Status_Map::iterator i_VEdgeToLocationStatus = m_VEdgeToLocationStatus.find(const_cast<VEdge2D*>(edge));
    if (i_VEdgeToLocationStatus != m_VEdgeToLocationStatus.end())
        return i_VEdgeToLocationStatus->second;
    else
        return VDEntity_Location_Status::UNKNOWN_LOCATION;
}



rg_Circle2D PolygonVD2D::refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle( const VVertex2D * const vertex )
{
    rg_INT numDisks, numLines, numPoints;
    Generator2D* parentGens[3];
    // parent generators in CCW orientation
    count_number_of_occurred_generator_types_and_get_parent_generators_at_VVertex(vertex, numDisks, numLines, numPoints, parentGens);

    rg_Circle2D refinedTangentCircle;

    switch (numDisks)
    {
    case 3: // numDisks == 3, numLines == 0, numPoints == 0
    {
        rg_Circle2D disks[3];
        disks[0] = ((DiskGenerator2D*)parentGens[0])->getDisk();
        disks[1] = ((DiskGenerator2D*)parentGens[1])->getDisk();
        disks[2] = ((DiskGenerator2D*)parentGens[2])->getDisk();
        rg_Circle2D nullCircle;
        rg_Circle2D::makeCircumcircle( disks[0], disks[1], disks[2], refinedTangentCircle, nullCircle );
    }
    break;

    case 2:
    {
        // numDisks == 2, numLines == 1, numPoints == 0
        if (numLines == 1)
        {
            rg_Circle2D disks[2];
            rg_Line2D lineSeg3;
            rg_INT diskIndex = 0;
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                    disks[diskIndex++] = ((DiskGenerator2D*)parentGens[i])->getDisk();
                if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                    ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg3);
            }

            Ellipse2D ellipse[2] = { Ellipse2D(disks[0].getCenterPt(), disks[0].getRadius(), disks[0].getRadius(), 0.0), Ellipse2D(disks[1].getCenterPt(), disks[1].getRadius(), disks[1].getRadius(), 0.0) };
            refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(ellipse[0], ellipse[1], lineSeg3, vertex->getCircumcircle());
        }
        // numDisks == 2, numLines == 0, numPoints == 1
        else
        {
            rg_Circle2D disks[2];
            rg_Point2D point3;
            rg_INT diskIndex = 0;
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                    disks[diskIndex++] = ((DiskGenerator2D*)parentGens[i])->getDisk();
                if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                    ((VertexGenerator2D*)parentGens[i])->get_geometry(point3);
            }
            rg_Circle2D nullCircle;
            rg_Circle2D::makeCircumcircle( disks[0], disks[1], rg_Circle2D(point3, 0.0), refinedTangentCircle, nullCircle );
        }
    }
    break;
    
    case 1:
    {
        switch (numLines)
        {
        case 2: // numDisks == 1, numLines == 2, numPoints == 0
        {
            rg_Circle2D disk1;
            rg_Line2D lineSeg[2];
            rg_INT lineSegIndex = 0;
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                    disk1 = ((DiskGenerator2D*)parentGens[i])->getDisk();
                if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                    ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg[lineSegIndex++]);
            }

            refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(Ellipse2D(disk1.getCenterPt(), disk1.getRadius(), disk1.getRadius(), 0.0), lineSeg[0], lineSeg[1], vertex->getCircumcircle());
        }
        break;

        case 1: // numDisks == 1, numLines == 1, numPoints == 1
        {
            EdgeGenerator2D*     edgeGen = rg_NULL;
            VertexGenerator2D* vertexGen = rg_NULL;

            rg_Circle2D disk1;
            rg_Line2D  lineSeg2;
            rg_Point2D point3;
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType()) {
                    disk1 = ((DiskGenerator2D*)parentGens[i])->getDisk();
                }
                if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                {
                    edgeGen = (EdgeGenerator2D*)parentGens[i];
                    ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg2);
                }
                if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                {
                    vertexGen = (VertexGenerator2D*)parentGens[i];
                    ((VertexGenerator2D*)parentGens[i])->get_geometry(point3);
                }
            }

            if (   vertexGen->get_next_edge_generator() != edgeGen && vertexGen->get_previous_edge_generator() != edgeGen ) {
                refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(Ellipse2D(disk1.getCenterPt(), disk1.getRadius(), disk1.getRadius(), 0.0), lineSeg2, point3, vertex->getCircumcircle());
            }
            else
            {
                // search VVertex iteratively on the bisector line between reflex vertex and line segment
                rg_Line2D perpendicularBisectorBetweenReflexVertexGenAndEdgeGen(point3, point3 + lineSeg2.getNormalVector());

                // VVetex is little bit away from perpendicular bisector between reflex vertex and edge (line segment) because this VVertex is not yet refined					
                rg_Point2D seedPoint = vertex->getLocation();
                // Therefore, we project seedPoint onto the perpendicular bisector and then go into the search process
                rg_Point2D footPrintOfSeedPointOntoPerpendicularBisector;
                perpendicularBisectorBetweenReflexVertexGenAndEdgeGen.compute_perpendicular_footprint_of_point_onto_entire_line(seedPoint, footPrintOfSeedPointOntoPerpendicularBisector);

                rg_Point2D refinedVVertexLocation = rg_GeoFunc::compute_a_point_on_line_whose_distance_with_anchor_point_on_line_is_equal_to_distance_between_the_point_and_its_footprint_of_ellipse(footPrintOfSeedPointOntoPerpendicularBisector, perpendicularBisectorBetweenReflexVertexGenAndEdgeGen, point3, Ellipse2D(disk1.getCenterPt(), disk1.getRadius(), disk1.getRadius(), 0.0), 10e-3, 10);
                refinedTangentCircle.setCircle(refinedVVertexLocation, refinedVVertexLocation.distance(point3));

                //return seedVVertexOfFirstSonVD->getCircumcircle();
            }
        }
        break;
        case 0: // numDisks == 1, numLines == 0, numPoints == 2
        {
            rg_Circle2D disk1;
            rg_Point2D point[2];
            rg_INT pointIndex = 0;
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                    disk1 = ((DiskGenerator2D*)parentGens[i])->getDisk();
                if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                    ((VertexGenerator2D*)parentGens[i])->get_geometry(point[pointIndex++]);
            }
            rg_Circle2D nullCircle;
            rg_Circle2D::makeCircumcircle(  disk1, rg_Circle2D(point[0], 0.0), rg_Circle2D(point[1], 0.0), refinedTangentCircle, nullCircle );
        }
        break;
        }
    }
    break;
    case 0:
    {
        switch (numLines)
        {
        case 3: // numDisks == 0, numLines == 3, numPoints == 0
        {
            rg_Line2D lineSeg[3];
            ((EdgeGenerator2D*)parentGens[0])->get_geometry(lineSeg[0]);
            ((EdgeGenerator2D*)parentGens[1])->get_geometry(lineSeg[1]);
            ((EdgeGenerator2D*)parentGens[2])->get_geometry(lineSeg[2]);
            refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(lineSeg[0], lineSeg[1], lineSeg[2], vertex->getCircumcircle());
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
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                {
                    edgeGen[lineSegindex] = (EdgeGenerator2D*)parentGens[i];
                    ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg[lineSegindex++]);
                }
                if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                {
                    ((VertexGenerator2D*)parentGens[i])->get_geometry(point3);
                    vertexGen = (VertexGenerator2D*)parentGens[i];
                }
            }

            bool bReflexVertexBoundsEdge = true;
            rg_Line2D perpendicularBisector;
            if (vertexGen->get_next_edge_generator() == edgeGen[0] || vertexGen->get_previous_edge_generator() == edgeGen[0])  
            {
                perpendicularBisector.setSP(point3);
                perpendicularBisector.setEP(point3 + lineSeg[0].getNormalVector());
            }
            else if (vertexGen->get_next_edge_generator() == edgeGen[1] || vertexGen->get_previous_edge_generator() == edgeGen[1])
            {
                perpendicularBisector.setSP(point3);
                perpendicularBisector.setEP(point3 + lineSeg[1].getNormalVector());
            }
            else
            {
                bReflexVertexBoundsEdge = false;
            }

            // this code should be rewritten as follows:
            // precise method: compute intersection between both "bisector line" of two lineSeg generators 
            // and "bisector parabola" of lineSeg and point generators

            // CASE II
            if (bReflexVertexBoundsEdge)
            {
                rg_Line2D bisectorLine = rg_GeoFunc::compute_bisector_line_between_two_line_segments(lineSeg[0].get_reversed_line2D(), lineSeg[1]);
                bool bTwoLinesAreParallel = false;
                rg_Point2D refinedCenter = bisectorLine.compute_intersection_with_line(perpendicularBisector, bTwoLinesAreParallel);
                rg_REAL    refinedRadius = refinedCenter.distance(point3);
                refinedTangentCircle.setCircle(refinedCenter, refinedRadius);
            }
            // CASE I
            else
                refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(lineSeg[0], lineSeg[1], point3, vertex->getCircumcircle());
        }
        break;
        case 1: // numDisks == 0, numLines == 1, numPoints == 2
        {
            // THERE ARE THREE CASES 
            // CASE I: ADD THE CASE THAT LINE SEPERATES TWO POINTS: INFEASIBLE CASE

            EdgeGenerator2D*   edgeGen = rg_NULL;
            VertexGenerator2D* vertexGen[2] = { rg_NULL, rg_NULL };

            rg_Line2D  lineSeg;
            rg_INT pointIndex = 0;
            rg_Point2D point[2];
            for (rg_INT i = 0; i < 3; i++)
            {
                if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                {
                    edgeGen = (EdgeGenerator2D*)parentGens[i];
                    ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg);
                }
                if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                {
                    vertexGen[pointIndex] = (VertexGenerator2D*)parentGens[i];
                    ((VertexGenerator2D*)parentGens[i])->get_geometry(point[pointIndex++]);
                }
            }

            bool bfirstReflexVertexBoundsEdge = false;
            bool bSecondReflexVertexBoundsEdge = false;
            rg_Line2D perpendicularBisectorBetweenReflexVertexAndEdge;
            if (vertexGen[0]->get_next_edge_generator() == edgeGen || vertexGen[0]->get_previous_edge_generator() == edgeGen) 
            {
                perpendicularBisectorBetweenReflexVertexAndEdge.setSP(point[0]);
                perpendicularBisectorBetweenReflexVertexAndEdge.setEP(point[0] + lineSeg.getNormalVector());
                bfirstReflexVertexBoundsEdge = true;
            }
            if (vertexGen[1]->get_next_edge_generator() == edgeGen || vertexGen[1]->get_previous_edge_generator() == edgeGen)
            {
                perpendicularBisectorBetweenReflexVertexAndEdge.setSP(point[1]);
                perpendicularBisectorBetweenReflexVertexAndEdge.setEP(point[1] + lineSeg.getNormalVector());
                bSecondReflexVertexBoundsEdge = true;
            }

            // CASE II: ONE VERTEX BOUNDS EDGE(LINESEGMENT) ON THE EDGE
            if (bfirstReflexVertexBoundsEdge || bSecondReflexVertexBoundsEdge)
            {
                rg_Line2D perpendicularBisectorBetweenTwoReflexVertices = rg_GeoFunc::compute_bisector_line_between_two_points(point[0], point[1]);
                bool bTwoLinesAreParllel = false;
                rg_Point2D refinedCenter = perpendicularBisectorBetweenReflexVertexAndEdge.compute_intersection_with_line(perpendicularBisectorBetweenTwoReflexVertices, bTwoLinesAreParllel);

                rg_Point2D edgeDirVec = lineSeg.evaluateVector();
                rg_Point2D perpendicularBisectorTowardPolygonInside = refinedCenter - perpendicularBisectorBetweenReflexVertexAndEdge.getSP();
                if (edgeDirVec.operator*(perpendicularBisectorTowardPolygonInside) > 0.0)
                {
                    rg_REAL    refinedRadius = refinedCenter.distance(point[0]);
                    refinedTangentCircle.setCircle(refinedCenter, refinedRadius);
                }
                else
                {
                    // MAYBE THIS CASE WILL NOT BE ENCOUNTERED ... REMOVED THIS BLOCK
                    //AfxMessageBox(_T("Orthogonalize VEdge by projection"));

                    if(bfirstReflexVertexBoundsEdge)
                        orthogonalize_VEdge_by_projecting_VVertex(const_cast<VVertex2D*>(vertex), edgeGen, point[0]);
                    else
                        orthogonalize_VEdge_by_projecting_VVertex(const_cast<VVertex2D*>(vertex), edgeGen, point[1]);

                    refinedTangentCircle = vertex->getCircumcircle();
                }
            }
            // CASE III: TWO VERTICES ARE OFF THE EDGE(LINESEGMENT) IN THE SAME HALFSPACE
            else
            {
                //refinedTangentCircle = rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(lineSeg, point[0], point[1], vertex->getCircumcircle());
                int numTangentCircles = 0;
                rg_Circle2D tangentCircle[2];
                numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( rg_Circle2D( point[0], 0.0 ), rg_Circle2D( point[1], 0.0 ), lineSeg, tangentCircle[0], tangentCircle[1] );

                if ( numTangentCircles == 1 ) {
                    refinedTangentCircle = tangentCircle[0];
                }
                else {
                    if ( tangentCircle[0].getRadius() > tangentCircle[1].getRadius() ) {
                        rg_Circle2D tempCircle = tangentCircle[0];
                        tangentCircle[0] = tangentCircle[1];
                        tangentCircle[1] = tempCircle;
                    }

                    if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                        refinedTangentCircle = tangentCircle[0];
                    }
                    else {
                        refinedTangentCircle = tangentCircle[1];
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
            ((VertexGenerator2D*)parentGens[0])->get_geometry(point[0]);
            ((VertexGenerator2D*)parentGens[1])->get_geometry(point[1]);
            ((VertexGenerator2D*)parentGens[2])->get_geometry(point[2]);
            rg_Circle2D nullCircle;
            rg_Circle2D::makeCircumcircle( rg_Circle2D(point[0], 0.0), rg_Circle2D(point[1], 0.0), rg_Circle2D(point[2], 0.0), refinedTangentCircle, nullCircle );
        }
        break;
        }
    }
    break;
    }
    return refinedTangentCircle;
}



void PolygonVD2D::count_number_of_occurred_generator_types_and_get_parent_generators_at_VVertex( const VVertex2D * const vertex, int & numberOfDisks, int & numberOfLineSegs, int & numberOfPoints, Generator2D * threeParentGenerators[] )
{
    numberOfDisks = numberOfLineSegs = numberOfPoints = 0;
    list<Generator2D*> diskGens;
    vertex->getDefining3Generators(diskGens);

    int index_gen = 0;
    for ( list<Generator2D*>::iterator i_gen = diskGens.begin(); i_gen != diskGens.end(); ++i_gen ) {
        Generator2D* parentGen = *i_gen;
        threeParentGenerators[index_gen++] = parentGen;
    }

    for ( int i = 0; i < 3; ++i ) {
        Generator2D::Generator_Type type = threeParentGenerators[i]->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::DISK_G:
        {
            ++numberOfDisks;
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            ++numberOfLineSegs;
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            ++numberOfPoints;
        }
        break;

        default:
            break;
        }
    }
}



void PolygonVD2D::orthogonalize_VEdge_by_projecting_VVertex( VVertex2D * vertex, EdgeGenerator2D * edgeGen, const rg_Point2D & passingPt )
{
    rg_Line2D PEdgeline;
    edgeGen->get_geometry(PEdgeline);
    rg_Line2D perpendicularLine = PEdgeline.make_perpendicular_line(passingPt);
    rg_REAL parameter;
    rg_Point2D projectedLocation = perpendicularLine.project(vertex->getLocation(), parameter);
    vertex->setLocation(projectedLocation);
}



void PolygonVD2D::adjust_topology_of_interior_of_polygon_by_edge_flipping()
{    
    for (list<Generator2D*>::iterator i_generators = m_generators.begin(); i_generators != m_generators.end(); ++i_generators) {
        Generator2D* currGenerator = *i_generators;
        Generator2D::Generator_Type type = currGenerator->getType();

        if (type == Generator2D::Generator_Type::VERTEX_G) 
        {
            VertexGenerator2D* currVertexGen = (VertexGenerator2D*)currGenerator;

            if (currVertexGen->vertex_type() != VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE)
            {
                continue;
            }
            flatten_topology_of_VFace_boundary_of_reflex_vertex_gen_with_possible_flipping_out_of_VEdges(currVertexGen);
        }
        else if ( type == Generator2D::Generator_Type::EDGE_G ) {
            EdgeGenerator2D* currEdgeGen = (EdgeGenerator2D*)currGenerator;
            flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges(currEdgeGen);
        }
    }
}



void PolygonVD2D::flatten_topology_of_VFace_boundary_of_reflex_vertex_gen_with_possible_flipping_out_of_VEdges( VertexGenerator2D * vertexGenerator )
{
    EdgeGenerator2D* prevEdgeGen = vertexGenerator->get_previous_edge_generator();
    EdgeGenerator2D* nextEdgeGen = vertexGenerator->get_next_edge_generator();

    VEdge2D* edge_between_prevEdgeGen_N_vertexGen = rg_NULL;
    VEdge2D* edge_between_vertexGen_N_nextEdgeGen = rg_NULL;
    find_two_VEdges_between_vertex_generator_N_incident_edge_generators(vertexGenerator, prevEdgeGen, nextEdgeGen, edge_between_prevEdgeGen_N_vertexGen, edge_between_vertexGen_N_nextEdgeGen);

    VFace2D* childVFaceOfReflexVertexGen = getGeneratorWhichHasThisID(vertexGenerator->getID())->getOuterFace();

    // orthogonalize both VEdges incident to relocated VVertex by projections
    VVertex2D* vertexToBeProjected;
    if (edge_between_vertexGen_N_nextEdgeGen->getLeftFace() == childVFaceOfReflexVertexGen)
    {
        vertexToBeProjected = edge_between_vertexGen_N_nextEdgeGen->getEndVertex();
    }
    else
    {
        vertexToBeProjected = edge_between_vertexGen_N_nextEdgeGen->getStartVertex();
    }

    if (edge_between_prevEdgeGen_N_vertexGen->getLeftFace() == childVFaceOfReflexVertexGen)
    {
        vertexToBeProjected = edge_between_prevEdgeGen_N_vertexGen->getStartVertex();
    }
    else
    {
        vertexToBeProjected = edge_between_prevEdgeGen_N_vertexGen->getEndVertex();
    }

    list<VEdge2D*> boundaryEdgesOfVertexGen;
    childVFaceOfReflexVertexGen->getBoundaryVEdges(boundaryEdgesOfVertexGen);

    VEdge2D* startVEdge     = NULL;
    VEdge2D* currentVEdge   = NULL;
    VEdge2D* lastVEdgeInCCW = NULL;
    VEdge2D* nextVEdge      = NULL;
   

    if ( vertexGenerator->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
        startVEdge      = getCCWNextEdgeOnFace( edge_between_vertexGen_N_nextEdgeGen, childVFaceOfReflexVertexGen );
        currentVEdge    = startVEdge;
        lastVEdgeInCCW  = edge_between_prevEdgeGen_N_vertexGen;
        nextVEdge       = getCCWNextEdgeOnFace( currentVEdge, childVFaceOfReflexVertexGen );
    }
    else {
        startVEdge      = getCCWNextEdgeOnFace( edge_between_prevEdgeGen_N_vertexGen, childVFaceOfReflexVertexGen );
        currentVEdge    = startVEdge;
        lastVEdgeInCCW  = edge_between_vertexGen_N_nextEdgeGen;
        nextVEdge       = getCCWNextEdgeOnFace( currentVEdge, childVFaceOfReflexVertexGen );
    }

    rg_INT numberOfBoundaryVEdges = boundaryEdgesOfVertexGen.size();
    rg_Point2D vertexCoord = vertexGenerator->get_point();

    do
    {
        rg_Point2D consecutiveThreeVertexLocationsInCCW[2];

        if (currentVEdge->getLeftFace() == childVFaceOfReflexVertexGen) {
            consecutiveThreeVertexLocationsInCCW[0] = currentVEdge->getStartVertex()->getLocation();
            consecutiveThreeVertexLocationsInCCW[1] = currentVEdge->getEndVertex()->getLocation();
        }
        else {
            consecutiveThreeVertexLocationsInCCW[0] = currentVEdge->getEndVertex()->getLocation();
            consecutiveThreeVertexLocationsInCCW[1] = currentVEdge->getStartVertex()->getLocation();
        }

        rg_Point2D vec0 = consecutiveThreeVertexLocationsInCCW[0] - vertexCoord;
        rg_Point2D vec1 = consecutiveThreeVertexLocationsInCCW[1] - vertexCoord;

        if (rg_NEG(vec0.operator*(vec1)))
        {
            currentVEdge->flip();
            --numberOfBoundaryVEdges;
            VVertex2D* startVVertex = currentVEdge->getStartVertex();
            VVertex2D* endVVertex   = currentVEdge->getEndVertex();

            // refine VVertex coordinate
            rg_Circle2D refinedTangentCircle[2];
            refinedTangentCircle[0] = refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle(startVVertex);
            refinedTangentCircle[1] = refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle(endVVertex);
            startVVertex->setCircumcircle(refinedTangentCircle[0]);
            endVVertex->setCircumcircle(refinedTangentCircle[1]);
        }

        currentVEdge = nextVEdge;
        nextVEdge = getCCWNextEdgeOnFace(currentVEdge, childVFaceOfReflexVertexGen);

    } while (currentVEdge != lastVEdgeInCCW && numberOfBoundaryVEdges > 3);
}



void PolygonVD2D::flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges( EdgeGenerator2D * edgeGenerator )
{
    list<VEdge2D*> VEdgeChainCCW;
    find_interior_VEdge_chain_of_this_PEdge_CCW(edgeGenerator, VEdgeChainCCW);

    if (VEdgeChainCCW.size() == 2)
        return;

    VFace2D* VFaceOfEdgeGen = getGeneratorWhichHasThisID(edgeGenerator->getID())->getOuterFace();

    rg_Line2D edgeGenLineSeg;
    edgeGenerator->get_geometry(edgeGenLineSeg);
    edgeGenLineSeg = edgeGenLineSeg.get_reversed_line2D();
    rg_REAL currParameter = 0.0;
    rg_REAL prevParameter = 0.0;

    VEdge2D* prevVEdge = rg_NULL;
    list<VEdge2D*>::iterator i_VEdge = VEdgeChainCCW.begin();
    prevVEdge = *i_VEdge;
    ++i_VEdge;

    while (i_VEdge != VEdgeChainCCW.end())
    {		
        VEdge2D* currentVEdge = *i_VEdge;
        rg_Point2D pointToBeProjected;
        if (currentVEdge->getLeftFace() == VFaceOfEdgeGen)
            pointToBeProjected = currentVEdge->getEndVertex()->getLocation();
        else
            pointToBeProjected = currentVEdge->getStartVertex()->getLocation();

        edgeGenLineSeg.project(pointToBeProjected, currParameter);
        if (rg_GT(currParameter, 1.0) || rg_LT(currParameter, 0.0) || rg_LT(currParameter, prevParameter))
        {
            currentVEdge->flip();

            refine_VVertices_coordinates_of_flipped_VEdge(currentVEdge);

            if (prevVEdge->getLeftFace() == VFaceOfEdgeGen) {
                edgeGenLineSeg.project(prevVEdge->getEndVertex()->getLocation(), prevParameter);
            }
            else { 
                edgeGenLineSeg.project(prevVEdge->getStartVertex()->getLocation(), prevParameter);
            }
        }
        else
        {
            prevParameter = currParameter;
            prevVEdge = currentVEdge;
        }

        ++i_VEdge;
    }
}



void PolygonVD2D::flatten_topology_of_VFace_boundary_of_edge_gen_with_possible_flipping_out_of_VEdges_exterior( EdgeGenerator2D* edgeGenerator )
{
    list<VEdge2D*> VEdgeChainCCW;
    find_exterior_VEdge_chain_of_this_PEdge_CCW( edgeGenerator, VEdgeChainCCW );

    if ( VEdgeChainCCW.size() == 2 )
        return;

    VFace2D* VFaceOfEdgeGen = getGeneratorWhichHasThisID( edgeGenerator->getID() )->getOuterFace();

    rg_Line2D edgeGenLineSeg;
    edgeGenerator->get_geometry( edgeGenLineSeg );
    rg_REAL currParameter = 0.0;
    rg_REAL prevParameter = 0.0;

    VEdge2D* prevVEdge = rg_NULL;
    list<VEdge2D*>::iterator i_VEdge = VEdgeChainCCW.begin();
    prevVEdge = *i_VEdge;
    ++i_VEdge;

    while ( i_VEdge != VEdgeChainCCW.end() )
    {
        VEdge2D* currentVEdge = *i_VEdge;
        rg_Point2D pointToBeProjected;
        if ( currentVEdge->getLeftFace() == VFaceOfEdgeGen )
            pointToBeProjected = currentVEdge->getEndVertex()->getLocation();
        else
            pointToBeProjected = currentVEdge->getStartVertex()->getLocation();

        edgeGenLineSeg.project( pointToBeProjected, currParameter );
        if ( rg_GT( currParameter, 1.0 ) || rg_LT( currParameter, 0.0 ) || rg_LT( currParameter, prevParameter ) )
        {
            currentVEdge->flip();

            refine_VVertices_coordinates_of_flipped_VEdge( currentVEdge );

            if ( prevVEdge->getLeftFace() == VFaceOfEdgeGen ) {
                edgeGenLineSeg.project( prevVEdge->getEndVertex()->getLocation(), prevParameter );
            }
            else {
                edgeGenLineSeg.project( prevVEdge->getStartVertex()->getLocation(), prevParameter );
            }
        }
        else
        {
            prevParameter = currParameter;
            prevVEdge = currentVEdge;
        }

        ++i_VEdge;
    }
}



void PolygonVD2D::find_interior_VEdge_chain_of_this_PEdge_CCW( EdgeGenerator2D * edgeGenerator, list<VEdge2D*>& VEdgeChainCCW )
{
    Generator2D* prevGen = edgeGenerator->get_start_vertex_generator();
    Generator2D* nextGen = edgeGenerator->get_end_vertex_generator();

    if ( edgeGenerator->get_start_vertex_generator()->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
        prevGen = edgeGenerator->get_start_vertex_generator();
    }
    else {
        prevGen = edgeGenerator->get_start_vertex_generator()->get_previous_edge_generator();
    }

    if( edgeGenerator->get_end_vertex_generator()->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
        nextGen = edgeGenerator->get_end_vertex_generator();
    }
    else {
        nextGen = edgeGenerator->get_end_vertex_generator()->get_next_edge_generator();
    }

    VFace2D* prevVFace      = getGeneratorWhichHasThisID(prevGen->getID())->getOuterFace();
    VFace2D* nextVFace      = getGeneratorWhichHasThisID(nextGen->getID())->getOuterFace();
    VFace2D* edgeGenVFace   = getGeneratorWhichHasThisID(edgeGenerator->getID())->getOuterFace();
    
    list<VEdge2D*> boundaryEdgesOfEdgeGen;
    edgeGenVFace->getBoundaryVEdges(boundaryEdgesOfEdgeGen);

    VEdge2D* startVEdge = rg_NULL;
    list<VEdge2D*>::iterator i_edge;
    for (i_edge = boundaryEdgesOfEdgeGen.begin(); i_edge != boundaryEdgesOfEdgeGen.end(); ++i_edge)
    {
        VEdge2D* currBoundaryEdge = *i_edge;
        if (currBoundaryEdge->getLeftFace() == nextVFace || currBoundaryEdge->getRightFace() == nextVFace)
        {
            startVEdge = currBoundaryEdge;
            break;
        }
    }

    VEdge2D* currVEdge = startVEdge;
    while (currVEdge->getLeftFace() != prevVFace && currVEdge->getRightFace() != prevVFace)
    {
        VEdgeChainCCW.push_back(currVEdge);
        currVEdge = getCCWNextEdgeOnFace(currVEdge, edgeGenVFace);
    }
    VEdgeChainCCW.push_back(currVEdge);
}



void PolygonVD2D::find_exterior_VEdge_chain_of_this_PEdge_CCW( EdgeGenerator2D* edgeGenerator, list<VEdge2D*>& VEdgeChainCCW )
{
    Generator2D* prevGen = edgeGenerator->get_start_vertex_generator();
    Generator2D* nextGen = edgeGenerator->get_end_vertex_generator();

    if ( edgeGenerator->get_start_vertex_generator()->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
        prevGen = edgeGenerator->get_start_vertex_generator()->get_previous_edge_generator();
    }
    else {
        prevGen = edgeGenerator->get_start_vertex_generator();
    }

    if ( edgeGenerator->get_end_vertex_generator()->vertex_type() == VertexGenerator2D::Vertex_Type::REFLEX_FROM_POLYGON_INSIDE ) {
        nextGen = edgeGenerator->get_end_vertex_generator()->get_next_edge_generator();
    }
    else {
        nextGen = edgeGenerator->get_end_vertex_generator();
    }

    VFace2D* prevVFace = getGeneratorWhichHasThisID( prevGen->getID() )->getOuterFace();
    VFace2D* nextVFace = getGeneratorWhichHasThisID( nextGen->getID() )->getOuterFace();
    VFace2D* edgeGenVFace = getGeneratorWhichHasThisID( edgeGenerator->getID() )->getOuterFace();

    list<VEdge2D*> boundaryEdgesOfEdgeGen;
    edgeGenVFace->getBoundaryVEdges( boundaryEdgesOfEdgeGen );

    VEdge2D* startVEdge = rg_NULL;
    list<VEdge2D*>::iterator i_edge;
    for ( i_edge = boundaryEdgesOfEdgeGen.begin(); i_edge != boundaryEdgesOfEdgeGen.end(); ++i_edge )
    {
        VEdge2D* currBoundaryEdge = *i_edge;
        if ( currBoundaryEdge->getLeftFace() == prevVFace || currBoundaryEdge->getRightFace() == prevVFace )
        {
            startVEdge = currBoundaryEdge;
            break;
        }
    }

    VEdge2D* currVEdge = startVEdge;
    while ( currVEdge->getLeftFace() != nextVFace && currVEdge->getRightFace() != nextVFace )
    {
        VEdgeChainCCW.push_back( currVEdge );
        currVEdge = getCCWNextEdgeOnFace( currVEdge, edgeGenVFace );
    }
    VEdgeChainCCW.push_back( currVEdge );
}



void PolygonVD2D::refine_VVertices_coordinates_of_flipped_VEdge( VEdge2D * edge )
{
    VVertex2D* startVVertex = edge->getStartVertex();
    VVertex2D* endVVertex = edge->getEndVertex();
    rg_Circle2D refinedTangentCircle[2];
    bool tangentCircleLocatedPolygonOutward = false;
    refinedTangentCircle[0] = refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle(startVVertex);
    refinedTangentCircle[1] = refine_location_of_VVertex_by_finding_maximum_empty_tangent_circle(endVVertex);
    startVVertex->setCircumcircle(refinedTangentCircle[0]);
    endVVertex->setCircumcircle(refinedTangentCircle[1]);
}



rg_RQBzCurve2D PolygonVD2D::createRQBzCurveOfParabola( const rg_Circle2D& disk, const rg_Line2D& line )
{
    double radius_of_focus = disk.getRadius();
    rg_Point2D focus = disk.getCenterPt();
    rg_Line2D  directrix = line;
    rg_Point2D footprint_focus;
    directrix.compute_perpendicular_footprint_of_point_onto_entire_line(focus, footprint_focus);
    rg_Point2D vec_to_directrix = (footprint_focus - focus).getUnitVector();

    rg_Point2D dirVec = -vec_to_directrix;
    directrix.setSP(directrix.getSP() + vec_to_directrix * radius_of_focus);
    directrix.setEP(directrix.getEP() + vec_to_directrix * radius_of_focus);

    rg_Point2D footprint_passingPt = (directrix.getSP() + directrix.getEP()) / 2.0;

    rg_Point2D SP, PP, EP;

    double x_s = directrix.getSP().getX();
    double y_s = directrix.getSP().getY();
    double x_v = dirVec.getX();
    double y_v = dirVec.getY();
    double x_d = focus.getX();
    double y_d = focus.getY();

    double x_s_d = x_s - x_d;
    double y_s_d = y_s - y_d;

    // double a = x_v * x_v + y_v * y_v - 1.0;  // a is always 0 because x_v and y_v are elements of 
    double b = x_s_d * x_v + y_s_d * y_v;
    double c = x_s_d * x_s_d + y_s_d * y_s_d;

    double d = -c / b / 2.0;
    SP = directrix.getSP() + d * dirVec;

    double x_p = footprint_passingPt.getX();
    double y_p = footprint_passingPt.getY();
    double x_p_d = x_p - x_d;
    double y_p_d = y_p - y_d;

    double b_p = x_p_d * x_v + y_p_d * y_v;
    double c_p = x_p_d * x_p_d + y_p_d * y_p_d;

    double d_p = -c_p / b_p / 2.0;
    PP = footprint_passingPt + d_p * dirVec;

    double x_e = directrix.getEP().getX();
    double y_e = directrix.getEP().getY();
    double x_e_d = x_e - x_d;
    double y_e_d = y_e - y_d;

    double b_e = x_e_d * x_v + y_e_d * y_v;
    double c_e = x_e_d * x_e_d + y_e_d * y_e_d;

    double d_e = -c_e / b_e / 2.0;
    EP = directrix.getEP() + d_e * dirVec;

    rg_Point2D tangentVec_SP, tangentVec_EP;

    rg_Point2D vec1_s = (SP - focus).getUnitVector();
    //rg_Point2D vec1_s = (focus - SP).getUnitVector();
    rg_Point2D vec2_s = dirVec;
    if( rg_ZERO( (vec1_s+vec2_s).magnitude() ) )
        tangentVec_SP =  rg_Point2D( -vec1_s.getY(), vec1_s.getX() );
    else
        tangentVec_SP = (vec1_s + vec2_s).getUnitVector();


    rg_Point2D vec1_e = (EP - focus).getUnitVector();
    //rg_Point2D vec1_e = (focus - EP).getUnitVector();
    rg_Point2D vec2_e = dirVec;
    if( rg_ZERO( (vec1_e+vec2_e).magnitude() ) )
        tangentVec_EP =  rg_Point2D( -vec1_e.getY(), vec1_e.getX() );
    else
        tangentVec_EP = (vec1_e + vec2_e).getUnitVector();

    rg_RQBzCurve2D curve;
    curve.makeRQBezier(SP, tangentVec_SP, EP, tangentVec_EP, PP);

    return curve;
}



DiskGenerator2D* PolygonVD2D::createGeneratorOfInsertedDisk( const rg_Circle2D & disk )
{
    int lastID = m_generators.back()->getID();
    DiskGenerator2D* newGenerator = new DiskGenerator2D( lastID + 1, disk );
    newGenerator->setUserData( newGenerator );
    m_generators.push_back( newGenerator );

    return newGenerator;
}



VFace2D * PolygonVD2D::findAnchorCell( const rg_Point2D & pt )
{
    VFace2D*        anchorCell          = m_VFaces.back();
    Generator2D*    anchorGen           = (Generator2D*)anchorCell->getGenerator();
    double          shortestDistance    = 0.0;
    
    Generator2D::Generator_Type currType = anchorGen->getType();

    switch ( currType )
    {
    case Generator2D::Generator_Type::VERTEX_G:
    {
        shortestDistance = ((VertexGenerator2D*)anchorGen)->get_point().distance(pt);
    }
    break;

    case Generator2D::Generator_Type::EDGE_G:
    {
        rg_Line2D lineSeg;
        ((EdgeGenerator2D*)anchorGen)->get_geometry(lineSeg);
        shortestDistance = lineSeg.getDistance(pt);
    }
    break;

    case Generator2D::Generator_Type::DISK_G:
    {
        shortestDistance = ((DiskGenerator2D*)anchorGen)->getDisk().distance(pt);
    }
    break;

    default:
        break;
    }

    



    VEdge2D*            beginningEdge       = anchorCell->getFirstVEdge();
    VEdge2D*            currBoundaryEdge    = beginningEdge;
    VFace2D*            currOppositeFace    = NULL;
    Generator2D*        neighborGenerator   = NULL;
    double              distanceToNeighbor  = DBL_MAX;


    do {
        currOppositeFace    = getFaceOfOppositeSide( anchorCell, currBoundaryEdge );
        
        Generator2D::Generator_Type type = Generator2D::Generator_Type::UNKNOWN_G;
        if( currOppositeFace->getGenerator() != NULL ) {
            neighborGenerator   = (Generator2D*)currOppositeFace->getGenerator();
            type = neighborGenerator->getType();

            switch ( type )
            {
            case Generator2D::Generator_Type::VERTEX_G:
            {
                distanceToNeighbor = ((VertexGenerator2D*)neighborGenerator)->get_point().distance(pt);
            }
            break;

            case Generator2D::Generator_Type::EDGE_G:
            {
                rg_Line2D lineSeg;
                ((EdgeGenerator2D*)neighborGenerator)->get_geometry(lineSeg);
                distanceToNeighbor = lineSeg.getDistance(pt);
            }
            break;

            case Generator2D::Generator_Type::DISK_G:
            {
                distanceToNeighbor = ((DiskGenerator2D*)neighborGenerator)->getDisk().distance(pt);
            }
            break;

            default:
                break;
            }
        }

        if ( currType == Generator2D::Generator_Type::EDGE_G && type == Generator2D::Generator_Type::VERTEX_G ) {
            if ( ( (EdgeGenerator2D*)anchorCell->getGenerator() )->get_start_vertex_generator() == neighborGenerator
                || ( (EdgeGenerator2D*)anchorCell->getGenerator() )->get_end_vertex_generator() == neighborGenerator ) 
            {
                if ( rg_EQ( distanceToNeighbor, shortestDistance ) ) {
                    anchorCell = currOppositeFace;
                    shortestDistance = distanceToNeighbor;
                    beginningEdge = currBoundaryEdge;
                    currType = type;
                }
            }
        }
        else {
            if ( distanceToNeighbor < shortestDistance ) {
                anchorCell = currOppositeFace;
                shortestDistance = distanceToNeighbor;
                beginningEdge = currBoundaryEdge;
                currType = type;
            }
        }

        currBoundaryEdge = getCCWNextEdgeOnFace( currBoundaryEdge, anchorCell );

    } while ( currBoundaryEdge != beginningEdge );
    

    return anchorCell;
}



void PolygonVD2D::makeNewCellAndConnectToCurrVD( DiskGenerator2D * newGenerator, list<VVertex2D*>& redVVertices, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& newVVertices, list<VEdge2D*>& newVEdges )
{
    list<VEdge2D*> crossingEdges;
    findCrossingEdges( redVVertices, crossingEdges );

    reconfigureByConvertingBlueVVerticesToWhite( blueVVertices );

    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( newGenerator, crossingEdges, newVVertices, newVEdges );
}



VVertex2D * PolygonVD2D::findSeedRedVertex( DiskGenerator2D * const newGenerator, Generator2D * const closestGenerator )
{
    list<VVertex2D*> bndVerticesOfAnchorCell;
    closestGenerator->getOuterFace()->getBoundaryVVertices( bndVerticesOfAnchorCell );

    double      minMUValue = DBL_MAX;
    VVertex2D*  firstRedVertex = NULL;

    for( list<VVertex2D*>::iterator i_vtx = bndVerticesOfAnchorCell.begin() ; i_vtx != bndVerticesOfAnchorCell.end() ; ++i_vtx ) {
        VVertex2D* currVVertex = *i_vtx;

        if ( get_location_status_of_Vvertex( currVVertex ) == ON_POLYGON_BOUNDARY ) {
            continue;
        }

        double currMUValue = computeMUValue ( currVVertex, newGenerator );
        if( currMUValue < minMUValue ) {
            minMUValue      = currMUValue;
            firstRedVertex  = currVVertex;
        }
    }

    if( minMUValue >= 0.0) {
        list<VEdge2D*> boundaryVEdges;
        closestGenerator->getOuterFace()->getBoundaryVEdges( boundaryVEdges );

        for( list<VEdge2D*>::iterator i_edge = boundaryVEdges.begin() ; i_edge != boundaryVEdges.end() ; ++i_edge ) {
            VEdge2D* currVEdge = *i_edge;

            if( isAnomalizingEdge( currVEdge, newGenerator ) ) {
                VVertex2D* fictitiousVVertex = splitVEdgeAtFictitiousVVertex( currVEdge );
                firstRedVertex = fictitiousVVertex;
                break;
            }
        }
    }

    return firstRedVertex;
}



VVertex2D* PolygonVD2D::splitVEdgeAtFictitiousVVertex( VEdge2D* const currVEdge )
{
    VVertex2D* fictitiousVVertex = createVertex(m_VVertices.back()->getID()+1);
    fictitiousVVertex->setStatus( YELLOW_V );
    fictitiousVVertex->setFirstVEdge( currVEdge );
    fictitiousVVertex->setFictitious();

    VEdge2D* dividedVEdge = createEdge(m_VEdges.back()->getID() + 1);
    m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>(dividedVEdge, INSIDE_POLYGON));    

    if( currVEdge->getLeftHand()->getLeftLeg() == currVEdge ) {
        currVEdge->getLeftHand()->setLeftLeg( dividedVEdge );
    }
    else {
        currVEdge->getLeftHand()->setRightHand( dividedVEdge );
    }

    if( currVEdge->getRightHand()->getRightLeg() == currVEdge ) {
        currVEdge->getRightHand()->setRightLeg( dividedVEdge );
    }
    else {
        currVEdge->getRightHand()->setLeftHand( dividedVEdge );
    }

    dividedVEdge->setLeftHand(      currVEdge->getLeftHand() );
    dividedVEdge->setRightHand(     currVEdge->getRightHand() );
    dividedVEdge->setLeftLeg(       currVEdge );
    dividedVEdge->setRightLeg(      currVEdge );
    dividedVEdge->setStartVertex(   fictitiousVVertex );
    dividedVEdge->setEndVertex(     currVEdge->getEndVertex() );
    dividedVEdge->setLeftFace(      currVEdge->getLeftFace() );
    dividedVEdge->setRightFace(     currVEdge->getRightFace() );

    currVEdge->setEndVertex(        fictitiousVVertex );
    currVEdge->setLeftHand(         dividedVEdge );
    currVEdge->setRightHand(        dividedVEdge );

    dividedVEdge->getEndVertex()->setFirstVEdge( dividedVEdge );

    return fictitiousVVertex;
}



bool PolygonVD2D::isAnomalizingEdge_for_DiskGenerator( VEdge2D * const incidentEdge, Generator2D * const newGenerator )
{
    if( incidentEdge->getStartVertex()->isFictitious() || incidentEdge->getEndVertex()->isFictitious() ) {
        return false;
    }

    Generator2D* leftGen  = (Generator2D*)incidentEdge->getLeftFace()->getGenerator();
    Generator2D* rightGen = (Generator2D*)incidentEdge->getRightFace()->getGenerator();

    Generator2D::Generator_Type leftType  = leftGen->getType();
    Generator2D::Generator_Type rightType = rightGen->getType();

    bool this_edge_is_an_anomalyEdge = false;

    switch ( leftType ) {
    case Generator2D::Generator_Type::DISK_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            this_edge_is_an_anomalyEdge = VoronoiDiagram2DC::isAnomalizingEdge(incidentEdge, newGenerator);
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;
            double radius_of_disk = diskGen->getDisk().getRadius();

            rg_Point2D focus = diskGen->getDisk().getCenterPt();
            rg_Line2D  directrix = edgeGen->get_geometry();
            rg_Point2D footprint_focus;
            directrix.compute_perpendicular_footprint_of_point_onto_entire_line(focus, footprint_focus);
            rg_Point2D vec_to_directrix = footprint_focus - focus;
            vec_to_directrix = vec_to_directrix.getUnitVector();
            directrix.setSP(directrix.getSP() + vec_to_directrix * radius_of_disk);
            directrix.setEP(directrix.getEP() + vec_to_directrix * radius_of_disk);
            Parabola2D parabola(focus, directrix);

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola(incidentEdge, parabola, newGenerator, radius_of_disk);
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola(incidentEdge, newGenerator->getDisk(), diskGen->getDisk(), rg_Circle2D(vertexGen->get_point(), 0.0));
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::EDGE_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;
            double radius_of_disk = diskGen->getDisk().getRadius();

            rg_Point2D focus = diskGen->getDisk().getCenterPt();
            rg_Line2D  directrix = edgeGen->get_geometry();
            rg_Point2D footprint_focus;
            directrix.compute_perpendicular_footprint_of_point_onto_entire_line(focus, footprint_focus);
            rg_Point2D vec_to_directrix = footprint_focus - focus;
            vec_to_directrix = vec_to_directrix.getUnitVector();
            directrix.setSP(directrix.getSP() + vec_to_directrix * radius_of_disk);
            directrix.setEP(directrix.getEP() + vec_to_directrix * radius_of_disk);
            Parabola2D parabola(focus, directrix);

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola(incidentEdge, parabola, newGenerator, radius_of_disk);
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            rg_Circle2D circumcircle_sp = incidentEdge->getStartVertex()->getCircumcircle();
            rg_Circle2D circumcircle_ep = incidentEdge->getEndVertex()->getCircumcircle();

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_line(rg_Line2D(circumcircle_sp.getCenterPt(), circumcircle_ep.getCenterPt()), circumcircle_sp, circumcircle_ep, newGenerator);
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)leftGen;
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)rightGen;

            if ( edgeGen->get_start_vertex_generator() == vertexGen || edgeGen->get_end_vertex_generator() == vertexGen )
                return false;
            
            rg_Point2D focus = vertexGen->get_point();
            rg_Line2D  directrix = edgeGen->get_geometry();
            Parabola2D parabola(focus, directrix);

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola(incidentEdge, parabola, newGenerator, 0.0);
        }
        break;

        default:
            break;
        }
    }
    break;

    case Generator2D::Generator_Type::VERTEX_G:
    {
        switch ( rightType ) {
        case Generator2D::Generator_Type::DISK_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            DiskGenerator2D* diskGen = (DiskGenerator2D*)rightGen;

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_hyperbola(incidentEdge, newGenerator->getDisk(), diskGen->getDisk(), rg_Circle2D(vertexGen->get_point(), 0.0));
        }
        break;

        case Generator2D::Generator_Type::EDGE_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)leftGen;
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)rightGen;

            if ( edgeGen->get_start_vertex_generator() == vertexGen || edgeGen->get_end_vertex_generator() == vertexGen )
                return false;

            rg_Point2D focus = vertexGen->get_point();
            rg_Line2D  directrix = edgeGen->get_geometry();
            Parabola2D parabola(focus, directrix);

            this_edge_is_an_anomalyEdge = isAnomalizingEdge_parabola(incidentEdge, parabola, newGenerator, 0.0);
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            // this_edge_is_an_anomalyEdge = false;
        }
        break;

        default:
            break;
        }
    }
    break;

    default:
        break;
    }

    return this_edge_is_an_anomalyEdge;
}



bool PolygonVD2D::isAnomalizingEdge_parabola( VEdge2D * const incidentEdge, const Parabola2D & parabola, Generator2D * const newGenerator, const double & radiusOfFocus )
{
    rg_Point2D focus = parabola.get_focus();
    rg_Line2D  directrix = parabola.get_directrix();

    rg_Point2D footprint_focus;
    directrix.compute_perpendicular_footprint_of_point_onto_entire_line(focus, footprint_focus);

    rg_Point2D translatePt = (focus + footprint_focus) / 2.0;
    double angle = parabola.get_rotation_angle_from_positive_X_axis_for_axis_of_symmetry();

    rg_TMatrix2D transformMat;
    transformMat.translate(-translatePt);
    transformMat.rotate(-angle);

    rg_Point2D transformedSP = transformMat * incidentEdge->getStartVertex()->getLocation();
    rg_Point2D transformedEP = transformMat * incidentEdge->getEndVertex()->getLocation();
    rg_Point2D transformedPt = transformMat * newGenerator->getDisk().getCenterPt();
    double x_1 = transformedPt.getX();
    double y_1 = transformedPt.getY();
    double r_1 = newGenerator->getDisk().getRadius();

    double p   = directrix.getDistance(focus) / 2.0;
    double t   = p + r_1 - radiusOfFocus;
    double c   = x_1 * x_1 + y_1 * y_1 - t * t;
    double b   = x_1;
    double a   = 1.0 - (y_1 + t) / 2.0 / p;

    double value_in_squared = b * b -  a * c;

    if ( value_in_squared <= 0.0 )
        return false;

    if ( rg_ZERO(a, resNeg6) )
        return false;

    double squared_ = sqrt( value_in_squared ); 
    double x_intersectionPt_1 = ( b - squared_ ) / a;
    double x_intersectionPt_2 = ( b + squared_ ) / a;
    double x_sp = transformedSP.getX();
    double x_ep = transformedEP.getX();

    if ( x_sp > x_ep ) {
        double x_temp = x_sp;
        x_sp = x_ep;
        x_ep = x_temp;
    }

    if (   x_intersectionPt_1 > x_sp && x_intersectionPt_1 < x_ep
        && x_intersectionPt_2 > x_sp && x_intersectionPt_2 < x_ep ) {
        return true;
    }
    else {
        return false;
    }
}

bool PolygonVD2D::isAnomalizingEdge_hyperbola( VEdge2D * const incidentEdge, const rg_Circle2D & disk1, const rg_Circle2D & disk2, const rg_Circle2D & disk3 )
{
    rg_Circle2D firstCircle;
    rg_Circle2D secondCircle;

    int numCircumcircles = rg_Circle2D::makeCircumcircle( disk1, disk2, disk3, firstCircle, secondCircle );

    rg_Point2D center0;
    rg_Point2D CWBoundary;
    rg_Point2D CCWBoundary;

    rg_Point2D firstLocation;
    rg_Point2D secondLocation;

    rg_Point3D CWBoundary3D;
    rg_Point3D CCWBoundary3D;
    rg_Point3D firstLocation3D;
    rg_Point3D secondLocation3D;

    rg_Point3D CW_First_Test;
    rg_Point3D CCW_First_Test;
    rg_Point3D CW_Second_Test;
    rg_Point3D CCW_Second_Test;

    bool this_edge_is_anomaly_edge = false;

    switch (numCircumcircles) {

    case 0:
    case 1:
    {
        // this_edge_is_anomaly_edge = false;
    }
    break;
        
    case 2:
    {
        if( incidentEdge->getLeftFace()->getGenerator()->getDisk().getRadius() < incidentEdge->getRightFace()->getGenerator()->getDisk().getRadius() ) {
            center0     = incidentEdge->getRightFace()->getGenerator()->getDisk().getCenterPt();
            CWBoundary  = incidentEdge->getEndVertex()->getLocation()   - center0;
            CCWBoundary = incidentEdge->getStartVertex()->getLocation() - center0;
        }
        else {
            center0     = incidentEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt() ;
            CWBoundary  = incidentEdge->getStartVertex()->getLocation() - center0;
            CCWBoundary = incidentEdge->getEndVertex()->getLocation()   - center0;
        }

        firstLocation   = firstCircle.getCenterPt();
        secondLocation  = secondCircle.getCenterPt();
        firstLocation   = firstLocation  - center0;
        secondLocation  = secondLocation - center0;

        CWBoundary3D     = rg_Point3D( CWBoundary );
        CCWBoundary3D    = rg_Point3D( CCWBoundary );
        firstLocation3D  = rg_Point3D( firstLocation );
        secondLocation3D = rg_Point3D( secondLocation );

        CW_First_Test   = CWBoundary3D.crossProduct( firstLocation3D );
        CCW_First_Test  = CCWBoundary3D.crossProduct( firstLocation3D );
        CW_Second_Test  = CWBoundary3D.crossProduct( secondLocation3D );
        CCW_Second_Test = CCWBoundary3D.crossProduct( secondLocation3D );

        if( CW_First_Test.getZ() < 0 ) {
            break;
        }
        if( CW_Second_Test.getZ() < 0 ) {
            break;
        }
        if( CCW_First_Test.getZ() > 0 ) {
            break;
        }
        if( CCW_Second_Test.getZ() > 0 ) {
            break;
        }

        this_edge_is_anomaly_edge = true;
    }
    break;

    default:
        break;
    }

    return this_edge_is_anomaly_edge;
}



bool PolygonVD2D::isAnomalizingEdge_line( const rg_Line2D & lineSeg, const rg_Circle2D& circumcircle_sp, const rg_Circle2D& circumcircle_ep, Generator2D * const newGenerator )
{
    double coeff_x, coeff_y, coeff_const;
    lineSeg.get_coefficients_of_implicit_form_of_line_equation(coeff_x, coeff_y, coeff_const);
    double angle_from_positive_X_axis = 0.0;
    if (rg_ZERO(coeff_x))
        angle_from_positive_X_axis = rg_PI / 2.0;
    else if (rg_ZERO(coeff_y))
        angle_from_positive_X_axis = 0.0;
    else {
        rg_Point2D lineDirVec(-coeff_y, coeff_x);
        rg_Point2D xAxisDirVec(1.0, 0.0);
        //angle_from_positive_X_axis = rg_GeoFunc::calculateAngle(lineDirVec, xAxisDirVec);
        angle_from_positive_X_axis = angleFromVec1toVec2(xAxisDirVec, lineDirVec);
    }

    rg_TMatrix2D transformMat;
    transformMat.translate(-lineSeg.getSP());
    transformMat.rotate(-angle_from_positive_X_axis);

    rg_Point2D transformed_newPt = transformMat * newGenerator->getDisk().getCenterPt();
    double x_2 = lineSeg.getLength();
    double x_3 = transformed_newPt.getX();
    double y_3 = transformed_newPt.getY();
    double r_3 = newGenerator->getDisk().getRadius();
    double d_1 = circumcircle_sp.getRadius();
    double d_2 = circumcircle_ep.getRadius();

    double t_1 = d_1 + r_3;
    double t_2 = d_2 - d_1;

    double a   = x_2 * x_2 - t_2 * t_2;
    double b   = -x_2 * x_3 - t_1 * t_2;
    double c   = x_3 * x_3 + y_3 * y_3 - t_1 * t_1;

    double inside_square = b * b - a * c;

    if ( inside_square <= 0 )
        return false;

    double squared = sqrt(inside_square);
    double i_1 = (- b - squared) / a;
    double i_2 = ( -b + squared) / a;

    if (   i_1 > 0 && i_1 < 1
        && i_2 > 0 && i_2 < 1 ) {
        return true;
    }
    else
        return false;
}



void PolygonVD2D::markLocationStatusAtNewEntities( list<VVertex2D*>& newVertices, list<VEdge2D*>& newEdges )
{
    for ( list<VVertex2D*>::iterator i_vtx = newVertices.begin(); i_vtx != newVertices.end(); ++i_vtx ) {
        m_VVertexToLocationStatus.insert(pair<VVertex2D*, VDEntity_Location_Status>((*i_vtx), INSIDE_POLYGON));
    }

    for ( list<VEdge2D*>::iterator i_edge = newEdges.begin(); i_edge != newEdges.end(); ++i_edge ) {
        m_VEdgeToLocationStatus.insert(pair<VEdge2D*, VDEntity_Location_Status>((*i_edge), INSIDE_POLYGON));
    }
}



void PolygonVD2D::markLocationStatusAtNewEntities_outside( list<VVertex2D*>& newVertices, list<VEdge2D*>& newEdges )
{
    for ( list<VVertex2D*>::iterator i_vtx = newVertices.begin(); i_vtx != newVertices.end(); ++i_vtx ) {
        m_VVertexToLocationStatus.insert( pair<VVertex2D*, VDEntity_Location_Status>( ( *i_vtx ), OUTSIDE_POLYGON ) );
    }

    for ( list<VEdge2D*>::iterator i_edge = newEdges.begin(); i_edge != newEdges.end(); ++i_edge ) {
        m_VEdgeToLocationStatus.insert( pair<VEdge2D*, VDEntity_Location_Status>( ( *i_edge ), OUTSIDE_POLYGON ) );
    }
}



void PolygonVD2D::computeCoordOfNewVVertices( list<VVertex2D*>& newVVertices )
{
    computeCoordOfNewVertices( newVVertices );
    return;
#ifdef DEBUG_VERTEX
    ofstream fout_debug_vertex("vertex_debugging.txt");
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
        count_number_of_occurred_generator_types_and_get_parent_generators_at_VVertex(vertex, numDisks, numLines, numPoints, parentGens);

#ifdef CHECK_COMP_TIME
        endTime = clock();
        t_check_parents = t_check_parents + endTime - startTime;
#endif

        //rg_Circle2D refinedTangentCircle;
        rg_Circle2D circumcircle_of_vertex;
        rg_Circle2D tangentCircle[2];

        switch (numDisks)
        {
        case 3: // numDisks == 3, numLines == 0, numPoints == 0
        {
#ifdef CHECK_COMP_TIME
            ++num_DDD;
            startTime = clock();
#endif
            rg_Circle2D disks[3];
            disks[0] = ((DiskGenerator2D*)parentGens[0])->getDisk();
            disks[1] = ((DiskGenerator2D*)parentGens[1])->getDisk();
            disks[2] = ((DiskGenerator2D*)parentGens[2])->getDisk();

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
            if (numLines == 1)
            {
#ifdef CHECK_COMP_TIME
                ++num_DDL;
                startTime = clock();
#endif
                rg_Circle2D disks[2];
                rg_Line2D lineSeg3;
                rg_INT diskIndex = 0;
                for (rg_INT i = 0; i < 3; i++)
                {
                    if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                        disks[diskIndex++] = ((DiskGenerator2D*)parentGens[i])->getDisk();
                    if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                        ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg3);
                }

                int numTangentCircles = 0;
                rg_Circle2D tangentCircle[2];
                numTangentCircles = computeTangentCircles_of_two_disks_and_a_line(disks[0], disks[1], lineSeg3, tangentCircle[0], tangentCircle[1]);

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

                /*
                double      radius_of_focus = disks[0].getRadius();
                rg_Point2D  focus           = disks[0].getCenterPt();
                rg_Line2D   directrix       = lineSeg3;
                rg_Point2D  normalVec       = directrix.getNormalVector().getUnitVector();
                
                directrix.setSP(directrix.getSP() - normalVec * radius_of_focus);
                directrix.setEP(directrix.getEP() - normalVec * radius_of_focus);
                Parabola2D parabola(focus, directrix);

                computeCoordOfNewVVertex_on_parabola(parabola, disks[1], radius_of_focus, tangentCircle[0], tangentCircle[1]);
                
                rg_Point2D footprint_debug[3];
                lineSeg3.compute_footprint_of_point_onto_line_segment(tangentCircle[0].getCenterPt(), footprint_debug[0]);
                disks[0].compute_perpendicular_footprint_of_point_onto_circle(tangentCircle[0].getCenterPt(), footprint_debug[1]);
                disks[1].compute_perpendicular_footprint_of_point_onto_circle(tangentCircle[0].getCenterPt(), footprint_debug[2]);

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
                */



                /////////////////// numerical from initial center /////////////////////
                /*
                rg_Point2D targetCenter;

                {
                    rg_Point2D init_tangentPt[3];
                    lineSeg3.compute_footprint_of_point_onto_line_segment(disks[0].getCenterPt(), init_tangentPt[0]);
                    lineSeg3.compute_footprint_of_point_onto_line_segment(disks[1].getCenterPt(), init_tangentPt[1]);
                    init_tangentPt[0] = (init_tangentPt[0] + init_tangentPt[1]) / 2.0;
                    disks[0].compute_perpendicular_footprint_of_point_onto_circle(init_tangentPt[0], init_tangentPt[1]);
                    disks[1].compute_perpendicular_footprint_of_point_onto_circle(init_tangentPt[0], init_tangentPt[2]);

                    targetCenter = (init_tangentPt[0] + init_tangentPt[1] + init_tangentPt[2]) / 3.0;
                }

                rg_Point2D currentCenter;
                rg_Circle2D refinedTangentCircle;
                rg_INT numIterations = 0;
                do
                {
                    currentCenter = targetCenter;

                    rg_Point2D footprint[3];
                    lineSeg3.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[0]);
                    disks[0].compute_perpendicular_footprint_of_point_onto_circle(currentCenter, footprint[1]);
                    disks[1].compute_perpendicular_footprint_of_point_onto_circle(currentCenter, footprint[2]);

                    refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
                    targetCenter = refinedTangentCircle.getCenterPt();
                    ++numIterations;
                } while (targetCenter.distance(currentCenter) >= 0.001 && numIterations <= 20);

                double radius = refinedTangentCircle.getCenterPt().distance(disks[0].getCenterPt()) - disks[0].getRadius();
                refinedTangentCircle.setRadius(radius);

                circumcircle_of_vertex = refinedTangentCircle;
                */



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
                for (rg_INT i = 0; i < 3; i++)
                {
                    if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                        disks[diskIndex++] = ((DiskGenerator2D*)parentGens[i])->getDisk();
                    if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                        ((VertexGenerator2D*)parentGens[i])->get_geometry(point3);
                }

                int numTangentCircles = rg_Circle2D::makeCircumcircle( disks[0], disks[1], rg_Circle2D(point3, 0.0), tangentCircle[0], tangentCircle[1] );
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
            switch (numLines)
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
                for (rg_INT i = 0; i < 3; i++)
                {
                    if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                        disk1 = ((DiskGenerator2D*)parentGens[i])->getDisk();
                    if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                        ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg[lineSegIndex++]);
                }


                VEdge2D* edge_of_two_polygon_edges = NULL;
                list<VEdge2D*> incidentEdges;
                vertex->getIncident3VEdges(incidentEdges);
                for ( list<VEdge2D*>::iterator i_edge = incidentEdges.begin(); i_edge != incidentEdges.end(); ++i_edge ) {
                    VEdge2D* currEdge = *i_edge;

                    Generator2D* leftParent  = (Generator2D*)currEdge->getLeftFace()->getGenerator()->getUserData();
                    Generator2D* rightParent = (Generator2D*)currEdge->getRightFace()->getGenerator()->getUserData();

                    if ( leftParent->getType() == Generator2D::Generator_Type::EDGE_G && rightParent ->getType() == Generator2D::Generator_Type::EDGE_G ) {
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
                rg_Line2D geometry_edge(SP, EP);

                rg_Point2D footprint_SP, footprint_EP;
                lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( SP, footprint_SP );
                lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line( EP, footprint_EP );

                double radius_SP = SP.distance( footprint_SP );
                double radius_EP = EP.distance( footprint_EP );

                
                


                computeCoordOfNewVVertex_on_line_of_two_polygon_edges( geometry_edge, rg_Circle2D(SP, radius_SP), rg_Circle2D(EP, radius_EP), disk1, tangentCircle[0], tangentCircle[1] );
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
                lineSeg[0].compute_perpendicular_footprint_of_point_onto_entire_line(circumcircle_of_vertex.getCenterPt(), foorprint_debug[0]);
                lineSeg[1].compute_perpendicular_footprint_of_point_onto_entire_line(circumcircle_of_vertex.getCenterPt(), foorprint_debug[1]);
                insertedDisk = disk1;
#endif
            }
            break;

            case 1: // numDisks == 1, numLines == 1, numPoints == 1
            {
#ifdef CHECK_COMP_TIME
                startTime = clock();
#endif
                EdgeGenerator2D*     edgeGen = rg_NULL;
                VertexGenerator2D* vertexGen = rg_NULL;

                rg_Circle2D disk1;
                rg_Line2D  lineSeg2;
                rg_Point2D point3;
                for (rg_INT i = 0; i < 3; i++)
                {
                    if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType()) {
                        disk1 = ((DiskGenerator2D*)parentGens[i])->getDisk();
                    }
                    if (Generator2D::Generator_Type::EDGE_G == parentGens[i]->getType())
                    {
                        edgeGen = (EdgeGenerator2D*)parentGens[i];
                        ((EdgeGenerator2D*)parentGens[i])->get_geometry(lineSeg2);
                    }
                    if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                    {
                        vertexGen = (VertexGenerator2D*)parentGens[i];
                        ((VertexGenerator2D*)parentGens[i])->get_geometry(point3);
                    }
                }
#ifdef CHECK_COMP_TIME
                endTime = clock();
                double localCompTime = endTime - startTime;
#endif
              

                if ( vertexGen->get_next_edge_generator() != edgeGen && vertexGen->get_previous_edge_generator() != edgeGen ) {
#ifdef CHECK_COMP_TIME
                    num_DDL++;
                    startTime = clock();
#endif
                    int numTangentCircles = 0;
                    rg_Circle2D tangentCircle[2];
                    numTangentCircles = computeTangentCircles_of_two_disks_and_a_line(rg_Circle2D(point3, 0.0), disk1, lineSeg2, tangentCircle[0], tangentCircle[1]);
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
                    lineSeg2.compute_perpendicular_footprint_of_point_onto_entire_line(circumcircle_of_vertex.getCenterPt(), foorprint_debug[0]);
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
                for (rg_INT i = 0; i < 3; i++)
                {
                    if (Generator2D::Generator_Type::DISK_G == parentGens[i]->getType())
                        disk1 = ((DiskGenerator2D*)parentGens[i])->getDisk();
                    if (Generator2D::Generator_Type::VERTEX_G == parentGens[i]->getType())
                        ((VertexGenerator2D*)parentGens[i])->get_geometry(point[pointIndex++]);
                }

                int numTangentCircles = rg_Circle2D::makeCircumcircle( disk1, rg_Circle2D(point[0], 0.0), rg_Circle2D(point[1], 0.0), tangentCircle[0], tangentCircle[1] );
                
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
        vertex->setCircumcircle(circumcircle_of_vertex);
        vertex->setStatus( WHITE_V );
    }
#ifdef DEBUG_VERTEX
    fout_debug_vertex.close();
#endif
}



void PolygonVD2D::computeCoordOfNewVVertices_outside( list<VVertex2D*>& newVVertices )
{
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
        count_number_of_occurred_generator_types_and_get_parent_generators_at_VVertex( vertex, numDisks, numLines, numPoints, parentGens );

#ifdef CHECK_COMP_TIME
        endTime = clock();
        t_check_parents = t_check_parents + endTime - startTime;
#endif

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

                lineSeg3 = lineSeg3.get_reversed_line2D();
                int numTangentCircles = 0;
                rg_Circle2D tangentCircle[2];
                numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( disks[0], disks[1], lineSeg3, tangentCircle[0], tangentCircle[1] );

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

                /*
                double      radius_of_focus = disks[0].getRadius();
                rg_Point2D  focus           = disks[0].getCenterPt();
                rg_Line2D   directrix       = lineSeg3;
                rg_Point2D  normalVec       = directrix.getNormalVector().getUnitVector();

                directrix.setSP(directrix.getSP() - normalVec * radius_of_focus);
                directrix.setEP(directrix.getEP() - normalVec * radius_of_focus);
                Parabola2D parabola(focus, directrix);

                computeCoordOfNewVVertex_on_parabola(parabola, disks[1], radius_of_focus, tangentCircle[0], tangentCircle[1]);

                rg_Point2D footprint_debug[3];
                lineSeg3.compute_footprint_of_point_onto_line_segment(tangentCircle[0].getCenterPt(), footprint_debug[0]);
                disks[0].compute_perpendicular_footprint_of_point_onto_circle(tangentCircle[0].getCenterPt(), footprint_debug[1]);
                disks[1].compute_perpendicular_footprint_of_point_onto_circle(tangentCircle[0].getCenterPt(), footprint_debug[2]);

                if ( this_circumcircle_has_right_orientation( vertex, tangentCircle[0] ) ) {
                    circumcircle_of_vertex = tangentCircle[0];
                }
                else {
                    circumcircle_of_vertex = tangentCircle[1];
                }
                */



                /////////////////// numerical from initial center /////////////////////
                /*
                rg_Point2D targetCenter;

                {
                    rg_Point2D init_tangentPt[3];
                    lineSeg3.compute_footprint_of_point_onto_line_segment(disks[0].getCenterPt(), init_tangentPt[0]);
                    lineSeg3.compute_footprint_of_point_onto_line_segment(disks[1].getCenterPt(), init_tangentPt[1]);
                    init_tangentPt[0] = (init_tangentPt[0] + init_tangentPt[1]) / 2.0;
                    disks[0].compute_perpendicular_footprint_of_point_onto_circle(init_tangentPt[0], init_tangentPt[1]);
                    disks[1].compute_perpendicular_footprint_of_point_onto_circle(init_tangentPt[0], init_tangentPt[2]);

                    targetCenter = (init_tangentPt[0] + init_tangentPt[1] + init_tangentPt[2]) / 3.0;
                }

                rg_Point2D currentCenter;
                rg_Circle2D refinedTangentCircle;
                rg_INT numIterations = 0;
                do
                {
                    currentCenter = targetCenter;

                    rg_Point2D footprint[3];
                    lineSeg3.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[0]);
                    disks[0].compute_perpendicular_footprint_of_point_onto_circle(currentCenter, footprint[1]);
                    disks[1].compute_perpendicular_footprint_of_point_onto_circle(currentCenter, footprint[2]);

                    refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
                    targetCenter = refinedTangentCircle.getCenterPt();
                    ++numIterations;
                } while (targetCenter.distance(currentCenter) >= 0.001 && numIterations <= 20);

                double radius = refinedTangentCircle.getCenterPt().distance(disks[0].getCenterPt()) - disks[0].getRadius();
                refinedTangentCircle.setRadius(radius);

                circumcircle_of_vertex = refinedTangentCircle;
                */



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
                    numTangentCircles = computeTangentCircles_of_two_disks_and_a_line( rg_Circle2D( point3, 0.0 ), disk1, lineSeg2, tangentCircle[0], tangentCircle[1] );
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



void PolygonVD2D::computeCoordOfNewVVertex_on_parabola( const Parabola2D & parabola, const rg_Circle2D & disk, const double & radiusOfFocus, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    rg_Point2D focus = parabola.get_focus();
    rg_Line2D  directrix = parabola.get_directrix();

    rg_Point2D footprint_focus;
    directrix.compute_perpendicular_footprint_of_point_onto_entire_line(focus, footprint_focus);

    rg_Point2D translatePt = (focus + footprint_focus) / 2.0;
    double angle = parabola.get_rotation_angle_from_positive_X_axis_for_axis_of_symmetry();

    rg_TMatrix2D transformMat;
    transformMat.translate(-translatePt);
    transformMat.rotate(-angle);

    rg_Point2D transformedPt = transformMat * disk.getCenterPt();
    double x_1 = transformedPt.getX();
    double y_1 = transformedPt.getY();
    double r_1 = disk.getRadius();

    double p   = directrix.getDistance(focus) / 2.0;
    double t   = p + r_1 - radiusOfFocus;
    double c   = x_1 * x_1 + y_1 * y_1 - t * t;
    double b   = -x_1;
    double a   = 1.0 - (y_1 + t) / 2.0 / p;

    double x_intersectionPt_1 = 0.0;
    double x_intersectionPt_2 = 0.0;

    if ( rg_ZERO( a, resNeg6 ) ) {
        x_intersectionPt_1 = c / -b / 2.0;
        x_intersectionPt_2 = x_intersectionPt_1;
    }
    else {
        double value_in_squared = b * b -  a * c;

        double squared_ = sqrt( value_in_squared );
        x_intersectionPt_1 = ( -b - squared_ ) / a;
        x_intersectionPt_2 = ( -b + squared_ ) / a;
    }

    rg_TMatrix2D backtransformMat;
    backtransformMat.rotate(angle);
    backtransformMat.translate(translatePt);    

    double y_intersectionPt_1 = x_intersectionPt_1 * x_intersectionPt_1 / 4.0 / p;
    double y_intersectionPt_2 = x_intersectionPt_2 * x_intersectionPt_2 / 4.0 / p;

    double r_intersectionPt_1 = y_intersectionPt_1 + p - radiusOfFocus;
    double r_intersectionPt_2 = y_intersectionPt_2 + p - radiusOfFocus;

    rg_Point2D centerPt[2] = {rg_Point2D(x_intersectionPt_1, y_intersectionPt_1), rg_Point2D(x_intersectionPt_2, y_intersectionPt_2)};
    centerPt[0] = backtransformMat * centerPt[0];
    centerPt[1] = backtransformMat * centerPt[1];

    tangentCircle1.setCenterPt(centerPt[0]);
    tangentCircle2.setCenterPt(centerPt[1]);
    tangentCircle1.setRadius(r_intersectionPt_1);
    tangentCircle2.setRadius(r_intersectionPt_2);
}



int PolygonVD2D::computeCoordOfNewVVertex_on_line_of_two_polygon_edges( const rg_Line2D& geometry_of_edge, const rg_Circle2D& circumcircle_sp, const rg_Circle2D& circumcircle_ep, const rg_Circle2D& disk, rg_Circle2D& tangentCircle1, rg_Circle2D& tangentCircle2 )
{
    rg_Point2D lineDirVec = geometry_of_edge.evaluateVector();
    lineDirVec = lineDirVec.getUnitVector();


    //rg_Point2D xAxisDirVec(1.0, 0.0);
    //double angle_from_positive_X_axis = angleFromVec1toVec2(xAxisDirVec, lineDirVec);

    //rg_TMatrix2D transformMat;
    //transformMat.translate(-geometry_of_edge.getSP());
    //transformMat.rotate(-angle_from_positive_X_axis);

    rg_Line2D bisector( circumcircle_sp.getCenterPt(), circumcircle_ep.getCenterPt() );
    rg_Point2D footprint_of_center_to_bisector;
    bisector.compute_perpendicular_footprint_of_point_onto_entire_line( disk.getCenterPt(), footprint_of_center_to_bisector );

    rg_Point2D SP = circumcircle_sp.getCenterPt();
    double x_1 = SP.distance( footprint_of_center_to_bisector );
    if ( !bisector.does_contain_in_line_segment( footprint_of_center_to_bisector ) ) {
        double distance_from_EP = circumcircle_ep.getCenterPt().distance( footprint_of_center_to_bisector );
        if ( distance_from_EP > x_1 ) {
            x_1 = -x_1;
        }
    }
    double y_1 = abs( bisector.signed_distance( disk.getCenterPt() ) );
    double r_1 = disk.getRadius();
    double clearance_sp = circumcircle_sp.getRadius() + r_1;
    double clearance_ep = circumcircle_ep.getRadius() + r_1;
    double diff = clearance_ep - clearance_sp;

    double a = diff * diff - 1.0;
    double half_b = diff * clearance_sp + x_1;
    double c = clearance_sp * clearance_sp - x_1 * x_1 - y_1 * y_1;
    double bb_ac = half_b * half_b - a * c;

    if ( bb_ac < 0.0 ) {
        return 0;
    }

    if ( rg_ZERO( a ) ) {
        double t = -c / half_b / 2.0;
        tangentCircle1.setCenterPt( SP + t * lineDirVec );
        tangentCircle1.setRadius( clearance_sp + t * diff - r_1 );

        return 1;
    }

    double      squareRoot = sqrt( bb_ac );
    if ( rg_ZERO( squareRoot ) ) {
        double  t = -half_b / a;
        tangentCircle1.setCenterPt( SP + t * lineDirVec );
        tangentCircle1.setRadius( clearance_sp + t * diff - r_1 );

        return 1;
    }

    double t_1 = ( -half_b + squareRoot ) / a;
    double t_2 = ( -half_b - squareRoot ) / a;
    tangentCircle1.setCenterPt( SP + t_1 * lineDirVec );
    tangentCircle1.setRadius( clearance_sp + t_1 * diff - r_1 );
    tangentCircle2.setCenterPt( SP + t_2 * lineDirVec );
    tangentCircle2.setRadius( clearance_sp + t_2 * diff - r_1 );

    return 2;





    /*
    rg_Point2D transformed_newPt = transformMat * disk.getCenterPt();
    double x_2 = geometry_of_edge.getLength();
    double x_3 = transformed_newPt.getX();
    double y_3 = transformed_newPt.getY();
    double r_3 = disk.getRadius();
    double d_1 = circumcircle_sp.getRadius();
    double d_2 = circumcircle_ep.getRadius();

    double t_1 = d_1 + r_3;
    double t_2 = d_2 - d_1;

    double a   = 1.0 - t_2 * t_2;
    double b   = -x_3 - t_1 * t_2;
    double c   = x_3 * x_3 + y_3 * y_3 - t_1 * t_1;

    double inside_square = b * b - a * c;

    if ( inside_square < 0.0 )
        return 0;

    double squared = sqrt(inside_square);
    double i_1 = (- b - squared) / a;
    double i_2 = ( -b + squared) / a;

    //double tangentCircle1_x = x_2 * i_1;
    //double tangentCircle2_x = x_2 * i_2;

    //rg_Point2D tangentCircle1_center( x_2 * i_1, 0.0 );
    //rg_Point2D tangentCircle2_center( x_2 * i_2, 0.0 );

    //rg_TMatrix2D backtransformMat;
    //backtransformMat.rotate(angle_from_positive_X_axis);
    //backtransformMat.translate(geometry_of_edge.getSP());

    //tangentCircle1_center = backtransformMat * tangentCircle1_center;
    //tangentCircle2_center = backtransformMat * tangentCircle2_center;

    rg_Point2D SP = circumcircle_sp.getCenterPt();
    rg_Point2D tangentCircle1_center = SP + lineDirVec * i_1;
    rg_Point2D tangentCircle2_center = SP + lineDirVec * i_2;

    double tangentCircle1_radius = d_1 + ( d_2 - d_1 ) * i_1;
    double tangentCircle2_radius = d_1 + ( d_2 - d_1 ) * i_2;

    tangentCircle1.setCenterPt(tangentCircle1_center);
    tangentCircle2.setCenterPt(tangentCircle2_center);
    tangentCircle1.setRadius(tangentCircle1_radius);
    tangentCircle2.setRadius(tangentCircle2_radius);

    return 2;
    */
}



void PolygonVD2D::computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex( const rg_Point2D& SP, const rg_Point2D& dirVec, const rg_Circle2D& disk, rg_Circle2D& tangentCircle )
{
    double x_s = SP.getX();
    double y_s = SP.getY();
    double x_v = dirVec.getX();
    double y_v = dirVec.getY();
    double x_d = disk.getCenterPt().getX();
    double y_d = disk.getCenterPt().getY();
    double r_d = disk.getRadius();

    double x_s_d = x_s - x_d;
    double y_s_d = y_s - y_d;

    // double a = x_v * x_v + y_v * y_v - 1.0;  // a is always 0 because x_v and y_v are elements of 
    double b = x_s_d * x_v + y_s_d * y_v - r_d;
    double c = x_s_d * x_s_d + y_s_d * y_s_d - r_d * r_d;

    double d = -c / b / 2.0;
    rg_Point2D tangentCircle_center = SP + d * dirVec;
    tangentCircle.setCenterPt( tangentCircle_center );
    tangentCircle.setRadius(d);
}



bool PolygonVD2D::findIntersectionPtWithinThisArrange( const rg_Line2D& line, const rg_Circle2D& disk_base, const rg_Circle2D& disk_another, const rg_Point2D& SP, const rg_Point2D& EP, rg_Circle2D& tangentCircle )
{
    rg_Circle2D tangentCircle_LP_base;
    rg_Circle2D tangentCircle_RP_base;
    rg_Circle2D tangentCircle_MP_base;
    rg_Circle2D tangentCircle_LP_another;
    rg_Circle2D tangentCircle_RP_another;
    rg_Circle2D tangentCircle_MP_another;

    rg_Point2D LP = SP;
    rg_Point2D RP = EP;
    rg_Point2D MP;

    rg_Point2D dirVec = line.getNormalVector().getUnitVector();
    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex(LP, dirVec, disk_base,    tangentCircle_LP_base);
    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex(LP, dirVec, disk_another, tangentCircle_LP_another);
    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex(RP, dirVec, disk_base,    tangentCircle_RP_base);
    computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex(RP, dirVec, disk_another, tangentCircle_RP_another);

    bool b_base_has_bigger_tangentCircle_LP = false;
    bool b_base_has_bigger_tangentCircle_MP = false;
    bool b_base_has_bigger_tangentCircle_RP = false;

    double radius_LP_base       = tangentCircle_LP_base.getRadius();
    double radius_MP_base;
    double radius_RP_base       = tangentCircle_RP_base.getRadius();
    double radius_LP_another    = tangentCircle_LP_another.getRadius();
    double radius_MP_another;
    double radius_RP_another    = tangentCircle_RP_another.getRadius();

    if ( rg_EQ( radius_LP_base, radius_LP_another ) ) {
        tangentCircle = tangentCircle_LP_base;
        return true;
    }
    if ( rg_EQ( radius_RP_base, radius_RP_another ) ) {
        tangentCircle = tangentCircle_RP_base;
        return true;
    }

    if ( radius_LP_base > radius_LP_another )
        b_base_has_bigger_tangentCircle_LP = true;
    if ( radius_RP_base > radius_RP_another )
        b_base_has_bigger_tangentCircle_RP = true;


    if ( b_base_has_bigger_tangentCircle_LP == b_base_has_bigger_tangentCircle_RP )
        return false;



    int MAX_ITERATION = 100;
    int i = 0;
    while ( i < MAX_ITERATION ) {
        MP = (LP + RP) / 2.0;
        computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex(MP, dirVec, disk_base,    tangentCircle_MP_base);
        computeCoordOfNewVVertex_on_line_of_incident_polygon_edge_N_polygon_vertex(MP, dirVec, disk_another, tangentCircle_MP_another);

        radius_MP_base      = tangentCircle_MP_base.getRadius();
        radius_MP_another   = tangentCircle_MP_another.getRadius();

        if ( rg_EQ( radius_MP_base, radius_MP_another ) ) {
            tangentCircle = tangentCircle_MP_base;
            return true;
        }

        if ( radius_MP_base > radius_MP_another )
            b_base_has_bigger_tangentCircle_MP = true;
        else {
            b_base_has_bigger_tangentCircle_MP = false;
        }

        if ( b_base_has_bigger_tangentCircle_LP == b_base_has_bigger_tangentCircle_MP ) {
            LP = MP;
            tangentCircle_LP_base       = tangentCircle_MP_base;
            tangentCircle_LP_another    = tangentCircle_MP_another;
            radius_LP_base              = radius_MP_base;
            radius_LP_another           = radius_MP_another;
        }
        else {
            RP = MP;
            tangentCircle_RP_base       = tangentCircle_MP_base;
            tangentCircle_RP_another    = tangentCircle_MP_another;
            radius_RP_base              = radius_MP_base;
            radius_RP_another           = radius_MP_another;
        }

        ++i;
    }

    tangentCircle = tangentCircle_MP_base;
    return true;
}



bool PolygonVD2D::this_circumcircle_has_right_orientation( const VVertex2D * const vertex, const rg_Circle2D & circumcircle )
{
    rg_Point2D vec_to_footprint[3];

    int index_footprint = 0;
    int index_generator = 0;


    list<Generator2D*> generators;
    vertex->getDefining3Generators(generators);
    for ( list<Generator2D*>::iterator i_gen = generators.begin(); i_gen != generators.end(); ++i_gen, ++index_footprint ) {
        Generator2D* currGen = *i_gen;
        Generator2D::Generator_Type type = currGen->getType();

        switch ( type )
        {
        case Generator2D::Generator_Type::EDGE_G:
        {
            EdgeGenerator2D* edgeGen = (EdgeGenerator2D*)currGen;
            rg_Line2D lineSeg = edgeGen->get_geometry();
            //lineSeg.compute_footprint_of_point_onto_line_segment(circumcircle.getCenterPt(), vec_to_footprint[index_footprint]);
            lineSeg.compute_perpendicular_footprint_of_point_onto_entire_line(circumcircle.getCenterPt(), vec_to_footprint[index_footprint]);

            for ( list<Generator2D*>::iterator j_gen = generators.begin(); j_gen != generators.end(); ++j_gen ) {
                Generator2D* currGen_to_check_adjacency = *j_gen;

                if ( currGen_to_check_adjacency == edgeGen ) {
                    continue;
                }

                if ( currGen_to_check_adjacency == edgeGen->get_start_vertex_generator() || currGen_to_check_adjacency == edgeGen->get_end_vertex_generator() ) {
                    vec_to_footprint[index_footprint] = ( lineSeg.getSP() + lineSeg.getEP() ) / 2.0;
                    break;
                }
            }
            vec_to_footprint[index_footprint] = ( vec_to_footprint[index_footprint] - circumcircle.getCenterPt() ).getUnitVector();
        }
        break;

        case Generator2D::Generator_Type::VERTEX_G:
        {
            VertexGenerator2D* vertexGen = (VertexGenerator2D*)currGen;
            vec_to_footprint[index_footprint] = vertexGen->get_point();
            vec_to_footprint[index_footprint] = ( vec_to_footprint[index_footprint] - circumcircle.getCenterPt() ).getUnitVector();
        }
        break;

        case Generator2D::Generator_Type::DISK_G:
        {
            DiskGenerator2D* diskGen = (DiskGenerator2D*)currGen;
            vec_to_footprint[index_footprint] = diskGen->getDisk().getCenterPt();
            vec_to_footprint[index_footprint] = ( vec_to_footprint[index_footprint] - circumcircle.getCenterPt() ).getUnitVector();

            double distanceBetweenCenters = diskGen->getDisk().getCenterPt().distance( circumcircle.getCenterPt() );
            if ( distanceBetweenCenters < diskGen->getDisk().getRadius() ) {
                vec_to_footprint[index_footprint] = -vec_to_footprint[index_footprint];
            }
        }
        break;

        default:
            break;
        }

        //vec_to_footprint[index_footprint] = vec_to_footprint[index_footprint] - circumcircle.getCenterPt();
    }

    double angle_1_to_2 = angleFromVec1toVec2( vec_to_footprint[0], vec_to_footprint[1] );
    double angle_1_to_3 = angleFromVec1toVec2( vec_to_footprint[0], vec_to_footprint[2] );

    if ( angle_1_to_2 < angle_1_to_3 ) {
        return true;
    }
    else {
        return false;
    }
}



rg_Point2D PolygonVD2D::getTangentVector( const rg_Point2D & pt, const rg_Point2D & leftCenter, const rg_Point2D & rightCenter )
{
    rg_Point2D vec1, vec2;
    vec1 = (leftCenter - pt).getUnitVector();
    vec2 = (rightCenter - pt).getUnitVector();

    if( rg_ZERO( (vec1+vec2).magnitude() ) )
    {
        return rg_Point2D( -vec1.getY(), vec1.getX() );
    }
    else
    {
        return vec1 + vec2;
    }
}

rg_Point2D PolygonVD2D::getPassingPt( const rg_Circle2D & leftDisk, const rg_Circle2D & rightDisk )
{
    rg_Circle2D c1, c2;	//c1.r > c2.r
    if( leftDisk.getRadius() > rightDisk.getRadius() )
    {
        c1 = leftDisk;
        c2 = rightDisk;
    }
    else
    {
        c2 = leftDisk;
        c1 = rightDisk;
    }

    rg_Point2D cp = ( c1.getCenterPt() + c2.getCenterPt() ) / 2.;
    rg_Point2D vector = (cp - c2.getCenterPt()).getUnitVector();

    rg_REAL A = c2.getCenterPt().getX() - c1.getCenterPt().getX();
    rg_REAL B = c2.getCenterPt().getY() - c1.getCenterPt().getY();
    rg_REAL r = c1.getRadius() - c2.getRadius();

    rg_REAL t = (r*r - A*A - B*B)/(A*vector.getX()+B*vector.getY() - r)/2.;

    //Is t always positive?
    // --> t should be nonnegative.
    rg_Point2D passingPt = c2.getCenterPt() + t * vector;

    return passingPt;
}

