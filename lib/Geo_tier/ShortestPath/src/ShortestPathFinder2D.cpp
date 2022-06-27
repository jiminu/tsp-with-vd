#include "ShortestPathFinder2D.h"
#include "VoronoiDiagram2DC.h"
#include "EntityAccessiblePriorityQ.h"
#include "ConstForShortestPathFinder.h"
#include "rg_Line.h"
#include "TPTNode.h"

#include "VoronoiDiagramCIC.h"
#include "rg_BoundingBox2D.h"

#include <unordered_set>
using namespace std;

//constructor
//
ShortestPathFinder2D::ShortestPathFinder2D()
{
    m_PathTotalDistanceByBFS           = -1.0;
    m_PathTotalDistanceByAllPathSearch = -1.0;
    m_PathTotalDistanceByDFS           = -1.0;
}

ShortestPathFinder2D::ShortestPathFinder2D(const ShortestPathFinder2D& SPF)
{
    copy_from(SPF);
}

ShortestPathFinder2D::~ShortestPathFinder2D()
{
}

void ShortestPathFinder2D::copy_from(const ShortestPathFinder2D& SPF)
{

}

//setter
//
void ShortestPathFinder2D::set_obstacles(const list<rg_Circle2D>& obstacles)
{
    m_Obstacles.clear();
    m_Obstacles = obstacles;
}


void ShortestPathFinder2D::set_input_argments(const rg_Circle2D& startPt,
                                              const rg_Circle2D& endPt, 
                                              const rg_Circle2D& probe, 
                                              const list<rg_Circle2D>& obstacles)
{
    set_start_pt(startPt);
    set_end_pt(endPt);
    set_probe(probe); 
    set_obstacles(obstacles);
}


void ShortestPathFinder2D::set_input_argments(const rg_Circle2D & startPt, const rg_Circle2D & endPt, const rg_Circle2D & probe, const list<rg_Circle2D>& obstacles, const double & probeMaxSpeed, const double & travelTime)
{
    set_start_pt(startPt);
    set_end_pt(endPt);
    set_probe(probe);
    set_obstacles(obstacles);
    set_probe_max_speed(probeMaxSpeed);
    set_max_travel_time(travelTime);
}


//getter
//
void ShortestPathFinder2D::get_obstacles(list<rg_Circle2D*>& obstacles)
{
    for (list<rg_Circle2D>::iterator it_obstacle = m_Obstacles.begin(); it_obstacle != m_Obstacles.end(); it_obstacle++)
    {
        rg_Circle2D* currObstacle = &(*it_obstacle);
        obstacles.push_back(currObstacle);
    }
}

void ShortestPathFinder2D::construct_voronoi_diagram(VoronoiDiagram2DC& VD2DC)
{
    VD2DC.constructVoronoiDiagramWithoutSorting(m_Obstacles);
}


void ShortestPathFinder2D::construct_quasi_triangulation_from_VD(QuasiTriangulation2D& QT2D, VoronoiDiagram2DC& VD2DC)
{
    QT2D.construct(VD2DC);
}


void ShortestPathFinder2D::construct_beta_shape_from_QT(BetaUniverse2D& betaUniverse, QuasiTriangulation2D& QT2D)
{
    betaUniverse.construct(QT2D);
}



//function
//

void ShortestPathFinder2D::compute_optimal_path(const rg_Circle2D& startPt, const rg_Circle2D& endPt, const rg_Circle2D& probe, const list<rg_Circle2D>& obstacles)
{
    set_input_argments(startPt, endPt, probe, obstacles);

    //2. construct Voronoi diaram and quasi-triangulation and beta-shape
    VoronoiDiagram2DC VD2DC;

    construct_voronoi_diagram(VD2DC);
    construct_quasi_triangulation_from_VD(m_QuasiTriangulation, VD2DC);
    construct_beta_shape_from_QT(m_BetaUniverse, m_QuasiTriangulation);

    find_connected_components_of_beta_shape(m_BetaUniverse, m_BetaComponents, m_BetaComponentFinder);

    // 3. mark qt edge whether in or out of ellipse
    compute_ellipse_filter(m_StartPt, m_EndPt, m_EllipseFilter);
    mark_QT_vertex_whether_in_or_out_of_ellipse(m_BetaUniverse, m_EllipseFilter, m_QTVertexValidation);
    mark_QT_edge_whether_in_or_out_of_ellipse(m_BetaUniverse, m_QTVertexValidation, m_QTEdgeValidation);


    compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS_ver1();
}


bool ShortestPathFinder2D::compute_geodesic_path(const rg_Circle2D& startPt, 
                                                 const rg_Circle2D& endPt, 
                                                 const rg_Circle2D& probe, 
                                                 const list<rg_Circle2D>& obstacles /**/)
{
    bool isThereaAGeodesicPath = false;

    if (startPt == endPt)
    {
        return false;
    }

    initialize_statistics();

    //1. set input
    ___start_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);

    compute_optimal_path(startPt, endPt, probe, obstacles);

    ___end_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);


    ___start_clock(SPF_TIME_FIND_PATH_BY_BFS);

    isThereaAGeodesicPath = compute_geodesic_path_based_on_VD(startPt, endPt, probe, obstacles);

    ___end_clock(SPF_TIME_FIND_PATH_BY_BFS);


    finalize_statistics();

    return isThereaAGeodesicPath;

    /*
    if (startPt == endPt)
    {
        return;
    }

    initialize_statistics();

    //1. set input
    set_input_argments(startPt, endPt, probe, obstacles);

    ___start_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);

    //2. construct Voronoi diaram and quasi-triangulation and beta-shape
    VoronoiDiagram2DC VD2DC;

    construct_voronoi_diagram(VD2DC);
    construct_quasi_triangulation_from_VD(m_QuasiTriangulation, VD2DC);
    construct_beta_shape_from_QT(m_BetaUniverse, m_QuasiTriangulation);
    
    find_connected_components_of_beta_shape(m_BetaUniverse, m_BetaComponents, m_BetaComponentFinder);

    ___end_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);




    ___start_clock(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE);

    // 3. mark qt edge whether in or out of ellipse
#ifdef ELLIPSE_FILTER
    compute_ellipse_filter(m_StartPt, m_EndPt, m_EllipseFilter);
    mark_QT_vertex_whether_in_or_out_of_ellipse(m_BetaUniverse, m_EllipseFilter, m_QTVertexValidation);
    mark_QT_edge_whether_in_or_out_of_ellipse(m_BetaUniverse, m_QTVertexValidation, m_QTEdgeValidation);
#endif
    ___end_clock(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE);




    ___start_clock(SPF_TIME_FIND_PATH_BY_BFS);

    //4. compute geodesic path by BFS
#ifdef ELLIPSE_FILTER
    //compute_geodesic_path_by_BFS_with_filtered_QT_edges();
    compute_geodesic_path_by_BFS_with_filtered_QT_edges_ver1();
#else
    compute_geodesic_path_by_BFS(m_BetaUniverse, m_PathByBFS);
#endif
    
    ___end_clock(SPF_TIME_FIND_PATH_BY_BFS);




    ___start_clock(SPF_TIME_FIND_PATH_BY_DFS);

#ifdef ELLIPSE_FILTER
    //compute_geodesic_path_by_DFS_with_filtered_QT_edges();
    compute_geodesic_path_by_DFS_with_filtered_QT_edges_ver1();
#else
    compute_geodesic_path_by_DFS(m_BetaUniverse, m_PathByBFS);
#endif

    ___end_clock(SPF_TIME_FIND_PATH_BY_DFS);




    ___start_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);

#ifdef ELLIPSE_FILTER
    //compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS();
    compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS_ver1();
#else
    compute_geodesic_path_by_all_path_search(m_BetaUniverse, m_PathByAllPathSearch);
#endif
    ___end_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);



    finalize_statistics();

    */
}



void ShortestPathFinder2D::compute_geodesic_path(const rg_Circle2D& startPt, 
                                                 const rg_Circle2D& endPt, 
                                                 const rg_Circle2D& probe, 
                                                 const double& travelTime, 
                                                 const double& probeMaxSpeed, 
                                                 const list<rg_Circle2D>& obstacles)
{
    if (startPt == endPt)
    {
        return;
    }

    initialize_statistics();

    //1. set input
    set_input_argments(startPt, endPt, probe, obstacles, probeMaxSpeed, travelTime);

    ___start_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);

    //2. construct Voronoi diaram and quasi-triangulation and beta-shape
    VoronoiDiagram2DC VD2DC;

    construct_voronoi_diagram(VD2DC);
    construct_quasi_triangulation_from_VD(m_QuasiTriangulation, VD2DC);
    construct_beta_shape_from_QT(m_BetaUniverse, m_QuasiTriangulation);
    
    find_connected_components_of_beta_shape(m_BetaUniverse, m_BetaComponents, m_BetaComponentFinder);

    ___end_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);


    ___start_clock(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE);

    // 3. mark qt edge whether in or out of ellipse
#ifdef ELLIPSE_FILTER
    compute_box_filter(m_StartPt, m_EndPt, m_MaxSpeedOfProbe, m_TravelTime, m_LineSegmentForBoxHeight, m_HalfLineSegmentForBoxWidth);
    mark_QT_vertex_whether_in_or_out_of_box(m_BetaUniverse, m_LineSegmentForBoxHeight, m_HalfLineSegmentForBoxWidth, m_QTVertexValidation);
    mark_QT_edge_whether_in_or_out_of_box(m_BetaUniverse, m_LineSegmentForBoxHeight, m_HalfLineSegmentForBoxWidth, m_QTVertexValidation, m_QTEdgeValidation);
    //compute_ellipse_filter(m_StartPt, m_EndPt, m_EllipseFilter);
    //mark_QT_vertex_whether_in_or_out_of_ellipse(m_BetaUniverse, m_EllipseFilter, m_QTVertexValidation);
    //mark_QT_edge_whether_in_or_out_of_ellipse(m_BetaUniverse, m_QTVertexValidatiovin, m_QTEdgeValidation);
#endif
    ___end_clock(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE);



    ___start_clock(SPF_TIME_FIND_PATH_BY_BFS);

    //4. compute geodesic path by BFS
#ifdef ELLIPSE_FILTER
    compute_geodesic_path_by_BFS_with_filtered_QT_edges();
#else
    compute_geodesic_path_by_BFS(m_BetaUniverse, m_PathByBFS);
#endif
    
    ___end_clock(SPF_TIME_FIND_PATH_BY_BFS);




    ___start_clock(SPF_TIME_FIND_PATH_BY_DFS);

#ifdef ELLIPSE_FILTER
    //compute_geodesic_path_by_DFS_with_filtered_QT_edges();
    compute_geodesic_path_by_DFS_with_filtered_QT_edges_ver1();
#else
    compute_geodesic_path_by_DFS(m_BetaUniverse, m_PathByBFS);
#endif

    ___end_clock(SPF_TIME_FIND_PATH_BY_DFS);




    ___start_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);

#ifdef ELLIPSE_FILTER
    compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS();
#else
    compute_geodesic_path_by_all_path_search(m_BetaUniverse, m_PathByAllPathSearch);
#endif
    ___end_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);



    finalize_statistics();
}



void ShortestPathFinder2D::compute_geodesic_path(const rg_Circle2D& startPt,
                                                 const rg_Circle2D& endPt,
                                                 const rg_Circle2D& probe,
                                                 const list<rg_Circle2D>& obstacles,
                                                 VoronoiDiagram2DC& VD2DC)
{

    if (startPt == endPt)
    {
        return;
    }


    initialize_statistics();


    //1. set input
    set_input_argments(startPt, endPt, probe, obstacles);



    ___start_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);

    VoronoiDiagram2DC VD = VD2DC;
    construct_quasi_triangulation_from_VD(m_QuasiTriangulation, VD);
    construct_beta_shape_from_QT(m_BetaUniverse, m_QuasiTriangulation);

    find_connected_components_of_beta_shape(m_BetaUniverse, m_BetaComponents, m_BetaComponentFinder);

    ___end_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);



    ___start_clock(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE);

    // 3. mark qt face whether in or out of ellipse
#ifdef ELLIPSE_FILTER
    compute_ellipse_filter(m_StartPt, m_EndPt, m_EllipseFilter);
    mark_QT_vertex_whether_in_or_out_of_ellipse(m_BetaUniverse, m_EllipseFilter, m_QTVertexValidation);
    mark_QT_edge_whether_in_or_out_of_ellipse(m_BetaUniverse, m_QTVertexValidation, m_QTEdgeValidation);
#endif

    ___end_clock(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE);



    ___start_clock(SPF_TIME_FIND_PATH_BY_BFS);

    //4. compute geodesic path by BFS
#ifdef ELLIPSE_FILTER
    compute_geodesic_path_by_BFS_with_filtered_QT_edges();
#else
    compute_geodesic_path_by_BFS(m_BetaUniverse, m_PathByBFS);
#endif

    ___end_clock(SPF_TIME_FIND_PATH_BY_BFS);




    ___start_clock(SPF_TIME_FIND_PATH_BY_DFS);

    //compute_geodesic_path_by_DFS(m_BetaUniverse, m_PathByDFS);

    ___end_clock(SPF_TIME_FIND_PATH_BY_DFS);




    ___start_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);

#ifdef ELLIPSE_FILTER
    compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS();
    //compute_geodesic_path_by_all_path_search_with_filtered_QT_edges(m_BetaUniverse, m_IsThisQTEdgeInsideOfEllipse, m_PathByAllPathSearch);
#else
    compute_geodesic_path_by_all_path_search(m_BetaUniverse, m_PathByAllPathSearch);
#endif
    ___end_clock(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH);



    finalize_statistics();

}




bool ShortestPathFinder2D::compute_geodesic_path_based_on_VD(const rg_Circle2D& startPt,
                                const rg_Circle2D& endPt,
                                const rg_Circle2D& probe,
                                const list<rg_Circle2D>& obstacles)
{
    set_input_argments(startPt, endPt, probe, obstacles);

    bool isThereAPath = false;

    rg_Line2D lineFromStartToEnd(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());

    if ( !is_there_straight_line_path_between_two_points(lineFromStartToEnd)) {
 
        switch (m_Obstacles.size())
        {
        case 1:
        {

            isThereAPath = find_shortest_path_when_only_one_disk_is_intersected_with_straight_line(
                m_Obstacles.front(), m_StartPt, m_EndPt, m_SolutionPoolGraphByBFS, m_PathByBFS, m_PathTotalDistanceByBFS);

            break;
        }

        default:
        {
            VoronoiDiagramCIC VD_of_obstacles;
            construct_voronoi_diagram_of_CIC(m_StartPt, m_EndPt, m_Obstacles, VD_of_obstacles);

            //   topological path by Dijkstra algorithm for whole V-edges
            isThereAPath = find_topogical_path_by_Dijkstra_in_VD(
                m_StartPt, m_EndPt, m_Probe, VD_of_obstacles, m_TopologicalPathByBFS_InVD);

            //   topological path by BF search
            //find_topogical_path_by_breadth_first_search_in_VD(
            //    m_StartPt, m_EndPt, m_Probe, VD_of_obstacles, m_TopologicalPathByBFS_InVD);

            ////////////////////////////////////////////////////////////
            //
            //  find shortest path using topological path in VD
            //
            if (isThereAPath)
            {
                isThereAPath = find_shortest_path_from_topological_path_using_VD(
                    VD_of_obstacles, m_TopologicalPathByBFS_InVD, m_SolutionPoolGraphByBFS, m_PathByBFS, m_PathTotalDistanceByBFS);
            }

            break;
        }
        }
    }
    else {
        rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(lineFromStartToEnd.getSP(), lineFromStartToEnd.getEP());

        find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(
                lineEqFromSourceToDestination,
                m_SolutionPoolGraphByBFS,
                m_PathByBFS, m_PathTotalDistanceByBFS);

        isThereAPath = true;
    }

    return isThereAPath;
}



bool ShortestPathFinder2D::find_topogical_path_by_Dijkstra_in_VD(
                const rg_Circle2D& startPt,
                const rg_Circle2D& endPt,
                const rg_Circle2D& probe,
                VoronoiDiagramCIC& VD_of_obstacles,
                list<VEdge2D*>& topologicalPath) const
{
    ShortestPathSolutionPoolGraph graphForDijkstra;
    map<VVertex2D*, VertexForSolutionPoolGraph*> vertexMap;
    map<VEdge2D*, EdgeForSolutionPoolGraph*>     edgeMap;

    list<VEdge2D*>  edgesInVD;
    VD_of_obstacles.getVoronoiEdges(edgesInVD);

    for (list<VEdge2D*>::iterator i_edge = edgesInVD.begin(); i_edge != edgesInVD.end(); ++i_edge) {
        VEdge2D*   currEdge = *i_edge;
        if (!currEdge->is_passable(0.0)) {
            continue;
        }

        VVertex2D* vertex[2] = { currEdge->getStartVertex(), currEdge->getEndVertex() };

        VertexForSolutionPoolGraph* vertexInGraph[2] = { rg_NULL, rg_NULL };
        for (int i = 0; i < 2; ++i) {
            map<VVertex2D*, VertexForSolutionPoolGraph*>::iterator it_vtx = vertexMap.find(vertex[i]);
            if (it_vtx == vertexMap.end()) {
                vertexInGraph[i] = graphForDijkstra.create_vertex(VertexForSolutionPoolGraph());
                vertexMap.insert(make_pair(vertex[i], vertexInGraph[i]));
            }
            else {
                vertexInGraph[i] = it_vtx->second;
            }
        }
        EdgeForSolutionPoolGraph*   edgeInGraph = graphForDijkstra.create_edge(EdgeForSolutionPoolGraph());
        edgeMap.insert(make_pair(currEdge, edgeInGraph));

        edgeInGraph->set_start_vtx(vertexInGraph[0]);
        edgeInGraph->set_end_vtx(vertexInGraph[1]);
        edgeInGraph->set_edge_length(vertex[0]->getLocation().distance(vertex[1]->getLocation()));
        edgeInGraph->set_user_data(currEdge);

        vertexInGraph[0]->add_line_edge(edgeInGraph);
        vertexInGraph[1]->add_line_edge(edgeInGraph);
    }


    //VVertex2D* startVtx = VD_of_obstacles.findClosestVVertexToPoint(startPt.getCenterPt());
    //VVertex2D* endVtx   = VD_of_obstacles.findClosestVVertexToPoint(endPt.getCenterPt());

    VVertex2D* startVtx = VD_of_obstacles.findVisibleClosestVVertexToPoint(startPt.getCenterPt());
    VVertex2D* endVtx   = VD_of_obstacles.findVisibleClosestVVertexToPoint(endPt.getCenterPt());
    if (startVtx == endVtx) {
        list<VVertex2D*> adjacentVVertices;
        startVtx->getAdjacent3VVertices(adjacentVVertices);

        double minDistance = DBL_MAX;
        for (list<VVertex2D*>::iterator i_vtx = adjacentVVertices.begin(); i_vtx != adjacentVVertices.end(); ++i_vtx) {
            VVertex2D* currVtx = (*i_vtx);

            double dist = currVtx->getCircumcircle().getCenterPt().distance(endPt.getCenterPt());
            if (dist < minDistance) {
                endVtx = currVtx;
                minDistance = dist;
            }
        }
    }


    bool      isThereATopologicalPath = true;

    VertexForSolutionPoolGraph* startVtxInGraph = NULL;
    VertexForSolutionPoolGraph* endVtxInGraph = NULL;

    map<VVertex2D*, VertexForSolutionPoolGraph*>::const_iterator it_StartVtxInGraph = vertexMap.find(startVtx);

    if (it_StartVtxInGraph != vertexMap.end())
    {
        startVtxInGraph = it_StartVtxInGraph->second;

        graphForDijkstra.set_start_vtx(startVtxInGraph);
    }
    else
    {
        isThereATopologicalPath = false;
    }

    map<VVertex2D*, VertexForSolutionPoolGraph*>::const_iterator it_EndVtxInGraph = vertexMap.find(endVtx);

    if (it_EndVtxInGraph != vertexMap.end())
    {
        endVtxInGraph  = it_EndVtxInGraph->second;

        graphForDijkstra.set_end_vtx(endVtxInGraph);

    }
    else
    {
        isThereATopologicalPath = false;

    }

    list<EdgeForSolutionPoolGraph*> pathByDijkstra;

    if (isThereATopologicalPath)
    {
        find_geodesic_path_by_dijkstra(graphForDijkstra);

        double    totalPathDistance;
        isThereATopologicalPath = back_trace_to_start_vtx_and_find_geodesic_path(startVtxInGraph, endVtxInGraph, pathByDijkstra, totalPathDistance);


        if (isThereATopologicalPath)
        {
            for (list<EdgeForSolutionPoolGraph*>::iterator it_edge = pathByDijkstra.begin(); it_edge != pathByDijkstra.end(); ++it_edge) {
                topologicalPath.push_back((VEdge2D*)((*it_edge)->user_data()));
            }
        }
    }


    return isThereATopologicalPath;
}



void ShortestPathFinder2D::find_topogical_path_by_breadth_first_search_in_VD(
            const rg_Circle2D& startPt,
            const rg_Circle2D& endPt,
            const rg_Circle2D& probe,
            VoronoiDiagramCIC& VD_of_obstacles,
            list<VEdge2D*>&    topologicalPath)
{
    ShortestPathSolutionPoolGraph graphForBFS;
    map<VVertex2D*, VertexForSolutionPoolGraph*> vertexMap;
    map<VEdge2D*, EdgeForSolutionPoolGraph*>     edgeMap;


    VVertex2D* startVtx = VD_of_obstacles.findClosestVVertexToPoint(startPt.getCenterPt());
    VVertex2D* endVtx   = VD_of_obstacles.findClosestVVertexToPoint(endPt.getCenterPt());
    if (startVtx == endVtx) {
        list<VVertex2D*> adjacentVVertices;
        startVtx->getAdjacent3VVertices(adjacentVVertices);

        double minDistance = DBL_MAX;
        for (list<VVertex2D*>::iterator i_vtx = adjacentVVertices.begin(); i_vtx != adjacentVVertices.end(); ++i_vtx) {
            VVertex2D* currVtx = (*i_vtx);

            double dist = currVtx->getCircumcircle().getCenterPt().distance(endPt.getCenterPt());
            if (dist < minDistance) {
                endVtx = currVtx;
                minDistance = dist;
            }
        }
    }


    VertexForSolutionPoolGraph* startVtxInGraph = graphForBFS.create_vertex(VertexForSolutionPoolGraph());
    startVtxInGraph->set_user_data(startVtx);
    startVtxInGraph->set_accumulate_length_from_source_vtx(0.0);

    VertexForSolutionPoolGraph* endVtxInGraph = graphForBFS.create_vertex(VertexForSolutionPoolGraph());
    endVtxInGraph->set_user_data(endVtx);
    endVtxInGraph->set_accumulate_length_from_source_vtx(0.0);

    vertexMap.insert(make_pair(startVtx, startVtxInGraph));
    vertexMap.insert(make_pair(endVtx,   endVtxInGraph));


    list< pair<VVertex2D*, VEdge2D*> > queue;
    list<VEdge2D*> edgesIncidentToStart;
    startVtx->getIncident3VEdges(edgesIncidentToStart);
    for (list<VEdge2D*>::iterator i_edge = edgesIncidentToStart.begin(); i_edge != edgesIncidentToStart.end(); ++i_edge) {
        queue.push_back(make_pair(startVtx, (*i_edge)));
    }


    while (!queue.empty()) {
        pair<VVertex2D*, VEdge2D*> currPath = queue.front();
        queue.pop_front();

        VVertex2D* prevVtx = currPath.first;
        VVertex2D* currVtx = currPath.second->getOppositVVertex(prevVtx);
        VEdge2D*   currEdge = currPath.second;

        if (edgeMap.find(currEdge) != edgeMap.end() || !currEdge->is_passable( 0.0 ) ) {
            continue;
        }


        VertexForSolutionPoolGraph* prevVtxInGraph = vertexMap.find(prevVtx)->second;

        EdgeForSolutionPoolGraph* edgeInGraph = graphForBFS.create_edge(EdgeForSolutionPoolGraph());
        edgeInGraph->set_start_vtx(   prevVtxInGraph);
        edgeInGraph->set_edge_length( prevVtx->getLocation().distance(currVtx->getLocation()) );
        edgeInGraph->set_user_data(currEdge);

        edgeMap.insert(make_pair(currEdge, edgeInGraph));


        map<VVertex2D*, VertexForSolutionPoolGraph*>::iterator it_currVtx = vertexMap.find(currVtx);
        if (it_currVtx == vertexMap.end()) {
            VertexForSolutionPoolGraph* currVtxInGraph = graphForBFS.create_vertex(VertexForSolutionPoolGraph());
            currVtxInGraph->set_user_data(currVtx);
            double pathLength = prevVtxInGraph->get_accumulate_length_from_source_vtx() + edgeInGraph->get_edge_length();
            currVtxInGraph->set_accumulate_length_from_source_vtx(pathLength);

            currVtxInGraph->set_prev_edge(edgeInGraph);
            edgeInGraph->set_end_vtx(currVtxInGraph);

            vertexMap.insert(make_pair(currVtx, currVtxInGraph));

            list<VEdge2D*> edgesIncidentToCurr;
            currVtx->getIncident3VEdges(edgesIncidentToCurr);
            for (list<VEdge2D*>::iterator i_edge = edgesIncidentToCurr.begin(); i_edge != edgesIncidentToCurr.end(); ++i_edge) {
                VEdge2D* nextEdge = (*i_edge);
                if ( nextEdge != currEdge) {
                    queue.push_back(make_pair(currVtx, nextEdge));
                }
            }
        }
        else {
            VertexForSolutionPoolGraph* currVtxInGraph = it_currVtx->second;
            double pathLengthByCurrEdge = prevVtxInGraph->get_accumulate_length_from_source_vtx() + edgeInGraph->get_edge_length();
            edgeInGraph->set_end_vtx(currVtxInGraph);


            if (currVtxInGraph->get_prev_edge() == rg_NULL){
                //  This case is for only end vertex.
                currVtxInGraph->set_accumulate_length_from_source_vtx(pathLengthByCurrEdge);
                currVtxInGraph->set_prev_edge(edgeInGraph);
            }
            else {
                double pathLengthByAccumulated = currVtxInGraph->get_accumulate_length_from_source_vtx();
                if (pathLengthByCurrEdge < pathLengthByAccumulated) {
                    currVtxInGraph->set_accumulate_length_from_source_vtx(pathLengthByCurrEdge);
                    currVtxInGraph->set_prev_edge(edgeInGraph);
                }
            }
        }
    }


    VertexForSolutionPoolGraph* currVtxInBackTracing = endVtxInGraph;
    do {
        EdgeForSolutionPoolGraph* edgeInBackTracing = currVtxInBackTracing->get_prev_edge();
        topologicalPath.push_front((VEdge2D*)edgeInBackTracing->user_data());
        currVtxInBackTracing = edgeInBackTracing->get_opposite_vtx(currVtxInBackTracing);
    } while (currVtxInBackTracing != startVtxInGraph);
    
}



bool ShortestPathFinder2D::find_shortest_path_when_only_one_disk_is_intersected_with_straight_line(const rg_Circle2D& obstacle, const rg_Circle2D& startPt, const rg_Circle2D& endPt, ShortestPathSolutionPoolGraph& solutionPoolGraph, list<EdgeForSolutionPoolGraph*>& geodesicPath, double& totalPathDistance)
{
    bool isThereAShortestPath = false;

    // 1. solution pool graph
    make_graph_of_solution_pool(obstacle, startPt, endPt, solutionPoolGraph);

    //2. dijkstra algorithm in solution pool graph
    find_geodesic_path_by_dijkstra(solutionPoolGraph);

    //3. back tracking
    isThereAShortestPath = back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(),
                                                                          solutionPoolGraph.get_end_vtx(),
                                                                          geodesicPath,
                                                                          totalPathDistance);

    return isThereAShortestPath;
}



void ShortestPathFinder2D::make_graph_of_solution_pool(const rg_Circle2D& obstacle, const rg_Circle2D& startPt, const rg_Circle2D& endPt, ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
     // 1. make grpah of solution pool
    rg_Circle2D sourcePoint[2] = { startPt, endPt };
    VertexForSolutionPoolGraph* vertexOfSource[2] = { solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph()),
                                                      solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph()) };
    solutionPoolGraph.set_start_vtx(vertexOfSource[0]);
    solutionPoolGraph.set_end_vtx(vertexOfSource[1]);


    vector<VertexForSolutionPoolGraph*> tangentPtsOnObstacle;

    for (int i = 0; i < 2; ++i)
    {
        vector<rg_ImplicitEquation> tangentLine;
        if (obstacle.getRadius() == 0.0) {
            //  this case does not happen.
            //tangentLine.resize(1);
        }
        else {
            tangentLine.resize(2);
            tangentLine[0].setDegree(1);
            tangentLine[1].setDegree(1);
            make_exterior_tangent_lines_of_two_circles(sourcePoint[i], obstacle, tangentLine[0], tangentLine[1]);
        }

        for (int i_line = 0; i_line < tangentLine.size(); ++i_line)
        {
            rg_Point2D tangentPoint[2];
            tangentPoint[0] = sourcePoint[i].getCenterPt();
            tangentPoint[1] = compute_tangent_point_between_line_and_circle(obstacle, tangentLine[i_line]);


            VertexForSolutionPoolGraph* vertexOnObstacle[2] = { vertexOfSource[i], rg_NULL };
            EdgeForSolutionPoolGraph*   edgeForTangentLine = rg_NULL;

            insert_tangent_line_into_graph(tangentPoint[0], tangentPoint[1], tangentLine[i_line],
                vertexOnObstacle[0], vertexOnObstacle[1], edgeForTangentLine,
                solutionPoolGraph);

            tangentPtsOnObstacle.push_back(vertexOnObstacle[1]);
        }

    }



    // 2. make graph of arcs
    int numTangentPoints = tangentPtsOnObstacle.size();

    vector<VertexForSolutionPoolGraph*> tangentPoint(1);
    tangentPoint.insert(tangentPoint.begin(), tangentPtsOnObstacle.begin(), tangentPtsOnObstacle.end());

    //  compute angles of tangent points and sort them
    for (int i = 0; i < numTangentPoints; ++i) 
    {
        rg_Point2D  vecToTangentPt = tangentPoint[i]->get_tangent_point_to_disk() - obstacle.getCenterPt();
        double      angle          = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToTangentPt);
        tangentPoint[i]->set_angle(angle);
    }
    sort(tangentPoint.begin(), --tangentPoint.end(), VertexForSolutionPoolGraph::AngleLess());
    tangentPoint[numTangentPoints] = tangentPoint[0];

    //  find arc as solution
    for (int i = 0; i < numTangentPoints; ++i) 
    {
        //  make arc from tangentPoint[i] to tangentPoint[i+1]
        EdgeForSolutionPoolGraph* edgeForArc = rg_NULL;
        insert_arc_into_graph(obstacle, tangentPoint[i], tangentPoint[i + 1],
                              edgeForArc, solutionPoolGraph);
    }


}



bool ShortestPathFinder2D::find_shortest_path_from_topological_path_using_VD(
                    VoronoiDiagramCIC& VD_of_obstacles,
                    const list<VEdge2D*>& topologicalPath,
                    ShortestPathSolutionPoolGraph & solutionPoolGraph,
                    list<EdgeForSolutionPoolGraph*>& geodesicPath,
                    double & totalPathDistacne) const
{
    bool isThereAShortestPath = false;

    //1. Make obstacle pool 
    vector<VFace2D*> obstaclePool;
    make_obstacle_pool(topologicalPath, obstaclePool);


    //2. generate solution pool
    make_graph_of_solution_pool(VD_of_obstacles, obstaclePool, solutionPoolGraph);

    //3. dijkstra algorithm in solution pool graph
    find_geodesic_path_by_dijkstra(solutionPoolGraph);

    //4. back tracking
    isThereAShortestPath = back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(),
                                                                          solutionPoolGraph.get_end_vtx(),
                                                                          geodesicPath,
                                                                          totalPathDistacne);
    return isThereAShortestPath;
}



void    ShortestPathFinder2D::make_obstacle_pool(
                const list<VEdge2D*>&   topologicalPath, 
                vector<VFace2D*>&       obstaclePool) const
{
    set<VFace2D*> facesForObstaclePool;
    for (list<VEdge2D*>::const_iterator i_edge = topologicalPath.begin(); i_edge != topologicalPath.end(); ++i_edge) {
        VEdge2D* currEdge = *i_edge;

        VFace2D* faceIncidentToCurrEdge[4] = { 
                    currEdge->getLeftFace(), currEdge->getRightFace(), 
                    currEdge->getMateFace(currEdge->getStartVertex()), currEdge->getMateFace(currEdge->getEndVertex()) };

        for (int i = 0; i < 4; ++i) {
            if (faceIncidentToCurrEdge[i]->getID() != -1) {
                facesForObstaclePool.insert(faceIncidentToCurrEdge[i]);
            }
        }
    }

    obstaclePool.insert(obstaclePool.begin(), facesForObstaclePool.begin(), facesForObstaclePool.end());
}



void    ShortestPathFinder2D::make_graph_of_solution_pool(
                VoronoiDiagramCIC&              VD_of_obstacles,
                const vector<VFace2D*>&         obstaclePool,
                ShortestPathSolutionPoolGraph & solutionPoolGraph) const
{
    unsigned int numObstacles = obstaclePool.size();

    //1. Generate solution pool
    unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>> mapObstacleToTangentPoints;
    for (int i = 0; i < numObstacles; ++i) {
        mapObstacleToTangentPoints.insert(make_pair(obstaclePool[i], list<VertexForSolutionPoolGraph*>()));
    }

    
    //  1. Insert tangent lines between two obstacles into the solution pool
    make_tangent_lines_between_two_obstacles_and_insert_them_into_solution_pool(
        VD_of_obstacles, obstaclePool, solutionPoolGraph, mapObstacleToTangentPoints);


    //  2. Insert tangent lines from start and end to obstacles into the solution pool
    make_tangent_lines_between_start_end_points_and_obstacles_and_insert_them_into_solution_pool(
        VD_of_obstacles, obstaclePool, solutionPoolGraph, mapObstacleToTangentPoints);


    //  3. Insert arcs into the solution pool
    make_arcs_on_obstacles_and_insert_them_into_solution_pool(
        obstaclePool, solutionPoolGraph, mapObstacleToTangentPoints);

}



void ShortestPathFinder2D::make_tangent_lines_between_two_obstacles_and_insert_them_into_solution_pool(
                VoronoiDiagramCIC& VD_of_obstacles,
                const vector<VFace2D*>& obstaclePool,
                ShortestPathSolutionPoolGraph& solutionPoolGraph,
                unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const
{
    unsigned int numObstacles = obstaclePool.size();

    list<const VFace2D*> wholeObstacles;
    VD_of_obstacles.getVoronoiFaces(wholeObstacles);

    

    VFace2D*    obstacle[2] = { rg_NULL, rg_NULL };
    rg_Circle2D obstacle_geometry[2];
    for (int i = 0; i < (numObstacles - 1); ++i) {
        obstacle[0] = obstaclePool[i];
        obstacle_geometry[0] = obstaclePool[i]->getGenerator()->getDisk();

        for (int j = i + 1; j < numObstacles; ++j) {
            obstacle[1] = obstaclePool[j];
            obstacle_geometry[1] = obstaclePool[j]->getGenerator()->getDisk();

            vector<rg_ImplicitEquation> tangentLine;
            if (obstacle_geometry[0].getRadius() == 0.0 && obstacle_geometry[1].getRadius() == 0.0) {
                //  this case does not happen.
                //tangentLine.resize(1);
            }
            else if (obstacle_geometry[0].getRadius() == 0.0) {
                tangentLine.resize(2);
                tangentLine[0].setDegree(1);
                tangentLine[1].setDegree(1);
                make_exterior_tangent_lines_of_two_circles(obstacle_geometry[0], obstacle_geometry[1], tangentLine[0], tangentLine[1]);

            }
            else if (obstacle_geometry[1].getRadius() == 0.0) {
                tangentLine.resize(2);
                tangentLine[0].setDegree(1);
                tangentLine[1].setDegree(1);
                make_exterior_tangent_lines_of_two_circles(obstacle_geometry[0], obstacle_geometry[1], tangentLine[0], tangentLine[1]);
            }
            else {
                if (obstacle_geometry[0].isIntersectWith(obstacle_geometry[1])) {
                    tangentLine.resize(2);
                    tangentLine[0].setDegree(1);
                    tangentLine[1].setDegree(1);
                    make_exterior_tangent_lines_of_two_circles(obstacle_geometry[0], obstacle_geometry[1], tangentLine[0], tangentLine[1]);
                }
                else {
                    tangentLine.resize(4);
                    tangentLine[0].setDegree(1);
                    tangentLine[1].setDegree(1);
                    tangentLine[2].setDegree(1);
                    tangentLine[3].setDegree(1);
                    make_exterior_tangent_lines_of_two_circles(obstacle_geometry[0], obstacle_geometry[1], tangentLine[0], tangentLine[1]);
                    make_interior_tangent_lines_of_two_circles(obstacle_geometry[0], obstacle_geometry[1], tangentLine[2], tangentLine[3]);
                }
            }



            for (int i_line = 0; i_line < tangentLine.size(); ++i_line) {
                rg_Point2D tangentPoint[2];
                tangentPoint[0] = compute_tangent_point_between_line_and_circle(obstacle_geometry[0], tangentLine[i_line]);
                tangentPoint[1] = compute_tangent_point_between_line_and_circle(obstacle_geometry[1], tangentLine[i_line]);

                rg_Line2D tangentLineSement(tangentPoint[0], tangentPoint[1]);
                bool      isThisLineIntersectedWithObstacle = false;

                for (list<const VFace2D*>::iterator it_obstacle = wholeObstacles.begin(); it_obstacle != wholeObstacles.end(); ++it_obstacle) {
                    VFace2D* currObstacle = (VFace2D*)(*it_obstacle);
                    if (currObstacle->getID() == -1 || currObstacle->isInfinite() ) continue;

                    if (currObstacle == obstaclePool[i] || currObstacle == obstaclePool[j]) {
                        continue;
                    }

                    rg_Circle2D otherDisk = currObstacle->getGenerator()->getDisk();
                    if (otherDisk.hasIntersectionWith(tangentLineSement)) {
                        isThisLineIntersectedWithObstacle = true;
                        break;
                    }
                }


                if (isThisLineIntersectedWithObstacle) {
                    continue;
                }
                else {
                    VertexForSolutionPoolGraph* vertexOnObstacle[2] = { rg_NULL, rg_NULL };
                    EdgeForSolutionPoolGraph*   edgeForTangentLine = rg_NULL;

                    insert_tangent_line_into_graph(tangentPoint[0], tangentPoint[1], tangentLine[i_line],
                        vertexOnObstacle[0], vertexOnObstacle[1], edgeForTangentLine,
                        solutionPoolGraph);

                    mapObstacleToTangentPoints.find(obstacle[0])->second.push_back(vertexOnObstacle[0]);
                    mapObstacleToTangentPoints.find(obstacle[1])->second.push_back(vertexOnObstacle[1]);
                }
            }
        }
    }
}



void ShortestPathFinder2D::make_tangent_lines_between_start_end_points_and_obstacles_and_insert_them_into_solution_pool(
                VoronoiDiagramCIC& VD_of_obstacles,
                const vector<VFace2D*>& obstaclePool,
                ShortestPathSolutionPoolGraph& solutionPoolGraph,
                unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const
{
    unsigned int numObstacles = obstaclePool.size();

    list<const VFace2D*> wholeObstacles;
    VD_of_obstacles.getVoronoiFaces(wholeObstacles);


    rg_Circle2D sourcePoint[2] = { m_StartPt, m_EndPt };
    VertexForSolutionPoolGraph* vertexOfSource[2] = { solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph()),
                                                      solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph()) };
    solutionPoolGraph.set_start_vtx(vertexOfSource[0]);
    solutionPoolGraph.set_end_vtx(  vertexOfSource[1]);


    for (int i = 0; i < 2; ++i) {

        for (int j = 0; j < numObstacles; ++j) {
            VFace2D*    obstacle          = obstaclePool[j];
            rg_Circle2D obstacle_geometry = obstaclePool[j]->getGenerator()->getDisk();

            vector<rg_ImplicitEquation> tangentLine;
            if (obstacle_geometry.getRadius() == 0.0) {
                //  this case does not happen.
                //tangentLine.resize(1);
            }
            else { 
                tangentLine.resize(2);
                tangentLine[0].setDegree(1);
                tangentLine[1].setDegree(1);
                make_exterior_tangent_lines_of_two_circles(sourcePoint[i], obstacle_geometry, tangentLine[0], tangentLine[1]);
            }



            for (int i_line = 0; i_line < tangentLine.size(); ++i_line) {
                rg_Point2D tangentPoint[2];
                tangentPoint[0] = sourcePoint[i].getCenterPt();
                tangentPoint[1] = compute_tangent_point_between_line_and_circle(obstacle_geometry, tangentLine[i_line]);

                rg_Line2D tangentLineSement(tangentPoint[0], tangentPoint[1]);
                bool      isThisLineIntersectedWithObstacle = false;


                for (list<const VFace2D*>::iterator it_obstacle = wholeObstacles.begin(); it_obstacle != wholeObstacles.end(); ++it_obstacle) {
                    VFace2D* currObstacle = (VFace2D*)(*it_obstacle);
                    if (currObstacle->getID() == -1 || currObstacle->isInfinite()) continue;

                    if (currObstacle == obstaclePool[j]) {
                        continue;
                    }

                    rg_Circle2D otherDisk = currObstacle->getGenerator()->getDisk();
                    if (otherDisk.hasIntersectionWith(tangentLineSement)) {
                        isThisLineIntersectedWithObstacle = true;
                        break;
                    }
                }


                if (isThisLineIntersectedWithObstacle) {
                    continue;
                }
                else {
                    VertexForSolutionPoolGraph* vertexOnObstacle[2] = { vertexOfSource[i], rg_NULL };
                    EdgeForSolutionPoolGraph*   edgeForTangentLine = rg_NULL;

                    insert_tangent_line_into_graph(tangentPoint[0],     tangentPoint[1],     tangentLine[i_line],
                                                   vertexOnObstacle[0], vertexOnObstacle[1], edgeForTangentLine,
                                                   solutionPoolGraph);

                    mapObstacleToTangentPoints.find(obstacle)->second.push_back(vertexOnObstacle[1]);
                }
            }
        }
    }
}



void ShortestPathFinder2D::insert_tangent_line_into_graph(
    const rg_Point2D& ptOnObstacle1,
    const rg_Point2D& ptOnObstacle2,
    const rg_ImplicitEquation& tangentLineEq,
    VertexForSolutionPoolGraph*& vertexOnObstacle1,
    VertexForSolutionPoolGraph*& vertexOnObstacle2,
    EdgeForSolutionPoolGraph*& edgeForTangentLine,
    ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    if (vertexOnObstacle1 == rg_NULL) {
        vertexOnObstacle1 = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
    }

    if (vertexOnObstacle2 == rg_NULL) {
        vertexOnObstacle2 = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
    }
    edgeForTangentLine = solutionPoolGraph.create_edge(EdgeForSolutionPoolGraph());

    //2. Set two vertices

    vertexOnObstacle1->set_tangent_point_to_disk(ptOnObstacle1);
    vertexOnObstacle1->add_line_edge(edgeForTangentLine);
    vertexOnObstacle1->set_accumulate_length_from_source_vtx(DBL_MAX);
    vertexOnObstacle1->set_angle(-1.0);

    vertexOnObstacle2->set_tangent_point_to_disk(ptOnObstacle2);
    vertexOnObstacle2->add_line_edge(edgeForTangentLine);
    vertexOnObstacle2->set_accumulate_length_from_source_vtx(DBL_MAX);
    vertexOnObstacle2->set_angle(-1.0);

    //2. Set a tangent line edge
    double tangentLineSementLength = ptOnObstacle1.distance(ptOnObstacle2);

    edgeForTangentLine->set_start_vtx(vertexOnObstacle1);
    edgeForTangentLine->set_end_vtx(vertexOnObstacle2);
    edgeForTangentLine->set_edge_type(LINE_SEGMENT_SPG);
    edgeForTangentLine->set_edge_length(tangentLineSementLength);
    edgeForTangentLine->set_implicit_equation(tangentLineEq);
}

/*

void ShortestPathFinder2D::make_arcs_on_obstacles_and_insert_them_into_solution_pool(
                const vector<VFace2D*>& obstaclePool,
                ShortestPathSolutionPoolGraph& solutionPoolGraph,
                unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const
{
    unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>::iterator i_obstacle;
    for (i_obstacle = mapObstacleToTangentPoints.begin(); i_obstacle != mapObstacleToTangentPoints.end(); ++i_obstacle) {
        VFace2D*     obstacle = i_obstacle->first;

        int numTangentPoints = i_obstacle->second.size();
        if (numTangentPoints <= 1) {
            continue;
        }

        rg_Circle2D  obstacle_geometry = obstacle->getGenerator()->getDisk();
        rg_Point2D   center = obstacle_geometry.getCenterPt();

        vector<VertexForSolutionPoolGraph*> tangentPoint(1);
        tangentPoint.insert(tangentPoint.begin(), i_obstacle->second.begin(), i_obstacle->second.end());
        
        //  compute angles of tangent points and sort them
        for (int i = 0; i < numTangentPoints; ++i) {
            rg_Point2D  vecToTangentPt  = tangentPoint[i]->get_tangent_point_to_disk() - center;
            double      angle           = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToTangentPt);
            tangentPoint[i]->set_angle(angle);
        }
        sort(tangentPoint.begin(), --tangentPoint.end(), VertexForSolutionPoolGraph::AngleLess());
        tangentPoint[numTangentPoints] = tangentPoint[0];

        vector<double> angleToTangentPoints(numTangentPoints + 1);
        for (int i = 0; i < numTangentPoints; ++i) {
            angleToTangentPoints[i] = tangentPoint[i]->get_angle();
        }
        angleToTangentPoints[numTangentPoints] = tangentPoint[0]->get_angle() + 2*rg_PI;

        //  find first neighbors of obstacle
        list<VFace2D*> firstNeighbors;
        obstacle->getAdjacentVFaces(firstNeighbors);

        vector<double>  angleToFirstNeighborIntersectingObstacle;
        for (list<VFace2D*>::iterator i_neighbor = firstNeighbors.begin(); i_neighbor != firstNeighbors.end(); ++i_neighbor) {
            rg_Circle2D  neighbor_geometry = (*i_neighbor)->getGenerator()->getDisk();
            if (obstacle_geometry.isIntersectWith(neighbor_geometry)) {
                rg_Point2D  vecToNeighbor = neighbor_geometry.getCenterPt() - center;
                double      angle = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToNeighbor);
                if (angle < angleToTangentPoints[0]) {
                    angle += 2 * rg_PI;
                }
                angleToFirstNeighborIntersectingObstacle.push_back(angle);
            }
        }
        sort(angleToFirstNeighborIntersectingObstacle.begin(), angleToFirstNeighborIntersectingObstacle.end());


        //  find arc as solution
        int numIntersections = angleToFirstNeighborIntersectingObstacle.size();
        for (int i = 0; i < numTangentPoints; ++i) {
            bool    isArcAvailable = true;
            for (int j = 0; j < numIntersections; ++j) {
                if (   angleToTangentPoints[i] < angleToFirstNeighborIntersectingObstacle[j]
                    && angleToTangentPoints[i+1] > angleToFirstNeighborIntersectingObstacle[j]) {
                    isArcAvailable = false;
                    break;
                }
            }

            if (isArcAvailable == false) {
                continue;
            }


            //  make arc from tangentPoint[i] to tangentPoint[i+1]
            EdgeForSolutionPoolGraph* edgeForArc = rg_NULL;
            insert_arc_into_graph( obstacle_geometry, tangentPoint[i], tangentPoint[i + 1],
                                   edgeForArc, solutionPoolGraph);

        }
    }

}
*/


void ShortestPathFinder2D::make_arcs_on_obstacles_and_insert_them_into_solution_pool(
                                                                            const vector<VFace2D*>& obstaclePool,
                                                                            ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                            unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const
{
    unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>::iterator i_obstacle;
    for (i_obstacle = mapObstacleToTangentPoints.begin(); i_obstacle != mapObstacleToTangentPoints.end(); ++i_obstacle) {
        VFace2D*     obstacle          = i_obstacle->first;
        rg_Circle2D  obstacle_geometry = obstacle->getGenerator()->getDisk();
        rg_Point2D   center = obstacle_geometry.getCenterPt();

        if (i_obstacle->second.size() <= 1)
        {
            continue;
        }

        vector<VertexForSolutionPoolGraph*> tangentPoints;
        make_sorted_unique_set_of_vertices(i_obstacle->second, tangentPoints, obstacle);

        if (tangentPoints.size() <= 1)
        {
            continue;
        }
        tangentPoints.push_back(tangentPoints[0]);


        int numOfTangentPoints = tangentPoints.size();
       
        vector<double> angleToTangentPoints(numOfTangentPoints);
        for (int i = 0; i < (numOfTangentPoints-1); ++i) {
            angleToTangentPoints[i] = tangentPoints[i]->get_angle();
        }
        angleToTangentPoints[numOfTangentPoints - 1] = tangentPoints[numOfTangentPoints - 1]->get_angle() + 2 * rg_PI;



        //  find first neighbors of obstacle
        list<VFace2D*> firstNeighbors;
        obstacle->getAdjacentVFaces(firstNeighbors);

        vector<double>  angleToFirstNeighborIntersectingObstacle;
        for (list<VFace2D*>::iterator i_neighbor = firstNeighbors.begin(); i_neighbor != firstNeighbors.end(); ++i_neighbor) {
            rg_Circle2D  neighbor_geometry = (*i_neighbor)->getGenerator()->getDisk();
            if (obstacle_geometry.isIntersectWith(neighbor_geometry)) {
                rg_Point2D  vecToNeighbor = neighbor_geometry.getCenterPt() - center;
                double      angle = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToNeighbor);
                if (angle < angleToTangentPoints[0]) {
                    angle += 2 * rg_PI;
                }
                angleToFirstNeighborIntersectingObstacle.push_back(angle);
            }
        }
        sort(angleToFirstNeighborIntersectingObstacle.begin(), angleToFirstNeighborIntersectingObstacle.end());


        //  find arc as solution
        int numIntersections = angleToFirstNeighborIntersectingObstacle.size();
        for (int i = 0; i < (numOfTangentPoints-1); ++i) {
            bool    isArcAvailable = true;
            for (int j = 0; j < numIntersections; ++j) {
                if (angleToTangentPoints[i] < angleToFirstNeighborIntersectingObstacle[j]
                    && angleToTangentPoints[i + 1] > angleToFirstNeighborIntersectingObstacle[j]) {
                    isArcAvailable = false;
                    break;
                }
            }

            if (isArcAvailable == false) {
                continue;
            }


            //  make arc from tangentPoint[i] to tangentPoint[i+1]
            EdgeForSolutionPoolGraph* edgeForArc = rg_NULL;
            insert_arc_into_graph(obstacle_geometry, tangentPoints[i], tangentPoints[i + 1],
                edgeForArc, solutionPoolGraph);

        }
    }

}


void ShortestPathFinder2D::make_sorted_unique_set_of_vertices(const list<VertexForSolutionPoolGraph*>& vertices, 
                                                              vector<VertexForSolutionPoolGraph*>& uniqueVertices,
                                                              VFace2D* obstacle) const
{
    vector<VertexForSolutionPoolGraph*> copiedVertices;
    copiedVertices.insert(copiedVertices.end(), vertices.begin(), vertices.end());

    rg_Point2D   center = obstacle->getGenerator()->getDisk().getCenterPt();

    //  compute angles of tangent points and sort them
    for (vector<VertexForSolutionPoolGraph*>::const_iterator it_Vertex = copiedVertices.begin(); it_Vertex != copiedVertices.end(); ++it_Vertex)
    {
        VertexForSolutionPoolGraph* currVertex = *it_Vertex;

        rg_Point2D  vecToTangentPt = currVertex->get_tangent_point_to_disk() - center;
        double      angle = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToTangentPt);
        currVertex->set_angle(angle);
    }

    sort(copiedVertices.begin(), copiedVertices.end(), VertexForSolutionPoolGraph::AngleLess());
    copiedVertices.push_back(copiedVertices[0]);


    int  numOfVertices = copiedVertices.size();
    int  firstSamePtIndex = 0;

    for (int i = 0; i < (numOfVertices - 1) ; ++i)
    {
        VertexForSolutionPoolGraph* currVertex = copiedVertices[i];
        VertexForSolutionPoolGraph* nextVertex = copiedVertices[i+1];

        if (!VertexForSolutionPoolGraph::are_two_same_point(currVertex, nextVertex))
        {
            VertexForSolutionPoolGraph* remainingVtx = copiedVertices[firstSamePtIndex];

            for (int j = firstSamePtIndex + 1; j <= i; ++j)
            {
                VertexForSolutionPoolGraph* deletingVtx = copiedVertices[j];
                
                if (remainingVtx == deletingVtx)
                {
                    continue;
                }

                list<EdgeForSolutionPoolGraph*> tangentLineOfDeletingVtx;
                deletingVtx->get_line_edges(tangentLineOfDeletingVtx);

                for (list<EdgeForSolutionPoolGraph*>::const_iterator it_Edge = tangentLineOfDeletingVtx.begin(); it_Edge != tangentLineOfDeletingVtx.end(); ++it_Edge)
                {
                    EdgeForSolutionPoolGraph* currEdge = *it_Edge;

                    remainingVtx->add_line_edge(currEdge);
                    
                    //VertexForSolutionPoolGraph::are_two_same_point(currVertex, nextVertex)
                    if (VertexForSolutionPoolGraph::are_two_same_point(currEdge->get_start_vtx(), remainingVtx))
                    {
                        currEdge->set_start_vtx(remainingVtx);
                    }
                    else //currEdge->get_start_vtx() == remainingVtx
                    {
                        currEdge->set_end_vtx(remainingVtx);
                    }
                }

                deletingVtx->clear();
            }

            uniqueVertices.push_back(remainingVtx);

            firstSamePtIndex = i + 1;
        }
    }

    if (firstSamePtIndex != (numOfVertices - 1))
    {
        VertexForSolutionPoolGraph* remainingVtx = copiedVertices[firstSamePtIndex];

        for (int j = firstSamePtIndex + 1; j < numOfVertices; ++j)
        {
            VertexForSolutionPoolGraph* deletingVtx = copiedVertices[j];

            if (remainingVtx == deletingVtx)
            {
                continue;
            }

            list<EdgeForSolutionPoolGraph*> tangentLineOfDeletingVtx;
            deletingVtx->get_line_edges(tangentLineOfDeletingVtx);

            for (list<EdgeForSolutionPoolGraph*>::const_iterator it_Edge = tangentLineOfDeletingVtx.begin(); it_Edge != tangentLineOfDeletingVtx.end(); ++it_Edge)
            {
                remainingVtx->add_line_edge(*it_Edge);
            }

            deletingVtx->clear();
        }

        uniqueVertices.push_back(remainingVtx);
    }
}


void    ShortestPathFinder2D::insert_arc_into_graph(
                        const rg_Circle2D&              obstacle,
                        VertexForSolutionPoolGraph*     vertexOfStartPt,
                        VertexForSolutionPoolGraph*     vertexOfEndPt,
                        EdgeForSolutionPoolGraph*&      edgeForArc,
                        ShortestPathSolutionPoolGraph&  solutionPoolGraph) const
{
    edgeForArc = solutionPoolGraph.create_edge(EdgeForSolutionPoolGraph());

    rg_ImplicitEquation arcInplicitEquation = make_implicit_equation_of_arc_edge(obstacle, 0.0);
    double              arcLength = compute_arc_segment_length_on_obstacle_boundary(
                                        obstacle, vertexOfStartPt->get_tangent_point_to_disk(), vertexOfEndPt->get_tangent_point_to_disk());
    
    edgeForArc->set_start_vtx(vertexOfStartPt);
    edgeForArc->set_end_vtx(vertexOfEndPt);
    edgeForArc->set_edge_type(ARC_SPG);
    edgeForArc->set_edge_length(arcLength);
    edgeForArc->set_implicit_equation(arcInplicitEquation);

    vertexOfStartPt->add_arc_edge(edgeForArc);
    vertexOfEndPt->add_arc_edge(edgeForArc);
}


bool ShortestPathFinder2D::is_there_straight_line_path_between_two_points(const rg_Line2D& straightLinePath) const
{
    bool isThereStraightLinePath = true;

    for (list<rg_Circle2D>::const_iterator it_Obstacle = m_Obstacles.begin(); it_Obstacle != m_Obstacles.end(); ++it_Obstacle) {
        rg_Circle2D obstacle = (*it_Obstacle);
        if (obstacle.hasIntersectionWith(straightLinePath)) {
            isThereStraightLinePath = false;
            break;
        }
    }

    return isThereStraightLinePath;
}


void ShortestPathFinder2D::construct_voronoi_diagram_of_CIC(const rg_Circle2D& startPt, const rg_Circle2D& endPt, const list<rg_Circle2D>& obstacles, VoronoiDiagramCIC& VD_of_obstacles) const
{
    //  construct VD for CIC
    //  1) make container for VD_CIC
    list<rg_Circle2D> obstaclesForCIC = obstacles;
    rg_BoundingBox2D  boundingBoxOfObstacles;
    for (list<rg_Circle2D>::const_iterator i_obstacle = obstacles.begin(); i_obstacle != obstacles.end(); ++i_obstacle) {
        boundingBoxOfObstacles.updateBoxByAddingCircle(*i_obstacle);
    }
    boundingBoxOfObstacles.updateBoxByAddingCircle(startPt);
    boundingBoxOfObstacles.updateBoxByAddingCircle(endPt);
    rg_Circle2D container(boundingBoxOfObstacles.getCenterPt(), boundingBoxOfObstacles.evaluateLongestLength()/2.0);

    obstaclesForCIC.push_front(container);

    //  2) construct VD
    VD_of_obstacles.constructVoronoiDiagramCIC(obstaclesForCIC);

    //  3) set radius interval of edge
    list<VEdge2D*> edges_In_VD_of_obstacles;
    VD_of_obstacles.getVoronoiEdges(edges_In_VD_of_obstacles);
    for (list<VEdge2D*>::iterator i_vedge = edges_In_VD_of_obstacles.begin(); i_vedge != edges_In_VD_of_obstacles.end(); ++i_vedge) {
        VEdge2D* currVEdge = *i_vedge;
        currVEdge->compute_radius_interval_of_tangent_circles();
    }
}

void ShortestPathFinder2D::compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS()
{
     //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  m_BetaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = m_BetaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = m_BetaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {
                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                m_BetaUniverse,
                                                                                                m_SolutionPoolGraphByAllPathSearch,
                                                                                                m_PathByAllPathSearch,
                                                                                                m_PathTotalDistanceByAllPathSearch);
            }
            else
            {
                //1. find one topological path
                list<EdgeBU2D*> topologicalPathInEllipseFilter;
                find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, topologicalPathInEllipseFilter, m_QTEdgeValidation, m_BetaUniverse);

                switch (!topologicalPathInEllipseFilter.empty()) //TOPOLOGICAL PATH EXIST
                {
                case true:
                {
                    //2. find geodesic path
                    convert_filltered_QT_edges_to_topological_path(m_QTEdgeValidation, m_TopologicalPathByAllPathSearch);
                    find_one_geodesic_path_from_QT_topological_path_for_all_path_search(m_TopologicalPathByAllPathSearch, 
                                                                                        m_BetaUniverse, 
                                                                                        m_SolutionPoolGraphByAllPathSearch, 
                                                                                        m_PathByAllPathSearch, 
                                                                                        m_PathTotalDistanceByAllPathSearch);
                    break;
                }
                
                case false:
                {     
                    switch(solverForBarrierProblemType)
                    { 
                    case TPT_METHOD:
                    {
                        unordered_map<EdgeBU2D*, bool> isThisQTEdgeInsideOfExpandedEllipse;
                        compute_several_geodesic_paths_by_topological_path_tree_with_filtered_QT_edges(qtStartFace, qtEndFace, m_QTEdgeValidation, isThisQTEdgeInsideOfExpandedEllipse, m_BetaUniverse);
                        sort_index_of_distances_in_non_decreasing_order(m_PathTotalDistancesByTPT, m_IndexOfBFSsWithIncreasingOrderOfDistance);

                        break;
                    }

                    case EXPAND_PATH_ALONG_BARRIER_METHOD:
                    {
                        unordered_map<EdgeBU2D*, bool> expandedQTEdgeValidation;
                        expand_candidate_QT_edges_by_including_on_a_barrier(m_QTEdgeValidation, expandedQTEdgeValidation, m_BetaUniverse);

                        //find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByAllPathSearch, isThisQTEdgeInsideOfExpandedEllipse, betaUniverse);

                        convert_filltered_QT_edges_to_topological_path(expandedQTEdgeValidation, m_TopologicalPathByAllPathSearch);
                        find_one_geodesic_path_from_QT_topological_path_for_all_path_search(m_TopologicalPathByAllPathSearch,
                                                                                            m_BetaUniverse, 
                                                                                            m_SolutionPoolGraphByAllPathSearch, 
                                                                                            m_PathByAllPathSearch, 
                                                                                            m_PathTotalDistanceByAllPathSearch);

                        
                        break;
                    }

                    default:
                        break;
                    }
                    //compute_several_geodesic_paths_by_topological_path_tree(qtStartFace, qtEndFace, betaUniverse);
                    //sort_index_of_distances_in_non_decreasing_order(m_PathTotalDistancesByTPT, m_IndexOfBFSsWithIncreasingOrderOfDistance);
                    break;
                }

                default:
                    break;
                }
            }
          
            break;
        }

        case false:
        {
            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByAllPathSearch,
                                                                                                                           m_PathByAllPathSearch, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }

}

void ShortestPathFinder2D::compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS_ver1()
{
      //1. check intersection btw pq line and obstacles
    rg_Line2D lineFromStartToEnd(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());

    switch (is_this_line_segement_intersected_with_one_of_ALL_obstacles(lineFromStartToEnd))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = m_BetaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = m_BetaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {
                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                m_BetaUniverse,
                                                                                                m_SolutionPoolGraphByAllPathSearch,
                                                                                                m_PathByAllPathSearch,
                                                                                                m_PathTotalDistanceByAllPathSearch);
            }
            else
            {
                //1. find one topological path
                list<EdgeBU2D*> topologicalPathInEllipseFilter;
                find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, topologicalPathInEllipseFilter, m_QTEdgeValidation, m_BetaUniverse);

                switch (!topologicalPathInEllipseFilter.empty()) //TOPOLOGICAL PATH EXIST
                {
                case true:
                {
                    //2. find geodesic path
                    convert_filltered_QT_edges_to_topological_path(m_QTEdgeValidation, m_TopologicalPathByAllPathSearch);
                    find_one_geodesic_path_from_QT_topological_path_for_all_path_search(m_TopologicalPathByAllPathSearch, 
                                                                                        m_BetaUniverse, 
                                                                                        m_SolutionPoolGraphByAllPathSearch, 
                                                                                        m_PathByAllPathSearch, 
                                                                                        m_PathTotalDistanceByAllPathSearch);
                    break;
                }
                
                case false:
                {     
                    switch(solverForBarrierProblemType)
                    { 
                    case TPT_METHOD:
                    {
                        unordered_map<EdgeBU2D*, bool> isThisQTEdgeInsideOfExpandedEllipse;
                        compute_several_geodesic_paths_by_topological_path_tree_with_filtered_QT_edges(qtStartFace, qtEndFace, m_QTEdgeValidation, isThisQTEdgeInsideOfExpandedEllipse, m_BetaUniverse);
                        sort_index_of_distances_in_non_decreasing_order(m_PathTotalDistancesByTPT, m_IndexOfBFSsWithIncreasingOrderOfDistance);

                        break;
                    }

                    case EXPAND_PATH_ALONG_BARRIER_METHOD:
                    {
                        unordered_map<EdgeBU2D*, bool> expandedQTEdgeValidation;
                        expand_candidate_QT_edges_by_including_on_a_barrier(m_QTEdgeValidation, expandedQTEdgeValidation, m_BetaUniverse);

                        //find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByAllPathSearch, isThisQTEdgeInsideOfExpandedEllipse, betaUniverse);

                        convert_filltered_QT_edges_to_topological_path(expandedQTEdgeValidation, m_TopologicalPathByAllPathSearch);
                        find_one_geodesic_path_from_QT_topological_path_for_all_path_search(m_TopologicalPathByAllPathSearch,
                                                                                            m_BetaUniverse, 
                                                                                            m_SolutionPoolGraphByAllPathSearch, 
                                                                                            m_PathByAllPathSearch, 
                                                                                            m_PathTotalDistanceByAllPathSearch);

                        
                        break;
                    }

                    default:
                        break;
                    }
                    //compute_several_geodesic_paths_by_topological_path_tree(qtStartFace, qtEndFace, betaUniverse);
                    //sort_index_of_distances_in_non_decreasing_order(m_PathTotalDistancesByTPT, m_IndexOfBFSsWithIncreasingOrderOfDistance);
                    break;
                }

                default:
                    break;
                }
            }
          
            break;
        }

        case false:
        {
            rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(lineFromStartToEnd.getSP(), lineFromStartToEnd.getEP());
            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByAllPathSearch,
                                                                                                                           m_PathByAllPathSearch, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }
}

void ShortestPathFinder2D::compute_geodesic_path_by_DFS_with_filtered_QT_edges_ver1()
{
    //1. check intersection btw pq line and obstacles
    rg_Line2D lineFromStartToEnd(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());

    switch (is_this_line_segement_intersected_with_one_of_ALL_obstacles(lineFromStartToEnd))
    {

    case true:
    {
        // 1. Get start and end face of QT
        FaceBU2D* qtStartFace = m_BetaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
        FaceBU2D* qtEndFace   = m_BetaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

        if (qtStartFace == qtEndFace)
        {

            find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace,
                                                                                            m_BetaUniverse,
                                                                                            m_SolutionPoolGraphByDFS,
                                                                                            m_PathByDFS,
                                                                                            m_PathTotalDistanceByDFS);
        }
        else
        {
            //3. find one topological path
            find_one_topological_path_from_QT_by_DFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByDFS, m_QTEdgeValidation, m_BetaUniverse);

            if (m_TopologicalPathByDFS.size() == 0)
            {
                unordered_map<EdgeBU2D*, bool> expandedQTEdgeValidation;
                expand_candidate_QT_edges_by_including_on_a_barrier(m_QTEdgeValidation, expandedQTEdgeValidation, m_BetaUniverse);
                find_one_topological_path_from_QT_by_DFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByDFS, expandedQTEdgeValidation, m_BetaUniverse);

                //find_one_topological_path_from_QT_by_DFS(qtStartFace, qtEndFace, m_TopologicalPathByDFS, betaUniverse);
            }

            //4. find geodesic path
            find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByDFS,
                                                            m_BetaUniverse,
                                                            m_SolutionPoolGraphByDFS,
                                                            m_PathByDFS,
                                                            m_PathTotalDistanceByDFS);
        }

        break;
    }

    case false:
    {
        rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(lineFromStartToEnd.getSP(), lineFromStartToEnd.getEP());
        find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                       m_SolutionPoolGraphByDFS,
                                                                                                                       m_PathByDFS, m_PathTotalDistanceByDFS);
        break;
    }

    default:
        break;
    }
}

void ShortestPathFinder2D::compute_geodesic_path_by_BFS_with_filtered_QT_edges_ver1()
{
    //1. check intersection btw pq line and obstacles
    rg_Line2D lineFromStartToEnd(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());

    switch (is_this_line_segement_intersected_with_one_of_ALL_obstacles(lineFromStartToEnd))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = m_BetaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = m_BetaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {

                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                m_BetaUniverse,
                                                                                                m_SolutionPoolGraphByBFS,
                                                                                                m_PathByBFS,
                                                                                                m_PathTotalDistanceByBFS);
            }
            else
            {          
                //3. find one topological path
                find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByBFS, m_QTEdgeValidation, m_BetaUniverse);

                if (m_TopologicalPathByBFS.size() == 0)
                {
                    unordered_map<EdgeBU2D*, bool> expandedQTEdgeValidation;
                    expand_candidate_QT_edges_by_including_on_a_barrier(m_QTEdgeValidation, expandedQTEdgeValidation, m_BetaUniverse);
                    find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByBFS, expandedQTEdgeValidation, m_BetaUniverse);

                    //find_one_topological_path_from_QT_by_BFS(qtStartFace, qtEndFace, m_TopologicalPathByBFS, betaUniverse);
                }

                //4. find geodesic path
                find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByBFS, 
                                                                m_BetaUniverse, 
                                                                m_SolutionPoolGraphByBFS, 
                                                                m_PathByBFS, 
                                                                m_PathTotalDistanceByBFS);
            }
          
            break;
        }

        case false:
        {
            rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(lineFromStartToEnd.getSP(), lineFromStartToEnd.getEP());
            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByBFS,
                                                                                                                           m_PathByBFS, m_PathTotalDistanceByBFS);
            break;
        }

    default:
        break;
    }
}



void ShortestPathFinder2D::compute_geodesic_path_by_DFS_with_filtered_QT_edges()
{
    //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());

    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination,
                                                                                  m_StartPt.getCenterPt(),
                                                                                  m_EndPt.getCenterPt(),
                                                                                  m_BetaUniverse))
    {
    case true:
    {
        // 1. Get start and end face of QT
        FaceBU2D* qtStartFace = m_BetaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
        FaceBU2D* qtEndFace   = m_BetaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

        if (qtStartFace == qtEndFace)
        {

            find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace,
                                                                                            m_BetaUniverse,
                                                                                            m_SolutionPoolGraphByDFS,
                                                                                            m_PathByDFS,
                                                                                            m_PathTotalDistanceByDFS);
        }
        else
        {
            //3. find one topological path
            find_one_topological_path_from_QT_by_DFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByDFS, m_QTEdgeValidation, m_BetaUniverse);

            if (m_TopologicalPathByDFS.size() == 0)
            {
                unordered_map<EdgeBU2D*, bool> expandedQTEdgeValidation;
                expand_candidate_QT_edges_by_including_on_a_barrier(m_QTEdgeValidation, expandedQTEdgeValidation, m_BetaUniverse);
                find_one_topological_path_from_QT_by_DFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByDFS, expandedQTEdgeValidation, m_BetaUniverse);

                //find_one_topological_path_from_QT_by_DFS(qtStartFace, qtEndFace, m_TopologicalPathByDFS, betaUniverse);
            }

            //4. find geodesic path
            find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByDFS,
                                                            m_BetaUniverse,
                                                            m_SolutionPoolGraphByDFS,
                                                            m_PathByDFS,
                                                            m_PathTotalDistanceByDFS);
        }

        break;
    }

    case false:
    {

        find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                       m_SolutionPoolGraphByDFS,
                                                                                                                       m_PathByDFS, m_PathTotalDistanceByDFS);
        break;
    }

    default:
        break;
    }
}



void ShortestPathFinder2D::compute_geodesic_path_by_BFS(BetaUniverse2D& betaUniverse, list<EdgeForSolutionPoolGraph*>& geodesicPath)
{
    //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  betaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = betaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = betaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {

                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                betaUniverse,
                                                                                                m_SolutionPoolGraphByBFS,
                                                                                                geodesicPath,
                                                                                                m_PathTotalDistanceByBFS);
            }
            else
            {          
                //3. find one topological path
                find_one_topological_path_from_QT_by_BFS(qtStartFace, qtEndFace, m_TopologicalPathByBFS, betaUniverse);

                //4. find geodesic path
                find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByBFS, 
                                                                betaUniverse, 
                                                                m_SolutionPoolGraphByBFS, 
                                                                geodesicPath, 
                                                                m_PathTotalDistanceByBFS);
            }
          
            break;
        }

        case false:
        {

            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByBFS,
                                                                                                                           geodesicPath, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }

}

void ShortestPathFinder2D::compute_several_geodesic_paths_by_topological_path_tree(FaceBU2D* qtStartFace, FaceBU2D* qtEndFace, BetaUniverse2D & betaUniverse)
{
    m_TopologicalPathsByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);
    m_SolutionPoolGraphsByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);
    m_PathsByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);
    m_PathTotalDistancesByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);


    compute_topological_path_tree(qtStartFace, qtEndFace, betaUniverse, m_TopologicalPathTree);

    find_several_topological_path_from_QT_by_topological_path_tree(qtStartFace, qtEndFace, m_TopologicalPathsByTPT, m_TopologicalPathTree, betaUniverse);

    find_several_geodesic_paths_from_QT_topological_paths(m_TopologicalPathsByTPT, 
                                                          betaUniverse, 
                                                          m_SolutionPoolGraphsByTPT, 
                                                          m_PathsByTPT, 
                                                          m_PathTotalDistancesByTPT);
}

void ShortestPathFinder2D::compute_several_geodesic_paths_by_topological_path_tree_with_filtered_QT_edges(FaceBU2D * qtStartFace, FaceBU2D * qtEndFace, const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfExpandedEllipse, BetaUniverse2D & betaUniverse)
{
    m_TopologicalPathsByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);
    m_SolutionPoolGraphsByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);
    m_PathsByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);
    m_PathTotalDistancesByTPT.resize(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE);

    expand_candidate_QT_edges_by_including_on_a_barrier(isThisQTEdgeInsideOfEllipse, isThisQTEdgeInsideOfExpandedEllipse, betaUniverse);

    compute_topological_path_tree_with_filtered_QT_edges(qtStartFace, qtEndFace, betaUniverse, m_TopologicalPathTree, isThisQTEdgeInsideOfExpandedEllipse);

    find_several_topological_path_from_QT_by_topological_path_tree(qtStartFace, qtEndFace, m_TopologicalPathsByTPT, m_TopologicalPathTree, betaUniverse);

    find_several_geodesic_paths_from_QT_topological_paths(m_TopologicalPathsByTPT, 
                                                          betaUniverse, 
                                                          m_SolutionPoolGraphsByTPT, 
                                                          m_PathsByTPT, 
                                                          m_PathTotalDistancesByTPT);
}



void ShortestPathFinder2D::compute_geodesic_path_by_BFS_with_filtered_QT_edges()
{
      //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  m_BetaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = m_BetaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = m_BetaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {

                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                m_BetaUniverse,
                                                                                                m_SolutionPoolGraphByBFS,
                                                                                                m_PathByBFS,
                                                                                                m_PathTotalDistanceByBFS);
            }
            else
            {          
                //3. find one topological path
                find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByBFS, m_QTEdgeValidation, m_BetaUniverse);

                if (m_TopologicalPathByBFS.size() == 0)
                {
                    unordered_map<EdgeBU2D*, bool> expandedQTEdgeValidation;
                    expand_candidate_QT_edges_by_including_on_a_barrier(m_QTEdgeValidation, expandedQTEdgeValidation, m_BetaUniverse);
                    find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByBFS, expandedQTEdgeValidation, m_BetaUniverse);

                    //find_one_topological_path_from_QT_by_BFS(qtStartFace, qtEndFace, m_TopologicalPathByBFS, betaUniverse);
                }

                //4. find geodesic path
                find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByBFS, 
                                                                m_BetaUniverse, 
                                                                m_SolutionPoolGraphByBFS, 
                                                                m_PathByBFS, 
                                                                m_PathTotalDistanceByBFS);
            }
          
            break;
        }

        case false:
        {

            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByBFS,
                                                                                                                           m_PathByBFS, m_PathTotalDistanceByBFS);
            break;
        }

    default:
        break;
    }
}


void ShortestPathFinder2D::compute_geodesic_path_by_DFS(BetaUniverse2D& betaUniverse, list<EdgeForSolutionPoolGraph*>& geodesicPath)
{
       //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  betaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = betaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = betaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {

                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                betaUniverse,
                                                                                                m_SolutionPoolGraphByDFS,
                                                                                                geodesicPath,
                                                                                                m_PathTotalDistanceByDFS);
            }
            else
            {          
                //3. find one topological path
                find_one_topological_path_from_QT_by_DFS(qtStartFace, qtEndFace, m_TopologicalPathByDFS, betaUniverse);

                //4. find geodesic path
                find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByDFS, 
                                                                betaUniverse, 
                                                                m_SolutionPoolGraphByDFS, 
                                                                geodesicPath, 
                                                                m_PathTotalDistanceByDFS);
            }
          
            break;
        }

        case false:
        {

            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByDFS,
                                                                                                                           geodesicPath, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }
}


void ShortestPathFinder2D::compute_geodesic_path_by_all_path_search(BetaUniverse2D& betaUniverse, list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath)
{
        //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  betaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = betaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = betaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {
                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                betaUniverse,
                                                                                                m_SolutionPoolGraphByAllPathSearch,
                                                                                                optimalGeodesicPath,
                                                                                                m_PathTotalDistanceByAllPathSearch);
            }
            else
            {
                compute_topological_path_tree(qtStartFace, qtEndFace, betaUniverse, m_TopologicalPathTree);

                find_optimal_path_from_many_topological_ones(m_TopologicalPathByAllPathSearch,
                                                             m_SolutionPoolGraphByAllPathSearch,
                                                             m_PathByAllPathSearch,
                                                             m_PathTotalDistanceByAllPathSearch,
                                                             m_TopologicalPathTree,
                                                             betaUniverse);

                /* RECURSIVE FUNCTION

                list<EdgeBU2D*> currTopologicalPath;

                // 1. Initialize for breadth-first search
                initialize_visitation_tag_of_QT_face(betaUniverse);
                qtStartFace->isVisited(true);

                // 2. Start from maximum three edge
                scan_neighbor_QT_face_and_step_into_one_QT_face(qtStartFace, 
                                                                qtEndFace,
                                                                betaUniverse, 
                                                                currTopologicalPath, 
                                                                m_TopologicalPathByAllPathSearch, 
                                                                m_SolutionPoolGraphByAllPathSearch,
                                                                optimalGeodesicPath, 
                                                                m_PathTotalDistanceByAllPathSearch);

                */
            }
          
            break;
        }

        case false:
        {

            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByAllPathSearch,
                                                                                                                           optimalGeodesicPath, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }
}



void ShortestPathFinder2D::compute_geodesic_path_by_all_path_search_with_filtered_QT_edges(BetaUniverse2D& betaUniverse, const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath)
{
    //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  betaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = betaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = betaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {
                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                betaUniverse,
                                                                                                m_SolutionPoolGraphByAllPathSearch,
                                                                                                optimalGeodesicPath,
                                                                                                m_PathTotalDistanceByAllPathSearch);
            }
            else
            {
                //1. find one topological path
                find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByBFS, isThisQTEdgeInsideOfEllipse, betaUniverse);

                switch (is_topological_path_by_BFS_computed())
                {
                case true:
                {
                    //2. find geodesic path
                    convert_filltered_QT_edges_to_topological_path(isThisQTEdgeInsideOfEllipse, m_TopologicalPathByAllPathSearch);
                    find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByAllPathSearch, 
                                                                    betaUniverse, 
                                                                    m_SolutionPoolGraphByAllPathSearch, 
                                                                    optimalGeodesicPath, 
                                                                    m_PathTotalDistanceByAllPathSearch);

                    //for drawing
                   find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByBFS,
                                                                   betaUniverse, 
                                                                   m_SolutionPoolGraphByBFS, 
                                                                   m_PathByBFS, 
                                                                   m_PathTotalDistanceByBFS);

                    break;
                }
                
                case false:
                {              
                    //2. find geodesic path
                    find_one_topological_path_from_QT_by_BFS(qtStartFace, qtEndFace, m_TopologicalPathByBFS, betaUniverse);

                    find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByBFS,
                                                                    betaUniverse, 
                                                                    m_SolutionPoolGraphByBFS, 
                                                                    m_PathByBFS, 
                                                                    m_PathTotalDistanceByBFS);
                    
                    //for drawing
                     m_TopologicalPathByAllPathSearch   = m_TopologicalPathByBFS;
                     m_SolutionPoolGraphByAllPathSearch = m_SolutionPoolGraphByBFS;
                     m_PathTotalDistanceByAllPathSearch = m_PathTotalDistanceByBFS;
                     optimalGeodesicPath                = m_PathByBFS;

                    break;
                }

                default:
                    break;
                }
 

                /*
                convert_filltered_QT_edges_to_topological_path(isThisQTEdgeInsideOfEllipse, m_TopologicalPathByAllPathSearch);
                find_one_geodesic_path_from_QT_topological_path(m_TopologicalPathByAllPathSearch, 
                                                                betaUniverse, 
                                                                m_SolutionPoolGraphByAllPathSearch, 
                                                                optimalGeodesicPath, 
                                                                m_PathTotalDistanceByAllPathSearch);
                */

                /*
                compute_topological_path_tree_with_filtered_QT_edges(qtStartFace, qtEndFace, betaUniverse, m_TopologicalPathTree, isThisQTEdgeInsideOfEllipse);

                find_optimal_path_from_many_topological_ones(m_TopologicalPathByAllPathSearch,
                                                             m_SolutionPoolGraphByAllPathSearch,
                                                             m_PathByAllPathSearch,
                                                             m_PathTotalDistanceByAllPathSearch,
                                                             m_TopologicalPathTree,
                                                             betaUniverse);
                */
            }
          
            break;
        }

        case false:
        {

            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByAllPathSearch,
                                                                                                                           optimalGeodesicPath, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }
}




void ShortestPathFinder2D::compute_topological_path_tree(FaceBU2D* qtStartFace, 
                                                         FaceBU2D* qtEndFace, 
                                                         const BetaUniverse2D& betaUniverse, 
                                                         TopologicalPathTree& topologicalPathTree) const
{
    if ((!are_some_bounding_edges_of_this_QT_face_exterior(qtStartFace))
        || (!are_some_bounding_edges_of_this_QT_face_exterior(qtEndFace)))
    {
        return;
    }

    list<pair<FaceBU2D*, TPTNode*>> queueOfFace;
    queueOfFace.push_back(make_pair(qtStartFace, (TPTNode*)NULL));

    while (!queueOfFace.empty())
    {
        pair<FaceBU2D*, TPTNode*> currFaceNParentNode = queueOfFace.front();
        queueOfFace.pop_front();

        FaceBU2D* currFace  = currFaceNParentNode.first;
        TPTNode* parentNode = currFaceNParentNode.second;

        TPTNode* currNode = topologicalPathTree.create_node(TPTNode(currFace, parentNode));

        if (parentNode != NULL)
        {
            parentNode->add_a_child(currNode);
        }

        if (currFace == qtEndFace)
        {
            topologicalPathTree.add_leaf_node(currNode);
            continue;
        }

        //find next face and insert into queue
        rg_dList<EdgeBU2D*> boundingEdges;
        currFace->getBoundingEdges(boundingEdges);

        boundingEdges.reset4Loop();
        while (boundingEdges.setNext4Loop())
        {
            EdgeBU2D* boundingEdge = boundingEdges.getEntity();
            FaceBU2D* adjacentFace = get_opposite_QT_face_of_Edge(boundingEdge, currFace);

            if (is_this_QT_edge_in_an_exterior_state(boundingEdge)
            && (!topologicalPathTree.does_this_node_have_ancestor_having_this_face(currNode, adjacentFace)))
            {
                // insert into queue
                queueOfFace.push_back(make_pair(adjacentFace, currNode));
            }
        }
    }
}


void ShortestPathFinder2D::compute_topological_path_tree_with_filtered_QT_edges(FaceBU2D* qtStartFace, FaceBU2D* qtEndFace,
                                                                                const BetaUniverse2D& betaUniverse, 
                                                                                TopologicalPathTree& topologicalPathTree, 
                                                                                const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const
{
     if ((!are_some_bounding_edges_of_this_QT_face_exterior(qtStartFace))
        || (!are_some_bounding_edges_of_this_QT_face_exterior(qtEndFace)))
    {
        return;
    }

    list<pair<FaceBU2D*, TPTNode*>> queueOfFace;
    queueOfFace.push_back(make_pair(qtStartFace, (TPTNode*)NULL));

    while (!queueOfFace.empty())
    {
        pair<FaceBU2D*, TPTNode*> currFaceNParentNode = queueOfFace.front();
        queueOfFace.pop_front();

        FaceBU2D* currFace  = currFaceNParentNode.first;
        TPTNode* parentNode = currFaceNParentNode.second;

        TPTNode* currNode = topologicalPathTree.create_node(TPTNode(currFace, parentNode));

        if (parentNode != NULL)
        {
            parentNode->add_a_child(currNode);
        }

        if (currFace == qtEndFace)
        {
            topologicalPathTree.add_leaf_node(currNode);
            continue;
        }

        //find next face and insert into queue
        rg_dList<EdgeBU2D*> boundingEdges;
        currFace->getBoundingEdges(boundingEdges);

        boundingEdges.reset4Loop();
        while (boundingEdges.setNext4Loop())
        {
            EdgeBU2D* boundingEdge = boundingEdges.getEntity();
            FaceBU2D* adjacentFace = get_opposite_QT_face_of_Edge(boundingEdge, currFace);

            if (isThisQTEdgeInsideOfEllipse.at(boundingEdge)
            && is_this_QT_edge_in_an_exterior_state(boundingEdge)
            && (!topologicalPathTree.does_this_node_have_ancestor_having_this_face(currNode, adjacentFace)))
            {
                // insert into queue
                queueOfFace.push_back(make_pair(adjacentFace, currNode));
            }
        }
    }
}




void ShortestPathFinder2D::find_optimal_path_from_many_topological_ones(list<EdgeBU2D*>& optimalTopologicalPath, 
                                                                        ShortestPathSolutionPoolGraph& optimalSolutionPoolGraph, 
                                                                        list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, 
                                                                        double& optimalDistance, 
                                                                        TopologicalPathTree& topologicalPathTree,
                                                                        BetaUniverse2D& betaUniverse)
{
    list<TPTNode*> leafNodes;
    topologicalPathTree.get_leaf_nodes(leafNodes);

    for (list<TPTNode*>::const_iterator it_LeafNode = leafNodes.begin();
        it_LeafNode != leafNodes.end();
        it_LeafNode++)
    {
        TPTNode* currLeafNode = *it_LeafNode;

        //1. get topological qt-face path
        list<FaceBU2D*> currTopologicalQTFacePath;
        topologicalPathTree.get_topological_path_from_this_leaf_node(currLeafNode, currTopologicalQTFacePath);


        //2. convert qt-face path to qt-edge path
        list<EdgeBU2D*> currTopologicalQTEdgePath;
        convert_QT_face_path_to_QT_edge_path(currTopologicalQTFacePath, currTopologicalQTEdgePath);


        //3. compute geodesic from qt-edge path
        ShortestPathSolutionPoolGraph currSolutionGraph;
        list<EdgeForSolutionPoolGraph*> currGeodesicPath;
        double currGeodesicDistance;

        find_one_geodesic_path_from_QT_topological_path(currTopologicalQTEdgePath, 
                                                        betaUniverse, 
                                                        currSolutionGraph,
                                                        currGeodesicPath,
                                                        currGeodesicDistance);
        if (currGeodesicDistance < optimalDistance)
        {
            optimalTopologicalPath.clear();
            optimalGeodesicPath.clear();

            optimalDistance          = currGeodesicDistance;
            optimalTopologicalPath   = currTopologicalQTEdgePath;
            optimalSolutionPoolGraph = currSolutionGraph;
            optimalGeodesicPath      = currGeodesicPath;

            change_optimal_path_in_old_graph_to_new_one(currSolutionGraph, optimalSolutionPoolGraph, optimalGeodesicPath);
        }
    }
}


void ShortestPathFinder2D::convert_QT_face_path_to_QT_edge_path(const list<FaceBU2D*>& qtFacePath, list<EdgeBU2D*>& qtEdgePath) const
{
    for (list<FaceBU2D*>::const_iterator it_CurrFace = qtFacePath.begin();
         it_CurrFace != qtFacePath.end();
         it_CurrFace++)
    {
        list<FaceBU2D*>::const_iterator it_NextFace = it_CurrFace;
        it_NextFace++;

        if (it_NextFace == qtFacePath.end())
        {
            break;
        }
        else
        {
            FaceBU2D* currFace = *it_CurrFace;
            FaceBU2D* nextFace = *it_NextFace;

            EdgeBU2D* sharingQTEdge = get_sharing_QT_edge_of_two_QT_faces(currFace, nextFace);
            qtEdgePath.push_back(sharingQTEdge);
        }
    }
}


void ShortestPathFinder2D::convert_QT_edge_path_to_QT_faces(const list<EdgeBU2D*>& qtEdgePath, unordered_set<FaceBU2D*>& qtFacePath) const
{
    for (list<EdgeBU2D*>::const_iterator it_CurrEdge = qtEdgePath.begin();
        it_CurrEdge != qtEdgePath.end();
        it_CurrEdge++)
    {
        FaceBU2D* leftFace  = (*it_CurrEdge)->getLeftFace();
        FaceBU2D* rightFace = (*it_CurrEdge)->getRightFace();

        qtFacePath.insert(leftFace);
        qtFacePath.insert(rightFace);
    }
}


void ShortestPathFinder2D::find_a_geodesic_path_by_stepping_into_the_graph(EdgeBU2D* qtStartEdge, 
                                                                           FaceBU2D* qtStartFace, 
                                                                           FaceBU2D* qtEndFace, 
                                                                           BetaUniverse2D& betaUniverse, 
                                                                           list<EdgeBU2D*>& currTopologicalPath, 
                                                                           list<EdgeBU2D*>& optimalTopologicalPath, 
                                                                           ShortestPathSolutionPoolGraph& optimalSolutionPoolGraph, 
                                                                           list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, 
                                                                           double& optimalDistance)
{
    qtStartFace->isVisited(true);
    currTopologicalPath.push_back(qtStartEdge);

    if (qtStartFace == qtEndFace)
    {
        ShortestPathSolutionPoolGraph   solutionPoolGraph;
        list<EdgeForSolutionPoolGraph*> currGeodesicPath;

        double totalDistance = DBL_MAX;

        find_one_geodesic_path_from_QT_topological_path(currTopologicalPath,
                                                        betaUniverse, 
                                                        solutionPoolGraph,
                                                        currGeodesicPath,
                                                        totalDistance);

        if (totalDistance < optimalDistance)
        {
            optimalTopologicalPath.clear();
            optimalGeodesicPath.clear();

            optimalDistance          = totalDistance;
            optimalTopologicalPath   = currTopologicalPath;
            optimalSolutionPoolGraph = solutionPoolGraph;
            optimalGeodesicPath      = currGeodesicPath;

            change_optimal_path_in_old_graph_to_new_one(solutionPoolGraph, optimalSolutionPoolGraph, optimalGeodesicPath);
        }
    }                                                   
    else
    {
        scan_neighbor_QT_face_and_step_into_one_QT_face(qtStartFace, 
                                                        qtEndFace,
                                                        betaUniverse,
                                                        currTopologicalPath, 
                                                        optimalTopologicalPath,
                                                        optimalSolutionPoolGraph,
                                                        optimalGeodesicPath, 
                                                        optimalDistance);
    }


    currTopologicalPath.pop_back();
    qtStartFace->isVisited(false);
}



void ShortestPathFinder2D::scan_neighbor_QT_face_and_step_into_one_QT_face(FaceBU2D* qtStartFace, 
                                                                           FaceBU2D* qtEndFace, 
                                                                           BetaUniverse2D& betaUniverse, 
                                                                           list<EdgeBU2D*>& currTopologicalPath, 
                                                                           list<EdgeBU2D*>& optimalTopologicalPath, 
                                                                           ShortestPathSolutionPoolGraph& optimalSolutionPoolGraph, 
                                                                           list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, 
                                                                           double& optimalDistance) 
{
     list<EdgeBU2D*> boundingQTEdges;
     insert_three_bounding_QT_edge_into_container(qtStartFace, boundingQTEdges);
     
     for (list<EdgeBU2D*>::const_iterator it_QTEdge = boundingQTEdges.begin();
         it_QTEdge != boundingQTEdges.end();
         it_QTEdge++)
     {
         EdgeBU2D* currQTEdge = *it_QTEdge;
     
         if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
         {
             FaceBU2D* nextQTFace = get_opposite_QT_face_of_Edge(currQTEdge, qtStartFace);
     
             if (!nextQTFace->isVisited())
             {
                 find_a_geodesic_path_by_stepping_into_the_graph(currQTEdge,
                                                                 nextQTFace, 
                                                                 qtEndFace,
                                                                 betaUniverse,
                                                                 currTopologicalPath,
                                                                 optimalTopologicalPath,
                                                                 optimalSolutionPoolGraph,
                                                                 optimalGeodesicPath,
                                                                 optimalDistance);
             }
             
         }
     }
}


void ShortestPathFinder2D::find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(FaceBU2D* qtFace, 
                                                                                                           BetaUniverse2D& betaUniverse, 
                                                                                                           ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                                                           list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                                                           double& totalPathDistacne) 
{
     //1. Make obstacle pool 
     vector<VertexBU2D*> obstaclePool;
     make_obstacle_pool(obstaclePool, qtFace);

     //2. generate solution pool
     generate_solution_pool_graph(obstaclePool, betaUniverse, solutionPoolGraph);

     //3. dijkstra algorithm in solution pool graph
     find_geodesic_path_by_dijkstra(solutionPoolGraph);

     //4. back tracking
     back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(),
                                                       solutionPoolGraph.get_end_vtx(),
                                                       geodesicPath,
                                                       totalPathDistacne);
}


void ShortestPathFinder2D::find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(const rg_ImplicitEquation& lineEqFromSourceToDestination, 
                                                                                                                                         ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                                                                                         list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath,
                                                                                                                                         double& totalPathDistance) const
{
     VertexForSolutionPoolGraph* sourceVtx      = NULL;
     VertexForSolutionPoolGraph* destinationVtx = NULL;
     EdgeForSolutionPoolGraph*   lineEdge       = NULL;
     
     create_two_vertices_and_a_tangent_line_edge(sourceVtx, destinationVtx, lineEdge, solutionPoolGraph);
     set_two_vertices_and_a_tangent_line_edge(sourceVtx, 
                                              destinationVtx, 
                                              lineEdge, 
                                              make_pair(m_StartPt.getCenterPt(), m_EndPt.getCenterPt()), 
                                              lineEqFromSourceToDestination);

     totalPathDistance = lineEdge->get_edge_length();
     optimalGeodesicPath.push_back(lineEdge);
}



rg_Circle2D ShortestPathFinder2D::find_next_stop_point_when_both_of_QT_faces_from_source_N_destination_are_same(FaceBU2D * qtFace, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, list<EdgeForSolutionPoolGraph*>& geodesicPath, double & totalPathDistacne)
{


    return rg_Circle2D();
}



rg_Circle2D ShortestPathFinder2D::find_next_stop_point_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(const rg_ImplicitEquation & lineEqFromSourceToDestination, ShortestPathSolutionPoolGraph & solutionPoolGraph, list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, double & totalPathDistance) const
{


    return rg_Circle2D();
}

rg_Circle2D ShortestPathFinder2D::find_next_stop_point_in_general_case(BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, double & totalPathDistance) 
{
    rg_Circle2D nextStopPoint;

    //1. compute ellipse filter
    unordered_map<VertexBU2D*, bool> isThisQTVertexInsideOfExpandedEllipse;

    compute_ellipse_filter(m_StartPt, m_EndPt, m_EllipseFilter);
    mark_QT_vertex_whether_in_or_out_of_ellipse(m_BetaUniverse, m_EllipseFilter, m_QTVertexValidation);
    expand_candidate_QT_vertices_by_including_on_a_barrier(m_QTVertexValidation, isThisQTVertexInsideOfExpandedEllipse, betaUniverse);

    //2. Make obstacle pool 
    vector<VertexBU2D*> obstaclePool;
    make_obstacle_pool(obstaclePool, isThisQTVertexInsideOfExpandedEllipse);

    //3. generate solution pool
    generate_visible_solution_pool_graph_from_current_position(m_StartPt.getCenterPt(), obstaclePool, betaUniverse, solutionPoolGraph);

    //4. dijkstra algorithm in solution pool graph
    


    return nextStopPoint;
}



void ShortestPathFinder2D::find_one_topological_path_from_QT_by_BFS(FaceBU2D* qtStartFace, 
                                                                    FaceBU2D* qtEndFace, 
                                                                    list<EdgeBU2D*>& topologicalPath, 
                                                                    BetaUniverse2D& betaUniverse) const
{
    // 1. Initialize for breadth-first search
    initialize_visitation_tag_of_QT_face(betaUniverse);

    // 2. Breadth-first search
    EdgeBU2D* gatewayQtEdgeOfEndFace = NULL;
    unordered_map<EdgeBU2D*, EdgeBU2D*> currToPrevQTEdgeMapper;

    compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_BFS(qtStartFace, qtEndFace, gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper);

    // 3. Back-tracing
    back_trace_to_QT_start_face_and_find_topological_path(gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, topologicalPath);
}



void ShortestPathFinder2D::find_several_topological_path_from_QT_by_topological_path_tree(FaceBU2D * qtStartFace, FaceBU2D * qtEndFace, vector<list<EdgeBU2D*>>& topologicalPaths, TopologicalPathTree& topologicalPathTree, BetaUniverse2D & betaUniverse) const
{
    vector<list<FaceBU2D*>> topologicalQTFacePaths;

    topologicalPathTree.get_topological_paths_in_non_decreasing_order_by(INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE, topologicalQTFacePaths);

    for (int i = 0; i < topologicalQTFacePaths.size(); i++)
    {
        list<EdgeBU2D*> currQTEdgeTopologicalPath;

        list<FaceBU2D*>& currQTFaceTopologicalPath = topologicalQTFacePaths.at(i);
        convert_QT_face_path_to_QT_edge_path(currQTFaceTopologicalPath, currQTEdgeTopologicalPath);
        topologicalPaths.at(i) = currQTEdgeTopologicalPath;
    }


    /*
    // 1. Initialize for breadth-first search
    initialize_visitation_tag_of_QT_face(betaUniverse);

    // 2. Breadth-first search
    vector<EdgeBU2D*> gatewayQtEdgesOfEndFace;
    unordered_map<EdgeBU2D*, EdgeBU2D*> currToPrevQTEdgeMapper;

    compute_gateway_QT_edges_of_end_face_and_its_prev_edges_by_BFS(qtStartFace, qtEndFace, gatewayQtEdgesOfEndFace, currToPrevQTEdgeMapper);

    // 3. Back-tracing
    back_trace_to_QT_start_face_and_find_topological_paths(gatewayQtEdgesOfEndFace, currToPrevQTEdgeMapper, topologicalPaths);
    */
}



void ShortestPathFinder2D::find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(FaceBU2D* qtStartFace, 
                                                                                           FaceBU2D* qtEndFace, 
                                                                                           list<EdgeBU2D*>& topologicalPath, 
                                                                                           const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, 
                                                                                           BetaUniverse2D& betaUniverse) const
{
    // 1. Initialize for breadth-first search
    initialize_visitation_tag_of_QT_face(betaUniverse);

    // 2. Breadth-first search
    EdgeBU2D* gatewayQtEdgeOfEndFace = NULL;
    unordered_map<EdgeBU2D*, EdgeBU2D*> currToPrevQTEdgeMapper;

    compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, isThisQTEdgeInsideOfEllipse);

    // 3. Back-tracing
    back_trace_to_QT_start_face_and_find_topological_path(gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, topologicalPath);

}

void ShortestPathFinder2D::find_one_topological_path_from_QT_by_DFS_with_filtered_QT_edges(FaceBU2D* qtStartFace, 
                                                                                           FaceBU2D* qtEndFace, 
                                                                                           list<EdgeBU2D*>& topologicalPath, 
                                                                                           const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, 
                                                                                           BetaUniverse2D& betaUniverse) const
{
    // 1. Initialize for breadth-first search
    initialize_visitation_tag_of_QT_face(betaUniverse);

    // 2. Breadth-first search
    EdgeBU2D* gatewayQtEdgeOfEndFace = NULL;
    unordered_map<EdgeBU2D*, EdgeBU2D*> currToPrevQTEdgeMapper;

    compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_DFS_with_filtered_QT_edges(qtStartFace, qtEndFace, gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, isThisQTEdgeInsideOfEllipse);
    //compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_adjusted_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, isThisQTEdgeInsideOfEllipse);
    // 3. Back-tracing
    back_trace_to_QT_start_face_and_find_topological_path(gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, topologicalPath);

}



void ShortestPathFinder2D::find_one_topological_path_from_QT_by_DFS(FaceBU2D* qtStartFace, FaceBU2D* qtEndFace, list<EdgeBU2D*>& topologicalPath, BetaUniverse2D& betaUniverse) const
{
    // 1. Initialize for breadth-first search
    initialize_visitation_tag_of_QT_face(betaUniverse);

    // 2. Breadth-first search
    EdgeBU2D* gatewayQtEdgeOfEndFace = NULL;
    unordered_map<EdgeBU2D*, EdgeBU2D*> currToPrevQTEdgeMapper;

    compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_DFS(qtStartFace, qtEndFace, gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper);

    // 3. Back-tracing
    back_trace_to_QT_start_face_and_find_topological_path(gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, topologicalPath);
}


void ShortestPathFinder2D::initialize_visitation_tag_of_QT_face(BetaUniverse2D& betaUniverse) const
{
    rg_dList<FaceBU2D>& allQTFaces  = betaUniverse.getFaces();

    allQTFaces.reset4Loop();
    while (allQTFaces.setNext4Loop())
    {
        FaceBU2D* currQTFace = allQTFaces.getpEntity();
        currQTFace->isVisited(false);
    }
}


void ShortestPathFinder2D::insert_three_bounding_QT_edge_into_container(FaceBU2D* qtFace, list <EdgeBU2D*>& stackOfQTEdges) const
{
    rg_dList<EdgeBU2D*> boundingEdgesForQTStartFace;
    qtFace->getBoundingEdges(boundingEdgesForQTStartFace);

    boundingEdgesForQTStartFace.reset4Loop();
    while (boundingEdgesForQTStartFace.setNext4Loop())
    {
        EdgeBU2D* currQTEdge = boundingEdgesForQTStartFace.getEntity();

        stackOfQTEdges.push_back(currQTEdge);
    }
}


void ShortestPathFinder2D::insert_three_bounding_QT_edge_into_container(FaceBU2D * qtFace, EntityAccessiblePriorityQ<EdgeBU2D*>& priorityQ) const
{
    rg_dList<EdgeBU2D*> boundingEdgesForQTStartFace;
    qtFace->getBoundingEdges(boundingEdgesForQTStartFace);

    boundingEdgesForQTStartFace.reset4Loop();
    while (boundingEdgesForQTStartFace.setNext4Loop())
    {
        EdgeBU2D* currQTEdge = boundingEdgesForQTStartFace.getEntity();

        insert_QT_edge_into_PQ(currQTEdge, priorityQ);
    }
}


void ShortestPathFinder2D::insert_QT_edge_into_PQ(EdgeBU2D * qtEdge, EntityAccessiblePriorityQ<EdgeBU2D*>& priorityQ) const
{
    if (qtEdge->isVirtual())
    {
        priorityQ.push(make_pair(qtEdge, 2));
    }
    else
    {
        priorityQ.push(make_pair(qtEdge, 1));
    }
}


void ShortestPathFinder2D::compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_BFS(FaceBU2D* qtStartFace, 
                                                                                         FaceBU2D* qtEndFace, 
                                                                                         EdgeBU2D*& gatewayQtEdgeOfEndFace,
                                                                                         unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const
{
    list<EdgeBU2D*> queueOfQTEdges;

    insert_three_bounding_QT_edge_into_container(qtStartFace, queueOfQTEdges);

    qtStartFace->isVisited(true);

    while (!queueOfQTEdges.empty())
    {
        EdgeBU2D* currQTEdge = queueOfQTEdges.front();
        queueOfQTEdges.pop_front();

        if (qtEndFace->isIncidentTo(currQTEdge) && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            gatewayQtEdgeOfEndFace = currQTEdge;
            break;
        }

        if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            FaceBU2D* currRightQTFace = currQTEdge->getRightFace();

            if ((!is_this_QT_face_visited(currRightQTFace)))
            {
                currRightQTFace->isVisited(true);

                queueOfQTEdges.push_back(currQTEdge->getRightHand());
                queueOfQTEdges.push_back(currQTEdge->getRightLeg());

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightLeg(), currQTEdge));
            }

            FaceBU2D* currLeftQTFace = currQTEdge->getLeftFace();

            if ((!is_this_QT_face_visited(currLeftQTFace)))
            {
                currLeftQTFace->isVisited(true);

                queueOfQTEdges.push_back(currQTEdge->getLeftHand());
                queueOfQTEdges.push_back(currQTEdge->getLeftLeg());

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftLeg(), currQTEdge));
            }
        }
    }
}



void ShortestPathFinder2D::compute_gateway_QT_edges_of_end_face_and_its_prev_edges_by_BFS(FaceBU2D * qtStartFace, 
                                                                                          FaceBU2D * qtEndFace, 
                                                                                          vector<EdgeBU2D*>& gatewayQtEdgesOfEndFace, 
                                                                                          unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const
{
    int numOfGatewayQTEdge = 0;

    list<EdgeBU2D*> queueOfQTEdges;

    insert_three_bounding_QT_edge_into_container(qtStartFace, queueOfQTEdges);

    qtStartFace->isVisited(true);

    while (!queueOfQTEdges.empty())
    {
        EdgeBU2D* currQTEdge = queueOfQTEdges.front();
        queueOfQTEdges.pop_front();

        if (qtEndFace->isIncidentTo(currQTEdge) && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            gatewayQtEdgesOfEndFace.push_back(currQTEdge);
            numOfGatewayQTEdge++;
        }

        if (numOfGatewayQTEdge >= INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE)
        {
            break;
        }

        if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            FaceBU2D* currRightQTFace = currQTEdge->getRightFace();

            if ((!is_this_QT_face_visited(currRightQTFace)))
            {
                currRightQTFace->isVisited(true);

                queueOfQTEdges.push_back(currQTEdge->getRightHand());
                queueOfQTEdges.push_back(currQTEdge->getRightLeg());

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightLeg(), currQTEdge));
            }

            FaceBU2D* currLeftQTFace = currQTEdge->getLeftFace();

            if ((!is_this_QT_face_visited(currLeftQTFace)))
            {
                currLeftQTFace->isVisited(true);

                queueOfQTEdges.push_back(currQTEdge->getLeftHand());
                queueOfQTEdges.push_back(currQTEdge->getLeftLeg());

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftLeg(), currQTEdge));
            }
        }
    }

}

void ShortestPathFinder2D::compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_BFS_with_filtered_QT_edges(FaceBU2D* qtStartFace, 
                                                                                         FaceBU2D* qtEndFace, 
                                                                                         EdgeBU2D*& gatewayQtEdgeOfEndFace, 
                                                                                         unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper, 
                                                                                         const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const
{
    list<EdgeBU2D*> queueOfQTEdges;

    insert_three_bounding_QT_edge_into_container(qtStartFace, queueOfQTEdges);

    qtStartFace->isVisited(true);

    while (!queueOfQTEdges.empty())
    {
        EdgeBU2D* currQTEdge = queueOfQTEdges.front();
        queueOfQTEdges.pop_front();

        if (isThisQTEdgeInsideOfEllipse.at(currQTEdge)
            && qtEndFace->isIncidentTo(currQTEdge) 
            && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            gatewayQtEdgeOfEndFace = currQTEdge;
            break;
        }

        if (isThisQTEdgeInsideOfEllipse.at(currQTEdge) 
            && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            FaceBU2D* currRightQTFace = currQTEdge->getRightFace();

            if ((!is_this_QT_face_visited(currRightQTFace)))
            {
                currRightQTFace->isVisited(true);

                queueOfQTEdges.push_back(currQTEdge->getRightHand());
                queueOfQTEdges.push_back(currQTEdge->getRightLeg());

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightLeg(), currQTEdge));
            }

            FaceBU2D* currLeftQTFace = currQTEdge->getLeftFace();

            if ((!is_this_QT_face_visited(currLeftQTFace)))
            {
                currLeftQTFace->isVisited(true);

                queueOfQTEdges.push_back(currQTEdge->getLeftHand());
                queueOfQTEdges.push_back(currQTEdge->getLeftLeg());

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftLeg(), currQTEdge));
            }
        }
    }
}



void ShortestPathFinder2D::compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_adjusted_BFS_with_filtered_QT_edges(FaceBU2D * qtStartFace, FaceBU2D * qtEndFace, EdgeBU2D *& gatewayQtEdgeOfEndFace, unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper, const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const
{
    EntityAccessiblePriorityQ<EdgeBU2D*> queueOfQTEdges;

    insert_three_bounding_QT_edge_into_container(qtStartFace, queueOfQTEdges);

    qtStartFace->isVisited(true);

    while (!queueOfQTEdges.empty())
    {
        EdgeBU2D* currQTEdge = queueOfQTEdges.pop();

        if (isThisQTEdgeInsideOfEllipse.at(currQTEdge)
            && qtEndFace->isIncidentTo(currQTEdge)
            && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            gatewayQtEdgeOfEndFace = currQTEdge;
            break;
        }

        if (isThisQTEdgeInsideOfEllipse.at(currQTEdge)
            && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            FaceBU2D* currRightQTFace = currQTEdge->getRightFace();

            if ((!is_this_QT_face_visited(currRightQTFace)))
            {
                currRightQTFace->isVisited(true);

                insert_QT_edge_into_PQ(currQTEdge->getRightHand(), queueOfQTEdges);
                insert_QT_edge_into_PQ(currQTEdge->getRightLeg(), queueOfQTEdges);

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getRightLeg(), currQTEdge));
            }

            FaceBU2D* currLeftQTFace = currQTEdge->getLeftFace();

            if ((!is_this_QT_face_visited(currLeftQTFace)))
            {
                currLeftQTFace->isVisited(true);

                insert_QT_edge_into_PQ(currQTEdge->getLeftHand(), queueOfQTEdges);
                insert_QT_edge_into_PQ(currQTEdge->getLeftLeg(), queueOfQTEdges);

                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftHand(), currQTEdge));
                currToPrevQTEdgeMapper.insert(pair<EdgeBU2D*, EdgeBU2D*>(currQTEdge->getLeftLeg(), currQTEdge));
            }
        }
    }
}



void ShortestPathFinder2D::compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_DFS_with_filtered_QT_edges(FaceBU2D* qtStartFace, 
                                                                                         FaceBU2D* qtEndFace, 
                                                                                         EdgeBU2D*& gatewayQtEdgeOfEndFace, 
                                                                                         unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper, 
                                                                                         const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const
{
    list<EdgeBU2D*> stackOfQTEdges;
    insert_three_bounding_QT_edge_into_container(qtStartFace, stackOfQTEdges);

    qtStartFace->isVisited(true);

    unordered_map<FaceBU2D*, EdgeBU2D*> currFaceToPrevQTEdgeMapper;
    currFaceToPrevQTEdgeMapper.insert(make_pair(qtStartFace, (EdgeBU2D*)NULL));

    while (!stackOfQTEdges.empty())
    {
        EdgeBU2D* currQTEdge = stackOfQTEdges.back();
        stackOfQTEdges.pop_back();

        if (isThisQTEdgeInsideOfEllipse.at(currQTEdge) 
            && qtEndFace->isIncidentTo(currQTEdge)
            && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            EdgeBU2D* prevQTEdge = NULL;

            if (currQTEdge->getRightFace() == qtEndFace)
            {
                prevQTEdge = currFaceToPrevQTEdgeMapper.at(currQTEdge->getLeftFace());
            }
            else // (currQTEdge->getLeftFace() == qtEndFace)
            {
                prevQTEdge = currFaceToPrevQTEdgeMapper.at(currQTEdge->getRightFace());
            }

            currToPrevQTEdgeMapper.insert(make_pair(currQTEdge, prevQTEdge));

            gatewayQtEdgeOfEndFace = currQTEdge;
            break;
        }

        if (isThisQTEdgeInsideOfEllipse.at(currQTEdge) && is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            FaceBU2D* currRightQTFace = currQTEdge->getRightFace();
            FaceBU2D* currLeftQTFace = currQTEdge->getLeftFace();

            if ((!is_this_QT_face_visited(currRightQTFace)) && (is_this_QT_face_visited(currLeftQTFace)))
            {
                currRightQTFace->isVisited(true);
                currFaceToPrevQTEdgeMapper.insert(make_pair(currRightQTFace, currQTEdge));

                stackOfQTEdges.push_back(currQTEdge->getRightHand());
                stackOfQTEdges.push_back(currQTEdge->getRightLeg());

                EdgeBU2D* prevQTEdge = currFaceToPrevQTEdgeMapper.at(currLeftQTFace);

                if (prevQTEdge != NULL)
                {
                    currToPrevQTEdgeMapper.insert(make_pair(currQTEdge, prevQTEdge));
                }
            }

            if ((!is_this_QT_face_visited(currLeftQTFace)) && (is_this_QT_face_visited(currRightQTFace)))
            {
                currLeftQTFace->isVisited(true);
                currFaceToPrevQTEdgeMapper.insert(make_pair(currLeftQTFace, currQTEdge));

                stackOfQTEdges.push_back(currQTEdge->getLeftHand());
                stackOfQTEdges.push_back(currQTEdge->getLeftLeg());

                EdgeBU2D* prevQTEdge = currFaceToPrevQTEdgeMapper.at(currRightQTFace);

                if (prevQTEdge != NULL)
                {
                    currToPrevQTEdgeMapper.insert(make_pair(currQTEdge, prevQTEdge));
                }
            }

        }
    }
}

void ShortestPathFinder2D::compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_DFS(FaceBU2D* qtStartFace, FaceBU2D* qtEndFace, EdgeBU2D*& gatewayQtEdgeOfEndFace, unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const
{
    list<EdgeBU2D*> stackOfQTEdges;
    insert_three_bounding_QT_edge_into_container(qtStartFace, stackOfQTEdges);

    qtStartFace->isVisited(true);

    unordered_map<FaceBU2D*, EdgeBU2D*> currFaceToPrevQTEdgeMapper;
    currFaceToPrevQTEdgeMapper.insert(make_pair(qtStartFace, (EdgeBU2D*) NULL));

    while (!stackOfQTEdges.empty())
    {
        EdgeBU2D* currQTEdge = stackOfQTEdges.back();
        stackOfQTEdges.pop_back();

        if (qtEndFace->isIncidentTo(currQTEdge))
        {
            EdgeBU2D* prevQTEdge = NULL;

            if (currQTEdge->getRightFace() == qtEndFace)
            {
                prevQTEdge = currFaceToPrevQTEdgeMapper.at(currQTEdge->getLeftFace());
            }
            else // (currQTEdge->getLeftFace() == qtEndFace)
            {
                prevQTEdge = currFaceToPrevQTEdgeMapper.at(currQTEdge->getRightFace());
            }
    
            currToPrevQTEdgeMapper.insert(make_pair(currQTEdge, prevQTEdge));

            gatewayQtEdgeOfEndFace = currQTEdge;
            break;
        }

        if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            FaceBU2D* currRightQTFace = currQTEdge->getRightFace();
            FaceBU2D* currLeftQTFace  = currQTEdge->getLeftFace();

            if ((!is_this_QT_face_visited(currRightQTFace)) && (is_this_QT_face_visited(currLeftQTFace)))
            {
                currRightQTFace->isVisited(true);
                currFaceToPrevQTEdgeMapper.insert(make_pair(currRightQTFace, currQTEdge));

                stackOfQTEdges.push_back(currQTEdge->getRightHand());
                stackOfQTEdges.push_back(currQTEdge->getRightLeg());

                EdgeBU2D* prevQTEdge = currFaceToPrevQTEdgeMapper.at(currLeftQTFace);
               
                if (prevQTEdge != NULL)
                {
                    currToPrevQTEdgeMapper.insert(make_pair(currQTEdge, prevQTEdge));
                }
            }

            if ((!is_this_QT_face_visited(currLeftQTFace)) && (is_this_QT_face_visited(currRightQTFace)))
            {
                currLeftQTFace->isVisited(true);
                currFaceToPrevQTEdgeMapper.insert(make_pair(currLeftQTFace, currQTEdge));

                stackOfQTEdges.push_back(currQTEdge->getLeftHand());
                stackOfQTEdges.push_back(currQTEdge->getLeftLeg());

                EdgeBU2D* prevQTEdge = currFaceToPrevQTEdgeMapper.at(currRightQTFace);
                
                if (prevQTEdge != NULL)
                {
                    currToPrevQTEdgeMapper.insert(make_pair(currQTEdge, prevQTEdge));
                }
            }

        }
    }
}


void ShortestPathFinder2D::compute_gateway_QT_edge_of_end_face_and_its_prev_edges(FaceBU2D* qtStartFace, 
                                                                                  FaceBU2D* qtEndFace, 
                                                                                  EdgeBU2D*& gatewayQtEdgeOfEndFace, 
                                                                                  unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const
{
    list<FaceBU2D*> stackOfQTFaces;
    stackOfQTFaces.push_back(qtStartFace);

    FaceBU2D* prevQTFace = NULL;

    while (!stackOfQTFaces.empty())
    {
        FaceBU2D* currQTFace = stackOfQTFaces.back();
        stackOfQTFaces.pop_back();

        EdgeBU2D* bridgeQTEdge = get_sharing_QT_edge_of_two_QT_faces(prevQTFace, currQTFace);


        if (currQTFace == qtEndFace)
        {
            rg_dList<EdgeBU2D*> boundingQTEdges;
            currQTFace->getBoundingEdges(boundingQTEdges);

            boundingQTEdges.reset4Loop();
            while (boundingQTEdges.setNext4Loop())
            {
                EdgeBU2D* currQTEdge = boundingQTEdges.getEntity();
                
                if (currQTEdge->getRightFace()->isVisited() || currQTEdge->getLeftFace()->isVisited())
                {
                    gatewayQtEdgeOfEndFace = currQTEdge;
                    break;
                }
            }

            break;
        }

        if (currQTFace->isVisited())
        {
            continue;
        }

        currQTFace->isVisited(true);

        rg_dList<EdgeBU2D*> boundingQTEdges;
        currQTFace->getBoundingEdges(boundingQTEdges);

        boundingQTEdges.reset4Loop();
        while (boundingQTEdges.setNext4Loop())
        {
            EdgeBU2D* currQTEdge = boundingQTEdges.getEntity();

            if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
            {
                if (currQTEdge->getRightFace() == currQTFace)
                {
                    stackOfQTFaces.push_back(currQTEdge->getLeftFace());
                }
                else //(currQTEdge->getLeftFace() == currQTFace)
                {
                    stackOfQTFaces.push_back(currQTEdge->getRightFace());
                }
            }
        }

        prevQTFace = currQTFace;
    }
}




EdgeBU2D* ShortestPathFinder2D::get_sharing_QT_edge_of_two_QT_faces(FaceBU2D* qtFace1, FaceBU2D* qtFace2) const
{
    EdgeBU2D* sharingQTEdge = NULL;

    rg_dList<EdgeBU2D*> boundingEdgesOfQTFace1;

    qtFace1->getBoundingEdges(boundingEdgesOfQTFace1);

    boundingEdgesOfQTFace1.reset4Loop();

    while (boundingEdgesOfQTFace1.setNext4Loop())
    {
        EdgeBU2D* currQTEdge = boundingEdgesOfQTFace1.getEntity();

        if ((currQTEdge->getRightFace() == qtFace2) || (currQTEdge->getLeftFace() == qtFace2) )
        {
            sharingQTEdge = currQTEdge;
            break;
        }
    }

    return sharingQTEdge;
}



void ShortestPathFinder2D::back_trace_to_QT_start_face_and_find_topological_path(EdgeBU2D* gatewayQtEdgeOfEndFace, 
                                                                                 unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper, 
                                                                                 list<EdgeBU2D*>& topologicalPath) const
{
    EdgeBU2D* currQTEdge = gatewayQtEdgeOfEndFace;

    while (currQTEdge != NULL)
    {
        topologicalPath.push_front(currQTEdge);
        
        unordered_map<EdgeBU2D*, EdgeBU2D*>::iterator it_CurrQTEdge = currToPrevQTEdgeMapper.find(currQTEdge);

        if (it_CurrQTEdge != currToPrevQTEdgeMapper.end())
        {
            currQTEdge = (*it_CurrQTEdge).second;
        }
        else
        {
            currQTEdge = NULL;
        }
    }
}


void ShortestPathFinder2D::back_trace_to_QT_start_face_and_find_topological_paths(vector<EdgeBU2D*>& gatewayQtEdgesOfEndFace, 
                                                                                  unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper, 
                                                                                  vector<list<EdgeBU2D*>>& topologicalPaths) const
{
    for (int i = 0; i < gatewayQtEdgesOfEndFace.size(); i++)
    {
        list<EdgeBU2D*> topologicalPath;

        EdgeBU2D* gatewayQtEdgeOfEndFace = gatewayQtEdgesOfEndFace.at(i);

        back_trace_to_QT_start_face_and_find_topological_path(gatewayQtEdgeOfEndFace, currToPrevQTEdgeMapper, topologicalPath);

        topologicalPaths.at(i) = topologicalPath;
    }
}



void ShortestPathFinder2D::find_one_geodesic_path_from_QT_topological_path(const list<EdgeBU2D*>& topologicalPath, 
                                                                           BetaUniverse2D& betaUniverse, 
                                                                           ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                           list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                           double& totalPathDistacne) 
{
    if (topologicalPath.size() == 0)
    {
        return;
    }

#ifdef OBSTACLE_PAIR_ON_QT_EDGE

    vector<pair<VertexBU2D*, VertexBU2D*>> obstaclePairs;
    make_obstacle_pool(obstaclePairs, topologicalPath);

    generate_solution_pool_graph(obstaclePairs, betaUniverse, solutionPoolGraph);

#else // OBSTACLE_PAIR_ON_QT_EDGE
    //1. Make obstacle pool 
    vector<VertexBU2D*> obstaclePool;
    make_obstacle_pool(obstaclePool, topologicalPath);

    //2. generate solution pool
    //generate_solution_pool_graph(obstaclePool, betaUniverse, solutionPoolGraph);
    generate_solution_pool_graph_with_topological_path(obstaclePool, topologicalPath, betaUniverse, solutionPoolGraph);

#endif
    //3. dijkstra algorithm in solution pool graph
    find_geodesic_path_by_dijkstra(solutionPoolGraph);
    
    //4. back tracking
    if (solutionPoolGraph.get_end_vtx()->get_prev_edge() == NULL)
    {
        return;
    }

    back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(), 
                                                      solutionPoolGraph.get_end_vtx(),
                                                      geodesicPath, 
                                                      totalPathDistacne);
}


void ShortestPathFinder2D::find_one_geodesic_path_from_QT_topological_path_for_all_path_search(const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, list<EdgeForSolutionPoolGraph*>& geodesicPath, double & totalPathDistacne)
{
      if (topologicalPath.size() == 0)
    {
        return;
    }

#ifdef OBSTACLE_PAIR_ON_QT_EDGE

    vector<pair<VertexBU2D*, VertexBU2D*>> obstaclePairs;
    make_obstacle_pool(obstaclePairs, topologicalPath);

    generate_solution_pool_graph(obstaclePairs, betaUniverse, solutionPoolGraph);

#else // OBSTACLE_PAIR_ON_QT_EDGE
    //1. Make obstacle pool 
    vector<VertexBU2D*> obstaclePool;
    make_obstacle_pool(obstaclePool, topologicalPath);

    //2. generate solution pool
    generate_solution_pool_graph(obstaclePool, betaUniverse, solutionPoolGraph);
    //generate_solution_pool_graph_with_topological_path(obstaclePool, topologicalPath, betaUniverse, solutionPoolGraph);

#endif
    //3. dijkstra algorithm in solution pool graph
    find_geodesic_path_by_dijkstra(solutionPoolGraph);
    
    //4. back tracking
    back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(), 
                                                   solutionPoolGraph.get_end_vtx(),
                                                   geodesicPath, 
                                                   totalPathDistacne);
}


void ShortestPathFinder2D::find_one_geodesic_path_from_candidate_obstacles(const unordered_map<VertexBU2D*, bool>& isThisQTVertexInsideOfEllipse, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, list<EdgeForSolutionPoolGraph*>& geodesicPath, double & totalPathDistacne)
{
    if (isThisQTVertexInsideOfEllipse.size() == 0)
    {
        return;
    }

    //1. Make obstacle pool 
    vector<VertexBU2D*> obstaclePool;
    make_obstacle_pool(obstaclePool, isThisQTVertexInsideOfEllipse);

    //2. generate solution pool
    generate_solution_pool_graph(obstaclePool, betaUniverse, solutionPoolGraph);

    //3. dijkstra algorithm in solution pool graph
    find_geodesic_path_by_dijkstra(solutionPoolGraph);
    
    //4. back tracking
    back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(), 
                                                   solutionPoolGraph.get_end_vtx(),
                                                   geodesicPath, 
                                                   totalPathDistacne);
}




void ShortestPathFinder2D::find_several_geodesic_paths_from_QT_topological_paths(const vector<list<EdgeBU2D*>>& topologicalPaths, 
                                                                                 BetaUniverse2D & betaUniverse, 
                                                                                 vector<ShortestPathSolutionPoolGraph>& solutionPoolGraphs, 
                                                                                 vector<list<EdgeForSolutionPoolGraph*>>& geodesicPaths, 
                                                                                 vector<double>& totalPathDistances)
{
    for (int i = 0; i < m_TopologicalPathsByTPT.size(); i++)
    {
        ShortestPathSolutionPoolGraph*  solutionPool = create_solution_pool_of_BFS(ShortestPathSolutionPoolGraph(), i);
        list<EdgeForSolutionPoolGraph*> geodesicPath;
        
        double                          totalDistance = DBL_MAX;

        find_one_geodesic_path_from_QT_topological_path(topologicalPaths.at(i),
                                                        betaUniverse,
                                                        *solutionPool,
                                                        geodesicPath,
                                                        totalDistance);

        geodesicPaths.at(i)      = geodesicPath;
        totalPathDistances.at(i) = totalDistance;
    }
}


void ShortestPathFinder2D::find_one_geodesic_path_from_filltered_disks_by_ellipse(const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, 
                                                                                  BetaUniverse2D& betaUniverse, 
                                                                                  ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                                  list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                                  double& totalPathDistance) 
{
    //1. Make obstacle pool 
    vector<VertexBU2D*> obstaclePool;
    make_obstacle_pool(obstaclePool, isThisQTEdgeInsideOfEllipse);

    //2. generate solution pool
    generate_solution_pool_graph(obstaclePool, betaUniverse, solutionPoolGraph);

    //3. dijkstra algorithm in solution pool graph
    find_geodesic_path_by_dijkstra(solutionPoolGraph);

    //4. back tracking
    back_trace_to_start_vtx_and_find_geodesic_path(solutionPoolGraph.get_start_vtx(),
                                                      solutionPoolGraph.get_end_vtx(),
                                                      geodesicPath,
                                                      totalPathDistance);
}


void ShortestPathFinder2D::get_all_QT_vertices_from_topological_path(set<VertexBU2D*>& qtVertice, const list<EdgeBU2D*>& qtEdge) const
{
    for (list<EdgeBU2D*>::const_iterator it_QTEdge = qtEdge.begin(); it_QTEdge != qtEdge.end(); it_QTEdge++)
    {
        EdgeBU2D* currQTEdge = *it_QTEdge;

        qtVertice.insert(currQTEdge->getStartVertex());
        qtVertice.insert(currQTEdge->getEndVertex());
    }
}



void ShortestPathFinder2D::generate_solution_pool_graph(const vector<VertexBU2D*>& obstaclePool, 
                                                        BetaUniverse2D& betaUniverse, 
                                                        ShortestPathSolutionPoolGraph& solutionPoolGraph)
{
    //1. Generate solution pool
    unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>> mapperForQTVtxToVtxOnDisk;
    generate_graph_for_obstacles_by_using_QT_and_BS(obstaclePool, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);

    //2 Add start and end vtx into the graph
    insert_start_end_point_into_the_graph(m_StartPt, m_EndPt, obstaclePool, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}



void ShortestPathFinder2D::generate_solution_pool_graph_with_topological_path(const vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph)
{
    m_QTEdgesOfTopologicalPath.insert(topologicalPath.begin(), topologicalPath.end());

    //1. Generate solution pool
    unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>> mapperForQTVtxToVtxOnDisk;
    generate_graph_for_obstacles_by_using_QT_and_BS_with_topological_path(obstaclePool, topologicalPath, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);

    //2 Add start and end vtx into the graph
    insert_start_end_point_into_the_graph_with_topological_path(m_StartPt, m_EndPt, obstaclePool, topologicalPath, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);

    m_QTEdgesOfTopologicalPath.clear();
}


void ShortestPathFinder2D::generate_solution_pool_graph(const vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, BetaUniverse2D& betaUniverse, ShortestPathSolutionPoolGraph& solutionPoolGraph)
{
    //1. Generate solution pool
    unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>> mapperForQTVtxToVtxOnDisk;
    generate_graph_for_obstacles_by_using_QT_and_BS(obstaclePairs, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);

    //2 Add start and end vtx into the graph
    vector<VertexBU2D*> obstaclePool;
    make_set_of_obstacles(obstaclePairs, obstaclePool);
    insert_start_end_point_into_the_graph(m_StartPt, m_EndPt, obstaclePool, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}


void ShortestPathFinder2D::generate_visible_solution_pool_graph_from_current_position(const rg_Point2D & currPosition, const vector<VertexBU2D*>& obstaclePool, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph) const
{
    //1. Generate solution pool
    unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>> mapperForQTVtxToVtxOnDisk;

    VertexForSolutionPoolGraph* startVtx = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());

    startVtx->set_tangent_point_to_disk(currPosition);
    startVtx->set_accumulate_length_from_source_vtx(DBL_MAX);

    solutionPoolGraph.set_start_vtx(startVtx);

    insert_a_point_into_the_graph(startVtx, obstaclePool, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}


void ShortestPathFinder2D::make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, const unordered_map<VertexBU2D*, bool>& isThisQTVertexInsideOfEllipse) const
{

    for (unordered_map<VertexBU2D*, bool>::const_iterator it_Vertex = isThisQTVertexInsideOfEllipse.begin();
        it_Vertex != isThisQTVertexInsideOfEllipse.end();
        it_Vertex++)
    {
        const pair<VertexBU2D*, bool>& pairOfVertex = *it_Vertex;

        VertexBU2D* currVertex = pairOfVertex.first;
        bool       isInEllipse = pairOfVertex.second;

        if (isInEllipse)
        {
            obstaclePool.push_back(currVertex);
        }
    }
}

void ShortestPathFinder2D::make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, 
                                              const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const
{
    for (unordered_map<EdgeBU2D*, bool>::const_iterator it_Edge = isThisQTEdgeInsideOfEllipse.begin();
        it_Edge != isThisQTEdgeInsideOfEllipse.end();
        it_Edge++)
    {
        const pair<EdgeBU2D*, bool>& pairOfEdge = *it_Edge;

        EdgeBU2D* currEdge    = pairOfEdge.first;
        bool      isInEllipse = pairOfEdge.second;

        if (isInEllipse)
        {
            if (!currEdge->getStartVertex()->isVirtual())
            {
                obstaclePool.push_back(currEdge->getStartVertex());
            }

            if (!currEdge->getEndVertex()->isVirtual())
            {
                obstaclePool.push_back(currEdge->getEndVertex());
            }
        }
    }
}


void ShortestPathFinder2D::make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& qtEdge) const
{

    set<VertexBU2D*> obstacleSet;

    for (list<EdgeBU2D*>::const_iterator it_QTEdge = qtEdge.begin();
        it_QTEdge != qtEdge.end();
        it_QTEdge++)
    {
        EdgeBU2D* currQTEdge = *it_QTEdge;

        switch (currQTEdge->isVirtual())
        {
        case true: //virtual edge should include neighbor disks
        {
            list<VertexBU2D*> neighborQTVertex;
            if (currQTEdge->getStartVertex()->isVirtual())
            {
                get_neighbor_obstacles_interseted_with_this_obstacle(neighborQTVertex, currQTEdge->getEndVertex(), m_BetaComponentFinder);
                obstacleSet.insert(currQTEdge->getEndVertex());
            }
            else //currQTEdge->getEndVertex()->isVirtual()
            {
                get_neighbor_obstacles_interseted_with_this_obstacle(neighborQTVertex, currQTEdge->getStartVertex(), m_BetaComponentFinder);
                obstacleSet.insert(currQTEdge->getStartVertex());
            }

            for (list<VertexBU2D*>::const_iterator it_NeighborVertex = neighborQTVertex.begin();
                it_NeighborVertex != neighborQTVertex.end();
                it_NeighborVertex++)
            {
                VertexBU2D* currNeighborVertex = *it_NeighborVertex;

                if (!currNeighborVertex->isVirtual())
                {
                    obstacleSet.insert(currNeighborVertex);
                }
            }
        }

        default:
            break;
        }

        VertexBU2D* currQTStartVertex = (currQTEdge)->getStartVertex();
        VertexBU2D* currQTEndVertex = (currQTEdge)->getEndVertex();
        VertexBU2D* currQTRightVertex = (currQTEdge)->getVertexOfRightHand();
        VertexBU2D* currQTLeftVertex = (currQTEdge)->getVertexOfLeftHand();

        if (!currQTStartVertex->isVirtual())
        {
            obstacleSet.insert(currQTStartVertex);
        }

        if (!currQTEndVertex->isVirtual())
        {
            obstacleSet.insert(currQTEndVertex);
        }

        if (!currQTRightVertex->isVirtual())
        {
            obstacleSet.insert(currQTRightVertex);
        }

        if (!currQTLeftVertex->isVirtual())
        {
            obstacleSet.insert(currQTLeftVertex);
        }
    }

    obstaclePool.insert(obstaclePool.end(), obstacleSet.begin(), obstacleSet.end());

#ifdef _DEBUG
    sort(obstaclePool.begin(), obstaclePool.end(), ShortestPathFinder2D::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order);
#endif


    /*
    set<VertexBU2D*> obstacleSet;

    for (list<EdgeBU2D*>::const_iterator it_QTEdge = qtEdge.begin();
        it_QTEdge != qtEdge.end();
        it_QTEdge++)
    {
        VertexBU2D* currQTStartVertex = (*it_QTEdge)->getStartVertex();
        VertexBU2D* currQTEndVertex = (*it_QTEdge)->getEndVertex();
        VertexBU2D* currQTRightVertex = (*it_QTEdge)->getVertexOfRightHand();
        VertexBU2D* currQTLeftVertex = (*it_QTEdge)->getVertexOfLeftHand();


        if (!currQTStartVertex->isVirtual())
        {
            obstacleSet.insert(currQTStartVertex);
        }

        if (!currQTEndVertex->isVirtual())
        {
            obstacleSet.insert(currQTEndVertex);
        }

        if (!currQTRightVertex->isVirtual())
        {
            obstacleSet.insert(currQTRightVertex);
        }

        if (!currQTLeftVertex->isVirtual())
        {
            obstacleSet.insert(currQTLeftVertex);
        }
    }

    obstaclePool.insert(obstaclePool.end(), obstacleSet.begin(), obstacleSet.end());

#ifdef _DEBUG
    sort(obstaclePool.begin(), obstaclePool.end(), ShortestPathFinder2D::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order);
#endif
    */


    /*
    set<VertexBU2D*> obstacleSet;

    for (list<EdgeBU2D*>::const_iterator it_QTEdge = qtEdge.begin();
        it_QTEdge != qtEdge.end();
        it_QTEdge++)
    {
        VertexBU2D* currQTStartVertex = (*it_QTEdge)->getStartVertex();
        VertexBU2D* currQTEndVertex   = (*it_QTEdge)->getEndVertex();

        if (!currQTStartVertex->isVirtual())
        {
            obstacleSet.insert(currQTStartVertex);
        }

        if (!currQTEndVertex->isVirtual())
        {
            obstacleSet.insert(currQTEndVertex);
        }
    }

    obstaclePool.insert(obstaclePool.end(), obstacleSet.begin(), obstacleSet.end());

    */
}



void ShortestPathFinder2D::make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, FaceBU2D* qtFace) const
{
    set<VertexBU2D*> obstcleSet;

    rg_dList<VertexBU2D*> boundingVertices;
    qtFace->getBoundingVertices(boundingVertices);

    boundingVertices.reset4Loop();
    while (boundingVertices.setNext4Loop())
    {
        VertexBU2D* currQTVertex = boundingVertices.getEntity();
        
        rg_dList<VertexBU2D*> currNItsNeighbor;
        currQTVertex->getNeighborVertices(currNItsNeighbor);
        currNItsNeighbor.add(currQTVertex);

        currNItsNeighbor.reset4Loop();
        while (currNItsNeighbor.setNext4Loop())
        {
            VertexBU2D* neighborQTVertex = currNItsNeighbor.getEntity();

            if (!neighborQTVertex->isVirtual())
            {
                obstcleSet.insert(neighborQTVertex);
            }
        }     
    }

    obstaclePool.insert(obstaclePool.end(), obstcleSet.begin(), obstcleSet.end());

#ifdef _DEBUG
    sort(obstaclePool.begin(), obstaclePool.end(), ShortestPathFinder2D::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order);
#endif
}



void ShortestPathFinder2D::convert_filltered_QT_edges_to_topological_path(const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, list<EdgeBU2D*>& topologicalPath) const
{
    for (unordered_map<EdgeBU2D*, bool>::const_iterator it_EdgeNBoolean = isThisQTEdgeInsideOfEllipse.begin();
        it_EdgeNBoolean != isThisQTEdgeInsideOfEllipse.end();
        it_EdgeNBoolean++)
    {
        const pair<EdgeBU2D*, bool>& currPair = *it_EdgeNBoolean;

        EdgeBU2D* currEdge = currPair.first;
        bool      currFilteredStatus = currPair.second;

        if (currFilteredStatus)
        {
            topologicalPath.push_back(currEdge);
        }
    }
}


void ShortestPathFinder2D::make_obstacle_pool(vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, const list<EdgeBU2D*>& qtEdge) const
{
    set<pair<VertexBU2D*, VertexBU2D*>> obstaclePairSet;

    for (list<EdgeBU2D*>::const_iterator it_QTEdge = qtEdge.begin();
        it_QTEdge != qtEdge.end();
        it_QTEdge++)
    {
        EdgeBU2D* currQTEdge            = (*it_QTEdge);
        EdgeBU2D* rightHandOfCurrQTEdge = currQTEdge->getRightHand();
        EdgeBU2D* leftHandOfCurrQTEdge  = currQTEdge->getLeftHand();
        EdgeBU2D* rightLegOfCurrQTEdge  = currQTEdge->getRightLeg();
        EdgeBU2D* leftLegOfCurrQTEdge   = currQTEdge->getLeftLeg();

        if (!currQTEdge->isVirtual())
        {
            if (currQTEdge->getStartVertex() < currQTEdge->getEndVertex())
            {
                obstaclePairSet.insert(make_pair(currQTEdge->getStartVertex(), currQTEdge->getEndVertex()));
            }
            else
            {
                obstaclePairSet.insert(make_pair(currQTEdge->getEndVertex(), currQTEdge->getStartVertex()));
            }
        }

        if (!rightHandOfCurrQTEdge->isVirtual())
        {
            if (rightHandOfCurrQTEdge->getStartVertex() < rightHandOfCurrQTEdge->getEndVertex())
            {
                obstaclePairSet.insert(make_pair(rightHandOfCurrQTEdge->getStartVertex(), rightHandOfCurrQTEdge->getEndVertex()));
            }
            else
            {
                obstaclePairSet.insert(make_pair(rightHandOfCurrQTEdge->getEndVertex(), rightHandOfCurrQTEdge->getStartVertex()));
            }
        }

        if (!leftHandOfCurrQTEdge->isVirtual())
        {
            if (leftHandOfCurrQTEdge->getStartVertex() < leftHandOfCurrQTEdge->getEndVertex())
            {
                obstaclePairSet.insert(make_pair(leftHandOfCurrQTEdge->getStartVertex(), leftHandOfCurrQTEdge->getEndVertex()));
            }
            else
            {
                obstaclePairSet.insert(make_pair(leftHandOfCurrQTEdge->getEndVertex(), leftHandOfCurrQTEdge->getStartVertex()));
            }
        }

        if (!rightLegOfCurrQTEdge->isVirtual())
        {
            if (rightLegOfCurrQTEdge->getStartVertex() < rightLegOfCurrQTEdge->getEndVertex())
            {
                obstaclePairSet.insert(make_pair(rightLegOfCurrQTEdge->getStartVertex(), rightLegOfCurrQTEdge->getEndVertex()));
            }
            else
            {
                obstaclePairSet.insert(make_pair(rightLegOfCurrQTEdge->getEndVertex(), rightLegOfCurrQTEdge->getStartVertex()));
            }
        }

        if (!leftLegOfCurrQTEdge->isVirtual())
        {
            if (leftLegOfCurrQTEdge->getStartVertex() < leftLegOfCurrQTEdge->getEndVertex())
            {
                obstaclePairSet.insert(make_pair(leftLegOfCurrQTEdge->getStartVertex(), leftLegOfCurrQTEdge->getEndVertex()));
            }
            else
            {
                obstaclePairSet.insert(make_pair(leftLegOfCurrQTEdge->getEndVertex(), leftLegOfCurrQTEdge->getStartVertex()));
            }
        }
    }

    obstaclePairs.insert(obstaclePairs.end(), obstaclePairSet.begin(), obstaclePairSet.end());
}



void ShortestPathFinder2D::generate_graph_for_obstacles_by_using_QT_and_BS(const vector<VertexBU2D*>& obstaclePool,
                                                                           BetaUniverse2D& betaUniverse, 
                                                                           ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                           unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk)
{
    if ((obstaclePool.size() == 0) || (obstaclePool.size() == 1))
    {
        return;
    }

    for (int i = 0; i < (obstaclePool.size() - 1); i++)
    {
        for (int j = i + 1; j < obstaclePool.size(); j++)
        {
            VertexBU2D* currObs1 = obstaclePool.at(i);
            VertexBU2D* currObs2 = obstaclePool.at(j);

            generate_graph_incrementally_by_adding_two_obstacles_into_the_graph(currObs1,
                                                                                currObs2,
                                                                                m_BetaComponentFinder,
                                                                                betaUniverse,
                                                                                solutionPoolGraph, 
                                                                                mapperForQTVtxToVtxOnDisk);

        }
    }
}


void ShortestPathFinder2D::generate_graph_for_obstacles_by_using_QT_and_BS_with_topological_path(const vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk)
{
    if ((obstaclePool.size() == 0) || (obstaclePool.size() == 1))
    {
        return;
    }

    for (int i = 0; i < (obstaclePool.size() - 1); i++)
    {
        for (int j = i + 1; j < obstaclePool.size(); j++)
        {
            if ((i == 1) && (j == 5))
            {
                cout << endl;
            }
            
            VertexBU2D* currObs1 = obstaclePool.at(i);
            VertexBU2D* currObs2 = obstaclePool.at(j);

            generate_graph_incrementally_by_adding_two_obstacles_into_the_graph_with_topological_path(currObs1,
                                                                                currObs2,
                                                                                topologicalPath,
                                                                                m_BetaComponentFinder,
                                                                                betaUniverse,
                                                                                solutionPoolGraph, 
                                                                                mapperForQTVtxToVtxOnDisk);

        }
    }
}


void ShortestPathFinder2D::generate_graph_for_obstacles_by_using_QT_and_BS(const vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, BetaUniverse2D& betaUniverse, ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk)
{
    for (int i = 0; i < obstaclePairs.size(); i++)
    {
        VertexBU2D* currObs1 = obstaclePairs.at(i).first;
        VertexBU2D* currObs2 = obstaclePairs.at(i).second;
        
        generate_graph_incrementally_by_adding_two_obstacles_into_the_graph(currObs1,
                                                                            currObs2,
                                                                            m_BetaComponentFinder,
                                                                            betaUniverse,
                                                                            solutionPoolGraph, 
                                                                            mapperForQTVtxToVtxOnDisk);

    }
}





void ShortestPathFinder2D::make_graph_of_tangent_lines(VertexBU2D* obs1, 
                                                       VertexBU2D* obs2, 
                                                       const BetaComponentFinder& betaComponentFinder, 
                                                       BetaUniverse2D& betaUniverse, 
                                                       ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                       vector<VertexForSolutionPoolGraph*>& verticesOnObstacle1,
                                                       vector<VertexForSolutionPoolGraph*>& verticesOnObstacle2)
{
    //1. find tangent line
    vector<rg_ImplicitEquation> tangentLineEqs;   
    find_tangent_line_equations_of_two_obstacles(obs1, obs2, betaComponentFinder, tangentLineEqs);
   
    //2. compute tangent point and filter the tangent line by intersection with other disks using QT
    vector<rg_ImplicitEquation> filteredTangentLineEqs;
    vector<pair<rg_Point2D, rg_Point2D>> listOfTwoPtsOnObstacle1And2;

    filter_tangent_line_equations_and_points_by_disks_vers1(tangentLineEqs,
                                             filteredTangentLineEqs, 
                                             listOfTwoPtsOnObstacle1And2, 
                                             betaUniverse, 
                                             obs1, obs2);

    //3. make verices and edges for graph
    make_vertices_N_tangent_line_edges(filteredTangentLineEqs, 
                                       listOfTwoPtsOnObstacle1And2, 
                                       verticesOnObstacle1, 
                                       verticesOnObstacle2, 
                                       solutionPoolGraph);
}



void ShortestPathFinder2D::make_graph_of_tangent_lines_with_topological_path(VertexBU2D * obs1, 
                                                                             VertexBU2D * obs2, 
                                                                             const list<EdgeBU2D*>& topologicalPath, 
                                                                             const BetaComponentFinder & betaComponentFinder, 
                                                                             BetaUniverse2D & betaUniverse,
                                                                             ShortestPathSolutionPoolGraph & solutionPoolGraph, 
                                                                             vector<VertexForSolutionPoolGraph*>& verticesOnObstacle1, 
                                                                             vector<VertexForSolutionPoolGraph*>& verticesOnObstacle2)
{
     //1. find tangent line
    vector<rg_ImplicitEquation> tangentLineEqs;   
    find_tangent_line_equations_of_two_obstacles(obs1, obs2, betaComponentFinder, tangentLineEqs);
   
    //2. compute tangent point and filter the tangent line by intersection with other disks using QT
    vector<rg_ImplicitEquation> filteredTangentLineEqs;
    vector<pair<rg_Point2D, rg_Point2D>> listOfTwoPtsOnObstacle1And2;

    /*
    filter_tangent_line_equations_and_points_by_disks(tangentLineEqs, 
                                                      filteredTangentLineEqs, 
                                                      listOfTwoPtsOnObstacle1And2, 
                                                      betaUniverse, 
                                                      obs1, obs2);
    */
    
    filter_tangent_line_equations_and_points_by_disks_vers1(tangentLineEqs,
                                                            filteredTangentLineEqs, 
                                                            listOfTwoPtsOnObstacle1And2, 
                                                            betaUniverse, 
                                                            obs1, obs2);
                                                          

    vector<rg_ImplicitEquation> filteredTangentLineEqsByTopoloticalPath;
    vector<pair<rg_Point2D, rg_Point2D>> filteredListOfTwoPtsOnObstacle1And2ByTopoloticalPath;

    filter_tangent_line_equations_and_points_by_topological_path(filteredTangentLineEqs, 
                                                                 listOfTwoPtsOnObstacle1And2, 
                                                                 filteredTangentLineEqsByTopoloticalPath,
                                                                 filteredListOfTwoPtsOnObstacle1And2ByTopoloticalPath,
                                                                 betaUniverse,
                                                                 topologicalPath);

    //3. make verices and edges for graph
    make_vertices_N_tangent_line_edges(filteredTangentLineEqsByTopoloticalPath,
                                       filteredListOfTwoPtsOnObstacle1And2ByTopoloticalPath, 
                                       verticesOnObstacle1, 
                                       verticesOnObstacle2, 
                                       solutionPoolGraph);
}

void ShortestPathFinder2D::make_graph_of_tangent_lines(VertexBU2D* obstacle, 
                                                      VertexForSolutionPoolGraph* pointVtx, 
                                                      BetaUniverse2D& betaUniverse,
                                                      ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                      vector<VertexForSolutionPoolGraph*>& verticesOnObstacle) const
{
    //1. find tangent line
    vector<rg_ImplicitEquation> tangentLineEqs;   
    find_tangent_line_equations_of_point_and_obstacle(pointVtx, obstacle, tangentLineEqs);
   
    //2. compute tangent point and filter the tangent line by intersection with other disks using QT
    vector<rg_ImplicitEquation> filteredTangentLineEqs;
    vector<rg_Point2D> listOfPtsOnObstacle;

    filter_tangent_line_equations_and_points_by_disks_vers1(tangentLineEqs,
                                             filteredTangentLineEqs, 
                                             listOfPtsOnObstacle,
                                             betaUniverse, 
                                             pointVtx, obstacle);


    //3. make verices and edges for graph
    make_vertices_N_tangent_line_edges(filteredTangentLineEqs, 
                                       listOfPtsOnObstacle, 
                                       pointVtx,
                                       verticesOnObstacle,
                                       solutionPoolGraph);
}


void ShortestPathFinder2D::make_graph_of_tangent_lines_with_topological_path(VertexBU2D * obstacle, VertexForSolutionPoolGraph * pointVtx, const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, vector<VertexForSolutionPoolGraph*>& verticesOnObstacle) const
{ 
    //1. find tangent line
    vector<rg_ImplicitEquation> tangentLineEqs;   
    find_tangent_line_equations_of_point_and_obstacle(pointVtx, obstacle, tangentLineEqs);
   
    //2. compute tangent point and filter the tangent line by intersection with other disks using QT
    vector<rg_ImplicitEquation> filteredTangentLineEqs;
    vector<rg_Point2D> listOfPtsOnObstacle;

    /*
    filter_tangent_line_equations_and_points_by_disks(tangentLineEqs, 
                                             filteredTangentLineEqs, 
                                             listOfPtsOnObstacle,
                                             betaUniverse, 
                                             pointVtx, obstacle);
      */                                     
    
    filter_tangent_line_equations_and_points_by_disks_vers1(tangentLineEqs,
                                             filteredTangentLineEqs, 
                                             listOfPtsOnObstacle,
                                             betaUniverse, 
                                             pointVtx, obstacle);
                                           
    vector<rg_ImplicitEquation> filteredTangentLineEqsByTopoloticalPath;
    vector<rg_Point2D> filteredListOfPtsOnObstacleByTopoloticalPath;

    filter_tangent_line_equations_and_points_by_topological_path(filteredTangentLineEqs, 
                                                                 listOfPtsOnObstacle, 
                                                                 filteredTangentLineEqsByTopoloticalPath,
                                                                 filteredListOfPtsOnObstacleByTopoloticalPath,
                                                                 pointVtx,
                                                                 betaUniverse,
                                                                 topologicalPath);
    //3. make verices and edges for graph
    make_vertices_N_tangent_line_edges(filteredTangentLineEqsByTopoloticalPath,
                                       filteredListOfPtsOnObstacleByTopoloticalPath,
                                       pointVtx,
                                       verticesOnObstacle,
                                       solutionPoolGraph);
}



void ShortestPathFinder2D::add_new_arc_edges_made_by_new_vertices_into_the_graph(VertexBU2D* obstacle, 
                                                                                 const vector<VertexForSolutionPoolGraph*>& newVerticesOnObstacle, 
                                                                                 const BetaComponentFinder& betaComponentFinder,
                                                                                 ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                                 unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    for (int i = 0; i < newVerticesOnObstacle.size(); i++)
    {
        VertexForSolutionPoolGraph* newVertexOnObstacle = newVerticesOnObstacle.at(i);

        if (does_this_obstacle_have_already_maden_vetices_of_tangent_points(obstacle, mapperForQTVtxToVtxOnDisk)) // 
        {
            const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex = mapperForQTVtxToVtxOnDisk.at(obstacle);

            list<pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>> vertexPairsForArcEdges;
            filter_vertices_for_making_arc_edges(obstacle, 
                                                 newVertexOnObstacle, 
                                                 verticesOnObstacleBeforeMakingNewVertex, 
                                                 betaComponentFinder, 
                                                 vertexPairsForArcEdges);

            make_arc_edges(vertexPairsForArcEdges, obstacle, solutionPoolGraph);
        }
        else
        {
            list<pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>> vertexPairsForArcEdges;

            filter_vertices_for_making_arc_edges_when_NONE_verices_on_obstacle(obstacle, 
                                                                               newVertexOnObstacle,
                                                                               betaComponentFinder, 
                                                                               vertexPairsForArcEdges);

            make_arc_edges(vertexPairsForArcEdges, obstacle, solutionPoolGraph);
        }
        
        insert_new_vertex_into_set_of_vertices_on_obstacle(newVertexOnObstacle, obstacle, mapperForQTVtxToVtxOnDisk);
    }
}



void ShortestPathFinder2D::insert_into_list_of_tangent_vertices_on_disk_by_angle(VertexBU2D* obstacle, VertexForSolutionPoolGraph * newVertexOnObstacle, list<VertexForSolutionPoolGraph*>& verticesOnObstacle) const
{
    list<VertexForSolutionPoolGraph*>::iterator it_CurrVertexOnObstacle = verticesOnObstacle.begin();

    while (it_CurrVertexOnObstacle != verticesOnObstacle.end())
    {
        VertexForSolutionPoolGraph* currVertex = *it_CurrVertexOnObstacle;

        if (get_angle_of_this_vertex_with_regard_to_this_obstacle(newVertexOnObstacle, obstacle) 
          < get_angle_of_this_vertex_with_regard_to_this_obstacle(currVertex, obstacle))
        {
            verticesOnObstacle.insert(it_CurrVertexOnObstacle, newVertexOnObstacle);
            break;
        }

        ++it_CurrVertexOnObstacle;
    }

    if (it_CurrVertexOnObstacle == verticesOnObstacle.end())
    {
        verticesOnObstacle.insert(it_CurrVertexOnObstacle, newVertexOnObstacle);
    }
}



bool ShortestPathFinder2D::are_these_two_obstacles_in_same_component(VertexBU2D* obstacle1, 
                                                                   VertexBU2D* obstacle2, 
                                                                   const BetaComponentFinder& betaComponentFinder) const
{
    BetaComponent2D* betaComponent1 = betaComponentFinder.findIncludedComponent(obstacle1);
    BetaComponent2D* betaComponent2 = betaComponentFinder.findIncludedComponent(obstacle2);

    if (betaComponent1 == betaComponent2)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void ShortestPathFinder2D::find_tangent_line_equations_of_two_obstacles_in_same_component(const rg_Circle2D& obstacleDisk1, const rg_Circle2D& obstacleDisk2, vector<rg_ImplicitEquation>& tangentLineEqs) const
{
    rg_ImplicitEquation tangentLineEq1(1);
    rg_ImplicitEquation tangentLineEq2(1);

    if (obstacleDisk1.isIncludedIn(obstacleDisk2) || obstacleDisk2.isIncludedIn(obstacleDisk1))
    {
        return;
    }

    make_exterior_tangent_lines_of_two_circles(obstacleDisk1, obstacleDisk2, tangentLineEq1, tangentLineEq2);

    tangentLineEqs.push_back(tangentLineEq1); 
    tangentLineEqs.push_back(tangentLineEq2);
}


void ShortestPathFinder2D::find_tangent_line_equations_of_two_obstacles_in_diff_component(const rg_Circle2D& obstacleDisk1, const rg_Circle2D& obstacleDisk2, vector<rg_ImplicitEquation>& tangentLineEqs) const
{
    rg_ImplicitEquation tangentLineEq1(1);
    rg_ImplicitEquation tangentLineEq2(1);
    rg_ImplicitEquation tangentLineEq3(1);
    rg_ImplicitEquation tangentLineEq4(1);

    if (obstacleDisk1.isIncludedIn(obstacleDisk2) || obstacleDisk2.isIncludedIn(obstacleDisk1))
    {
        return;
    }

    make_exterior_tangent_lines_of_two_circles(obstacleDisk1, obstacleDisk2, tangentLineEq1, tangentLineEq2);
    make_interior_tangent_lines_of_two_circles(obstacleDisk1, obstacleDisk2, tangentLineEq3, tangentLineEq4);

    tangentLineEqs.push_back(tangentLineEq1);
    tangentLineEqs.push_back(tangentLineEq2);
    tangentLineEqs.push_back(tangentLineEq3);
    tangentLineEqs.push_back(tangentLineEq4);
}




void ShortestPathFinder2D::filter_tangent_line_equations_and_points_by_disks(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                    vector<rg_ImplicitEquation>& filletedTangentLineEqs, 
                                                                    vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2, 
                                                                    BetaUniverse2D& betaUniverse, 
                                                                    VertexBU2D* obs1, VertexBU2D* obs2) const
{
    rg_Circle2D obstacleDisk1(obs1->getCoord(), obs1->getCircle().getRadius() + m_Probe.getRadius());
    rg_Circle2D obstacleDisk2(obs2->getCoord(), obs2->getCircle().getRadius() + m_Probe.getRadius());

    for (int i = 0; i < tangentLineEqs.size(); i++)
    {

        if (i == 2)
        {
            cout << endl;
        }

        rg_Point2D tangentPoint1 = compute_tangent_point_between_line_and_circle(obstacleDisk1, tangentLineEqs[i]);
        rg_Point2D tangentPoint2 = compute_tangent_point_between_line_and_circle(obstacleDisk2, tangentLineEqs[i]);

        unordered_set<VertexBU2D*> visitedObstacles;
        visitedObstacles.insert(obs1);
        visitedObstacles.insert(obs2);

        if (!is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(tangentLineEqs[i], 
                                                                                   tangentPoint1, tangentPoint2, 
                                                                                   betaUniverse, 
                                                                                   obs1, obs2,
                                                                                   visitedObstacles))
        {
            filletedTangentLineEqs.push_back(tangentLineEqs[i]);
            listOfTwoPtsOnObstacle1And2.push_back(make_pair(tangentPoint1, tangentPoint2));
        }
        /*
        if (!is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(tangentLineEqs[i], 
                                                                                   tangentPoint1, tangentPoint2, 
                                                                                   betaUniverse, 
                                                                                   visitedObstacles))
        {
            filletedTangentLineEqs.push_back(tangentLineEqs[i]);
            listOfTwoPtsOnObstacle1And2.push_back(make_pair(tangentPoint1, tangentPoint2));
        }
        */
    }
}

void ShortestPathFinder2D::filter_tangent_line_equations_and_points_by_disks_vers1(const vector<rg_ImplicitEquation>& tangentLineEqs, vector<rg_ImplicitEquation>& filletedTangentLineEqs, vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2, BetaUniverse2D & betaUniverse, VertexBU2D * obs1, VertexBU2D * obs2) const
{
    list<VertexBU2D*> allObstaclesWithoutTwoObstacles;

    rg_dList<VertexBU2D>& allObstacles = betaUniverse.getVertices();
    allObstacles.reset4Loop();
    while (allObstacles.setNext4Loop())
    {
        VertexBU2D* currObstacle = allObstacles.getpEntity();

        if (currObstacle == obs1)
        {
            continue;
        }
        else if (currObstacle == obs2)
        {
            continue;
        }

        allObstaclesWithoutTwoObstacles.push_back(currObstacle);
    }



    rg_Circle2D obstacleDisk1(obs1->getCoord(), obs1->getCircle().getRadius() + m_Probe.getRadius());
    rg_Circle2D obstacleDisk2(obs2->getCoord(), obs2->getCircle().getRadius() + m_Probe.getRadius());

    for (int i = 0; i < tangentLineEqs.size(); i++)
    {
        rg_Point2D tangentPoint1 = compute_tangent_point_between_line_and_circle(obstacleDisk1, tangentLineEqs[i]);
        rg_Point2D tangentPoint2 = compute_tangent_point_between_line_and_circle(obstacleDisk2, tangentLineEqs[i]);

        rg_Line2D tangentLineSementForTwoPts(tangentPoint1, tangentPoint2);
      

        if (!is_this_line_segement_intersected_with_one_of_obstacles(tangentLineSementForTwoPts, allObstaclesWithoutTwoObstacles))
        {
            filletedTangentLineEqs.push_back(tangentLineEqs[i]);
            listOfTwoPtsOnObstacle1And2.push_back(make_pair(tangentPoint1, tangentPoint2));
        }
    }
}


void ShortestPathFinder2D::filter_tangent_line_equations_and_points_by_topological_path(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                                        const vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2, 
                                                                                        vector<rg_ImplicitEquation>& filletedTangentLineEqs, 
                                                                                        vector<pair<rg_Point2D, rg_Point2D>>& filletedListOfTwoPtsOnObstacle1And2, 
                                                                                        BetaUniverse2D & betaUniverse, 
                                                                                        const list<EdgeBU2D*>& topologicalPath) const
{
    unordered_set<EdgeBU2D*> boundaryQTEdgesOfTP;

    //1. insert all qt-edges into the set
    for (list<EdgeBU2D*>::const_iterator it_QTEdge = topologicalPath.begin();
         it_QTEdge != topologicalPath.end();
         ++it_QTEdge)
    {
        EdgeBU2D* currQTEdge = *it_QTEdge;

        boundaryQTEdgesOfTP.insert(currQTEdge);
        boundaryQTEdgesOfTP.insert(currQTEdge->getRightHand());
        boundaryQTEdgesOfTP.insert(currQTEdge->getLeftHand());
        boundaryQTEdgesOfTP.insert(currQTEdge->getRightLeg());
        boundaryQTEdgesOfTP.insert(currQTEdge->getLeftLeg());
    }

    //2. erase all qt-edges of qt-path
    for (list<EdgeBU2D*>::const_iterator it_Edge = topologicalPath.begin(); 
         it_Edge != topologicalPath.end(); 
         ++it_Edge)
    {
        EdgeBU2D* currQTEdge = *it_Edge;
        boundaryQTEdgesOfTP.erase(currQTEdge);
    }

    //boundaryQTEdgesOfTP.erase(topologicalPath.begin(), topologicalPath.end());

    //3. check intersection btw the remaining qt-edges and tangent line
    for (int i = 0 ; i < listOfTwoPtsOnObstacle1And2.size(); ++i)
    {
        const pair<rg_Point2D, rg_Point2D>& twoPoints = listOfTwoPtsOnObstacle1And2.at(i);
        const rg_ImplicitEquation tangentLineEq       = tangentLineEqs.at(i);

        if (!are_tangent_line_segment_intersected_with_one_of_QT_edges(tangentLineEq, twoPoints.first, twoPoints.second, boundaryQTEdgesOfTP))
        {
            unordered_set<FaceBU2D*> qtFacesOfPath;
            convert_QT_edge_path_to_QT_faces(topologicalPath, qtFacesOfPath);

            if (do_any_of_QT_faces_include_this_pt(twoPoints.first, qtFacesOfPath, betaUniverse) 
             && do_any_of_QT_faces_include_this_pt(twoPoints.second, qtFacesOfPath, betaUniverse))
            {
                filletedTangentLineEqs.push_back(tangentLineEq);
                filletedListOfTwoPtsOnObstacle1And2.push_back(twoPoints);
            }
        }
        
    }
}



void ShortestPathFinder2D::filter_tangent_line_equations_and_points_by_topological_path(const vector<rg_ImplicitEquation>& tangentLineEqs, const vector<rg_Point2D>& listOfPtsOnObstacle, vector<rg_ImplicitEquation>& filletedTangentLineEqs, vector<rg_Point2D>& filletedListOfPtsOnObstacle, VertexForSolutionPoolGraph * pointVtx, BetaUniverse2D & betaUniverse, const list<EdgeBU2D*>& topologicalPath) const
{
     unordered_set<EdgeBU2D*> boundaryQTEdgesOfTP;

    //1. insert all qt-edges into the set
    for (list<EdgeBU2D*>::const_iterator it_QTEdge = topologicalPath.begin();
         it_QTEdge != topologicalPath.end();
         ++it_QTEdge)
    {
        EdgeBU2D* currQTEdge = *it_QTEdge;

        boundaryQTEdgesOfTP.insert(currQTEdge);
        boundaryQTEdgesOfTP.insert(currQTEdge->getRightHand());
        boundaryQTEdgesOfTP.insert(currQTEdge->getLeftHand());
        boundaryQTEdgesOfTP.insert(currQTEdge->getRightLeg());
        boundaryQTEdgesOfTP.insert(currQTEdge->getLeftLeg());
    }

    //2. erase all qt-edges of qt-path
    //boundaryQTEdgesOfTP.erase(topologicalPath.begin(), topologicalPath.end());
    for (list<EdgeBU2D*>::const_iterator it_Edge = topologicalPath.begin();
        it_Edge != topologicalPath.end();
        ++it_Edge)
    {
        EdgeBU2D* currQTEdge = *it_Edge;
        boundaryQTEdgesOfTP.erase(currQTEdge);
    }

    //3. check intersection btw the remaining qt-edges and tangent line
    for (int i = 0 ; i < listOfPtsOnObstacle.size(); ++i)
    {
        const rg_Point2D point = listOfPtsOnObstacle.at(i);
        const rg_ImplicitEquation tangentLineEq = tangentLineEqs.at(i);

        if (!are_tangent_line_segment_intersected_with_one_of_QT_edges(tangentLineEq, pointVtx->get_tangent_point_to_disk(), point, boundaryQTEdgesOfTP))
        {
            unordered_set<FaceBU2D*> qtFacesOfPath;
            convert_QT_edge_path_to_QT_faces(topologicalPath, qtFacesOfPath);

            if (do_any_of_QT_faces_include_this_pt(point, qtFacesOfPath, betaUniverse))
            {
                filletedTangentLineEqs.push_back(tangentLineEq);
                filletedListOfPtsOnObstacle.push_back(point);
            }
        }
    }
}


void ShortestPathFinder2D::filter_tangent_line_equations_and_points_by_disks(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                    vector<rg_ImplicitEquation>& filletedTangentLineEqs, 
                                                                    vector<rg_Point2D>& listOfPtsOnObstacle,
                                                                    BetaUniverse2D& betaUniverse, 
                                                                    VertexForSolutionPoolGraph* pointVtx, 
                                                                    VertexBU2D* obstacle) const
{
    rg_Circle2D obstacleDisk(obstacle->getCoord(), obstacle->getCircle().getRadius() + m_Probe.getRadius());

    for (int i = 0; i < tangentLineEqs.size(); i++)
    {
        rg_Point2D tangentPointOnObstacle = compute_tangent_point_between_line_and_circle(obstacleDisk, tangentLineEqs[i]);

        unordered_set<VertexBU2D*> visitedObstacles;
        visitedObstacles.insert(obstacle);

        if (!is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(tangentLineEqs[i], 
                                                                                   pointVtx->get_tangent_point_to_disk(), 
                                                                                   tangentPointOnObstacle, 
                                                                                   betaUniverse,
                                                                                   obstacle,
                                                                                   visitedObstacles))
        {
            filletedTangentLineEqs.push_back(tangentLineEqs[i]);
            listOfPtsOnObstacle.push_back(tangentPointOnObstacle);
        }

        /*
        if (!is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(tangentLineEqs[i], 
                                                                                   pointVtx->get_tangent_point_to_disk(), 
                                                                                   tangentPointOnObstacle, 
                                                                                   betaUniverse, 
                                                                                   visitedObstacles))
        {
            filletedTangentLineEqs.push_back(tangentLineEqs[i]);
            listOfPtsOnObstacle.push_back(tangentPointOnObstacle);
        }
        */
    }
}

void ShortestPathFinder2D::filter_tangent_line_equations_and_points_by_disks_vers1(const vector<rg_ImplicitEquation>& tangentLineEqs, vector<rg_ImplicitEquation>& filletedTangentLineEqs, vector<rg_Point2D>& listOfPtsOnObstacle, BetaUniverse2D & betaUniverse, VertexForSolutionPoolGraph * pointVtx, VertexBU2D * obstacle) const
{
    rg_Circle2D obstacleDisk(obstacle->getCoord(), obstacle->getCircle().getRadius() + m_Probe.getRadius());

    for (int i = 0; i < tangentLineEqs.size(); i++)
    {
        rg_Point2D tangentPointOnObstacle = compute_tangent_point_between_line_and_circle(obstacleDisk, tangentLineEqs[i]);

        rg_Line2D tangentLineSementForTwoPts(pointVtx->get_tangent_point_to_disk(), tangentPointOnObstacle);
        
        list<VertexBU2D*> allObstaclesWithoutTwoObstacles;

        rg_dList<VertexBU2D>& allObstacles = betaUniverse.getVertices();
        allObstacles.reset4Loop();
        while (allObstacles.setNext4Loop())
        {
            VertexBU2D* currObstacle = allObstacles.getpEntity();

            if (currObstacle == obstacle)
            {
                continue;
            }

            allObstaclesWithoutTwoObstacles.push_back(currObstacle);
        }

        if (!is_this_line_segement_intersected_with_one_of_obstacles(tangentLineSementForTwoPts, allObstaclesWithoutTwoObstacles))
        {
            filletedTangentLineEqs.push_back(tangentLineEqs[i]);
            listOfPtsOnObstacle.push_back(tangentPointOnObstacle);
        }
    }
}



bool ShortestPathFinder2D::is_this_line_segement_intersected_with_one_of_obstacles(const rg_Line2D & lineSegment, const list<VertexBU2D*>& obstacles) const
{
    bool isThisLineSegmentIntersectedWithOneOfObstacles = false;

    for(list<VertexBU2D*>::const_iterator it_Obstacle = obstacles.begin(); it_Obstacle != obstacles.end(); ++it_Obstacle)
    {
        VertexBU2D* currObstacle = *it_Obstacle;

        if (currObstacle->isVirtual())
        {
            continue;
        }

        rg_Circle2D currObstacleCircle(currObstacle->getCircle());

        if (currObstacleCircle.hasIntersectionWith(lineSegment))
        {
            isThisLineSegmentIntersectedWithOneOfObstacles = true;
            break;
        }
    }

    return isThisLineSegmentIntersectedWithOneOfObstacles;
}

bool ShortestPathFinder2D::is_this_line_segement_intersected_with_one_of_ALL_obstacles(const rg_Line2D& lineSegment)
{
    list<VertexBU2D*> allObstacleList;

    rg_dList<VertexBU2D>& allObstacles = m_BetaUniverse.getVertices();

    allObstacles.reset4Loop();
    while (allObstacles.setNext4Loop())
    {
        VertexBU2D* currObstacle = allObstacles.getpEntity();
        allObstacleList.push_back(currObstacle);
    }

    return is_this_line_segement_intersected_with_one_of_obstacles(lineSegment, allObstacleList);
}

bool ShortestPathFinder2D::is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                                                 const rg_Point2D& startPtOfTangentLine, 
                                                                                                 const rg_Point2D& endPtOfTangentLine, 
                                                                                                 BetaUniverse2D& betaUniverse, 
                                                                                                 FaceBU2D* startQTFace, FaceBU2D* endQTFace, 
                                                                                                 unordered_set<VertexBU2D*>& visitedObstacles) const
{
     bool isThisLineIntersectedWithAnotherObstacle = false;
     
     FaceBU2D* currQTFace = startQTFace;
     FaceBU2D* prevQTFace = NULL;
     
     if (currQTFace == endQTFace)
     {
        
         isThisLineIntersectedWithAnotherObstacle = are_this_line_intersected_with_one_of_three_disks_of_QT_face(lineEq,
                                                                                                                 currQTFace,
                                                                                                                 startPtOfTangentLine,
                                                                                                                 endPtOfTangentLine,
                                                                                                                 visitedObstacles);
     }
     else
     {
         while (currQTFace != endQTFace)
         {
             //0. SPECIAL CASE - RED N GREEN QT-FACE // CHECK BY FULL DISKS
             if ((currQTFace->computeSignedArea() < 0) || currQTFace->isAnomalyInBetaComplex(m_Probe.getRadius()))
             {
                 isThisLineIntersectedWithAnotherObstacle = check_intersection_between_line_segment_N_all_disks(lineEq,
                                                                                                                startPtOfTangentLine,
                                                                                                                endPtOfTangentLine,
                                                                                                                betaUniverse,
                                                                                                                visitedObstacles);
                 break;
             }


             //1. check intersection with three disks 
             if (are_this_line_intersected_with_one_of_three_disks_of_QT_face(lineEq,
                                                                              currQTFace,
                                                                              startPtOfTangentLine,
                                                                              endPtOfTangentLine,
                                                                              visitedObstacles))
             {
                 isThisLineIntersectedWithAnotherObstacle = true;
                 break;
             }

             //2. propagate to one neighbor face
             propagate_to_one_of_neighbor_QT_faces_passed_through_by_tangent_line(currQTFace, prevQTFace, startPtOfTangentLine, endPtOfTangentLine, lineEq);
         
         }

        if (currQTFace == endQTFace)
        {
        
         isThisLineIntersectedWithAnotherObstacle = are_this_line_intersected_with_one_of_three_disks_of_QT_face(lineEq,
                                                                                                                 currQTFace,
                                                                                                                 startPtOfTangentLine,
                                                                                                                 endPtOfTangentLine,
                                                                                                                 visitedObstacles);
        }
     }

     return isThisLineIntersectedWithAnotherObstacle;
}




void ShortestPathFinder2D::expand_candidate_QT_edges_by_including_on_a_barrier(const unordered_map<EdgeBU2D*, bool>& QTEdgeValidation, unordered_map<EdgeBU2D*, bool>& expandedQTEdgeValidation, BetaUniverse2D & betaUniverse)
{
    //1. find all component in ellipse
    unordered_set<BetaComponent2D*> betaComponentSet;

    for (unordered_map<EdgeBU2D*, bool>::const_iterator it_Pair = QTEdgeValidation.begin();
        it_Pair != QTEdgeValidation.end();
        it_Pair++)
    {
        const pair<EdgeBU2D*, bool>& pair = *it_Pair;

        EdgeBU2D* currQTEdge  = pair.first;
        bool      isInEllipse = pair.second;

        if (isInEllipse)
        {
            VertexBU2D* currQTStartVertex = currQTEdge->getStartVertex();
            VertexBU2D* currQTEndVertex   = currQTEdge->getEndVertex();

            BetaComponent2D* foundComponent = NULL;

            foundComponent = m_BetaComponentFinder.findIncludedComponent(currQTStartVertex);
            
            if (foundComponent != NULL)
            {
                betaComponentSet.insert(foundComponent);
            }

            foundComponent = m_BetaComponentFinder.findIncludedComponent(currQTEndVertex);

            if (foundComponent != NULL)
            {
                betaComponentSet.insert(foundComponent);
            }
        }
    }

    //2. add qt-edges which have an vertices of the component
    const rg_dList<EdgeBU2D>& allQTEdges = betaUniverse.getEdges();

    allQTEdges.reset4Loop();
    while (allQTEdges.setNext4Loop())
    {
        EdgeBU2D* currQTEdge = allQTEdges.getpEntity();

        bool isInPrevEllipse = QTEdgeValidation.at(currQTEdge);

        if (isInPrevEllipse)
        {
            expandedQTEdgeValidation.insert(make_pair(currQTEdge, true));
            continue;
        }
        
        if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            VertexBU2D* currQTStartVertex = currQTEdge->getStartVertex();
            VertexBU2D* currQTEndVertex   = currQTEdge->getEndVertex();

            BetaComponent2D* foundComponent = NULL;

            foundComponent = m_BetaComponentFinder.findIncludedComponent(currQTStartVertex);

            if (betaComponentSet.find(foundComponent) != betaComponentSet.end())
            {
                expandedQTEdgeValidation.insert(make_pair(currQTEdge, true));
                continue;
            }

            foundComponent = m_BetaComponentFinder.findIncludedComponent(currQTEndVertex);

            if (betaComponentSet.find(foundComponent) != betaComponentSet.end())
            {
                expandedQTEdgeValidation.insert(make_pair(currQTEdge, true));
                continue;
            }
        }
       
        expandedQTEdgeValidation.insert(make_pair(currQTEdge, false));
    }
}


void ShortestPathFinder2D::expand_candidate_QT_vertices_by_including_on_a_barrier(const unordered_map<VertexBU2D*, bool>& QTVertexValidation, unordered_map<VertexBU2D*, bool>& expandedQTVertexValidation, BetaUniverse2D & betaUniverse)
{
    //1. find all component in ellipse
    unordered_set<BetaComponent2D*> betaComponentSet;

    for (unordered_map<VertexBU2D*, bool>::const_iterator it_Pair = QTVertexValidation.begin();
        it_Pair != QTVertexValidation.end();
        it_Pair++)
    {
        const pair<VertexBU2D*, bool>& pair = *it_Pair;

        VertexBU2D* currQTVertex = pair.first;
        bool       isInEllipse   = pair.second;

        if (isInEllipse)
        {
            BetaComponent2D* foundComponent = NULL;

            foundComponent = m_BetaComponentFinder.findIncludedComponent(currQTVertex);

            if (foundComponent != NULL)
            {
                betaComponentSet.insert(foundComponent);
            }
        }
    }

    //2. add qt-edges which have an vertices of the component
    const rg_dList<VertexBU2D>& allQTVertices = betaUniverse.getVertices();

    allQTVertices.reset4Loop();
    while (allQTVertices.setNext4Loop())
    {
        VertexBU2D* currQTVertex = allQTVertices.getpEntity();

        bool isInPrevEllipse = QTVertexValidation.at(currQTVertex);

        if (isInPrevEllipse)
        {
            expandedQTVertexValidation.insert(make_pair(currQTVertex, true));
            continue;
        }

        BetaComponent2D* foundComponent = NULL;

        foundComponent = m_BetaComponentFinder.findIncludedComponent(currQTVertex);

        if (betaComponentSet.find(foundComponent) != betaComponentSet.end())
        {
            expandedQTVertexValidation.insert(make_pair(currQTVertex, true));
            continue;
        }

        expandedQTVertexValidation.insert(make_pair(currQTVertex, false));
    }
}




bool ShortestPathFinder2D::are_this_line_intersected_with_these_disks(const rg_ImplicitEquation& tangentLineEq, 
                                                                     rg_dList<VertexBU2D*>& obstacles, 
                                                                     const rg_Point2D& startPt, 
                                                                     const rg_Point2D& endPt,
                                                                     unordered_set<VertexBU2D*>& visitedObstacles) const
{
    bool is_this_tangent_line_intersected_with_one_disk = false;

    obstacles.reset4Loop();
    while (obstacles.setNext4Loop())
    {
        VertexBU2D* currQTVertex = obstacles.getEntity();

        rg_dList<VertexBU2D*> qtCurrAndNeighborVertices;
        currQTVertex->getNeighborVertices(qtCurrAndNeighborVertices);

        qtCurrAndNeighborVertices.add(currQTVertex);

        qtCurrAndNeighborVertices.reset4Loop();
        while (qtCurrAndNeighborVertices.setNext4Loop())
        {
            VertexBU2D* currQTVertexOfNeighbors = qtCurrAndNeighborVertices.getEntity();

            if (currQTVertexOfNeighbors->isVirtual())
            {
                continue;
            }

            // if already check the intersection btw disk and a line
            if (visitedObstacles.find(currQTVertexOfNeighbors) != visitedObstacles.end())
            {
                continue;
            }

            //0. visit obstacle 
            visitedObstacles.insert(currQTVertexOfNeighbors);

            //1 .check distance vs radius 
            rg_Circle2D currObstacle(currQTVertexOfNeighbors->getCoord(), currQTVertexOfNeighbors->getCircle().getRadius() + get_probe().getRadius());

            rg_Line<rg_Point2D> line(startPt, endPt);
            double distanceFromCircleCenterToTangentLine = line.getDistance(currObstacle.getCenterPt());

            if (!rg_LT(distanceFromCircleCenterToTangentLine, currObstacle.getRadius(),resNeg3))
            {
                continue;
            }

            //2. compute intersection point and check intersection point in the interval

            // Circie center : (a, b), radius : r
            // Line : ex + fy + g = 0
            // 
            // 1. if f = 0, x = - (g / e), y = b +- sqrt(r^2 - (g/e + a)^2)
            // 2. if e = 0, x = a +- sqrt(r^2 - (g/f + b)^2), y = - (g / f)
            // 3. else, x = (af^2 - eA +- sqrt(r^2B - (A + ea)^2)) / B, y = - (e / f)x - (g / f)
            // ref. A = bf + g, B = f^2 + e^2
            double a = currObstacle.getCenterPt().getX();
            double b = currObstacle.getCenterPt().getY();
            double r = currObstacle.getRadius();
            double e = tangentLineEq.getCoeff(1, 0);
            double f = tangentLineEq.getCoeff(0, 1); 
            double g = tangentLineEq.getCoeff(0, 0);

            double x1, x2, y1, y2;

            if (rg_ZERO(f))
            {
                // 1. if f = 0, x = - (g / e), y = b +- sqrt(r^2 - (g/e + a)^2)

                double sqrtVal = sqrt(pow(r, 2) - pow(g / e + a, 2));

                y1 = b + sqrtVal;
                y2 = b - sqrtVal;

                is_this_tangent_line_intersected_with_one_disk = is_this_val_in_this_interval(y1, startPt.getY(), endPt.getY())
                                                                || is_this_val_in_this_interval(y2, startPt.getY(), endPt.getY())
                                                                || is_this_val_in_this_interval(startPt.getY(), y1, y2)
                                                                || is_this_val_in_this_interval(endPt.getY(), y1, y2);
            }
            else if (rg_ZERO(e))
            {
                // 2. if e = 0, x = a +- sqrt(r^2 - (g/f + b)^2), y = - (g / f)

                double sqrtVal = sqrt(pow(r, 2) - pow(g / f + b, 2));

                x1 = a + sqrtVal;
                x2 = a - sqrtVal;

                is_this_tangent_line_intersected_with_one_disk = is_this_val_in_this_interval(x1, startPt.getX(), endPt.getX())
                                                                || is_this_val_in_this_interval(x2, startPt.getX(), endPt.getX())
                                                                || is_this_val_in_this_interval(startPt.getX(), x1, x2)
                                                                || is_this_val_in_this_interval(endPt.getX(), x1, x2);
            }
            else // choose x or y value. It works.
            {
                // 3. else, x = (af^2 - eA +- f * sqrt(r^2B - (A + ea)^2)) / B, y = - (e / f)x - (g / f)
                // ref. A = bf + g, B = f^2 + e^2

                double A = b*f + g;
                double B = pow(f, 2) + pow(e, 2);
                double sqrtVal = sqrt(pow(r, 2)*B - pow(A + e*a, 2));

                x1 = (a*pow(f, 2) - e*A + f * sqrtVal) / B;
                x2 = (a*pow(f, 2) - e*A - f * sqrtVal) / B;

                is_this_tangent_line_intersected_with_one_disk = is_this_val_in_this_interval(x1, startPt.getX(), endPt.getX())  
                                                                || is_this_val_in_this_interval(x2, startPt.getX(), endPt.getX()) 
                                                                || is_this_val_in_this_interval(startPt.getX(), x1, x2)
                                                                || is_this_val_in_this_interval(endPt.getX(), x1, x2);
            }

            if (is_this_tangent_line_intersected_with_one_disk)
            {
                break;
            }
        }

        if (is_this_tangent_line_intersected_with_one_disk)
        {
            break;
        }
    }

    return is_this_tangent_line_intersected_with_one_disk;
}



bool ShortestPathFinder2D::check_intersection_between_line_segment_N_all_disks(const rg_ImplicitEquation& lineEq, 
                                                                               const rg_Point2D& startPtOfTangentLine, 
                                                                               const rg_Point2D& endPtOfTangentLine, 
                                                                               BetaUniverse2D& betaUniverse, 
                                                                               unordered_set<VertexBU2D*>& visitedObstacles) const
{
    rg_dList<VertexBU2D*> verticesPointers;
    rg_dList<VertexBU2D>& vertices = betaUniverse.getVertices();

    vertices.reset4Loop();
    while (vertices.setNext4Loop())
    {
        VertexBU2D* currVertex = vertices.getpEntity();
        verticesPointers.add(currVertex);
    }

    return are_this_line_intersected_with_these_disks(lineEq, verticesPointers, startPtOfTangentLine, endPtOfTangentLine, visitedObstacles);
}




void ShortestPathFinder2D::propagate_to_one_of_neighbor_QT_faces_passed_through_by_tangent_line(FaceBU2D*& currQTFace,
                                                                                                FaceBU2D*& prevQTFace,
                                                                                                const rg_Point2D& startPtOfTangentLine,
                                                                                                const rg_Point2D& endPtOfTangentLine,
                                                                                                const rg_ImplicitEquation& tangentLineEq) const
{
    switch (currQTFace->isVirtual())
    {
    case true:
    {
        //1. get a general edge
        EdgeBU2D* generalQTEdge = NULL;
        find_one_general_QT_edge_bounding_QT_Face(generalQTEdge, currQTFace);

        //2. check whether tangent line go into inner qt
        if (are_tangent_line_segment_N_QT_edge_intersected(tangentLineEq, startPtOfTangentLine, endPtOfTangentLine, generalQTEdge))
        {
            set_prev_and_curr_QT_face_for_loop(currQTFace, prevQTFace, generalQTEdge);
        }
        else
        {
            set_prev_and_curr_virtual_QT_face_for_loop(currQTFace, prevQTFace, generalQTEdge, endPtOfTangentLine);
        }

        break;
    }

    case false:
    {
        rg_dList<EdgeBU2D*> boundingQTEdges;
        currQTFace->getBoundingEdges(boundingQTEdges);

        boundingQTEdges.reset4Loop();
        while (boundingQTEdges.setNext4Loop())
        {
            EdgeBU2D* currQTEdge = boundingQTEdges.getEntity();

            if ((currQTEdge->getRightFace() == prevQTFace) || (currQTEdge->getLeftFace() == prevQTFace))
            {
                continue;
            }

            if (are_tangent_line_segment_N_QT_edge_intersected(tangentLineEq, startPtOfTangentLine, endPtOfTangentLine, currQTEdge))
            {
                set_prev_and_curr_QT_face_for_loop(currQTFace, prevQTFace, currQTEdge);
                break;
            }
        }

        break;
    }
    }
}


void ShortestPathFinder2D::find_one_general_QT_edge_bounding_QT_Face(EdgeBU2D*& generalQTEdge, FaceBU2D* qtFace) const
{
    rg_dList<EdgeBU2D*> boundingQTEdges;
    qtFace->getBoundingEdges(boundingQTEdges);

    boundingQTEdges.reset4Loop();
    while (boundingQTEdges.setNext4Loop())
    {
        EdgeBU2D* currBoundingQTEdge = boundingQTEdges.getEntity();

        if (!currBoundingQTEdge->isVirtual())
        {
            generalQTEdge = currBoundingQTEdge;
            break;
        }
    }
}


void ShortestPathFinder2D::set_prev_and_curr_virtual_QT_face_for_loop(FaceBU2D*& currQTFace,
                                                                      FaceBU2D*& prevQTFace, 
                                                                      EdgeBU2D* currQTEdge, 
                                                                      const rg_Point2D& endPt) const
{
    rg_dList<FaceBU2D*> incidentQTFaces;

    rg_Point2D obstacleCenterOfStartVtx = currQTEdge->getStartVertex()->getCoord();
    rg_Point2D obstacleCenterOfEndVtx   = currQTEdge->getEndVertex()->getCoord();

    if (endPt.distance(obstacleCenterOfStartVtx) <= endPt.distance(obstacleCenterOfEndVtx))
    {

        currQTEdge->getStartVertex()->getIncidentFaces(incidentQTFaces);
    }
    else
    {
        currQTEdge->getEndVertex()->getIncidentFaces(incidentQTFaces);
    }


    incidentQTFaces.reset4Loop();
    while (incidentQTFaces.setNext4Loop())
    {
        FaceBU2D* currQTFaceOfLoop = incidentQTFaces.getEntity();

        if (currQTFaceOfLoop == currQTFace) // always two virtual face at the boundary vtx
        {
            continue;
        }

        if (currQTFaceOfLoop->isVirtual())
        {
            prevQTFace = currQTFace;
            currQTFace = currQTFaceOfLoop;
            break;
        }
    }
}


bool ShortestPathFinder2D::are_tangent_line_segment_N_QT_edge_intersected(const rg_ImplicitEquation& tangentLineEq,
                                                                          const rg_Point2D& startPt,
                                                                          const rg_Point2D& endPt,
                                                                          EdgeBU2D* qtEdge) const
{
    rg_Point2D qtPoint1;
    rg_Point2D qtPoint2;

    if (qtEdge->getStartVertex()->isVirtual())
    {
        qtPoint1 = qtEdge->computeASamplePointOnThisQTEdge();
    }
    else
    {
        qtPoint1 = qtEdge->getStartVertex()->getCoord();
    }


    if (qtEdge->getEndVertex()->isVirtual())
    {
        qtPoint2 = qtEdge->computeASamplePointOnThisQTEdge();
    }
    else
    {
        qtPoint2 = qtEdge->getEndVertex()->getCoord();
    }

    rg_ImplicitEquation qtEdgeLineEq = make_line_equation(qtPoint1, qtPoint2);

    // Line1. a1x + b1y + c1 = 0
    // Line2. a2x + b2y + c2 = 0
    // x = (b1c2 - b2c1) / (A),
    // y = (a2c1 - a1c2) / (A),
    // A = a1b2 - a2b1

    double a1 = tangentLineEq.getCoeff(1, 0);
    double b1 = tangentLineEq.getCoeff(0, 1);
    double c1 = tangentLineEq.getCoeff(0, 0);

    double a2 = qtEdgeLineEq.getCoeff(1, 0);
    double b2 = qtEdgeLineEq.getCoeff(0, 1);
    double c2 = qtEdgeLineEq.getCoeff(0, 0);

    double A = a1*b2 - a2*b1;

    if (rg_ZERO(A))
    {
        return false;
    }

    double x = (b1*c2 - b2*c1) / A;
    double y = (a2*c1 - a1*c2) / A;

    double x_qt_start = qtPoint1.getX();
    double x_qt_end = qtPoint2.getX();

    double y_qt_start = qtPoint1.getY();
    double y_qt_end = qtPoint2.getY();


    bool is_this_tangent_line_intersected_with_qt_line = false;

    if (rg_EQ(x_qt_start, x_qt_end, resNeg3) || rg_EQ(startPt.getX(), endPt.getX(), resNeg3))
    {

        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(y, startPt.getY(), endPt.getY())
                                                     && is_this_val_in_this_interval(y, y_qt_start, y_qt_end);
    }
    else if (rg_EQ(y_qt_start, y_qt_end, resNeg3) || rg_EQ(startPt.getY(), endPt.getY(), resNeg3))
    {
        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(x, startPt.getX(), endPt.getX())
                                                     && is_this_val_in_this_interval(x, x_qt_start, x_qt_end);
    }
    else
    {
        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(x, startPt.getX(), endPt.getX())
                                                     && is_this_val_in_this_interval(x, x_qt_start, x_qt_end);
    }


    /*
    if (rg_EQ(x_qt_start, x_qt_end, resNeg3) && rg_EQ(startPt.getX(), endPt.getX(), resNeg3))
    {
        double y = (a2*c1 - a1*c2) / A;
        double y_qt_start = qtPoint1.getY();
        double y_qt_end   = qtPoint2.getY();

        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(y, startPt.getY(), endPt.getY())
                                                     && is_this_val_in_this_interval(y, y_qt_start, y_qt_end);
    }
    else if (rg_EQ(x_qt_start, x_qt_end, resNeg3) && !rg_EQ(startPt.getX(), endPt.getX(), resNeg3))
    {
        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(x, startPt.getX(), endPt.getX());
    }
    else if (!rg_EQ(x_qt_start, x_qt_end, resNeg3) && rg_EQ(startPt.getX(), endPt.getX(), resNeg3))
    {
        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(x, x_qt_start, x_qt_end);
    }
    else
    {
        is_this_tangent_line_intersected_with_qt_line = is_this_val_in_this_interval(x, startPt.getX(), endPt.getX()) 
                                                     && is_this_val_in_this_interval(x, x_qt_start, x_qt_end);
    }
    */

    return is_this_tangent_line_intersected_with_qt_line;
}


bool ShortestPathFinder2D::are_tangent_line_segment_intersected_with_one_of_QT_edges(const rg_ImplicitEquation & tangentLineEq, const rg_Point2D & startPt, const rg_Point2D & endPt, const unordered_set<EdgeBU2D*>& qtEdges) const
{
    bool areTheseTwoLineIntersected = false;

    for (unordered_set<EdgeBU2D*>::const_iterator it_QTEdge = qtEdges.begin();
        it_QTEdge != qtEdges.end();
        ++it_QTEdge)
    {
        EdgeBU2D* currQTEdge = *it_QTEdge;

        if (are_tangent_line_segment_N_QT_edge_intersected(tangentLineEq, startPt, endPt, currQTEdge))
        {
            areTheseTwoLineIntersected = true;
            break;
        }
    }

    return areTheseTwoLineIntersected;
}



bool ShortestPathFinder2D::are_tangent_line_segment_intersected_with_one_of_QT_edges(const rg_ImplicitEquation & tangentLineEq, const rg_Point2D & startPt, const rg_Point2D & endPt, const list<EdgeBU2D*>& qtEdges) const
{
    bool areTheseTwoLineIntersected = false;

    for (list<EdgeBU2D*>::const_iterator it_QTEdge = qtEdges.begin();
        it_QTEdge != qtEdges.end();
        ++it_QTEdge)
    {
        EdgeBU2D* currQTEdge = *it_QTEdge;

        if (are_tangent_line_segment_N_QT_edge_intersected(tangentLineEq, startPt, endPt, currQTEdge))
        {
            areTheseTwoLineIntersected = true;
            break;
        }
    }

    return areTheseTwoLineIntersected;
}




void ShortestPathFinder2D::make_vertices_N_tangent_line_edges(const vector<rg_ImplicitEquation>& filteredTangentLineEqs, 
                                                              const vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2,
                                                              vector<VertexForSolutionPoolGraph*>& verticesOnObstacle1,
                                                              vector<VertexForSolutionPoolGraph*>& verticesOnObstacle2,
                                                              ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    //0. Check whether there are tangent lines
    if (filteredTangentLineEqs.size() == 0)
    {
        return;
    }

    // REPEAT the followings
    for (int i = 0; i < listOfTwoPtsOnObstacle1And2.size(); i++)
    {
        //1. Make two vertices and a tangent line edge
        VertexForSolutionPoolGraph* vtxOnObstacle1OnMaking = NULL; 
        VertexForSolutionPoolGraph* vtxOnObstacle2OnMaking = NULL;
        EdgeForSolutionPoolGraph* tangentLineEdge = NULL;

        create_two_vertices_and_a_tangent_line_edge(vtxOnObstacle1OnMaking, vtxOnObstacle2OnMaking, tangentLineEdge, solutionPoolGraph);

        //2. Set two vertices aqnd a tangent line edge
        set_two_vertices_and_a_tangent_line_edge(vtxOnObstacle1OnMaking, 
                                                 vtxOnObstacle2OnMaking, 
                                                 tangentLineEdge, 
                                                 listOfTwoPtsOnObstacle1And2.at(i),
                                                 filteredTangentLineEqs.at(i));

        //3. Insert into container..
        verticesOnObstacle1.push_back(vtxOnObstacle1OnMaking);
        verticesOnObstacle2.push_back(vtxOnObstacle2OnMaking);
    }
}


void ShortestPathFinder2D::make_vertices_N_tangent_line_edges(const vector<rg_ImplicitEquation>& filteredTangentLineEqs,
                                                              const vector<rg_Point2D>& listOfPtsOnObstacle, 
                                                              VertexForSolutionPoolGraph* pointVtx,
                                                              vector<VertexForSolutionPoolGraph*>& verticesOnObstacle,
                                                              ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
     //0. Check whether there are tangent lines
    if (filteredTangentLineEqs.size() == 0)
    {
        return;
    }

    // REPEAT the followings
    for (int i = 0; i < listOfPtsOnObstacle.size(); i++)
    {
        //1. Make a vertex and a tangent line edge
        VertexForSolutionPoolGraph* vtxOnObstacleOnMaking = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
        EdgeForSolutionPoolGraph* tangentLineEdge         = solutionPoolGraph.create_edge(EdgeForSolutionPoolGraph());
        
        //2. Set two vertices aqnd a tangent line edge
        set_two_vertices_and_a_tangent_line_edge(pointVtx,
                                                 vtxOnObstacleOnMaking, 
                                                 tangentLineEdge, 
                                                 pointVtx->get_tangent_point_to_disk(),
                                                 listOfPtsOnObstacle.at(i),
                                                 filteredTangentLineEqs.at(i));

        //3. Insert into container..
        verticesOnObstacle.push_back(vtxOnObstacleOnMaking);
    }
}




void ShortestPathFinder2D::set_two_vertices_and_a_tangent_line_edge(VertexForSolutionPoolGraph* vtxOfPoint, 
                                                                    VertexForSolutionPoolGraph* vtxOnObstacle, 
                                                                    EdgeForSolutionPoolGraph* edge, 
                                                                    const rg_Point2D& point,
                                                                    const rg_Point2D& ptOnObstacle,
                                                                    const rg_ImplicitEquation& implicitEquationOfTangentLine) const
{
    //2. Set two vertices
    vtxOfPoint->add_line_edge(edge);
    
    vtxOnObstacle->set_tangent_point_to_disk(ptOnObstacle);
    vtxOnObstacle->add_line_edge(edge);
    vtxOnObstacle->set_accumulate_length_from_source_vtx(DBL_MAX);

    //2. Set a tangent line edge
    double tangentLineSementLength = point.distance(ptOnObstacle);

    edge->set_start_vtx(vtxOfPoint);
    edge->set_end_vtx(vtxOnObstacle);
    edge->set_edge_type(LINE_SEGMENT_SPG);
    edge->set_edge_length(tangentLineSementLength);
    edge->set_implicit_equation(implicitEquationOfTangentLine);
}



EdgeForSolutionPoolGraph * ShortestPathFinder2D::find_arc_edge_having_these_two_vertices(VertexForSolutionPoolGraph * vertex1, VertexForSolutionPoolGraph * vertex2) const
{
    if (vertex1 == NULL)
    {
        return NULL;
    }
    else if (vertex2 == NULL)
    {
        return NULL;
    }

    EdgeForSolutionPoolGraph * targetArcEdge = NULL;

    list<EdgeForSolutionPoolGraph*> arcEdges;
    vertex1->get_arc_edges(arcEdges);

    for (list<EdgeForSolutionPoolGraph*>::const_iterator it_ArcEdge = arcEdges.begin();
         it_ArcEdge != arcEdges.end();
         ++it_ArcEdge)
    {
        EdgeForSolutionPoolGraph* currArcEdge = *it_ArcEdge;

        if ((currArcEdge->get_start_vtx() == vertex1) && (currArcEdge->get_end_vtx() == vertex2))
        {
            targetArcEdge = currArcEdge;
            break;
        }
        else if ((currArcEdge->get_start_vtx() == vertex2) && (currArcEdge->get_end_vtx() == vertex1))
        {
            targetArcEdge = currArcEdge;
            break;
        }
    }

    return targetArcEdge;
}

rg_ImplicitEquation ShortestPathFinder2D::make_implicit_equation_of_arc_edge(VertexBU2D* obstacle, const double& probeRadius) const
{
    // Circle center : (a, b) , Radii : r
    // x^2 -2ax + y^2 -2by + (a^2 + b^2 - r^2) = 0

    rg_Point2D obstacleCenter = obstacle->getCircle().getCenterPt();

    double a = obstacleCenter.getX();
    double b = obstacleCenter.getY();
    double r = obstacle->getCircle().getRadius() + probeRadius;

    rg_ImplicitEquation arcImplicitEq(2);
    arcImplicitEq.setCoeff(2, 0, 1.0);
    arcImplicitEq.setCoeff(1, 0, -2.0 * a);
    arcImplicitEq.setCoeff(0, 2, 1.0);
    arcImplicitEq.setCoeff(0, 1, -2.0 * b);
    arcImplicitEq.setCoeff(0, 0, pow(a,2) + pow(b,2) - pow(r,2));

    return arcImplicitEq;
}


rg_ImplicitEquation ShortestPathFinder2D::make_implicit_equation_of_arc_edge(const rg_Circle2D& obstacle, const double& probeRadius) const
{
    rg_Point2D obstacleCenter = obstacle.getCenterPt();

    double a = obstacleCenter.getX();
    double b = obstacleCenter.getY();
    double r = obstacle.getRadius() + probeRadius;

    rg_ImplicitEquation arcImplicitEq(2);
    arcImplicitEq.setCoeff(2, 0, 1.0);
    arcImplicitEq.setCoeff(1, 0, -2.0 * a);
    arcImplicitEq.setCoeff(0, 2, 1.0);
    arcImplicitEq.setCoeff(0, 1, -2.0 * b);
    arcImplicitEq.setCoeff(0, 0, pow(a, 2) + pow(b, 2) - pow(r, 2));

    return arcImplicitEq;
}


double ShortestPathFinder2D::compute_arc_segment_length_on_obstacle_boundary(VertexBU2D* obstacle, 
                                                                             const rg_Point2D& arcStartPt, 
                                                                             const rg_Point2D& arcEndPt) const
{
    if (arcStartPt == arcEndPt)
    {
        double     obstacleRadii = obstacle->getCircle().getRadius();
        
        return (2*rg_PI*obstacleRadii);
    }
    else
    {
        rg_Point2D obstacleCenter = obstacle->getCircle().getCenterPt();
        double     obstacleRadii = obstacle->getCircle().getRadius();

        rg_Point2D vec1 = arcStartPt - obstacleCenter;
        rg_Point2D vec2 = arcEndPt - obstacleCenter;

        double theta = angleFromVec1toVec2(vec1, vec2); //always less than PI

        return (theta*obstacleRadii);
    }
}


double ShortestPathFinder2D::compute_arc_segment_length_on_obstacle_boundary(const rg_Circle2D& obstacle,
                                                                            const rg_Point2D& arcStartPt,
                                                                            const rg_Point2D& arcEndPt) const
{/*
    if (arcStartPt == arcEndPt)
    {
        double     obstacleRadii = obstacle.getRadius();

        return 0.0;// (2 * rg_PI*obstacleRadii);
    }
    else
    {
        rg_Point2D obstacleCenter = obstacle.getCenterPt();
        double     obstacleRadii = obstacle.getRadius();

        rg_Point2D vec1 = arcStartPt - obstacleCenter;
        rg_Point2D vec2 = arcEndPt - obstacleCenter;

        double theta = angleFromVec1toVec2(vec1, vec2); //always less than PI

        return (theta*obstacleRadii);
    }
    */
    rg_Point2D obstacleCenter = obstacle.getCenterPt();
    double     obstacleRadii = obstacle.getRadius();

    rg_Point2D vec1 = arcStartPt - obstacleCenter;
    rg_Point2D vec2 = arcEndPt - obstacleCenter;

    double theta = angleFromVec1toVec2(vec1, vec2); //always less than PI

    return (theta*obstacleRadii);
}



rg_ImplicitEquation ShortestPathFinder2D::make_line_equation(const rg_Point2D& pt1, const rg_Point2D& pt2) const
{
    //line equation : (y1-y2)x + (x2-x1)y + (x1y2 - x2y1) = 0

    rg_ImplicitEquation lineEq(1);
    lineEq.setCoeff(1, 1, 0.0);
    lineEq.setCoeff(1, 0, pt1.getY() - pt2.getY());
    lineEq.setCoeff(0, 1, pt2.getX() - pt1.getX());
    lineEq.setCoeff(0, 0, pt1.getX()*pt2.getY() - pt2.getX()*pt1.getY());

    return lineEq;
}



bool ShortestPathFinder2D::are_these_two_lines_intersected_in_this_interval(const rg_ImplicitEquation& lineEq1,
                                                                            const rg_ImplicitEquation& lineEq2, 
                                                                            double x_s, 
                                                                            double x_e) const
{
    // Line1. a1x + b1y + c1 = 0
    // Line2. a2x + b2y + c2 = 0
    // x = (b1c2 - b2c1) / (A),
    // y = (a2c1 - a1c2) / (A),
    // A = a1b2 - a2b1

    double a1 = lineEq1.getCoeff(1, 0);
    double b1 = lineEq1.getCoeff(0, 1);
    double c1 = lineEq1.getCoeff(0, 0);

    double a2 = lineEq2.getCoeff(1, 0);
    double b2 = lineEq2.getCoeff(0, 1);
    double c2 = lineEq2.getCoeff(0, 0);

    double A = a1*b2 - a2*b1;

    if (rg_ZERO(A))
    {
        return false;
    }

    double x = (b1*c2 - b2*c1) / A;

    return (((x >= x_s) && (x <= x_e)) || ((x <= x_s) && (x >= x_e)));
}




void ShortestPathFinder2D::insert_a_point_into_the_graph(VertexForSolutionPoolGraph* pointVtx, 
                                                         const vector<VertexBU2D*>& obstaclePool, 
                                                         BetaUniverse2D& betaUniverse, 
                                                         ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                         unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    if (obstaclePool.size() == 0)
    {
        return;
    }

    for (int i = 0; i < obstaclePool.size(); i++)
    {
        VertexBU2D* currObstacle = obstaclePool.at(i);

        generate_graph_incrementally_by_adding_one_obstacle_N_point_into_the_graph(currObstacle,
                                                                                   pointVtx,
                                                                                   m_BetaComponentFinder,
                                                                                   betaUniverse,
                                                                                   solutionPoolGraph, 
                                                                                   mapperForQTVtxToVtxOnDisk);
    }

}


void ShortestPathFinder2D::insert_a_point_into_the_graph_with_topological_path(VertexForSolutionPoolGraph * pointVtx, const vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    if (obstaclePool.size() == 0)
    {
        return;
    }

    for (int i = 0; i < obstaclePool.size(); i++)
    {
        VertexBU2D* currObstacle = obstaclePool.at(i);

        generate_graph_incrementally_by_adding_one_obstacle_N_point_into_the_graph_with_topological_path(currObstacle,
                                                                                   pointVtx,
                                                                                   topologicalPath,
                                                                                   m_BetaComponentFinder,
                                                                                   betaUniverse,
                                                                                   solutionPoolGraph, 
                                                                                   mapperForQTVtxToVtxOnDisk);
    }

}



void ShortestPathFinder2D::make_exterior_tangent_lines_of_two_circles(const rg_Circle2D& circle1,
											                          const rg_Circle2D& circle2,
											                          rg_ImplicitEquation& result1,
											                          rg_ImplicitEquation& result2)
{
    rg_Circle2D w1, w2;

	//the radius of w2 is less than that of w1
	if( circle1.getRadius() < circle2.getRadius() )
	{
		w1 = circle2;
		w2 = circle1;
	}
	else
	{
		w1 = circle1;
		w2 = circle2;
	}

	rg_Point2D c2cVector = w1.getCenterPt() - w2.getCenterPt();
	rg_REAL r = w1.getRadius() - w2.getRadius();
	rg_REAL length = c2cVector.magnitude();
	rg_REAL sine = r / length;
	rg_REAL cosine = sqrt( length*length - r*r ) / length;

	//rotate theta  /  -theta
	rg_Point2D normal1( c2cVector.getX() * cosine - c2cVector.getY() * sine ,
					 c2cVector.getX() * sine + c2cVector.getY() * cosine );
	rg_Point2D normal2( c2cVector.getX() * cosine + c2cVector.getY() * sine ,
					 c2cVector.getY() * cosine - c2cVector.getX() * sine );
	normal1 = normal1.getUnitVector();
	normal2 = normal2.getUnitVector();
	
	//rotate -PI/2  /  PI/2
	normal1.setPoint( normal1.getY() , -1 * normal1.getX() );
	normal2.setPoint( -1 * normal2.getY() , normal2.getX() );

	result1.setCoeff(1, 0, normal1.getX());
	result1.setCoeff(0, 1, normal1.getY());
	result1.setCoeff(0, 0, -1 * normal1.getX() * w2.getCenterPt().getX() 
						   - normal1.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius());

	result2.setCoeff(1, 0, normal2.getX());
    result2.setCoeff(0, 1, normal2.getY());
    result2.setCoeff(0, 0, -1 * normal2.getX() * w2.getCenterPt().getX()
    - normal2.getY() * w2.getCenterPt().getY()
    + w2.getRadius());

}


void ShortestPathFinder2D::make_interior_tangent_lines_of_two_circles(const rg_Circle2D& circle1,
    const rg_Circle2D& circle2,
    rg_ImplicitEquation& result1,
    rg_ImplicitEquation& result2)
{
    if (circle1.isIntersectWith(circle2))
    {
        result1.setDegree(0);
        result2.setDegree(0);

        return;
    }

    rg_Circle2D w1, w2;

    //the radius of w2 is less than that of w1
    if (circle1.getRadius() < circle2.getRadius())
    {
        w1 = circle2;
        w2 = circle1;
    }
    else
    {
        w1 = circle1;
        w2 = circle2;
    }

    rg_Point2D c2cVector = w1.getCenterPt() - w2.getCenterPt();
    rg_REAL r = w1.getRadius() + w2.getRadius();
    rg_REAL length = c2cVector.magnitude();
    rg_REAL sine = r / length;
    rg_REAL cosine = sqrt(length*length - r*r) / length;

    //rotate theta  /  -theta
    rg_Point2D normal1(c2cVector.getX() * cosine - c2cVector.getY() * sine,
        c2cVector.getX() * sine + c2cVector.getY() * cosine);
    rg_Point2D normal2(c2cVector.getX() * cosine + c2cVector.getY() * sine,
        c2cVector.getY() * cosine - c2cVector.getX() * sine);
    normal1 = normal1.getUnitVector();
    normal2 = normal2.getUnitVector();

    //rotate -PI/2  /  PI/2
    normal1.setPoint(normal1.getY(), -1 * normal1.getX());
    normal2.setPoint(-1 * normal2.getY(), normal2.getX());

    result1.setCoeff(1, 0, normal1.getX());
    result1.setCoeff(0, 1, normal1.getY());
    result1.setCoeff(0, 0, -1 * normal1.getX() * w2.getCenterPt().getX()
        - normal1.getY() * w2.getCenterPt().getY()
        - w2.getRadius());

    result2.setCoeff(1, 0, normal2.getX());
    result2.setCoeff(0, 1, normal2.getY());
    result2.setCoeff(0, 0, -1 * normal2.getX() * w2.getCenterPt().getX()
        - normal2.getY() * w2.getCenterPt().getY()
        - w2.getRadius());
}


rg_Point2D ShortestPathFinder2D::compute_tangent_point_between_line_and_circle(const rg_Circle2D& circle, const rg_ImplicitEquation& line)
{
    rg_Point2D tangentPt;

    double x = circle.getCenterPt().getX();
    double y = circle.getCenterPt().getY();
    double a = line.getCoeff(1, 0);
    double b = line.getCoeff(0, 1);
    double c = line.getCoeff(0, 0);
    double a2_b2 = pow(a, 2) + pow(b, 2);

    if (rg_EQ(a2_b2, 1.0))
    {
        tangentPt.setX(-a*line.evaluateImpEquation(x, y) + x);
        tangentPt.setY(-b*line.evaluateImpEquation(x, y) + y);
    }
    else
    {
        tangentPt.setX(-a*line.evaluateImpEquation(x, y) / a2_b2 + x);
        tangentPt.setY(-b*line.evaluateImpEquation(x, y) / a2_b2 + y);
    }

    return tangentPt;
}



void ShortestPathFinder2D::sort_index_of_distances_in_non_decreasing_order(const vector<double>& distances, vector<int>& indexOfDistancesWithIncreasingOrder) const
{
    //1. make pairs of distance and index
    vector<pair<double, int>> distanceNIndexPairs;

    for (int i = 0; i < distances.size(); i++)
    {
        distanceNIndexPairs.push_back(make_pair(distances.at(i), i));
    }

    //2. sort the pairs
    sort(distanceNIndexPairs.begin(), distanceNIndexPairs.end(), ShortestPathFinder2D::compare_pair_of_distance_N_index_in_non_decreasing_order);

    //3. copy the index
    for (int i = 0; i < distanceNIndexPairs.size(); i++)
    {
        int index = distanceNIndexPairs.at(i).second;
        indexOfDistancesWithIncreasingOrder.push_back(index);
    }
}



bool ShortestPathFinder2D::compare_pair_of_distance_N_index_in_non_decreasing_order(const pair<double, int>& pair1, const pair<double, int>& pair2)
{
    if (pair1.first < pair2.first)
    {
        return true;
    }
    else
    {
        return false;
    }
}


ShortestPathSolutionPoolGraph * ShortestPathFinder2D::create_solution_pool_of_BFS(const ShortestPathSolutionPoolGraph & solutionPoolGraph, const int& index)
{
    m_SolutionPoolGraphsByTPT.at(index) = solutionPoolGraph;

    return (&m_SolutionPoolGraphsByTPT.at(index));
}


void ShortestPathFinder2D::find_prev_N_next_vertices_of_this_new_vertex_without_intersection_with_obstacles(list<VertexForSolutionPoolGraph*>& prevVertices,
                                                                                                            list<VertexForSolutionPoolGraph*>& nextVertices,
                                                                                                            VertexForSolutionPoolGraph* newVertexOnObstacle, 
                                                                                                            VertexBU2D* obstacle, 
                                                                                                            const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex, 
                                                                                                            const list<VertexBU2D*>& neighborIntersectedObstacles) const
{
    list<double> anglesOfNeighborIntersectedObstacles;
    get_angle_of_neighbor_intersected_obstacles(anglesOfNeighborIntersectedObstacles, neighborIntersectedObstacles, obstacle);

    double angleOfNewVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(newVertexOnObstacle, obstacle);

    for (list<VertexForSolutionPoolGraph*>::const_iterator it_VertexOnObstacle = verticesOnObstacleBeforeMakingNewVertex.begin();
        it_VertexOnObstacle != verticesOnObstacleBeforeMakingNewVertex.end();
        it_VertexOnObstacle++)
    {
        VertexForSolutionPoolGraph* oldVertexOnObstacle = (*it_VertexOnObstacle);

        int    intervalOption = 0;
        double angleOfCurrVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(oldVertexOnObstacle, obstacle);

        if (do_any_of_two_intervals_include_none_of_obstacles(angleOfCurrVertexOnObstacle,
            angleOfNewVertexOnObstacle,
            anglesOfNeighborIntersectedObstacles,
            intervalOption))
        {

            switch (intervalOption)
            {
            case 1:
            {
                      prevVertices.push_back(oldVertexOnObstacle);
                      break;
            }

            case 2:
            {
                      nextVertices.push_back(oldVertexOnObstacle);
                      break;
            }

            default:
                break;
            }
        }
    }
}



void ShortestPathFinder2D::find_nearest_two_vertices_including_this_new_vertex(VertexForSolutionPoolGraph*& prevVertex, VertexForSolutionPoolGraph*& nextVertex, const list<VertexForSolutionPoolGraph*>& prevVertices, const list<VertexForSolutionPoolGraph*>& nextVertices, VertexForSolutionPoolGraph* newVertexOnObstacle, VertexBU2D* obstacle) const
{

}



void ShortestPathFinder2D::filter_vertices_for_making_arc_edges(VertexBU2D* obstacle, 
                                                                VertexForSolutionPoolGraph* newVertexOnObstacle,
                                                                const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex,
                                                                const BetaComponentFinder& betaComponentFinder,
                                                                list<pair<VertexForSolutionPoolGraph*,VertexForSolutionPoolGraph*>>& vertexPairsForArcEdges) const
{

    list<VertexBU2D*> neighborIntersectedObstacles;
    get_neighbor_obstacles_interseted_with_this_obstacle(neighborIntersectedObstacles, obstacle, betaComponentFinder);

    switch (neighborIntersectedObstacles.size())
    {
    case 0:
    {
        VertexForSolutionPoolGraph* prevOldVertex = NULL;
        VertexForSolutionPoolGraph* nextOldVertex = NULL;

        find_nearest_two_vertices_including_this_new_vertex(prevOldVertex, 
                                                            nextOldVertex, 
                                                            newVertexOnObstacle, 
                                                            verticesOnObstacleBeforeMakingNewVertex, 
                                                            obstacle);


        if (m_QTEdgesOfTopologicalPath.size() == 0)
        {
            vertexPairsForArcEdges.push_back(make_pair(prevOldVertex, newVertexOnObstacle));
            vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, nextOldVertex));
        }
        else
        {
            double angleOfPrevVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(prevOldVertex, obstacle);
            double angleOfNextVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(nextOldVertex, obstacle);
            double angleOfNewVertexOnObstacle  = get_angle_of_this_vertex_with_regard_to_this_obstacle(newVertexOnObstacle, obstacle);
            
            list<double> anglesOfNeighborQTEdgesNotInTopologicalpath;
            get_angle_of_neighbor_qt_edges_NOT_in_topological_path(anglesOfNeighborQTEdgesNotInTopologicalpath, obstacle);
            
            if(do_any_of_two_intervals_include_none_of_obstacles(angleOfPrevVertexOnObstacle,
                                                                 angleOfNewVertexOnObstacle,
                                                                 anglesOfNeighborQTEdgesNotInTopologicalpath))
            {
                vertexPairsForArcEdges.push_back(make_pair(prevOldVertex, newVertexOnObstacle));
            }
            
            if(do_any_of_two_intervals_include_none_of_obstacles(angleOfNewVertexOnObstacle,
                                                                 angleOfNextVertexOnObstacle,
                                                                 anglesOfNeighborQTEdgesNotInTopologicalpath))
            {
                vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, nextOldVertex));
            }
        }
        
        break;
    }

    default: // size > 0
    {
        VertexForSolutionPoolGraph* prevOldVertex = NULL;
        VertexForSolutionPoolGraph* nextOldVertex = NULL;

        find_nearest_two_vertices_including_this_new_vertex(prevOldVertex, 
                                                            nextOldVertex, 
                                                            newVertexOnObstacle, 
                                                            verticesOnObstacleBeforeMakingNewVertex, 
                                                            obstacle);

        double angleOfPrevVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(prevOldVertex, obstacle);
        double angleOfNextVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(nextOldVertex, obstacle);
        double angleOfNewVertexOnObstacle  = get_angle_of_this_vertex_with_regard_to_this_obstacle(newVertexOnObstacle, obstacle);

        list<double> anglesOfNeighborIntersectedObstacles;
        get_angle_of_neighbor_intersected_obstacles(anglesOfNeighborIntersectedObstacles, neighborIntersectedObstacles, obstacle);

        if (m_QTEdgesOfTopologicalPath.size() == 0)
        {

            if (do_any_of_two_intervals_include_none_of_obstacles(angleOfPrevVertexOnObstacle,
                angleOfNewVertexOnObstacle,
                anglesOfNeighborIntersectedObstacles))
            {
                vertexPairsForArcEdges.push_back(make_pair(prevOldVertex, newVertexOnObstacle));
            }

            if (do_any_of_two_intervals_include_none_of_obstacles(angleOfNewVertexOnObstacle,
                angleOfNextVertexOnObstacle,
                anglesOfNeighborIntersectedObstacles))
            {
                vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, nextOldVertex));
            }
        }
        else
        {
            list<double> anglesOfNeighborQTEdgesNotInTopologicalpath;
            get_angle_of_neighbor_qt_edges_NOT_in_topological_path(anglesOfNeighborQTEdgesNotInTopologicalpath, obstacle);

            if (do_any_of_two_intervals_include_none_of_obstacles(angleOfPrevVertexOnObstacle,
                angleOfNewVertexOnObstacle,
                anglesOfNeighborIntersectedObstacles))
            {
                if (do_any_of_two_intervals_include_none_of_obstacles(angleOfPrevVertexOnObstacle,
                    angleOfNewVertexOnObstacle,
                    anglesOfNeighborQTEdgesNotInTopologicalpath))
                {

                    vertexPairsForArcEdges.push_back(make_pair(prevOldVertex, newVertexOnObstacle));
                }
            }

            if (do_any_of_two_intervals_include_none_of_obstacles(angleOfNewVertexOnObstacle,
                angleOfNextVertexOnObstacle,
                anglesOfNeighborIntersectedObstacles))
            {
                if (do_any_of_two_intervals_include_none_of_obstacles(angleOfNewVertexOnObstacle,
                    angleOfNextVertexOnObstacle,
                    anglesOfNeighborQTEdgesNotInTopologicalpath))
                {
                    vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, nextOldVertex));
                }
            }
        }

        break;
    }
    }


    /*
    list<VertexBU2D*> neighborIntersectedObstacles;
    get_neighbor_obstacles_interseted_with_this_obstacle(neighborIntersectedObstacles, obstacle, betaComponentFinder);

    switch (neighborIntersectedObstacles.size())
    {
    case 0:
    {
              for (list<VertexForSolutionPoolGraph*>::const_iterator it_VertexOnObstacle = verticesOnObstacleBeforeMakingNewVertex.begin();
                   it_VertexOnObstacle != verticesOnObstacleBeforeMakingNewVertex.end();
                   it_VertexOnObstacle++)
              {
                  VertexForSolutionPoolGraph* oldVertexOnObstacle = (*it_VertexOnObstacle);

                  rg_Point2D center(obstacle->getCoord());
                  rg_Point2D newVertexPoint(newVertexOnObstacle->get_tangent_point_to_disk());
                  rg_Point2D oldVertexPoint(oldVertexOnObstacle->get_tangent_point_to_disk());

                  rg_Point2D vecCenterToNew = newVertexPoint - center;
                  rg_Point2D vecCenterToOld = oldVertexPoint - center;

                  double theta = angleFromVec1toVec2(vecCenterToNew, vecCenterToOld);

                  if (theta > rg_PI)
                  {
                      vertexPairsForArcEdges.push_back(make_pair(oldVertexOnObstacle, newVertexOnObstacle));
                  }
                  else
                  {
                      vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, oldVertexOnObstacle));
                  }
              }
              
              break;
    }

    default: // size > 0
    {

               list<double> anglesOfNeighborIntersectedObstacles;
               get_angle_of_neighbor_intersected_obstacles(anglesOfNeighborIntersectedObstacles, neighborIntersectedObstacles, obstacle);

               double angleOfNewVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(newVertexOnObstacle, obstacle);

               for (list<VertexForSolutionPoolGraph*>::const_iterator it_VertexOnObstacle = verticesOnObstacleBeforeMakingNewVertex.begin();
                   it_VertexOnObstacle != verticesOnObstacleBeforeMakingNewVertex.end();
                   it_VertexOnObstacle++)
               {
                   VertexForSolutionPoolGraph* oldVertexOnObstacle = (*it_VertexOnObstacle);

                   int    intervalOption = 0;
                   double angleOfCurrVertexOnObstacle = get_angle_of_this_vertex_with_regard_to_this_obstacle(oldVertexOnObstacle, obstacle);

                   if (do_any_of_two_intervals_include_none_of_obstacles(angleOfCurrVertexOnObstacle,
                                                                         angleOfNewVertexOnObstacle,
                                                                         anglesOfNeighborIntersectedObstacles,
                                                                         intervalOption))
                   {

                       switch (intervalOption)
                       {
                       case 1:
                       {
                                 vertexPairsForArcEdges.push_back(make_pair(oldVertexOnObstacle, newVertexOnObstacle));
                                 break;
                       }

                       case 2:
                       {
                                 vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, oldVertexOnObstacle));
                                 break;
                       }

                       default:
                           break;
                       }
                   }
               }

               break;
    }
    }
 */
}

void ShortestPathFinder2D::filter_vertices_for_making_arc_edges_when_NONE_verices_on_obstacle(VertexBU2D * obstacle, VertexForSolutionPoolGraph * newVertexOnObstacle, const BetaComponentFinder & betaComponentFinder, list<pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>>& vertexPairsForArcEdges) const
{
    list<VertexBU2D*> neighborIntersectedObstacles;
    get_neighbor_obstacles_interseted_with_this_obstacle(neighborIntersectedObstacles, obstacle, betaComponentFinder);

    switch (neighborIntersectedObstacles.size())
    {
    case 0:
    {
        if (m_QTEdgesOfTopologicalPath.size() == 0)
        {
            vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, newVertexOnObstacle));
        }
        else
        {
            list<double> anglesOfNeighborQTEdgesNotInTopologicalpath;
            get_angle_of_neighbor_qt_edges_NOT_in_topological_path(anglesOfNeighborQTEdgesNotInTopologicalpath, obstacle);

            if (anglesOfNeighborQTEdgesNotInTopologicalpath.size() == 0)
            {
                vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, newVertexOnObstacle));
            }
        }
        
        break;
    }

    default: // size > 0
    {
        break;
    }
    }
}



void ShortestPathFinder2D::find_nearest_two_vertices_including_this_new_vertex(VertexForSolutionPoolGraph*& prevVertex, 
                                                                               VertexForSolutionPoolGraph*& nextVertex, 
                                                                               VertexForSolutionPoolGraph* newVertexOnObstacle,
                                                                               const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex,
                                                                               VertexBU2D* obstacle) const
{
    if (verticesOnObstacleBeforeMakingNewVertex.size() == 0)
    {
        return;
    }

    list<VertexForSolutionPoolGraph*>::const_iterator it_CurrVertexOnObstacle = verticesOnObstacleBeforeMakingNewVertex.begin();

    while (it_CurrVertexOnObstacle != verticesOnObstacleBeforeMakingNewVertex.end())
    {
        VertexForSolutionPoolGraph* currVertex = *it_CurrVertexOnObstacle;

        if (get_angle_of_this_vertex_with_regard_to_this_obstacle(newVertexOnObstacle, obstacle)
          < get_angle_of_this_vertex_with_regard_to_this_obstacle(currVertex, obstacle))
        {
            if (it_CurrVertexOnObstacle == verticesOnObstacleBeforeMakingNewVertex.begin())
            {
                nextVertex = *it_CurrVertexOnObstacle;
                prevVertex = *(--verticesOnObstacleBeforeMakingNewVertex.end());
            }
            else
            {
                nextVertex = *it_CurrVertexOnObstacle;
                prevVertex = *(--it_CurrVertexOnObstacle);
            }

            break;
        }

        ++it_CurrVertexOnObstacle;
    }

    if (it_CurrVertexOnObstacle == verticesOnObstacleBeforeMakingNewVertex.end())
    {
        nextVertex = *verticesOnObstacleBeforeMakingNewVertex.begin();
        prevVertex = *(--verticesOnObstacleBeforeMakingNewVertex.end());
    }
}



void ShortestPathFinder2D::get_neighbor_obstacles_interseted_with_this_obstacle(list<VertexBU2D*>& neighborIntersectedObstacles,
                                                                                VertexBU2D* obstacle, 
                                                                                const BetaComponentFinder& betaComponentFinder) const
{
    BetaComponent2D* betaComponentOfObstacle = betaComponentFinder.findIncludedComponent(obstacle);
    list<EdgeBU2D*>  edgesInSameComponent    = betaComponentOfObstacle->getEdges();

    for (list<EdgeBU2D*>::const_iterator it_EdgeInSameComponent = edgesInSameComponent.begin();
        it_EdgeInSameComponent != edgesInSameComponent.end();
        it_EdgeInSameComponent++)
    {
        EdgeBU2D* currEdgeInSameComponent = *it_EdgeInSameComponent;

        if (currEdgeInSameComponent->getStartVertex() == obstacle)
        {
            neighborIntersectedObstacles.push_back(currEdgeInSameComponent->getEndVertex());
        }
        else if (currEdgeInSameComponent->getEndVertex() == obstacle)
        {
            neighborIntersectedObstacles.push_back(currEdgeInSameComponent->getStartVertex());
        }
    }
    

    /*
    BetaComponent2D* betaComponentOfObstacle = betaComponentFinder.findIncludedComponent(obstacle);
    list<VertexBU2D*> verticesInSameComponent = betaComponentOfObstacle->getVertices();

    rg_dList<VertexBU2D*> neighborVetices;
    obstacle->getNeighborVertices(neighborVetices);

    neighborVetices.reset4Loop();
    while (neighborVetices.setNext4Loop())
    {
        VertexBU2D* currVertex = neighborVetices.getEntity();

        list<VertexBU2D*>::iterator it_FoundVertex = find(verticesInSameComponent.begin(), verticesInSameComponent.end(), currVertex);

        if (it_FoundVertex != verticesInSameComponent.end())
        {
            neighborIntersectedObstacles.push_back(currVertex);
        }
    }
    */
}



void ShortestPathFinder2D::get_angle_of_neighbor_intersected_obstacles(list<double>& anglesOfNeighborIntersectedObstacles, 
                                                                       const list<VertexBU2D*>& neighborIntersectedObstacles, 
                                                                       VertexBU2D* obstacle) const
{
    for (list<VertexBU2D*>::const_iterator it_NIObstacle = neighborIntersectedObstacles.begin();
        it_NIObstacle != neighborIntersectedObstacles.end();
        it_NIObstacle++)
    {
        VertexBU2D* currNIObstace = (*it_NIObstacle);

        double angle = get_angle_of_start_to_end_point(obstacle->getCoord(), currNIObstace->getCoord());

        anglesOfNeighborIntersectedObstacles.push_back(angle);
    }
}


void ShortestPathFinder2D::get_angle_of_neighbor_qt_edges_NOT_in_topological_path(list<double>& anglesOfNeighborQTEdges, VertexBU2D * obstacle) const
{
    rg_dList<EdgeBU2D*> incidentQTEdges;
    obstacle->getIncidentEdges(incidentQTEdges);

    incidentQTEdges.reset4Loop();
    while (incidentQTEdges.setNext4Loop())
    {
        EdgeBU2D* qtEdge = incidentQTEdges.getEntity();

        if (m_QTEdgesOfTopologicalPath.find(qtEdge) == m_QTEdgesOfTopologicalPath.end())
        {
            double angle = 0.0;

            if (qtEdge->getStartVertex() == obstacle)
            {
                if (qtEdge->getEndVertex()->isVirtual())
                {
                    angle = get_angle_of_start_to_end_point(qtEdge->getStartVertex()->getCoord(), qtEdge->computeASamplePointOnThisQTEdge());
                }
                else
                {
                    angle = get_angle_of_start_to_end_point(qtEdge->getStartVertex()->getCoord(), qtEdge->getEndVertex()->getCoord());
                }
            }
            else // qtEdge->getEndVertex() == obstacle
            {
                if (qtEdge->getStartVertex()->isVirtual())
                {
                    angle = get_angle_of_start_to_end_point(qtEdge->getEndVertex()->getCoord(), qtEdge->computeASamplePointOnThisQTEdge());
                }
                else
                {
                    angle = get_angle_of_start_to_end_point(qtEdge->getEndVertex()->getCoord(), qtEdge->getStartVertex()->getCoord());
                }
            }

            anglesOfNeighborQTEdges.push_back(angle);
        }
    }
}


void ShortestPathFinder2D::make_arc_edges(const list<pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>>& vertexPairsForArcEdges, 
                                          VertexBU2D* obstacle,
                                          ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    /* ORDER OF INPUT   
        vertexPairsForArcEdges.push_back(make_pair(prevOldVertex, newVertexOnObstacle));
        vertexPairsForArcEdges.push_back(make_pair(newVertexOnObstacle, nextOldVertex));
     
    */
    switch (vertexPairsForArcEdges.size())
    {
    case 2:
    {
        pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*> vertexPair1ForArcEdge(*vertexPairsForArcEdges.begin());
        pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*> vertexPair2ForArcEdge(*(++vertexPairsForArcEdges.begin()));

        VertexForSolutionPoolGraph* prevVertex = vertexPair1ForArcEdge.first;
        VertexForSolutionPoolGraph* newVertex  = vertexPair1ForArcEdge.second;
        VertexForSolutionPoolGraph* nextVertex = vertexPair2ForArcEdge.second;

        EdgeForSolutionPoolGraph* newArcEdgeOnObstacle = NULL;
        create_an_arc_edge(newArcEdgeOnObstacle, solutionPoolGraph);
        set_two_arc_edges(newArcEdgeOnObstacle, prevVertex, newVertex, nextVertex, obstacle, solutionPoolGraph);

        break;
    }

    case 1:
    {
        pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*> vertexPairForArcEdge(*vertexPairsForArcEdges.begin());

        EdgeForSolutionPoolGraph* arcEdgeOnObstacle = NULL;
        create_an_arc_edge(arcEdgeOnObstacle, solutionPoolGraph);
        set_an_arc_edge(arcEdgeOnObstacle, vertexPairForArcEdge.first, vertexPairForArcEdge.second, obstacle, solutionPoolGraph);

        break;
    }

    default:
        break;
    }

}



double ShortestPathFinder2D::get_angle_of_this_vertex_with_regard_to_this_obstacle(VertexForSolutionPoolGraph* vertexOnObstacle, VertexBU2D* obstacle) const
{
    double angle = vertexOnObstacle->get_angle();

    if (rg_NEG(angle))
    {
        rg_Point2D centerOfObstacle    = obstacle->getCoord();
        rg_Point2D tangentPtOnObstacle = vertexOnObstacle->get_tangent_point_to_disk();

        angle = get_angle_of_start_to_end_point(centerOfObstacle, tangentPtOnObstacle);

        vertexOnObstacle->set_angle(angle);
    }

    return angle;
}


bool ShortestPathFinder2D::do_this_interval_include_none_of_obstacles(const double& startAngle,
                                                                      const double& endAngle, 
                                                                      const list<double>& anglesOfNeighborIntersectedObstacles) const
{
    if (rg_EQ(startAngle, endAngle))
    {
        return true;
    }

    bool   doThisIntervalIncludeNoneOfObstacles = true;

    if (startAngle < endAngle)
    {
        for (list<double>::const_iterator it_AngleOfNIObstacle = anglesOfNeighborIntersectedObstacles.begin();
            it_AngleOfNIObstacle != anglesOfNeighborIntersectedObstacles.end();
            it_AngleOfNIObstacle++)
        {
            double currAngleOfNIObstacle = *it_AngleOfNIObstacle;

            bool isThisAngleInThisInterval = is_this_val_in_this_interval(currAngleOfNIObstacle, 
                                                                          startAngle, 
                                                                          endAngle);

            if (isThisAngleInThisInterval)
            {
                doThisIntervalIncludeNoneOfObstacles = false;
                break;
            }
        }
    }
    else  //startAngle > endAngle
    {
        for (list<double>::const_iterator it_AngleOfNIObstacle = anglesOfNeighborIntersectedObstacles.begin();
             it_AngleOfNIObstacle != anglesOfNeighborIntersectedObstacles.end();
             it_AngleOfNIObstacle++)
        {
            double currAngleOfNIObstacle = *it_AngleOfNIObstacle;

            bool isThisAngleInThisInterval = is_this_val_in_this_interval(currAngleOfNIObstacle, 
                                                                          startAngle, 
                                                                          2*rg_PI) 
                                           ||is_this_val_in_this_interval(currAngleOfNIObstacle, 
                                                                          0.0, 
                                                                          endAngle) ;

            if (isThisAngleInThisInterval)
            {
                doThisIntervalIncludeNoneOfObstacles = false;
                break;
            }
        }
    }

    return doThisIntervalIncludeNoneOfObstacles;
}



void ShortestPathFinder2D::find_geodesic_path_by_dijkstra(const ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    VertexForSolutionPoolGraph* startVtx = solutionPoolGraph.get_start_vtx();

    EntityAccessiblePriorityQ<VertexForSolutionPoolGraph*> priorityQ;

    list<VertexForSolutionPoolGraph*> allVertices;
    solutionPoolGraph.get_vertices(allVertices);
    
    // 1. initialize - vertex
    for (list<VertexForSolutionPoolGraph*>::iterator it_Vertex = allVertices.begin();
        it_Vertex != allVertices.end();
        it_Vertex++)
    {
        VertexForSolutionPoolGraph* currVtx = (*it_Vertex);

        if (currVtx == startVtx)
        {
            currVtx->set_accumulate_length_from_source_vtx(0.0);
        }
        else
        {
            currVtx->set_accumulate_length_from_source_vtx(DBL_MAX);
        }
        
        currVtx->set_prev_edge(NULL);

        priorityQ.push(currVtx, currVtx->get_accumulate_length_from_source_vtx());
    }


    // 2. REPEAT the follwings
    while (!priorityQ.empty())
    {
        VertexForSolutionPoolGraph* minVertex = priorityQ.pop();

        list<EdgeForSolutionPoolGraph*> neighborEdges;
        minVertex->get_all_edges(neighborEdges);
 
        // TANGENT LINE N ARC LINE
        for (list<EdgeForSolutionPoolGraph*>::iterator it_NeighborEdge = neighborEdges.begin();
            it_NeighborEdge != neighborEdges.end();
            it_NeighborEdge++)
        {
            EdgeForSolutionPoolGraph*    currEdge   = *it_NeighborEdge;
            VertexForSolutionPoolGraph*  nextVertex = currEdge->get_opposite_vtx(minVertex);

            double currLength = minVertex->get_accumulate_length_from_source_vtx() + currEdge->get_edge_length();

            if (currLength < nextVertex->get_accumulate_length_from_source_vtx())
            {
                nextVertex->set_accumulate_length_from_source_vtx(currLength);
                nextVertex->set_prev_edge(currEdge);

                priorityQ.replaceKey(nextVertex, currLength);
            }
        } 
    }
}



bool ShortestPathFinder2D::back_trace_to_start_vtx_and_find_geodesic_path(VertexForSolutionPoolGraph* startVtx, 
                                                                             VertexForSolutionPoolGraph* endVtx, 
                                                                             list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                             double& totalPathDistance) const
{
    bool isThereAPath = true;

    list<EdgeForSolutionPoolGraph*> path;

    VertexForSolutionPoolGraph* currVtx   = endVtx;
    EdgeForSolutionPoolGraph*  currEdge   = currVtx->get_prev_edge();

    while (currVtx != startVtx)
    {
        if (currEdge == NULL)
        {
            isThereAPath = false;
            break;
        }

        path.push_front(currEdge);

        currVtx  = currEdge->get_opposite_vtx(currVtx);
        currEdge = currVtx->get_prev_edge();
    }

    if (isThereAPath)
    {
        geodesicPath      = path;
        totalPathDistance = endVtx->get_accumulate_length_from_source_vtx();
    }

    return isThereAPath;
}


void ShortestPathFinder2D::find_connected_components_of_beta_shape(const BetaUniverse2D& betaUniverse, list<BetaComponent2D>& betaComponents, BetaComponentFinder& betaComponentFinder)
{
    m_BetaUniverse.findConnectedComponents(m_Probe.getRadius(), betaComponents, betaComponentFinder);

}


void ShortestPathFinder2D::compute_ellipse_filter(const rg_Circle2D& startPt, const rg_Circle2D& endPt, Ellipse2D& ellipseFilter) const
{
    // input : start point(s), end point (e)
    // center: (e + s) /2.0, distance d = (e - s)
    // angle : from x-axis to vector(e - s)
    // semi-major axis length a = rg_PI*d/2.0
    // semi-minor axis length b = sqrt( a^2 - (d/2)^2)
    //1. center 
    rg_Point2D center = (startPt.getCenterPt() + endPt.getCenterPt()) / 2.0;
    
    rg_Point2D vecFromStartToEnd = (endPt.getCenterPt() - startPt.getCenterPt());
    rg_Point2D vecForXAxis       = rg_Point2D(1.0, 0.0);

    double distance = vecFromStartToEnd.magnitude();
    double angle    = angleFromVec1toVec2(vecForXAxis, vecFromStartToEnd);
    double semiMajorAxisLength = rg_PI * distance / 4.0 ;
    double semiMinorAxisLength = sqrt(pow(semiMajorAxisLength, 2) - pow((distance / 2.0), 2));

    ellipseFilter = Ellipse2D(center, semiMajorAxisLength, semiMinorAxisLength, angle);
}


void ShortestPathFinder2D::mark_QT_vertex_whether_in_or_out_of_ellipse(BetaUniverse2D & betaUniverse, Ellipse2D & elipseFilter, unordered_map<VertexBU2D*, bool>& QTVertexValidation)
{
    if (QTVertexValidation.empty())
    {
        //1. GET ALL QT VERTEX AND CHECK WHETHER IN OR OUT OF ELLIPSE
        const rg_dList<VertexBU2D>& allQTVertices = betaUniverse.getVertices();

        allQTVertices.reset4Loop();
        while (allQTVertices.setNext4Loop())
        {
            VertexBU2D* currVertex = allQTVertices.getpEntity();

            if (currVertex->isVirtual())
            {
                QTVertexValidation.insert(make_pair(currVertex, false));
            }
            else
            {
                QTVertexValidation.insert(make_pair(currVertex, is_this_QT_vertex_inside_of_ellipse(currVertex, elipseFilter)));
            }
        }
    }
    else
    {
        for (unordered_map<VertexBU2D*, bool>::iterator it_VertexPair = QTVertexValidation.begin(); it_VertexPair != QTVertexValidation.end(); ++it_VertexPair)
        {
            VertexBU2D* currVertex = (*it_VertexPair).first;

            if (QTVertexValidation.at(currVertex) == true)
            {
                if (!is_this_QT_vertex_inside_of_ellipse(currVertex, elipseFilter))
                {
                    QTVertexValidation.at(currVertex) = false;
                }
            }
        }
    }

}


void ShortestPathFinder2D::mark_QT_edge_whether_in_or_out_of_ellipse(BetaUniverse2D & betaUniverse, const unordered_map<VertexBU2D*, bool>& QTVertexValidation, unordered_map<EdgeBU2D*, bool>& QTEdgeValidation)
{
    if (QTEdgeValidation.empty())
    {
        //2. CHECK EACH QT EDGE
        const rg_dList<EdgeBU2D>& allQTEdges = betaUniverse.getEdges();

        allQTEdges.reset4Loop();
        while (allQTEdges.setNext4Loop())
        {
            EdgeBU2D* currEdge = allQTEdges.getpEntity();

            QTEdgeValidation.insert(make_pair(currEdge, is_this_QT_edge_inside_of_ellipse(currEdge, QTVertexValidation)));
        }
    }
    else
    {
        for (unordered_map<EdgeBU2D*, bool>::iterator it_EdgePair = QTEdgeValidation.begin(); it_EdgePair != QTEdgeValidation.end(); ++it_EdgePair)
        {
            EdgeBU2D* currEdge = (*it_EdgePair).first;

            if (QTEdgeValidation.at(currEdge) == true)
            {
                if (!is_this_QT_edge_inside_of_ellipse(currEdge, QTVertexValidation))
                {
                    QTEdgeValidation.at(currEdge) = false;
                }
            }
        }

    }
   
}


void ShortestPathFinder2D::mark_QT_vertex_whether_in_or_out_of_box(BetaUniverse2D & betaUniverse, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth, unordered_map<VertexBU2D*, bool>& QTVertexValidation)
{
    if (QTVertexValidation.empty())
    {
        //1. GET ALL QT VERTEX AND CHECK WHETHER IN OR OUT OF ELLIPSE
        const rg_dList<VertexBU2D>& allQTVertices = betaUniverse.getVertices();

        allQTVertices.reset4Loop();
        while (allQTVertices.setNext4Loop())
        {
            VertexBU2D* currVertex = allQTVertices.getpEntity();

            if (currVertex->isVirtual())
            {
                QTVertexValidation.insert(make_pair(currVertex, false));
            }
            else
            {
                QTVertexValidation.insert(make_pair(currVertex, is_this_QT_vertex_inside_of_box(currVertex, lineSegmentForBoxHeight, halfLineSegmentForBoxWidth)));
            }
        }
    }
    else
    {
        for (unordered_map<VertexBU2D*, bool>::iterator it_VertexPair = QTVertexValidation.begin(); it_VertexPair != QTVertexValidation.end(); ++it_VertexPair)
        {
            VertexBU2D* currVertex = (*it_VertexPair).first;

            if (QTVertexValidation.at(currVertex) == true)
            {
                if (!is_this_QT_vertex_inside_of_box(currVertex, lineSegmentForBoxHeight, halfLineSegmentForBoxWidth))
                {
                    QTVertexValidation.at(currVertex) = false;
                }
            }
        }
    }

}


void ShortestPathFinder2D::mark_QT_edge_whether_in_or_out_of_box(BetaUniverse2D& betaUniverse, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth, const unordered_map<VertexBU2D*, bool>& QTVertexValidation, unordered_map<EdgeBU2D*, bool>& QTEdgeValidation)
{
    if (QTEdgeValidation.empty())
    {
        //2. CHECK EACH QT EDGE
        const rg_dList<EdgeBU2D>& allQTEdges = betaUniverse.getEdges();

        allQTEdges.reset4Loop();
        while (allQTEdges.setNext4Loop())
        {
            EdgeBU2D* currEdge = allQTEdges.getpEntity();

            QTEdgeValidation.insert(make_pair(currEdge, is_this_QT_edge_inside_of_box(currEdge, lineSegmentForBoxHeight, halfLineSegmentForBoxWidth, QTVertexValidation)));
        }
    }
    else
    {
        for (unordered_map<EdgeBU2D*, bool>::iterator it_EdgePair = QTEdgeValidation.begin(); it_EdgePair != QTEdgeValidation.end(); ++it_EdgePair)
        {
            EdgeBU2D* currEdge = (*it_EdgePair).first;

            if (QTEdgeValidation.at(currEdge) == true)
            {
                if (!is_this_QT_edge_inside_of_box(currEdge, lineSegmentForBoxHeight, halfLineSegmentForBoxWidth, QTVertexValidation))
                {
                    QTEdgeValidation.at(currEdge) = false;
                }
            }
        }

    }

}


void ShortestPathFinder2D::compute_box_filter(const rg_Circle2D & startPt, const rg_Circle2D & endPt, const double & probeMaxSpeed, const double & travelTime, rg_Line2D & lineSegmentForBoxHeight, double& halfLineSegmentForBoxWidth)
{
    lineSegmentForBoxHeight    = rg_Line2D(startPt.getCenterPt(), endPt.getCenterPt());
    halfLineSegmentForBoxWidth = probeMaxSpeed*travelTime;
}


void ShortestPathFinder2D::initialize_statistics()
{
    m_TimeStatistics.reset(SPF_TIME_SIZE_PHASE_ONE);
}


void ShortestPathFinder2D::finalize_statistics()
{
    m_TimeStatistics.setTime(SPF_TIME_TOTAL, m_TimeStatistics.time(SPF_TIME_CONSTRUCT_VD_FAMILY)
                                           + m_TimeStatistics.time(SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE)
                                           + m_TimeStatistics.time(SPF_TIME_FIND_PATH_BY_BFS)
                                           + m_TimeStatistics.time(SPF_TIME_FIND_PATH_BY_DFS)
                                           + m_TimeStatistics.time(SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH));

}


void ShortestPathFinder2D::change_optimal_path_in_old_graph_to_new_one(const ShortestPathSolutionPoolGraph& oldGraph, const ShortestPathSolutionPoolGraph& newGraph, list<EdgeForSolutionPoolGraph*>& geodesicPath) const
{
    unordered_map<EdgeForSolutionPoolGraph*, EdgeForSolutionPoolGraph*> oldToNewEdgeMapper;

    list<EdgeForSolutionPoolGraph*> oldEdges;
    list<EdgeForSolutionPoolGraph*> newEdges;

    oldGraph.get_edges(oldEdges);
    newGraph.get_edges(newEdges);

    list<EdgeForSolutionPoolGraph*>::const_iterator it_OldEdge = oldEdges.begin();
    list<EdgeForSolutionPoolGraph*>::const_iterator it_NewEdge = newEdges.begin();

    for (; it_OldEdge != oldEdges.end(); it_OldEdge++, it_NewEdge++)
    {
        EdgeForSolutionPoolGraph* oldEdge = *it_OldEdge;
        EdgeForSolutionPoolGraph* newEdge = *it_NewEdge;

        oldToNewEdgeMapper.insert(make_pair(oldEdge, newEdge));
    }


    list<EdgeForSolutionPoolGraph*> newGeodesic;

    for (list<EdgeForSolutionPoolGraph*>::const_iterator it_Path = geodesicPath.begin();
        it_Path != geodesicPath.end();
        it_Path++)
    {
        EdgeForSolutionPoolGraph* oldPathEdge = *it_Path;
        EdgeForSolutionPoolGraph* newPathEdge = oldToNewEdgeMapper.at(oldPathEdge);

        newGeodesic.push_back(newPathEdge);
    }

    geodesicPath.clear();
    geodesicPath = newGeodesic;
}


void ShortestPathFinder2D::make_mapper_from_QT_obstacle_pair_to_its_tangent_lines(VertexBU2D* obstacle1, VertexBU2D* obstacle2, 
                                                                                  const vector<rg_ImplicitEquation>& tangentLines, 
                                                                                  unordered_map<pair<VertexBU2D*, VertexBU2D*>, vector<rg_ImplicitEquation>, pair_hash>& obstaclePairNItsTangentLines) const
{
    if (obstacle1 < obstacle2)
    {
        pair<VertexBU2D*, VertexBU2D*> pairOfObstacle = make_pair(obstacle1, obstacle2);
        obstaclePairNItsTangentLines.insert(make_pair(pairOfObstacle, tangentLines));
    }
    else
    {
        pair<VertexBU2D*, VertexBU2D*> pairOfObstacle = make_pair(obstacle2, obstacle1);
        obstaclePairNItsTangentLines.insert(make_pair(pairOfObstacle, tangentLines));
    }
}


bool ShortestPathFinder2D::are_the_tangent_lines_of_this_obstacle_pair_already_computed(VertexBU2D* obstacle1, VertexBU2D* obstacle2) const
{
    pair<VertexBU2D*, VertexBU2D*> pairOfObstacle;

    if (obstacle1 < obstacle2)
    {
       pairOfObstacle = make_pair(obstacle1, obstacle2);
    }
    else
    {
       pairOfObstacle = make_pair(obstacle2, obstacle1);
    }

    bool areTheTangentLinesOfThisObstaclePairAlreadyComputed = true;

    if (m_ObstaclePairNItsTangentLines.find(pairOfObstacle) == m_ObstaclePairNItsTangentLines.end())
    {
        areTheTangentLinesOfThisObstaclePairAlreadyComputed = false;
    }

    return areTheTangentLinesOfThisObstaclePairAlreadyComputed;
}


void ShortestPathFinder2D::get_tangent_lines_of_this_obstacle_pair_from_already_computed(VertexBU2D* obstacle1, VertexBU2D* obstacle2, 
                                                                                         vector<rg_ImplicitEquation>& tangentLines) const
{
    if (are_the_tangent_lines_of_this_obstacle_pair_already_computed(obstacle1, obstacle2))
    {
        pair<VertexBU2D*, VertexBU2D*> pairOfObstacle;

        if (obstacle1 < obstacle2)
        {
            pairOfObstacle = make_pair(obstacle1, obstacle2);
        }
        else
        {
            pairOfObstacle = make_pair(obstacle2, obstacle1);
        }

        tangentLines = m_ObstaclePairNItsTangentLines.at(pairOfObstacle);
    }
}


void ShortestPathFinder2D::make_set_of_obstacles(const vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, vector<VertexBU2D*>& obstaclePool) const
{
    set<VertexBU2D*> obstacleSet;

    for (int i = 0; i < obstaclePairs.size(); i++)
    {
        VertexBU2D* obstacle1 = obstaclePairs.at(i).first;
        VertexBU2D* obstacle2 = obstaclePairs.at(i).second;

        obstacleSet.insert(obstacle1);
        obstacleSet.insert(obstacle2);
    }

    obstaclePool.insert(obstaclePool.end(), obstacleSet.begin(), obstacleSet.end());
}


bool ShortestPathFinder2D::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order(const VertexBU2D* obstacle1, const VertexBU2D* obstacle2)
{
    bool isLargerThan = false;

    rg_Circle2D disk1(obstacle1->getCircle());
    rg_Circle2D disk2(obstacle2->getCircle());

    if (disk1.getRadius() > disk2.getRadius())
    {
        isLargerThan = true;
    }
    else if (disk1.getRadius() == disk2.getRadius())
    {

        if (disk1.getCenterPt().getX() > disk2.getCenterPt().getX())
        {
            isLargerThan = true;
        }
        else if (disk1.getCenterPt().getX() == disk2.getCenterPt().getX())
        {


            if (disk1.getCenterPt().getY() > disk2.getCenterPt().getY())
            {
                isLargerThan = true;
            }
        }
    }

    return isLargerThan;
}



bool ShortestPathFinder2D::are_some_bounding_edges_of_this_QT_face_exterior(FaceBU2D* qtFace) const
{
    bool areSomeBoundingEdgesOfThisQTFaceExterior = false;

    rg_dList<EdgeBU2D*> boundingQTEdges;
    qtFace->getBoundingEdges(boundingQTEdges);

    boundingQTEdges.reset4Loop();
    while (boundingQTEdges.setNext4Loop())
    {
        EdgeBU2D* currQTEdge = boundingQTEdges.getEntity();

        if (is_this_QT_edge_in_an_exterior_state(currQTEdge))
        {
            areSomeBoundingEdgesOfThisQTFaceExterior = true;
            break;
        }
    }

    return areSomeBoundingEdgesOfThisQTFaceExterior;
}




/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                   TEMP. CODE                                                           */
/*                                                                                                                        */
/*                                                                                                                        */
/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/



rg_Circle2D ShortestPathFinder2D::find_next_stop_point_by_greedy_method(const rg_Circle2D & targetPt, 
                                                                        const rg_Circle2D & currProbe, 
                                                                        const double & speed, 
                                                                        const double & travelTime, 
                                                                        const list<rg_Circle2D>& obstacles)
{
    rg_Circle2D startPt(currProbe.getCenterPt(), 0.0);

    if (startPt.getCenterPt() == targetPt.getCenterPt())
    {
        return rg_Circle2D();
    }

    initialize_statistics();


    //1. set input
    set_input_argments(startPt, targetPt, currProbe, obstacles, speed, travelTime);


    ___start_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);

    //2. construct Voronoi diaram and quasi-triangulation and beta-shape
    VoronoiDiagram2DC VD2DC;

    construct_voronoi_diagram(VD2DC);
    construct_quasi_triangulation_from_VD(m_QuasiTriangulation, VD2DC);
    construct_beta_shape_from_QT(m_BetaUniverse, m_QuasiTriangulation);

    find_connected_components_of_beta_shape(m_BetaUniverse, m_BetaComponents, m_BetaComponentFinder);

    ___end_clock(SPF_TIME_CONSTRUCT_VD_FAMILY);

    
    rg_Circle2D nextStopPoint = find_next_stop_point(m_BetaUniverse, m_PathByAllPathSearch);

    finalize_statistics();

    return nextStopPoint;
}


rg_Circle2D ShortestPathFinder2D::find_next_stop_point(BetaUniverse2D & betaUniverse, list<EdgeForSolutionPoolGraph*>& geodesicPath)
{
    rg_Circle2D nextStopPoint;

    //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  betaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = betaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = betaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {
                nextStopPoint = find_next_stop_point_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace,
                                                                                                              betaUniverse,
                                                                                                              m_SolutionPoolGraphByAllPathSearch,
                                                                                                              m_PathByAllPathSearch,
                                                                                                              m_PathTotalDistanceByAllPathSearch);
            }
            else
            {
                nextStopPoint = find_next_stop_point_in_general_case(betaUniverse,
                                                                     m_SolutionPoolGraphByBFS,
                                                                     geodesicPath,
                                                                     m_PathTotalDistanceByBFS);
            }
          
            break;
        }

        case false:
        {
            nextStopPoint = find_next_stop_point_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                                         m_SolutionPoolGraphByAllPathSearch,
                                                                                                                                         m_PathByAllPathSearch, 
                                                                                                                                         m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }


    return nextStopPoint;
}


void ShortestPathFinder2D::compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS(BetaUniverse2D & betaUniverse, const unordered_map<VertexBU2D*, bool>& isThisQTVertexInsideOfEllipse, const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse)
{
       //1. check intersection btw pq line and obstacles
    rg_ImplicitEquation lineEqFromSourceToDestination = make_line_equation(m_StartPt.getCenterPt(), m_EndPt.getCenterPt());
    
    switch (is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEqFromSourceToDestination, 
                                                                                  m_StartPt.getCenterPt(), 
                                                                                  m_EndPt.getCenterPt(), 
                                                                                  betaUniverse))
    {
        case true:
        {
            // 1. Get start and end face of QT
            FaceBU2D* qtStartFace = betaUniverse.findFaceContainingInputPoint(m_StartPt.getCenterPt());
            FaceBU2D* qtEndFace   = betaUniverse.findFaceContainingInputPoint(m_EndPt.getCenterPt());

            if (qtStartFace == qtEndFace)
            {
                find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(qtStartFace, 
                                                                                                betaUniverse,
                                                                                                m_SolutionPoolGraphByAllPathSearch,
                                                                                                m_PathByAllPathSearch,
                                                                                                m_PathTotalDistanceByAllPathSearch);
            }
            else
            {
                //1. find one topological path
                list<EdgeBU2D*> topologicalPathInEllipseFilter;
                find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, topologicalPathInEllipseFilter, isThisQTEdgeInsideOfEllipse, betaUniverse);

                switch (!topologicalPathInEllipseFilter.empty()) //TOPOLOGICAL PATH EXIST
                {
                case true:
                {
                    //2. find geodesic path
                    find_one_geodesic_path_from_candidate_obstacles(isThisQTVertexInsideOfEllipse,
                                                                    betaUniverse, 
                                                                    m_SolutionPoolGraphByAllPathSearch, 
                                                                    m_PathByAllPathSearch, 
                                                                    m_PathTotalDistanceByAllPathSearch);
                    break;
                }
                
                case false:
                {     
                    switch(solverForBarrierProblemType)
                    { 
                    case TPT_METHOD:
                    {
                        unordered_map<EdgeBU2D*, bool> isThisQTEdgeInsideOfExpandedEllipse;
                        compute_several_geodesic_paths_by_topological_path_tree_with_filtered_QT_edges(qtStartFace, qtEndFace, isThisQTEdgeInsideOfEllipse, isThisQTEdgeInsideOfExpandedEllipse, betaUniverse);
                        sort_index_of_distances_in_non_decreasing_order(m_PathTotalDistancesByTPT, m_IndexOfBFSsWithIncreasingOrderOfDistance);

                        break;
                    }

                    case EXPAND_PATH_ALONG_BARRIER_METHOD:
                    {
                        unordered_map<VertexBU2D*, bool> isThisQTVertexInsideOfExpandedEllipse;
                        expand_candidate_QT_vertices_by_including_on_a_barrier(isThisQTVertexInsideOfEllipse, isThisQTVertexInsideOfExpandedEllipse, betaUniverse);

                        //find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(qtStartFace, qtEndFace, m_TopologicalPathByAllPathSearch, isThisQTEdgeInsideOfExpandedEllipse, betaUniverse);

                        //2. find geodesic path
                        find_one_geodesic_path_from_candidate_obstacles(isThisQTVertexInsideOfExpandedEllipse,
                                                                        betaUniverse, 
                                                                        m_SolutionPoolGraphByAllPathSearch, 
                                                                        m_PathByAllPathSearch, 
                                                                        m_PathTotalDistanceByAllPathSearch);

                        
                        break;
                    }

                    default:
                        break;
                    }
                    //compute_several_geodesic_paths_by_topological_path_tree(qtStartFace, qtEndFace, betaUniverse);
                    //sort_index_of_distances_in_non_decreasing_order(m_PathTotalDistancesByTPT, 
                    
                    break;
                }

                default:
                    break;
                }
            }
          
            break;
        }

        case false:
        {
            find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(lineEqFromSourceToDestination,
                                                                                                                           m_SolutionPoolGraphByAllPathSearch,
                                                                                                                           m_PathByAllPathSearch, m_PathTotalDistanceByAllPathSearch);
            break;
        }

    default:
        break;
    }


}



void ShortestPathFinder2D::mark_QT_edge_whether_in_or_out_of_ellipse(BetaUniverse2D& betaUniverse, 
                                                                     Ellipse2D& elipseFilter, 
                                                                     unordered_map<EdgeBU2D*, bool>& QTEdgeValidation)
{
    unordered_map<VertexBU2D*, bool> QTVertexValidation;

    //1. GET ALL QT VERTEX AND CHECK WHETHER IN OR OUT OF ELLIPSE
    const rg_dList<VertexBU2D>& allQTVertices = betaUniverse.getVertices();

    allQTVertices.reset4Loop();
    while (allQTVertices.setNext4Loop())
    {
        VertexBU2D* currVertex = allQTVertices.getpEntity();

        if (currVertex->isVirtual())
        {
            QTVertexValidation.insert(make_pair(currVertex, false));
        }
        else
        {
            QTVertexValidation.insert(make_pair(currVertex, is_this_QT_vertex_inside_of_ellipse(currVertex, elipseFilter)));
        }
    }

    //2. CHECK EACH QT EDGE
    const rg_dList<EdgeBU2D>& allQTEdges = betaUniverse.getEdges();

    allQTEdges.reset4Loop();
    while (allQTEdges.setNext4Loop())
    {
        EdgeBU2D* currEdge = allQTEdges.getpEntity();

        QTEdgeValidation.insert(make_pair(currEdge, is_this_QT_edge_inside_of_ellipse(currEdge, QTVertexValidation)));
    }
}