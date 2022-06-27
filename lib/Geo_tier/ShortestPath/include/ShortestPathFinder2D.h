#ifndef _SHORTEST_PATH_FINDER_2D_
#define _SHORTEST_PATH_FINDER_2D_

//#define OBSTACLE_PAIR_ON_QT_EDGE
#define  ELLIPSE_FILTER

#include "rg_Circle2D.H"
#include "rg_Point2D.h"
#include "ShortestPathSolutionPoolGraph.h"
#include "EdgeForSolutionPoolGraph.h"
#include "VertexForSolutionPoolGraph.h"
#include "QuasiTriangulation2D.h"
#include "BetaUniverse2D.h"
#include "ConstForBetaComplex.h"
#include "TimeStatistics.h"
#include "TopologicalPathTree.h"
#include "EntityAccessiblePriorityQ.h"
#include "Ellipse2D.H"
#include "rg_Line2D.h"

#include <list>
#include <queue>
#include <unordered_map>
#include <unordered_set>
using namespace std;

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};



class VoronoiDiagramCIC;


class ShortestPathFinder2D
{
private:
    //input
    rg_Circle2D       m_StartPt;
    rg_Circle2D       m_EndPt;
    list<rg_Circle2D> m_Obstacles;
    rg_Circle2D       m_Probe;

    double            m_MaxSpeedOfProbe;
    double            m_TravelTime;

    //Output.
    list<EdgeForSolutionPoolGraph*> m_PathByAllPathSearch;
    double                          m_PathTotalDistanceByAllPathSearch;

    list<EdgeForSolutionPoolGraph*> m_PathByBFS;
    double                          m_PathTotalDistanceByBFS;

    list<EdgeForSolutionPoolGraph*> m_PathByDFS;
    double                          m_PathTotalDistanceByDFS;

    //Tool 1. ELLIPSE FILTER & BOX FILTER
    Ellipse2D                        m_EllipseFilter;
    rg_Line2D                        m_LineSegmentForBoxHeight;
    double                           m_HalfLineSegmentForBoxWidth;
    
    unordered_map<EdgeBU2D*, bool>   m_QTEdgeValidation;
    unordered_map<VertexBU2D*, bool> m_QTVertexValidation;


    //Tool 2. QT
    QuasiTriangulation2D             m_QuasiTriangulation;
    BetaUniverse2D                   m_BetaUniverse;
    list<BetaComponent2D>            m_BetaComponents;
    BetaComponentFinder              m_BetaComponentFinder;

    //Tool 3. SPEED UP 
    unordered_set<EdgeBU2D*>                                                              m_QTEdgesOfTopologicalPath;
    unordered_map<pair<VertexBU2D*, VertexBU2D*>, vector<rg_ImplicitEquation>, pair_hash> m_ObstaclePairNItsTangentLines;
 

    //Tool 4. STATISTICS
    TimeStatistics                m_TimeStatistics;

    //Temporary Member Data 
    ShortestPathSolutionPoolGraph   m_SolutionPoolGraphByAllPathSearch;
    list<EdgeBU2D*>                 m_TopologicalPathByAllPathSearch;

    ShortestPathSolutionPoolGraph   m_SolutionPoolGraphByBFS;
    list<EdgeBU2D*>                 m_TopologicalPathByBFS;
    list<VEdge2D*>                  m_TopologicalPathByBFS_InVD;

    ShortestPathSolutionPoolGraph   m_SolutionPoolGraphByDFS;
    list<EdgeBU2D*>                 m_TopologicalPathByDFS;
    list<VEdge2D*>                  m_TopologicalPathByDijkstra_InVD;


    //TOPOLOGICAL PATH TREE
    vector<list<EdgeForSolutionPoolGraph*>> m_PathsByTPT;
    vector<double>                          m_PathTotalDistancesByTPT;
    vector<int>                             m_IndexOfBFSsWithIncreasingOrderOfDistance;

    TopologicalPathTree                   m_TopologicalPathTree;
    vector<ShortestPathSolutionPoolGraph> m_SolutionPoolGraphsByTPT;
    vector<list<EdgeBU2D*>>               m_TopologicalPathsByTPT;

public:
    //constructor
    ShortestPathFinder2D();
    ShortestPathFinder2D(const ShortestPathFinder2D& SPF);

    //destructor;
    ~ShortestPathFinder2D();

    //setter
    inline void set_start_pt(const rg_Circle2D& startPt)  { m_StartPt = startPt; };
    inline void set_end_pt(const rg_Circle2D& endPt)      { m_EndPt = endPt; };
    inline void set_probe(const rg_Circle2D& probe)       { m_Probe   = probe;   };
    inline void set_ellipse(const Ellipse2D& ellipse)     { m_EllipseFilter = ellipse; };
    inline void set_probe_max_speed(const double& probeMaxSpeed)      { m_MaxSpeedOfProbe = probeMaxSpeed; };
    inline void set_max_travel_time(const double& travelTime) { m_TravelTime = travelTime; };

    void set_obstacles(const list<rg_Circle2D>& obstacles);
    
    void set_input_argments(const rg_Circle2D& startPt,
                            const rg_Circle2D& endPt,
                            const rg_Circle2D& probe,
                            const list<rg_Circle2D>& obstacles);

    void set_input_argments(const rg_Circle2D& startPt,
                            const rg_Circle2D& endPt,
                            const rg_Circle2D& probe,
                            const list<rg_Circle2D>& obstacles, 
                            const double& probeMaxSpeed,
                            const double& travelTime);

    //getter
    inline const rg_Circle2D&  get_start_pt() const { return m_StartPt; };
    inline const rg_Circle2D&  get_end_pt()   const { return m_EndPt;   };
    inline const rg_Circle2D&  get_probe()    const { return m_Probe;   };
    inline       Ellipse2D&    get_ellipse()        { return m_EllipseFilter; };
    void                       get_obstacles(list<rg_Circle2D*>& obstacles);
    
    inline void    get_geodesic_path(list<EdgeForSolutionPoolGraph*>& path) { path = m_PathByAllPathSearch; };
    inline void    get_geodesic_path_by_BFS(list<EdgeForSolutionPoolGraph*>& path) { path = m_PathByBFS; };
    inline void    get_geodesic_path_by_DFS(list<EdgeForSolutionPoolGraph*>& path) { path = m_PathByDFS; };
    inline void    get_geodesic_path_by_BFS(list<EdgeForSolutionPoolGraph*>& path, const int& index);
    inline void    get_geodesic_path_by_BFS_with_priority(list<EdgeForSolutionPoolGraph*>& path, const int& priority);

    inline double  get_path_total_distance() const { return m_PathTotalDistanceByAllPathSearch; };
    inline double  get_path_total_distance_by_BFS() const { return m_PathTotalDistanceByBFS; };
    inline double  get_path_total_distance_by_DFS() const { return m_PathTotalDistanceByDFS; };
    inline double  get_path_total_distance_by_BFS(const int& index) const;
    inline double  get_path_total_distance_by_BFS_with_priority(const int& priority) const;

    inline const QuasiTriangulation2D& get_quasi_triangluation() const { return m_QuasiTriangulation;};
    inline const BetaUniverse2D&       get_beta_universe2D()     const { return m_BetaUniverse;};
    inline void get_topological_path(list<EdgeBU2D*>& topologicalPath){ topologicalPath = m_TopologicalPathByAllPathSearch; };
    inline void get_topological_path_by_BFS(list<EdgeBU2D*>& topologicalPath){ topologicalPath = m_TopologicalPathByBFS; };
    inline void get_topological_path_by_DFS(list<EdgeBU2D*>& topologicalPath) { topologicalPath = m_TopologicalPathByDFS; };
    inline void get_topological_path_by_BFS(list<EdgeBU2D*>& topologicalPath, const int& index);
    inline void get_topological_path_by_BFS_with_priority(list<EdgeBU2D*>& topologicalPath, const int& priority);

    inline ShortestPathSolutionPoolGraph& get_solution_pool_graph()        { return m_SolutionPoolGraphByAllPathSearch; };
    inline ShortestPathSolutionPoolGraph& get_solution_pool_graph_by_BFS() { return m_SolutionPoolGraphByBFS; };
    inline ShortestPathSolutionPoolGraph& get_solution_pool_graph_by_DFS() { return m_SolutionPoolGraphByDFS; };
    inline ShortestPathSolutionPoolGraph& get_solution_pool_graph_by_BFS(const int& index);
    inline ShortestPathSolutionPoolGraph& get_solution_pool_graph_by_BFS_with_priority(const int& priority);

    inline int get_num_of_BFS_path() const { return m_PathTotalDistancesByTPT.size(); };
    inline unordered_map<VertexBU2D*, bool>& get_is_this_QT_vertex_inside_of_ellipse() { return m_QTVertexValidation; };

    inline const TimeStatistics&  get_time_statistics()  const { return m_TimeStatistics; };

    //Main function
    bool compute_geodesic_path(const rg_Circle2D& startPt,
                               const rg_Circle2D& endPt,
                               const rg_Circle2D& probe,
                               const list<rg_Circle2D>& obstacles);
    
    void compute_geodesic_path(const rg_Circle2D& startPt,
                               const rg_Circle2D& endPt,
                               const rg_Circle2D& probe,
                               const double& travelTime,
                               const double& probeMaxSpeed,
                               const list<rg_Circle2D>& obstacles);

    void compute_geodesic_path(const rg_Circle2D& startPt,
                               const rg_Circle2D& endPt,
                               const rg_Circle2D& probe,
                               const list<rg_Circle2D>& obstacle,
                               VoronoiDiagram2DC& VD2DC);

    void compute_optimal_path(const rg_Circle2D& startPt,
                               const rg_Circle2D& endPt,
                               const rg_Circle2D& probe,
                               const list<rg_Circle2D>& obstacles);

    bool compute_geodesic_path_based_on_VD( const rg_Circle2D& startPt,
                                            const rg_Circle2D& endPt,
                                            const rg_Circle2D& probe,
                                            const list<rg_Circle2D>& obstacles);
private:
    bool is_there_straight_line_path_between_two_points(const rg_Line2D& straightLinePath) const;
    void construct_voronoi_diagram_of_CIC(const rg_Circle2D& startPt, const rg_Circle2D& endPt, const list<rg_Circle2D>& obstacles, VoronoiDiagramCIC& VD_of_obstacles) const;
    bool find_topogical_path_by_Dijkstra_in_VD(
                const rg_Circle2D& startPt,
                const rg_Circle2D& endPt,
                const rg_Circle2D& probe,
                VoronoiDiagramCIC& VD_of_obstacles,
                list<VEdge2D*>& topologicalPath) const;

    void find_topogical_path_by_breadth_first_search_in_VD(
                const rg_Circle2D& startPt,
                const rg_Circle2D& endPt,
                const rg_Circle2D& probe,
                VoronoiDiagramCIC& VD_of_obstacles,
                list<VEdge2D*>&    topologicalPath);

    bool find_shortest_path_when_only_one_disk_is_intersected_with_straight_line(const rg_Circle2D& obstacle, 
                                                                                 const rg_Circle2D& startPt,
                                                                                 const rg_Circle2D& endPt,
                                                                                 ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                                 list<EdgeForSolutionPoolGraph*>& geodesicPath,
                                                                                 double& totalPathDistance);
        void make_graph_of_solution_pool(const rg_Circle2D& obstacle,
                                         const rg_Circle2D& startPt,
                                         const rg_Circle2D& endPt,
                                         ShortestPathSolutionPoolGraph& solutionPoolGraph) const;

    bool find_shortest_path_from_topological_path_using_VD(
                VoronoiDiagramCIC& VD_of_obstacles,
                const list<VEdge2D*>& topologicalPath,
                ShortestPathSolutionPoolGraph & solutionPoolGraph,
                list<EdgeForSolutionPoolGraph*>& geodesicPath,
                double & totalPathDistacne) const;

        void    make_obstacle_pool(const list<VEdge2D*>& topologicalPath, vector<VFace2D*>& obstaclePool) const;
        void    make_graph_of_solution_pool(
                        VoronoiDiagramCIC&              VD_of_obstacles,
                        const vector<VFace2D*>&         obstaclePool, 
                        ShortestPathSolutionPoolGraph & solutionPoolGraph) const;
        void make_tangent_lines_between_two_obstacles_and_insert_them_into_solution_pool(
                        VoronoiDiagramCIC& VD_of_obstacles,
                        const vector<VFace2D*>& obstaclePool,
                        ShortestPathSolutionPoolGraph& solutionPoolGraph,
                        unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const;
        void make_tangent_lines_between_start_end_points_and_obstacles_and_insert_them_into_solution_pool(
                        VoronoiDiagramCIC& VD_of_obstacles,
                        const vector<VFace2D*>& obstaclePool,
                        ShortestPathSolutionPoolGraph& solutionPoolGraph,
                        unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const;

        void insert_tangent_line_into_graph( 
                        const rg_Point2D& ptOnObstacle1, 
                        const rg_Point2D& ptOnObstacle2, 
                        const rg_ImplicitEquation& tangentLineEq,
                        VertexForSolutionPoolGraph*& vertexOnObstacle1,
                        VertexForSolutionPoolGraph*& vertexOnObstacle2,
                        EdgeForSolutionPoolGraph*& edgeForTangentLine,
                        ShortestPathSolutionPoolGraph& solutionPoolGraph ) const;

        void make_arcs_on_obstacles_and_insert_them_into_solution_pool(
                        const vector<VFace2D*>& obstaclePool,
                        ShortestPathSolutionPoolGraph& solutionPoolGraph,
                        unordered_map<VFace2D*, list<VertexForSolutionPoolGraph*>>& mapObstacleToTangentPoints) const;
        
        void make_sorted_unique_set_of_vertices(const list<VertexForSolutionPoolGraph*>& vertices, 
                                                vector<VertexForSolutionPoolGraph*>& uniqueVertices,
                                                VFace2D* obstacle) const;

        void    insert_arc_into_graph(
                        const rg_Circle2D&              obstacle,
                        VertexForSolutionPoolGraph*     vertexOfStartPt,
                        VertexForSolutionPoolGraph*     vertexOfEndPt,
                        EdgeForSolutionPoolGraph*&      edgeForArc,
                        ShortestPathSolutionPoolGraph&  solutionPoolGraph) const ;

    //Inner function
private:
    //0. PREPROCESS
    void construct_voronoi_diagram(VoronoiDiagram2DC& VD2DC);
    void construct_quasi_triangulation_from_VD(QuasiTriangulation2D& QT2D, VoronoiDiagram2DC& VD2DC);
    void construct_beta_shape_from_QT(BetaUniverse2D& betaUniverse, QuasiTriangulation2D& QT2D);
    void find_connected_components_of_beta_shape(const BetaUniverse2D& betaUniverse, list<BetaComponent2D>& betaComponents, BetaComponentFinder& betaComponentFinder);
    void compute_ellipse_filter(const rg_Circle2D& startPt, const rg_Circle2D& endPt, Ellipse2D& ellipseFilter) const;
    void mark_QT_vertex_whether_in_or_out_of_ellipse(BetaUniverse2D& betaUniverse, Ellipse2D& elipseFilter, unordered_map<VertexBU2D*, bool>& QTVertexValidation);
    void mark_QT_edge_whether_in_or_out_of_ellipse(BetaUniverse2D& betaUniverse, const unordered_map<VertexBU2D*, bool>& QTVertexValidation, unordered_map<EdgeBU2D*, bool>& QTEdgeValidation);
    
    void mark_QT_vertex_whether_in_or_out_of_box(BetaUniverse2D& betaUniverse, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth, unordered_map<VertexBU2D*, bool>& QTVertexValidation);
    void mark_QT_edge_whether_in_or_out_of_box(BetaUniverse2D& betaUniverse, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth, const unordered_map<VertexBU2D*, bool>& QTVertexValidation, unordered_map<EdgeBU2D*, bool>& QTEdgeValidation);
    
    void compute_box_filter(const rg_Circle2D& startPt, const rg_Circle2D& endPt, const double& probeMaxSpeed, const double& travelTime , rg_Line2D& lineSegmentForBoxHeight, double& halfLineSegmentForBoxWidth);

    inline bool is_this_QT_vertex_inside_of_ellipse(VertexBU2D* qtVertex, Ellipse2D& elipseFilter) const;
    inline bool is_this_QT_edge_inside_of_ellipse(EdgeBU2D* qtEdge, const unordered_map<VertexBU2D*, bool>& QTVertexValidation) const;
    
    inline bool is_this_QT_vertex_inside_of_box(VertexBU2D* qtVertex, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth) const;
    inline bool is_this_QT_edge_inside_of_box(EdgeBU2D* qtEdge, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth, const unordered_map<VertexBU2D*, bool>& QTVertexValidation) const;


    //1. COMPUTE GEODESIC
    void compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS();
    void compute_geodesic_path_by_DFS_with_filtered_QT_edges();
    void compute_geodesic_path_by_BFS_with_filtered_QT_edges();
    
    void compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS_ver1();
    void compute_geodesic_path_by_DFS_with_filtered_QT_edges_ver1();
    void compute_geodesic_path_by_BFS_with_filtered_QT_edges_ver1();


    //2. LINE AND OBSTACLE INTESECTION CHECKING
    bool is_this_line_segement_intersected_with_one_of_obstacles(const rg_Line2D& lineSegment, const list<VertexBU2D*>& obstacles) const;
    bool is_this_line_segement_intersected_with_one_of_ALL_obstacles(const rg_Line2D& lineSegment);

    inline bool is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                            const rg_Point2D& startPtOfTangentLine, 
                                                                            const rg_Point2D& endPtOfTangentLine,
                                                                            BetaUniverse2D& betaUniverse) const;

    inline bool is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                            const rg_Point2D& startPtOfTangentLineFromPoint, 
                                                                            const rg_Point2D& endPtOfTangentLineFromObstacle,
                                                                            BetaUniverse2D& betaUniverse,
                                                                            unordered_set<VertexBU2D*>& visitedObstacles)const;
   
    inline bool is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq,
                                                                            const rg_Point2D& startPtOfTangentLineFromPoint, 
                                                                            const rg_Point2D& endPtOfTangentLineFromObstacle,
                                                                            BetaUniverse2D& betaUniverse,
                                                                            VertexBU2D* obstacle,
                                                                            unordered_set<VertexBU2D*>& visitedObstacles)const;

    inline bool is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                            const rg_Point2D& startPtOfTangentLine, 
                                                                            const rg_Point2D& endPtOfTangentLine,
                                                                            BetaUniverse2D& betaUniverse,
                                                                            VertexBU2D* obs1, 
                                                                            VertexBU2D* obs2,
                                                                            unordered_set<VertexBU2D*>& visitedObstacles)const;
    bool is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                            const rg_Point2D& startPtOfTangentLine, 
                                                                            const rg_Point2D& endPtOfTangentLine,
                                                                            BetaUniverse2D& betaUniverse,
                                                                            FaceBU2D* startQTFace, 
                                                                            FaceBU2D* endQTFace,
                                                                            unordered_set<VertexBU2D*>& visitedObstacles)const;


    //3. FILTERING 
    void expand_candidate_QT_edges_by_including_on_a_barrier(const unordered_map<EdgeBU2D*, bool>& QTEdgeValidation,
                                                             unordered_map<EdgeBU2D*, bool>& expandedQTEdgeValidation,
                                                             BetaUniverse2D& betaUniverse);

    void expand_candidate_QT_vertices_by_including_on_a_barrier(const unordered_map<VertexBU2D*, bool>& QTVertexValidation,
                                                             unordered_map<VertexBU2D*, bool>& expandedQTVertexValidation,
                                                             BetaUniverse2D& betaUniverse);




    void compute_geodesic_path_by_BFS(BetaUniverse2D& betaUniverse, 
                                      list<EdgeForSolutionPoolGraph*>& geodesicPath);

    void compute_several_geodesic_paths_by_topological_path_tree(FaceBU2D* qtStartFace,FaceBU2D* qtEndFace, BetaUniverse2D& betaUniverse);

    void compute_several_geodesic_paths_by_topological_path_tree_with_filtered_QT_edges(FaceBU2D* qtStartFace, 
                                                                                        FaceBU2D* qtEndFace, 
                                                                                        const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse,
                                                                                        unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfExpandedEllipse,
                                                                                        BetaUniverse2D& betaUniverse);




    void compute_geodesic_path_by_DFS(BetaUniverse2D& betaUniverse, 
                                      list<EdgeForSolutionPoolGraph*>& geodesicPath);
   
    void compute_geodesic_path_by_all_path_search(BetaUniverse2D& betaUniverse, 
                                                  list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath);

    void compute_geodesic_path_by_all_path_search_with_filtered_QT_edges(BetaUniverse2D& betaUniverse,
                                                                         const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse,
                                                                         list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath);


    void compute_topological_path_tree(FaceBU2D* qtStartFace,
                                        FaceBU2D* qtEndFace,
                                        const BetaUniverse2D& betaUniverse,
                                        TopologicalPathTree& topologicalPathTree) const;

    void compute_topological_path_tree_with_filtered_QT_edges(FaceBU2D* qtStartFace,
                                                              FaceBU2D* qtEndFace,
                                                              const BetaUniverse2D& betaUniverse,
                                                              TopologicalPathTree& topologicalPathTree,
                                                              const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const;


    void find_one_geodesic_path_from_filltered_disks_by_ellipse(const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse,
                                                                BetaUniverse2D& betaUniverse, 
                                                                ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                double& totalPathDistance);

    void find_optimal_path_from_many_topological_ones(list<EdgeBU2D*>& optimalTopologicalPath,
                                                      ShortestPathSolutionPoolGraph& optimalSolutionPoolGraph,
                                                      list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, 
                                                      double& optimalDistance,
                                                      TopologicalPathTree& topologicalPathTree,
                                                      BetaUniverse2D& betaUniverse);

    void convert_QT_face_path_to_QT_edge_path(const list<FaceBU2D*>& qtFacePath, list<EdgeBU2D*>& qtEdgePath) const;
    void convert_QT_edge_path_to_QT_faces(const list<EdgeBU2D*>& qtEdgePath, unordered_set<FaceBU2D*>& qtFacePath) const;
    

    void find_a_geodesic_path_by_stepping_into_the_graph(EdgeBU2D* qtStartEdge, 
                                                            FaceBU2D* qtStartFace, 
                                                            FaceBU2D* qtEndFace, 
                                                            BetaUniverse2D& betaUniverse, 
                                                            list<EdgeBU2D*>& currTopologicalPath, 
                                                            list<EdgeBU2D*>& optimalTopologicalPath, 
                                                            ShortestPathSolutionPoolGraph& optimalSolutionPoolGraph, 
                                                            list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, 
                                                            double& optimalDistance) ;
    void scan_neighbor_QT_face_and_step_into_one_QT_face(FaceBU2D* qtStartFace, 
                                                            FaceBU2D* qtEndFace, 
                                                            BetaUniverse2D& betaUniverse, 
                                                            list<EdgeBU2D*>& currTopologicalPath, 
                                                            list<EdgeBU2D*>& optimalTopologicalPath, 
                                                            ShortestPathSolutionPoolGraph& optimalSolutionPoolGraph, 
                                                            list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath, 
                                                            double& optimalDistance) ;

    void find_one_geodesic_path_when_both_of_QT_faces_from_source_N_destination_are_same(FaceBU2D* qtFace,
                                                                                         BetaUniverse2D& betaUniverse, ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                                         list<EdgeForSolutionPoolGraph*>& geodesicPath, double& totalPathDistacne);
     void find_one_geodesic_path_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(const rg_ImplicitEquation& lineEqFromSourceToDestination,
                                                                                                                        ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                                                                        list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath,
                                                                                                                        double& totalPathDistance) const;

     rg_Circle2D find_next_stop_point_when_both_of_QT_faces_from_source_N_destination_are_same(FaceBU2D* qtFace,
                                                                                              BetaUniverse2D& betaUniverse, 
                                                                                              ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                                              list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                                              double& totalPathDistacne);


       
     rg_Circle2D find_next_stop_point_when_straight_line_from_source_to_destination_having_no_intersection_with_any_obstacles(const rg_ImplicitEquation& lineEqFromSourceToDestination,
                                                                                                                             ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                                                                                             list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath,
                                                                                                                             double& totalPathDistance) const;
     rg_Circle2D find_next_stop_point_in_general_case(BetaUniverse2D& betaUniverse, 
                                                     ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                     list<EdgeForSolutionPoolGraph*>& optimalGeodesicPath,
                                                     double& totalPathDistance) ;


    void find_one_topological_path_from_QT_by_BFS(FaceBU2D* qtStartFace, 
                                                  FaceBU2D* qtEndFace, 
                                                  list<EdgeBU2D*>& topologicalPath, 
                                                  BetaUniverse2D& betaUniverse) const;

    void find_several_topological_path_from_QT_by_topological_path_tree(FaceBU2D* qtStartFace, 
                                                                        FaceBU2D* qtEndFace, 
                                                                        vector<list<EdgeBU2D*>>& topologicalPaths, 
                                                                        TopologicalPathTree& topologicalPathTree,
                                                                        BetaUniverse2D& betaUniverse) const;

    void find_one_topological_path_from_QT_by_BFS_with_filtered_QT_edges(FaceBU2D* qtStartFace,
                                                                         FaceBU2D* qtEndFace, 
                                                                         list<EdgeBU2D*>& topologicalPath,
                                                                         const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse,
                                                                         BetaUniverse2D& betaUniverse) const;

    void find_one_topological_path_from_QT_by_DFS(FaceBU2D* qtStartFace, FaceBU2D* qtEndFace, list<EdgeBU2D*>& topologicalPath, BetaUniverse2D& betaUniverse) const;
    
    void find_one_topological_path_from_QT_by_DFS_with_filtered_QT_edges(FaceBU2D* qtStartFace,
                                                                         FaceBU2D* qtEndFace, 
                                                                         list<EdgeBU2D*>& topologicalPath,
                                                                         const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse,
                                                                         BetaUniverse2D& betaUniverse) const;

    void initialize_visitation_tag_of_QT_face(BetaUniverse2D& betaUniverse) const;
    void insert_three_bounding_QT_edge_into_container(FaceBU2D* qtFace, list <EdgeBU2D*>& stackOfQTEdges) const;
    void insert_three_bounding_QT_edge_into_container(FaceBU2D* qtFace, EntityAccessiblePriorityQ<EdgeBU2D*>& priorityQ) const;
    void insert_QT_edge_into_PQ(EdgeBU2D* qtEdge, EntityAccessiblePriorityQ<EdgeBU2D*>& priorityQ) const;

    void compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_BFS(FaceBU2D* qtStartFace, 
                                                                       FaceBU2D* qtEndFace, 
                                                                       EdgeBU2D*& gatewayQtEdgeOfEndFace,
                                                                       unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const;
    
    void compute_gateway_QT_edges_of_end_face_and_its_prev_edges_by_BFS(FaceBU2D* qtStartFace, 
                                                                        FaceBU2D* qtEndFace, 
                                                                        vector<EdgeBU2D*>& gatewayQtEdgesOfEndFace,
                                                                        unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const;

    void compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_BFS_with_filtered_QT_edges(FaceBU2D* qtStartFace,
                                                                                              FaceBU2D* qtEndFace, 
                                                                                              EdgeBU2D*& gatewayQtEdgeOfEndFace,
                                                                                              unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper,
                                                                                              const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const;
    void compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_adjusted_BFS_with_filtered_QT_edges(FaceBU2D* qtStartFace,
                                                                                                       FaceBU2D* qtEndFace, 
                                                                                                       EdgeBU2D*& gatewayQtEdgeOfEndFace,
                                                                                                       unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper,
                                                                                                       const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const;
    
    void compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_DFS(FaceBU2D* qtStartFace, 
                                                                       FaceBU2D* qtEndFace, 
                                                                       EdgeBU2D*& gatewayQtEdgeOfEndFace,
                                                                       unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const;
    
    void compute_gateway_QT_edge_of_end_face_and_its_prev_edges_by_DFS_with_filtered_QT_edges(FaceBU2D* qtStartFace,
                                                                                              FaceBU2D* qtEndFace, 
                                                                                              EdgeBU2D*& gatewayQtEdgeOfEndFace,
                                                                                              unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper,
                                                                                              const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const;


    void compute_gateway_QT_edge_of_end_face_and_its_prev_edges(FaceBU2D* qtStartFace, 
                                                                FaceBU2D* qtEndFace, 
                                                                EdgeBU2D*& gatewayQtEdgeOfEndFace, 
                                                                unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper) const;
    
    EdgeBU2D* get_sharing_QT_edge_of_two_QT_faces(FaceBU2D* qtFace1, FaceBU2D* qtFace2) const;
    
    void back_trace_to_QT_start_face_and_find_topological_path(EdgeBU2D* gatewayQtEdgeOfEndFace, 
                                                               unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper,
                                                               list<EdgeBU2D*>& topologicalPath) const;
    
    void back_trace_to_QT_start_face_and_find_topological_paths(vector<EdgeBU2D*>& gatewayQtEdgesOfEndFace, 
                                                                unordered_map<EdgeBU2D*, EdgeBU2D*>& currToPrevQTEdgeMapper,
                                                                vector<list<EdgeBU2D*>>& topologicalPaths) const;

    void find_one_geodesic_path_from_QT_topological_path(const list<EdgeBU2D*>& topologicalPath, 
                                                         BetaUniverse2D& betaUniverse, 
                                                         ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                         list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                         double& totalPathDistacne);

    void find_one_geodesic_path_from_QT_topological_path_for_all_path_search(const list<EdgeBU2D*>& topologicalPath, 
                                                                             BetaUniverse2D& betaUniverse, 
                                                                             ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                             list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                                             double& totalPathDistacne);

    void find_one_geodesic_path_from_candidate_obstacles(const unordered_map<VertexBU2D*, bool>& isThisQTVertexInsideOfEllipse,
                                                         BetaUniverse2D& betaUniverse, 
                                                         ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                         list<EdgeForSolutionPoolGraph*>& geodesicPath, 
                                                         double& totalPathDistacne);

    void find_several_geodesic_paths_from_QT_topological_paths(const vector<list<EdgeBU2D*>>& topologicalPaths, 
                                                                BetaUniverse2D& betaUniverse, 
                                                                vector<ShortestPathSolutionPoolGraph>& solutionPoolGraphs, 
                                                                vector<list<EdgeForSolutionPoolGraph*>>& geodesicPaths, 
                                                                vector<double>& totalPathDistances);

    void get_all_QT_vertices_from_topological_path(set<VertexBU2D*>& qtVertice, const list<EdgeBU2D*>& qtEdge) const;
    void generate_solution_pool_graph(const vector<VertexBU2D*>& obstaclePool, BetaUniverse2D& betaUniverse, ShortestPathSolutionPoolGraph& solutionPoolGraph);
    void generate_solution_pool_graph_with_topological_path(const vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D& betaUniverse,  ShortestPathSolutionPoolGraph& solutionPoolGraph);

    void generate_solution_pool_graph(const vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, BetaUniverse2D& betaUniverse, ShortestPathSolutionPoolGraph& solutionPoolGraph);
    void generate_visible_solution_pool_graph_from_current_position(const rg_Point2D& currPosition, const vector<VertexBU2D*>& obstaclePool, BetaUniverse2D& betaUniverse, ShortestPathSolutionPoolGraph& solutionPoolGraph) const;
  


    void make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, const unordered_map<VertexBU2D*, bool>& isThisQTVertexInsideOfEllipse) const;
    void make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse) const;
    void make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& qtEdge) const;
    void make_obstacle_pool(vector<VertexBU2D*>& obstaclePool, FaceBU2D* qtFace) const;
    void make_obstacle_pool(vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, const list<EdgeBU2D*>& qtEdge) const;

    void convert_filltered_QT_edges_to_topological_path(const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse, 
                                                        list<EdgeBU2D*>& topologicalPath) const;

    void make_mapper_from_QT_obstacle_pair_to_its_tangent_lines(VertexBU2D* obstacle1, VertexBU2D* obstacle2, const vector<rg_ImplicitEquation>& tangentLines,
                                                                unordered_map<pair<VertexBU2D*, VertexBU2D*>, vector<rg_ImplicitEquation>, pair_hash>& obstaclePairNItsTangentLines) const;
    bool are_the_tangent_lines_of_this_obstacle_pair_already_computed(VertexBU2D* obstacle1, VertexBU2D* obstacle2) const;
    void get_tangent_lines_of_this_obstacle_pair_from_already_computed(VertexBU2D* obstacle1, VertexBU2D* obstacle2, vector<rg_ImplicitEquation>& tangentLines) const;
            
    void make_set_of_obstacles(const vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs, vector<VertexBU2D*>& obstaclePool) const;
    void generate_graph_for_obstacles_by_using_QT_and_BS(const vector<VertexBU2D*>& obstaclePool,
                                                            BetaUniverse2D& betaUniverse, 
                                                            ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) ;
    void generate_graph_for_obstacles_by_using_QT_and_BS_with_topological_path(const vector<VertexBU2D*>& obstaclePool,
                                                                               const list<EdgeBU2D*>& topologicalPath,
                                                                               BetaUniverse2D& betaUniverse, 
                                                                               ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) ;



    void generate_graph_for_obstacles_by_using_QT_and_BS(const vector<pair<VertexBU2D*, VertexBU2D*>>& obstaclePairs,
                                                            BetaUniverse2D& betaUniverse, 
                                                            ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) ;

    inline void generate_graph_incrementally_by_adding_two_obstacles_into_the_graph(VertexBU2D* obs1,
                                                                                    VertexBU2D* obs2,
                                                                                    const BetaComponentFinder& betaComponentFinder,
                                                                                    BetaUniverse2D& betaUniverse, 
                                                                                    ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) ;
    inline void generate_graph_incrementally_by_adding_two_obstacles_into_the_graph_with_topological_path(VertexBU2D* obs1,
                                                                                                          VertexBU2D* obs2,
                                                                                                          const list<EdgeBU2D*>& topologicalPath,
                                                                                                          const BetaComponentFinder& betaComponentFinder,
                                                                                                          BetaUniverse2D& betaUniverse, 
                                                                                                          ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) ;

    void make_graph_of_tangent_lines(VertexBU2D* obs1, 
                                     VertexBU2D* obs2, 
                                     const BetaComponentFinder& betaComponentFinder,
                                     BetaUniverse2D& betaUniverse,
                                     ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                     vector<VertexForSolutionPoolGraph*>& verticesOnObstacle1,
                                     vector<VertexForSolutionPoolGraph*>& verticesOnObstacle2);

    void make_graph_of_tangent_lines_with_topological_path(VertexBU2D* obs1, 
                                                           VertexBU2D* obs2, 
                                                           const list<EdgeBU2D*>& topologicalPath,
                                                           const BetaComponentFinder& betaComponentFinder,
                                                           BetaUniverse2D& betaUniverse,
                                                           ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                           vector<VertexForSolutionPoolGraph*>& verticesOnObstacle1,
                                                           vector<VertexForSolutionPoolGraph*>& verticesOnObstacle2);
    inline void find_tangent_line_equations_of_two_obstacles(VertexBU2D* obs1, VertexBU2D* obs2, 
                                                                const BetaComponentFinder& betaComponentFinder, 
                                                                vector<rg_ImplicitEquation>& tangentLineEqs);
    rg_ImplicitEquation make_line_equation(const rg_Point2D& pt1, const rg_Point2D& pt2) const;
    bool are_these_two_lines_intersected_in_this_interval(const rg_ImplicitEquation& lineEq1, const rg_ImplicitEquation& lineEq2, double x_s, double x_e) const;
    bool are_these_two_obstacles_in_same_component(VertexBU2D* obstacle1, VertexBU2D* obstacle2, const BetaComponentFinder& betaComponentFinder) const;
    void find_tangent_line_equations_of_two_obstacles_in_same_component(const rg_Circle2D& obstacleDisk1, const rg_Circle2D& obstacleDisk2, vector<rg_ImplicitEquation>& tangentLineEqs) const;
    void find_tangent_line_equations_of_two_obstacles_in_diff_component(const rg_Circle2D& obstacleDisk1, const rg_Circle2D& obstacleDisk2, vector<rg_ImplicitEquation>& tangentLineEqs) const;
  
    void filter_tangent_line_equations_and_points_by_disks(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                           vector<rg_ImplicitEquation>& filletedTangentLineEqs,
                                                           vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2,
                                                           BetaUniverse2D& betaUniverse,
                                                           VertexBU2D* obs1, VertexBU2D* obs2) const;
    void filter_tangent_line_equations_and_points_by_disks_vers1(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                 vector<rg_ImplicitEquation>& filletedTangentLineEqs,
                                                                 vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2,
                                                                 BetaUniverse2D& betaUniverse,
                                                                 VertexBU2D* obs1, VertexBU2D* obs2) const;
    void filter_tangent_line_equations_and_points_by_disks(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                           vector<rg_ImplicitEquation>& filletedTangentLineEqs, 
                                                           vector<rg_Point2D>& listOfPtsOnObstacle,
                                                           BetaUniverse2D& betaUniverse, 
                                                           VertexForSolutionPoolGraph* pointVtx, VertexBU2D* obstacle) const;
    void filter_tangent_line_equations_and_points_by_disks_vers1(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                 vector<rg_ImplicitEquation>& filletedTangentLineEqs, 
                                                                 vector<rg_Point2D>& listOfPtsOnObstacle,
                                                                 BetaUniverse2D& betaUniverse, 
                                                                 VertexForSolutionPoolGraph* pointVtx, VertexBU2D* obstacle) const;

    void filter_tangent_line_equations_and_points_by_topological_path(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                      const vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2,
                                                                      vector<rg_ImplicitEquation>& filletedTangentLineEqs,
                                                                      vector<pair<rg_Point2D, rg_Point2D>>& filletedListOfTwoPtsOnObstacle1And2,
                                                                      BetaUniverse2D& betaUniverse,
                                                                      const list<EdgeBU2D*>& topologicalPath) const;
    void filter_tangent_line_equations_and_points_by_topological_path(const vector<rg_ImplicitEquation>& tangentLineEqs, 
                                                                      const vector<rg_Point2D>& listOfPtsOnObstacle,
                                                                      vector<rg_ImplicitEquation>& filletedTangentLineEqs,
                                                                      vector<rg_Point2D>& filletedListOfPtsOnObstacle,
                                                                      VertexForSolutionPoolGraph* pointVtx,
                                                                      BetaUniverse2D& betaUniverse,
                                                                      const list<EdgeBU2D*>& topologicalPath) const;
    
    inline bool are_this_line_intersected_with_one_of_three_disks_of_QT_face(const rg_ImplicitEquation& tangentLineEq, 
                                                                      FaceBU2D* qtFace, 
                                                                      const rg_Point2D& startPt, 
                                                                      const rg_Point2D& endPt,
                                                                      unordered_set<VertexBU2D*>& visitedObstacles) const;
    
    bool are_this_line_intersected_with_these_disks(const rg_ImplicitEquation& tangentLineEq, 
                                                    rg_dList<VertexBU2D*>& obstacles, 
                                                    const rg_Point2D& startPt, 
                                                    const rg_Point2D& endPt,
                                                    unordered_set<VertexBU2D*>& visitedObstacles) const;
    bool check_intersection_between_line_segment_N_all_disks(const rg_ImplicitEquation& lineEq, 
                                                             const rg_Point2D& startPtOfTangentLine,
                                                             const rg_Point2D& endPtOfTangentLine, 
                                                             BetaUniverse2D& betaUniverse, 
                                                             unordered_set<VertexBU2D*>& visitedObstacles) const;
    inline bool is_this_val_in_this_interval(double val, double startInterval, double endInterval) const { return (((val >= startInterval) && (val <= endInterval)) 
                                                                                                                || ((val <= startInterval) && (val >= endInterval))) ; };
    void propagate_to_one_of_neighbor_QT_faces_passed_through_by_tangent_line(FaceBU2D*& currQTFace, 
                                                                            FaceBU2D*& prevQTFace,
                                                                            const rg_Point2D& startPtOfTangentLine,
                                                                            const rg_Point2D& endPtOfTangentLine,
                                                                            const rg_ImplicitEquation& tangentLineEq) const;
    void find_one_general_QT_edge_bounding_QT_Face(EdgeBU2D*& generalQTEdge, FaceBU2D* qtFace) const;
    bool are_tangent_line_segment_N_QT_edge_intersected(const rg_ImplicitEquation& tangentLineEq, 
                                                        const rg_Point2D& startPt,
                                                        const rg_Point2D& endPt, 
                                                        EdgeBU2D* qtEdge) const;
     bool are_tangent_line_segment_intersected_with_one_of_QT_edges(const rg_ImplicitEquation& tangentLineEq, 
                                                                   const rg_Point2D& startPt,
                                                                   const rg_Point2D& endPt, 
                                                                   const list<EdgeBU2D*>& qtEdges) const;

    bool are_tangent_line_segment_intersected_with_one_of_QT_edges(const rg_ImplicitEquation& tangentLineEq, 
                                                                   const rg_Point2D& startPt,
                                                                   const rg_Point2D& endPt, 
                                                                   const unordered_set<EdgeBU2D*>& qtEdges) const;

    inline void set_prev_and_curr_QT_face_for_loop(FaceBU2D*& currQTFace, FaceBU2D*& prevQTFace, EdgeBU2D* currQTEdge) const;
    void set_prev_and_curr_virtual_QT_face_for_loop(FaceBU2D*& currQTFace, FaceBU2D*& prevQTFace, EdgeBU2D* currQTEdge, const rg_Point2D& endPt) const;
    void make_vertices_N_tangent_line_edges(const vector<rg_ImplicitEquation>& filteredTangentLineEqs,
                                            const vector<pair<rg_Point2D, rg_Point2D>>& listOfTwoPtsOnObstacle1And2,
                                            vector<VertexForSolutionPoolGraph*>& verticesOnObstacle1,
                                            vector<VertexForSolutionPoolGraph*>& verticesOnObstacle2,
                                            ShortestPathSolutionPoolGraph& solutionPoolGraph) const;
    inline void create_two_vertices_and_a_tangent_line_edge(VertexForSolutionPoolGraph*& vertex1, 
                                                            VertexForSolutionPoolGraph*& vertex2, 
                                                            EdgeForSolutionPoolGraph*& edge, 
                                                            ShortestPathSolutionPoolGraph& solutionPoolGraph) const;
    inline void set_two_vertices_and_a_tangent_line_edge(VertexForSolutionPoolGraph* vertex1,
                                                        VertexForSolutionPoolGraph* vertex2,
                                                        EdgeForSolutionPoolGraph* edge,
                                                        const pair<rg_Point2D, rg_Point2D>& twoPtsOnObstacle1And2,
                                                        const rg_ImplicitEquation& implicitEquationOfTangentLine) const;
    void add_new_arc_edges_made_by_new_vertices_into_the_graph(VertexBU2D* obstacle, 
                                                                const vector<VertexForSolutionPoolGraph*>& newVerticesOnObstacle,
                                                                const BetaComponentFinder& betaComponentFinder,
                                                                ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;
    void insert_into_list_of_tangent_vertices_on_disk_by_angle(VertexBU2D* obstacle, VertexForSolutionPoolGraph* newVertexOnObstacle, list<VertexForSolutionPoolGraph*>& verticesOnObstacle) const;
    inline bool does_this_obstacle_have_already_maden_vetices_of_tangent_points(VertexBU2D* obstacle, const unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;
    inline void insert_new_vertex_into_set_of_vertices_on_obstacle(VertexForSolutionPoolGraph* newVertexOnObstacle, VertexBU2D* obstacle, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;

    void filter_vertices_for_making_arc_edges(VertexBU2D* obstacle, 
                                              VertexForSolutionPoolGraph* newVertexOnObstacle,
                                              const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex,
                                              const BetaComponentFinder& betaComponentFinder,
                                              list<pair<VertexForSolutionPoolGraph*,VertexForSolutionPoolGraph*>>& vertexPairsForArcEdges) const;
    void filter_vertices_for_making_arc_edges_when_NONE_verices_on_obstacle(VertexBU2D* obstacle,
                                                                            VertexForSolutionPoolGraph* newVertexOnObstacle,
                                                                            const BetaComponentFinder& betaComponentFinder,
                                                                            list<pair<VertexForSolutionPoolGraph*,VertexForSolutionPoolGraph*>>& vertexPairsForArcEdges) const;
    
    void find_nearest_two_vertices_including_this_new_vertex(VertexForSolutionPoolGraph*& prevVertex, 
                                                                VertexForSolutionPoolGraph*& nextVertex, 
                                                                VertexForSolutionPoolGraph* newVertexOnObstacle,
                                                                const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex,
                                                                VertexBU2D* obstacle) const;
                    
    void find_prev_N_next_vertices_of_this_new_vertex_without_intersection_with_obstacles(list<VertexForSolutionPoolGraph*>& prevVertices,
                                                                                            list<VertexForSolutionPoolGraph*>& nextVertices, 
                                                                                            VertexForSolutionPoolGraph* newVertexOnObstacle,
                                                                                            VertexBU2D* obstacle,
                                                                                            const list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex,
                                                                                            const list<VertexBU2D*>& neighborIntersectedObstacles) const;
    void find_nearest_two_vertices_including_this_new_vertex(VertexForSolutionPoolGraph*& prevVertex, 
                                                                VertexForSolutionPoolGraph*& nextVertex,
                                                                const list<VertexForSolutionPoolGraph*>& prevVertices,
                                                                const list<VertexForSolutionPoolGraph*>& nextVertices, 
                                                                VertexForSolutionPoolGraph* newVertexOnObstacle,
                                                                VertexBU2D* obstacle) const;
    void get_neighbor_obstacles_interseted_with_this_obstacle(list<VertexBU2D*>& neighborIntersectedObstacles, 
                                                                VertexBU2D* obstacle, 
                                                                const BetaComponentFinder& betaComponentFinder) const;
    void get_angle_of_neighbor_intersected_obstacles(list<double>& anglesOfNeighborIntersectedObstacles, 
                                                     const list<VertexBU2D*>& neighborIntersectedObstacles,
                                                     VertexBU2D* obstacle) const;

    void get_angle_of_neighbor_qt_edges_NOT_in_topological_path(list<double>& anglesOfNeighborQTEdges, 
                                                            VertexBU2D* obstacle) const;

    double get_angle_of_this_vertex_with_regard_to_this_obstacle(VertexForSolutionPoolGraph* vertexOnObstacle, VertexBU2D* obstacle) const;

    inline double get_angle_of_start_to_end_point(const rg_Point2D& startPt, const rg_Point2D& endPt) const;

    inline bool do_any_of_two_intervals_include_none_of_obstacles(const double& angle1, 
                                                                  const double& angle2, 
                                                                  const list<double>& anglesOfNeighborIntersectedObstacles, int& intervalOption) const;
   
    inline bool do_any_of_two_intervals_include_none_of_obstacles(const double& angle1, 
                                                                  const double& angle2, 
                                                                  const list<double>& anglesOfNeighborIntersectedObstacles) const;
   
    inline bool do_any_of_QT_faces_include_this_pt(const rg_Point2D& pt, const unordered_set<FaceBU2D*>& qtFaces, BetaUniverse2D& betaUniverse) const;


    bool do_this_interval_include_none_of_obstacles(const double& startAngle,
                                                    const double& endAngle, 
                                                    const list<double>& anglesOfNeighborIntersectedObstacles) const;

    void make_arc_edges(const list<pair<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>>& vertexPairsForArcEdges,
                        VertexBU2D* obstacle,
                        ShortestPathSolutionPoolGraph& solutionPoolGraph) const;

    inline void create_an_arc_edge(EdgeForSolutionPoolGraph*& edge, ShortestPathSolutionPoolGraph& solutionPoolGraph) const;
    inline void set_an_arc_edge(EdgeForSolutionPoolGraph* arcEdge,
                                VertexForSolutionPoolGraph* startVertex, VertexForSolutionPoolGraph* endVertex, 
                                VertexBU2D* obstacle, 
                                ShortestPathSolutionPoolGraph& solutionPoolGraph) const;

    inline void set_two_arc_edges(EdgeForSolutionPoolGraph* newArcEdge,
                                  VertexForSolutionPoolGraph* prevVertex, VertexForSolutionPoolGraph* newVertex, VertexForSolutionPoolGraph* nextVertex,
                                  VertexBU2D* obstacle, 
                                  ShortestPathSolutionPoolGraph& solutionPoolGraph) const;

    EdgeForSolutionPoolGraph* find_arc_edge_having_these_two_vertices(VertexForSolutionPoolGraph* vertex1, VertexForSolutionPoolGraph* vertex2) const;

    rg_ImplicitEquation make_implicit_equation_of_arc_edge(VertexBU2D* obstacle, const double& probeRadius) const;
    rg_ImplicitEquation make_implicit_equation_of_arc_edge(const rg_Circle2D& obstacle, const double& probeRadius) const;
    double compute_arc_segment_length_on_obstacle_boundary(VertexBU2D* obstacle, const rg_Point2D& arcStartPt, const rg_Point2D& arcEndPt) const;
    double compute_arc_segment_length_on_obstacle_boundary(const rg_Circle2D& obstacle, const rg_Point2D& arcStartPt, const rg_Point2D& arcEndPt) const;

    void insert_start_end_point_into_the_graph(const rg_Circle2D& startPt,
                                                const rg_Circle2D& endPt,
                                                const vector<VertexBU2D*>& obstaclePool,
                                                BetaUniverse2D& betaUniverse,
                                                ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;

     void insert_start_end_point_into_the_graph_with_topological_path(const rg_Circle2D& startPt, 
                                                                      const rg_Circle2D& endPt,
                                                                      const vector<VertexBU2D*>& obstaclePool,
                                                                      const list<EdgeBU2D*>& topologicalPath,
                                                                      BetaUniverse2D& betaUniverse,
                                                                      ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;


    inline void generate_graph_incrementally_by_adding_one_obstacle_N_point_into_the_graph(VertexBU2D* obstacle,
                                                                                           VertexForSolutionPoolGraph* pointVtx,
                                                                                           const BetaComponentFinder& betaComponentFinder,
                                                                                           BetaUniverse2D& betaUniverse, 
                                                                                           ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const ;
    inline void generate_graph_incrementally_by_adding_one_obstacle_N_point_into_the_graph_with_topological_path(VertexBU2D* obstacle,
                                                                                                                 VertexForSolutionPoolGraph* pointVtx,
                                                                                                                 const list<EdgeBU2D*>& topologicalPath,
                                                                                                                 const BetaComponentFinder& betaComponentFinder,
                                                                                                                 BetaUniverse2D& betaUniverse, 
                                                                                                                 ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const ;

    void make_graph_of_tangent_lines(VertexBU2D* obstacle, 
                                    VertexForSolutionPoolGraph* pointVtx, 
                                    BetaUniverse2D& betaUniverse,
                                    ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                    vector<VertexForSolutionPoolGraph*>& verticesOnObstacle) const;

    void make_graph_of_tangent_lines_with_topological_path(VertexBU2D* obstacle, 
                                                           VertexForSolutionPoolGraph* pointVtx, 
                                                           const list<EdgeBU2D*>& topologicalPath,
                                                           BetaUniverse2D& betaUniverse,
                                                           ShortestPathSolutionPoolGraph& solutionPoolGraph,
                                                           vector<VertexForSolutionPoolGraph*>& verticesOnObstacle) const;


    void insert_a_point_into_the_graph(VertexForSolutionPoolGraph* pointVtx,
                                        const vector<VertexBU2D*>& obstaclePool,
                                        BetaUniverse2D& betaUniverse,
                                        ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;
    void insert_a_point_into_the_graph_with_topological_path(VertexForSolutionPoolGraph* pointVtx,
                                                              const vector<VertexBU2D*>& obstaclePool,
                                                              const list<EdgeBU2D*>& topologicalPath,
                                                              BetaUniverse2D& betaUniverse,
                                                              ShortestPathSolutionPoolGraph& solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const;
    inline void find_tangent_line_equations_of_point_and_obstacle(VertexForSolutionPoolGraph* pointVtx,
                                                                  VertexBU2D* obstacle,
                                                                  vector<rg_ImplicitEquation>& tangentLineEqs) const;

    void make_vertices_N_tangent_line_edges(const vector<rg_ImplicitEquation>& filteredTangentLineEqs,
                                            const vector<rg_Point2D>& listOfPtsOnObstacle,
                                            VertexForSolutionPoolGraph* pointVtx,
                                            vector<VertexForSolutionPoolGraph*>& verticesOnObstacle,
                                            ShortestPathSolutionPoolGraph& solutionPoolGraph) const;
    void set_two_vertices_and_a_tangent_line_edge(VertexForSolutionPoolGraph* vtxOfPoint,
                                                    VertexForSolutionPoolGraph* vtxOnObstacle,
                                                    EdgeForSolutionPoolGraph* edge,
                                                    const rg_Point2D& point,
                                                    const rg_Point2D& ptOnObstacle,
                                                    const rg_ImplicitEquation& implicitEquationOfTangentLine) const;
        
    void find_geodesic_path_by_dijkstra(const ShortestPathSolutionPoolGraph& solutionPoolGraph) const;
    bool back_trace_to_start_vtx_and_find_geodesic_path(VertexForSolutionPoolGraph* startVtx, VertexForSolutionPoolGraph* endVtx, list<EdgeForSolutionPoolGraph*>& geodesicPath, double& totalPathDistance) const;

    inline bool is_this_QT_edge_in_an_exterior_state(EdgeBU2D* qtEdge) const;
    inline bool is_this_QT_face_in_an_exterior_state(FaceBU2D* qtFace) const;
    bool are_some_bounding_edges_of_this_QT_face_exterior(FaceBU2D* qtFace) const;
    bool is_this_QT_face_visited(FaceBU2D* qtFace) const;
     
    FaceBU2D* get_opposite_QT_face_of_Edge(EdgeBU2D* qtEdge, FaceBU2D* qtFace) const;

    void change_optimal_path_in_old_graph_to_new_one(const ShortestPathSolutionPoolGraph& oldGraph, 
                                                    const ShortestPathSolutionPoolGraph& newGraph,
                                                    list<EdgeForSolutionPoolGraph*>& geodesicPath) const;
    static bool compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order (const VertexBU2D* obstacle1, const VertexBU2D* obstacle2);

    //Query for Shortest path problem
    inline bool is_topological_path_by_BFS_computed();

private:
    void copy_from(const ShortestPathFinder2D& SPF);
    

    static void make_exterior_tangent_lines_of_two_circles(  const rg_Circle2D& circle1,
											                 const rg_Circle2D& circle2,
											                 rg_ImplicitEquation& result1,
											                 rg_ImplicitEquation& result2);

    static void make_interior_tangent_lines_of_two_circles(  const rg_Circle2D& circle1,
											                 const rg_Circle2D& circle2,
											                 rg_ImplicitEquation& result1,
											                 rg_ImplicitEquation& result2);

    static rg_Point2D compute_tangent_point_between_line_and_circle(const rg_Circle2D& circle,
                                                                    const rg_ImplicitEquation& line);

    //For BFS
    void sort_index_of_distances_in_non_decreasing_order(const vector<double>& distances, vector<int>& indexOfDistancesWithIncreasingOrder) const;
    static bool compare_pair_of_distance_N_index_in_non_decreasing_order(const pair<double, int>& pair1, const pair<double, int>& pair2);
    
    ShortestPathSolutionPoolGraph* create_solution_pool_of_BFS(const ShortestPathSolutionPoolGraph& solutionPoolGraph, const int& index);

    //For statistics
    void initialize_statistics();
    void finalize_statistics();

    inline void ___start_clock(const unsigned int& it){ m_TimeStatistics.start_clock(it); };
    inline void ___end_clock(const unsigned int& it)  { m_TimeStatistics.end_clock(it); };




/******************************************************************************************************/
/******************************************************************************************************/
/*                                                                                                    */
/*                                        TEMP. CODE                                                  */
/*                                                                                                    */
/******************************************************************************************************/
/******************************************************************************************************/
 private:  
     
    rg_Circle2D find_next_stop_point_by_greedy_method(const rg_Circle2D& targetPt,
                                                      const rg_Circle2D& currProbe,
                                                      const double& speed,
                                                      const double& travelTime,
                                                      const list<rg_Circle2D>& obstacles);
    
    rg_Circle2D find_next_stop_point(BetaUniverse2D& betaUniverse, list<EdgeForSolutionPoolGraph*>& geodesicPath);

    void compute_geodesic_path_by_all_path_search_with_filtered_QT_edges_N_BFS(BetaUniverse2D& betaUniverse,
                                                                               const unordered_map<VertexBU2D*, bool>& isThisQTVertexInsideOfEllipse,
                                                                               const unordered_map<EdgeBU2D*, bool>& isThisQTEdgeInsideOfEllipse);

    void mark_QT_edge_whether_in_or_out_of_ellipse(BetaUniverse2D& betaUniverse, Ellipse2D& elipseFilter, unordered_map<EdgeBU2D*, bool>& QTEdgeValidation);

};


inline bool ShortestPathFinder2D::is_this_QT_vertex_inside_of_ellipse(VertexBU2D* qtVertex, 
                                                                      Ellipse2D& elipseFilter) const
{
    return elipseFilter.does_contain(qtVertex->getCoord());
}



inline bool ShortestPathFinder2D::is_this_QT_edge_inside_of_ellipse(EdgeBU2D* qtEdge,
                                                                    const unordered_map<VertexBU2D*, bool>& QTVertexValidation) const
{
    bool isStartVtxInsideOfEllipse = QTVertexValidation.at(qtEdge->getStartVertex());
    bool isEndVtxInsideOfEllipse = QTVertexValidation.at(qtEdge->getEndVertex());

    if (isStartVtxInsideOfEllipse || isEndVtxInsideOfEllipse)
    //if (isStartVtxInsideOfEllipse && isEndVtxInsideOfEllipse)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool ShortestPathFinder2D::is_this_QT_vertex_inside_of_box(VertexBU2D* qtVertex, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth) const
{
    rg_Point2D targetPt(qtVertex->getCoord());

    double parameter = 0.0;
    lineSegmentForBoxHeight.project(targetPt, parameter);

    if (parameter < 0)
    {
        return false;
    }
    else if (parameter > 1)
    {
        return false;
    }
    else //(parameter >=0 && parameter <= 1)
    {
        if (lineSegmentForBoxHeight.getDistance(targetPt) < halfLineSegmentForBoxWidth)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


inline bool ShortestPathFinder2D::is_this_QT_edge_inside_of_box(EdgeBU2D* qtEdge, const rg_Line2D & lineSegmentForBoxHeight, const double& halfLineSegmentForBoxWidth, const unordered_map<VertexBU2D*, bool>& QTVertexValidation) const
{
    bool isStartVtxInsideOfBox = QTVertexValidation.at(qtEdge->getStartVertex());
    bool isEndVtxInsideOfBox   = QTVertexValidation.at(qtEdge->getEndVertex());

    if (isStartVtxInsideOfBox || isEndVtxInsideOfBox)
        //if (isStartVtxInsideOfEllipse && isEndVtxInsideOfEllipse)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline void ShortestPathFinder2D::get_geodesic_path_by_BFS(list<EdgeForSolutionPoolGraph*>& path, const int & index)
{
    if ((index >= 0) && (index < m_PathsByTPT.size()))
    {
        path = m_PathsByTPT[index];
    }
}


inline void ShortestPathFinder2D::get_geodesic_path_by_BFS_with_priority(list<EdgeForSolutionPoolGraph*>& path, const int & priority)
{
    if ((priority >= 0) && (priority < m_IndexOfBFSsWithIncreasingOrderOfDistance.size()))
    {
        path = m_PathsByTPT[m_IndexOfBFSsWithIncreasingOrderOfDistance[priority]];
    }
}


inline double ShortestPathFinder2D::get_path_total_distance_by_BFS(const int & index) const
{
    if ((index >= 0) && (index < m_PathTotalDistancesByTPT.size()))
    { 
        return m_PathTotalDistancesByTPT[index]; 
    } 
    else
    {
        return DBL_MAX;
    }
}


inline double ShortestPathFinder2D::get_path_total_distance_by_BFS_with_priority(const int & priority) const
{
    if ((priority >= 0) && (priority < m_IndexOfBFSsWithIncreasingOrderOfDistance.size()))
    {
        return m_PathTotalDistancesByTPT[m_IndexOfBFSsWithIncreasingOrderOfDistance[priority]];
    }
    else
    {
        return DBL_MAX;
    }
}


inline void ShortestPathFinder2D::get_topological_path_by_BFS(list<EdgeBU2D*>& topologicalPath, const int & index)
{
    if ((index >= 0) && (index < m_TopologicalPathsByTPT.size()))
    { 
        topologicalPath = m_TopologicalPathsByTPT[index]; 
    } 
}


inline void ShortestPathFinder2D::get_topological_path_by_BFS_with_priority(list<EdgeBU2D*>& topologicalPath, const int & priority)
{
    if ((priority >= 0) && (priority < m_IndexOfBFSsWithIncreasingOrderOfDistance.size()))
    {
        topologicalPath = m_TopologicalPathsByTPT[m_IndexOfBFSsWithIncreasingOrderOfDistance[priority]];
    }
}


inline ShortestPathSolutionPoolGraph & ShortestPathFinder2D::get_solution_pool_graph_by_BFS(const int & index)
{
    if ((index >= 0) && (index < m_SolutionPoolGraphsByTPT.size()))
    { 
        return m_SolutionPoolGraphsByTPT[index]; 
    }
    else
    {
        return m_SolutionPoolGraphsByTPT[0];
    }
}


inline ShortestPathSolutionPoolGraph & ShortestPathFinder2D::get_solution_pool_graph_by_BFS_with_priority(const int & priority)
{
    if ((priority >= 0) && (priority < m_IndexOfBFSsWithIncreasingOrderOfDistance.size()))
    {
        return m_SolutionPoolGraphsByTPT[m_IndexOfBFSsWithIncreasingOrderOfDistance[priority]];
    }
    else
    {
        return m_SolutionPoolGraphsByTPT[0];
    }
}


inline bool ShortestPathFinder2D::is_topological_path_by_BFS_computed()
{
    return !m_TopologicalPathsByTPT.empty();
}


inline bool ShortestPathFinder2D::are_this_line_intersected_with_one_of_three_disks_of_QT_face(const rg_ImplicitEquation& tangentLineEq, 
                                                                                        FaceBU2D* qtFace, 
                                                                                        const rg_Point2D& startPt, 
                                                                                        const rg_Point2D& endPt,
                                                                                        unordered_set<VertexBU2D*>& visitedObstacles) const
{
    rg_dList<VertexBU2D*> qtVertices;
    qtFace->getBoundingVertices(qtVertices);

    return are_this_line_intersected_with_these_disks(tangentLineEq, qtVertices, startPt, endPt, visitedObstacles);
}




inline FaceBU2D* ShortestPathFinder2D::get_opposite_QT_face_of_Edge(EdgeBU2D* qtEdge, FaceBU2D* qtFace) const
{
    if (qtEdge->getRightFace() == qtFace)
    {
        return qtEdge->getLeftFace();
    }
    else if (qtEdge->getLeftFace() == qtFace)
    {
        return qtEdge->getRightFace();
    }
    else
    {
        return NULL;
    }
}


inline bool ShortestPathFinder2D::is_this_QT_edge_in_an_exterior_state(EdgeBU2D* qtEdge) const
{
    if (qtEdge->getBoundingState(get_probe().getRadius()) == EXTERIOR_SIMPLEX)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool ShortestPathFinder2D::is_this_QT_face_in_an_exterior_state(FaceBU2D* qtFace) const
{
    if (qtFace->getBoundingState(get_probe().getRadius()) == EXTERIOR_SIMPLEX)
    {
        return true;
    }
    else
    {
        return false;
    }
}



inline bool ShortestPathFinder2D::is_this_QT_face_visited(FaceBU2D* qtFace) const
{
    if (qtFace->isVisited())
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline void ShortestPathFinder2D::set_prev_and_curr_QT_face_for_loop(FaceBU2D*& currQTFace, FaceBU2D*& prevQTFace, EdgeBU2D* currQTEdge) const
{
    prevQTFace = currQTFace;

    if (currQTEdge->getRightFace() == currQTFace)
    {
        currQTFace = currQTEdge->getLeftFace();
    }
    else
    {
        currQTFace = currQTEdge->getRightFace();
    }
}


inline bool ShortestPathFinder2D::is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq,
                                                                                                        const rg_Point2D& startPtOfTangentLine, 
                                                                                                        const rg_Point2D& endPtOfTangentLine,
                                                                                                        BetaUniverse2D& betaUniverse) const
{
     unordered_set<VertexBU2D*> visitedObstacles;
    
     return is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEq, startPtOfTangentLine, endPtOfTangentLine, betaUniverse, visitedObstacles);
}


inline bool ShortestPathFinder2D::is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                                                        const rg_Point2D& startPtOfTangentLine, 
                                                                                                        const rg_Point2D& endPtOfTangentLine,
                                                                                                        BetaUniverse2D& betaUniverse,
                                                                                                        unordered_set<VertexBU2D*>& visitedObstacles)const
{
    FaceBU2D*  startQTFace = betaUniverse.findFaceContainingInputPoint(startPtOfTangentLine);
    FaceBU2D*  endQTFace   = betaUniverse.findFaceContainingInputPoint(endPtOfTangentLine);

    return is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEq, 
                                                                                 startPtOfTangentLine, 
                                                                                 endPtOfTangentLine,
                                                                                 betaUniverse, 
                                                                                 startQTFace, 
                                                                                 endQTFace, 
                                                                                 visitedObstacles);
}

inline bool ShortestPathFinder2D::is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq, 
                                                                                                 const rg_Point2D& startPtOfTangentLine, 
                                                                                                 const rg_Point2D& endPtOfTangentLine, 
                                                                                                 BetaUniverse2D& betaUniverse, 
                                                                                                 VertexBU2D* obs1, VertexBU2D* obs2, 
                                                                                                 unordered_set<VertexBU2D*>& visitedObstacles) const
{
    FaceBU2D*  startQTFace = betaUniverse.findFaceContainingInputPointWithNearestVtx(startPtOfTangentLine, obs1);
    FaceBU2D*  endQTFace  = betaUniverse.findFaceContainingInputPointWithNearestVtx(endPtOfTangentLine, obs2);

    return is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEq, 
                                                                                 startPtOfTangentLine, 
                                                                                 endPtOfTangentLine, 
                                                                                 betaUniverse, 
                                                                                 startQTFace, 
                                                                                 endQTFace, 
                                                                                 visitedObstacles);

}



inline bool ShortestPathFinder2D::is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(const rg_ImplicitEquation& lineEq,
                                                                                                 const rg_Point2D& startPtOfTangentLineFromPoint, 
                                                                                                 const rg_Point2D& endPtOfTangentLineFromObstacle, 
                                                                                                 BetaUniverse2D& betaUniverse, 
                                                                                                 VertexBU2D* obstacle, 
                                                                                                 unordered_set<VertexBU2D*>& visitedObstacles) const
{
     FaceBU2D*  startQTFace = betaUniverse.findFaceContainingInputPoint(startPtOfTangentLineFromPoint);
     FaceBU2D*  endQTFace   = betaUniverse.findFaceContainingInputPointWithNearestVtx(endPtOfTangentLineFromObstacle, obstacle);

     return is_a_line_segment_of_these_two_points_intersected_with_some_obstacles(lineEq, 
                                                                                  startPtOfTangentLineFromPoint, 
                                                                                  endPtOfTangentLineFromObstacle, 
                                                                                  betaUniverse, 
                                                                                  startQTFace, 
                                                                                  endQTFace,
                                                                                  visitedObstacles);
}



inline void ShortestPathFinder2D::generate_graph_incrementally_by_adding_two_obstacles_into_the_graph(VertexBU2D* obs1,
                                                                                                      VertexBU2D* obs2,
                                                                                                      const BetaComponentFinder& betaComponentFinder,
                                                                                                      BetaUniverse2D& betaUniverse, 
                                                                                                      ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                                                      unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk)
{

    //1. make graph of tangent lines
    vector<VertexForSolutionPoolGraph*> verticesOnObstacle1;
    vector<VertexForSolutionPoolGraph*> verticesOnObstacle2;

    make_graph_of_tangent_lines(obs1, obs2, betaComponentFinder, betaUniverse, solutionPoolGraph, verticesOnObstacle1, verticesOnObstacle2);
    
    //2. add_new_arcs_made_by_tangent_point
    add_new_arc_edges_made_by_new_vertices_into_the_graph(obs1, verticesOnObstacle1, betaComponentFinder, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
    add_new_arc_edges_made_by_new_vertices_into_the_graph(obs2, verticesOnObstacle2, betaComponentFinder, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}


inline void ShortestPathFinder2D::generate_graph_incrementally_by_adding_two_obstacles_into_the_graph_with_topological_path(VertexBU2D * obs1, VertexBU2D * obs2, const list<EdgeBU2D*>& topologicalPath, const BetaComponentFinder & betaComponentFinder, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk)
{
    //1. make graph of tangent lines
    vector<VertexForSolutionPoolGraph*> verticesOnObstacle1;
    vector<VertexForSolutionPoolGraph*> verticesOnObstacle2;

    make_graph_of_tangent_lines_with_topological_path(obs1, obs2, topologicalPath ,betaComponentFinder, betaUniverse, solutionPoolGraph, verticesOnObstacle1, verticesOnObstacle2);
    
    //2. add_new_arcs_made_by_tangent_point
    add_new_arc_edges_made_by_new_vertices_into_the_graph(obs1, verticesOnObstacle1, betaComponentFinder, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
    add_new_arc_edges_made_by_new_vertices_into_the_graph(obs2, verticesOnObstacle2, betaComponentFinder, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}


inline void ShortestPathFinder2D::insert_start_end_point_into_the_graph(const rg_Circle2D& startPt, 
                                                                        const rg_Circle2D& endPt, 
                                                                        const vector<VertexBU2D*>& obstaclePool,
                                                                        BetaUniverse2D& betaUniverse, 
                                                                        ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                        unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    VertexForSolutionPoolGraph* startVtx = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
    VertexForSolutionPoolGraph* endVtx = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());

    startVtx->set_tangent_point_to_disk(startPt.getCenterPt());
    startVtx->set_accumulate_length_from_source_vtx(DBL_MAX);

    endVtx->set_tangent_point_to_disk(endPt.getCenterPt());
    endVtx->set_accumulate_length_from_source_vtx(DBL_MAX);

    solutionPoolGraph.set_start_vtx(startVtx);
    solutionPoolGraph.set_end_vtx(endVtx);

    insert_a_point_into_the_graph(startVtx, obstaclePool, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
    insert_a_point_into_the_graph(endVtx, obstaclePool, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}



inline void ShortestPathFinder2D::insert_start_end_point_into_the_graph_with_topological_path(const rg_Circle2D & startPt, const rg_Circle2D & endPt, const vector<VertexBU2D*>& obstaclePool, const list<EdgeBU2D*>& topologicalPath, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    VertexForSolutionPoolGraph* startVtx = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
    VertexForSolutionPoolGraph* endVtx = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());

    startVtx->set_tangent_point_to_disk(startPt.getCenterPt());
    startVtx->set_accumulate_length_from_source_vtx(DBL_MAX);

    endVtx->set_tangent_point_to_disk(endPt.getCenterPt());
    endVtx->set_accumulate_length_from_source_vtx(DBL_MAX);

    solutionPoolGraph.set_start_vtx(startVtx);
    solutionPoolGraph.set_end_vtx(endVtx);

    insert_a_point_into_the_graph_with_topological_path(startVtx, obstaclePool, topologicalPath, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
    insert_a_point_into_the_graph_with_topological_path(endVtx, obstaclePool, topologicalPath, betaUniverse, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}


inline void ShortestPathFinder2D::generate_graph_incrementally_by_adding_one_obstacle_N_point_into_the_graph(VertexBU2D* obstacle, 
                                                                                                             VertexForSolutionPoolGraph* pointVtx, 
                                                                                                             const BetaComponentFinder& betaComponentFinder, 
                                                                                                             BetaUniverse2D& betaUniverse, 
                                                                                                             ShortestPathSolutionPoolGraph& solutionPoolGraph, 
                                                                                                             unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    //1. make graph of tangent lines
    vector<VertexForSolutionPoolGraph*> verticesOnObstacle;

    make_graph_of_tangent_lines(obstacle, pointVtx, betaUniverse, solutionPoolGraph, verticesOnObstacle);

    //2. add_new_arcs_made_by_tangent_point
    add_new_arc_edges_made_by_new_vertices_into_the_graph(obstacle, verticesOnObstacle, betaComponentFinder, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}


inline void ShortestPathFinder2D::generate_graph_incrementally_by_adding_one_obstacle_N_point_into_the_graph_with_topological_path(VertexBU2D * obstacle, VertexForSolutionPoolGraph * pointVtx, const list<EdgeBU2D*>& topologicalPath, const BetaComponentFinder & betaComponentFinder, BetaUniverse2D & betaUniverse, ShortestPathSolutionPoolGraph & solutionPoolGraph, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    //1. make graph of tangent lines
    vector<VertexForSolutionPoolGraph*> verticesOnObstacle;

    make_graph_of_tangent_lines_with_topological_path(obstacle, pointVtx, topologicalPath, betaUniverse, solutionPoolGraph, verticesOnObstacle);

    //2. add_new_arcs_made_by_tangent_point
    add_new_arc_edges_made_by_new_vertices_into_the_graph(obstacle, verticesOnObstacle, betaComponentFinder, solutionPoolGraph, mapperForQTVtxToVtxOnDisk);
}



inline void ShortestPathFinder2D::find_tangent_line_equations_of_two_obstacles(VertexBU2D* obs1, VertexBU2D* obs2,
                                                                               const BetaComponentFinder& betaComponentFinder, 
                                                                               vector<rg_ImplicitEquation>& tangentLineEqs)
{
    if (are_the_tangent_lines_of_this_obstacle_pair_already_computed(obs1, obs2))
    {
        get_tangent_lines_of_this_obstacle_pair_from_already_computed(obs1, obs2, tangentLineEqs);
    }
    else
    {
        //1. find tangent line
        rg_Circle2D obstacleDisk1(obs1->getCoord(), obs1->getCircle().getRadius() + m_Probe.getRadius());
        rg_Circle2D obstacleDisk2(obs2->getCoord(), obs2->getCircle().getRadius() + m_Probe.getRadius());

        if (are_these_two_obstacles_in_same_component(obs1, obs2, betaComponentFinder))
        {
            find_tangent_line_equations_of_two_obstacles_in_same_component(obstacleDisk1, obstacleDisk2, tangentLineEqs);
        }
        else
        {
            find_tangent_line_equations_of_two_obstacles_in_diff_component(obstacleDisk1, obstacleDisk2, tangentLineEqs);
        }

        make_mapper_from_QT_obstacle_pair_to_its_tangent_lines(obs1, obs2, tangentLineEqs, m_ObstaclePairNItsTangentLines);
    }
}


inline void ShortestPathFinder2D::find_tangent_line_equations_of_point_and_obstacle(VertexForSolutionPoolGraph* pointVtx, VertexBU2D* obstacle, vector<rg_ImplicitEquation>& tangentLineEqs) const
{
    rg_Circle2D point(pointVtx->get_tangent_point_to_disk(), 0.0);
    rg_Circle2D obstacleDisk(obstacle->getCoord(), obstacle->getCircle().getRadius() + m_Probe.getRadius());

    find_tangent_line_equations_of_two_obstacles_in_same_component(point, obstacleDisk, tangentLineEqs);
}



inline bool ShortestPathFinder2D::does_this_obstacle_have_already_maden_vetices_of_tangent_points(VertexBU2D* obstacle, const unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>::const_iterator it_Obstacle = mapperForQTVtxToVtxOnDisk.find(obstacle);

    if (it_Obstacle == mapperForQTVtxToVtxOnDisk.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}



inline void ShortestPathFinder2D::insert_new_vertex_into_set_of_vertices_on_obstacle(VertexForSolutionPoolGraph* newVertexOnObstacle, VertexBU2D* obstacle, unordered_map<VertexBU2D*, list<VertexForSolutionPoolGraph*>>& mapperForQTVtxToVtxOnDisk) const
{
    if (does_this_obstacle_have_already_maden_vetices_of_tangent_points(obstacle, mapperForQTVtxToVtxOnDisk))
    {
        list<VertexForSolutionPoolGraph*>& verticesOnObstacleBeforeMakingNewVertex = mapperForQTVtxToVtxOnDisk.at(obstacle);

        insert_into_list_of_tangent_vertices_on_disk_by_angle(obstacle, newVertexOnObstacle, verticesOnObstacleBeforeMakingNewVertex);
        //verticesOnObstacleBeforeMakingNewVertex.push_back(newVertexOnObstacle);
    }
    else
    {
        list<VertexForSolutionPoolGraph*> verticesOnObstacle;
        verticesOnObstacle.push_back(newVertexOnObstacle);

        mapperForQTVtxToVtxOnDisk.insert(make_pair(obstacle, verticesOnObstacle));
    }
}



inline void ShortestPathFinder2D::create_two_vertices_and_a_tangent_line_edge(VertexForSolutionPoolGraph*& vertex1, 
                                                                              VertexForSolutionPoolGraph*& vertex2, 
                                                                              EdgeForSolutionPoolGraph*& edge,
                                                                              ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    vertex1 = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
    vertex2 = solutionPoolGraph.create_vertex(VertexForSolutionPoolGraph());
    edge    = solutionPoolGraph.create_edge(EdgeForSolutionPoolGraph());
}


inline void ShortestPathFinder2D::set_two_vertices_and_a_tangent_line_edge(VertexForSolutionPoolGraph* vertex1, 
                                                                           VertexForSolutionPoolGraph* vertex2, 
                                                                           EdgeForSolutionPoolGraph* edge, 
                                                                           const pair<rg_Point2D, rg_Point2D>& twoPtsOnObstacle1And2, 
                                                                           const rg_ImplicitEquation& implicitEquationOfTangentLine) const
{
    //2. Set two vertices
    const rg_Point2D& ptOnObstacle1 = twoPtsOnObstacle1And2.first;
    const rg_Point2D& ptOnObstacle2 = twoPtsOnObstacle1And2.second;

    vertex1->set_tangent_point_to_disk(ptOnObstacle1);
    vertex1->add_line_edge(edge);
    vertex1->set_accumulate_length_from_source_vtx(DBL_MAX);
    vertex1->set_angle(-1.0);

    vertex2->set_tangent_point_to_disk(ptOnObstacle2);
    vertex2->add_line_edge(edge);
    vertex2->set_accumulate_length_from_source_vtx(DBL_MAX);
    vertex2->set_angle(-1.0);

    //2. Set a tangent line edge
    double tangentLineSementLength = ptOnObstacle1.distance(ptOnObstacle2);

    edge->set_start_vtx(vertex1);
    edge->set_end_vtx(vertex2);
    edge->set_edge_type(LINE_SEGMENT_SPG);
    edge->set_edge_length(tangentLineSementLength);
    edge->set_implicit_equation(implicitEquationOfTangentLine);
}


inline void ShortestPathFinder2D::create_an_arc_edge(EdgeForSolutionPoolGraph*& edge, ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    edge = solutionPoolGraph.create_edge(EdgeForSolutionPoolGraph());
}


inline void ShortestPathFinder2D::set_an_arc_edge(EdgeForSolutionPoolGraph* arcEdge,
                                                  VertexForSolutionPoolGraph* startVertex, VertexForSolutionPoolGraph* endVertex, 
                                                  VertexBU2D* obstacle, ShortestPathSolutionPoolGraph& solutionPoolGraph) const
{
    rg_Point2D arcStartPoint = startVertex->get_tangent_point_to_disk();
    rg_Point2D arcEndPoint   = endVertex->get_tangent_point_to_disk();

    rg_ImplicitEquation arcInplicitEquation = make_implicit_equation_of_arc_edge(obstacle, m_Probe.getRadius());
    double arcLength                        = compute_arc_segment_length_on_obstacle_boundary(obstacle, arcStartPoint, arcEndPoint);

    arcEdge->set_start_vtx(startVertex);
    arcEdge->set_end_vtx(endVertex);
    arcEdge->set_edge_type(ARC_SPG);
    arcEdge->set_edge_length(arcLength);
    arcEdge->set_implicit_equation(arcInplicitEquation);

    if (startVertex == endVertex)
    {
        startVertex->add_arc_edge(arcEdge);
    }
    else
    {
        startVertex->add_arc_edge(arcEdge);
        endVertex->add_arc_edge(arcEdge);
    }
}



inline void ShortestPathFinder2D::set_two_arc_edges(EdgeForSolutionPoolGraph * newArcEdge, VertexForSolutionPoolGraph * prevVertex, VertexForSolutionPoolGraph * newVertex, VertexForSolutionPoolGraph * nextVertex, VertexBU2D * obstacle, ShortestPathSolutionPoolGraph & solutionPoolGraph) const
{
    rg_Point2D arcPrevPoint  = prevVertex->get_tangent_point_to_disk();
    rg_Point2D arcNewPoint   = newVertex->get_tangent_point_to_disk();
    rg_Point2D arcNextPoint  = nextVertex->get_tangent_point_to_disk();

    rg_ImplicitEquation arcInplicitEquation = make_implicit_equation_of_arc_edge(obstacle, m_Probe.getRadius());
    double arcLengthOfPrevToNew = compute_arc_segment_length_on_obstacle_boundary(obstacle, arcPrevPoint, arcNewPoint);
    double arcLengthOfNewToNext = compute_arc_segment_length_on_obstacle_boundary(obstacle, arcNewPoint, arcNextPoint);

    //1. SET OLD ARC EDGE
    EdgeForSolutionPoolGraph * oldArcEdge = find_arc_edge_having_these_two_vertices(prevVertex, nextVertex);

    oldArcEdge->set_start_vtx(prevVertex);
    oldArcEdge->set_end_vtx(newVertex);
    oldArcEdge->set_edge_type(ARC_SPG);
    oldArcEdge->set_edge_length(arcLengthOfPrevToNew);
    oldArcEdge->set_implicit_equation(arcInplicitEquation);

    //prevVertex->add_arc_edge(oldArcEdge);
    newVertex->add_arc_edge(oldArcEdge);


    //2. SET NEW ARC EDGE
    newArcEdge->set_start_vtx(newVertex);
    newArcEdge->set_end_vtx(nextVertex);
    newArcEdge->set_edge_type(ARC_SPG);
    newArcEdge->set_edge_length(arcLengthOfNewToNext);
    newArcEdge->set_implicit_equation(arcInplicitEquation);

    newVertex->add_arc_edge(newArcEdge);
    nextVertex->add_arc_edge(newArcEdge);
    nextVertex->erase_arc_edge(oldArcEdge); //erase
}


inline double ShortestPathFinder2D::get_angle_of_start_to_end_point(const rg_Point2D& startPt, const rg_Point2D& endPt) const
{
    rg_Point2D vecFromStartToEndPt = (endPt - startPt);

    double angle = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecFromStartToEndPt);

    return angle;
}

inline bool ShortestPathFinder2D::do_any_of_two_intervals_include_none_of_obstacles(const double& angle1, 
                                                                                    const double& angle2, 
                                                                                    const list<double>& anglesOfNeighborIntersectedObstacles, 
                                                                                    int& intervalOption) const
{
    bool doAnyOfTwoItervalsIncludeNoneOfObstacles = true;

    if (do_this_interval_include_none_of_obstacles(angle1, angle2, anglesOfNeighborIntersectedObstacles))
    {
        doAnyOfTwoItervalsIncludeNoneOfObstacles = true;
        intervalOption = 1;
    }
    else if (do_this_interval_include_none_of_obstacles(angle2, angle1, anglesOfNeighborIntersectedObstacles))
    {
        doAnyOfTwoItervalsIncludeNoneOfObstacles = true;
        intervalOption = 2;
    }
    else 
    {
        doAnyOfTwoItervalsIncludeNoneOfObstacles = false;
        intervalOption = 0;
    }

    /*
    doAnyOfTwoItervalsIncludeNoneOfObstacles = do_this_interval_include_none_of_obstacles(angle1, 
                                                                                          angle2, 
                                                                                          anglesOfNeighborIntersectedObstacles)
                                            || do_this_interval_include_none_of_obstacles(angle2, 
                                                                                          angle1, 
                                                                                          anglesOfNeighborIntersectedObstacles);

                                                                                          */
    return doAnyOfTwoItervalsIncludeNoneOfObstacles;
}

inline bool ShortestPathFinder2D::do_any_of_two_intervals_include_none_of_obstacles(const double & angle1, const double & angle2, const list<double>& anglesOfNeighborIntersectedObstacles) const
{
    bool doAnyOfTwoItervalsIncludeNoneOfObstacles = true;

    if (do_this_interval_include_none_of_obstacles(angle1, angle2, anglesOfNeighborIntersectedObstacles))
    {
        doAnyOfTwoItervalsIncludeNoneOfObstacles = true;
    }

    else 
    {
        doAnyOfTwoItervalsIncludeNoneOfObstacles = false;
    }         

    return doAnyOfTwoItervalsIncludeNoneOfObstacles;
}


inline bool ShortestPathFinder2D::do_any_of_QT_faces_include_this_pt(const rg_Point2D & pt, const unordered_set<FaceBU2D*>& qtFaces, BetaUniverse2D& betaUniverse) const
{
    bool doContainPt = false;

    for (unordered_set<FaceBU2D*>::const_iterator it_QTFace = qtFaces.begin();
         it_QTFace != qtFaces.end();
         ++it_QTFace)
    {
        FaceBU2D* qtFace = *it_QTFace;

        if (qtFace->isVirtual())
        {
            if (betaUniverse.thisFaceContainsInputPointForVirtualQTFace(qtFace, pt))
            {
                doContainPt = true;
                break;
            }
        }
        else
        {
            if (betaUniverse.thisFaceContainsInputPoint(qtFace, pt))
            {
                doContainPt = true;
                break;
            }
        }
    }

    return doContainPt;
}


#endif