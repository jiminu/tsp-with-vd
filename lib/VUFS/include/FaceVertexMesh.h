#ifndef FACEVERTEXMESH_H
#define FACEVERTEXMESH_H

#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "rg_Point3D.h"

#include "PolygonMeshModel.h"

#include <list>
#include <map>
#include <set>
#include <vector>
using namespace std;


class FaceVertexMesh
{
public:
    class Face;
    class Edge;
    class Vertex;

private:
    list<Face*>     m_faces;
    list<Vertex*>   m_vertices;

    unsigned int    m_number_of_faces;
    unsigned int    m_number_of_vertices;

public:
    FaceVertexMesh();
    FaceVertexMesh(const FaceVertexMesh& FVmesh);
    ~FaceVertexMesh();

    FaceVertexMesh& operator =(const FaceVertexMesh& FVmesh);

    void    clear();

    unsigned int    number_of_faces() const;
    unsigned int    number_of_vertices() const;

    const list<Face*>&     get_faces() const;
    const list<Vertex*>&   get_vertices() const;

    Face*           create_face();
    Vertex*         create_vertex();

    void            remove_face(const Face* const face);
    void            remove_vertex(const Vertex* const vertex);


    Vertex*         find(const rg_Point3D& point) const;


    bool            convert_to(PolygonMeshModel& polygonMeshModel) const;


    void            check_face_with_same_vertex( list<Face*>&   facesWithSameVertex ) const;
    void            check_non_manifold_vertices( list<Vertex*>& nonManifoldVertices ) const;

private:
    void            make_map_from_FaceVertexMesh_to_PolygonMesh(
                                PolygonMeshModel&                                        polygonMeshModel,
                                map<FaceVertexMesh::Face*, PolygonMeshModel::Face*>&     faceMap,
                                map<FaceVertexMesh::Vertex*, PolygonMeshModel::Vertex*>& vertexMap  ) const;

    void            construct_shells_of_PolygonMesh(
                                PolygonMeshModel&                                        polygonMeshModel,
                                map<FaceVertexMesh::Face*, PolygonMeshModel::Face*>&     faceMap,
                                map<FaceVertexMesh::Vertex*, PolygonMeshModel::Vertex*>& vertexMap  ) const;
    void            construct_bodies_of_PolygonMesh(
                                PolygonMeshModel&                                        polygonMeshModel ) const; 
};




///////////////////////////////////////////////////////////////////////////////
//
// class FaceVertexMesh::Face 
//
class FaceVertexMesh::Face : public TopologicalEntity
{
private:
    typedef FaceVertexMesh::Edge    Edge;
    typedef FaceVertexMesh::Vertex  Vertex;

    Vertex*         m_bounding_vertex[3];
    rg_Point3D      m_normal;


public:
    Face();
    Face(const rg_INT& ID);
    Face(const Face& face);
    virtual ~Face();

    Face& operator =(const Face& face);

    Vertex*         bounding_vertex(const int& pos) const;
    rg_Point3D      normal() const;
    void            set_bounding_vertex(const int& pos, const Vertex* const vertex);
    void            set_normal(const rg_Point3D& normal);

    bool            is_incident_to(const Vertex* const vertex) const;
    bool            is_incident_to(const Edge& edge) const;

    bool            has_same_vertex() const;

    Edge            bounding_edge(const int& pos) const;
    Face*           adjacent_face(const int& pos) const;

    int             find_edge_position(const Edge& edge) const;
    int             find_vertex_position(const Vertex* const vertex) const;

    unsigned int    number_of_bounding_vertices() const;

    unsigned int    find_bounding_vertices( list<Vertex*>& boundingVertices ) const;
    unsigned int    find_bounding_edges(    list<Edge>&    boundingEdges    ) const;
    unsigned int    find_adjacent_faces(    list<Face*>&   adjacentFaces    ) const;

    bool            find_edges_incident_to_vertex(const Vertex* const vertex, Edge& inEdge, Edge& outEdge) const;

};




///////////////////////////////////////////////////////////////////////////////
//
// class FaceVertexMesh::Edge 
//
class FaceVertexMesh::Edge
{
private:
    typedef FaceVertexMesh::Face    Face;
    typedef FaceVertexMesh::Vertex  Vertex;

    Face*       m_face;
    Vertex*     m_start_vertex;
    Vertex*     m_end_vertex;

public:
    Edge();
    Edge(const Face* const face, const Vertex* const startVertex, const Vertex* const endVertex );
    Edge(const Edge& edge);
    ~Edge();

    Edge& operator =(const Edge& edge);

    Face*        face() const;
    Vertex*      start_vertex() const;
    Vertex*      end_vertex() const;
                 
    void         set_face(         const Face*   const face );
    void         set_start_vertex( const Vertex* const vertex );
    void         set_end_vertex(   const Vertex* const vertex );
    void         set( const Face*   const face,
                      const Vertex* const startVertex, 
                      const Vertex* const endVertex );
                 
    Face*        find_twin_face() const;

    unsigned int find_incident_faces( list<Face*>& incidentFaces ) const;
};


///////////////////////////////////////////////////////////////////////////////
//
// class FaceVertexMesh::Vertex 
//
class FaceVertexMesh::Vertex : public TopologicalEntity
{
private:
    typedef FaceVertexMesh::Face Face;

    rg_Point3D  m_coordinate;
    set<Face*>  m_incident_faces;

public:
    Vertex();
    Vertex(const rg_INT& ID);
    Vertex(const Vertex& vertex);
    virtual ~Vertex();

    Vertex& operator =(const Vertex& vertex);

    rg_Point3D      coordinate() const;

    void            set_coordinate( const rg_Point3D& coordinate );
    void            add_incident_face( const Face* const face );
    void            remove_incident_face( const Face* const face );

    bool            has_same_coordinate(const rg_Point3D& point, const rg_REAL& res=rg_SYSTEM_RES);

    bool            is_equal(const Vertex* const vertex, const rg_REAL& res=rg_SYSTEM_RES);
    bool            is_member_of( const Face* const face ) const;

    unsigned int    number_of_incident_faces() const;

    unsigned int    find_incident_faces(list<Face*>& incidentFaces) const;


};



#endif

