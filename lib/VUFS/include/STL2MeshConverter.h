#ifndef STL2MESHCONVERTER_H
#define STL2MESHCONVERTER_H

#include "rg_Triangle.h"
#include "FaceVertexMesh.h"
#include "PolygonMeshModel.h"

#include "AxisAlignedBox.h"

#include <string>
#include <list>
using namespace std;


class STL2MeshConverter
{
public:
    class Triangle;

private:
    list< Triangle > m_triangles;
    AxisAlignedBox   m_boundingBox;

    set< Triangle* > m_triangles_with_same_vertex;


public:
    STL2MeshConverter();
    ~STL2MeshConverter();

    bool    read(const string& stl_filename);

    bool    convert_to(FaceVertexMesh& faceVertexMeshModel) const;
    bool    convert_to(PolygonMeshModel& polygonMeshModel) const;

    bool    make_mesh_model_from_ply_file( const string& ply_filename, PolygonMeshModel& polygonMeshModel );

    static  bool write(const string& stl_filename, const PolygonMeshModel& model, const bool& ASCIIOrNot=true);

private:
    bool    is_available_file(    const string& stl_filename );
    bool    is_ASCII_STL_file(    const string& stl_filename );
    bool    read_ASCII_STL_file(  const string& stl_filename );
    bool    read_binary_STL_file( const string& stl_filename );

    void    check_triangles_with_same_vertex();

    AxisAlignedBox  get_bounding_box() const;

    static  bool write_ASCII_file(const string& stl_filename, const PolygonMeshModel& model);
    static  bool write_binary_file(const string& stl_filename, const PolygonMeshModel& model);

};



class STL2MeshConverter::Triangle : public rg_Triangle {
private:
    rg_Point3D m_normal;

public:
    Triangle();
    Triangle(const rg_Point3D& normal, const rg_Point3D& vertex1, const rg_Point3D& vertex2, const rg_Point3D& vertex3);
    Triangle(const Triangle& triangle);
    ~Triangle();

    Triangle&   operator =(const Triangle& triangle);

    rg_Point3D  normal() const;
    void        set_normal(const rg_Point3D& normal);

    void        set(const rg_Point3D& normal, const rg_Point3D& vertex1, const rg_Point3D& vertex2, const rg_Point3D& vertex3 );


    bool        are_vertices_with_same_coordinates() const;
};


#endif

