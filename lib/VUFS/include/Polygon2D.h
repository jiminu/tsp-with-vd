#ifndef POLYGON2D_H
#define POLYGON2D_H

/***************************************************************************************/
/* 1. Summary                                                                          */
/*   This class represents a polygon in 2-dimesional space.                            */
/*   It can include many polygonal holes in a big exterior polygon boundary.           */
/*   The first list element in m_boundaryVertices is for exterior polygon boundary.    */
/*   The boundary vertices are in CCW order.                                           */
/*                                                                                     */
/* 2. History                                                                          */
/*   2020.03.04. Polygon2D using one array(Ver-A) is merged into this class by CYSONG. */
/***************************************************************************************/


#include "rg_Const.h"
#include "rg_BoundingBox2D.h"
#include "rg_Point3D.h"
class rg_Point2D;
#include <string>
#include <list>
using namespace std;

class Polygon2D
{
private:
	list< list<rg_Point2D> > m_boundaryVertices; // first(front) list is for exterior boundary.
	rg_INT                   m_numberOfVertices;
	rg_BoundingBox2D         m_boundingBox;

public:
    //Constructor and destructor
	Polygon2D();
    Polygon2D(const list<rg_Point2D>& exteriorBoundaryVertices); //ADDED BY CYSONG [FEB05, 20]
    Polygon2D(const rg_Point2D* exteriorBoundaryVertices, const rg_INT& numVertices);
    Polygon2D(const list< list<rg_Point2D> >& boundaryVertices);
	Polygon2D(const Polygon2D& polygon);
    ~Polygon2D();

    //Getter and setter
    inline void get_bounding_box(rg_BoundingBox2D& boundingBox) const;
    inline void get_boundary_vertices(list< list<rg_Point2D> >& boundaryVertices) const;
    void get_boundary_vertices(list<rg_Point2D>& boundaryVertices) const;
    rg_INT get_exterior_boundary_vertices(list<rg_Point2D>& exteriorBoundaryVertices) const;
    rg_INT get_exterior_boundary_vertices(rg_Point2D*& exteriorBoundaryVertices) const;

    inline list< list<rg_Point2D> >&       boundary_vertices();
    inline const list< list<rg_Point2D> >& boundary_vertices() const;
    inline rg_BoundingBox2D&               bounding_box();
    inline const rg_BoundingBox2D&         bounding_box() const;

    void set_exterior_boundary_vertices(const rg_Point2D* boundaryVertex, const rg_INT& numVertices);
    void set_exterior_boundary_vertices(const list<rg_Point2D>& exteriorBoundaryVertices); //ADDED BY CYSONG [FEB05, 20]
    void set_boundary_vertices(const list< list<rg_Point2D> >& boundaryVertices);

    //Operator
    Polygon2D& operator=(const Polygon2D& polygon);

    //Function
    void    clear();

    rg_REAL find_minimum_edge_length()    const;
    rg_REAL compute_length_of_perimeter() const;
    rg_REAL compute_area()                const;

    void split_all_edges_which_are_longer_than(const rg_REAL& minEdgeLengthToSplit, const rg_REAL& edgeLengthToSurvive);
    void split_all_edges_which_are_longer_than(const rg_REAL& minEdgeLengthToSplit);

    void merge_all_consecutive_edge_which_are_on_same_line();
    void perturb_vertices_which_has_same_coordinate();

    void create_polygon_as_rectangle(const rg_Point2D& minPt, const rg_Point2D& maxPt);
    void translate(const rg_Point2D& translationAmount);

    //Queries: it seems to be updated [commented by CYSONG, March 06,2020]
    bool contains(const rg_Point2D& queryPoint);
    bool this_disk_is_outside_of_polygon_exterior_boundary(const rg_Circle2D& disk, const rg_REAL& tolerance) const;
    bool this_disk_intersects_or_is_outside_of_polygon_exterior_boundary(const rg_Circle2D& disk, const rg_REAL& tolerance) const;
    bool this_polygon_is_convex(const list<rg_Point2D>& boundaryVertices, list<rg_Point2D*>& reflexVertices) const;
    inline bool this_polygon_vertex_is_reflex(const list<rg_Point2D>& boundaryVertices, list<rg_Point2D>::const_iterator& i_vertex) const;

    //File I/O: (August 23, 2020 by Joonghyun)
    bool read_polygon_file(const std::string& filename);
    bool write_polygon_file(const std::string& filename);
    static bool read_polygon_file(list<Polygon2D>& polygons,  const std::string& filename,       bool& bFirstPolygonIsContainer);
    static bool write_polygon_file(list<Polygon2D>& polygons, const std::string& filename, const bool& bFirstPolygonIsContainer);

private:
    static bool read_a_polygon(std::istringstream& istrstream, std::string& line, Polygon2D& polygon, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge);
    static void write_description_of_polygon_file_in_commented_lines(ofstream& fout);
    static bool read_a_container_polygon(std::istringstream& istrstream, std::string& line, Polygon2D& polygon, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge);
    static bool write_a_polygon(ofstream& fout, Polygon2D& polygon, rg_INDEX& polygonIndex, rg_INDEX& vertexIndex, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge);
    static bool write_a_container_polygon(ofstream& fout, Polygon2D& polygon, rg_INDEX& polygonIndex, rg_INDEX& vertexIndex, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge);

public:
    //Statistics
    inline rg_INT num_vertices() const;
    inline rg_INT num_edges()    const;
    inline rg_INT num_holes()    const;
    inline rg_INT num_shells()   const;

	void compute_mean_and_standard_deviation_of_edge_length(rg_REAL& meanOfEdgeLength, rg_REAL& standardDeviationOfEdgeLength) const;
        rg_REAL compute_mean_of_edge_length() const;
        rg_REAL compute_standard_deviation_of_edge_length(const rg_REAL& meanOfEdgeLength) const;
	void compute_mean_and_standard_deviation_of_edge_length_which_is_longer_than(const rg_REAL& miminumLengthToMeasure, rg_REAL& meanOfEdgeLength, rg_REAL& standardDeviationOfEdgeLength) const;
	    rg_REAL compute_mean_of_edge_length_which_is_longer_than(const rg_REAL& miminumLengthToMeasure) const;
	    rg_REAL compute_standard_deviation_of_edge_length_which_is_longer_than(const rg_REAL& miminumLengthToMeasure, const rg_REAL& meanOfEdgeLength) const;
	void compute_mean_and_standard_deviation_of_edge_length_which_is_shorter_than(const rg_REAL& maximumLengthToMeasure, rg_REAL& meanOfEdgeLength, rg_REAL& standardDeviationOfEdgeLength) const;
	    rg_REAL compute_mean_of_edge_length_which_is_shorter_than(const rg_REAL& maximumLengthToMeasure) const;
	    rg_REAL compute_standard_deviation_of_edge_length_which_is_shorter_than(const rg_REAL& maximumLengthToMeasure, const rg_REAL& meanOfEdgeLength) const;

};

inline void Polygon2D::get_boundary_vertices(list< list< rg_Point2D > >& boundaryVertices) const
{
	boundaryVertices = m_boundaryVertices;
}

inline void Polygon2D::get_bounding_box(rg_BoundingBox2D & boundingBox) const
{
	boundingBox = m_boundingBox;
}

inline list< list< rg_Point2D > >& Polygon2D::boundary_vertices()
{
	return m_boundaryVertices;
}

inline const list< list<rg_Point2D> >& Polygon2D::boundary_vertices() const
{
    return m_boundaryVertices;
}

inline rg_BoundingBox2D& Polygon2D::bounding_box()
{
    return m_boundingBox;
}

inline const rg_BoundingBox2D& Polygon2D::bounding_box() const 
{
    return m_boundingBox;
}

inline rg_INT Polygon2D::num_vertices() const
{
	return m_numberOfVertices;
}

inline rg_INT Polygon2D::num_edges() const
{
	return m_numberOfVertices; // numOfEdges is equal to numOfVertices
}

inline rg_INT Polygon2D::num_holes() const
{
	return m_boundaryVertices.size() - 1;
}

inline rg_INT Polygon2D::num_shells() const
{
    return m_boundaryVertices.size();
}

inline bool Polygon2D::this_polygon_vertex_is_reflex(const list<rg_Point2D>& boundaryVertices, list<rg_Point2D>::const_iterator& i_vertex) const
{
    rg_Point2D vecCurrent2Previous(0.0, 0.0);
    rg_Point2D vecCurrent2Next(0.0, 0.0);

    list<rg_Point2D>::const_iterator i_lastVertex = prev(boundaryVertices.end());

    if (i_vertex == boundaryVertices.begin())
    {
        vecCurrent2Previous = (*i_lastVertex)                   - (*boundaryVertices.begin());
        vecCurrent2Next     = (*next(boundaryVertices.begin())) - (*boundaryVertices.begin());
    }
    else if (i_vertex == i_lastVertex)
    {
        vecCurrent2Previous = (*prev(i_lastVertex))       - (*i_lastVertex);
        vecCurrent2Next     = (*boundaryVertices.begin()) - (*i_lastVertex);
    }
    else
    {
        vecCurrent2Previous = (*prev(i_vertex)) - (*i_vertex);
        vecCurrent2Next     = (*next(i_vertex)) - (*i_vertex);
    }

    if (vecCurrent2Next * vecCurrent2Previous >= 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}



#endif
