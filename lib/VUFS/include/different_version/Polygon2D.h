#ifndef POLYGON2D_H
#define POLYGON2D_H

#include "rg_Const.h"
#include "rg_BoundingBox2D.h"
class rg_Point2D;
#include <string>
#include <list>
using namespace std;

class Polygon2D
{
private:
	rg_Point2D* m_boundaryVertex;
	rg_INT      m_numVertices;
	rg_BoundingBox2D m_boundingBox;

public:
	Polygon2D();
	~Polygon2D();
	Polygon2D(list<rg_Point2D>& boundaryVertices);
	Polygon2D(const rg_Point2D* m_boundaryVertex, const rg_INT& numVertices);
	Polygon2D(const Polygon2D& polygon);
	
	
	void getBoundaryVertices(rg_Point2D*& boundaryVertex, rg_INT& numVertices) const;
	void getBoundaryVertices(list<rg_Point2D>& boundaryVertices) const;
	inline void getBoundingBox(rg_BoundingBox2D& boundingBox)
	{
		boundingBox = m_boundingBox;
	}
	inline rg_BoundingBox2D getBoundingBox() const { return m_boundingBox; }
	
	inline const rg_Point2D* boundary_vertices() const
	{
		return m_boundaryVertex;
	}
	inline const rg_INT num_vertices() const
	{
		return m_numVertices;
	}
	void boundary_vertices(list<rg_Point2D>& boundaryVertices) const;

	void setBoundaryVertices(const rg_Point2D* boundaryVertex, const rg_INT& numVertices);
	void setBoundaryVertices(list<rg_Point2D>& boundaryVertices);
	void update_bounding_box();


	void readPolygonFile(const char* fileName);



	rg_REAL findMinEdgeLength() const;
	bool isReflexVertex(const rg_INT& vertexArrayIndex);

	rg_REAL get_length_of_perimeter() const;

	rg_REAL get_area() const;

	void translate(const rg_Point2D& translationAmount);

	bool disk_protrudes_boundary(const rg_Circle2D& disk, const rg_REAL& tolerance);

	bool disk_is_outside_polygon(const rg_Circle2D& disk, const rg_REAL& tolerance);

	Polygon2D& operator=(const Polygon2D& polygon);

    void setRectangle( const rg_Point2D& minPt, const rg_Point2D& maxPt );
};

#endif
