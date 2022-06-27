#include "Polygon2D.h"
#include "rg_Point2D.h"
#include "rg_GeoFunc.h"
#include <fstream>

Polygon2D::Polygon2D()
{
	m_boundaryVertex = rg_NULL;
	m_numVertices = 0;
}

Polygon2D::~Polygon2D()
{
	m_numVertices = 0;
	if (m_boundaryVertex != rg_NULL)
	{
		delete[] m_boundaryVertex;
		m_boundaryVertex = rg_NULL;
	}
}


Polygon2D::Polygon2D(list<rg_Point2D>& boundaryVertices)
{
	setBoundaryVertices(boundaryVertices);
}


Polygon2D::Polygon2D(const rg_Point2D* boundaryVertex, const rg_INT& numVertices)
{
	setBoundaryVertices(boundaryVertex, numVertices);
}

Polygon2D::Polygon2D(const Polygon2D& polygon)
{
	setBoundaryVertices(polygon.m_boundaryVertex, polygon.m_numVertices);
	m_boundingBox = polygon.m_boundingBox;
}

void Polygon2D::getBoundaryVertices(rg_Point2D*& boundaryVertex, rg_INT& numVertices) const
{
	numVertices = m_numVertices;
	boundaryVertex = new rg_Point2D[numVertices];
	for (rg_INT i = 0; i < m_numVertices; i++)
		boundaryVertex[i] = m_boundaryVertex[i];
}

void Polygon2D::getBoundaryVertices(list<rg_Point2D>& boundaryVertices) const
{
	for (rg_INT i = 0; i < m_numVertices; i++)
		boundaryVertices.push_back(m_boundaryVertex[i]);
}

void Polygon2D::setBoundaryVertices(const rg_Point2D* boundaryVertex, const rg_INT& numVertices)
{
	if (numVertices <= 0 || boundaryVertex == rg_NULL)
		return;

	if (m_numVertices > 0)
	{
		if (m_boundaryVertex != rg_NULL)
		{
			delete[] m_boundaryVertex;
		}
	}

	m_numVertices = numVertices;
	m_boundaryVertex = new rg_Point2D[m_numVertices];
	for (rg_INT i = 0; i < m_numVertices; i++)
		m_boundaryVertex[i] = boundaryVertex[i];

	update_bounding_box();
}


void Polygon2D::setBoundaryVertices(list<rg_Point2D>& boundaryVertices)
{
	rg_INT numVertices = boundaryVertices.size();
	rg_Point2D* boundaryVertex = new rg_Point2D[numVertices];

	rg_INT i = 0;
	list<rg_Point2D>::iterator i_boundaryVertex = boundaryVertices.begin();
	while (i_boundaryVertex != boundaryVertices.end())
	{
		boundaryVertex[i++] = *i_boundaryVertex;
		++i_boundaryVertex;
	}

	setBoundaryVertices(boundaryVertex, numVertices);

	if (boundaryVertex != rg_NULL)
		delete[] boundaryVertex;
}

void Polygon2D::update_bounding_box()
{
	m_boundingBox.reset();
	for (rg_INT i = 0; i < m_numVertices; i++)
		m_boundingBox.contain(m_boundaryVertex[i]);
}

void Polygon2D::readPolygonFile(const char* fileName)
{
	ifstream fin;
	fin.open(fileName);

	char* seps = " \t\n";
	char buffer[200];
	char* context = NULL;

	fin.getline(buffer, 200);

#if defined(WIN32) || defined(WIN64)
	rg_INT numVertices = atoi(strtok_s(buffer, seps, &context));

	rg_Point2D* boundaryVertices = new rg_Point2D[numVertices];

	for (rg_INT i = 0; i < numVertices; i++)
	{
		fin.getline(buffer, 200);
		rg_REAL x = atof(strtok_s(buffer, seps, &context));
		rg_REAL y = atof(strtok_s(NULL, seps, &context));
		boundaryVertices[i].setPoint(x, y);

		rg_Point2D point(x, y);
		m_boundingBox.contain(point);
	}
#else
	rg_INT numVertices = atoi(strtok_r(buffer, seps, &context));

	rg_Point2D* boundaryVertices = new rg_Point2D[numVertices];

	for (rg_INT i = 0; i < numVertices; i++)
	{
		fin.getline(buffer, 200);
		rg_REAL x = atof(strtok_r(buffer, seps, &context));
		rg_REAL y = atof(strtok_r(NULL, seps, &context));
		boundaryVertices[i].setPoint(x, y);

		rg_Point2D point(x, y);
		m_boundingBox.contain(point);
	}
#endif	
	
	rg_Point2D minPt = m_boundingBox.getMinPt();
	rg_Point2D maxPt = m_boundingBox.getMaxPt();
	rg_REAL deltaX = minPt.getX();
	rg_REAL deltaY = minPt.getY();

	//for (rg_INT i = 0; i < numVertices; i++)
	//{
	//	boundaryVertices[i].setX(boundaryVertices[i].getX() - deltaX);
	//	boundaryVertices[i].setY(boundaryVertices[i].getY() - deltaY);
	//}

	//m_boundingBox.setMinPt(rg_Point2D(minPt.getX() - deltaX, minPt.getY() - deltaY));
	//m_boundingBox.setMaxPt(rg_Point2D(maxPt.getX() - deltaX, maxPt.getY() - deltaY));

	setBoundaryVertices(boundaryVertices, numVertices);

	//minPt = m_boundingBox.getMinPt();
	//maxPt = m_boundingBox.getMaxPt();
	//ofstream fout;
	//fout.open("D://Development//data//polygon//BBOX_info.txt");
	//fout << "  MinX  " << minPt.getX() << "  MinY  " << minPt.getY() << endl;
	//fout << "  MaxX  " << maxPt.getX() << "  MaxY  " << maxPt.getY() << endl;
	//fout.close();
}


void Polygon2D::boundary_vertices(list<rg_Point2D>& boundaryVertices) const
{
	for (int i = 0; i < m_numVertices; i++)
	{
		boundaryVertices.push_back(m_boundaryVertex[i]);
	}
}


rg_REAL Polygon2D::findMinEdgeLength() const
{
	//rg_INT i = 0;
	//rg_INT j = i + 1;
	//rg_REAL minEdgeLength = DBL_MAX;
	//while (j < m_numVertices)
	//{
	//	rg_REAL currentEdgeLength = m_boundaryVertex[i].distance(m_boundaryVertex[j]);
	//	if (currentEdgeLength < minEdgeLength)
	//		minEdgeLength = currentEdgeLength;
	//	++i;
	//	++j;
	//}

	rg_REAL minEdgeLength = DBL_MAX;
	for (int i = 0; i < m_numVertices-1; i++)
	{
		rg_REAL currentEdgeLength = m_boundaryVertex[i].distance(m_boundaryVertex[i+1]);
		if (currentEdgeLength < minEdgeLength)
			minEdgeLength = currentEdgeLength;
	}

	rg_REAL currentEdgeLength = m_boundaryVertex[m_numVertices - 1].distance(m_boundaryVertex[ 0 ]);
	if (currentEdgeLength < minEdgeLength)
		minEdgeLength = currentEdgeLength;

	return minEdgeLength;
}

bool Polygon2D::isReflexVertex(const rg_INT& vertexArrayIndex)
{
	rg_Point2D vecCurrent2Previous(0.0, 0.0);
	rg_Point2D vecCurrent2Next(0.0, 0.0);
	if (0 < vertexArrayIndex < m_numVertices - 1)
	{
		vecCurrent2Previous = m_boundaryVertex[vertexArrayIndex - 1] - m_boundaryVertex[vertexArrayIndex];
		vecCurrent2Next     = m_boundaryVertex[vertexArrayIndex + 1] - m_boundaryVertex[vertexArrayIndex];
	}
	else if (vertexArrayIndex == 0)
	{
		vecCurrent2Previous = m_boundaryVertex[m_numVertices - 1] - m_boundaryVertex[vertexArrayIndex];
		vecCurrent2Next     = m_boundaryVertex[vertexArrayIndex + 1] - m_boundaryVertex[vertexArrayIndex];
	}
	// vertexArrayIndex == m_numVertices - 1
	else
	{
		vecCurrent2Previous = m_boundaryVertex[m_numVertices - 1] - m_boundaryVertex[vertexArrayIndex];
		vecCurrent2Next     = m_boundaryVertex[ 0 ]               - m_boundaryVertex[vertexArrayIndex];
	}
	
	if (vecCurrent2Next * vecCurrent2Previous > 0)
		return false;
	else
		return true;
}

rg_REAL Polygon2D::get_length_of_perimeter() const
{
	rg_REAL perimeter = 0.0;
	for (int i = 0; i < m_numVertices - 1; i++)
	{
		rg_REAL currentEdgeLength = m_boundaryVertex[i].distance(m_boundaryVertex[i + 1]);
		perimeter += currentEdgeLength;
	}

	rg_REAL currentEdgeLength = m_boundaryVertex[m_numVertices - 1].distance(m_boundaryVertex[0]);
	perimeter += currentEdgeLength;

	return perimeter;
}


rg_REAL Polygon2D::get_area() const
{
	rg_Point2D anchorPoint(0.0, 0.0);
	for (rg_INT i = 0; i < m_numVertices; i++)
		anchorPoint += m_boundaryVertex[i];
	anchorPoint = anchorPoint / m_numVertices;

	rg_REAL area = 0.0;
	for (rg_INT i = 0; i < m_numVertices-1; i++)
	{
		area += rg_GeoFunc::computeSignedAreaOfTriangle(anchorPoint, m_boundaryVertex[i], m_boundaryVertex[i+1]);
	}
	area += rg_GeoFunc::computeSignedAreaOfTriangle(anchorPoint, m_boundaryVertex[m_numVertices - 1], m_boundaryVertex[0]);

	return area;
}

void Polygon2D::translate(const rg_Point2D& translationAmount)
{
	list<rg_Point2D> boundaryVertices;
	getBoundaryVertices(boundaryVertices);

	list<rg_Point2D>::iterator i_vertex = boundaryVertices.begin();
	for (; i_vertex != boundaryVertices.end(); ++i_vertex)
	{
		rg_Point2D currVertex = (*i_vertex);
		i_vertex->setPoint(currVertex + translationAmount);
	}

	setBoundaryVertices(boundaryVertices);
}

bool Polygon2D::disk_protrudes_boundary(const rg_Circle2D & disk, const rg_REAL& tolerance)
{
	rg_Point2D diskCenter = disk.getCenterPt();

	// Check whether disk center is either outside the polygon or on its boundary.
	for (rg_INT i = 0; i < m_numVertices-1; ++i)
	{
		rg_Line2D line1(m_boundaryVertex[i], m_boundaryVertex[i + 1]);
		if(rg_NPOS(line1.signed_distance(diskCenter), tolerance))
			return true;
	}

	rg_Line2D line1(m_boundaryVertex[m_numVertices-1], m_boundaryVertex[0]);
	if (rg_NPOS(line1.signed_distance(diskCenter), tolerance))
		return true;

	rg_REAL radius = disk.getRadius();

	// Test intersection between disk and each line (not line segment) which constitutes polygon boundary.
	for (rg_INT i = 0; i < m_numVertices - 1; ++i)
	{		
		rg_Line<rg_Point2D> line2(m_boundaryVertex[i], m_boundaryVertex[i + 1]);
		if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
			return true;
	}

	rg_Line<rg_Point2D> line2(m_boundaryVertex[m_numVertices - 1], m_boundaryVertex[0]);
	if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
		return true;

	return false;
}

bool Polygon2D::disk_is_outside_polygon(const rg_Circle2D & disk, const rg_REAL & tolerance)
{
	rg_Point2D diskCenter = disk.getCenterPt();

	// Check whether disk center is either inside the polygon or on its boundary.
	for (rg_INT i = 0; i < m_numVertices - 1; ++i)
	{
		rg_Line2D line1(m_boundaryVertex[i], m_boundaryVertex[i + 1]);
		if (rg_NNEG(line1.signed_distance(diskCenter), tolerance))
			return false;
	}

	rg_Line2D line1(m_boundaryVertex[m_numVertices - 1], m_boundaryVertex[0]);
	if (rg_NPOS(line1.signed_distance(diskCenter), tolerance))
		return false;

	rg_REAL radius = disk.getRadius();

	// Test intersection between disk and each line (not line segment) which constitutes polygon boundary.
	for (rg_INT i = 0; i < m_numVertices - 1; ++i)
	{
		rg_Line<rg_Point2D> line2(m_boundaryVertex[i], m_boundaryVertex[i + 1]);
		if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
			return false;
	}

	rg_Line<rg_Point2D> line2(m_boundaryVertex[m_numVertices - 1], m_boundaryVertex[0]);
	if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
		return false;

	return true;
}


Polygon2D& Polygon2D::operator=(const Polygon2D& polygon)
{
	if (this != &polygon)
	{
		setBoundaryVertices(polygon.m_boundaryVertex, polygon.m_numVertices);
		m_boundingBox = polygon.m_boundingBox;
	}
	return *this;
}



void Polygon2D::setRectangle( const rg_Point2D& minPt, const rg_Point2D& maxPt )
{
    m_numVertices = 4;
    
    if (m_boundaryVertex != NULL )
        delete[] m_boundaryVertex;


    double minX = minPt.getX();
    double minY = minPt.getY();
    double maxX = maxPt.getX();
    double maxY = maxPt.getY();

    m_boundaryVertex = new rg_Point2D[4];
    m_boundaryVertex[0] = rg_Point2D(minX, maxY);
    m_boundaryVertex[1] = rg_Point2D(minX, minY);
    m_boundaryVertex[2] = rg_Point2D(maxX, minY);
    m_boundaryVertex[3] = rg_Point2D(maxX, maxY);

    update_bounding_box();
}


