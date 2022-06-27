#include "rg_Const.h"

#include "rg_QTTetrahedron.h"
#include "rg_QTVertex.h"
#include "rg_QTEdge.h"
#include "rg_QTFace.h"
using namespace V::GeometryTier;


QTFace::QTFace()
: m_tetrahedron(rg_NULL), m_mateVertex(rg_NULL)
{
}



QTFace::QTFace(QTTetrahedron* tetrahedron, QTVertex* vertex)
: m_tetrahedron(tetrahedron), m_mateVertex(vertex)
{
}



QTFace::QTFace(const QTFace& face)
: m_tetrahedron(face.m_tetrahedron), m_mateVertex(face.m_mateVertex)
{
}



QTFace::~QTFace()
{
}




QTTetrahedron* QTFace::getTetrahedron() const
{
	return m_tetrahedron;
}



QTVertex* QTFace::getMateVertex() const
{
	return m_mateVertex;
}



QTFace QTFace::getTwinFace() const
{
    QTVertex* vertex[3] = {rg_NULL, rg_NULL, rg_NULL};
    inquireBoundingVerticesInThisWorld( vertex );

    QTTetrahedron* mateTetrahedron = m_tetrahedron->getMateTetrahedron(m_mateVertex);
    QTVertex**     neighborVertex  = mateTetrahedron->getVertices();

    QTVertex*      mateVertexOfMateTetrahedron = rg_NULL;
    for ( rg_INT i=0; i<4; i++ )  {
        if (    neighborVertex[i] != vertex[0] 
             && neighborVertex[i] != vertex[1] 
             && neighborVertex[i] != vertex[2] ) {
            mateVertexOfMateTetrahedron = neighborVertex[i];
            break;
        }
    }

    return QTFace(mateTetrahedron, mateVertexOfMateTetrahedron);
}



QTTetrahedron* QTFace::getIncidentTetrahedron() const
{
    return m_tetrahedron->getMateTetrahedron(m_mateVertex);
}



void QTFace::setTetrahedron(QTTetrahedron* tetrahedron)
{
	m_tetrahedron = tetrahedron;
}



void QTFace::setMateVertex(QTVertex* vertex)
{
	m_mateVertex = vertex;
}



void QTFace::setQTFace(QTTetrahedron* tetrahedron, QTVertex* vertex)
{
	m_tetrahedron = tetrahedron;
	m_mateVertex  = vertex;
}




QTFace& QTFace::operator =(const QTFace& face)
{
	if ( this == &face )
		return *this;

	m_tetrahedron = face.m_tetrahedron;
	m_mateVertex  = face.m_mateVertex;

	return *this;
}



rg_FLAG QTFace::operator==(const QTFace& face)
{
    if ( m_tetrahedron == face.m_tetrahedron && m_mateVertex == face.m_mateVertex )
        return rg_TRUE;
    else  {
        QTFace twinFace = face.getTwinFace();

        if (    twinFace.m_tetrahedron == face.m_tetrahedron 
             && twinFace.m_mateVertex  == face.m_mateVertex )
            return rg_TRUE;
        else 
            return rg_FALSE;
    }
}




rg_BOOL QTFace::isThere(const QTEdge& edge) const
{
    if ( m_tetrahedron->isThere( edge.getStartVertex(), edge.getEndVertex() ) ) {
        if (    m_mateVertex != edge.getStartVertex() 
             && m_mateVertex != edge.getEndVertex()   ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_FLAG QTFace::getEdgeOrientation(const QTEdge& edge) const
{
    if ( isThere( edge ) ) {
        rg_FLAG   edgeOrientation = rg_FALSE;

        QTVertex* startVtx = edge.getStartVertex();
        QTVertex* endVtx   = edge.getEndVertex();

        QTVertex* vertex[3] = {rg_NULL, rg_NULL, rg_NULL};
        inquireBoundingVerticesInThisWorld( vertex );

        if ( vertex[0] == startVtx ) {
            if ( vertex[1] == endVtx ) {
                edgeOrientation = rg_TRUE;
            }
        }
        else if ( vertex[1] == startVtx ) {
            if ( vertex[2] == endVtx ) {
                edgeOrientation = rg_TRUE;
            }
        }
        else { // ( vertex[2] == startVtx ) {
            if ( vertex[0] == endVtx ) {
                edgeOrientation = rg_TRUE;
            }
        }
        
        return edgeOrientation;
    }
    else {
        return rg_UNKNOWN;
    }
}


    /////  For Queres  /////
void QTFace::inquireBoundingVerticesInThisWorld(QTVertex** vertexArray) const
{
    if ( m_tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<3; i++ )  {
            vertexArray[i] = m_tetrahedron->getVertex(i);
        }
    }
    else  {
        QTVertex** vertices = m_tetrahedron->getVertices();

        if ( m_mateVertex == vertices[0] )  {
            vertexArray[0] = vertices[1];
            vertexArray[1] = vertices[2];
            vertexArray[2] = vertices[3];
        }
        else if  ( m_mateVertex == vertices[1] )  {
            vertexArray[0] = vertices[0];
            vertexArray[1] = vertices[3];
            vertexArray[2] = vertices[2];
        }
        else if  ( m_mateVertex == vertices[2] )  {
            vertexArray[0] = vertices[0];
            vertexArray[1] = vertices[1];
            vertexArray[2] = vertices[3];
        }
        else  {
            vertexArray[0] = vertices[0];
            vertexArray[1] = vertices[2];
            vertexArray[2] = vertices[1];
        }
    }
}



void QTFace::inquireBoundingVerticesInThisWorld(rg_dList<QTVertex*>& vertexList) const
{
    if ( m_tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<3; i++ )  
            vertexList.add( m_tetrahedron->getVertex(i) );
    }
    else  {
        QTVertex** vertices = m_tetrahedron->getVertices();

        if ( m_mateVertex == vertices[0] )  {
            vertexList.add(vertices[1]);
            vertexList.add(vertices[2]);
            vertexList.add(vertices[3]);
        }
        else if  ( m_mateVertex == vertices[1] )  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[3]);
            vertexList.add(vertices[2]);
        }
        else if  ( m_mateVertex == vertices[2] )  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[1]);
            vertexList.add(vertices[3]);
        }
        else  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[2]);
            vertexList.add(vertices[1]);
        }
    }
}



void QTFace::inquireBoundingEdgesInThisWorld(rg_dList<QTEdge>& edgeList) const
{
    rg_dList<QTVertex*> vertexList;
    inquireBoundingVerticesInThisWorld(vertexList);
    QTVertex** vertices = vertexList.getArray();

    edgeList.add( QTEdge(m_tetrahedron, vertices[0], vertices[1]) );
    edgeList.add( QTEdge(m_tetrahedron, vertices[1], vertices[2]) );
    edgeList.add( QTEdge(m_tetrahedron, vertices[2], vertices[0]) );

    delete [] vertices;
}



//void QTFace::inquireIncidentFacesInThisWorld(rg_dList<QTFace>& faceList)
//{
//}



void QTFace::inquireIncidentTetrahedraInThisWorld(rg_dList<QTTetrahedron*>& tetrahedronList) const
{
    tetrahedronList.add( m_tetrahedron );
    tetrahedronList.add( m_tetrahedron->getMateTetrahedron(m_mateVertex) );
}



