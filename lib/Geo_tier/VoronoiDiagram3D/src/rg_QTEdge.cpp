#include "rg_Const.h"

#include "rg_QTEdge.h"
#include "rg_QTTetrahedron.h"
#include "rg_QTVertex.h"
using namespace V::GeometryTier;

QTEdge::QTEdge()
: m_tetrahedron(rg_NULL), m_startVertex(rg_NULL), m_endVertex(rg_NULL)
{
}



QTEdge::QTEdge(QTTetrahedron* tetrahedron, QTVertex* stVertex, QTVertex* edVertex)
: m_tetrahedron(tetrahedron), m_startVertex(stVertex), m_endVertex(edVertex)
{
}



QTEdge::QTEdge(const QTEdge& edge)
: m_tetrahedron(edge.m_tetrahedron), m_startVertex(edge.m_startVertex), m_endVertex(edge.m_endVertex)
{
}



QTEdge::~QTEdge()
{
}




QTTetrahedron* QTEdge::getTetrahedron() const
{
	return m_tetrahedron;
}



QTVertex*      QTEdge::getStartVertex() const
{
	return m_startVertex;
}



QTVertex*      QTEdge::getEndVertex() const
{
	return m_endVertex;
}



void QTEdge::setTetrahedron(QTTetrahedron* tetrahedron)
{
	m_tetrahedron = tetrahedron;
}



void QTEdge::setStartVertex(QTVertex* stVertex)
{
	m_startVertex = stVertex;
}



void QTEdge::setEndVertex(QTVertex* edVertex)
{
	m_endVertex = edVertex;
}



void QTEdge::setEdge(QTTetrahedron* tetrahedron, QTVertex* stVertex, QTVertex* edVertex)
{
	m_tetrahedron = tetrahedron;
	m_startVertex = stVertex;
	m_endVertex   = edVertex;
}



QTEdge& QTEdge::operator =(const QTEdge& edge)
{
	if ( this == &edge )
		return *this;

	m_tetrahedron = edge.m_tetrahedron;
	m_startVertex = edge.m_startVertex;
	m_endVertex   = edge.m_endVertex;

	return *this;
}



rg_FLAG QTEdge::operator==(const QTEdge& edge)
{
    if ( m_startVertex == edge.m_startVertex && m_endVertex == edge.m_endVertex)
        return rg_TRUE;
    else if ( m_startVertex == edge.m_endVertex && m_endVertex == edge.m_startVertex)
        return rg_TRUE;
    else
        return rg_FALSE;
}



QTFace QTEdge::getCCWNextFace(const QTFace& currFace) const
{
    QTFace ccwNextFace;

    rg_FLAG edgeOrientation = currFace.getEdgeOrientation( *this );
    if ( edgeOrientation == rg_TRUE ) {
        QTVertex* vertexOfFace[3] = {rg_NULL, rg_NULL, rg_NULL};
        currFace.inquireBoundingVerticesInThisWorld(vertexOfFace);

        QTVertex* mateVertex = rg_NULL;
        for ( rg_INT i=0; i<3; i++ ) {
            if ( vertexOfFace[i] != m_startVertex && vertexOfFace[i] != m_endVertex ) {
                mateVertex = vertexOfFace[i];
                break;
            }
        }

        ccwNextFace.setQTFace( currFace.getTetrahedron(), mateVertex );
    }
    else if ( edgeOrientation == rg_FALSE ) {
        if ( currFace.getIncidentTetrahedron() != rg_NULL ) {
            QTFace twinFace = currFace.getTwinFace();

            QTVertex* vertexOfFace[3] = {rg_NULL, rg_NULL, rg_NULL};
            twinFace.inquireBoundingVerticesInThisWorld(vertexOfFace);

            QTVertex* mateVertex = rg_NULL;
            for ( rg_INT i=0; i<3; i++ ) {
                if ( vertexOfFace[i] != m_startVertex && vertexOfFace[i] != m_endVertex ) {
                    mateVertex = vertexOfFace[i];
                    break;
                }
            }

            ccwNextFace.setQTFace( twinFace.getTetrahedron(), mateVertex );
        }
    }
    else { //if ( edgeOrientation == rg_UNKNOWN ) {
    }

    return ccwNextFace;
}



QTFace QTEdge::getCWNextFace( const QTFace& currFace) const
{
    QTFace ccwNextFace;

    rg_FLAG edgeOrientation = currFace.getEdgeOrientation( *this );
    if ( edgeOrientation == rg_TRUE ) {
        if ( currFace.getIncidentTetrahedron() != rg_NULL ) {
            QTFace twinFace = currFace.getTwinFace();

            QTVertex* vertexOfFace[3] = {rg_NULL, rg_NULL, rg_NULL};
            twinFace.inquireBoundingVerticesInThisWorld(vertexOfFace);

            QTVertex* mateVertex = rg_NULL;
            for ( rg_INT i=0; i<3; i++ ) {
                if ( vertexOfFace[i] != m_startVertex && vertexOfFace[i] != m_endVertex ) {
                    mateVertex = vertexOfFace[i];
                    break;
                }
            }

            ccwNextFace.setQTFace( twinFace.getTetrahedron(), mateVertex );
        }
    }
    else if ( edgeOrientation == rg_FALSE ) {
        QTVertex* vertexOfFace[3] = {rg_NULL, rg_NULL, rg_NULL};
        currFace.inquireBoundingVerticesInThisWorld(vertexOfFace);

        QTVertex* mateVertex = rg_NULL;
        for ( rg_INT i=0; i<3; i++ ) {
            if ( vertexOfFace[i] != m_startVertex && vertexOfFace[i] != m_endVertex ) {
                mateVertex = vertexOfFace[i];
                break;
            }
        }

        ccwNextFace.setQTFace( currFace.getTetrahedron(), mateVertex );
    }
    else { //if ( edgeOrientation == rg_UNKNOWN ) {
    }

    return ccwNextFace;
}



QTTetrahedron* QTEdge::getCCWNextCell(QTTetrahedron* currCell) const
{
    QTTetrahedron* ccwNextCell = rg_NULL;

    if ( currCell == rg_NULL ) {
        if ( m_tetrahedron->isAnomaly() != DEGENERACY_ANOMALY )  {
            QTFace currFace(m_tetrahedron, m_startVertex);
    
            QTVertex* vertexOfFace[3] = {rg_NULL, rg_NULL, rg_NULL};
            currFace.inquireBoundingVerticesInThisWorld(vertexOfFace);

            if ( m_endVertex == vertexOfFace[0] ) {
                ccwNextCell = m_tetrahedron->getMateTetrahedron( vertexOfFace[2] );
            }
            else if ( m_endVertex == vertexOfFace[1] ) {
                ccwNextCell = m_tetrahedron->getMateTetrahedron( vertexOfFace[0] );
            }
            else if ( m_endVertex == vertexOfFace[2] ) {
                ccwNextCell = m_tetrahedron->getMateTetrahedron( vertexOfFace[1] );
            }
            else {
                //  ccwNextCell = rg_NULL;
            }
        }
    }
    else {
        QTVertex* mateVertex[2] = {rg_NULL, rg_NULL};
        rg_INT i_mateVtx = 0;
        rg_INT i=0;
        for ( i=0; i<4; i++ ) {
            QTVertex* vertex = currCell->getVertex(i);
            if ( vertex != m_startVertex && vertex != m_endVertex ) {
                mateVertex[i_mateVtx] = vertex;
                i_mateVtx++;
            }
        }

        QTFace faceIncidentToThisEdge( currCell, mateVertex[0] );
        if ( faceIncidentToThisEdge.getEdgeOrientation( *this ) == rg_FALSE ) {
            ccwNextCell = faceIncidentToThisEdge.getIncidentTetrahedron();
        }
        else {
            faceIncidentToThisEdge.setMateVertex( mateVertex[1] );
            ccwNextCell = faceIncidentToThisEdge.getIncidentTetrahedron();
        }
    }

    return ccwNextCell;
}



QTTetrahedron* QTEdge::getCWNextCell(QTTetrahedron* currCell) const
{
    QTTetrahedron* cwNextCell = rg_NULL;

    if ( currCell == rg_NULL ) {
        if ( m_tetrahedron->isAnomaly() != DEGENERACY_ANOMALY )  {
            QTFace currFace(m_tetrahedron, m_startVertex);
    
            QTVertex* vertexOfFace[3] = {rg_NULL, rg_NULL, rg_NULL};
            currFace.inquireBoundingVerticesInThisWorld(vertexOfFace);

            if ( m_endVertex == vertexOfFace[0] ) {
                cwNextCell = m_tetrahedron->getMateTetrahedron( vertexOfFace[1] );
            }
            else if ( m_endVertex == vertexOfFace[1] ) {
                cwNextCell = m_tetrahedron->getMateTetrahedron( vertexOfFace[2] );
            }
            else if ( m_endVertex == vertexOfFace[2] ) {
                cwNextCell = m_tetrahedron->getMateTetrahedron( vertexOfFace[0] );
            }
            else {
                //  cwNextCell = rg_NULL;
            }
        }
    }
    else {
        QTVertex* mateVertex[2] = {rg_NULL, rg_NULL};
        rg_INT i_mateVtx = 0;
        rg_INT i=0;
        for ( i=0; i<4; i++ ) {
            QTVertex* vertex = currCell->getVertex(i);
            if ( vertex != m_startVertex && vertex != m_endVertex ) {
                mateVertex[i_mateVtx] = vertex;
                i_mateVtx++;
            }
        }

        QTFace faceIncidentToThisEdge( currCell, mateVertex[0] );
        if ( faceIncidentToThisEdge.getEdgeOrientation( *this ) == rg_TRUE ) {
            cwNextCell = faceIncidentToThisEdge.getIncidentTetrahedron();
        }
        else {
            faceIncidentToThisEdge.setMateVertex( mateVertex[1] );
            cwNextCell = faceIncidentToThisEdge.getIncidentTetrahedron();
        }
    }

    return cwNextCell;
}



/////  For Queres  /////
void QTEdge::inquireBoundingVerticesInThisWorld(rg_dList<QTVertex*>& vertexList)
{
    vertexList.add(m_startVertex);
    vertexList.add(m_endVertex);
}



void QTEdge::inquireIncidentEdgesInThisWorld(rg_dList<QTEdge>& edgeList)
{
    if ( m_tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        QTVertex* thirdVertex = rg_NULL;
        for ( rg_INT i=0; i<3; i++ )  {
            QTVertex* currVertex = m_tetrahedron->getVertex(i);
            if ( currVertex != m_startVertex  && currVertex != m_endVertex )  {
                thirdVertex = currVertex;
                break;
            }
        }
        
        edgeList.add( QTEdge(m_tetrahedron, m_endVertex, thirdVertex) );
        edgeList.add( QTEdge(m_tetrahedron, thirdVertex, m_startVertex) );

        return;
    }

    
    QTTetrahedron* currTetrahedron = m_tetrahedron;
    do  {
        QTFace currFace(currTetrahedron, m_startVertex);
        rg_dList<QTVertex*> vertexList;
        currFace.inquireBoundingVerticesInThisWorld(vertexList);
        QTVertex** vertices = vertexList.getArray();

        if ( m_endVertex == vertices[0] )  {
            edgeList.add( QTEdge(currTetrahedron, m_endVertex, vertices[1]) );
            edgeList.add( QTEdge(currTetrahedron, vertices[1], m_startVertex) );

            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
        }
        else if ( m_endVertex == vertices[1] )  {
            edgeList.add( QTEdge(currTetrahedron, m_endVertex, vertices[1]) );
            edgeList.add( QTEdge(currTetrahedron, vertices[1], m_startVertex) );

            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
        }
        else  {
            edgeList.add( QTEdge(currTetrahedron, m_endVertex, vertices[1]) );
            edgeList.add( QTEdge(currTetrahedron, vertices[1], m_startVertex) );
            
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
        }

        delete [] vertices;
    } while ( m_tetrahedron != currTetrahedron );
}



void QTEdge::inquireIncidentFacesInThisWorld(rg_dList<QTFace>& faceList)
{
    if ( m_tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        faceList.add( QTFace(m_tetrahedron, rg_NULL) );        
        return;
    }

    
    QTTetrahedron* currTetrahedron = m_tetrahedron;
    do  {
        QTFace currFace(currTetrahedron, m_startVertex);
        rg_dList<QTVertex*> vertexList;
        currFace.inquireBoundingVerticesInThisWorld(vertexList);
        QTVertex** vertices = vertexList.getArray();

        if ( m_endVertex == vertices[0] )  {
            faceList.add( QTFace(currTetrahedron, vertices[2]) );        
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
        }
        else if ( m_endVertex == vertices[1] )  {
            faceList.add( QTFace(currTetrahedron, vertices[0]) );        
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
        }
        else  {
            faceList.add( QTFace(currTetrahedron, vertices[1]) );        
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
        }

        delete [] vertices;
    } while ( m_tetrahedron != currTetrahedron );
}



void QTEdge::inquireIncidentTetrahedraInThisWorld(rg_dList<QTTetrahedron*>& tetrahedronList)
{
    if ( m_tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )
        return;

    
    QTTetrahedron* currTetrahedron = m_tetrahedron;

    do  {
        tetrahedronList.add(currTetrahedron);

        QTFace currFace(currTetrahedron, m_startVertex);
        rg_dList<QTVertex*> vertexList;
        currFace.inquireBoundingVerticesInThisWorld(vertexList);
        QTVertex** vertices = vertexList.getArray();

        if ( m_endVertex == vertices[0] )  
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
        else if ( m_endVertex == vertices[1] )  
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
        else  
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);

        delete [] vertices;
    } while ( m_tetrahedron != currTetrahedron );

}




