#include "rg_Const.h"
#include "rg_QTTetrahedron.h"
#include "rg_QTVertex.h"
#include "VDVertex.h"
using namespace V::GeometryTier;


QTTetrahedron::QTTetrahedron()
: m_emptyTangentSphere(), m_visited(rg_FALSE)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i] = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }

    m_vVertex = rg_NULL;
}



QTTetrahedron::QTTetrahedron(const rg_INT& id)
: TopologicalEntity(id), m_emptyTangentSphere(), m_visited(rg_FALSE)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i] = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }

    m_vVertex = rg_NULL;
}



QTTetrahedron::QTTetrahedron(const rg_INT& id, const Sphere& emptyTSphere)
: TopologicalEntity(id), m_emptyTangentSphere(emptyTSphere), m_visited(rg_FALSE)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i] = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }

    m_vVertex = rg_NULL;
}



QTTetrahedron::QTTetrahedron(const Sphere& emptyTSphere)
: m_emptyTangentSphere(emptyTSphere), m_visited(rg_FALSE)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i] = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }

    m_vVertex = rg_NULL;
}



QTTetrahedron::QTTetrahedron(const QTTetrahedron& tetrahedron)
: TopologicalEntity(tetrahedron.m_ID),  
  m_emptyTangentSphere(tetrahedron.m_emptyTangentSphere), m_visited(tetrahedron.m_visited)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i]   = tetrahedron.m_vertex[i];
        m_neighbor[i] = tetrahedron.m_neighbor[i];
    }

    m_vVertex = tetrahedron.m_vVertex;
}



QTTetrahedron::~QTTetrahedron()
{
    //if ( m_vVertex != rg_NULL ) {
    //    m_vVertex->disconnectQCell(this);
    //}
}




QTVertex* QTTetrahedron::getVertex(const rg_INDEX& i) const
{
    if ( i<0 || i>=4)
        return rg_NULL;
    else
        return m_vertex[i];
}



QTVertex** QTTetrahedron::getVertices()
{
    return m_vertex;
}



QTTetrahedron* QTTetrahedron::getNeighbor(const rg_INDEX& i) const
{
    if ( i<0 || i>=4)
        return rg_NULL;
    else
        return m_neighbor[i];
}



QTTetrahedron** QTTetrahedron::getNeighbors()
{
    return m_neighbor;
}



Sphere QTTetrahedron::getEmptyTangentSphere() const
{
    return m_emptyTangentSphere;
}




QTTetrahedron*  QTTetrahedron::getMateTetrahedron(const rg_INT& pos) const
{
    if ( pos<0 || pos>=4 )
        return rg_NULL;
    else
        return m_neighbor[pos];
}



QTTetrahedron*  QTTetrahedron::getMateTetrahedron(QTVertex* vertex) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( vertex == m_vertex[i] )
            return m_neighbor[i];
    }

    return rg_NULL;
}



rg_INT QTTetrahedron::getPosOfVertex(QTVertex* vertex) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( vertex == m_vertex[i] )
            return i;
    }

    return -1;
}



rg_INT QTTetrahedron::getMateVertexPosOfIncidentFace(QTTetrahedron* neighbor, const rg_INT& mateVertexPos)
{
	rg_INT i=0;
    QTVertex* incidentVertex[3];
    rg_INT    vIndex = 0;
    for ( i=0; i<4; i++ )  {
        if ( i == mateVertexPos )
            continue;

        incidentVertex[vIndex++] = neighbor->m_vertex[i];
    }

    vIndex = 0;
    rg_INT  thisMateVertexPos = -1;

    for ( i=0; i<4; i++ )  {
        rg_FLAG isIncident = rg_FALSE;
        for ( rg_INT j=0; j<3; j++)  {        
            if ( m_vertex[i] == incidentVertex[j] )  {
                isIncident = rg_TRUE;
                break;
            }
        }

        if ( isIncident == rg_FALSE )  {
            return i;
        }
    }

    return thisMateVertexPos;
}



rg_INT QTTetrahedron::getMatePosInNeighbor(QTTetrahedron* neighbor, const rg_INT& neighborPos) const
{
    if ( m_neighbor[neighborPos] != neighbor ) {
        return -1;
    }

    QTVertex* vertexOnIncidentFace[3];
    rg_INT    i_vtx = 0;
    rg_INT    i     = 0;
    for ( i=0; i<4; i++ )  {
        if ( i == neighborPos ) {
            continue;
        }

        vertexOnIncidentFace[i_vtx++] = m_vertex[i];
    }

    rg_INT matePosInNeighbor = -1;
    QTVertex** vertexOfNeighbor = neighbor->getVertices();
    for ( i=0; i<4; i++ ) {
        if (    vertexOfNeighbor[i] != vertexOnIncidentFace[0]
             && vertexOfNeighbor[i] != vertexOnIncidentFace[1]
             && vertexOfNeighbor[i] != vertexOnIncidentFace[2] ) {
            matePosInNeighbor = i;
            break;
        }
    }

    return matePosInNeighbor;
}


rg_INT QTTetrahedron::getPosOfNeighbor(QTTetrahedron* neighbor) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( neighbor == m_neighbor[i] )
            return i;
    }

    return -1;
}



QTTetrahedron* QTTetrahedron::getNeighborInMultiplicityAnomaly() const
{
    for ( rg_INT i=0; i<4; i++ )  {
        for ( rg_INT j=i+1; j<4; j++ )  {
            if ( m_neighbor[i] == m_neighbor[j] )  
                return m_neighbor[i];
        }
    }

    return rg_NULL;
}



rg_INT QTTetrahedron::getEdgesNotSharedWithNeighbor(QTTetrahedron* neighbor, QTEdge* edges) 
{
    if ( m_neighbor[0]==neighbor && m_neighbor[1]==neighbor )  {
        edges[0].setEdge(this, m_vertex[0], m_vertex[1]);
        return 1;
    }
    else if ( m_neighbor[0]==neighbor && m_neighbor[2]==neighbor )  {
        edges[0].setEdge(this, m_vertex[0], m_vertex[2]);
        return 1;
    }
    else if ( m_neighbor[0]==neighbor && m_neighbor[3]==neighbor )  {
        edges[0].setEdge(this, m_vertex[0], m_vertex[3]);
        return 1;
    }
    else if ( m_neighbor[1]==neighbor && m_neighbor[2]==neighbor )  {
        edges[0].setEdge(this, m_vertex[1], m_vertex[2]);
        return 1;
    }
    else if ( m_neighbor[1]==neighbor && m_neighbor[3]==neighbor )  {
        edges[0].setEdge(this, m_vertex[1], m_vertex[3]);
        return 1;
    }
    else if ( m_neighbor[2]==neighbor && m_neighbor[3]==neighbor )  {
        edges[0].setEdge(this, m_vertex[2], m_vertex[3]);
        return 1;
    }
    else if ( m_neighbor[0]==neighbor )  {
        edges[0].setEdge(this, m_vertex[0], m_vertex[1]);
        edges[1].setEdge(this, m_vertex[0], m_vertex[2]);
        edges[2].setEdge(this, m_vertex[0], m_vertex[3]);
        return 3;
    }
    else if ( m_neighbor[1]==neighbor )  {
        edges[0].setEdge(this, m_vertex[1], m_vertex[0]);
        edges[1].setEdge(this, m_vertex[1], m_vertex[2]);
        edges[2].setEdge(this, m_vertex[1], m_vertex[3]);
        return 3;
    }
    else if ( m_neighbor[2]==neighbor )  {
        edges[0].setEdge(this, m_vertex[2], m_vertex[0]);
        edges[1].setEdge(this, m_vertex[2], m_vertex[1]);
        edges[2].setEdge(this, m_vertex[2], m_vertex[3]);
        return 3;
    }
    else if ( m_neighbor[3]==neighbor )  {
        edges[0].setEdge(this, m_vertex[3], m_vertex[0]);
        edges[1].setEdge(this, m_vertex[3], m_vertex[1]);
        edges[2].setEdge(this, m_vertex[3], m_vertex[2]);
        return 3;
    }
    else {
        return 0;
    }
}



rg_BOOL QTTetrahedron::isVirtual() const
{
    if ( isIsolatedFace() ) {
        return rg_FALSE;
    }


    rg_BOOL isVirtualTetrahedron = rg_FALSE;
    for (rg_INT i=0; i<4; i++ )  {
        if ( m_vertex[i]->isVirtual() ) {
            isVirtualTetrahedron = rg_TRUE;
            break;
        }
    }

    return isVirtualTetrahedron;
}



rg_FLAG QTTetrahedron::isThere(QTVertex* vertex) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( vertex == m_vertex[i] )
            return rg_TRUE;
    }

    return rg_FALSE;
}



rg_FLAG QTTetrahedron::isThere(QTVertex* vertex1, QTVertex* vertex2) const
{
    if ( isThere(vertex1) && isThere(vertex2) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG QTTetrahedron::isThereThisNeighbor(QTTetrahedron* neighbor) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( neighbor == m_neighbor[i] )
            return rg_TRUE;
    }

    return rg_FALSE;
}



rg_FLAG QTTetrahedron::isInfiniteTetrahedron() const
{
    return isVirtual();
//     if ( m_vertex[3] == rg_NULL )
//         return rg_FALSE;
// 
//     if (    ( m_vertex[0]->getID() == 0 )
//          || ( m_vertex[1]->getID() == 0 )
//          || ( m_vertex[2]->getID() == 0 )
//          || ( m_vertex[3]->getID() == 0 ) )
//          return rg_TRUE;
//     else
//         return rg_FALSE;
}



rg_BOOL QTTetrahedron::isIsolatedFace() const
{
    if ( m_vertex[3] == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_FLAG QTTetrahedron::isAnomaly() const
{
    if ( m_vertex[3] == rg_NULL )
        return DEGENERACY_ANOMALY;

    /*
    if ( m_vertex[0]==m_vertex[1] || m_vertex[0]==m_vertex[2] || m_vertex[0]==m_vertex[3] )
        return MULTIPLICITY_ANOMALY;
    else if ( m_vertex[1]==m_vertex[2] || m_vertex[1]==m_vertex[3] )
        return MULTIPLICITY_ANOMALY;
    else if ( m_vertex[2]==m_vertex[3] )
        return MULTIPLICITY_ANOMALY;
    else 
        return rg_FALSE;
    */

    rg_INT multiplicity = 0;
    for ( rg_INT i=0; i<4; i++ )  {
        for ( rg_INT j=i+1; j<4; j++ )  {
            if ( m_neighbor[i] == m_neighbor[j] )  {
                multiplicity++;
            }
        }

        if ( multiplicity != 0 )
            break;
    }

    switch ( multiplicity )  {
    case 0:
        return rg_FALSE;
        break;
    case 1:
        return MULTIPLICITY_ANOMALY_BY_TWO_COMMON_FACE;
        break;
    case 2:
        return MULTIPLICITY_ANOMALY_BY_THREE_COMMON_FACE;
        break;
    case 3:
        return MULTIPLICITY_ANOMALY_BY_FOUR_COMMON_FACE;
        break;
    default:
        return rg_FALSE;
        break;        
    }
}



rg_INT QTTetrahedron::isMultiplicity(QTTetrahedron*& neighbor) const
{
    rg_INT multiplicity = 0;
    for ( rg_INT i=0; i<4; i++ )  {
        for ( rg_INT j=i+1; j<4; j++ )  {
            if ( m_neighbor[i] == m_neighbor[j] )  {
                multiplicity++;
                neighbor = m_neighbor[i];
            }
        }

        if ( multiplicity != 0 )
            return multiplicity+1;
    }

    return multiplicity;
}



void QTTetrahedron::setVertex(const rg_INDEX& i, QTVertex* vertex)
{
    if ( i<0 || i>=4)
        return;
    else
        m_vertex[i] = vertex;
}



void QTTetrahedron::setVertices(QTVertex* vertex1, QTVertex* vertex2, 
                                QTVertex* vertex3, QTVertex* vertex4)
{
    m_vertex[0] = vertex1;
    m_vertex[1] = vertex2;
    m_vertex[2] = vertex3;
    m_vertex[3] = vertex4;    
}



void QTTetrahedron::setVertices(QTVertex** vertices)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i] = vertices[i];
    }
}



void QTTetrahedron::setNeighbor(const rg_INDEX& i, QTTetrahedron* neighbor)
{
    if ( i<0 || i>=4)
        return;
    else
        m_neighbor[i] = neighbor;
}



void QTTetrahedron::setNeighbors(QTTetrahedron* neighbor1, QTTetrahedron* neighbor2, 
                                 QTTetrahedron* neighbor3, QTTetrahedron* neighbor4)
{
    m_neighbor[0] = neighbor1;
    m_neighbor[1] = neighbor2;
    m_neighbor[2] = neighbor3;
    m_neighbor[3] = neighbor4;
}



void QTTetrahedron::setNeighbors(QTTetrahedron** neighbors)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_neighbor[i] = neighbors[i];
    }
}



void QTTetrahedron::setEmptyTangentSphere(const Sphere& emptyTSphere)
{
    m_emptyTangentSphere = emptyTSphere;
}




QTTetrahedron& QTTetrahedron::operator =(const QTTetrahedron& tetrahedron)
{
    if ( this == &tetrahedron )
        return *this;

    m_ID      = tetrahedron.m_ID;
    m_visited = tetrahedron.m_visited;
    
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i]   = tetrahedron.m_vertex[i];
        m_neighbor[i] = tetrahedron.m_neighbor[i];
    }
    m_emptyTangentSphere = tetrahedron.m_emptyTangentSphere;

    return *this;
}



void QTTetrahedron::findCommonEdgesWithNeighbor(QTTetrahedron* neighbor, rg_dList<QTEdge>& edgeList) 
{
    rg_dList<QTFace> commonFaceList;
    findCommonFacesWithNeighbor(neighbor, commonFaceList);

    QTFace* currFace = rg_NULL;
    commonFaceList.reset4Loop();
    while ( commonFaceList.setNext4Loop() )  {
        currFace = commonFaceList.getpEntity();

        rg_dList<QTEdge> boundingEdges;
        currFace->inquireBoundingEdgesInThisWorld( boundingEdges );
        
        boundingEdges.reset4Loop();
        while ( boundingEdges.setNext4Loop() )  {
            QTEdge currEdge = boundingEdges.getEntity();

            edgeList.addWithoutSame( currEdge );
        }
    }
}



void QTTetrahedron::findCommonFacesWithNeighbor(QTTetrahedron* neighbor, rg_dList<QTFace>& faceList) 
{
    for ( rg_INT i=0; i<4; i++ )      {
        if ( m_neighbor[i] == neighbor )  {
            faceList.add( QTFace(this, m_vertex[i]) );
        }
    }    
}


void QTTetrahedron::inquireBoundingVerticesInThisWorld(rg_dList<QTVertex*>& vertexList)
{
    if ( isAnomaly() == DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<3; i++ )  
            vertexList.add( m_vertex[i] );
    }
    else  {
        for ( rg_INT i=0; i<4; i++ )  
            vertexList.add( m_vertex[i] );
    }
}



void QTTetrahedron::inquireBoundingEdgesInThisWorld(rg_dList<QTEdge>& edgeList)
{
    if ( isAnomaly() == DEGENERACY_ANOMALY )  {
        edgeList.add( QTEdge(this, m_vertex[0], m_vertex[1]) );
        edgeList.add( QTEdge(this, m_vertex[1], m_vertex[2]) );
        edgeList.add( QTEdge(this, m_vertex[2], m_vertex[0]) );
    }
    else  {
        edgeList.add( QTEdge(this, m_vertex[0], m_vertex[1]) );
        edgeList.add( QTEdge(this, m_vertex[0], m_vertex[2]) );
        edgeList.add( QTEdge(this, m_vertex[0], m_vertex[3]) );
        edgeList.add( QTEdge(this, m_vertex[1], m_vertex[2]) );
        edgeList.add( QTEdge(this, m_vertex[2], m_vertex[3]) );
        edgeList.add( QTEdge(this, m_vertex[3], m_vertex[1]) );
    }
}



void QTTetrahedron::inquireBoundingFacesInThisWorld(rg_dList<QTFace>& faceList)
{
    if ( isAnomaly() == DEGENERACY_ANOMALY )  {
        faceList.add( QTFace(this, m_vertex[3]) );        
    }
    else  {
        faceList.add( QTFace(this, m_vertex[0]) );
        faceList.add( QTFace(this, m_vertex[1]) );
        faceList.add( QTFace(this, m_vertex[2]) );
        faceList.add( QTFace(this, m_vertex[3]) );
    }
}



void QTTetrahedron::inquireIncidentTetrahedraInThisWorld(rg_dList<QTTetrahedron*>& tetrahedronList)
{
    if ( isAnomaly() != DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<4; i++ )  
            tetrahedronList.addWithoutSame( m_neighbor[i] );
    }
}



// function to connect with vertex in VD
void QTTetrahedron::connectVVertex(VDVertex* v_vtx)
{
    m_vVertex = v_vtx;
}

void QTTetrahedron::disconnectVVertex(VDVertex* v_vtx)
{
    if ( m_vVertex == v_vtx ) {
        m_vVertex = rg_NULL;
    }
}
