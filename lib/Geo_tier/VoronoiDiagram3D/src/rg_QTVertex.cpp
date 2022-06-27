#include "rg_Const.h"
#include "rg_QTVertex.h"
#include "rg_QTTetrahedron.h"
#include "VDCell.h"
using namespace V::GeometryTier;


QTVertex::QTVertex()
: m_firstTetrahedron(rg_NULL), m_ball(rg_NULL), m_visited(rg_FALSE)
{
    m_vCell = rg_NULL;
}



QTVertex::QTVertex(const rg_INT& id)
: TopologicalEntity(id), m_firstTetrahedron(rg_NULL), m_ball(rg_NULL), m_visited(rg_FALSE)
{
    m_vCell = rg_NULL;
}



QTVertex::QTVertex(const rg_INT& id, Ball* ball)
: TopologicalEntity(id), m_firstTetrahedron(rg_NULL), m_ball(ball), m_visited(rg_FALSE)
{
    m_vCell = rg_NULL;
}



QTVertex::QTVertex(const QTVertex& vertex)
: TopologicalEntity(vertex.m_ID), m_firstTetrahedron(vertex.m_firstTetrahedron), 
  m_ball(vertex.m_ball), m_visited(vertex.m_visited)
{
    m_vCell = vertex.m_vCell;
}



QTVertex::~QTVertex()
{
    m_firstTetrahedron = rg_NULL;
    m_ball             = rg_NULL;

    //if ( m_vCell != rg_NULL ) {
    //    m_vCell->disconnectQVertex(this);
    //}
}




QTTetrahedron* QTVertex::getFirstTetrahedron() const
{
    return m_firstTetrahedron;
}




Sphere QTVertex::getBall() const
{
    return m_ball->getGeometry();
}



void*  QTVertex::getProperty() const
{
    return m_ball->getProperty();
}



rg_INT QTVertex::getIDFromInput() const
{
    return m_ball->getIDFromInput();
}


    
Ball* QTVertex::getBallProperty()
{
    return m_ball;
}



rg_BOOL QTVertex::isVirtual() const
{
    if ( m_ball == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



void   QTVertex::setFirstTetrahedron(QTTetrahedron* tetrahedron)
{
    m_firstTetrahedron = tetrahedron;
}


    
void   QTVertex::setBall(Ball* ball)
{
    m_ball = ball;
}



// void   QTVertex::setBall(const Sphere& ball)
// {
//     m_ball = ball;
// }
// 
// 
// 
// void   QTVertex::setProperty(void* property)
// {
//     m_property = property;
// }




QTVertex& QTVertex::operator =(const QTVertex& vertex)
{
    if ( this == &vertex )
        return *this;

    m_ID               = vertex.m_ID;
    m_firstTetrahedron = vertex.m_firstTetrahedron;
    m_ball             = vertex.m_ball;
//     m_property         = vertex.m_property;
    m_visited          = vertex.m_visited;
    
    m_vCell = vertex.m_vCell;

    return *this;
}



void QTVertex::inquireIncidentTetrahedra(rg_dList<QTTetrahedron*>& incidentTetrahedra)
{
    incidentTetrahedra.add( m_firstTetrahedron );

    rg_INT i=0;
    rg_dNode<QTTetrahedron*>* currNode = incidentTetrahedra.getFirstpNode();    
    do  {
        QTTetrahedron** neighbor = currNode->getEntity()->getNeighbors();
        for ( i=0; i<4; i++ ) {
            if ( neighbor[i]->isThere(this) )
                incidentTetrahedra.addWithoutSame( neighbor[i] );
        }

        currNode = currNode->getNext();
    } while ( currNode != incidentTetrahedra.getFirstpNode() );
}



void QTVertex::connectVCell(VDCell* v_cell)
{
    m_vCell = v_cell;
}



void QTVertex::disconnectVCell(VDCell* v_cell)
{
    if ( m_vCell == v_cell ) {
        m_vCell = rg_NULL;
    }
}



