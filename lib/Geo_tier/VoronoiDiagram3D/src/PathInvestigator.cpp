#include "PathInvestigator.h"
using namespace V::GeometryTier;

#include "rg_RelativeOp.h"

#include <float.h>
//#include <time.h>
//#include <math.h>
//#include <stdlib.h> 

PathInvestigator::PathInvestigator()
: m_vertex(rg_NULL), m_currPathID(-1)
{
    for (int i=0; i<3; i++)
    {
        m_paths[i] = rg_NULL;
        m_orientation[i] = rg_FALSE;
    }
}

PathInvestigator::PathInvestigator(const PathInvestigator& anInvestigator)
: m_vertex(anInvestigator.m_vertex), m_currPathID( anInvestigator.m_currPathID )
{
    for (int i=0; i<3; i++)
    {
        m_paths[i] = anInvestigator.m_paths[i];
        m_orientation[i] = anInvestigator.m_orientation[i];
    }
}

PathInvestigator::~PathInvestigator()
{
}

VDVertex* PathInvestigator::getPropagatingVertex() const
{
    return m_vertex;
}

VDEdge*  PathInvestigator::getPath(const rg_INT& i) const
{
    if ( i<0 || i>2 )
        return rg_NULL;
    else
        return m_paths[i];
}

VDEdge** PathInvestigator::getPaths() 
{
    return m_paths;
}

rg_FLAG PathInvestigator::getOrientation(const rg_INT& i) const
{
    if ( i<0 || i>2 )
        return rg_NULL;
    else
        return m_orientation[i];
}

VDEdge*  PathInvestigator::getCurrPath() const
{
    if ( m_currPathID<0 || m_currPathID>2 )
        return rg_NULL;
    else
        return m_paths[ m_currPathID ];
}

rg_INT   PathInvestigator::getCurrPathID() const
{
     return m_currPathID;
}

rg_FLAG  PathInvestigator::getOrientationOfCurrPath() const
{
    return m_orientation[m_currPathID];
}


VDVertex* PathInvestigator::getNextPropagatingVertex() const
{
    if ( m_orientation[m_currPathID] )
        return m_paths[ m_currPathID ]->getEndVertex();
    else
        return m_paths[ m_currPathID ]->getStartVertex();
}

rg_FLAG PathInvestigator::isValidInvestigator() const
{
    if ( m_paths[0] == rg_NULL && m_paths[1] == rg_NULL && m_paths[2] == rg_NULL )
        return rg_FALSE;
    else 
        return rg_TRUE;
}


void PathInvestigator::setPropagatingVertex(VDVertex* pVertex)
{
    m_vertex = pVertex;
}

void PathInvestigator::setPath(const rg_INT& i, VDEdge* aPath, const rg_FLAG& orientation)
{
    if ( i<0 || i>2 )
        return;
    else
    {
        m_paths[ i ] = aPath; 
        m_orientation[i] = orientation;
    }
}

/*
void PathInvestigator::setPaths(VDEdge** paths)
{
    for (int i=0; i<3; i++)
        m_paths[i] = paths[i];
}

void PathInvestigator::setCurrPathID(const rg_INT& id)
{
    m_currPathID = id;
}

void PathInvestigator::setOrientationOfCurrPath(const rg_FLAG& orientation)
{
    m_orientationOfCurrPath = orientation;
}
*/
void PathInvestigator::setCurrPath(const rg_INT& id)
{
    if ( id<0 || id>2 )
        return;
    else
        m_currPathID = id;
}

/*
void PathInvestigator::setCurrPath(const rg_INT& id, VDEdge* aPath, const rg_FLAG& orientation)
{
    if ( id<0 || id>2 )
        return;
    else
    {
        m_currPathID = id;
        m_paths[ id ] = aPath; 
        m_orientation[id] = orientation;
    }
}
*/


//rg_FLAG PathInvestigator::makeNextInvestigator(const rg_Point3D& targetPt, PathInvestigator& nextInvestigator)
rg_FLAG PathInvestigator::makeNextInvestigator(const rg_Point3D& targetPt, rg_dList<VDVertex*>& verticesOnTrajectory, PathInvestigator& nextInvestigator)
{

    VDVertex* propagatingVertex = rg_NULL;
    if ( m_orientation[ m_currPathID ] == rg_TRUE )
        propagatingVertex = m_paths[ m_currPathID ]->getEndVertex();
    else
        propagatingVertex = m_paths[ m_currPathID ]->getStartVertex();

    VDVertex* currVertex = rg_NULL;
	verticesOnTrajectory.reset4Loop();
	while(verticesOnTrajectory.setNext4Loop())
	{
		currVertex = verticesOnTrajectory.getEntity();

        if ( propagatingVertex == currVertex )
            return rg_FALSE;
    }

    verticesOnTrajectory.addTail(propagatingVertex);
    nextInvestigator.m_vertex = propagatingVertex;
    rg_Point3D direction = targetPt - propagatingVertex->getPoint();


    rg_REAL angleOfEdges[4];
    rg_FLAG orientationOfEdges[4];
    rg_FLAG orderOfEdge[4];
	int i=0;
	for ( i=0; i<4; i++)
    {
        VDEdge* currEdge = propagatingVertex->getIncidentEdge(i);

        if ( currEdge->isBoundedEdge() && currEdge->isTangible() == FULLY_TANGIBLE && currEdge != m_paths[ m_currPathID ] )
        {
            if ( currEdge->getStartVertex() == propagatingVertex )
            {
                rg_Point3D edgeDirection = currEdge->getEndVertex()->getPoint() - currEdge->getStartVertex()->getPoint();
//                angleOfEdges[i]       = direction.angle(edgeDirection);
                angleOfEdges[i]       = currEdge->getEndVertex()->getPoint().distance(targetPt);
                orientationOfEdges[i] = rg_TRUE;
            }
            else
            {
                rg_Point3D edgeDirection = currEdge->getStartVertex()->getPoint() - currEdge->getEndVertex()->getPoint();
//                angleOfEdges[i] = direction.angle(edgeDirection);
                angleOfEdges[i]       = currEdge->getStartVertex()->getPoint().distance(targetPt);
                orientationOfEdges[i] = rg_FALSE;
            }
        }
        else
        {
            angleOfEdges[i] = DBL_MAX;
            orientationOfEdges[i] = -1;
        }
    }

    for ( i=0; i<4; i++ )
    {
        int order = 0;
        for ( int j=0; j<4; j++ )
        {
            if ( i==j )
                continue;

            if ( rg_LT( angleOfEdges[i], angleOfEdges[j] ) )
                continue;
            else
                order++;
        }
        orderOfEdge[i] = order;
    }

    for ( i=0; i<4; i++ )
    {
        if ( orderOfEdge[i] == 3 )
            continue;

        if ( orientationOfEdges[i] == -1 )
            continue;

        nextInvestigator.m_paths[ orderOfEdge[i] ]       = propagatingVertex->getIncidentEdge(i);
        nextInvestigator.m_orientation[ orderOfEdge[i] ] = orientationOfEdges[i];
    }    

    if ( nextInvestigator.m_paths[0] != NULL )
    {
        nextInvestigator.m_currPathID = 0;
        return rg_TRUE;
    }
    else
    {
        return rg_FALSE;
    }
}

void PathInvestigator::killCurrPath()
{
    m_paths[m_currPathID] = rg_NULL;
    m_orientation[m_currPathID] = -1;

    m_currPathID++;
}

PathInvestigator& PathInvestigator::operator =(const PathInvestigator& anInvestigator)
{
    if ( this == &anInvestigator )
        return *this;

    m_vertex     = anInvestigator.m_vertex;
    m_currPathID = anInvestigator.m_currPathID;

    for (int i=0; i<3; i++)
    {
        m_paths[i] = anInvestigator.m_paths[i];
        m_orientation[i] = anInvestigator.m_orientation[i];
    }


    return *this;
}

