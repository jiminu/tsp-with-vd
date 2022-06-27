#include "VDVertex.h"

#include "VDCell.h"
#include "VDEdge.h"
#include "VDPartialEdge.h"
#include "VDLoop.h"
#include "VDFace.h"


#include "rg_QTTetrahedron.h"
#include "BetaCell.h"
using namespace V::GeometryTier;


#include <set>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VDVertex::VDVertex()
: m_gates(rg_NULL), m_isOnInfinity(rg_UNKNOWN), m_radiusOfTangentSphere(0.0),
  m_ProbeTangibility( rg_UNKNOWN)
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++) {
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;

}
    

VDVertex::VDVertex( const rg_INT&     ID)
: TopologicalEntity(       ID ),
  m_gates(rg_NULL), m_isOnInfinity(rg_UNKNOWN), m_radiusOfTangentSphere(0.0),
  m_ProbeTangibility( rg_UNKNOWN)
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++) {
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;

}

VDVertex::VDVertex( const rg_FLAG&    onInfinity, 
                    const rg_Point3D& point, 
                    const rg_REAL&    radiusOfTS )
: m_gates(                 rg_NULL ), 
  m_isOnInfinity(          onInfinity ), 
  m_point(                 point ), 
  m_radiusOfTangentSphere( radiusOfTS ),
  m_ProbeTangibility( rg_UNKNOWN)
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)  {
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;
}


VDVertex::VDVertex(       Gate**      gates, 
                    const rg_FLAG&    onInfinity, 
                    const rg_Point3D& point, 
                    const rg_REAL&    radiusOfTS )
: m_isOnInfinity(          onInfinity ), 
  m_point(                 point ), 
  m_radiusOfTangentSphere( radiusOfTS ),
  m_ProbeTangibility( rg_UNKNOWN)
{
    m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        m_gates[i]         = gates[i];
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;
}


VDVertex::VDVertex( const rg_INT&     ID, 
                    const rg_FLAG&    onInfinity, 
                    const rg_Point3D& point, 
                    const rg_REAL&    radiusOfTS )
: TopologicalEntity(       ID ),
  m_gates(                 rg_NULL ), 
  m_isOnInfinity(          onInfinity ), 
  m_point(                 point ), 
  m_radiusOfTangentSphere( radiusOfTS ),
  m_ProbeTangibility( rg_UNKNOWN)
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++) {
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;
}


VDVertex::VDVertex( const rg_INT&     ID, 
                          Gate**      gates, 
                    const rg_FLAG&    onInfinity, 
                    const rg_Point3D& point, 
                    const rg_REAL&    radiusOfTS )
: TopologicalEntity(       ID ),
  m_isOnInfinity(          onInfinity ), 
  m_point(                 point ), 
  m_radiusOfTangentSphere( radiusOfTS ),
  m_ProbeTangibility( rg_UNKNOWN)
{
    m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        m_gates[i]         = gates[i];
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;
}


VDVertex::VDVertex( const TopologicalEntity& aTopoEntity, 
                    const rg_FLAG&           onInfinity, 
                    const rg_Point3D&        point, 
                    const rg_REAL&           radiusOfTS )
: TopologicalEntity(       aTopoEntity ),
  m_gates(                 rg_NULL ), 
  m_isOnInfinity(          onInfinity ), 
  m_point(                 point ), 
  m_radiusOfTangentSphere( radiusOfTS ),
  m_ProbeTangibility( rg_UNKNOWN)
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++) {
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;
}


VDVertex::VDVertex( const TopologicalEntity& aTopoEntity, 
                          Gate**             gates, 
                    const rg_FLAG&           onInfinity, 
                    const rg_Point3D&        point, 
                    const rg_REAL&           radiusOfTS )
: TopologicalEntity(       aTopoEntity ),
  m_isOnInfinity(          onInfinity ), 
  m_point(                 point ), 
  m_radiusOfTangentSphere( radiusOfTS ),
  m_ProbeTangibility( rg_UNKNOWN)
{
    if ( gates != rg_NULL )
    {
        m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
        for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            m_gates[i] = gates[i];
    }
    else
    {
        m_gates = rg_NULL;
    }
    
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++) {
        m_incidentEdges[i] = rg_NULL;
    }

    m_visited  = rg_FALSE;
    m_qCell    = rg_NULL;
    m_betaCell = rg_NULL;
}


VDVertex::VDVertex( const VDVertex& aVertex )
: TopologicalEntity(       aVertex ),
  m_isOnInfinity(          aVertex.m_isOnInfinity ), 
  m_point(                 aVertex.m_point ), 
  m_radiusOfTangentSphere( aVertex.m_radiusOfTangentSphere ),
  m_ProbeTangibility( aVertex.m_ProbeTangibility)
{
    if ( aVertex.m_gates != rg_NULL )
    {
        m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
        for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            m_gates[i] = aVertex.m_gates[i];
    }
    else
    {
        m_gates = rg_NULL;
    }

    m_visited  = aVertex.m_visited;
    m_qCell    = aVertex.m_qCell;
    m_betaCell = aVertex.m_betaCell;
}


VDVertex::~VDVertex()
{
    if ( m_gates != rg_NULL ) {
        delete [] m_gates;
    }

    //if ( m_qCell != rg_NULL ) {
    //    m_qCell->disconnectVVertex(this);
    //}
    //if ( m_betaCell != rg_NULL ) {
    //    m_betaCell->disconnectVVertex(this);
    //}
}






///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
VDEdge* VDVertex::getIncidentEdge( const rg_INT& i ) const
{
    if ( (i>-1) && (i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX) ) 
        return m_incidentEdges[i];
    else 
        return rg_NULL;
}

VDEdge** VDVertex::getAllIncidentEdges() 
{
    return m_incidentEdges;
}


Gate* VDVertex::getGate( const rg_INT& i ) const
{
    if ( m_gates == rg_NULL )
        return rg_NULL;

    if ( (i>-1) && (i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX) ) 
        return m_gates[i];
    else 
        return rg_NULL;
}

Gate** VDVertex::getAllGates() const
{
    return m_gates;
}


rg_FLAG VDVertex::isOnInfinity() const
{
    return m_isOnInfinity;
}


rg_Point3D VDVertex::getPoint() const
{
    return m_point;
}

rg_REAL VDVertex::getRadiusOfTangentSphere() const
{
    return m_radiusOfTangentSphere;
}

rg_FLAG VDVertex::isConnectedWithUnboundedEdge() const
{
    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        if ( !m_incidentEdges[i]->isBoundedEdge() )
            return rg_TRUE;
    }

    return rg_FALSE;
}

rg_FLAG VDVertex::isTangible() const
{
    return m_ProbeTangibility;
}

void    VDVertex::isTangible(const rg_FLAG& tangibility)
{
    m_ProbeTangibility = tangibility;
}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VDVertex::setIncidentEdge( const rg_INT& i, VDEdge* anEdge )
{
    m_incidentEdges[i] = anEdge;
}

void VDVertex::setAllIncidentEdge( VDEdge** incidentEdges )
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
        m_incidentEdges[i] = incidentEdges[i];
}

rg_FLAG VDVertex::setIncidentEdge( VDEdge* anEdge )
{
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        if ( m_incidentEdges[i] == rg_NULL )
        {
            m_incidentEdges[i] = anEdge;
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}


void VDVertex::setGate( const rg_INT& i, Gate* aGate )
{
    if ( m_gates == rg_NULL )  {
        m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
        for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            m_gates[i] = rg_NULL;
    }

    m_gates[i] = aGate;
}

void VDVertex::setAllGates( Gate** gates )
{
    if ( m_gates != rg_NULL )
        delete [] m_gates;

    m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
        m_gates[i] = gates[i];
}

rg_FLAG VDVertex::setGate( Gate* aGate )
{
    if ( m_gates == rg_NULL )  {
        m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
        for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            m_gates[i] = rg_NULL;
    }

    for (rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        if ( m_gates[i] == rg_NULL )
        {
            m_gates[i] = aGate;
            return rg_TRUE;
        }
    }

    return rg_FALSE;   
}

void VDVertex::isOnInfinity( const rg_FLAG& onInfinity )
{
    m_isOnInfinity = onInfinity;
}


void VDVertex::setPoint( const rg_Point3D& point )
{
    m_point = point;
}

void VDVertex::setRadiusOfTangentSphere( const rg_REAL& radiusOfTS )
{
    m_radiusOfTangentSphere = radiusOfTS;
}

void VDVertex::removeAllGates()
{
    if ( m_gates != rg_NULL )
        delete [] m_gates;

    m_gates = rg_NULL;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VDVertex& VDVertex::operator =(const VDVertex& aVertex )
{
    if ( this == &aVertex )
        return *this;

    m_ID = aVertex.m_ID;

    rg_INT i=0;
	for (i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
        m_incidentEdges[i] = aVertex.m_incidentEdges[i];


    if ( m_gates != rg_NULL )
        delete [] m_gates;

    if ( aVertex.m_gates != rg_NULL )
    {
        m_gates = new Gate*[NUM_INCIDENT_EDGES_OF_ONE_VERTEX];
        for (i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            m_gates[i] = aVertex.m_gates[i];
    }
    else
    {
        m_gates = rg_NULL;
    }

    m_isOnInfinity          = aVertex.m_isOnInfinity;
    m_point                 = aVertex.m_point;
    m_radiusOfTangentSphere = aVertex.m_radiusOfTangentSphere;

    m_visited  = aVertex.m_visited;
    m_qCell    = aVertex.m_qCell;
    m_betaCell = aVertex.m_betaCell;

    return *this;
}



///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..
void VDVertex::inquireAllCellsToDefineThisVertex(VDCell** cellToDefineVertex) const
{
    VDCell* cellsToDefineFirstEdge[3];
    VDCell* cellsToDefineSecondEdge[3];
    m_incidentEdges[0]->inquireIntoOnlyEdgeSharingCellsInCCW( cellsToDefineFirstEdge );
    m_incidentEdges[1]->inquireIntoOnlyEdgeSharingCellsInCCW( cellsToDefineSecondEdge );
    
	rg_INT i=0;
    for(i=0; i<3; i++)
    {
        cellToDefineVertex[i] = cellsToDefineFirstEdge[i];
    }

    rg_FLAG isInArray = rg_FALSE;
    for (i=0; i<3; i++)
    {
        isInArray = rg_FALSE;
        for (rg_INT j=0; j<3; j++)
        {
            if ( cellsToDefineSecondEdge[i] == cellToDefineVertex[j]  )
            {
                isInArray = rg_TRUE;
                break;
            }
        }

        if ( !isInArray )
        {
            cellToDefineVertex[3] = cellsToDefineSecondEdge[i];
            break;
        }
    }
}


//  topological quires..
void  VDVertex::inquireNeighborVertices(rg_dList<VDVertex*>& vertexList) const
{
    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++ )  {    
        if ( this == m_incidentEdges[i]->getStartVertex() )
            vertexList.add( m_incidentEdges[i]->getEndVertex() );
        else
            vertexList.add( m_incidentEdges[i]->getStartVertex() );
    }
}



void  VDVertex::inquireIncidentEdges(rg_dList<VDEdge*>& edgeList) const
{
    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++ )  {
        if ( m_incidentEdges[i] != rg_NULL ) {
            edgeList.add( m_incidentEdges[i] );
        }
    }
}



void  VDVertex::inquireIncidentFaces(rg_dList<VDFace*>& faceList) const
{
    m_incidentEdges[0]->inquireIncidentFaces(faceList);

    rg_dList<VDFace*> faceListOfSecondEdge;
    m_incidentEdges[1]->inquireIncidentFaces(faceListOfSecondEdge);

    faceListOfSecondEdge.reset4Loop();
    while( faceListOfSecondEdge.setNext4Loop() )  {
        faceList.addWithoutSame( faceListOfSecondEdge.getEntity() );
    }

    rg_dList<VDFace*> faceListOfThirdEdge;
    m_incidentEdges[2]->inquireIncidentFaces(faceListOfThirdEdge);
    faceListOfThirdEdge.reset4Loop();
    while( faceListOfThirdEdge.setNext4Loop() )  {
        faceList.addWithoutSame( faceListOfThirdEdge.getEntity() );
    }
}



void  VDVertex::inquireIncidentCells(rg_dList<VDCell*>& cellList) const 
{
    m_incidentEdges[0]->inquireIncidentCells( cellList );

    rg_dList<VDCell*> cellListOfSecondEdge;
    m_incidentEdges[1]->inquireIncidentCells( cellListOfSecondEdge );

    cellListOfSecondEdge.reset4Loop();
    while( cellListOfSecondEdge.setNext4Loop() )  {
        cellList.addWithoutSame( cellListOfSecondEdge.getEntity() );
    }
}


    
rg_BOOL VDVertex::isBoundingCell( VDCell* currCell ) const
{
    rg_dList<VDCell*> cellList;
    inquireIncidentCells( cellList );

    if ( cellList.isInList( currCell ) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL VDVertex::isIncidentTo(   VDEdge* currEdge ) const
{
    if ( currEdge->getStartVertex() == this ) {
        return rg_TRUE;
    }
    else if ( currEdge->getEndVertex() == this ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}


    
rg_BOOL VDVertex::hasTwinVertex() const
{
    rg_dList<VDVertex*> neighborVtx;    
    inquireNeighborVertices( neighborVtx );

    rg_BOOL isThereTwinVertex = rg_FALSE;
    VDVertex** neighborVtxArray = neighborVtx.getArray();

    for ( rg_INT i=0; i<4 && !isThereTwinVertex; i++ ) {
        for ( rg_INT j=i+1; j<4; j++ ) {
            if ( neighborVtxArray[i] == neighborVtxArray[j] ) {
                isThereTwinVertex = rg_TRUE;
                break;
            }
        }
    }

    delete [] neighborVtxArray;

    return isThereTwinVertex;
}



VDVertex* VDVertex::getTwinVertex() const
{
    rg_dList<VDVertex*> neighborVtx;    
    inquireNeighborVertices( neighborVtx );

    VDVertex*  twinVertex = rg_NULL;
    VDVertex** neighborVtxArray = neighborVtx.getArray();

    for ( rg_INT i=0; i<4; i++ ) {
        for ( rg_INT j=i+1; j<4; j++ ) {
            if ( neighborVtxArray[i] == neighborVtxArray[j] ) {
                twinVertex = neighborVtxArray[i];
                break;
            }
        }

        if ( twinVertex != rg_NULL ) {
            break;
        }
    }

    delete [] neighborVtxArray;

    return twinVertex;
}



void    VDVertex::getMateCellForIncidentEdges(rg_dList<VDCell*>& cellList) const
{
    for ( rg_INT i=0; i<4; i++ ) {
        VDCell* mateCell = rg_NULL;
        if ( this == m_incidentEdges[i]->getStartVertex() ) {
            mateCell = m_incidentEdges[i]->getMateCellOfStartVertex();
        }
        else {
            mateCell = m_incidentEdges[i]->getMateCellOfEndVertex();
        }

        cellList.add( mateCell );
    }
}




//rg_INT VDVertex::searchIncidentCells(rg_dList<VDCell*>& cellList) const
//{
//
//}



rg_INT VDVertex::getNumCellsToDefineThisVertex() const
{
    set<VDCell*> cellsDefineThisVertex;

    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++ )  {
        if ( m_incidentEdges[i] == rg_NULL ) continue;

        VDPartialEdge* startPrEdge = m_incidentEdges[i]->getPartialEdge();
        VDPartialEdge* currPrEdge  = startPrEdge;
        do {
            VDLoop* currLoop = currPrEdge->getLoop();
            if ( currLoop != rg_NULL ) {
                VDFace* currFace = currLoop->getFace();

                if ( currFace->getRightCell() != rg_NULL ) {
                    cellsDefineThisVertex.insert( currFace->getRightCell() );
                }
                if ( currFace->getLeftCell() != rg_NULL ) {
                    cellsDefineThisVertex.insert( currFace->getLeftCell() );
                }
            }

            currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
        } while ( currPrEdge != rg_NULL && currPrEdge != startPrEdge );
    }


    return cellsDefineThisVertex.size();
}



rg_INT VDVertex::getNumIncidentEdges() const
{
    rg_INT numEdges = 0;
    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++ )  {
        if ( m_incidentEdges[i] != rg_NULL ) {
            numEdges++;
        }
    }

    return numEdges;
}



// function to connect with cell in IWDS or eIWDS
void VDVertex::connectQCell(QTTetrahedron* q_cell)
{
    m_qCell = q_cell;
}



void VDVertex::disconnectQCell(QTTetrahedron* q_cell)
{
    if ( m_qCell == q_cell ) {
        m_qCell = rg_NULL;
    }
}



void VDVertex::connectBetaCell(BetaCell* b_cell)
{
    m_betaCell = b_cell;
}



void VDVertex::disconnectBetaCell(BetaCell* b_cell)
{
    if ( m_betaCell == b_cell ) {
        m_betaCell = rg_NULL;
    }
}


