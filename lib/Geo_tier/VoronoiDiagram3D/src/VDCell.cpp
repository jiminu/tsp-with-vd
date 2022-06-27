#include "VDCell.h"

#include "rg_BallGenerator.h"
#include "VDFace.h"

#include "rg_QTVertex.h"
#include "BetaVertex.h"
using namespace V::GeometryTier;


///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VDCell::VDCell()
: m_boundingFaces(), m_generator(rg_NULL), m_boundedness(rg_TRUE)
{
    m_qVertex    = rg_NULL;
    m_betaVertex = rg_NULL;
}

VDCell::VDCell(const rg_INT& ID)
: TopologicalEntity(ID), m_generator(rg_NULL), m_boundedness(rg_TRUE)
{
    m_qVertex    = rg_NULL;
    m_betaVertex = rg_NULL;
}

VDCell::VDCell(BallGenerator* generator)
: m_boundingFaces(), m_generator(generator), m_boundedness(rg_TRUE)
{
    m_qVertex    = rg_NULL;
    m_betaVertex = rg_NULL;
}

VDCell::VDCell(const rg_INT& ID, BallGenerator* generator)
: TopologicalEntity(ID), m_generator(generator), m_boundedness(rg_TRUE)
{
    m_qVertex    = rg_NULL;
    m_betaVertex = rg_NULL;
}

VDCell::VDCell(const TopologicalEntity& aTopoEntity, BallGenerator* generator)
: TopologicalEntity(aTopoEntity), m_generator(generator), m_boundedness(rg_TRUE)
{
    m_qVertex    = rg_NULL;
    m_betaVertex = rg_NULL;
}

VDCell::VDCell(const VDCell& aCell)
: TopologicalEntity(aCell), m_boundingFaces(aCell.m_boundingFaces), m_generator(aCell.m_generator), m_boundedness(aCell.m_boundedness)
{
    m_qVertex    = aCell.m_qVertex;
    m_betaVertex = aCell.m_betaVertex;
}


VDCell::~VDCell()
{
    //if ( m_qVertex != rg_NULL ) {
    //    m_qVertex->disconnectVCell(this);
    //}
    //if ( m_betaVertex != rg_NULL ) {
    //    m_betaVertex->disconnectVCell(this);
    //}
}

///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
BallGenerator* VDCell::getGenerator() const
{
    return m_generator;
}

rg_dList<VDFace*>* VDCell::getBoungindFaces() 
{
    return &m_boundingFaces;
}

rg_INT VDCell::getNumOfBoundingFaces() const
{
    return m_boundingFaces.getSize();
}

rg_FLAG VDCell::isBounded() const
{
    return m_boundedness;
}

void VDCell::isBounded( const rg_FLAG& boundedness)
{
    m_boundedness = boundedness;
}

///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VDCell::setGenerator(BallGenerator* generator)
{
    m_generator = generator;
    if ( m_generator != rg_NULL ) {
        m_generator->setCell( this );
    }
}

void VDCell::setBoundingFaces(const rg_dList<VDFace*>& boundingFaces)
{
    m_boundingFaces.duplicateList( boundingFaces );
}

void VDCell::addBoundingFace(VDFace* aBoundingFace)
{
    m_boundingFaces.addTail( aBoundingFace );
}

void V::GeometryTier::connectCellAndGenerator(VDCell* cell, BallGenerator* generator)
{
    cell->m_generator = generator;
    generator->setCell( cell );
}

///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VDCell& VDCell::operator =(const VDCell& aCell)
{
    if ( this == &aCell )
        return *this;

    m_ID        = aCell.m_ID;
    m_generator = aCell.m_generator;

    m_boundingFaces.duplicateList( aCell.m_boundingFaces );

    m_boundedness = aCell.m_boundedness;

    m_qVertex    = aCell.m_qVertex;
    m_betaVertex = aCell.m_betaVertex;

    return *this;
}


///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..
VDFace* VDCell::findFaceToShareWith(VDCell* anotherCell)
{
	VDFace* currentFace = rg_NULL;

	m_boundingFaces.reset4Loop();
	while( m_boundingFaces.setNext4Loop() )
	{
		currentFace = m_boundingFaces.getEntity();

        if ( currentFace->getLeftCell() == this )  
        {
            if ( currentFace->getRightCell() == anotherCell )  
                return currentFace;
        }
        else if ( currentFace->getLeftCell() == anotherCell )  
        {
            if ( currentFace->getRightCell() == this )  
                return currentFace;
        }
        else 
        {}
    }

    return rg_NULL;
}



rg_INT  VDCell::findFacesToShareWith(VDCell* anotherCell, rg_dList<VDFace*>& faceList) const
{
	m_boundingFaces.reset4Loop();
	while( m_boundingFaces.setNext4Loop() ) {
		VDFace* currFace = m_boundingFaces.getEntity();

        if ( currFace->getLeftCell() == anotherCell || currFace->getRightCell() == anotherCell ) {
            faceList.add( currFace );
        }
    }

    return faceList.getSize();
}


/*
///////////////////////////////////////////////////////////////////////////////
//
//  computations..
//
rg_REAL VDCell::computeVolume()
{
	VDFace* currFace = NULL;
	rg_dList<rg_Point3D>* faceMesh = NULL;
	rg_Point3D vec1, vec2, vec3;
	rg_Point3D center = ((BallGenerator*)m_generator)->getCenter();
	rg_REAL totalVolume = 0.;
	rg_REAL tetVolume = 0.;

	m_boundingFaces.reset4Loop();
	while(m_boundingFaces.setNext4Loop())
	{
		currFace = m_boundingFaces.getEntity();
		if(currFace->isOnInfinity())
			return rg_MAX_REAL;
		
		faceMesh = currFace->getFaceMesh();

		if(faceMesh == NULL)
		{
			currFace->makeMeshForVoronoiFace();
			faceMesh = currFace->getFaceMesh();

			//unbounded face case
			if(faceMesh == NULL)
				return rg_MAX_REAL;
		}

		faceMesh->reset4Loop();
		while(faceMesh->setNext4Loop())
		{
			//first threes are normals
			faceMesh->setNext4Loop();
			faceMesh->setNext4Loop();
			faceMesh->setNext4Loop();

			vec1 = faceMesh->getEntity() - center;
			faceMesh->setNext4Loop();
			vec2 = faceMesh->getEntity() - center;
			faceMesh->setNext4Loop();
			vec3 = faceMesh->getEntity() - center;
			
			tetVolume =  vec1.getX()*(vec2.getY()*vec3.getZ()-vec2.getZ()*vec3.getY())
					    -vec1.getY()*(vec2.getX()*vec3.getZ()-vec2.getZ()*vec3.getX())
					    +vec1.getZ()*(vec2.getX()*vec3.getY()-vec2.getY()*vec3.getX());
			
			if(tetVolume < 0)
				tetVolume = -tetVolume;

			totalVolume += tetVolume;
		}
	}
	return totalVolume/6.;
}



rg_REAL VDCell::computeArea()
{
	VDFace* currFace = NULL;
	rg_dList<rg_Point3D>* faceMesh = NULL;
	rg_Point3D vec1, vec2;
	rg_Point3D pt1, pt2, pt3;

	rg_REAL totalArea = 0.;
	rg_REAL triArea = 0.;

	m_boundingFaces.reset4Loop();
	while(m_boundingFaces.setNext4Loop())
	{
		currFace = m_boundingFaces.getEntity();
		if(currFace->isOnInfinity())
			return rg_MAX_REAL;
		
		faceMesh = currFace->getFaceMesh();

		if(faceMesh == NULL)
		{
			currFace->makeMeshForVoronoiFace();
			faceMesh = currFace->getFaceMesh();

			//unbounded face case
			if(faceMesh == NULL)
				return rg_MAX_REAL;
		}

		faceMesh->reset4Loop();
		while(faceMesh->setNext4Loop())
		{
			//first threes are normals
			faceMesh->setNext4Loop();
			faceMesh->setNext4Loop();
			faceMesh->setNext4Loop();

			pt1 = faceMesh->getEntity();
			faceMesh->setNext4Loop();
			pt2 = faceMesh->getEntity();
			faceMesh->setNext4Loop();
			pt3 = faceMesh->getEntity();
			
			vec1 = pt2 - pt1;
			vec2 = pt3 - pt2;
			triArea =  vec1.crossProduct(vec2).magnitude(); 
			
			totalArea += triArea;
		}
	}
	return totalArea/2.;
}
*/


///////////////////////////////////////////////////////////////////////////////
//
//  topological quires..
//
void  VDCell::inquireBoundingVertices(rg_dList<VDVertex*>& vertexList)
{
    VDFace* currFace = rg_NULL;

    m_boundingFaces.reset4Loop();
    while ( m_boundingFaces.setNext4Loop() )  {
        currFace = m_boundingFaces.getEntity();

        rg_dList<VDVertex*> vertexListOfCurrFace;
        currFace->inquireBoundingVertices(vertexListOfCurrFace);

        vertexListOfCurrFace.reset4Loop();
        while ( vertexListOfCurrFace.setNext4Loop() ) {
            vertexList.addWithoutSame( vertexListOfCurrFace.getEntity() );
        }
    }
}



void  VDCell::inquireBoundingEdges(rg_dList<VDEdge*>& edgeList)
{
    VDFace* currFace = rg_NULL;

    m_boundingFaces.reset4Loop();
    while ( m_boundingFaces.setNext4Loop() )  {
        currFace = m_boundingFaces.getEntity();

        rg_dList<VDEdge*> edgeListOfCurrFace;
        currFace->inquireBoundingEdges(edgeListOfCurrFace);

        edgeListOfCurrFace.reset4Loop();
        while ( edgeListOfCurrFace.setNext4Loop() ) {
            edgeList.addWithoutSame( edgeListOfCurrFace.getEntity() );
        }
    }
}



void  VDCell::inquireBoundingFaces(rg_dList<VDFace*>& faceList)
{
    faceList.duplicateList( m_boundingFaces );
}



void  VDCell::inquireIncidentCells(rg_dList<VDCell*>& cellList)
{
    VDFace* currFace = rg_NULL;

    m_boundingFaces.reset4Loop();
    while ( m_boundingFaces.setNext4Loop() )  {
        currFace = m_boundingFaces.getEntity();

        if ( this == currFace->getLeftCell() )
            cellList.add( currFace->getRightCell() );
        else 
            cellList.add( currFace->getLeftCell() );
    }
}




//  connect vertices in IWDS and eIWDS
void VDCell::connectQVertex(   QTVertex*   q_vtx)
{
    m_qVertex = q_vtx;
}



void VDCell::disconnectQVertex(QTVertex*   q_vtx)
{
    if ( m_qVertex == q_vtx ) {
        m_qVertex = rg_NULL;
    }
}



void VDCell::connectBetaVertex(BetaVertex* b_vtx)
{
    m_betaVertex = b_vtx;
}



void VDCell::disconnectBetaVertex(BetaVertex* b_vtx)
{
    if ( m_betaVertex == b_vtx ) {
        m_betaVertex = rg_NULL;
    }
}




