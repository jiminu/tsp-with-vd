#include "ManifoldizedBetaShape.h"
#include "MBSBody.h"
#include "rg_AugmentedBetaComplex.h"


ManifoldizedBetaShape::ManifoldizedBetaShape()
{
    m_betaUniverse = rg_NULL;
    m_betaValue = 0.0;
}


ManifoldizedBetaShape::~ManifoldizedBetaShape()
{
    removeAll();    
}

	
    //get function
BetaUniverse* ManifoldizedBetaShape::getBetaUniverse() const
{
    return m_betaUniverse;
}


rg_REAL ManifoldizedBetaShape::getBetaValue() const
{
    return m_betaValue;
}


rg_dList< MBSBody* >* ManifoldizedBetaShape::getBodies()
{
    return &m_bodies;
}


rg_INT ManifoldizedBetaShape::getNumOfBodies() const
{
    return m_bodies.getSize();
}


rg_INT ManifoldizedBetaShape::getNumOfShells() const
{
    rg_INT numOfShells = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        numOfShells++;
        numOfShells += currBody->getInteriorShells()->getSize();
    }

    return numOfShells;
}


rg_INT ManifoldizedBetaShape::getNumOfVertices() const
{
    rg_INT numOfVertices = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        numOfVertices += currBody->getExteriorShell()->getVertices()->getSize();

        rg_dList< MBSShell* >* currInteriorShells = currBody->getInteriorShells();
        currInteriorShells->reset4Loop();
        while ( currInteriorShells->setNext4Loop() ) {
             numOfVertices += currInteriorShells->getEntity()->getVertices()->getSize();
        }        
    }

    return numOfVertices;
}


rg_INT ManifoldizedBetaShape::getNumOfEdges() const
{
    rg_INT numOfEdges = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        numOfEdges += currBody->getExteriorShell()->getEdges()->getSize();

        rg_dList< MBSShell* >* currInteriorShells = currBody->getInteriorShells();
        currInteriorShells->reset4Loop();
        while ( currInteriorShells->setNext4Loop() ) {
             numOfEdges += currInteriorShells->getEntity()->getEdges()->getSize();
        }        
    }

    return numOfEdges;
}


rg_INT ManifoldizedBetaShape::getNumOfFaces() const
{
    rg_INT numOfFaces = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        numOfFaces += currBody->getExteriorShell()->getFaces()->getSize();

        rg_dList< MBSShell* >* currInteriorShells = currBody->getInteriorShells();
        currInteriorShells->reset4Loop();
        while ( currInteriorShells->setNext4Loop() ) {
             numOfFaces += currInteriorShells->getEntity()->getFaces()->getSize();
        }        
    }

    return numOfFaces;
}
    

//  Youngsong Cho - 2011-04-04
rg_INT ManifoldizedBetaShape::getNumOfHandlesByIgnoringInteriorVoids() const
{
    rg_INT numHandles = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        numHandles += currBody->getNumHandlesByIgnoringInteriorVoids();
    }

    return numHandles;
}




    
//set function
void ManifoldizedBetaShape::setBetaUniverse( BetaUniverse* betaUniverse )
{
    m_betaUniverse = betaUniverse;
}


void ManifoldizedBetaShape::setBetaValue( const rg_REAL& betaValue )
{
    m_betaValue = betaValue;
}

void ManifoldizedBetaShape::addBody( MBSBody* body )
{
    m_bodies.add( body );
}




//  Youngsong Cho - 2011-04-04
void ManifoldizedBetaShape::construct( BetaUniverse* betaUniverse, const rg_REAL& betaValue, const rg_INT& constructOption )
{
    AugmentedBetaComplex augmentedBetaComplex;
    augmentedBetaComplex.construct( betaValue, betaUniverse, constructOption );

	augmentedBetaComplex.extractBoundaryOfAugmentedBetaComplex( *this );
}





void ManifoldizedBetaShape::removeAll()
{
    m_betaUniverse = rg_NULL;
    m_betaValue = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        delete currBody;
        currBody = rg_NULL;
    }

    m_bodies.removeAll();
}



//rg_REAL --> rg_INT 현재는 debugging을 위해
rg_REAL ManifoldizedBetaShape::getNumOfHoles()
{
    rg_INT numOfHoles = 0;

    m_bodies.reset4Loop();
    while ( m_bodies.setNext4Loop() ) {
        MBSBody* currBody = m_bodies.getEntity();

        if ( currBody->isSingularBody() ) {
            continue;
        }

        MBSShell* currExteriorShell = currBody->getExteriorShell();
        rg_INT euler = currExteriorShell->getNumVertices() - currExteriorShell->getNumEdges() + currExteriorShell->getNumFaces();

        numOfHoles += 1.0 - euler/2.0;    //1.0, 2.0 --> 1, 2 현재는 debugging을 위해

        
        rg_dList< MBSShell* >* currInteriorShells = currBody->getInteriorShells();
        currInteriorShells->reset4Loop();
        while ( currInteriorShells->setNext4Loop() ) {
            MBSShell*  currInteriorShell = currInteriorShells->getEntity();

            euler = currInteriorShell->getNumVertices() - currInteriorShell->getNumEdges() + currInteriorShell->getNumFaces();

            numOfHoles += 1.0 - euler/2.0;    //1.0, 2.0 --> 1, 2 현재는 debugging을 위해
        }        
    }

    return numOfHoles;
}
