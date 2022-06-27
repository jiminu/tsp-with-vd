#include "MBSFace.h"
#include "MBSEdge.h"
#include "MBSVertex.h"
#include "BetaVertex.h"

MBSFace::MBSFace()
: TopologicalEntity(),
  m_firstEdge(rg_NULL),
  m_betaFace(rg_NULL),
  m_isArtificial(rg_FALSE),
  m_visited(rg_FALSE)
{
}



MBSFace::MBSFace(const rg_INT& ID)
: TopologicalEntity(ID),
  m_firstEdge(rg_NULL),
  m_betaFace(rg_NULL),
  m_isArtificial(rg_FALSE),
  m_visited(rg_FALSE)
{
}







MBSFace::MBSFace(const MBSFace& face)
: TopologicalEntity(face),
  m_firstEdge(face.m_firstEdge),
  m_betaFace(face.m_betaFace),
  m_isArtificial(face.m_isArtificial),
  m_visited(face.m_visited)
{    
}



MBSFace::~MBSFace()
{
}





MBSEdge* MBSFace::getFirstEdge() const
{
    return m_firstEdge;
}

BetaFace* MBSFace::getBetaFace() const
{
    return m_betaFace;
}


rg_BOOL MBSFace::isArtificial() const
{
    return m_isArtificial;
}


MBSShell* MBSFace::getShell() const
{
    return m_firstEdge->getStartVertex()->getShell();
}



void MBSFace::setFirstEdge( MBSEdge* firstEdge )
{
    m_firstEdge = firstEdge;
}


void MBSFace::setBetaFace( BetaFace* betaFace )
{
    m_betaFace = betaFace;
}


void MBSFace::isArtificial( const rg_BOOL& isArtificial )
{
    m_isArtificial = isArtificial;
}





MBSFace&             MBSFace::operator =(const MBSFace& face)
{
    if ( this == &face ) {
        return *this;
    }

    TopologicalEntity::operator =(face);
    m_firstEdge     = face.m_firstEdge;
    m_betaFace      = face.m_betaFace;
    m_isArtificial  = face.m_isArtificial;
    m_visited       = face.m_visited;

    return *this;
}


    



void MBSFace::searchBoundingVertices(rg_dList<MBSVertex*>& vertexList)  
{
    MBSEdge* startEdge = m_firstEdge;
    MBSEdge* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            vertexList.add( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            vertexList.add( currEdge->getStartVertex() );
            currEdge = currEdge->getRightLeg();
        }
//         else  {
//             return rg_FALSE;
//         }
    } while ( startEdge != currEdge );   

//     return rg_TRUE;
}



void MBSFace::searchBoundingEdges(rg_dList<MBSEdge*>& edgeList)  
{
    MBSEdge* startEdge = m_firstEdge;
    MBSEdge* currEdge  = startEdge;

    do  {
        edgeList.add( currEdge );

        if ( this == currEdge->getLeftFace() ) {
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            currEdge = currEdge->getRightLeg();
        }
//         else {
//             return rg_FALSE;
//         }

    } while ( startEdge != currEdge );   

//     return rg_TRUE;
}



rg_FLAG MBSFace::searchIncidentFaces(rg_dList<MBSFace*>& faceList)  
{
    /*
    if ( MBSFace::getBoundingState() == SINGULAR_SIMPLEX ) {
        return rg_FALSE;
    }
    */
    if (    m_firstEdge->getLeftFace() == rg_NULL 
        || m_firstEdge->getRightFace() == rg_NULL ) {  //위의 조건과 같다. boundingState를 보는 것 보다는 topology를 check하는 것이 더 낫다.
        return rg_FALSE;
    }

    MBSEdge* startEdge = m_firstEdge;
    MBSEdge* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            faceList.add( currEdge->getRightFace() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            faceList.add( currEdge->getLeftFace() );
            currEdge = currEdge->getRightLeg();
        }
        else  {
            return rg_FALSE;
        }
    } while ( startEdge != currEdge );   

    return rg_TRUE;
}

rg_FLAG MBSFace::isZeroVolumeFace()
{
    rg_dList< MBSVertex* > boundingVertices;
    this->searchBoundingVertices( boundingVertices );

    boundingVertices.reset4Loop();
    Ball* boundingSpheres[3] = { rg_NULL, rg_NULL, rg_NULL };
    rg_INT i;
    for ( i=0; i<3; i++ )  {
        boundingVertices.setNext4Loop();

        boundingSpheres[i] = boundingVertices.getEntity()->getOriginalBetaVertex()->getBallProperty();
    }

    if ( boundingSpheres[0] == boundingSpheres[1] 
      || boundingSpheres[0] == boundingSpheres[2] 
      || boundingSpheres[1] == boundingSpheres[2] 
      ) {                                             //area = 0
      return rg_TRUE;        
    }
    else {
        return rg_FALSE;
    }

}



