#include "MBSBody.h"
#include "MBSShell.h"


MBSBody::MBSBody()
: TopologicalEntity(),
  m_exteriorShell(rg_NULL),
  m_visited(rg_FALSE)
{
}



MBSBody::MBSBody(const rg_INT& ID)
: TopologicalEntity(ID),
  m_exteriorShell(rg_NULL),
  m_visited(rg_FALSE)
{
}



MBSBody::MBSBody(const MBSBody& body)
: TopologicalEntity(body),
  m_exteriorShell(body.m_exteriorShell),
  m_visited(body.m_visited)
{
    m_interiorShells.duplicateList( body.m_interiorShells );
}



MBSBody::~MBSBody()
{
    if ( m_exteriorShell != rg_NULL ) {
        delete m_exteriorShell;
        m_exteriorShell = rg_NULL;
    }


    rg_INT numInteriorShells = m_interiorShells.getSize();
    rg_dNode<MBSShell*>* shellNode = m_interiorShells.getFirstpNode();
    for ( rg_INT i=0; i<numInteriorShells; i++, shellNode=shellNode->getNext() ) {
        MBSShell* interiorShell = shellNode->getEntity();
        if ( interiorShell != rg_NULL ) {
            delete interiorShell;
            
            shellNode->setEntity( rg_NULL );
        }
    }
    m_interiorShells.removeAll();
}




MBSShell* MBSBody::getExteriorShell() const
{
    return m_exteriorShell;
}



rg_INT MBSBody::getNumShells() const
{
    rg_INT numShells = 0;

    if ( m_exteriorShell != rg_NULL ) {
        numShells++;
    }

    numShells += m_interiorShells.getSize();

    return numShells;
}



rg_INT MBSBody::getNumInteriorShells() const
{
    return m_interiorShells.getSize();
}



rg_dList<MBSShell*>* MBSBody::getInteriorShells() 
{
    return &m_interiorShells;
}




void MBSBody::setExteriorShell(MBSShell* exteriorShell)
{
    m_exteriorShell = exteriorShell;
}



void MBSBody::addInteriorShell(MBSShell* interiorShell)
{
    m_interiorShells.addWithoutSame( interiorShell );
}



void MBSBody::connectExteriorShell(MBSShell* exteriorShell)
{
    m_exteriorShell = exteriorShell;
    exteriorShell->setBody( this );
}



void MBSBody::connectInteriorShell(MBSShell* interiorShell)
{
    m_interiorShells.addWithoutSame( interiorShell );

    interiorShell->setBody( this );
}



MBSBody& MBSBody::operator =(const MBSBody& mbs_body)
{
    if ( this == &mbs_body ) {
        return *this;
    }

    TopologicalEntity::operator =(mbs_body);

    m_exteriorShell = mbs_body.m_exteriorShell;
    m_interiorShells.duplicateList( mbs_body.m_interiorShells );

    m_visited = mbs_body.m_visited;

    return *this;
}



void MBSBody::searchWholeShells(rg_dList<MBSShell*>& shellList)
{
    shellList.add( m_exteriorShell );
    shellList.appendTail( m_interiorShells );
}



rg_BOOL MBSBody::isSingularBody()
{
    if ( m_interiorShells.getSize() != 0 ) {
        return rg_FALSE;
    }

    if (    m_exteriorShell->getNumFaces() == 1
         || m_exteriorShell->getNumEdges() == 1
         || m_exteriorShell->getNumVertices() == 1 ) {
        return rg_TRUE;
    }
    

    return rg_FALSE;
}



//  Youngsong Cho - 2011-04-04
rg_INT MBSBody::getNumHandlesByIgnoringInteriorVoids() const
{
    rg_INT V = m_exteriorShell->getNumVertices();
    rg_INT E = m_exteriorShell->getNumEdges();
    rg_INT F = m_exteriorShell->getNumFaces();
    rg_INT eulerChar = V - E + F;
    rg_INT S = 1;
    rg_INT B = 0;

    rg_INT numHandles = S - (eulerChar + B)/2;

    return numHandles;
}
