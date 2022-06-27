//MBSBody.h
#ifndef _MBSBODY_H
#define _MBSBODY_H

#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"

class MBSShell;

class MBSBody : public TopologicalEntity
{
private:
    MBSShell*            m_exteriorShell;
    rg_dList<MBSShell*>  m_interiorShells;

    rg_BOOL              m_visited;

public:
    MBSBody();
    MBSBody(const rg_INT& ID);
    MBSBody(const MBSBody& body);
    ~MBSBody();

    MBSShell*            getExteriorShell() const;    
    rg_dList<MBSShell*>* getInteriorShells();
    rg_INT               getNumShells() const;
    rg_INT               getNumInteriorShells() const;

    inline rg_BOOL       isVisited() const { return m_visited; }
    inline void          isVisited(const rg_BOOL& visited) { m_visited = visited; }

    void                 setExteriorShell(MBSShell* exteriorShell);
    void                 addInteriorShell(MBSShell* interiorShell);
    
    MBSBody&             operator =(const MBSBody& o_body);
    
    void connectExteriorShell(MBSShell* exteriorShell);
    void connectInteriorShell(MBSShell* interiorShell);
    void searchWholeShells(rg_dList<MBSShell*>& shellList);

    rg_BOOL isSingularBody();

    //  Youngsong Cho - 2011-04-04
    rg_INT getNumHandlesByIgnoringInteriorVoids() const;
};

#endif 
