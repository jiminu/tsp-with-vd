#ifndef _MBSSHELL_H
#define _MBSSHELL_H

#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"

class MBSBody;
class MBSFace;
class MBSEdge;
class MBSVertex;

class MBSShell : public TopologicalEntity
{
private:
    MBSBody*             m_body;

    rg_dList<MBSFace*>   m_faces;
    rg_dList<MBSEdge*>   m_edges;
    rg_dList<MBSVertex*> m_vertices;

    rg_BOOL              m_visited;

public:
    MBSShell();
    MBSShell(const rg_INT& ID);
    MBSShell(const rg_INT& ID, MBSBody* body);
    MBSShell(const MBSShell& o_shell);
    ~MBSShell();

    MBSBody*             getBody() const;
    
    rg_INT              getNumFaces() const;
    rg_INT              getNumEdges() const;
    rg_INT              getNumVertices() const;

    rg_dList<MBSFace*>*   getFaces();
    rg_dList<MBSEdge*>*   getEdges();
    rg_dList<MBSVertex*>* getVertices();

    inline rg_BOOL      isVisited() const { return m_visited; }
    inline void         isVisited(const rg_BOOL& visited) { m_visited = visited; }

    void                setBody(MBSBody* body);
    void                addFace(MBSFace* face);
    void                addEdge(MBSEdge* edge);
    void                addVertex(MBSVertex* vertex);

    MBSShell&           operator =(const MBSShell& shell);
};

#endif
