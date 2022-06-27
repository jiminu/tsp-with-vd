#ifndef _T3DGATE_H
#define _T3DGATE_H

#include "rg_Const.h"
#include "T3DVertex.h"
#include "T3DTetrahedron.h"

namespace V {

namespace GeometryTier {


class T3DGate
{
private:
    T3DTetrahedron* m_tetrahedron;
    rg_INT          m_mateVertexPos;
    T3DVertex*      m_vertex[3];
    T3DVertex*      m_vertexForNextTetrahedron[3];
    rg_FLAG         m_isValidVertex[3];
    rg_INT          m_indexOfPriorEdge;

public:
    rg_INT          m_status;

    
    T3DGate();
    T3DGate(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos);
    T3DGate(const T3DGate& gate);
    ~T3DGate();

    T3DTetrahedron* getTetrahedron() const;
    rg_INT          getMateVertexPos() const;
    T3DVertex*      getVertex(const rg_INT& i) const;
    T3DVertex**     getVertices();
    T3DVertex*      getVertexForNextTetrahedron(const rg_INT& i) const;
    T3DVertex**     getVerticesForNextTetrahedron();
    T3DTetrahedron* getIncidentNeighbor() const;
    rg_INT          getIndexOfPriorEdge() const;
    void            getVerticesOfPriorEdge(T3DVertex*& startVertex, T3DVertex*& endVertex);

    rg_FLAG isValidVertexForNextTetrahedron(const rg_INT& i) const;
    rg_FLAG compareGateByVertices( const T3DGate& gate );
    rg_FLAG compareGateByVertices( T3DVertex** vertexOnGate );

    void setTetrahedron(T3DTetrahedron* tetrahedron);
    void setMateVertexPos(const rg_INT& pos);
    void setVertexForNextTetrahedron(const rg_INT& i, T3DVertex* vertex);
    void setGate(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos);
    void setIndexOfPriorEdge(const rg_INT& i);

    T3DGate& operator =(const T3DGate& gate);
};

} // namespace GeometryTier

} // namespace V


#endif

