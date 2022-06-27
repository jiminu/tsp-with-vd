#ifndef _T3DTETRAHEDRON_H
#define _T3DTETRAHEDRON_H

#include "TopologicalEntity.h"
#include "T3DFace.h"

namespace V {

namespace GeometryTier {


class T3DVertex;

class T3DTetrahedron : public TopologicalEntity
{
private:
	T3DVertex*      m_vertex[4];
	T3DTetrahedron* m_neighbor[4];

    rg_FLAG         m_check;

public:
	T3DTetrahedron();
	T3DTetrahedron(const rg_INT& ID);
    T3DTetrahedron(const rg_INT& ID,
                   T3DVertex* vertex1, T3DVertex* vertex2, 
                   T3DVertex* vertex3, T3DVertex* vertex4);
	T3DTetrahedron(const T3DTetrahedron& tetrahedron);
	~T3DTetrahedron();

	T3DVertex*       getVertex(const rg_INDEX& i) const;
	T3DVertex**      getVertices();
    T3DTetrahedron*  getNeighbor(const rg_INDEX& i) const;
    T3DTetrahedron** getNeighbors();
    T3DTetrahedron*  getMateTetrahedron(const rg_INT& pos) const;
    T3DTetrahedron*  getMateTetrahedron(T3DVertex* vertex) const;
    T3DVertex*       getMateVertex(T3DTetrahedron* neighbor) const;
    rg_INT           getMateVertexPos(T3DTetrahedron* neighbor) const;

    rg_INT           getPosOfVertex(T3DVertex* vertex) const;
    rg_INT           getMateVertexPosOfIncidentFace(T3DTetrahedron* neighbor, const rg_INT& mateVertexPos);
    void             getVerticesOnFace(const rg_INT& mateVertexPos, T3DVertex** vertices);
    T3DFace          getFace(const rg_INT& mateVertexPos);
    
    rg_INT  locateNeighborPos(T3DTetrahedron* neighbor);
    rg_FLAG locateFace(T3DFace& faceToLocate, const rg_INT& vertexID1, const rg_INT& vertexID2, const rg_INT& vertexID3);

    rg_FLAG isThereThisVertex(T3DVertex* vertex) const;
    rg_FLAG isThereThisEdge(T3DVertex* vertex1, T3DVertex* vertex2) const;
    rg_FLAG isThereThisNeighbor(T3DTetrahedron* neighbor) const;
    rg_FLAG isOnBoundaryOfConvexHull() const;

    rg_FLAG isChecked() const;
    void    setCheck(const rg_FLAG& check);

    void setVertex(const rg_INDEX& i, T3DVertex* vertex);
    void setVertices(T3DVertex* vertex1, T3DVertex* vertex2, 
                     T3DVertex* vertex3, T3DVertex* vertex4);
    void setVertices(T3DVertex** vertices);
    void setNeighbor(const rg_INDEX& i, T3DTetrahedron* neighbor);
    void setNeighbor(T3DVertex* vertex, T3DTetrahedron* neighbor);
    void setNeighbors(T3DTetrahedron* neighbor1, T3DTetrahedron* neighbor2, 
                      T3DTetrahedron* neighbor3, T3DTetrahedron* neighbor4);
    void setNeighbors(T3DTetrahedron** neighbors);

    T3DTetrahedron& operator =(const T3DTetrahedron& tetrahedron);

};

} // namespace GeometryTier

} // namespace V


#endif

