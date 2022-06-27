#ifndef T3DFACE_H
#define T3DFACE_H

#include "rg_Const.h"

namespace V {

namespace GeometryTier {


class T3DTetrahedron;
class T3DVertex;

class T3DFace
{
private:
    T3DTetrahedron* m_tetrahedron;
    rg_INT          m_mateVertexPos;
    
public:
    T3DFace();
    T3DFace(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos);
    T3DFace(const T3DFace& face);
    ~T3DFace();

    T3DTetrahedron* getTetrahedron() const;
    rg_INT          getMateVertexPos() const;
    T3DTetrahedron* getIncidentNeighbor() const;
    void            getBoundingVertices(T3DVertex** vertices) const;

    void setTetrahedron(T3DTetrahedron* tetrahedron);
    void setMateVertexPos(const rg_INT& mateVertexPos);
    void setFace(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos);

    rg_FLAG compareFaceByVertices( const T3DFace& face );

    T3DFace& operator =(const T3DFace& face);
};

} // namespace GeometryTier

} // namespace V


#endif

