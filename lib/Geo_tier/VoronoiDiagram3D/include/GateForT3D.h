#ifndef _GATEFORT3D_H
#define _GATEFORT3D_H

#include "T3DVertex.h"
#include "T3DTetrahedron.h"

namespace V {

namespace GeometryTier {


class GateForT3D
{
private:
    T3DVertex* m_vertex;
    T3DVertex* m_mateVertex;
    T3DTetrahedron* m_tetrahedron;

public:
    GateForT3D();
    GateForT3D(T3DVertex* vertex, T3DVertex* mateVertex, T3DTetrahedron* tetrahedron);
    GateForT3D(const GateForT3D& gate);
    ~GateForT3D();

    T3DVertex* getVertex() const;
    T3DVertex* getMateVertex() const;
    T3DTetrahedron* getTetrahedron() const;

    void setVertex(T3DVertex* vertex);
    void setMateVertex(T3DVertex* mateVertex);
    void setTetrahedron(T3DTetrahedron* tetrahedron);

    GateForT3D& operator =(const GateForT3D& gate);

};

} // namespace GeometryTier

} // namespace V


#endif


