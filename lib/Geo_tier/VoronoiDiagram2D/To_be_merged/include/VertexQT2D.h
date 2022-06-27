#ifndef _VERTEXQT2D_H
#define _VERTEXQT2D_H

#include "VertexSCDS.h"
#include "Disc.h"

namespace BULL2D {
namespace GeometryTier {


class VertexQT2D : public VertexSCDS
{
private:
    Disc* m_disc;

public:
    VertexQT2D();
    VertexQT2D(const rg_INT& ID, Disc* disc);
    VertexQT2D(const VertexQT2D& vertex);
    ~VertexQT2D();
    
    Disc*       getDisc();
    rg_Circle2D getDiscGeometry() const;

    rg_BOOL     isVirtual() const;

    void        setDisc(Disc* disc);

    VertexQT2D& operator =(const VertexQT2D& vertex);
};

} // namespace GeometryTier
} // namespace BULL2D


#endif

