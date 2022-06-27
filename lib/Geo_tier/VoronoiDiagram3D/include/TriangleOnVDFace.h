#ifndef _TRIANGLEONVDFACE_H
#define _TRIANGLEONVDFACE_H

#include "Triangle.h"

namespace V {

namespace GeometryTier {


class TriangleOnVDFace : public Triangle
{
private:
    rg_REAL m_distanceToGenerator;

    rg_REAL m_distance[3];

public:
    TriangleOnVDFace();
    TriangleOnVDFace(const TriangleOnVDFace& aTriangle);
    ~TriangleOnVDFace();

    rg_REAL getDistanceToGenerator() const;
    rg_REAL getDistance(const rg_INT& i) const;

    void    setDistanceToGenerator(const rg_REAL& distance);
    void    setDistance(const rg_INT& i, const rg_REAL& distance);

    TriangleOnVDFace& operator =(const TriangleOnVDFace& aTriangle);
};

} // namespace GeometryTier

} // namespace V


#endif
