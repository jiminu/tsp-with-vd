#ifndef _FACEQT2D_H
#define _FACEQT2D_H

#include "FaceSCDS.h"
#include "rg_Circle2D.h"


namespace V {
namespace GeometryTier {


class FaceQT2D : public FaceSCDS
{
private:
    rg_Circle2D m_emptyTangentCircle;

public:
    FaceQT2D();
    FaceQT2D(const rg_INT& ID, const rg_Circle2D& tangentCircle);
    FaceQT2D(const FaceQT2D& face);
    ~FaceQT2D();

    rg_Circle2D getEmptyTangentCircle() const;

    rg_BOOL     isVirtual() const;

    void        setEmptyTangentCircle(const rg_Circle2D& tangentCircle);

    FaceQT2D&   operator =(const FaceQT2D& face);

    ///////////////////////////////////////////////////////////////////////////

    rg_BOOL isAnomaly() const;
};


} // GeometryTier
} // V


#endif

