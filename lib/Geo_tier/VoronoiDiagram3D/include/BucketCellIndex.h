#ifndef _BUCKETCELLINDEX_H_
#define _BUCKETCELLINDEX_H_

#include "rg_Const.h"

namespace V {

namespace GeometryTier {


class BucketCellIndex
{
public:
    rg_INT m_i;
    rg_INT m_j;
    rg_INT m_k;

    BucketCellIndex();
    BucketCellIndex(const rg_INT& i, const rg_INT& j, const rg_INT& k);
    BucketCellIndex(const BucketCellIndex& cellIndex);
    ~BucketCellIndex();

    inline rg_INT getI() const { return m_i;};
    inline rg_INT getJ() const { return m_j;};
    inline rg_INT getK() const { return m_k;};

    void   setI(const rg_INT& i);
    void   setJ(const rg_INT& j);
    void   setK(const rg_INT& k);
    void   setCellIndex(const rg_INT& i, const rg_INT& j, const rg_INT& k);

    BucketCellIndex& operator =(const BucketCellIndex& cellIndex);
};

} // namespace GeometryTier

} // namespace V


#endif

