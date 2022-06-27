#ifndef BUCKET_CELL_FOR_DISK_H_
#define BUCKET_CELL_FOR_DISK_H_

#include "Generator2D.h"

namespace BULL2D{
namespace GeometryTier {


class BucketCellForDisk
{
private:
    Generator2D*                m_containingGenerator;

public:
    BucketCellForDisk(void);
    BucketCellForDisk( Generator2D* const containingGenerator );
    BucketCellForDisk( const BucketCellForDisk& bucketCell );
    ~BucketCellForDisk(void);

    inline Generator2D*     getContainingGenerator() const { return m_containingGenerator; }
    inline void             setContainingGenerator( Generator2D* generator ) { m_containingGenerator = generator; }
};

} // namespace GeometryTier
} // namespace BULL2D

#endif

