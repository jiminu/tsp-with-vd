#ifndef BUCKET_CELL_FOR_DISK_H_
#define BUCKET_CELL_FOR_DISK_H_

#include "Generator2D.h"

namespace V {
namespace GeometryTier {

class BucketCellForDisk
{
private:
    Generator2D*                m_containingGenerator;
    //bool                        m_isChecked;

public:
    BucketCellForDisk(void);
    BucketCellForDisk( Generator2D* const containingGenerator );
    BucketCellForDisk( const BucketCellForDisk& bucketCell );
    ~BucketCellForDisk(void);

    inline Generator2D*     getContainingGenerator() const { return m_containingGenerator; }
    inline void             setContainingGenerator( Generator2D* generator ) { m_containingGenerator = generator; }
};

} // GeometryTier
} // V

#endif

