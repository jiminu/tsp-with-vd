#ifndef BUCKET_FOR_DISKS_H_
#define BUCKET_FOR_DISKS_H_

#include "Generator2D.h"
#include "BucketCellForDisk.h"
#include "rg_BoundingBox2D.h"

namespace V {
namespace GeometryTier {


class BucketForDisks
{
private:
    BucketCellForDisk**     m_bucketCell;

    int                     m_numBucketCellOfXDir;
    int                     m_numBucketCellOfYDir;

    double                  m_xDirLengthOfBucketCell;
    double                  m_yDirLengthOfBucketCell;

    rg_BoundingBox2D        m_boundingBox;

public:
    BucketForDisks(void);
	BucketForDisks(const BucketForDisks& BFD); // cysong
    BucketForDisks( const list<Generator2D*>& diskSet );
    ~BucketForDisks(void);
	
    inline int          getNumBucketCellOfXDir() const { return m_numBucketCellOfXDir; }
    inline int          getNumBucketCellOfYDir() const { return m_numBucketCellOfYDir; }
    inline double       getxDirLengthOfBucketCell() const { return m_xDirLengthOfBucketCell; }
    inline double       getyDirLengthOfBucketCell() const { return m_yDirLengthOfBucketCell; }
    inline rg_Point2D   getMinPt() const { return m_boundingBox.getMinPt(); }
    inline rg_Point2D   getMaxPt() const { return m_boundingBox.getMaxPt(); }
	inline const rg_BoundingBox2D& getBoundingBox() const { return m_boundingBox; };

    void                getBucketCellIndexForPoint( const rg_Point2D& pt, int& xIndex, int& yIndex ) const;


    void                setBucketCell( const list<Generator2D*>& diskSet );
    void                setBoundingBox( const list<Generator2D*>& diskSet );


    void                addGenerator( Generator2D* generator );
    Generator2D*        findGoodCloseGenerator( const rg_Point2D& pt ) const;
    void                destroyBucket();  //cysong : bucket을 destroy할시에  cell 갯수 초기화
	void                copyBucketFrom(const BucketForDisks& BFD);
	
	//operators
	BucketForDisks&     operator=(const BucketForDisks& BFD);
};

} // GeometryTier
} // V

#endif
