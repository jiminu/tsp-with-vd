#ifndef _BUCKET_H_
#define _BUCKET_H_

#include "rg_Const.h"
#include "rg_dList.h"
#include "VDCell.h"
#include "Sphere.h"
#include "rg_BallGenerator.h"
#include "BucketCellIndex.h"


#include <fstream>
using namespace std;


namespace V {

namespace GeometryTier {


const rg_FLAG INVALID_CELL   = 0;
const rg_FLAG ON_PLANE       = 1;
const rg_FLAG IN_PLUS_PLANE  = 2;
const rg_FLAG IN_MINUS_PLANE = 3;


class Bucket
{
private:
	rg_dList< BallGenerator* >*** m_bucketTable;
    rg_FLAG***                    m_bucketMark;   //  by Youngsong Cho (2005. 7. 22)
    rg_INT   m_numOfXs;
	rg_INT   m_numOfYs;
	rg_INT   m_numOfZs;

	//bounding box;
	rg_Point3D m_minPointOfBoundingBox;	// POINT(min_X, min_Y, min_Z)
	rg_Point3D m_maxPointOfBoundingBox;	// POINT(max_X, max_Y, max_Z)

	//length of diagonal of a bucket
	rg_REAL m_diagonalLength;

    rg_REAL m_unitSize;

public:
	Bucket();
	Bucket(rg_dList< BallGenerator >& ballList);
	~Bucket();

	rg_Point3D getMinPointOfBoundingBox() const;
	rg_Point3D getMaxPointOfBoundingBox() const;
	rg_REAL    getDiagonalLength() const;

    void setMinPointOfBoundingBox(const rg_Point3D& minPt);
    void setMaxPointOfBoundingBox(const rg_Point3D& maxPt);

	void setBucket(rg_dList< BallGenerator >& ballList);
	void constructBucket( const rg_Point3D&         minPt, 
                          const rg_Point3D&         maxPt, 
                                rg_dList< BallGenerator >& ballList );

	void constructBoundingBox(rg_dList< BallGenerator >& ballList);

	rg_dList< BallGenerator* >* getCellsByMask(const Sphere& sphere);
	
	void addCellIntoBucket(BallGenerator* aBall);
	void getBucketIndex(const rg_Point3D& pt, rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ);

    
    rg_FLAG isThereIntersectingBallWithSphere(const Sphere& sphere);

//private:
    void removeBucketTable(); 


    //  by Youngsong Cho (2005. 7. 22)
	void constructBucket( const rg_REAL&            loadFactor,
                          const rg_Point3D&         minPt, 
                          const rg_Point3D&         maxPt, 
                          const rg_dList< BallGenerator >& ballList );
	void constructBucketBySizeOfElement( const rg_REAL&            sizeOfCube,
                                         const rg_Point3D&         minPt, 
                                         const rg_Point3D&         maxPt, 
                                         const rg_dList< BallGenerator >& ballList );
	    void addBallIntoBucketConsistingOfCubes(BallGenerator* aBall);

    inline rg_FLAG getMark(const rg_INT& i, const rg_INT& j, const rg_INT& k) const {return m_bucketMark[i][j][k];};
    inline rg_FLAG getMark(const BucketCellIndex& cellIndex) const {return m_bucketMark[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k];};
    inline void setMark(const rg_INT& i, const rg_INT& j, const rg_INT& k, const rg_FLAG& mark) 
                      { m_bucketMark[i][j][k] = mark;};
    inline void setMark(const BucketCellIndex& cellIndex, const rg_FLAG& mark)  
                      { m_bucketMark[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k] = mark;};

    rg_FLAG isVisited(const BucketCellIndex& cellIndex) const;

    rg_FLAG setMarkAfterCheck(const rg_INT& i, const rg_INT& j, const rg_INT& k, const rg_FLAG& mark);
    rg_FLAG setMarkAfterCheck(const BucketCellIndex& cellIndex, const rg_FLAG& mark);  

    rg_REAL     getXLengthOfCell() const;
    rg_REAL     getYLengthOfCell() const;
    rg_REAL     getZLengthOfCell() const;
    rg_Point3D  getLengthOfCell() const;

    rg_FLAG  getPositionOfCellWRTPlane(const BucketCellIndex& cell,
                                       const rg_REAL&         distance,
                                       const rg_Point3D&      cellLength);

    rg_INT                      getNumBalls(const rg_INT& indexX, const rg_INT& indexY, const rg_INT& indexZ);
    rg_dList< BallGenerator* >* getBalls(const BucketCellIndex& cellIndex);
    void                        getBalls( rg_dList< BallGenerator* >& balls, 
                                          const BucketCellIndex&      cellIndex );
    rg_dList< BallGenerator* >* getBalls( const BucketCellIndex& cellIndex, rg_REAL* plane);
    void                        getBalls( rg_dList< BallGenerator* >& balls, 
                                          const BucketCellIndex& cellIndex, rg_REAL* plane);
    rg_dList< BallGenerator* >* getBalls(rg_REAL* plane);

    rg_INT getCellsByMask(rg_dList< BallGenerator* >& resultList, const Sphere& sphere);
    void getBallsIntersectingSphere(const Sphere& sphere, rg_dList< BallGenerator* >& resultList);
    void getBallsBoundedByPlaneAndSphere(rg_dList< BallGenerator* >& resultList, 
                                         rg_REAL*                    plane, 
                                         const Sphere&               sphere);

	void getJustifiedBucketIndex( const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, 
                                  rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ);
	void getJustifiedBucketIndex( const rg_Point3D& pt,  
                                  rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ);

    void reportBucket(ofstream& fout);
    void getBalls(ofstream& fout, rg_REAL* plane);

    rg_FLAG isConstructed() const;
    rg_INT  getNumElementInX() const;
    rg_INT  getNumElementInY() const;
    rg_INT  getNumElementInZ() const;

    rg_INT  getNumEmptyElements();
};

} // namespace GeometryTier

} // namespace V


#endif


