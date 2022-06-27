#ifndef BUCKET_FOR_DISK_H_
#define BUCKET_FOR_DISK_H_

#include "BucketCellForDisks.h"
#include "rg_BoundingBox2D.h"
#include "rg_Circle2D.h"
#include "rg_Line2D.h"
#include "rg_GeoFunc.h"
#include<vector>
#include <list>
using namespace std;

class BucketForDisk
{
public:
    enum CellDivisionOptions {
        BY_LARGEST_RADIUS,
        BY_SMALLEST_RAIDUS,
        BY_FIX_NUM_100_N_100,
        NUM_OF_DISKS_SAME_AS_NUM_OF_CELLS
    };

    static constexpr double RES_OF_POINT_IS_CONTAINED_IN_THIS_BUCKET = 1.0e-5;
    static constexpr double RES_OF_LINE_SEGMENT                      = 1.0e-5;


private:
	//BucketCellForDisks**    m_bucketCell;

    vector < vector<BucketCellForDisks*> > m_bucketCell;
    int                     m_numBucketCellOfXDir;
    int                     m_numBucketCellOfYDir;

    double                  m_xDirLengthOfBucketCell;
    double                  m_yDirLengthOfBucketCell;

    rg_BoundingBox2D        m_boundingBox;
    list<rg_Circle2D*>      m_DisksInBucket;

public:
    BucketForDisk();
	BucketForDisk(const BucketForDisk& BFD); 
    BucketForDisk(const list<rg_Circle2D>& disks, const CellDivisionOptions& divisionOptions = BY_LARGEST_RADIUS, const bool& DisksAreArrangedInBucket = true);
    BucketForDisk(const list<rg_Circle2D*>& disks, const CellDivisionOptions& divisionOptions = BY_LARGEST_RADIUS, const bool& DisksAreArrangedInBucket = true);
    ~BucketForDisk(void);
	
	//inline BucketCellForDisks**		getBucketCell() const { return m_bucketCell; }
    inline const vector < vector<BucketCellForDisks*> >&		getBucketCell() const { return m_bucketCell; }

    inline int						getNumBucketCellOfXDir() const { return m_numBucketCellOfXDir; }
    inline int						getNumBucketCellOfYDir() const { return m_numBucketCellOfYDir; }
    inline double					getxDirLengthOfBucketCell() const { return m_xDirLengthOfBucketCell; }
    inline double					getyDirLengthOfBucketCell() const { return m_yDirLengthOfBucketCell; }
    inline rg_Point2D				getMinPt() const { return m_boundingBox.getMinPt(); }
    inline rg_Point2D				getMaxPt() const { return m_boundingBox.getMaxPt(); }
	inline const rg_BoundingBox2D&	getBoundingBox() const { return m_boundingBox; };
    inline void                     getDisksInBucket(list<rg_Circle2D*>& disksInBucket) const { disksInBucket = m_DisksInBucket; };

    void generateBucketForDisks(const list<rg_Circle2D*>& disks, const CellDivisionOptions& divisionOptions, const bool& DisksAreArrangedInBucket);

    inline void getBucketCellIndexForPoint( const rg_Point2D& pt, int& xIndex, int& yIndex ) const;
    BucketCellForDisks* getBucketCellIncludingPoint(const rg_Point2D& pt) const;

	BucketCellForDisks* getBucketCellContainingThisDisk(const rg_Circle2D* circle) const;
	void getMax9BucketCellsIncludingNAroundThisCell(BucketCellForDisks* cell, list<BucketCellForDisks*>& max9Cells) const;
	inline void getMaxNMinPtsOfThisCell(const int& xIndexOfCell, const int& yIndexOfCell, rg_Point2D& maxPt, rg_Point2D& minPt) const;
	inline BucketCellForDisks* getBucketCellAccordingIndex(const int& xIndexOfCell, const int& yIndexOfCell) const;
    
    double findMaxRadius(const list<rg_Circle2D*>& disks) const;
    double findMinRadius(const list<rg_Circle2D*>& disks) const;
    void setBucketCellSize(const double& cellSize); 
    void setBucketCellNumbers(const int& numX, const int& numY);  

    void setBoundingBox(const list<rg_Circle2D*>& disks); 
    void setBoundingBox(const rg_Point2D& maxPt, const rg_Point2D& minPt);
    void addDisksIntoBucket(const list<rg_Circle2D*>& disks);
    void addDiskIntoBucket(rg_Circle2D* disk);
    void destroyBucket();  //cysong : bucket을 destroy할시에  cell 갯수 초기화

    void  clear();
    void  copyBucketFrom(const BucketForDisk& BFD);
	
	//operators
	BucketForDisk&     operator=(const BucketForDisk& BFD);

    //CYSONG[Oct17, 19]: function added
    void findAllCellsIntersectedWithThisLineSegment(const rg_Line2D& lineSegment, list<BucketCellForDisks*>& intersectedBucketCells) const;
    void getMax9BucketCellsIncludingNAroundThisCells(const list<BucketCellForDisks*>& inputCells, list<BucketCellForDisks*>& cellsAroundNIncludingInputCells) const;
    void extractAllDisksIncludedInTheseCells(const list<BucketCellForDisks*>& bucketCells, list<rg_Circle2D*>& disks) const;
    inline pair<bool, rg_Point2D> thisLineSementIntersectsThisBox(const rg_Line2D& lineSegment, const rg_Point2D& minPtOfCell, const rg_Point2D& maxPtOfCell) const;
    inline bool this_value_is_in_this_interval(const double& value, const double& lowerVal, const double& upperVal, const double& res) const;
    inline void compute_four_boundary_line_segments(const rg_Point2D& minPt, const rg_Point2D& maxPt, vector<rg_Line2D>& fourBoundaries_buttom_right_top_left_order) const;
    inline bool two_line_segments_are_intersected(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2, const double& res) const;
    inline bool this_point_is_on_line_segment(const rg_Line2D& lineSegment, const rg_Point2D& endPointOfTheOtherLineSegment, const double& res) const;

    //CYSONG[Mar31, 20] 
    bool there_is_a_pair_of_disks_intersected(const double& tolerance = rg_MATH_RES) const;
};



inline void BucketForDisk::addDiskIntoBucket(rg_Circle2D* disk)
{
    int Xindex = 0;
    int Yindex = 0;

    getBucketCellIndexForPoint(disk->getCenterPt(), Xindex, Yindex);
    m_bucketCell[Xindex][Yindex]->insert_disk(disk);
    m_DisksInBucket.push_back(disk);
}

inline BucketCellForDisks* BucketForDisk::getBucketCellAccordingIndex(const int& xIndexOfCell, const int& yIndexOfCell) const
{
    BucketCellForDisks* targetBucketCell = const_cast<BucketCellForDisks*>(m_bucketCell[xIndexOfCell][yIndexOfCell]);

    return targetBucketCell;
}


inline pair<bool, rg_Point2D> BucketForDisk::thisLineSementIntersectsThisBox(const rg_Line2D& lineSegment, const rg_Point2D& minPtOfCell, const rg_Point2D& maxPtOfCell) const
{
     rg_Point2D boundaryPts[4] = { rg_Point2D(minPtOfCell),
                                   rg_Point2D(maxPtOfCell.getX(), minPtOfCell.getY()),
                                   rg_Point2D(maxPtOfCell),
                                   rg_Point2D(minPtOfCell.getX(), maxPtOfCell.getY()) };

    for (int i = 0; i < 4; ++i)
    {
        double parameterOfIntersectionForLineSeg1  = DBL_MAX;
        double parameterOfIntersectionForLineSeg2  = DBL_MAX;
        bool  bTwoLinesAreParallel = false;

        rg_Point2D startPointOfIntersection = rg_GeoFunc::compute_intersection_between_two_lines(
            boundaryPts[0], boundaryPts[1], lineSegment.getSP(), lineSegment.getEP(),
            parameterOfIntersectionForLineSeg1, parameterOfIntersectionForLineSeg2,
            bTwoLinesAreParallel);

        if (this_value_is_in_this_interval(parameterOfIntersectionForLineSeg1, 0.0, 1.0, resNeg5)
         && this_value_is_in_this_interval(parameterOfIntersectionForLineSeg2, 0.0, 1.0, resNeg5))
        {
            return make_pair(true, startPointOfIntersection);
        }
    }
  
    return make_pair(false, rg_Point2D(DBL_MAX, DBL_MAX));
}


inline bool BucketForDisk::this_value_is_in_this_interval(
    const double& value, const double& lowerVal, const double& upperVal, const double& res) const
{
    if (rg_GE(value, lowerVal, res) && rg_LE(value, upperVal, res))
    {
        return true;
    }
    else
    {
        return false;
    }
}



inline void BucketForDisk::getMaxNMinPtsOfThisCell(const int& xIndexOfCell, const int& yIndexOfCell, rg_Point2D& maxPt, rg_Point2D& minPt) const
{
    rg_Point2D bucketMinPt = m_boundingBox.getMinPt();

    double minPtX = bucketMinPt.getX() + (xIndexOfCell)*m_xDirLengthOfBucketCell;
    double minPtY = bucketMinPt.getY() + (yIndexOfCell)*m_yDirLengthOfBucketCell;

    double maxPtX = bucketMinPt.getX() + (xIndexOfCell + 1)*m_xDirLengthOfBucketCell;
    double maxPtY = bucketMinPt.getY() + (yIndexOfCell + 1)*m_yDirLengthOfBucketCell;

    minPt.setPoint(minPtX, minPtY);
    maxPt.setPoint(maxPtX, maxPtY);
}


inline void BucketForDisk::getBucketCellIndexForPoint(const rg_Point2D& pt, int& xIndex, int& yIndex) const
{
    xIndex = (rg_INT)((pt.getX() - m_boundingBox.getMinPt().getX()) / m_xDirLengthOfBucketCell);
    yIndex = (rg_INT)((pt.getY() - m_boundingBox.getMinPt().getY()) / m_yDirLengthOfBucketCell);

    if (xIndex < 0)
    {
        xIndex = 0;
    }
    if (xIndex > m_numBucketCellOfXDir - 1)
    {
        xIndex = m_numBucketCellOfXDir - 1;
    }

    if (yIndex < 0)
    {
        yIndex = 0;
    }
    if (yIndex > m_numBucketCellOfYDir - 1)
    {
        yIndex = m_numBucketCellOfYDir - 1;
    }
}


inline void BucketForDisk::compute_four_boundary_line_segments(const rg_Point2D& minPt, const rg_Point2D& maxPt, vector<rg_Line2D>& fourBoundaries_buttom_right_top_left_order) const
{
    rg_Point2D boundaryPts[4] = { rg_Point2D(minPt),
                                  rg_Point2D(maxPt.getX(), minPt.getY()),
                                  rg_Point2D(maxPt),
                                  rg_Point2D(minPt.getX(), maxPt.getY()) };

    fourBoundaries_buttom_right_top_left_order.resize(4);
    fourBoundaries_buttom_right_top_left_order[0] = rg_Line2D(boundaryPts[0], boundaryPts[1]);
    fourBoundaries_buttom_right_top_left_order[1] = rg_Line2D(boundaryPts[1], boundaryPts[2]);
    fourBoundaries_buttom_right_top_left_order[2] = rg_Line2D(boundaryPts[2], boundaryPts[3]);
    fourBoundaries_buttom_right_top_left_order[3] = rg_Line2D(boundaryPts[3], boundaryPts[0]); 
}


inline bool BucketForDisk::two_line_segments_are_intersected(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2, const double& res) const
{
    double orientationOfStartPtOfLine1 = lineSegment2.signed_distance(lineSegment1.getSP());
    double orientationOfEndPtOfLine1   = lineSegment2.signed_distance(lineSegment1.getEP());
    double orientationOfStartPtOfLine2 = lineSegment1.signed_distance(lineSegment2.getSP());
    double orientationOfEndPtOfLine2   = lineSegment1.signed_distance(lineSegment2.getEP());

    // 1. general case
    if (   abs(orientationOfStartPtOfLine1) > res
        && abs(orientationOfEndPtOfLine1)   > res
        && abs(orientationOfStartPtOfLine2) > res
        && abs(orientationOfEndPtOfLine2)   > res)
    {
        if ( orientationOfStartPtOfLine1 * orientationOfEndPtOfLine1 < 0.0
          && orientationOfStartPtOfLine2 * orientationOfEndPtOfLine2 < 0.0)
        {
            return true;
        }
    }
   
    //2. special case 
    if (   rg_ZERO(orientationOfStartPtOfLine1, res) 
        && this_point_is_on_line_segment(lineSegment2, lineSegment1.getSP(), res))
    {
        return true;
    }

    if (   rg_ZERO(orientationOfEndPtOfLine1, res)
        && this_point_is_on_line_segment(lineSegment2, lineSegment1.getEP(), res))
    {
        return true;
    }

    if (   rg_ZERO(orientationOfStartPtOfLine2, res)
        && this_point_is_on_line_segment(lineSegment1, lineSegment2.getSP(), res))
    {
        return true;
    }

    if (   rg_ZERO(orientationOfEndPtOfLine2, res)
        && this_point_is_on_line_segment(lineSegment1, lineSegment2.getEP(), res))
    {
        return true;
    }

    return false;
}


inline bool BucketForDisk::this_point_is_on_line_segment(const rg_Line2D& lineSegment, const rg_Point2D& endPointOfTheOtherLineSegment, const double& res) const
{
    rg_Point2D startPtOfLineSegment = lineSegment.getSP();
    rg_Point2D endPtOfLineSegment   = lineSegment.getEP();

    double maxX;
    double minX;
    double maxY;
    double minY;

    if (startPtOfLineSegment.getX() < endPtOfLineSegment.getX())
    {
        minX = startPtOfLineSegment.getX();
        maxX = endPtOfLineSegment.getX();
    }
    else
    {
        minX = endPtOfLineSegment.getX();
        maxX = startPtOfLineSegment.getX();
    }

    if (startPtOfLineSegment.getY() < endPtOfLineSegment.getY())
    {
        minY = startPtOfLineSegment.getY();
        maxY = endPtOfLineSegment.getY();
    }
    else
    {
        minY = endPtOfLineSegment.getY();
        maxY = startPtOfLineSegment.getY();
    }


    if (   rg_GE(endPointOfTheOtherLineSegment.getX(), minX, res)
        && rg_LE(endPointOfTheOtherLineSegment.getX(), maxX, res)
        && rg_GE(endPointOfTheOtherLineSegment.getY(), minY, res)
        && rg_LE(endPointOfTheOtherLineSegment.getY(), maxY, res))
    {
        return true;
    }
    else
    {
        return false;
    }
   
};


#endif
