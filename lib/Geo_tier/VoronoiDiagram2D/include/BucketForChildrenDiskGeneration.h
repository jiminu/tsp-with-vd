#ifndef BUCKET_FOR_CHILDREN_DISK_GENERATION_H
#define BUCKET_FOR_CHILDREN_DISK_GENERATION_H

#include "rg_Circle2D.h"
#include "rg_BoundingBox2D.h"
#include "Generator2D.h"

#include <list>
#include <set>
using namespace std;

#define ChildrenDisk pair<rg_Circle2D*, Generator2D*> 

class BucketForChildrenDisks;

class BucketCellForChildrenDisks
{
private:
    list<ChildrenDisk> m_containingDisks;
    

public:
    BucketCellForChildrenDisks();
    BucketCellForChildrenDisks( const BucketCellForChildrenDisks& bucketCell );
    ~BucketCellForChildrenDisks();


    list<ChildrenDisk>&    containingDisks();
    void            insert_disk( const ChildrenDisk& disk );
    void            delete_disk( const ChildrenDisk& disk );


    BucketCellForChildrenDisks& operator =( const BucketCellForChildrenDisks& bucketCell );
};



class BucketForChildrenDisks
{
private:
    BucketCellForChildrenDisks**    m_bucketCell;
    int                             m_numberOfDisks;

    int                             m_numberOfBucketCellOfXDir;
    int                             m_numberOfBucketCellOfYDir;

    double                          m_lengthOfBucketCellOfXDir;
    double                          m_lengthOfBucketCellOfYDir;

    rg_BoundingBox2D                m_boundingBoxOfPVerticesWithSmallEnlargement;

public:
    BucketForChildrenDisks();
    BucketForChildrenDisks(const list<ChildrenDisk>& disks);
    BucketForChildrenDisks(const BucketForChildrenDisks& bucket);
    ~BucketForChildrenDisks();

    int                             number_of_disks() const;
    int                             number_of_bucket_cells_of_x_direction() const;
    int                             number_of_bucket_cells_of_y_direction() const;
    rg_Point2D                      bottom_left_point_of_bounding_box() const;
    rg_Point2D                      top_right_point_of_bounding_box() const;

    BucketCellForChildrenDisks&     get_bucket_cell(const int& xIndex, const int& yIndex);
    BucketCellForChildrenDisks&     find_bucket_cell_containing_this_point(const rg_Point2D& point);
    void                            find_index_of_bucket_cell_containing_this_disk(const ChildrenDisk& disk, int& minXIndex, int& minYIndex, int& maxXIndex, int& maxYIndex);
    void                            find_index_of_bucket_cell_containing_this_disk( const rg_Circle2D& disk, int& minXIndex, int& minYIndex, int& maxXIndex, int& maxYIndex );
    void                            find_index_of_bucket_cell_containing_this_point(const rg_Point2D& point, int& xIndex, int& yIndex);
    void                            collect_candidate_small_disks_to_check_inclusion_by_big_disk( const ChildrenDisk& disk, list<ChildrenDisk>& candidateDisksToBeIncluded );
    void                            collect_candidate_disks_to_check_intersection_by_input_disk( const ChildrenDisk& disk, set<ChildrenDisk>& candidateDisksToBeIntersected );
    void                            collect_candidate_disks_to_check_intersection_by_input_disk( const rg_Circle2D& disk, set<ChildrenDisk>& candidateDisksToBeIntersected );

    void                            set_boundingBox(const list<ChildrenDisk>& disks);
    void                            set_bucket(const list<ChildrenDisk>& disks);

    void                            insertDiskInBucketCell(const ChildrenDisk& disk);
    void                            deleteDIskFromBucketCell(const ChildrenDisk& disk);

private:
    void                            destroyBucket();
    
};










inline BucketCellForChildrenDisks::BucketCellForChildrenDisks()
{

}



inline BucketCellForChildrenDisks::BucketCellForChildrenDisks( const BucketCellForChildrenDisks & bucketCell )
    : m_containingDisks(bucketCell.m_containingDisks)
{

}



inline BucketCellForChildrenDisks::~BucketCellForChildrenDisks()
{

}



inline list<ChildrenDisk>& BucketCellForChildrenDisks::containingDisks()
{
    return m_containingDisks;
}



inline void BucketCellForChildrenDisks::insert_disk( const ChildrenDisk& disk )
{
    m_containingDisks.push_back(disk);
}



inline void BucketCellForChildrenDisks::delete_disk( const ChildrenDisk& disk )
{
    list<ChildrenDisk>::iterator i_disk;
    for ( i_disk = m_containingDisks.begin(); i_disk != m_containingDisks.end(); ++i_disk ) {
        ChildrenDisk currDisk = *i_disk;

        if ( currDisk.first == disk.first ) {
            m_containingDisks.erase(i_disk);
            break;
        }
    }
}



inline BucketCellForChildrenDisks & BucketCellForChildrenDisks::operator=( const BucketCellForChildrenDisks & bucketCell )
{
    if ( this == &bucketCell ) {
        return *this;
    }

    m_containingDisks = bucketCell.m_containingDisks;

    return *this;
}










BucketForChildrenDisks::BucketForChildrenDisks()
    :   m_bucketCell(NULL),
        m_numberOfDisks(0),
        m_numberOfBucketCellOfXDir(0),
        m_numberOfBucketCellOfYDir(0),
        m_lengthOfBucketCellOfXDir(0.0),
        m_lengthOfBucketCellOfYDir(0.0)
{

}



BucketForChildrenDisks::BucketForChildrenDisks( const list<ChildrenDisk>& disks )
    :   m_bucketCell(NULL),
        m_numberOfDisks(0),
        m_numberOfBucketCellOfXDir(0),
        m_numberOfBucketCellOfYDir(0),
        m_lengthOfBucketCellOfXDir(0.0),
        m_lengthOfBucketCellOfYDir(0.0)
{
    set_bucket(disks);
}



BucketForChildrenDisks::BucketForChildrenDisks( const BucketForChildrenDisks& bucket )
    :   m_bucketCell(NULL),
        m_numberOfDisks(bucket.m_numberOfDisks),
        m_numberOfBucketCellOfXDir(bucket.m_numberOfBucketCellOfXDir),
        m_numberOfBucketCellOfYDir(bucket.m_numberOfBucketCellOfYDir),
        m_lengthOfBucketCellOfXDir(bucket.m_lengthOfBucketCellOfXDir),
        m_lengthOfBucketCellOfYDir(bucket.m_lengthOfBucketCellOfYDir)
{
    m_bucketCell = new BucketCellForChildrenDisks*[m_numberOfBucketCellOfXDir];
    for ( int i = 0; i < m_numberOfBucketCellOfXDir; ++i ) {
        m_bucketCell[i] = new BucketCellForChildrenDisks[m_numberOfBucketCellOfYDir];
    }

    for ( int i = 0; i < m_numberOfBucketCellOfXDir; ++i ) {
        for ( int j = 0; j < m_numberOfBucketCellOfYDir; ++j ) {
            m_bucketCell[i][j] = bucket.m_bucketCell[i][j];
        }
    }
}



BucketForChildrenDisks::~BucketForChildrenDisks( )
{
    destroyBucket();
}



inline int BucketForChildrenDisks::number_of_disks() const
{
    return m_numberOfDisks;
}



inline int BucketForChildrenDisks::number_of_bucket_cells_of_x_direction() const
{
    return m_numberOfBucketCellOfXDir;
}



inline int BucketForChildrenDisks::number_of_bucket_cells_of_y_direction() const
{
    return m_numberOfBucketCellOfYDir;
}



inline rg_Point2D BucketForChildrenDisks::bottom_left_point_of_bounding_box() const
{
    return m_boundingBoxOfPVerticesWithSmallEnlargement.getMinPt();
}



inline rg_Point2D BucketForChildrenDisks::top_right_point_of_bounding_box() const
{
    return m_boundingBoxOfPVerticesWithSmallEnlargement.getMaxPt();
}



inline BucketCellForChildrenDisks & BucketForChildrenDisks::get_bucket_cell( const int & xIndex, const int & yIndex )
{
    return m_bucketCell[xIndex][yIndex];
}



inline BucketCellForChildrenDisks & BucketForChildrenDisks::find_bucket_cell_containing_this_point( const rg_Point2D & point )
{
    int xIndex = 0;
    int yIndex = 0;

    find_index_of_bucket_cell_containing_this_point(point, xIndex, yIndex);

    return get_bucket_cell(xIndex, yIndex);
}



inline void BucketForChildrenDisks::find_index_of_bucket_cell_containing_this_disk(const ChildrenDisk& disk, int& minXIndex, int& minYIndex, int& maxXIndex, int& maxYIndex)
{
    rg_Point2D minPtOfBoundingBoxOfInputDisk( disk.first->getCenterPt().getX() - disk.first->getRadius(), disk.first->getCenterPt().getY() - disk.first->getRadius() );
    rg_Point2D maxPtOfBoundingBoxOfInputDisk( disk.first->getCenterPt().getX() + disk.first->getRadius(), disk.first->getCenterPt().getY() + disk.first->getRadius() );

    find_index_of_bucket_cell_containing_this_point( minPtOfBoundingBoxOfInputDisk, minXIndex, minYIndex );
    find_index_of_bucket_cell_containing_this_point( maxPtOfBoundingBoxOfInputDisk, maxXIndex, maxYIndex );
}



inline void BucketForChildrenDisks::find_index_of_bucket_cell_containing_this_disk( const rg_Circle2D& disk, int& minXIndex, int& minYIndex, int& maxXIndex, int& maxYIndex )
{
    rg_Point2D minPtOfBoundingBoxOfInputDisk( disk.getCenterPt().getX() - disk.getRadius(), disk.getCenterPt().getY() - disk.getRadius() );
    rg_Point2D maxPtOfBoundingBoxOfInputDisk( disk.getCenterPt().getX() + disk.getRadius(), disk.getCenterPt().getY() + disk.getRadius() );

    find_index_of_bucket_cell_containing_this_point( minPtOfBoundingBoxOfInputDisk, minXIndex, minYIndex );
    find_index_of_bucket_cell_containing_this_point( maxPtOfBoundingBoxOfInputDisk, maxXIndex, maxYIndex );
}



inline void BucketForChildrenDisks::find_index_of_bucket_cell_containing_this_point(const rg_Point2D& point, int& xIndex, int& yIndex)
{
    xIndex = (int)( (point.getX() - m_boundingBoxOfPVerticesWithSmallEnlargement.getMinPt().getX()) / m_lengthOfBucketCellOfXDir );
    yIndex = (int)( (point.getY() - m_boundingBoxOfPVerticesWithSmallEnlargement.getMinPt().getY()) / m_lengthOfBucketCellOfYDir );

    if( xIndex < 0 )
        xIndex = 0;
    if( xIndex > m_numberOfBucketCellOfXDir - 1 )
        xIndex = m_numberOfBucketCellOfXDir - 1;

    if( yIndex < 0 )
        yIndex = 0;
    if( yIndex > m_numberOfBucketCellOfYDir - 1 )
        yIndex = m_numberOfBucketCellOfYDir - 1;
}



inline void BucketForChildrenDisks::collect_candidate_small_disks_to_check_inclusion_by_big_disk( const ChildrenDisk& disk, list<ChildrenDisk>& candidateDisksToBeIncluded )
{
    int minXIndex, minYIndex, maxXIndex, maxYIndex;
    find_index_of_bucket_cell_containing_this_disk(disk, minXIndex, minYIndex, maxXIndex, maxYIndex);

    for ( int i_x = minXIndex; i_x <= maxXIndex; ++i_x ) {
        for ( int i_y = minYIndex; i_y <= maxYIndex; ++i_y ) {
            list<ChildrenDisk> containingDisks = get_bucket_cell(i_x, i_y).containingDisks();

            candidateDisksToBeIncluded.insert(candidateDisksToBeIncluded.end(), containingDisks.begin(), containingDisks.end());
        }
    }
}



inline void BucketForChildrenDisks::collect_candidate_disks_to_check_intersection_by_input_disk( const ChildrenDisk& disk, set<ChildrenDisk>& candidateDisksToBeIntersected )
{
    int minXIndex, minYIndex, maxXIndex, maxYIndex;
    find_index_of_bucket_cell_containing_this_disk(disk, minXIndex, minYIndex, maxXIndex, maxYIndex);

    for ( int i_x = minXIndex; i_x <= maxXIndex; ++i_x ) {
        for ( int i_y = minYIndex; i_y <= maxYIndex; ++i_y ) {
            list<ChildrenDisk> containingDisks = get_bucket_cell(i_x, i_y).containingDisks();

            candidateDisksToBeIntersected.insert(containingDisks.begin(), containingDisks.end());
                //.insert(candidateDisksToBeIncluded.end(), containingDisks.begin(), containingDisks.end());
        }
    }
}



inline void BucketForChildrenDisks::collect_candidate_disks_to_check_intersection_by_input_disk( const rg_Circle2D& disk, set<ChildrenDisk>& candidateDisksToBeIntersected )
{
    int minXIndex, minYIndex, maxXIndex, maxYIndex;
    find_index_of_bucket_cell_containing_this_disk( disk, minXIndex, minYIndex, maxXIndex, maxYIndex );

    for ( int i_x = minXIndex; i_x <= maxXIndex; ++i_x ) {
        for ( int i_y = minYIndex; i_y <= maxYIndex; ++i_y ) {
            list<ChildrenDisk> containingDisks = get_bucket_cell( i_x, i_y ).containingDisks();

            candidateDisksToBeIntersected.insert( containingDisks.begin(), containingDisks.end() );
            //.insert(candidateDisksToBeIncluded.end(), containingDisks.begin(), containingDisks.end());
        }
    }
}



inline void BucketForChildrenDisks::set_boundingBox( const list<ChildrenDisk>& disks )
{
    list<ChildrenDisk>::const_iterator i_disk;
    for ( i_disk = disks.begin(); i_disk != disks.end(); ++i_disk ) {
        const ChildrenDisk currDisk = *i_disk;
        m_boundingBoxOfPVerticesWithSmallEnlargement.updateBoxByAddingCircle(rg_Circle2D(currDisk.first->getCenterPt(), 0.0));
        //m_boundingBoxOfPVerticesWithSmallEnlargement.updateBoxByAddingCircle(currDisk);
    }
}



inline void BucketForChildrenDisks::set_bucket( const list<ChildrenDisk>& disks )
{
    set_boundingBox(disks);

    m_numberOfDisks = disks.size();

    double xDirLengthOfBoundingBox    = m_boundingBoxOfPVerticesWithSmallEnlargement.evaluateXLength();
    double yDirLengthOfBoundingBox    = m_boundingBoxOfPVerticesWithSmallEnlargement.evaluateYLength();
    double sumOfLengthOfXAxis_N_YAxis = xDirLengthOfBoundingBox + yDirLengthOfBoundingBox;

    m_numberOfBucketCellOfXDir = (int)( (double)m_numberOfDisks / 2.0 * ( xDirLengthOfBoundingBox / sumOfLengthOfXAxis_N_YAxis ) );
    m_numberOfBucketCellOfYDir = (int)( (double)m_numberOfDisks / 2.0 * ( yDirLengthOfBoundingBox / sumOfLengthOfXAxis_N_YAxis ) );

    m_lengthOfBucketCellOfXDir = xDirLengthOfBoundingBox / (double)m_numberOfBucketCellOfXDir;
    m_lengthOfBucketCellOfYDir = yDirLengthOfBoundingBox / (double)m_numberOfBucketCellOfYDir;

    if( m_lengthOfBucketCellOfXDir > m_lengthOfBucketCellOfYDir ) {
        m_lengthOfBucketCellOfYDir = m_lengthOfBucketCellOfXDir;
    }
    else {
        m_lengthOfBucketCellOfXDir = m_lengthOfBucketCellOfYDir;
    }


    m_bucketCell = new BucketCellForChildrenDisks*[m_numberOfBucketCellOfXDir];
    for ( int i = 0; i < m_numberOfBucketCellOfXDir; ++i ) {
        m_bucketCell[i] = new BucketCellForChildrenDisks[m_numberOfBucketCellOfYDir];
    }


    // loop for disk
    list<ChildrenDisk>::const_iterator i_disk;
    for ( i_disk = disks.begin(); i_disk != disks.end(); ++i_disk ) {
        const ChildrenDisk currDisk = *i_disk;
        insertDiskInBucketCell(currDisk);
    }    
}



inline void BucketForChildrenDisks::insertDiskInBucketCell( const ChildrenDisk& disk )
{
    int minXIndex, minYIndex, maxXIndex, maxYIndex;
    find_index_of_bucket_cell_containing_this_disk( disk, minXIndex, minYIndex, maxXIndex, maxYIndex );

    // loop for xAxis
    for ( int i_x = minXIndex; i_x <= maxXIndex; ++i_x ) {

        // loop for xAxis
        for ( int i_y = minYIndex; i_y <= maxYIndex; ++i_y ) {
            m_bucketCell[i_x][i_y].insert_disk(disk);
        }
    }
}



inline void BucketForChildrenDisks::deleteDIskFromBucketCell( const ChildrenDisk& disk )
{
    int minXIndex, minYIndex, maxXIndex, maxYIndex;
    find_index_of_bucket_cell_containing_this_disk( disk, minXIndex, minYIndex, maxXIndex, maxYIndex );

    // loop for xAxis
    for ( int i_x = minXIndex; i_x <= maxXIndex; ++i_x ) {

        // loop for xAxis
        for ( int i_y = minYIndex; i_y <= maxYIndex; ++i_y ) {
            m_bucketCell[i_x][i_y].delete_disk(disk);
        }
    }
}



inline void BucketForChildrenDisks::destroyBucket()
{
    if( m_bucketCell != NULL ) {
        for( int i = 0 ; i < m_numberOfBucketCellOfXDir ; i++ ) {
            if( m_bucketCell[i] != NULL ) {
                delete [] m_bucketCell[i];
            }
        }

        delete [] m_bucketCell;
    }


    m_bucketCell = NULL;

    m_numberOfBucketCellOfXDir = 0;
    m_numberOfBucketCellOfYDir = 0;

    m_lengthOfBucketCellOfXDir = 0.0;;
    m_lengthOfBucketCellOfYDir = 0.0;

    m_boundingBoxOfPVerticesWithSmallEnlargement = rg_BoundingBox2D();
}



#endif


