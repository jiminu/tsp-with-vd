#include"BucketForDisk.h"
#include "rg_Point2D.h"
#include "rg_GeoFunc.h"
#include <unordered_set>
#include <stack>
using namespace std;

BucketForDisk::BucketForDisk()
: m_numBucketCellOfXDir(0),
  m_numBucketCellOfYDir(0),
  m_xDirLengthOfBucketCell(0.0),
  m_yDirLengthOfBucketCell(0.0)
  //m_bucketCell(NULL)
{
}


BucketForDisk::BucketForDisk(const BucketForDisk& BFD)
{
    copyBucketFrom(BFD);
}


BucketForDisk::BucketForDisk(const list<rg_Circle2D>& disks, const CellDivisionOptions& divisionOptions /* = BY_LARGEST_RADIUS*/, const bool& DisksAreArrangedInBucket /* = true*/)//const list<Generator2D>& diskSet___zzfeng
{
    list<rg_Circle2D*> diskPtrs;
    for (list<rg_Circle2D>::const_iterator it_Disk = disks.begin(); it_Disk != disks.end(); ++it_Disk)
    {
        diskPtrs.push_back(const_cast<rg_Circle2D*>(&(*it_Disk)));
    }

    generateBucketForDisks(diskPtrs, divisionOptions, DisksAreArrangedInBucket);
}


BucketForDisk::BucketForDisk(const list<rg_Circle2D*>& disks, const CellDivisionOptions& divisionOptions /*= BY_LARGEST_RADIUS*/, const bool& DisksAreArrangedInBucket /* = true*/)
{
    generateBucketForDisks(disks, divisionOptions, DisksAreArrangedInBucket);
}


BucketForDisk::~BucketForDisk(void)
{
    destroyBucket();
}


double BucketForDisk::findMaxRadius(const list<rg_Circle2D*>& disks) const
{
    double maxRadius = 0.0;
    for (list<rg_Circle2D*>::const_iterator it_Disk = disks.begin(); it_Disk != disks.end(); ++it_Disk)
    {
        double currDiskRadius = (*it_Disk)->getRadius();

        if (currDiskRadius > maxRadius)
        {
            maxRadius = currDiskRadius;
        }
    }
    
    return maxRadius;
}

double BucketForDisk::findMinRadius(const list<rg_Circle2D*>& disks) const
{
    double minRadius = DBL_MAX;
    for (list<rg_Circle2D*>::const_iterator it_Disk = disks.begin(); it_Disk != disks.end(); ++it_Disk)
    {
        double currDiskRadius = (*it_Disk)->getRadius();

        if (currDiskRadius < minRadius)
        {
            minRadius = currDiskRadius;
        }
    }

    return minRadius;
}


void BucketForDisk::setBucketCellSize(const double& cellSize)
{
    double xDirLengthOfBoundingBox = m_boundingBox.getMaxPt().getX() - m_boundingBox.getMinPt().getX();
    double yDirLengthOfBoundingBox = m_boundingBox.getMaxPt().getY() - m_boundingBox.getMinPt().getY();

	m_xDirLengthOfBucketCell = cellSize;
	m_yDirLengthOfBucketCell = cellSize;

	m_numBucketCellOfXDir = ceil(xDirLengthOfBoundingBox / m_xDirLengthOfBucketCell);
	m_numBucketCellOfYDir = ceil(yDirLengthOfBoundingBox / m_yDirLengthOfBucketCell);
	
	double maxPtX = m_boundingBox.getMinPt().getX() + m_xDirLengthOfBucketCell * m_numBucketCellOfXDir;
	double maxPtY = m_boundingBox.getMinPt().getY() + m_yDirLengthOfBucketCell * m_numBucketCellOfYDir;
    m_boundingBox.setMaxPt(rg_Point2D(maxPtX, maxPtY));

	//int j = 0;
    m_bucketCell.resize(m_numBucketCellOfXDir);
    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        m_bucketCell[i].resize(m_numBucketCellOfYDir, NULL);
    }
    /*
    m_bucketCell = new BucketCellForDisks*[m_numBucketCellOfXDir];
    for( int i = 0 ; i < m_numBucketCellOfXDir ; ++i)
    {
        m_bucketCell[i] = new BucketCellForDisks[m_numBucketCellOfYDir];
    }
    */
    
    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        for (int j = 0; j < m_numBucketCellOfYDir; ++j)
        {
            m_bucketCell[i][j] = new BucketCellForDisks();
            m_bucketCell[i][j]->set_Xindex(i);
            m_bucketCell[i][j]->set_Yindex(j);
        }
    }
}


void BucketForDisk::setBucketCellNumbers(const int& numX, const int& numY)
{
    double xDirLengthOfBoundingBox = m_boundingBox.getMaxPt().getX() - m_boundingBox.getMinPt().getX();
    double yDirLengthOfBoundingBox = m_boundingBox.getMaxPt().getY() - m_boundingBox.getMinPt().getY();

    m_xDirLengthOfBucketCell = xDirLengthOfBoundingBox/numX;
    m_yDirLengthOfBucketCell = yDirLengthOfBoundingBox/numY;

    m_numBucketCellOfXDir = numX;
    m_numBucketCellOfYDir = numY;

    //int j = 0;
    m_bucketCell.resize(m_numBucketCellOfXDir);
    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        m_bucketCell[i].resize(m_numBucketCellOfYDir);
    }
    /*
        m_bucketCell = new BucketCellForDisks*[m_numBucketCellOfXDir];
        for (int i = 0; i < m_numBucketCellOfXDir; ++i)
        {
            m_bucketCell[i] = new BucketCellForDisks[m_numBucketCellOfYDir];
        }
    */


    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        for (int j = 0; j < m_numBucketCellOfYDir; ++j)
        {
            m_bucketCell[i][j] = new BucketCellForDisks();
            m_bucketCell[i][j]->set_Xindex(i);
            m_bucketCell[i][j]->set_Yindex(j);
        }
    }
}

void BucketForDisk::setBoundingBox(const list<rg_Circle2D*>& disks)
{
    for( list<rg_Circle2D*>::const_iterator it_Disk = disks.begin() ; it_Disk != disks.end() ; ++it_Disk )
	{
		const rg_Circle2D* currDisk = *it_Disk;    
        m_boundingBox.updateBoxByAddingCircle( *currDisk );
    }
}


void BucketForDisk::setBoundingBox(const rg_Point2D& maxPt, const rg_Point2D& minPt)
{
    m_boundingBox.setMaxPt(maxPt);
    m_boundingBox.setMinPt(minPt);
}


void BucketForDisk::addDisksIntoBucket(const list<rg_Circle2D*>& disks)
{
    for (list<rg_Circle2D*>::const_iterator it_Disk = disks.begin(); it_Disk != disks.end(); ++it_Disk)
    {
        rg_Circle2D* currDisk = const_cast<rg_Circle2D*>(*it_Disk);
        addDiskIntoBucket(currDisk);
    }
}

void BucketForDisk::generateBucketForDisks(const list<rg_Circle2D*>& disks, const CellDivisionOptions& divisionOptions, const bool& DisksAreArrangedInBucket)
{
    setBoundingBox(disks);

    double cellSize = 0.5;
    switch (divisionOptions)
    {
        case BY_LARGEST_RADIUS:
        {
            double maxRaidus = findMaxRadius(disks);
            cellSize = 2 * maxRaidus + 0.5;
            setBucketCellSize(cellSize);

            break;
        }

        case BY_SMALLEST_RAIDUS:
        {
            double minRadius = findMinRadius(disks);
            cellSize = 2 * minRadius + 0.5;
            setBucketCellSize(cellSize);

            break;
        }

        case BY_FIX_NUM_100_N_100:
        {
            setBucketCellNumbers(100, 100);
            break;
        }

        case NUM_OF_DISKS_SAME_AS_NUM_OF_CELLS:
        {
            int numOfCellsInOneAxis = ceil(sqrt(disks.size()));
            setBucketCellNumbers(numOfCellsInOneAxis, numOfCellsInOneAxis);
            break;
        }
    }

    if (DisksAreArrangedInBucket)
    {
        addDisksIntoBucket(disks);
    }
}




BucketCellForDisks* BucketForDisk::getBucketCellIncludingPoint(const rg_Point2D& pt) const
{
    if (!m_boundingBox.doesContain(pt, RES_OF_POINT_IS_CONTAINED_IN_THIS_BUCKET))
    {
        return NULL;
    }

    int Xindex, Yindex;
    BucketCellForDisks* targetBucketCell = NULL;

    getBucketCellIndexForPoint(pt, Xindex, Yindex);
    targetBucketCell = getBucketCellAccordingIndex(Xindex, Yindex);

    return targetBucketCell;
}


BucketCellForDisks* BucketForDisk::getBucketCellContainingThisDisk(const rg_Circle2D* circle) const
{
	int Xindex, Yindex;
	BucketCellForDisks* targetBucketCell = NULL;

	getBucketCellIndexForPoint(circle->getCircle().getCenterPt(), Xindex, Yindex);
	targetBucketCell = getBucketCellAccordingIndex(Xindex, Yindex);

	//possibleBucketCell¼öÁ¤
	bool targetBucketCellFound = false;

	const list<rg_Circle2D*>& dynamicDiskslist = targetBucketCell->get_contained_disks();

	for (list<rg_Circle2D*>::const_iterator it_Disk = dynamicDiskslist.begin(); it_Disk != dynamicDiskslist.end(); it_Disk++)
	{
		rg_Circle2D* currDynamicDisk = const_cast<rg_Circle2D*>(*it_Disk);

		if (circle == currDynamicDisk)
		{
			targetBucketCellFound = true;
		}

	}

	if (targetBucketCellFound == false)
	{
		list<BucketCellForDisks*> max9Cells;
		getMax9BucketCellsIncludingNAroundThisCell(targetBucketCell, max9Cells);

		for (list<BucketCellForDisks*>::const_iterator it_bucketCell = max9Cells.begin(); it_bucketCell != max9Cells.end(); it_bucketCell++)
		{
			BucketCellForDisks* currBucketCell = const_cast<BucketCellForDisks*>(*it_bucketCell);

			const list<rg_Circle2D*>&  disksInCurrBucket = currBucketCell->get_contained_disks();

			for (list<rg_Circle2D*>::const_iterator it_Disk = disksInCurrBucket.begin(); it_Disk != disksInCurrBucket.end(); it_Disk++)
			{
				rg_Circle2D* currDynamicDisk = const_cast<rg_Circle2D*>(*it_Disk);

				if (circle == currDynamicDisk)
				{
					targetBucketCell      = currBucketCell;
					targetBucketCellFound = true;
					break;
				}
			}

			if (targetBucketCellFound == true)
			{
				break;
			}
		}
	}

	return targetBucketCell;
}


void BucketForDisk::getMax9BucketCellsIncludingNAroundThisCell(BucketCellForDisks* cell, list<BucketCellForDisks*>& max9Cells) const
{
	int Xindex = cell->get_Xindex();
	int Yindex = cell->get_Yindex();
	
    if (m_numBucketCellOfXDir == 1)
    {
        if (m_numBucketCellOfYDir == 1)
        {
            max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
        }
        else if (m_numBucketCellOfYDir == 2)
        {
            max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
            max9Cells.push_back(getBucketCellAccordingIndex(0, 1));
        }
        else //(m_numBucketCellOfYDir >= 3)
        {
            if (Yindex == 0)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(0, 1));
            }
            else if (Yindex == m_numBucketCellOfYDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex));
            }
            else
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex + 1));
            }
        }
    }
    else if (m_numBucketCellOfXDir == 2)
    {
        if (m_numBucketCellOfYDir == 1)
        {
            max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
            max9Cells.push_back(getBucketCellAccordingIndex(1, 0));
        }
        else if (m_numBucketCellOfYDir == 2)
        {
            max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
            max9Cells.push_back(getBucketCellAccordingIndex(0, 1));
            max9Cells.push_back(getBucketCellAccordingIndex(1, 0));
            max9Cells.push_back(getBucketCellAccordingIndex(1, 1));
        }
        else //(m_numBucketCellOfYDir >= 3)
        {
            if (Yindex == 0)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(0, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 1));
            }
            else if (Yindex == m_numBucketCellOfYDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex));
            }
            else
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex + 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex + 1));
            }
        }
    }
    else // m_numBucketCellOfXDir >= 3
    {
        if (m_numBucketCellOfYDir == 1)
        {
            if (Xindex == 0)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 0));
            }
            else if (Xindex == m_numBucketCellOfXDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex,     0));
            }
            else
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex,     0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, 0));
            }
        }
        else if (m_numBucketCellOfYDir == 2)
        {
            if (Xindex == 0)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(0, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 1));

            }
            else if (Xindex == m_numBucketCellOfXDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 1));
            }
            else
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, 1));
            }
        }
        else // m_numBucketCellOfYDir >= 3
        {
            if (Xindex == 0 && Yindex == 0)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(0, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, 1));
            }
            else if (Xindex == m_numBucketCellOfXDir - 1 && Yindex == 0)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 1));
            }
            else if (Xindex == 0 && Yindex == m_numBucketCellOfYDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex));
            }
            else if (Xindex == m_numBucketCellOfXDir - 1 && Yindex == m_numBucketCellOfYDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex));
            }
            else if (0 < Xindex < m_numBucketCellOfXDir - 1 && Yindex == 0)//1
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, 0));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, 1));
            }
            else if (Xindex == 0 && 0 < Yindex < m_numBucketCellOfYDir - 1)//2
            {
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(0, Yindex + 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(1, Yindex + 1));
            }
            else if (Xindex == m_numBucketCellOfXDir - 1 && 0 < Yindex < m_numBucketCellOfYDir - 1)//3
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex + 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex + 1));
            }
            else if (0 < Xindex < m_numBucketCellOfXDir - 1 && Yindex == m_numBucketCellOfYDir - 1)//4
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, Yindex - 1));
            }
            else if (0 < Xindex < m_numBucketCellOfXDir - 1 && 0 < Yindex < m_numBucketCellOfYDir - 1)
            {
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex - 1, Yindex + 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex, Yindex + 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, Yindex - 1));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, Yindex));
                max9Cells.push_back(getBucketCellAccordingIndex(Xindex + 1, Yindex + 1));
            }
        }
    }
}



void BucketForDisk::destroyBucket()
{
    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        for (int j = 0; j < m_numBucketCellOfYDir; ++j)
        {
            delete m_bucketCell[i][j];
        }
    }

    m_bucketCell.resize(0);


    /*
        for( int i = 0 ; i < m_numBucketCellOfXDir ; i++ ) {
            if( m_bucketCell[i] != NULL ) {
                delete [] m_bucketCell[i];
            }
        }

        if( m_bucketCell != NULL ) {
            delete [] m_bucketCell;
        }

        m_bucketCell = NULL;
    */

}


void BucketForDisk::clear()
{
    destroyBucket();
    //m_bucketCell          = NULL;
    m_numBucketCellOfXDir = 0;
    m_numBucketCellOfYDir = 0;
    m_xDirLengthOfBucketCell = 0.0;
    m_yDirLengthOfBucketCell = 0.0;
}


void BucketForDisk::copyBucketFrom(const BucketForDisk& BFD)
{
    m_numBucketCellOfXDir		= BFD.m_numBucketCellOfXDir;
	m_numBucketCellOfYDir		= BFD.m_numBucketCellOfYDir;
	m_xDirLengthOfBucketCell	= BFD.m_xDirLengthOfBucketCell;
	m_yDirLengthOfBucketCell	= BFD.m_yDirLengthOfBucketCell;
	m_boundingBox				= BFD.m_boundingBox;
    m_DisksInBucket             = BFD.m_DisksInBucket;

    m_bucketCell.resize(m_numBucketCellOfXDir);
    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        m_bucketCell[i].resize(m_numBucketCellOfYDir);
    }

    /*
        m_bucketCell = new BucketCellForDisks*[m_numBucketCellOfXDir];
        for (int i = 0; i < m_numBucketCellOfXDir; ++i)
        {
            m_bucketCell[i] = new BucketCellForDisks[m_numBucketCellOfYDir];
        }
    */

    for (int i = 0; i < m_numBucketCellOfXDir; ++i)
    {
        for (int j = 0; j < m_numBucketCellOfYDir; ++j)
        {
            m_bucketCell[i][j] = BFD.m_bucketCell[i][j];
        }
    }
}

BucketForDisk& BucketForDisk:: operator=(const BucketForDisk& BFD)
{
	if (this != &BFD)
	{
        clear();
        copyBucketFrom(BFD);
	}

	return *this;
}



void BucketForDisk::findAllCellsIntersectedWithThisLineSegment(
    const rg_Line2D& lineSegment, list<BucketCellForDisks*>& intersectedBucketCells) const
{
    //1. find start cell
    BucketCellForDisks* startCell = NULL;

    if (m_boundingBox.doesContain(lineSegment.getSP(), RES_OF_POINT_IS_CONTAINED_IN_THIS_BUCKET))
    {
        startCell = getBucketCellIncludingPoint(lineSegment.getSP());
    }
    else if (m_boundingBox.doesContain(lineSegment.getEP(), RES_OF_POINT_IS_CONTAINED_IN_THIS_BUCKET))
    {
        startCell = getBucketCellIncludingPoint(lineSegment.getEP());
    }
    else
    {
        pair<bool, rg_Point2D> intersectionPointBtwLineSegementNBoundingBox = 
            thisLineSementIntersectsThisBox(lineSegment, m_boundingBox.getMinPt(), m_boundingBox.getMaxPt());

        bool thisLineSegmentIntersectsBoundingBox = intersectionPointBtwLineSegementNBoundingBox.first;

        if (thisLineSegmentIntersectsBoundingBox)
        {
            startCell = getBucketCellIncludingPoint(intersectionPointBtwLineSegementNBoundingBox.second);
        }
    }

    //2. propagate cells
    unordered_set<BucketCellForDisks*> visitedCells;
    stack<BucketCellForDisks*>         intersectedCellsInStack;

    intersectedCellsInStack.push(startCell);

    while (!intersectedCellsInStack.empty())
    {
        BucketCellForDisks* currCell = intersectedCellsInStack.top();
        intersectedCellsInStack.pop();

        // currCell is already visited.
        if (visitedCells.find(currCell) != visitedCells.end())
        {
            continue;
        }
        else
        {
            visitedCells.insert(currCell);
        }

        rg_Point2D maxPtOfCell;
        rg_Point2D minPtOfCell; 
        getMaxNMinPtsOfThisCell(currCell->get_Xindex(), currCell->get_Yindex(), maxPtOfCell, minPtOfCell);
        
        vector<rg_Line2D> fourBoundaries_buttom_right_top_left_order;
        compute_four_boundary_line_segments(minPtOfCell, maxPtOfCell, fourBoundaries_buttom_right_top_left_order);

        for (int i = 0; i < 4; ++i)
        {
            const rg_Line2D& currBoundary = fourBoundaries_buttom_right_top_left_order[i];

            if (two_line_segments_are_intersected(currBoundary, lineSegment, RES_OF_LINE_SEGMENT))
            {
                switch (i)
                {
                    case 0:
                    { 
                        if (currCell->get_Yindex() > 0)
                        {
                            BucketCellForDisks* nextBucketCell = getBucketCellAccordingIndex(
                                currCell->get_Xindex(), currCell->get_Yindex() - 1);
                            intersectedCellsInStack.push(nextBucketCell);
                        }

                        break;
                    }

                    case 1:
                    {
                        if (currCell->get_Xindex() < (m_numBucketCellOfXDir - 1))
                        {
                            BucketCellForDisks* nextBucketCell = getBucketCellAccordingIndex(
                                currCell->get_Xindex() + 1, currCell->get_Yindex());
                            intersectedCellsInStack.push(nextBucketCell);
                        }

                        break;
                    }

                    case 2:
                    {
                        if (currCell->get_Yindex() < (m_numBucketCellOfYDir - 1))
                        {
                            BucketCellForDisks* nextBucketCell = getBucketCellAccordingIndex(
                                currCell->get_Xindex(), currCell->get_Yindex()+1);
                            intersectedCellsInStack.push(nextBucketCell);
                        }

                        break;
                    }

                    case 3:
                    {
                        if (currCell->get_Xindex() > 0)
                        {
                            BucketCellForDisks* nextBucketCell = getBucketCellAccordingIndex(
                                currCell->get_Xindex() - 1, currCell->get_Yindex());
                            intersectedCellsInStack.push(nextBucketCell);
                        }

                        break;
                    }
                }
            }
        }
    }

    intersectedBucketCells.insert(intersectedBucketCells.end(), visitedCells.begin(), visitedCells.end());
}


void BucketForDisk::getMax9BucketCellsIncludingNAroundThisCells(const list<BucketCellForDisks*>& inputCells, list<BucketCellForDisks*>& cellsAroundNIncludingInputCells) const
{
    unordered_set<BucketCellForDisks*> allCells;
    for (list<BucketCellForDisks*>::const_iterator it_Cell = inputCells.begin();
         it_Cell != inputCells.end();
         ++it_Cell)
    {
        BucketCellForDisks* currCell = *it_Cell;
        list<BucketCellForDisks*> allAroundCellsIncludingCurrCell;
        getMax9BucketCellsIncludingNAroundThisCell(currCell, allAroundCellsIncludingCurrCell);
        allCells.insert(allAroundCellsIncludingCurrCell.begin(), allAroundCellsIncludingCurrCell.end());
    }

    cellsAroundNIncludingInputCells.insert(cellsAroundNIncludingInputCells.end(), allCells.begin(), allCells.end());
}


void BucketForDisk::extractAllDisksIncludedInTheseCells(const list<BucketCellForDisks*>& bucketCells, list<rg_Circle2D*>& disks) const
{
    for (list<BucketCellForDisks*>::const_iterator it_Cell = bucketCells.begin();
         it_Cell != bucketCells.end();
         ++it_Cell)
    {
        const list<rg_Circle2D*>& containedDisks = (*it_Cell)->get_contained_disks();
        disks.insert(disks.end(), containedDisks.begin(), containedDisks.end());
    }
}


bool BucketForDisk::there_is_a_pair_of_disks_intersected(const double& tolerance /*= rg_MATH_RES*/) const
{
    bool thereExistIntersection = false;

    for (list<rg_Circle2D*>::const_iterator it_Disk = m_DisksInBucket.begin(); it_Disk != m_DisksInBucket.end(); ++it_Disk)
    {
        rg_Circle2D* currDisk = *it_Disk;

        BucketCellForDisks* bucketCellIncludingThisPt = getBucketCellIncludingPoint(currDisk->getCenterPt());

        list<BucketCellForDisks*> bucketCells;
        getMax9BucketCellsIncludingNAroundThisCell(bucketCellIncludingThisPt, bucketCells);

        for (list<BucketCellForDisks*>::const_iterator it_Bucket = bucketCells.begin(); it_Bucket != bucketCells.end(); ++it_Bucket)
        {
            BucketCellForDisks* currBucket = const_cast<BucketCellForDisks*> (*it_Bucket);

            const list<rg_Circle2D*>& containedDisksInCurrBucket = currBucket->get_contained_disks();

            for (list<rg_Circle2D*>::const_iterator it_NeighborDiskInBucketCell = containedDisksInCurrBucket.begin(); it_NeighborDiskInBucketCell != containedDisksInCurrBucket.end(); ++it_NeighborDiskInBucketCell)
            {
                rg_Circle2D* neighborDiskInBucket = *it_NeighborDiskInBucketCell;

                if (currDisk == neighborDiskInBucket)
                {
                    continue;
                }

                if (currDisk->isIntersectWith(*neighborDiskInBucket, tolerance))
                {
                    thereExistIntersection = true;
                    break;
                }
            }

            if (thereExistIntersection)
            {
                break;
            }
        }
    }

    return thereExistIntersection;
}
