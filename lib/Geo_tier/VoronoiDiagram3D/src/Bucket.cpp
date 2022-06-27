#include "Bucket.h"
#include <float.h>
#include "rg_Generator.h"
#include "rg_BallGenerator.h"

#include "rg_RelativeOp.h"
#include "FunctionsForVoronoiDiagram3D.h"
using namespace V::GeometryTier;



Bucket::Bucket()
: m_bucketTable( rg_NULL ),
  m_bucketMark(  rg_NULL ),
  m_numOfXs(0),
  m_numOfYs(0),
  m_numOfZs(0),
  m_diagonalLength(0.0),
  m_minPointOfBoundingBox(rg_MAX_REAL, rg_MAX_REAL, rg_MAX_REAL), 
  m_maxPointOfBoundingBox(-rg_MAX_REAL, -rg_MAX_REAL, -rg_MAX_REAL)
{
}



Bucket::Bucket(rg_dList< BallGenerator >& ballList)
{
	constructBoundingBox( ballList );
	
	//bucket의 대략적인 size : 현재는 generator개수의 두배.
	rg_INT N = ballList.getSize() * 2;

	//size of bounding box
	rg_REAL xLength = m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX();
	rg_REAL yLength = m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY();
	rg_REAL zLength = m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ();

	m_numOfZs = (rg_INT)( pow( (zLength*zLength)/(xLength*yLength) * N, 1./3 ) + 1. );
	m_numOfYs = (rg_INT)( yLength/zLength * m_numOfZs + 1. );
	m_numOfXs = (rg_INT)( xLength/zLength * m_numOfZs + 1. );

	m_bucketTable = new rg_dList<BallGenerator*>** [m_numOfXs];
    m_bucketMark  = new rg_FLAG** [m_numOfXs];
	for(rg_INT i=0; i<m_numOfXs; i++)  {
		m_bucketTable[i] = new rg_dList<BallGenerator*>* [m_numOfYs];
		m_bucketMark[i]  = new rg_FLAG* [m_numOfYs];

        for(rg_INT j=0; j<m_numOfYs; j++)  {
			m_bucketTable[i][j] = new rg_dList<BallGenerator*> [m_numOfZs];
			m_bucketMark[i][j]  = new rg_FLAG [m_numOfZs];

            for ( rg_INT k=0; k<m_numOfZs; k++)  {
                m_bucketMark[i][j][k] = 0; 
            }
		}
	}

	//insert cell into bucket structure
	ballList.reset4Loop();
	while(ballList.setNext4Loop())
	{
		addCellIntoBucket( ballList.getpEntity() );
	}

	m_diagonalLength = sqrt(   pow(xLength/m_numOfXs, 2) 
                             + pow(yLength/m_numOfYs, 2) 
                             + pow(zLength/m_numOfZs, 2) );
}

Bucket::~Bucket()
{
    removeBucketTable();
}

rg_Point3D Bucket::getMinPointOfBoundingBox() const
{
	return m_minPointOfBoundingBox;
}

rg_Point3D Bucket::getMaxPointOfBoundingBox() const
{
	return m_maxPointOfBoundingBox;
}

rg_REAL Bucket::getDiagonalLength() const
{
	return m_diagonalLength;
}

void Bucket::setMinPointOfBoundingBox(const rg_Point3D& minPt)
{
    m_minPointOfBoundingBox = minPt;
}

void Bucket::setMaxPointOfBoundingBox(const rg_Point3D& maxPt)
{
    m_maxPointOfBoundingBox = maxPt;
}

//re-construct bucket by the given list of cells
void Bucket::constructBucket( const rg_Point3D&         minPt, 
                              const rg_Point3D&         maxPt, 
                                    rg_dList< BallGenerator >& ballList)
{
    m_minPointOfBoundingBox = minPt;
    m_maxPointOfBoundingBox = maxPt;

	if( m_bucketTable != rg_NULL )
        removeBucketTable();

	//bucket의 대략적인 size : 현재는 generator개수의 두배.
	rg_INT N = ballList.getSize() * 2;
//	rg_INT N = ballList.getSize() * 16;

	//size of bounding box
	rg_REAL xLength = m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX();
	rg_REAL yLength = m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY();
	rg_REAL zLength = m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ();

	m_numOfZs = (rg_INT)( pow( (zLength*zLength)/(xLength*yLength) * N, 1./3 ) + 1. );
	m_numOfYs = (rg_INT)( yLength/zLength * m_numOfZs + 1. );
	m_numOfXs = (rg_INT)( xLength/zLength * m_numOfZs + 1. );

	m_bucketTable = new rg_dList<BallGenerator*>** [m_numOfXs];
    m_bucketMark  = new rg_FLAG** [m_numOfXs];
	for(rg_INT i=0; i<m_numOfXs; i++)  {
		m_bucketTable[i] = new rg_dList<BallGenerator*>* [m_numOfYs];
		m_bucketMark[i]  = new rg_FLAG* [m_numOfYs];

        for(rg_INT j=0; j<m_numOfYs; j++)  {
			m_bucketTable[i][j] = new rg_dList<BallGenerator*> [m_numOfZs];
			m_bucketMark[i][j]  = new rg_FLAG [m_numOfZs];

            for ( rg_INT k=0; k<m_numOfZs; k++)  {
                m_bucketMark[i][j][k] = 0; 
            }
		}
	}

    //insert cell into bucket structure
	ballList.reset4Loop();
	while(ballList.setNext4Loop())
	{
		addCellIntoBucket( ballList.getpEntity() );
	}

	m_diagonalLength = sqrt( pow(xLength/m_numOfXs, 2) + pow(yLength/m_numOfYs, 2) + pow(zLength/m_numOfZs, 2) );
}


void Bucket::setBucket(rg_dList< BallGenerator >& ballList)
{
	if( m_bucketTable != rg_NULL )
        removeBucketTable();

    constructBoundingBox(ballList);
	
	//bucket의 대략적인 size : 현재는 generator개수의 두배.
	rg_INT N = ballList.getSize() * 2;

	//size of bounding box
	rg_REAL xLength = m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX();
	rg_REAL yLength = m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY();
	rg_REAL zLength = m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ();

	m_numOfZs = (rg_INT)( pow( (zLength*zLength)/(xLength*yLength) * N, 1./3 ) + 1. );
	m_numOfYs = (rg_INT)( yLength/zLength * m_numOfZs + 1. );
	m_numOfXs = (rg_INT)( xLength/zLength * m_numOfZs + 1. );


    m_bucketTable = new rg_dList<BallGenerator*>** [m_numOfXs];
    m_bucketMark  = new rg_FLAG** [m_numOfXs];
	for(rg_INT i=0; i<m_numOfXs; i++)  {
		m_bucketTable[i] = new rg_dList<BallGenerator*>* [m_numOfYs];
		m_bucketMark[i]  = new rg_FLAG* [m_numOfYs];

        for(rg_INT j=0; j<m_numOfYs; j++)  {
			m_bucketTable[i][j] = new rg_dList<BallGenerator*> [m_numOfZs];
			m_bucketMark[i][j]  = new rg_FLAG [m_numOfZs];

            for ( rg_INT k=0; k<m_numOfZs; k++)  {
                m_bucketMark[i][j][k] = 0; 
            }
		}
	}

	//insert cell into bucket structure
	ballList.reset4Loop();
	while(ballList.setNext4Loop())
	{
		addCellIntoBucket( ballList.getpEntity() );
	}

	m_diagonalLength = sqrt( pow(xLength/m_numOfXs, 2) + pow(yLength/m_numOfYs, 2) + pow(zLength/m_numOfZs, 2) );
}

void Bucket::constructBoundingBox(rg_dList< BallGenerator >& ballList)
{
    m_minPointOfBoundingBox.setPoint( rg_MAX_REAL,  rg_MAX_REAL,  rg_MAX_REAL);
    m_maxPointOfBoundingBox.setPoint(-rg_MAX_REAL, -rg_MAX_REAL, -rg_MAX_REAL);

    rg_REAL minX =  rg_MAX_REAL;
    rg_REAL minY =  rg_MAX_REAL;
    rg_REAL minZ =  rg_MAX_REAL;

    rg_REAL maxX = -rg_MAX_REAL;
    rg_REAL maxY = -rg_MAX_REAL;
    rg_REAL maxZ = -rg_MAX_REAL;

	rg_REAL x;
	rg_REAL y;
	rg_REAL z;
	rg_REAL radius;

    BallGenerator* currentGenerator = rg_NULL;

	ballList.reset4Loop();
	while( ballList.setNext4Loop() )
	{
        currentGenerator = ballList.getpEntity();

		x      = currentGenerator->getCenter().getX();
		y      = currentGenerator->getCenter().getY();
		z      = currentGenerator->getCenter().getZ();
		radius = currentGenerator->getRadius();

		if( rg_GT( x + radius, maxX ) )
			maxX = x + radius;
		if( rg_GT( y + radius, maxY ) )
			maxY = y + radius;
		if( rg_GT( z + radius, maxZ ) )
			maxZ = z + radius;

		if( rg_LT( x - radius, minX ) )
			minX = x - radius;
		if( rg_LT( y - radius, minY ) )
			minY = y - radius;
		if( rg_LT( z - radius, minZ ) )
			minZ = z - radius;
	}

    m_minPointOfBoundingBox.setPoint( minX, minY, minZ);
    m_maxPointOfBoundingBox.setPoint( maxX, maxY, maxZ);
}

//masking
//return generators (cells) without duplication
rg_dList< BallGenerator* >* Bucket::getCellsByMask(const Sphere& sphere)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D center = sphere.getCenter();
    rg_REAL    radius = sphere.getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	if( maxIndexX > m_numOfXs-1 )
		maxIndexX = m_numOfXs-1;
	if( maxIndexY > m_numOfYs-1 )
		maxIndexY = m_numOfYs-1;
	if( maxIndexZ > m_numOfZs-1 )
		maxIndexZ = m_numOfZs-1;

	if( minIndexX < 0 )
		minIndexX = 0;
	if( minIndexY < 0 )
		minIndexY = 0;
	if( minIndexZ < 0 )
		minIndexZ = 0;



	BallGenerator* currBall = NULL;
	rg_dList< BallGenerator* >* resultList = new rg_dList< BallGenerator* >;

	for(int i = minIndexX; i <=maxIndexX; i++)
	{
		for(int j = minIndexY; j <= maxIndexY; j++)
		{
			for(int k= minIndexZ; k <= maxIndexZ; k++)
			{
				m_bucketTable[i][j][k].reset4Loop();
				while( m_bucketTable[i][j][k].setNext4Loop() )
				{
					currBall = m_bucketTable[i][j][k].getEntity();
                    
					if( !(currBall->m_checkForBucket) )
					{
                        if ( currBall->isThereIntersectionWith(sphere) )
                        {
						    resultList->add(currBall);
						    currBall->m_checkForBucket = rg_TRUE;
                        }
                    }
				}
			}
		}
	}
    
	//resultList->reset4Loop();
	//while( resultList->setNext4Loop() )
	//{
	//	resultList->getEntity()->m_checkForBucket = rg_FALSE;
	//}

    return resultList;
    

    /*
	BallGenerator* currBall = NULL;
	rg_dList< BallGenerator* > cellsInMask;
	for(int i = minIndexX; i <=maxIndexX; i++)
	{
		for(int j = minIndexY; j <= maxIndexY; j++)
		{
			for(int k= minIndexZ; k <= maxIndexZ; k++)
			{
				m_bucketTable[i][j][k].reset4Loop();
				while( m_bucketTable[i][j][k].setNext4Loop() )
				{
					currBall = m_bucketTable[i][j][k].getEntity();
                    
					if( !(currBall->m_checkForBucket) )
					{
					    cellsInMask.add(currBall);

					    currBall->m_checkForBucket = rg_TRUE;
					}
				}
			}
		}
	}

	rg_dList< BallGenerator* >* resultList = new rg_dList< BallGenerator* >;
	cellsInMask.reset4Loop();
	while( cellsInMask.setNext4Loop() )
	{
		currBall = cellsInMask.getEntity();

        if ( currBall->isThereIntersectionWith(sphere) )
        {
			resultList->add(currBall);
        }
        else
        {
			currBall->m_checkForBucket = rg_FALSE;
        }
	}
	return resultList;
    */
}



rg_INT Bucket::getCellsByMask(rg_dList< BallGenerator* >& resultList, const Sphere& sphere)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D center = sphere.getCenter();
    rg_REAL    radius = sphere.getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	if( maxIndexX > m_numOfXs-1 )
		maxIndexX = m_numOfXs-1;
	if( maxIndexY > m_numOfYs-1 )
		maxIndexY = m_numOfYs-1;
	if( maxIndexZ > m_numOfZs-1 )
		maxIndexZ = m_numOfZs-1;

	if( minIndexX < 0 )
		minIndexX = 0;
	if( minIndexY < 0 )
		minIndexY = 0;
	if( minIndexZ < 0 )
		minIndexZ = 0;



	BallGenerator* currBall = NULL;

	for(int i = minIndexX; i <=maxIndexX; i++)
	{
		for(int j = minIndexY; j <= maxIndexY; j++)
		{
			for(int k= minIndexZ; k <= maxIndexZ; k++)
			{
				m_bucketTable[i][j][k].reset4Loop();
				while( m_bucketTable[i][j][k].setNext4Loop() )
				{
					currBall = m_bucketTable[i][j][k].getEntity();
                    
					if( !(currBall->m_checkForBucket) )
					{
                        if ( currBall->isThereIntersectionWith(sphere) )
                        {
						    resultList.add(currBall);
                            currBall->m_checkForBucket = rg_TRUE;
                        }
                    }
				}
			}
		}
	}

    rg_INT xs = maxIndexX - minIndexX + 1;
    rg_INT ys = maxIndexY - minIndexY + 1;
    rg_INT zs = maxIndexZ - minIndexZ + 1;

    return (xs*ys*zs);
}



void Bucket::getBallsIntersectingSphere(const Sphere& sphere, rg_dList< BallGenerator* >& resultList)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D center = sphere.getCenter();
    rg_REAL    radius = sphere.getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	if( maxIndexX > m_numOfXs-1 )
		maxIndexX = m_numOfXs-1;
	if( maxIndexY > m_numOfYs-1 )
		maxIndexY = m_numOfYs-1;
	if( maxIndexZ > m_numOfZs-1 )
		maxIndexZ = m_numOfZs-1;

	if( minIndexX < 0 )
		minIndexX = 0;
	if( minIndexY < 0 )
		minIndexY = 0;
	if( minIndexZ < 0 )
		minIndexZ = 0;



	BallGenerator* currBall = NULL;

	for(int i = minIndexX; i <=maxIndexX; i++)  {
		for(int j = minIndexY; j <= maxIndexY; j++)  {
			for(int k= minIndexZ; k <= maxIndexZ; k++)  {
				m_bucketTable[i][j][k].reset4Loop();
				
                while( m_bucketTable[i][j][k].setNext4Loop() )  {
					currBall = m_bucketTable[i][j][k].getEntity();
                    
					if( !(currBall->m_checkForBucket) )  {
                        if ( currBall->isThereIntersectionWith(sphere) )  {
						    resultList.add(currBall);
                            currBall->m_checkForBucket = rg_TRUE;
                        }
                    }
				}
			}
		}
	}

	resultList.reset4Loop();
	while( resultList.setNext4Loop() )  {
		resultList.getEntity()->m_checkForBucket = rg_FALSE;
	}

}



void Bucket::getBallsBoundedByPlaneAndSphere(rg_dList< BallGenerator* >& resultList, 
                                             rg_REAL*                    plane, 
                                             const Sphere&               sphere)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D center = sphere.getCenter();
    rg_REAL    radius = sphere.getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	if( maxIndexX > m_numOfXs-1 )
		maxIndexX = m_numOfXs-1;
	if( maxIndexY > m_numOfYs-1 )
		maxIndexY = m_numOfYs-1;
	if( maxIndexZ > m_numOfZs-1 )
		maxIndexZ = m_numOfZs-1;

	if( minIndexX < 0 )
		minIndexX = 0;
	if( minIndexY < 0 )
		minIndexY = 0;
	if( minIndexZ < 0 )
		minIndexZ = 0;


    rg_dList< BallGenerator* > ballsList;
	BallGenerator* currBall = NULL;

	for(int i = minIndexX; i <=maxIndexX; i++)  {
		for(int j = minIndexY; j <= maxIndexY; j++)  {
			for(int k= minIndexZ; k <= maxIndexZ; k++)  {
                ballsList.append( m_bucketTable[i][j][k] );
			}
		}
	}

	ballsList.reset4Loop();				
    while( ballsList.setNext4Loop() )  {
		currBall = ballsList.getEntity();
        
		if( !(currBall->m_checkForBucket) )  {
            rg_Point3D center   = currBall->getCenter();
            rg_REAL    distance = center.getX()*plane[0] + center.getY()*plane[1] + center.getZ()*plane[2] + plane[3];
            if ( distance + currBall->getRadius() >= 0.0 )  {

                if ( currBall->isThereIntersectionWith(sphere) )  {
				    resultList.add(currBall);
                    currBall->m_checkForBucket = rg_TRUE;
                }
            }
        }
	}

}



void Bucket::addCellIntoBucket(BallGenerator* aBall)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    {
        rg_Point3D center = aBall->getCenter();
        rg_REAL    radius = aBall->getRadius();

	    rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	    rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	    getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	    getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	    if( maxIndexX > m_numOfXs-1 )
		    maxIndexX = m_numOfXs-1;
	    if( maxIndexY > m_numOfYs-1 )
		    maxIndexY = m_numOfYs-1;
	    if( maxIndexZ > m_numOfZs-1 )
		    maxIndexZ = m_numOfZs-1;

	    if( minIndexX < 0 )
		    minIndexX = 0;
	    if( minIndexY < 0 )
		    minIndexY = 0;
	    if( minIndexZ < 0 )
		    minIndexZ = 0;

	    //해당 bucket에 입력
	    for(int i = minIndexX; i <=maxIndexX; i++)
	    {
		    for(int j = minIndexY; j <= maxIndexY; j++)
		    {
			    for(int k= minIndexZ; k <= maxIndexZ; k++)
			    {
				    m_bucketTable[i][j][k].add(aBall);
			    }
		    }
	    }
    }
}



void Bucket::getBucketIndex(const rg_Point3D& pt, rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ)
{
	indexX = (rg_INT) (   ( pt.getX() - m_minPointOfBoundingBox.getX() )
                        / ( m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX() )
                        * m_numOfXs );
	indexY = (rg_INT) (   ( pt.getY() - m_minPointOfBoundingBox.getY() )
                        / ( m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY() )
                        * m_numOfYs );
	indexZ = (rg_INT) (   ( pt.getZ() - m_minPointOfBoundingBox.getZ() )
                        / ( m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ() )
                        * m_numOfZs );
}



rg_FLAG Bucket::isThereIntersectingBallWithSphere(const Sphere& sphere)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D center = sphere.getCenter();
    rg_REAL    radius = sphere.getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	if( maxIndexX > m_numOfXs-1 )
		maxIndexX = m_numOfXs-1;
	if( maxIndexY > m_numOfYs-1 )
		maxIndexY = m_numOfYs-1;
	if( maxIndexZ > m_numOfZs-1 )
		maxIndexZ = m_numOfZs-1;

	if( minIndexX < 0 )
		minIndexX = 0;
	if( minIndexY < 0 )
		minIndexY = 0;
	if( minIndexZ < 0 )
		minIndexZ = 0;

	BallGenerator* currBall = NULL;
	for(int i = minIndexX; i <=maxIndexX; i++)
	{
		for(int j = minIndexY; j <= maxIndexY; j++)
		{
			for(int k= minIndexZ; k <= maxIndexZ; k++)
			{
                if ( m_bucketTable[i][j][k].getSize() == 0 )
                    continue;

				m_bucketTable[i][j][k].reset4Loop();
				while( m_bucketTable[i][j][k].setNext4Loop() )
				{
					currBall = m_bucketTable[i][j][k].getEntity();

                    if ( currBall->isThisGeneratorIntersectedWith(sphere) )
                        return rg_TRUE;
				}
			}
		}
	}

    return rg_FALSE;
}

/*
//This function is for the test to enhance the performance of the bucket data structure!!
void Bucket::outputStatistics()
{
	ofstream fout("stats.out");
	float aver = 0;
	for(int i=0; i<numOfX; i++)
	{
		for(int j=0; j<numOfY; j++)
		{
			for(int k=0; k<numOfZ; k++)
			{
				fout<<i<<' '<<j<<' '<<k<<":"<<bucket[i][j][k].getSize()<<endl;
				aver+=bucket[i][j][k].getSize();
			}
		}
	}

	fout<<"Average: "<<aver<<endl;
	fout<<"# of X: "<<numOfX<<endl;
	fout<<"# of Y: "<<numOfY<<endl;
	fout<<"# of Z: "<<numOfZ<<endl;

	fout.close();


}
*/

void Bucket::removeBucketTable()
{
	if( m_bucketTable != rg_NULL )
    {
	    for(rg_INT i=0; i<m_numOfXs; i++)
	    {
		    if( m_bucketTable[i] != rg_NULL )
		    {
			    for(rg_INT j=0; j<m_numOfYs; j++)
			    {
				    if( m_bucketTable[i][j] != rg_NULL )
					    delete [] m_bucketTable[i][j];
			    }
			    delete [] m_bucketTable[i];
		    }
	    }
	    delete [] m_bucketTable;

        m_bucketTable = rg_NULL;
    }

	if( m_bucketMark != rg_NULL )
    {
	    for(rg_INT i=0; i<m_numOfXs; i++)
	    {
		    if( m_bucketMark[i] != rg_NULL )
		    {
			    for(rg_INT j=0; j<m_numOfYs; j++)
			    {
				    if( m_bucketMark[i][j] != rg_NULL )
					    delete [] m_bucketMark[i][j];
			    }
			    delete [] m_bucketMark[i];
		    }
	    }
	    delete [] m_bucketMark;

        m_bucketMark = rg_NULL;
    }
}



//  by Youngsong Cho (2005. 7. 22)
void Bucket::constructBucket( const rg_REAL&            loadFactor,
                              const rg_Point3D&         minPt, 
                              const rg_Point3D&         maxPt, 
                              const rg_dList< BallGenerator >& ballList )
{    
    m_minPointOfBoundingBox = minPt;
    m_maxPointOfBoundingBox = maxPt;

	if( m_bucketTable != rg_NULL )
        removeBucketTable();

	//bucket의 대략적인 size : 현재는 generator개수의 두배.
	rg_INT N = (rg_INT) (ballList.getSize() * loadFactor);
    if ( N == 0 )
        N = 1;
//	rg_INT N = ballList.getSize() * 16;

	//size of bounding box
	rg_REAL xLength = m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX();
	rg_REAL yLength = m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY();
	rg_REAL zLength = m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ();

	m_numOfZs = (rg_INT)( pow( (zLength*zLength)/(xLength*yLength) * N, 1./3 ) + 1. );
	m_numOfYs = (rg_INT)( yLength/zLength * m_numOfZs + 1. );
	m_numOfXs = (rg_INT)( xLength/zLength * m_numOfZs + 1. );

	m_bucketTable = new rg_dList<BallGenerator*>** [m_numOfXs];
    m_bucketMark  = new rg_FLAG** [m_numOfXs];
	for(rg_INT i=0; i<m_numOfXs; i++)  {
		m_bucketTable[i] = new rg_dList<BallGenerator*>* [m_numOfYs];
		m_bucketMark[i]  = new rg_FLAG* [m_numOfYs];

        for(rg_INT j=0; j<m_numOfYs; j++)  {
			m_bucketTable[i][j] = new rg_dList<BallGenerator*> [m_numOfZs];
			m_bucketMark[i][j]  = new rg_FLAG [m_numOfZs];

            for ( rg_INT k=0; k<m_numOfZs; k++)  {
                m_bucketMark[i][j][k] = 0; 
            }
		}
	}

    //insert cell into bucket structure
	ballList.reset4Loop();
	while(ballList.setNext4Loop())
	{
		addCellIntoBucket( ballList.getpEntity() );
	}

	m_diagonalLength = sqrt( pow(xLength/m_numOfXs, 2) + pow(yLength/m_numOfYs, 2) + pow(zLength/m_numOfZs, 2) );
}



void Bucket::constructBucketBySizeOfElement( const rg_REAL&            sizeOfCube,
                                             const rg_Point3D&         minPt, 
                                             const rg_Point3D&         maxPt, 
                                             const rg_dList< BallGenerator >& ballList )
{
    m_unitSize = sizeOfCube;

    rg_REAL min[3] ={ minPt.getX(), minPt.getY(), minPt.getZ() };
    rg_REAL max[3] ={ maxPt.getX(), maxPt.getY(), maxPt.getZ() };

	rg_INT i=0;
    for ( i=0; i<3; i++)  {
        rg_REAL size        = max[i] - min[i];
        rg_INT  numElements = (rg_INT)(size/sizeOfCube);
        rg_REAL remnant     = size - (numElements*sizeOfCube);
        rg_REAL increment   = sizeOfCube - remnant;
        rg_REAL variation   = increment/2.;

        if ( remnant != 0.0 )  {
            max[i] += variation;
            min[i] -= variation;
        }
    }

    double numX = (max[0] - min[0])/sizeOfCube;
    double ceil_numX = ceil( numX );
    if ( numX < ceil_numX && rg_ZERO(numX - ceil_numX) ) {
        m_numOfXs = (rg_INT)(ceil(numX));
    }
    else {
        m_numOfXs = (rg_INT)(numX);
    }

    double numY = (max[1] - min[1])/sizeOfCube;
    double ceil_numY = ceil( numY );
    if ( numY < ceil_numY && rg_ZERO(numY - ceil_numY) ) {
        m_numOfYs = (rg_INT)(ceil(numY));
    }
    else {
        m_numOfYs = (rg_INT)(numY);
    }

    double numZ = (max[2] - min[2])/sizeOfCube;
    double ceil_numZ = ceil( numZ );
    if ( numZ < ceil_numZ && rg_ZERO(numZ - ceil_numZ) ) {
        m_numOfZs = (rg_INT)(ceil(numZ));
    }
    else {
        m_numOfZs = (rg_INT)(numZ);
    }

    m_minPointOfBoundingBox.setPoint( min[0], min[1], min[2] );
    m_maxPointOfBoundingBox.setPoint( max[0], max[1], max[2] );



	m_bucketTable = new rg_dList<BallGenerator*>** [m_numOfXs];
    m_bucketMark  = new rg_FLAG** [m_numOfXs];
	for( i=0; i<m_numOfXs; i++)  {
		m_bucketTable[i] = new rg_dList<BallGenerator*>* [m_numOfYs];
		m_bucketMark[i]  = new rg_FLAG* [m_numOfYs];

        for(rg_INT j=0; j<m_numOfYs; j++)  {
			m_bucketTable[i][j] = new rg_dList<BallGenerator*> [m_numOfZs];
			m_bucketMark[i][j]  = new rg_FLAG [m_numOfZs];

            for ( rg_INT k=0; k<m_numOfZs; k++)  {
                m_bucketMark[i][j][k] = 0; 
            }
		}
	}

    //insert cell into bucket structure
    BallGenerator* currBall = rg_NULL;
	ballList.reset4Loop();
	while(ballList.setNext4Loop())
	{
        currBall = ballList.getpEntity();

		addCellIntoBucket( currBall );
	}
}



/*
void Bucket::addBallIntoBucketConsistingOfCubes(BallGenerator* aBall)
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_INT centerIndexX = -1;
    rg_INT centerIndexY = -1;
    rg_INT centerIndexZ = -1;

    rg_Point3D center = aBall->getCenter();
    rg_REAL    radius = aBall->getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);


    rg_INT numDiffIndices = 0;
    if ( maxIndexX - minIndexX == 0 )
        numDiffIndices++;
    if ( maxIndexY - minIndexY == 0 )
        numDiffIndices++;
    if ( maxIndexZ - minIndexZ == 0 )
        numDiffIndices++;

    if ( numDiffIndices == 3 || numDiffIndices == 2 )  {
	    for(int i = minIndexX; i <=maxIndexX; i++)  {
		    for(int j = minIndexY; j <= maxIndexY; j++)  {
			    for(int k= minIndexZ; k <= maxIndexZ; k++)  {
				    m_bucketTable[i][j][k].add(aBall);
			    }
		    }
	    }
    } 
    else if ( numDiffIndices == 1 )  {

    }
    else  { // if ( numDiffIndices == 0 )

    }
}
*/



rg_FLAG Bucket::setMarkAfterCheck(const rg_INT& i, const rg_INT& j, const rg_INT& k, const rg_FLAG& mark)
{ 
	if( i < 0 || i > m_numOfXs-1 )
        return rg_FALSE;
	if( j < 0 || j > m_numOfYs-1 )
        return rg_FALSE;
	if( k < 0 || k > m_numOfZs-1 )
        return rg_FALSE;
  
    if ( m_bucketMark[i][j][k] == rg_TRUE )
        return rg_FALSE;

    m_bucketMark[i][j][k] = mark;

    return rg_TRUE;
}



rg_FLAG Bucket::isVisited(const BucketCellIndex& cellIndex) const
{
	if( cellIndex.m_i < 0 || cellIndex.m_i > m_numOfXs-1 )
        return rg_TRUE;
	if( cellIndex.m_j < 0 || cellIndex.m_j > m_numOfYs-1 )
        return rg_TRUE;
	if( cellIndex.m_k < 0 || cellIndex.m_k > m_numOfZs-1 )
        return rg_TRUE;
  
    return m_bucketMark[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k];
}




rg_FLAG Bucket::setMarkAfterCheck(const BucketCellIndex& cellIndex, const rg_FLAG& mark)  
{ 
	if( cellIndex.m_i < 0 || cellIndex.m_i > m_numOfXs-1 )
        return rg_FALSE;
	if( cellIndex.m_j < 0 || cellIndex.m_j > m_numOfYs-1 )
        return rg_FALSE;
	if( cellIndex.m_k < 0 || cellIndex.m_k > m_numOfZs-1 )
        return rg_FALSE;
  
    if ( m_bucketMark[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k] == rg_TRUE )
        return rg_FALSE;

    m_bucketMark[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k] = mark;

    return rg_TRUE;
}



rg_REAL  Bucket::getXLengthOfCell() const
{
    rg_REAL bucketXLength = m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX();
    rg_REAL xLength       = bucketXLength/m_numOfXs;

    return xLength;
}



rg_REAL  Bucket::getYLengthOfCell() const
{
    rg_REAL bucketYLength = m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY();
    rg_REAL yLength       = bucketYLength/m_numOfYs;

    return yLength;
}



rg_REAL  Bucket::getZLengthOfCell() const
{
    rg_REAL bucketZLength = m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ();
    rg_REAL zLength       = bucketZLength/m_numOfZs;

    return zLength;
}



rg_Point3D  Bucket::getLengthOfCell() const
{
    rg_Point3D length = m_maxPointOfBoundingBox - m_minPointOfBoundingBox;

    rg_REAL xLength   = length.getX()/m_numOfXs;
    rg_REAL yLength   = length.getY()/m_numOfYs;
    rg_REAL zLength   = length.getZ()/m_numOfZs;

    return rg_Point3D(xLength, yLength, zLength);
}




rg_FLAG  Bucket::getPositionOfCellWRTPlane(const BucketCellIndex& cell,
                                           const rg_REAL&         distance,
                                           const rg_Point3D&      cellLength)
{
    if ( cell.m_i < 0 || cell.m_j < 0 || cell.m_k < 0 )   
        return INVALID_CELL;
    if ( cell.m_i >= m_numOfXs || cell.m_j >= m_numOfYs || cell.m_k >= m_numOfZs)   
        return INVALID_CELL;
    if ( m_bucketMark[cell.m_i][cell.m_j][cell.m_k] == rg_TRUE )
        return INVALID_CELL;


    rg_REAL distBetPAndV[8];
    distBetPAndV[0]= distance;
    distBetPAndV[1]= distance + cellLength.getX();
    distBetPAndV[2]= distance + cellLength.getX() + cellLength.getY();
    distBetPAndV[3]= distance                     + cellLength.getY();
    distBetPAndV[4]= distBetPAndV[0] + cellLength.getZ();
    distBetPAndV[5]= distBetPAndV[1] + cellLength.getZ();
    distBetPAndV[6]= distBetPAndV[2] + cellLength.getZ();
    distBetPAndV[7]= distBetPAndV[3] + cellLength.getZ();

    rg_INT cellPosition = 0;
    if ( rg_GE(distBetPAndV[0], 0.0) )
        cellPosition += 1;
    if ( rg_GE(distBetPAndV[1], 0.0) )
        cellPosition += 2;
    if ( rg_GE(distBetPAndV[2], 0.0) )
        cellPosition += 4;
    if ( rg_GE(distBetPAndV[3], 0.0) )
        cellPosition += 8;
    if ( rg_GE(distBetPAndV[4], 0.0) )
        cellPosition += 16;
    if ( rg_GE(distBetPAndV[5], 0.0) )
        cellPosition += 32;
    if ( rg_GE(distBetPAndV[6], 0.0) )
        cellPosition += 64;
    if ( rg_GE(distBetPAndV[7], 0.0) )
        cellPosition += 128;

    if ( cellPosition == 0 )
        return IN_MINUS_PLANE;
    else if ( cellPosition == 255 )
        return IN_PLUS_PLANE;
    else
        return ON_PLANE;
}



rg_INT   Bucket::getNumBalls(const rg_INT& indexX, const rg_INT& indexY, const rg_INT& indexZ)
{
    return m_bucketTable[indexX][indexY][indexZ].getSize();
}



rg_dList< BallGenerator* >* Bucket::getBalls(const BucketCellIndex& cellIndex)
{
    rg_dList< BallGenerator* >* balls = new rg_dList< BallGenerator* >;

    BallGenerator* currBall = rg_NULL;
    m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].reset4Loop();
	while( m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].setNext4Loop() )
	{
	    currBall = m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].getEntity();

		if( currBall->m_checkForBucket )
		    continue;

	    balls->addTail( currBall );
        currBall->m_checkForBucket = rg_TRUE;
    }

    return balls;
}



void Bucket::getBalls( rg_dList< BallGenerator* >& balls, 
                       const BucketCellIndex&      cellIndex )
{
    BallGenerator* currBall = rg_NULL;
    m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].reset4Loop();
	while( m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].setNext4Loop() )
	{
	    currBall = m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].getEntity();

		if( currBall->m_checkForBucket )
		    continue;

        balls.addTail( currBall );
        currBall->m_checkForBucket = rg_TRUE;
    }
}




rg_dList< BallGenerator* >* Bucket::getBalls(const BucketCellIndex& cellIndex, rg_REAL* plane)
{
    rg_dList< BallGenerator* >* balls = new rg_dList< BallGenerator* >;

    BallGenerator* currBall = rg_NULL;
    m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].reset4Loop();
	while( m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].setNext4Loop() )
	{
	    currBall = m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].getEntity();

		if( currBall->m_checkForBucket )
		    continue;

        rg_Point3D center   = currBall->getCenter();
        rg_REAL    distance =   plane[0]*center.getX() + plane[1]*center.getY() 
                              + plane[2]*center.getZ() + plane[3] + currBall->getRadius();

        if ( rg_GE( distance, 0.0 ) )  {
            balls->addTail( currBall );
            currBall->m_checkForBucket = rg_TRUE;
        }
    }

    return balls;

}



void Bucket::getBalls( rg_dList< BallGenerator* >& balls, 
                       const BucketCellIndex& cellIndex, rg_REAL* plane)
{
    BallGenerator* currBall = rg_NULL;
    m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].reset4Loop();
	while( m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].setNext4Loop() )
	{
	    currBall = m_bucketTable[cellIndex.m_i][cellIndex.m_j][cellIndex.m_k].getEntity();

		if( currBall->m_checkForBucket )
		    continue;

        rg_Point3D center   = currBall->getCenter();
        rg_REAL    distance =   plane[0]*center.getX() + plane[1]*center.getY() 
                              + plane[2]*center.getZ() + plane[3] + currBall->getRadius();

        if ( rg_GE( distance, 0.0 ) )  {
            balls.addTail( currBall );
            currBall->m_checkForBucket = rg_TRUE;
        }
    }
}



rg_dList< BallGenerator* >* Bucket::getBalls(rg_REAL* plane)
{
    BucketCellIndex max;
    BucketCellIndex min;

    if ( rg_GE(plane[0], 0.0) && rg_GE(plane[1], 0.0) && rg_GE(plane[2], 0.0) )  {
		max.m_i = m_numOfXs-1;
		max.m_j = m_numOfYs-1;
		max.m_k = m_numOfZs-1;

        rg_REAL minX = -(plane[1]*m_maxPointOfBoundingBox.getY() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL minY = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL minZ = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[1]*m_maxPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex(minX, minY, minZ, min.m_i, min.m_j, min.m_k);
    }
    else if ( rg_GE(plane[0], 0.0) && rg_GE(plane[1], 0.0) && rg_LT(plane[2], 0.0) )  {
		max.m_i = m_numOfXs-1;
		max.m_j = m_numOfYs-1;
		min.m_k = 0;

        rg_REAL minX = -(plane[1]*m_maxPointOfBoundingBox.getY() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL minY = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL maxZ = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[1]*m_maxPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( minX, minY, maxZ, min.m_i, min.m_j, max.m_k );
    }
    else if ( rg_GE(plane[0], 0.0) && rg_LT(plane[1], 0.0) && rg_GE(plane[2], 0.0) )  {
		max.m_i = m_numOfXs-1;
		min.m_j = 0;
		max.m_k = m_numOfZs-1;

        rg_REAL minX = -(plane[1]*m_minPointOfBoundingBox.getY() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL maxY = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL minZ = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[1]*m_minPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( minX, maxY, minZ, min.m_i, max.m_j, min.m_k );
    }
    else if ( rg_GE(plane[0], 0.0) && rg_LT(plane[1], 0.0) && rg_LT(plane[2], 0.0) )  {
		max.m_i = m_numOfXs-1;
		min.m_j = 0;
		min.m_k = 0;

        rg_REAL minX = -(plane[1]*m_minPointOfBoundingBox.getY() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL maxY = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL maxZ = -(plane[0]*m_maxPointOfBoundingBox.getX() + plane[1]*m_minPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( minX, maxY, maxZ, min.m_i, max.m_j, max.m_k );
    }
    else if ( rg_LT(plane[0], 0.0) && rg_GE(plane[1], 0.0) && rg_GE(plane[2], 0.0) )  {
		min.m_i = 0;
		max.m_j = m_numOfYs-1;
		max.m_k = m_numOfZs-1;

        rg_REAL maxX = -(plane[1]*m_maxPointOfBoundingBox.getY() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL minY = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL minZ = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[1]*m_maxPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( maxX, minY, minZ, max.m_i, min.m_j, min.m_k );
    }
    else if ( rg_LT(plane[0], 0.0) && rg_GE(plane[1], 0.0) && rg_LT(plane[2], 0.0) )  {
		min.m_i = 0;
		max.m_j = m_numOfYs-1;
		min.m_k = 0;

        rg_REAL maxX = -(plane[1]*m_maxPointOfBoundingBox.getY() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL minY = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL maxZ = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[1]*m_maxPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( maxX, minY, maxZ, max.m_i, min.m_j, max.m_k );
    }
    else if ( rg_LT(plane[0], 0.0) && rg_LT(plane[1], 0.0) && rg_GE(plane[2], 0.0) )  {
		min.m_i = 0;
		min.m_j = 0;
		max.m_k = m_numOfZs-1;

        rg_REAL maxX = -(plane[1]*m_minPointOfBoundingBox.getY() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL maxY = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[2]*m_maxPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL minZ = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[1]*m_minPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( maxX, maxY, minZ, max.m_i, max.m_j, min.m_k);
    }
    else  {//if ( rg_LT(plane[0], 0.0) && rg_LT(plane[0], 0.0) && rg_LT(plane[0], 0.0) )  {
		min.m_i = 0;
		min.m_j = 0;
		min.m_k = 0;

        rg_REAL maxX = -(plane[1]*m_minPointOfBoundingBox.getY() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[0];
        rg_REAL maxY = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[2]*m_minPointOfBoundingBox.getZ() + plane[3])/plane[1];
        rg_REAL maxZ = -(plane[0]*m_minPointOfBoundingBox.getX() + plane[1]*m_minPointOfBoundingBox.getY() + plane[3])/plane[2];

        getJustifiedBucketIndex( maxX, maxY, maxZ, max.m_i, max.m_j, max.m_k);
    }


    rg_Point3D normalOfGatePlane(plane[0], plane[1], plane[2]);
	BallGenerator* currBall = NULL;
	rg_dList< BallGenerator* >* resultList = new rg_dList< BallGenerator* >;

	for(int i = min.m_i; i <=max.m_i; i++)  {
		for(int j = min.m_j; j <= max.m_j; j++)  {
			for(int k= min.m_k; k <= max.m_k; k++)  {

				m_bucketTable[i][j][k].reset4Loop();
				while( m_bucketTable[i][j][k].setNext4Loop() )
				{
					currBall = m_bucketTable[i][j][k].getEntity();
                    
					if( !(currBall->m_checkForBucket) )  {
                        rg_REAL excludingCondition = currBall->getCenter().innerProduct(normalOfGatePlane) 
                                                     + plane[3] + currBall->getRadius();
                        if ( rg_LT(excludingCondition, 0.0) )
                            continue;

						resultList->add(currBall);
						currBall->m_checkForBucket = rg_TRUE;
                    }
				}
			}
		}
	}

    return resultList;
}

void Bucket::getJustifiedBucketIndex( const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, 
                                      rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ)
{
	indexX = (rg_INT) (   ( x - m_minPointOfBoundingBox.getX() )
                        / ( m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX() )
                        * m_numOfXs );
	indexY = (rg_INT) (   ( y - m_minPointOfBoundingBox.getY() )
                        / ( m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY() )
                        * m_numOfYs );
	indexZ = (rg_INT) (   ( z - m_minPointOfBoundingBox.getZ() )
                        / ( m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ() )
                        * m_numOfZs );

	if( indexX > m_numOfXs-1 )
		indexX = m_numOfXs-1;
	if( indexY > m_numOfYs-1 )
		indexY = m_numOfYs-1;
	if( indexZ > m_numOfZs-1 )
		indexZ = m_numOfZs-1;

	if( indexX < 0 )
		indexX = 0;
	if( indexY < 0 )
		indexY = 0;
	if( indexZ < 0 )
		indexZ = 0;
}



void Bucket::getJustifiedBucketIndex( const rg_Point3D& pt,  
                                      rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ)
{
	indexX = (rg_INT) (   ( pt.getX() - m_minPointOfBoundingBox.getX() )
                        / ( m_maxPointOfBoundingBox.getX() - m_minPointOfBoundingBox.getX() )
                        * m_numOfXs );
	indexY = (rg_INT) (   ( pt.getY() - m_minPointOfBoundingBox.getY() )
                        / ( m_maxPointOfBoundingBox.getY() - m_minPointOfBoundingBox.getY() )
                        * m_numOfYs );
	indexZ = (rg_INT) (   ( pt.getZ() - m_minPointOfBoundingBox.getZ() )
                        / ( m_maxPointOfBoundingBox.getZ() - m_minPointOfBoundingBox.getZ() )
                        * m_numOfZs );

	if( indexX > m_numOfXs-1 )
		indexX = m_numOfXs-1;
	if( indexY > m_numOfYs-1 )
		indexY = m_numOfYs-1;
	if( indexZ > m_numOfZs-1 )
		indexZ = m_numOfZs-1;

	if( indexX < 0 )
		indexX = 0;
	if( indexY < 0 )
		indexY = 0;
	if( indexZ < 0 )
		indexZ = 0;
}



rg_FLAG Bucket::isConstructed() const
{
    if ( m_bucketTable != rg_NULL )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_INT  Bucket::getNumElementInX() const
{
    return m_numOfXs;
}



rg_INT  Bucket::getNumElementInY() const
{
    return m_numOfYs;
}



rg_INT  Bucket::getNumElementInZ() const
{
    return m_numOfZs;
}



rg_INT  Bucket::getNumEmptyElements()
{
    rg_INT numEmptyElements = 0;
    for ( rg_INT i=0; i<m_numOfXs; i++ )  {
        for ( rg_INT j=0; j<m_numOfYs; j++ )  {
            for ( rg_INT k=0; k<m_numOfZs; k++ )  {
                if ( m_bucketTable[i][j][k].getSize() == 0 )
                    numEmptyElements++;
            }
        }
    }

    return numEmptyElements;
}



void Bucket::reportBucket(ofstream& fout)
{
    /*
    rg_INT numVisited = 0;
    BallGenerator* currBall = NULL;
	for(int i = 0; i <m_numOfXs; i++)  {
		for(int j = 0; j <m_numOfYs; j++)  {
            for(int k= 0; k <m_numOfZs; k++)  {
                numVisited += m_bucketMark[i][j][k];
			}
		}
	}

    fout << numVisited << "\t";
    */
    fout << "# cells in X:" << "\t" << m_numOfXs << endl;
    fout << "# cells in Y:" << "\t" << m_numOfYs << endl;
    fout << "# cells in Z:" << "\t" << m_numOfZs << endl;
    fout << endl;

	BallGenerator* currBall = NULL;

	for(int i = 0; i <m_numOfXs; i++)  {
		for(int j = 0; j <m_numOfYs; j++)  {
            for(int k= 0; k <m_numOfZs; k++)  {
            
                fout << "[" << i << "]" << "[" << j << "]" << "[" << k << "]" << "\t";
                fout << m_bucketTable[i][j][k].getSize() << "\t";
                fout << m_bucketMark[i][j][k] << "\t";

				m_bucketTable[i][j][k].reset4Loop();
				while( m_bucketTable[i][j][k].setNext4Loop() )
				{
					currBall = m_bucketTable[i][j][k].getEntity();
                
                    fout << currBall->getCell()->getID() << "\t";
				}

                fout << endl;
			}
		}
	}
}

void Bucket::getBalls(ofstream& fout, rg_REAL* plane)
{
    rg_Point3D normal(plane[0], plane[1], plane[2]);

    rg_REAL minX = m_minPointOfBoundingBox.getX();
    rg_REAL minY = m_minPointOfBoundingBox.getY();
    rg_REAL minZ = m_minPointOfBoundingBox.getZ();
    rg_REAL maxX = m_maxPointOfBoundingBox.getX();
    rg_REAL maxY = m_maxPointOfBoundingBox.getY();
    rg_REAL maxZ = m_maxPointOfBoundingBox.getZ();

    rg_Point3D ptOnBucket[8];
    ptOnBucket[0].setPoint( minX, minY, minZ );
    ptOnBucket[1].setPoint( minX, minY, maxZ );
    ptOnBucket[2].setPoint( minX, maxY, minZ );
    ptOnBucket[3].setPoint( minX, maxY, maxZ );
    ptOnBucket[4].setPoint( maxX, minY, minZ );
    ptOnBucket[5].setPoint( maxX, minY, maxZ );
    ptOnBucket[6].setPoint( maxX, maxY, minZ );
    ptOnBucket[7].setPoint( maxX, maxY, maxZ );

    for (int i=0; i<8; i++)
    {
        rg_REAL distance = normal.innerProduct( ptOnBucket[i] ) + plane[3];
        fout << "pt" << i << "\t" << distance << endl;
    }
}


