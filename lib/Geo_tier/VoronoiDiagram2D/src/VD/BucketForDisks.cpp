#include "BucketForDisks.h"
using namespace V::GeometryTier;



BucketForDisks::BucketForDisks(void)
: m_numBucketCellOfXDir(0),
  m_numBucketCellOfYDir(0),
  m_xDirLengthOfBucketCell(0.0),
  m_yDirLengthOfBucketCell(0.0),
  m_bucketCell(NULL)
{
}



BucketForDisks::BucketForDisks( const list<Generator2D*>& diskSet )
{
    setBucketCell( diskSet );
}



BucketForDisks::BucketForDisks(const BucketForDisks& BFD)
{
	copyBucketFrom(BFD);
}



BucketForDisks::~BucketForDisks(void)
{
    destroyBucket();
}



BucketForDisks& BucketForDisks::operator=(const BucketForDisks& BFD)
{
	if (this == &BFD){
		return *this;
	}

	destroyBucket();
	copyBucketFrom(BFD);

	return *this;
}



void BucketForDisks::setBucketCell( const list<Generator2D*>& diskSet )
{
    setBoundingBox( diskSet );

    int     numOfDisks = diskSet.size();

    double xDirLengthOfBoundingBox = m_boundingBox.getMaxPt().getX() - m_boundingBox.getMinPt().getX();
    double yDirLengthOfBoundingBox = m_boundingBox.getMaxPt().getY() - m_boundingBox.getMinPt().getY();

    m_numBucketCellOfXDir = (int)( sqrt((double)numOfDisks * xDirLengthOfBoundingBox / yDirLengthOfBoundingBox) ) + 1;
    m_numBucketCellOfYDir = (int)( (double)m_numBucketCellOfXDir * yDirLengthOfBoundingBox / xDirLengthOfBoundingBox ) + 1;

    m_xDirLengthOfBucketCell = xDirLengthOfBoundingBox / (double)m_numBucketCellOfXDir;
    m_yDirLengthOfBucketCell = yDirLengthOfBoundingBox / (double)m_numBucketCellOfYDir;

    if( m_xDirLengthOfBucketCell > m_yDirLengthOfBucketCell ) {
        m_yDirLengthOfBucketCell = m_xDirLengthOfBucketCell;
    }
    else {
        m_xDirLengthOfBucketCell = m_yDirLengthOfBucketCell;
    }

    m_bucketCell = new BucketCellForDisk*[m_numBucketCellOfXDir];
    for( int i = 0 ; i < m_numBucketCellOfXDir ; i++ ) {
        m_bucketCell[i] = new BucketCellForDisk[m_numBucketCellOfYDir]; 
    }
}



void BucketForDisks::addGenerator( Generator2D* generator )
{
    double radiusOfGenerator    = generator->getDisk().getRadius();
    rg_Point2D minPtOfGenerator = generator->getDisk().getCenterPt() - rg_Point2D( radiusOfGenerator, radiusOfGenerator );
    rg_Point2D maxPtOfGenerator = generator->getDisk().getCenterPt() + rg_Point2D( radiusOfGenerator, radiusOfGenerator );

    int minXIndex = 0;
    int minYIndex = 0;
    int maxXIndex = 0;
    int maxYIndex = 0;

    getBucketCellIndexForPoint( minPtOfGenerator, minXIndex, minYIndex );
    getBucketCellIndexForPoint( maxPtOfGenerator, maxXIndex, maxYIndex );

    for( int i = minXIndex ; i <= maxXIndex ; i++ ) {
        for( int j = minYIndex ; j <= maxYIndex ; j++ ) {
            m_bucketCell[i][j].setContainingGenerator( generator );
        }
    }


}


/*
Generator2D* BucketForDisks::findGoodCloseGenerator( Generator2D* generator )
{
    rg_Point2D centerPtOfGenerator = generator->getDisk()->getCenterPt();

    int xIndexOfGenerator = 0;
    int yIndexOfGenerator = 0;

    //Generator2D* outputGenerator = NULL;

    getBucketCellIndexForPoint( centerPtOfGenerator, xIndexOfGenerator, yIndexOfGenerator );

    if( m_bucketCell[xIndexOfGenerator][yIndexOfGenerator].getContainingGenerator() != NULL ) {
        //outputGenerator = m_bucketCell[xIndexOfGenerator][yIndexOfGenerator].getContainingGenerator();
        return m_bucketCell[xIndexOfGenerator][yIndexOfGenerator].getContainingGenerator();
    }
    else {
        int i_shell = 1;

        bool isMinXIndexFeasible = true;
        bool isMinYIndexFeasible = true;
        bool isMaxXIndexFeasible = true;
        bool isMaxYIndexFeasible = true;

        int minXIndex = 0;
        int minYIndex = 0;
        int maxXIndex = 0;
        int maxYIndex = 0;


        while( (xIndexOfGenerator - i_shell >= 0) || (xIndexOfGenerator + i_shell < m_numBucketCellOfXDir) || (yIndexOfGenerator - i_shell >= 0) || (yIndexOfGenerator + i_shell < m_numBucketCellOfYDir)) {

            // MIN X INDEX
            if( xIndexOfGenerator - i_shell >= 0 ) {
                isMinXIndexFeasible = true;
                minXIndex           = xIndexOfGenerator - i_shell;
            }
            else {
                isMinXIndexFeasible = false;
                minXIndex           = 0;
            }

            // MIN Y INDEX
            if( yIndexOfGenerator - i_shell >= 0 ) {
                isMinYIndexFeasible = true;
                minYIndex           = yIndexOfGenerator - i_shell;
            }
            else {
                isMinYIndexFeasible = false;
                minYIndex           = 0;
            }

            // MAX X INDEX
            if( xIndexOfGenerator + i_shell < m_numBucketCellOfXDir ) {
                isMaxXIndexFeasible = true;
                maxXIndex           = xIndexOfGenerator + i_shell;
            }
            else {
                isMaxXIndexFeasible = false;
                maxXIndex           = m_numBucketCellOfXDir - 1;
            }

            // MAX Y INDEX
            if( yIndexOfGenerator + i_shell < m_numBucketCellOfYDir) {
                isMaxYIndexFeasible = true;
                maxYIndex           = yIndexOfGenerator + i_shell;
            }
            else {
                isMaxYIndexFeasible = false;
                maxYIndex           = m_numBucketCellOfYDir - 1;
            }

            //  SEARCH WAY
            //  - - - - |
            //  |       |
            //  |       |
            //  |       |
            //  | - - - -

            // LEFT LINE
            if( isMinXIndexFeasible ) {
                for( int y_index = minYIndex ; y_index < maxYIndex ; y_index++ ) {
                    if( m_bucketCell[minXIndex][y_index].getContainingGenerator() != NULL ) {
                        return m_bucketCell[minXIndex][y_index].getContainingGenerator();
                    }
                }

                if( !isMaxYIndexFeasible ) {
                    if( m_bucketCell[minXIndex][maxYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[minXIndex][maxYIndex].getContainingGenerator();
                    }
                }
            }

            // BOTTOM LINE
            if( isMinYIndexFeasible ) {
                for( int x_index = minXIndex + 1 ; x_index <= maxXIndex ; x_index++ ) {
                    if( m_bucketCell[x_index][minYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[x_index][minYIndex].getContainingGenerator();
                    }
                }

                if( !isMinXIndexFeasible ) {
                    if( m_bucketCell[minXIndex][minYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[minXIndex][minYIndex].getContainingGenerator();
                    }
                }
            }

            // RIGHT LINE
            if( isMaxXIndexFeasible ) {
                for( int y_index = minYIndex + 1 ; y_index <= maxYIndex ; y_index++ ) {
                    if( m_bucketCell[maxXIndex][y_index].getContainingGenerator() != NULL ) {
                        return m_bucketCell[maxXIndex][y_index].getContainingGenerator();
                    }
                }

                if( !isMinYIndexFeasible ) {
                    if( m_bucketCell[maxXIndex][minYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[maxXIndex][minYIndex].getContainingGenerator();
                    }
                }
            }

            // TOP LINE
            if( isMaxYIndexFeasible ) {
                for( int x_index = minXIndex ; x_index < maxXIndex ; x_index++ ) {
                    if( m_bucketCell[x_index][maxYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[x_index][maxYIndex].getContainingGenerator();
                    }
                }

                if( !isMaxXIndexFeasible ) {
                    if( m_bucketCell[maxXIndex][maxYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[maxXIndex][maxYIndex].getContainingGenerator();
                    }
                }
            }

        }
    }
    

    return NULL;
}
*/


Generator2D* BucketForDisks::findGoodCloseGenerator( const rg_Point2D& pt ) const
{
    if (m_bucketCell == nullptr) {
        return nullptr;
    }

    int xIndexOfGenerator = 0;
    int yIndexOfGenerator = 0;
    
    getBucketCellIndexForPoint( pt, xIndexOfGenerator, yIndexOfGenerator );

    if( m_bucketCell[xIndexOfGenerator][yIndexOfGenerator].getContainingGenerator() != NULL ) {
        //outputGenerator = m_bucketCell[xIndexOfGenerator][yIndexOfGenerator].getContainingGenerator();
        return m_bucketCell[xIndexOfGenerator][yIndexOfGenerator].getContainingGenerator();
    }
    else {
        int i_shell = 1;

        bool isMinXIndexFeasible = true;
        bool isMinYIndexFeasible = true;
        bool isMaxXIndexFeasible = true;
        bool isMaxYIndexFeasible = true;

        int minXIndex = 0;
        int minYIndex = 0;
        int maxXIndex = 0;
        int maxYIndex = 0;


        while( (xIndexOfGenerator - i_shell >= 0) || (xIndexOfGenerator + i_shell < m_numBucketCellOfXDir) || (yIndexOfGenerator - i_shell >= 0) || (yIndexOfGenerator + i_shell < m_numBucketCellOfYDir)) {

            // MIN X INDEX
            if( xIndexOfGenerator - i_shell >= 0 ) {
                isMinXIndexFeasible = true;
                minXIndex           = xIndexOfGenerator - i_shell;
            }
            else {
                isMinXIndexFeasible = false;
                minXIndex           = 0;
            }

            // MIN Y INDEX
            if( yIndexOfGenerator - i_shell >= 0 ) {
                isMinYIndexFeasible = true;
                minYIndex           = yIndexOfGenerator - i_shell;
            }
            else {
                isMinYIndexFeasible = false;
                minYIndex           = 0;
            }

            // MAX X INDEX
            if( xIndexOfGenerator + i_shell < m_numBucketCellOfXDir ) {
                isMaxXIndexFeasible = true;
                maxXIndex           = xIndexOfGenerator + i_shell;
            }
            else {
                isMaxXIndexFeasible = false;
                maxXIndex           = m_numBucketCellOfXDir - 1;
            }

            // MAX Y INDEX
            if( yIndexOfGenerator + i_shell < m_numBucketCellOfYDir) {
                isMaxYIndexFeasible = true;
                maxYIndex           = yIndexOfGenerator + i_shell;
            }
            else {
                isMaxYIndexFeasible = false;
                maxYIndex           = m_numBucketCellOfYDir - 1;
            }

            //  SEARCH WAY
            //  - - - - |
            //  |       |
            //  |       |
            //  |       |
            //  | - - - -

            // LEFT LINE
            if( isMinXIndexFeasible ) {
                for( int y_index = minYIndex ; y_index < maxYIndex ; y_index++ ) {
                    if( m_bucketCell[minXIndex][y_index].getContainingGenerator() != NULL ) {
                        return m_bucketCell[minXIndex][y_index].getContainingGenerator();
                    }
                }

                if( !isMaxYIndexFeasible ) {
                    if( m_bucketCell[minXIndex][maxYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[minXIndex][maxYIndex].getContainingGenerator();
                    }
                }
            }

            // BOTTOM LINE
            if( isMinYIndexFeasible ) {
                for( int x_index = minXIndex + 1 ; x_index <= maxXIndex ; x_index++ ) {
                    if( m_bucketCell[x_index][minYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[x_index][minYIndex].getContainingGenerator();
                    }
                }

                if( !isMinXIndexFeasible ) {
                    if( m_bucketCell[minXIndex][minYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[minXIndex][minYIndex].getContainingGenerator();
                    }
                }
            }

            // RIGHT LINE
            if( isMaxXIndexFeasible ) {
                for( int y_index = minYIndex + 1 ; y_index <= maxYIndex ; y_index++ ) {
                    if( m_bucketCell[maxXIndex][y_index].getContainingGenerator() != NULL ) {
                        return m_bucketCell[maxXIndex][y_index].getContainingGenerator();
                    }
                }

                if( !isMinYIndexFeasible ) {
                    if( m_bucketCell[maxXIndex][minYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[maxXIndex][minYIndex].getContainingGenerator();
                    }
                }
            }

            // TOP LINE
            if( isMaxYIndexFeasible ) {
                for( int x_index = minXIndex ; x_index < maxXIndex ; x_index++ ) {
                    if( m_bucketCell[x_index][maxYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[x_index][maxYIndex].getContainingGenerator();
                    }
                }

                if( !isMaxXIndexFeasible ) {
                    if( m_bucketCell[maxXIndex][maxYIndex].getContainingGenerator() != NULL ) {
                        return m_bucketCell[maxXIndex][maxYIndex].getContainingGenerator();
                    }
                }
            }
            i_shell++;


        }
    }


    return NULL;
}



void BucketForDisks::setBoundingBox( const list<Generator2D*>& diskSet )
{
    for( list<Generator2D*>::const_iterator it_generator = diskSet.begin() ; it_generator != diskSet.end() ; ++it_generator ) {
        rg_Circle2D currDisk = (*it_generator)->getDisk();
        m_boundingBox.updateBoxByAddingCircle( currDisk );
    }
}



void BucketForDisks::getBucketCellIndexForPoint( const rg_Point2D& pt, int& xIndex, int& yIndex ) const
{
    xIndex = (rg_INT) ( (pt.getX() - m_boundingBox.getMinPt().getX()) / m_xDirLengthOfBucketCell );
    yIndex = (rg_INT) ( (pt.getY() - m_boundingBox.getMinPt().getY()) / m_yDirLengthOfBucketCell );

    if(xIndex < 0)
        xIndex = 0;
    if(xIndex > m_numBucketCellOfXDir - 1)
        xIndex = m_numBucketCellOfXDir - 1;

    if(yIndex < 0)
        yIndex = 0;
    if(yIndex > m_numBucketCellOfYDir - 1)
        yIndex = m_numBucketCellOfYDir - 1;
}



void BucketForDisks::destroyBucket()
{
    for( int i = 0 ; i < m_numBucketCellOfXDir ; i++ ) {
        if( m_bucketCell[i] != NULL ) {
            delete [] m_bucketCell[i];
        }
    }

    if( m_bucketCell != NULL ) {
        delete [] m_bucketCell;
    }

	m_bucketCell = NULL;

	m_numBucketCellOfXDir = 0;
	m_numBucketCellOfYDir = 0;
	m_xDirLengthOfBucketCell = 0.0;
	m_yDirLengthOfBucketCell = 0.0;
}



void BucketForDisks::copyBucketFrom(const BucketForDisks& BFD)
{
	m_boundingBox = BFD.m_boundingBox;
	m_numBucketCellOfXDir = BFD.m_numBucketCellOfXDir;
	m_numBucketCellOfYDir = BFD.m_numBucketCellOfYDir;
	m_xDirLengthOfBucketCell = BFD.m_xDirLengthOfBucketCell;
	m_yDirLengthOfBucketCell = BFD.m_yDirLengthOfBucketCell;

	m_bucketCell = new BucketCellForDisk*[m_numBucketCellOfXDir];

	for (int i = 0; i < m_numBucketCellOfXDir; i++) {
		m_bucketCell[i] = new BucketCellForDisk[m_numBucketCellOfYDir];

		for (int j = 0; j < m_numBucketCellOfYDir; j++){
			m_bucketCell[i][j] = BucketCellForDisk(BFD.m_bucketCell[i][j]);
		}
	}
}


