#include "PCAAnalyzer.h"
#include "ConstForPCAAnalyzer.h"
#include "rg_Matrix.h"
using namespace V::GeometryTier;


PCAAnalyzer::PCAAnalyzer()
{
    m_dimension     = PCA_DEFAULT_DIMENSION;
    m_numOfRecords  = 0;
    m_inputRecords  = rg_NULL;
    m_outputRecords = rg_NULL;
    m_radiusOfAtoms = rg_NULL;

    m_eigenValues   = rg_NULL;
    m_eigenVectors  = rg_NULL;

}



PCAAnalyzer::PCAAnalyzer( const rg_INT& aDimension )
{
    m_dimension     = aDimension;
    m_numOfRecords  = 0;
    m_inputRecords  = rg_NULL;
    m_outputRecords = rg_NULL;
    m_radiusOfAtoms = rg_NULL;

    m_eigenValues   = rg_NULL;
    m_eigenVectors  = rg_NULL;
}



PCAAnalyzer::~PCAAnalyzer()
{
    for ( rg_INT i_dim=0; i_dim<m_dimension; i_dim++ ) {
        if( m_inputRecords[i_dim] != rg_NULL ) {
            delete [] m_inputRecords[i_dim];
        }
        if( m_outputRecords[i_dim] != rg_NULL ) {
            delete [] m_outputRecords[i_dim];
        }
    }

    if( m_inputRecords != rg_NULL ) {
        delete [] m_inputRecords;
    }
    if( m_outputRecords != rg_NULL ) {
        delete [] m_outputRecords;
    }		
    if( m_radiusOfAtoms != rg_NULL ) {
        delete [] m_radiusOfAtoms;
    }

    if( m_eigenValues != rg_NULL ) {
        delete [] m_eigenValues;
    }
    if( m_eigenVectors != rg_NULL ) {
        delete [] m_eigenVectors;
    }
}



rg_REAL PCAAnalyzer::getLengthOfLocalCoord( const rg_INT& aDimension )
{
    rg_REAL maxLength = 0.0;

    rg_REAL maxValue = -rg_MAX_REAL;
    rg_REAL minValue = rg_MAX_REAL;
    
    rg_INT idOfMaxValue = 0;
    rg_INT idOfMinValue = 0;

    for( rg_INT i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        if( m_outputRecords[aDimension][i_rec] > maxValue ) {
            maxValue = m_outputRecords[aDimension][i_rec];
            idOfMaxValue = i_rec;
        }
        if( m_outputRecords[aDimension][i_rec] < minValue ) {
            minValue = m_outputRecords[aDimension][i_rec];
            idOfMinValue = i_rec;
        }
    }

    maxLength = fabs(maxValue-minValue);
    
    return maxLength;
}



rg_REAL PCAAnalyzer::getLengthOfLocalCoordWithRadiusOfAtom( const rg_INT& aDimension )
{
    rg_REAL maxLength = 0.0;

    rg_REAL maxValue = -rg_MAX_REAL;
    rg_REAL minValue = rg_MAX_REAL;
    
    rg_INT idOfMaxValue = 0;
    rg_INT idOfMinValue = 0;

    for( rg_INT i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        if( m_outputRecords[aDimension][i_rec] > maxValue ) {
            maxValue = m_outputRecords[aDimension][i_rec];
            idOfMaxValue = i_rec;
        }
        if( m_outputRecords[aDimension][i_rec] < minValue ) {
            minValue = m_outputRecords[aDimension][i_rec];
            idOfMinValue = i_rec;
        }
    }

    maxLength = fabs(maxValue-minValue);
    
    if( m_radiusOfAtoms != rg_NULL ) {
        maxLength = maxLength + m_radiusOfAtoms[idOfMaxValue] + m_radiusOfAtoms[idOfMinValue];
    }

    return maxLength;
}



void PCAAnalyzer::getLengthsOfLocalCoordsForXYZ( rg_REAL* arrLengths )
{
    rg_REAL maxValueX = -rg_MAX_REAL;
    rg_REAL minValueX = rg_MAX_REAL;

    rg_REAL maxValueY = -rg_MAX_REAL;
    rg_REAL minValueY = rg_MAX_REAL;

    rg_REAL maxValueZ = -rg_MAX_REAL;
    rg_REAL minValueZ = rg_MAX_REAL;


    for( rg_INT i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        
        // X-COORD
        if( m_outputRecords[0][i_rec] > maxValueX ) {
            maxValueX = m_outputRecords[0][i_rec];
        }
        if( m_outputRecords[0][i_rec] < minValueX ) {
            minValueX = m_outputRecords[0][i_rec];
        }

        // Y-COORD
        if( m_outputRecords[1][i_rec] > maxValueY ) {
            maxValueY = m_outputRecords[1][i_rec];
        }
        if( m_outputRecords[1][i_rec] < minValueY ) {
            minValueY = m_outputRecords[1][i_rec];
        }

        // Z-COORD
        if( m_outputRecords[2][i_rec] > maxValueZ ) {
            maxValueZ = m_outputRecords[2][i_rec];
        }
        if( m_outputRecords[2][i_rec] < minValueZ ) {
            minValueZ = m_outputRecords[2][i_rec];
        }
    }

    arrLengths[0] = fabs(maxValueX-minValueX);
    arrLengths[1] = fabs(maxValueY-minValueY);
    arrLengths[2] = fabs(maxValueZ-minValueZ);
}



void PCAAnalyzer::getLengthsOfLocalCoordsWithRadiusOfAtomForXYZ( rg_REAL* arrLengths )
{
    rg_REAL maxValueX = -rg_MAX_REAL;
    rg_REAL minValueX =  rg_MAX_REAL;

    rg_REAL maxValueY = -rg_MAX_REAL;
    rg_REAL minValueY =  rg_MAX_REAL;

    rg_REAL maxValueZ = -rg_MAX_REAL;
    rg_REAL minValueZ =  rg_MAX_REAL;

    
    rg_INT idOfMaxValueX = 0;
    rg_INT idOfMinValueX = 0;

    rg_INT idOfMaxValueY = 0;
    rg_INT idOfMinValueY = 0;

    rg_INT idOfMaxValueZ = 0;
    rg_INT idOfMinValueZ = 0;


    for( rg_INT i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        
        // X-COORD
        if( m_outputRecords[0][i_rec] > maxValueX ) {
            maxValueX = m_outputRecords[0][i_rec];
            idOfMaxValueX = i_rec;
        }
        if( m_outputRecords[0][i_rec] < minValueX ) {
            minValueX = m_outputRecords[0][i_rec];
            idOfMinValueX = i_rec;
        }

        // Y-COORD
        if( m_outputRecords[1][i_rec] > maxValueY ) {
            maxValueY = m_outputRecords[1][i_rec];
            idOfMaxValueY = i_rec;
        }
        if( m_outputRecords[1][i_rec] < minValueY ) {
            minValueY = m_outputRecords[1][i_rec];
            idOfMinValueY = i_rec;
        }

        // Z-COORD
        if( m_outputRecords[2][i_rec] > maxValueZ ) {
            maxValueZ = m_outputRecords[2][i_rec];
            idOfMaxValueZ = i_rec;
        }
        if( m_outputRecords[2][i_rec] < minValueZ ) {
            minValueZ = m_outputRecords[2][i_rec];
            idOfMinValueZ = i_rec;
        }
    }

    rg_REAL addtionalAtomRadiusForX = 0.0;
    rg_REAL addtionalAtomRadiusForY = 0.0;
    rg_REAL addtionalAtomRadiusForZ = 0.0;

    if( m_radiusOfAtoms != rg_NULL ) {
        addtionalAtomRadiusForX = m_radiusOfAtoms[idOfMaxValueX] + m_radiusOfAtoms[idOfMinValueX];
        addtionalAtomRadiusForY = m_radiusOfAtoms[idOfMaxValueY] + m_radiusOfAtoms[idOfMinValueY];
        addtionalAtomRadiusForZ = m_radiusOfAtoms[idOfMaxValueZ] + m_radiusOfAtoms[idOfMinValueZ];
    }


    arrLengths[0] = fabs(maxValueX-minValueX) + addtionalAtomRadiusForX;
    arrLengths[1] = fabs(maxValueY-minValueY) + addtionalAtomRadiusForY;
    arrLengths[2] = fabs(maxValueZ-minValueZ) + addtionalAtomRadiusForZ;
}



void PCAAnalyzer::getObjectAlignedLocalAxesAndLocalCenter( const rg_FLAG& atomRadius, rg_Point3D* localAxisX, rg_Point3D* localAxisY, rg_Point3D* localAxisZ, rg_Point3D& localCenter )
{
    rg_REAL sumOfAtomCenters[3];
    rg_REAL meanOfAtomCenters[3]; // LOCAL CENTER OF MOLECULE
    
    // Calculate the mean of input records
    rg_INT i_dim = 0;
    rg_INT i_rec = 0;

    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        sumOfAtomCenters[i_dim]  = 0.0;

        for( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
            sumOfAtomCenters[i_dim] += m_inputRecords[i_dim][i_rec];
        }

        meanOfAtomCenters[i_dim] = sumOfAtomCenters[i_dim] / m_numOfRecords;
    }

    rg_REAL lengthsOfLocalCoords[3];

//     if( atomRadius == rg_TRUE ) {
//         getLengthsOfLocalCoordsWithRadiusOfAtomForXYZ( lengthsOfLocalCoords );
//     }
//     else {
        getLengthsOfLocalCoordsForXYZ( lengthsOfLocalCoords );
//    }

    localAxisX[0].setX( meanOfAtomCenters[0] + m_eigenVectors[0]*lengthsOfLocalCoords[0]/2 );
    localAxisX[0].setY( meanOfAtomCenters[1] + m_eigenVectors[1]*lengthsOfLocalCoords[0]/2 );
    localAxisX[0].setZ( meanOfAtomCenters[2] + m_eigenVectors[2]*lengthsOfLocalCoords[0]/2 );

    localAxisX[1].setX( meanOfAtomCenters[0] + m_eigenVectors[0]*(-1)*lengthsOfLocalCoords[0]/2 );
    localAxisX[1].setY( meanOfAtomCenters[1] + m_eigenVectors[1]*(-1)*lengthsOfLocalCoords[0]/2 );
    localAxisX[1].setZ( meanOfAtomCenters[2] + m_eigenVectors[2]*(-1)*lengthsOfLocalCoords[0]/2 );

    localAxisY[0].setX( meanOfAtomCenters[0] + m_eigenVectors[3]*lengthsOfLocalCoords[1]/2 );
    localAxisY[0].setY( meanOfAtomCenters[1] + m_eigenVectors[4]*lengthsOfLocalCoords[1]/2 );
    localAxisY[0].setZ( meanOfAtomCenters[2] + m_eigenVectors[5]*lengthsOfLocalCoords[1]/2 );

    localAxisY[1].setX( meanOfAtomCenters[0] + m_eigenVectors[3]*(-1)*lengthsOfLocalCoords[1]/2 );
    localAxisY[1].setY( meanOfAtomCenters[1] + m_eigenVectors[4]*(-1)*lengthsOfLocalCoords[1]/2 );
    localAxisY[1].setZ( meanOfAtomCenters[2] + m_eigenVectors[5]*(-1)*lengthsOfLocalCoords[1]/2 );

    localAxisZ[0].setX( meanOfAtomCenters[0] + m_eigenVectors[6]*lengthsOfLocalCoords[2]/2 );
    localAxisZ[0].setY( meanOfAtomCenters[1] + m_eigenVectors[7]*lengthsOfLocalCoords[2]/2 );
    localAxisZ[0].setZ( meanOfAtomCenters[2] + m_eigenVectors[8]*lengthsOfLocalCoords[2]/2 );

    localAxisZ[1].setX( meanOfAtomCenters[0] + m_eigenVectors[6]*(-1)*lengthsOfLocalCoords[2]/2 );
    localAxisZ[1].setY( meanOfAtomCenters[1] + m_eigenVectors[7]*(-1)*lengthsOfLocalCoords[2]/2 );
    localAxisZ[1].setZ( meanOfAtomCenters[2] + m_eigenVectors[8]*(-1)*lengthsOfLocalCoords[2]/2 );
  
    localCenter.setPoint( meanOfAtomCenters[0], meanOfAtomCenters[1], meanOfAtomCenters[2] );
}



void PCAAnalyzer::getObjectAlignedBoundingBoxForMolecule( const rg_FLAG& atomRadius, rg_Point3D* boxPoints )
{
    rg_Point3D eigenVector[3];
    eigenVector[0].setPoint( m_eigenVectors[0], m_eigenVectors[1], m_eigenVectors[2] );
    eigenVector[1].setPoint( m_eigenVectors[3], m_eigenVectors[4], m_eigenVectors[5] );
    eigenVector[2].setPoint( m_eigenVectors[6], m_eigenVectors[7], m_eigenVectors[8] );

    
    rg_Point3D localAxisX[2];
    rg_Point3D localAxisY[2];
    rg_Point3D localAxisZ[2]; 
    rg_Point3D localCenter;
    
    getObjectAlignedLocalAxesAndLocalCenter( atomRadius, localAxisX, localAxisY, localAxisZ, localCenter );
    
    // MOVE LOCAL AXES TO GLOBAL ORIGIN
    rg_Point3D translatedLocalAxisX[2];
    rg_Point3D translatedLocalAxisY[2];
    rg_Point3D translatedLocalAxisZ[2];

    translatedLocalAxisX[0] = localAxisX[0] - localCenter;
    translatedLocalAxisX[1] = localAxisX[1] - localCenter;

    translatedLocalAxisY[0] = localAxisY[0] - localCenter;
    translatedLocalAxisY[1] = localAxisY[1] - localCenter;

    translatedLocalAxisZ[0] = localAxisZ[0] - localCenter;
    translatedLocalAxisZ[1] = localAxisZ[1] - localCenter;

    rg_Point3D boxsOfTLocalAxes[8];

    // Connectivity of box points : 0-1, 1-2, 2-3, 3-0 
    //                              | |  | |  | |  | |
    //                              4-5, 5-6, 6-7, 7-8

    boxsOfTLocalAxes[0] = translatedLocalAxisX[0] + translatedLocalAxisY[0] + translatedLocalAxisZ[0];
    boxsOfTLocalAxes[1] = translatedLocalAxisX[0] + translatedLocalAxisY[1] + translatedLocalAxisZ[0];
    boxsOfTLocalAxes[2] = translatedLocalAxisX[0] + translatedLocalAxisY[1] + translatedLocalAxisZ[1];
    boxsOfTLocalAxes[3] = translatedLocalAxisX[0] + translatedLocalAxisY[0] + translatedLocalAxisZ[1];

    boxsOfTLocalAxes[4] = translatedLocalAxisX[1] + translatedLocalAxisY[0] + translatedLocalAxisZ[0];
    boxsOfTLocalAxes[5] = translatedLocalAxisX[1] + translatedLocalAxisY[1] + translatedLocalAxisZ[0];
    boxsOfTLocalAxes[6] = translatedLocalAxisX[1] + translatedLocalAxisY[1] + translatedLocalAxisZ[1];
    boxsOfTLocalAxes[7] = translatedLocalAxisX[1] + translatedLocalAxisY[0] + translatedLocalAxisZ[1];

    boxPoints[0] = boxsOfTLocalAxes[0] + localCenter;
    boxPoints[1] = boxsOfTLocalAxes[1] + localCenter;
    boxPoints[2] = boxsOfTLocalAxes[2] + localCenter;
    boxPoints[3] = boxsOfTLocalAxes[3] + localCenter;
    boxPoints[4] = boxsOfTLocalAxes[4] + localCenter;
    boxPoints[5] = boxsOfTLocalAxes[5] + localCenter;
    boxPoints[6] = boxsOfTLocalAxes[6] + localCenter;
    boxPoints[7] = boxsOfTLocalAxes[7] + localCenter;
}



void PCAAnalyzer::getObjectAlignedBoundingBoxForMoleculeB( const rg_FLAG& atomRadius, rg_Point3D* boxPoints )
{
    
    rg_REAL lengthsOfLocalCoords[3];

    if( atomRadius == rg_TRUE ) {
        getLengthsOfLocalCoordsWithRadiusOfAtomForXYZ( lengthsOfLocalCoords );
    }
    else {
        getLengthsOfLocalCoordsForXYZ( lengthsOfLocalCoords );
    }

    rg_REAL halfLengthX = lengthsOfLocalCoords[0]/2;
    rg_REAL halfLengthY = lengthsOfLocalCoords[1]/2;
    rg_REAL halfLengthZ = lengthsOfLocalCoords[2]/2;


    rg_Point3D boxsOfTLocalAxes[8];

	boxsOfTLocalAxes[0].setPoint(-halfLengthX, halfLengthY, halfLengthZ);
	boxsOfTLocalAxes[1].setPoint( halfLengthX, halfLengthY, halfLengthZ);
	boxsOfTLocalAxes[2].setPoint( halfLengthX,-halfLengthY, halfLengthZ);
	boxsOfTLocalAxes[3].setPoint(-halfLengthX,-halfLengthY, halfLengthZ);
	boxsOfTLocalAxes[4].setPoint(-halfLengthX, halfLengthY,-halfLengthZ);
	boxsOfTLocalAxes[5].setPoint( halfLengthX, halfLengthY,-halfLengthZ);
	boxsOfTLocalAxes[6].setPoint( halfLengthX,-halfLengthY,-halfLengthZ);
	boxsOfTLocalAxes[7].setPoint(-halfLengthX,-halfLengthY,-halfLengthZ);


    rg_Matrix eigenMatrix(m_dimension, m_dimension);

    rg_Matrix outputDataMatrix( m_dimension, 1 );
    rg_Matrix adjustedDataMatrix( m_dimension, 1 );


    rg_INT i_dim = 0;
    rg_INT j_dim = 0;
    rg_INT i_eigenVec = 0;

    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        for ( j_dim=0; j_dim<m_dimension; j_dim++ ) {
            i_eigenVec = j_dim + (i_dim * m_dimension);
            eigenMatrix.setElement( i_dim,  j_dim, m_eigenVectors[i_eigenVec] );
        }
    }

    rg_Matrix tEigenMatrix = eigenMatrix.inverse();


    // Calculate the mean of input records
    rg_REAL meanOfAtomCenterInXYZ[3]; // LOCAL CENTER OF MOLECULE
    rg_REAL sumOfAtomCenters[3];

    rg_INT  i_rec = 0;

    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        sumOfAtomCenters[i_dim]  = 0.0;
        for( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
            sumOfAtomCenters[i_dim] += m_inputRecords[i_dim][i_rec];
        }
        meanOfAtomCenterInXYZ[i_dim] = sumOfAtomCenters[i_dim] / m_numOfRecords;
    }

    rg_Point3D meanAtomCenter( meanOfAtomCenterInXYZ[0], meanOfAtomCenterInXYZ[1], meanOfAtomCenterInXYZ[2] );




    rg_Matrix  matBoxsOfTLocalAxes[8];
    rg_Point3D rotatedPointsOfBox[8];
    for( rg_INT i_pt=0; i_pt<8; i_pt++ ) {
        matBoxsOfTLocalAxes[i_pt].setSize(m_dimension, 1);
        matBoxsOfTLocalAxes[i_pt].setElement( 0, 0, boxsOfTLocalAxes[i_pt].getX() );
        matBoxsOfTLocalAxes[i_pt].setElement( 1, 0, boxsOfTLocalAxes[i_pt].getY() );
        matBoxsOfTLocalAxes[i_pt].setElement( 2, 0, boxsOfTLocalAxes[i_pt].getZ() );

        matBoxsOfTLocalAxes[i_pt] = tEigenMatrix*matBoxsOfTLocalAxes[i_pt];

        rotatedPointsOfBox[i_pt].setPoint( matBoxsOfTLocalAxes[i_pt].getElement(0, 0),
                                           matBoxsOfTLocalAxes[i_pt].getElement(1, 0),
                                           matBoxsOfTLocalAxes[i_pt].getElement(2, 0)  );

        boxPoints[i_pt] = rotatedPointsOfBox[i_pt] + meanAtomCenter;
    }
}



void PCAAnalyzer::setDimension( const rg_INT& aDimension )
{
    m_dimension     = aDimension;    
}



void PCAAnalyzer::setNumberOfRecords( const rg_INT& numOfRecrods )
{
    m_numOfRecords  = numOfRecrods;
}



void PCAAnalyzer::setDimensionAndNumberOfRecordsForMolecule( const rg_INT& aDimension, const rg_INT& numOfRecords )
{
    clear();

    m_dimension    = aDimension;
    m_numOfRecords = numOfRecords;

    m_inputRecords  = new rg_REAL* [m_dimension];
    m_outputRecords = new rg_REAL* [m_dimension];
    m_radiusOfAtoms = new rg_REAL  [m_numOfRecords];

    rg_INT i_dim = 0;
    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        m_inputRecords[i_dim]  = new rg_REAL[m_numOfRecords];
        m_outputRecords[i_dim] = new rg_REAL[m_numOfRecords];
    }

    m_eigenValues  = new rg_REAL [m_dimension];
    m_eigenVectors = new rg_REAL [m_dimension*m_dimension];
}



void PCAAnalyzer::setInputRecords( const rg_INT& aDimension, const rg_INT& numOfRecrods, rg_REAL** inputRecords )
{
    
}



void PCAAnalyzer::setInputRecords( list<rg_Point3D>* inputRecords )
{
    
}



void PCAAnalyzer::setInputRecords( rg_dList<rg_Point3D>* inputRecords )
{
    setDimensionAndNumberOfRecordsForMolecule( PCA_DEFAULT_DIMENSION, inputRecords->getSize() );

    rg_INT i_records = 0;

    inputRecords->reset4Loop();
    while ( inputRecords->setNext4Loop() ) {
        m_inputRecords[0][i_records] = inputRecords->getpEntity()->getX();
        m_inputRecords[1][i_records] = inputRecords->getpEntity()->getY();
        m_inputRecords[2][i_records] = inputRecords->getpEntity()->getZ();
        
        m_radiusOfAtoms[i_records] = 0.0;
        i_records++;
    }    
}



void PCAAnalyzer::setInputRecords( rg_Point3D* pointSet, const rg_INT& numOfRecords )
{
    setDimensionAndNumberOfRecordsForMolecule( PCA_DEFAULT_DIMENSION, numOfRecords );

    rg_INT i_records = 0;

	for( i_records=0; i_records<numOfRecords; i_records++ ) {

		m_inputRecords[0][i_records] = pointSet[i_records].getX();
		m_inputRecords[1][i_records] = pointSet[i_records].getY();
		m_inputRecords[2][i_records] = pointSet[i_records].getZ();

		m_radiusOfAtoms[i_records] = 0.0;
	}
}



// void PCAAnalyzer::setInputRecords( Particle* particleSet, const rg_INT& numOfRecords )
// {
//     setDimensionAndNumberOfRecordsForMolecule( PCA_DEFAULT_DIMENSION, numOfRecords );
// 
//     rg_INT i_records = 0;
// 
// 	for( i_records=0; i_records<numOfRecords; i_records++ ) {
// 
//         Atom* currAtom = (Atom*)(particleSet[i_records].getProperty());
//         rg_Point3D centerOfBall = currAtom->getpAtomBall()->getCenter();
// 
//         m_inputRecords[0][i_records] = centerOfBall.getX();
//         m_inputRecords[1][i_records] = centerOfBall.getY();
//         m_inputRecords[2][i_records] = centerOfBall.getZ();
//         
//         m_radiusOfAtoms[i_records] = currAtom->getpAtomBall()->getRadius();
// 	}    
// }



void PCAAnalyzer::setInputRecords( Molecule* aMolecule )
{
    rg_dList<Atom>* atoms = aMolecule->getAtoms();

    setDimensionAndNumberOfRecordsForMolecule( PCA_DEFAULT_DIMENSION, atoms->getSize() );

    rg_INT i_records = 0;

    atoms->reset4Loop();
    Sphere* atomBall = rg_NULL;
    while ( atoms->setNext4Loop() ) {
        atomBall = atoms->getpEntity()->getpAtomBall();
        rg_Point3D centerOfBall = atomBall->getCenter();

        m_inputRecords[0][i_records] = centerOfBall.getX();
        m_inputRecords[1][i_records] = centerOfBall.getY();
        m_inputRecords[2][i_records] = centerOfBall.getZ();
        
        m_radiusOfAtoms[i_records] = atomBall->getRadius();
        i_records++;
    }
}



void PCAAnalyzer::runPCAForMolecule()
{
    rg_INT numOfEntriesForCovariance = (m_dimension*(m_dimension+1))/2;
    rg_REAL* covariance = new rg_REAL[numOfEntriesForCovariance];  // Assemble covariance matrix as a semi-definite matrix.

    getCovarianceForXYZ( covariance );        // covariance[6] = {xx,xy,yy,xz,yz,zz};
    computeSymmetircEigenvectorsAndEigenvalues( covariance );

    rg_REAL** adjustedInput = getInputRecordsSubtractedByMean();

    setOutputRecords( adjustedInput );
    
    for ( rg_INT i_dim=0; i_dim<m_dimension; i_dim++ ) {
        delete [] adjustedInput[i_dim];
    }

    delete [] adjustedInput;
    delete [] covariance;
}



void PCAAnalyzer::getCovarianceForXYZ( rg_REAL*& covariance )
{
    // Matrix numbering:
    // ____x_y_z_
    // x | 0
    // y | 1 2
    // z | 3 4 5
    // covariance[6] = {xx,xy,yy,xz,yz,zz};

    rg_REAL sumXX = 0.0;
    rg_REAL sumXY = 0.0;
    rg_REAL sumYY = 0.0;
    rg_REAL sumZX = 0.0;
    rg_REAL sumYZ = 0.0;
    rg_REAL sumZZ = 0.0;
    
    rg_REAL** adjustedInput = getInputRecordsSubtractedByMean();
    
    rg_INT    i_dim = 0;
    rg_INT    i_rec = 0; 
    rg_INT    i_ent = 0;

    for ( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        sumXX += ( adjustedInput[0][i_rec]*adjustedInput[0][i_rec] );
        sumXY += ( adjustedInput[0][i_rec]*adjustedInput[1][i_rec] );
        sumYY += ( adjustedInput[1][i_rec]*adjustedInput[1][i_rec] );
        sumZX += ( adjustedInput[2][i_rec]*adjustedInput[0][i_rec] );
        sumYZ += ( adjustedInput[1][i_rec]*adjustedInput[2][i_rec] );
        sumZZ += ( adjustedInput[2][i_rec]*adjustedInput[2][i_rec] );
    }

    covariance[0] = sumXX / m_numOfRecords;
    covariance[1] = sumXY / m_numOfRecords;
    covariance[2] = sumYY / m_numOfRecords;
    covariance[3] = sumZX / m_numOfRecords;
    covariance[4] = sumYZ / m_numOfRecords;
    covariance[5] = sumZZ / m_numOfRecords;
    
    
    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        delete [] adjustedInput[i_dim];
    }

    delete [] adjustedInput;
}



void PCAAnalyzer::getCovariance( rg_REAL*& covariance ) // for flexible dimensions (to be modified...)
{
    rg_REAL x, y, z,
            sumX = 0., sumY = 0., sumZ = 0.,
            sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
            sumXY = 0., sumXZ = 0., sumYZ = 0., 
            xx, yy, zz, xy, xz, yz;

    rg_INT i_rec=0;

    for ( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        x = m_inputRecords[0][i_rec];
        y = m_inputRecords[1][i_rec];
        z = m_inputRecords[1][i_rec];
        sumX += x / m_numOfRecords;
        sumY += y / m_numOfRecords;
        sumZ += z / m_numOfRecords;
        sumX2 += x * x / m_numOfRecords;
        sumY2 += y * y / m_numOfRecords;
        sumZ2 += z * z / m_numOfRecords;
        sumXY += x * y / m_numOfRecords;
        sumXZ += x * z / m_numOfRecords;
        sumYZ += y * z / m_numOfRecords;
    }

    xx = sumX2 - sumX * sumX;
    yy = sumY2 - sumY * sumY;
    zz = sumZ2 - sumZ * sumZ;
    xy = sumXY - sumX * sumY;
    xz = sumXZ - sumX * sumZ;
    yz = sumYZ - sumY * sumZ;

    covariance[0] = xx;
    covariance[1] = xy;
    covariance[2] = yy;
    covariance[3] = xz;
    covariance[4] = yz;
    covariance[5] = zz;
}



void PCAAnalyzer::computeSymmetircEigenvectorsAndEigenvalues( const rg_REAL* covarianceMat, rg_REAL* eigenVectors, rg_REAL* eigenValues )
{
    const rg_INT MAX_ITER = 100;
    static const rg_REAL EPSILON = (rg_REAL)0.00001;

    // number of entries in covarianceMat
    rg_INT nn = (m_dimension*(m_dimension+1))/2;
  
    // copy matrix
    rg_REAL *a = new rg_REAL[nn];
    rg_INT ij;
    for(ij=0; ij<nn; ij++) {
        a[ij] = covarianceMat[ij];
    }
    // Fortran-porting
    a--;
  
    // init diagonalization matrix as the unit matrix
    rg_REAL *v = new rg_REAL[m_dimension*m_dimension];
    ij = 0;
    rg_INT i;
    for(i=0; i<m_dimension; i++)
    for(rg_INT j=0; j<m_dimension; j++) 
      if(i==j)
        v[ij++] = 1.0;
      else {
        v[ij++] = 0.0;
      }
    
    v--;
  
    // construct weight of the non diagonal terms 
    ij = 1;
    rg_REAL a_norm = 0.0;
    for(i=1; i<=m_dimension; i++) {
        for(rg_INT j=1; j<=i; j++) {
            if( i!=j ) {
                rg_REAL a_ij = a[ij];
                a_norm += a_ij * a_ij;
            }
            ij++;
        }
    }
  
    if(a_norm != 0.0) {
        rg_REAL a_normEPS = a_norm * EPSILON;
        rg_REAL thr = a_norm;

        // rotations
        rg_INT nb_iter = 0;
        while( thr > a_normEPS && nb_iter < MAX_ITER ) {
            nb_iter++;
            rg_REAL thr_nn = thr / nn;

            for(rg_INT l=1; l<m_dimension; l++) {
                for(rg_INT m=l+1; m<=m_dimension; m++) {
                  // construct sinx and cosx 
                    rg_INT lq = (l*l-l)/2;
                    rg_INT mq = (m*m-m)/2;

                    rg_INT lm = l + mq;
                    rg_REAL a_lm = a[lm];
                    rg_REAL a_lm_2 = a_lm * a_lm;

                    if(a_lm_2 < thr_nn)
                        continue;

                    rg_INT ll   = l + lq;
                    rg_INT mm   = m + mq;
                    rg_REAL a_ll = a[ll];
                    rg_REAL a_mm = a[mm];

                    rg_REAL delta = a_ll - a_mm;

                    rg_REAL x;
                    if(delta == 0.0) {
                        x = (rg_REAL) -3.14159265 / 4; 
                    }
                    else {
                        x = (rg_REAL)(- atan( (a_lm+a_lm) / delta ) / 2.0);
                    }

                    rg_REAL sinx    = sin(x);
                    rg_REAL cosx    = cos(x);
                    rg_REAL sinx_2  = sinx * sinx;
                    rg_REAL cosx_2  = cosx * cosx;
                    rg_REAL sincos  = sinx * cosx;

                    // rotate L and M columns 
                    rg_INT ilv = m_dimension*(l-1);
                    rg_INT imv = m_dimension*(m-1);

                    rg_INT i;
                    for( i=1; i<=m_dimension;i++ ) {
                        if( (i!=l) && (i!=m) ) {
                            rg_INT iq = (i*i-i)/2;

                            rg_INT im;
                            if( i<m ) {
                              im = i + mq; 
                            }
                            else {
                              im = m + iq;
                            }

                            rg_REAL a_im = a[im];

                            rg_INT il;
                            if( i<l ) {
                              il = i + lq; 
                            }
                            else {
                              il = l + iq;
                            }
                            rg_REAL a_il = a[il];

                            a[il] = a_il * cosx - a_im * sinx;
                            a[im] = a_il * sinx + a_im * cosx;
                        }

                        ilv++;
                        imv++;

                        rg_REAL v_ilv = v[ilv];
                        rg_REAL v_imv = v[imv];

                        v[ilv] = cosx * v_ilv - sinx * v_imv;
                        v[imv] = sinx * v_ilv + cosx * v_imv;
                    } 

                    x = a_lm * sincos; 
                    x += x;

                    a[ll] =  a_ll * cosx_2 + a_mm * sinx_2 - x;
                    a[mm] =  a_ll * sinx_2 + a_mm * cosx_2 + x;
                    a[lm] =  0.0;

                    thr = fabs(thr - a_lm_2);
                }
            }
        }
    }
  
    // convert indices and copy eigen values 
    a++;
    for(i=0; i<m_dimension; i++) {
        rg_INT k = i + (i*(i+1))/2;
        eigenValues[i] = a[k];
    }
    delete [] a;
  
    // sort eigen values and vectors 
    rg_INT *index = new rg_INT[m_dimension];
    for(i=0; i<m_dimension; i++) {
        index[i] = i;
    }
  
    for(i=0; i<(m_dimension-1); i++) {
        rg_REAL x = eigenValues[i];
        rg_INT k = i;

        for(rg_INT j=i+1; j<m_dimension; j++) {
            if( x < eigenValues[j] ) {
                k = j;
                x = eigenValues[j];
            }
        }

        eigenValues[k] = eigenValues[i];
        eigenValues[i] = x;

        rg_INT jj = index[k];
        index[k] = index[i];
        index[i] = jj;
    }


    // save eigen vectors 
    v++;
    ij = 0;
    for( rg_INT k=0; k<m_dimension; k++ )  {
        rg_INT ik = index[k]*m_dimension;
        for( rg_INT i=0; i<m_dimension; i++ ) {
          eigenVectors[ij++] = v[ik++];
        }
    }

    delete [] v;
    delete [] index;
}



void PCAAnalyzer::computeSymmetircEigenvectorsAndEigenvalues( const rg_REAL* covarianceMat )
{
    const rg_INT MAX_ITER = 100;
    static const rg_REAL EPSILON = (rg_REAL)0.00001;

    // number of entries in covarianceMat
    rg_INT nn = (m_dimension*(m_dimension+1))/2;
  
    // copy matrix
    rg_REAL *a = new rg_REAL[nn];
    rg_INT ij;
    for(ij=0; ij<nn; ij++) {
        a[ij] = covarianceMat[ij];
    }
    // Fortran-porting
    a--;
  
    // init diagonalization matrix as the unit matrix
    rg_REAL *v = new rg_REAL[m_dimension*m_dimension];
    ij = 0;
    rg_INT i;
    for(i=0; i<m_dimension; i++)
    for(rg_INT j=0; j<m_dimension; j++) 
      if(i==j)
        v[ij++] = 1.0;
      else {
        v[ij++] = 0.0;
      }
    
    v--;
  
    // construct weight of the non diagonal terms 
    ij = 1;
    rg_REAL a_norm = 0.0;
    for(i=1; i<=m_dimension; i++) {
        for(rg_INT j=1; j<=i; j++) {
            if( i!=j ) {
                rg_REAL a_ij = a[ij];
                a_norm += a_ij * a_ij;
            }
            ij++;
        }
    }
  
    if(a_norm != 0.0) {
        rg_REAL a_normEPS = a_norm * EPSILON;
        rg_REAL thr = a_norm;

        // rotations
        rg_INT nb_iter = 0;
        while( thr > a_normEPS && nb_iter < MAX_ITER ) {
            nb_iter++;
            rg_REAL thr_nn = thr / nn;

            for(rg_INT l=1; l<m_dimension; l++) {
                for(rg_INT m=l+1; m<=m_dimension; m++) {
                  // construct sinx and cosx 
                    rg_INT lq = (l*l-l)/2;
                    rg_INT mq = (m*m-m)/2;

                    rg_INT lm = l + mq;
                    rg_REAL a_lm = a[lm];
                    rg_REAL a_lm_2 = a_lm * a_lm;

                    if(a_lm_2 < thr_nn)
                        continue;

                    rg_INT ll   = l + lq;
                    rg_INT mm   = m + mq;
                    rg_REAL a_ll = a[ll];
                    rg_REAL a_mm = a[mm];

                    rg_REAL delta = a_ll - a_mm;

                    rg_REAL x;
                    if(delta == 0.0) {
                        x = (rg_REAL) -3.14159265 / 4; 
                    }
                    else {
                        x = (rg_REAL)(- atan( (a_lm+a_lm) / delta ) / 2.0);
                    }

                    rg_REAL sinx    = sin(x);
                    rg_REAL cosx    = cos(x);
                    rg_REAL sinx_2  = sinx * sinx;
                    rg_REAL cosx_2  = cosx * cosx;
                    rg_REAL sincos  = sinx * cosx;

                    // rotate L and M columns 
                    rg_INT ilv = m_dimension*(l-1);
                    rg_INT imv = m_dimension*(m-1);

                    rg_INT i;
                    for( i=1; i<=m_dimension;i++ ) {
                        if( (i!=l) && (i!=m) ) {
                            rg_INT iq = (i*i-i)/2;

                            rg_INT im;
                            if( i<m ) {
                              im = i + mq; 
                            }
                            else {
                              im = m + iq;
                            }

                            rg_REAL a_im = a[im];

                            rg_INT il;
                            if( i<l ) {
                              il = i + lq; 
                            }
                            else {
                              il = l + iq;
                            }
                            rg_REAL a_il = a[il];

                            a[il] = a_il * cosx - a_im * sinx;
                            a[im] = a_il * sinx + a_im * cosx;
                        }

                        ilv++;
                        imv++;

                        rg_REAL v_ilv = v[ilv];
                        rg_REAL v_imv = v[imv];

                        v[ilv] = cosx * v_ilv - sinx * v_imv;
                        v[imv] = sinx * v_ilv + cosx * v_imv;
                    } 

                    x = a_lm * sincos; 
                    x += x;

                    a[ll] =  a_ll * cosx_2 + a_mm * sinx_2 - x;
                    a[mm] =  a_ll * sinx_2 + a_mm * cosx_2 + x;
                    a[lm] =  0.0;

                    thr = fabs(thr - a_lm_2);
                }
            }
        }
    }
  
    // convert indices and copy eigen values 
    a++;
    for(i=0; i<m_dimension; i++) {
        rg_INT k = i + (i*(i+1))/2;
        m_eigenValues[i] = a[k];
    }
    delete [] a;
  
    // sort eigen values and vectors 
    rg_INT *index = new rg_INT[m_dimension];
    for(i=0; i<m_dimension; i++) {
        index[i] = i;
    }
  
    for(i=0; i<(m_dimension-1); i++) {
        rg_REAL x = m_eigenValues[i];
        rg_INT k = i;

        for(rg_INT j=i+1; j<m_dimension; j++) {
            if( x < m_eigenValues[j] ) {
                k = j;
                x = m_eigenValues[j];
            }
        }

        m_eigenValues[k] = m_eigenValues[i];
        m_eigenValues[i] = x;

        rg_INT jj = index[k];
        index[k] = index[i];
        index[i] = jj;
    }


    // save eigen vectors 
    v++;
    ij = 0;
    for( rg_INT k=0; k<m_dimension; k++ )  {
        rg_INT ik = index[k]*m_dimension;
        for( rg_INT i=0; i<m_dimension; i++ ) {
          m_eigenVectors[ij++] = v[ik++];
        }
    }

    delete [] v;
    delete [] index;
}



void PCAAnalyzer::clear()
{
    rg_INT i_dim=0;

    if( m_inputRecords != rg_NULL ) {
        for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
            if( m_inputRecords[i_dim] != rg_NULL ) {
                delete [] m_inputRecords[i_dim];
            }
        }
    }

    if( m_outputRecords != rg_NULL ) {
        for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
            if( m_outputRecords[i_dim] != rg_NULL ) {
                delete [] m_outputRecords[i_dim];
            }
        }
    }

    if( m_inputRecords != rg_NULL ) {
        delete [] m_inputRecords;
    }
    if( m_outputRecords != rg_NULL ) {
        delete [] m_outputRecords;
    }		
    if( m_radiusOfAtoms != rg_NULL ) {
        delete [] m_radiusOfAtoms;
    }

    if( m_eigenValues != rg_NULL ) {
        delete [] m_eigenValues;
    }
    if( m_eigenVectors != rg_NULL ) {
        delete [] m_eigenVectors;
    }
}



rg_REAL** PCAAnalyzer::getInputRecordsSubtractedByMean()
{
    rg_REAL** adjustedInput = new rg_REAL* [m_dimension];
    rg_INT i_dim = 0;
    rg_INT i_rec = 0;
    
    rg_REAL* sumOfRecords  = new rg_REAL [m_dimension];
    rg_REAL* meanOfRecords = new rg_REAL [m_dimension];
    
    // Calculate the mean of input records
    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        adjustedInput[i_dim] = new rg_REAL[m_numOfRecords];
        sumOfRecords[i_dim]  = 0.0;

        for( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
            sumOfRecords[i_dim] += m_inputRecords[i_dim][i_rec];
        }

        meanOfRecords[i_dim] = sumOfRecords[i_dim] / m_numOfRecords;
    }

    // set adjusted input records subtracted by mean.
    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        for( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
            adjustedInput[i_dim][i_rec] = m_inputRecords[i_dim][i_rec] - meanOfRecords[i_dim];
        }
    }

    delete [] meanOfRecords;
    delete [] sumOfRecords;

    return adjustedInput;
}



void PCAAnalyzer::setOutputRecords( rg_REAL** adjustedInput )
{
    rg_Matrix eigenMatrix(m_dimension, m_dimension);

    rg_Matrix outputDataMatrix( m_dimension, 1 );
    rg_Matrix adjustedDataMatrix( m_dimension, 1 );


    rg_INT i_dim = 0;
    rg_INT j_dim = 0;
    rg_INT i_eigenVec = 0;

    for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
        for ( j_dim=0; j_dim<m_dimension; j_dim++ ) {
            i_eigenVec = j_dim + (i_dim * m_dimension);
            eigenMatrix.setElement( i_dim,  j_dim, m_eigenVectors[i_eigenVec] );
        }
    }

    rg_INT i_rec = 0;

    for( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
        adjustedDataMatrix.setElement( 0, 0, adjustedInput[0][i_rec] );
        adjustedDataMatrix.setElement( 1, 0, adjustedInput[1][i_rec] );
        adjustedDataMatrix.setElement( 2, 0, adjustedInput[2][i_rec] );

        outputDataMatrix = eigenMatrix*adjustedDataMatrix;

        m_outputRecords[0][i_rec] = outputDataMatrix.getElement(0, 0);
        m_outputRecords[1][i_rec] = outputDataMatrix.getElement(1, 0);
        m_outputRecords[2][i_rec] = outputDataMatrix.getElement(2, 0);
    }


//     rg_INT i_dim = 0;
//     rg_INT i_rec = 0;
//     rg_INT j_dim = 0;
// 
//     rg_INT i_eigenVec = 0;
// 
//     for ( i_dim=0; i_dim<m_dimension; i_dim++ ) {
// 
//         for( i_rec=0; i_rec<m_numOfRecords; i_rec++ ) {
//             m_outputRecords[i_dim][i_rec] = 0.0;
//             
//             for ( j_dim=0; j_dim<m_dimension; j_dim++ ) {
//                 i_eigenVec = j_dim + (i_dim * m_dimension);
//                 m_outputRecords[i_dim][i_rec] += adjustedInput[j_dim][i_rec] * m_eigenVectors[i_eigenVec];
//             }
//         }
//     }
    
}



rg_REAL PCAAnalyzer::getPCValueWithRadiusOfAtom( const rg_INT& id )
{
    rg_REAL maxLengthXCoord = getLengthOfLocalCoordWithRadiusOfAtom( 0 );
    rg_REAL maxLengthYCoord = getLengthOfLocalCoordWithRadiusOfAtom( 1 );
    rg_REAL maxLengthZCoord = getLengthOfLocalCoordWithRadiusOfAtom( 2 );
    
    rg_REAL maxLength = 0.0;
    rg_REAL midLength = 0.0;
    rg_REAL minLength = 0.0;
    
    if( maxLengthXCoord > maxLengthYCoord && maxLengthXCoord > maxLengthZCoord )
        maxLength = maxLengthXCoord;
    else if( maxLengthYCoord > maxLengthXCoord && maxLengthYCoord > maxLengthZCoord )
        maxLength = maxLengthYCoord;
    else if( maxLengthZCoord > maxLengthXCoord && maxLengthZCoord > maxLengthYCoord )
        maxLength = maxLengthZCoord;
    
    if( ( maxLengthXCoord < maxLengthYCoord && maxLengthXCoord > maxLengthZCoord ) ||
        ( maxLengthXCoord > maxLengthYCoord && maxLengthXCoord < maxLengthZCoord )  )
        midLength = maxLengthXCoord;
    else if( ( maxLengthYCoord < maxLengthXCoord && maxLengthYCoord > maxLengthZCoord ) || 
        ( maxLengthYCoord > maxLengthXCoord && maxLengthYCoord < maxLengthZCoord )  )
        midLength = maxLengthYCoord;
    else if( ( maxLengthZCoord < maxLengthXCoord && maxLengthZCoord > maxLengthYCoord ) || 
        ( maxLengthZCoord > maxLengthXCoord && maxLengthZCoord < maxLengthYCoord )  )
        midLength = maxLengthZCoord;
    
    if( maxLengthXCoord < maxLengthYCoord && maxLengthXCoord < maxLengthZCoord )
        minLength = maxLengthXCoord;
    else if( maxLengthYCoord < maxLengthXCoord && maxLengthYCoord < maxLengthZCoord )
        minLength = maxLengthYCoord;
    else if( maxLengthZCoord < maxLengthXCoord && maxLengthZCoord < maxLengthYCoord )
        minLength = maxLengthZCoord;


    rg_REAL pcValueOnDemand = 0.0;

    switch ( id ) {
    case 1:
        pcValueOnDemand = maxLength;
        break;
    case 2:
        pcValueOnDemand = midLength;
        break;
    case 3:
        pcValueOnDemand = minLength;
        break;
    default:
        break;
    }

    return pcValueOnDemand;
}