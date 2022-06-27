//////////////////////////////////////////////////////////////////////////////////////////////////
//
//   PCAAnalyzer.h
//
//
// getObjectAlignedBoundingBoxForMolecule( const rg_FLAG& atomRadius, rg_Point3D* boxPoints );
// 
// EX )  rg_Point3D boxPoints[8];
//       getObjectAlignedBoundingBoxForMolecule( 1, boxPoints ) ;
//
// Connectivity of box points : 0-1, 1-2, 2-3, 3-0 
//                              | |  | |  | |  | |
//                              4-5, 5-6, 6-7, 7-4
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _PCAANALYZER_H
#define _PCAANALYZER_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_dList.h"
//#include "Particle.h"
#include "rg_Molecule.h"
#include "ConstForMolecule.h"

#include <list>
#include <iostream>
#include <math.h>
using namespace std;



namespace V {
namespace GeometryTier {



class PCAAnalyzer
{
private:
    rg_INT       m_dimension;
    rg_INT       m_numOfRecords;
    rg_REAL**    m_inputRecords;
    rg_REAL**    m_outputRecords;

    rg_REAL*     m_radiusOfAtoms;

    rg_REAL*     m_eigenValues;
    rg_REAL*     m_eigenVectors;


public:
    //  CONSTRUCTOR & DECONSTRUCTOR
    PCAAnalyzer();
    PCAAnalyzer( const rg_INT& aDimension );
    ~PCAAnalyzer();

    //  GET FUNCTION
    rg_REAL getLengthOfLocalCoord( const rg_INT& aDimension );
    rg_REAL getLengthOfLocalCoordWithRadiusOfAtom( const rg_INT& aDimension );

    void    getLengthsOfLocalCoordsForXYZ( rg_REAL* arrLengths );
    void    getLengthsOfLocalCoordsWithRadiusOfAtomForXYZ( rg_REAL* arrLengths );

    void    getObjectAlignedLocalAxesAndLocalCenter( const rg_FLAG& atomRadius, rg_Point3D* localAxisX, rg_Point3D* localAxisY, rg_Point3D* localAxisZ, rg_Point3D& localCenter );
    void    getObjectAlignedBoundingBoxForMolecule( const rg_FLAG& atomRadius, rg_Point3D* boxPoints );
    void    getObjectAlignedBoundingBoxForMoleculeB( const rg_FLAG& atomRadius, rg_Point3D* boxPoints );

    rg_REAL getPCValueWithRadiusOfAtom( const rg_INT& id );

    //  SET FUNCTION
    void    setDimension( const rg_INT& aDimension );
    void    setNumberOfRecords( const rg_INT& numOfRecrods );
    void    setInputRecords( const rg_INT& aDimension, const rg_INT& numOfRecrods, rg_REAL** inputRecords );
    void    setInputRecords( list<rg_Point3D>* inputRecords );
	void    setInputRecords( rg_dList<rg_Point3D>* inputRecords );
    void    setInputRecords( rg_Point3D* pointSet, const rg_INT& numOfRecords );
    //void    setInputRecords( Particle*   particleSet, const rg_INT& numOfRecords );
    void    setInputRecords( Molecule* aMolecule );


    //  COMPUTATION FUNCTION
    void    runPCAForMolecule();

    void    getCovarianceForXYZ( rg_REAL*& covariance );
    void    getCovariance( rg_REAL*& covariance );
	void    computeSymmetircEigenvectorsAndEigenvalues( const rg_REAL* covarianceMat, rg_REAL* eigenVectors, rg_REAL* eigenValues );

    void    clear();

private:

    void      setDimensionAndNumberOfRecordsForMolecule( const rg_INT& aDimension, const rg_INT& numOfRecords );
    
    void      computeSymmetircEigenvectorsAndEigenvalues( const rg_REAL* covarianceMat );
    rg_REAL** getInputRecordsSubtractedByMean();

	void      setOutputRecords( rg_REAL** adjustedInput );
    //void      setOutputRecords( const rg_REAL* eigenVectors, rg_REAL** adjustedInput );

};



} // namespace GeometryTier
} // namespace V


#endif

