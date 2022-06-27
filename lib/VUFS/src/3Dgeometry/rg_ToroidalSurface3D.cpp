#include "rg_ToroidalSurface3D.h"

// Constructors and destructors

rg_ToroidalSurface3D::rg_ToroidalSurface3D()
{
	minorRadius = majorRadius = 0.0;
	//rg_Point3D temp(0.0, 0.0, 0.0);
	//localOrigin = localZAxis = temp;
}

rg_ToroidalSurface3D::rg_ToroidalSurface3D(const rg_ToroidalSurface3D & sourceObj)
{
	minorRadius = sourceObj.minorRadius;
	majorRadius = sourceObj.majorRadius;
	localOrigin = sourceObj.localOrigin;
	localXAxis = sourceObj.localXAxis;
	localYAxis = sourceObj.localYAxis;
	localZAxis = sourceObj.localZAxis;
}

rg_ToroidalSurface3D::~rg_ToroidalSurface3D()
{
}

// Member functions

rg_REAL rg_ToroidalSurface3D::getMinorRadius() const
{
	return minorRadius;
}

rg_REAL rg_ToroidalSurface3D::getMajorRadius() const
{
	return majorRadius;
}

rg_Point3D rg_ToroidalSurface3D::getLocalOrigin() const
{
	return localOrigin;
}

rg_Point3D rg_ToroidalSurface3D::getLocalXAxis() const
{
	return localXAxis;
}

rg_Point3D rg_ToroidalSurface3D::getLocalYAxis() const
{
	return localYAxis;
}

rg_Point3D rg_ToroidalSurface3D::getLocalZAxis() const
{
	return localZAxis;
}

void rg_ToroidalSurface3D::setMinorRadius(const rg_REAL& minorR)
{
	minorRadius = minorR;
}

void rg_ToroidalSurface3D::setMajorRadius(const rg_REAL& majorR)
{
	majorRadius = majorR;
}

void rg_ToroidalSurface3D::setLocalOrigin(const rg_Point3D& localOrg)
{
	localOrigin = localOrg;
}

void rg_ToroidalSurface3D::setLocalXAxis(const rg_Point3D& localX)
{
	localXAxis = localX;
}

void rg_ToroidalSurface3D::setLocalYAxis(const rg_Point3D& localY)
{
	localYAxis = localY;
}

void rg_ToroidalSurface3D::setLocalZAxis(const rg_Point3D& localZ)
{
	localZAxis = localZ;
}

rg_Point3D** rg_ToroidalSurface3D::triangulateSurface(rg_INT& numOfRow, rg_INT& numOfCol, rg_REAL majorTol, rg_REAL minorTol)
{
	numOfRow = rg_INT(360 / minorTol); // v-direction
	numOfCol = rg_INT(360 / majorTol); // u-direction

	rg_Point3D** vertexList = new rg_Point3D* [numOfRow];

	for(rg_INT i = 0;i < numOfRow;i++)
		vertexList[ i ] = new rg_Point3D[numOfCol];

	rg_Point3D localXVector = localXAxis;
	rg_Point3D localYVector = localYAxis;
	rg_Point3D localZVector = localZAxis;

	rg_REAL majorAngle = 0.0; // parameter u
	rg_REAL minorAngle = 0.0; // parameter v

	rg_REAL PI = 3.1415926535;

	for(rg_INT indexOfRow = 0;indexOfRow < numOfRow;indexOfRow++)
	{
		for(rg_INT indexOfCol = 0;indexOfCol < numOfCol;indexOfCol++)
		{
            const rg_REAL cosMinorAngle=cos( minorAngle );
            const rg_REAL cosMajorAngle=cos(majorAngle );
            const rg_REAL sinMajorAngle=sin(majorAngle );
            const rg_REAL sinMinorAngle=sin(minorAngle );
			vertexList[indexOfRow][indexOfCol]
			= localOrigin + (majorRadius + minorRadius * cosMinorAngle) * ( cosMajorAngle* localXVector + sinMajorAngle * localYVector) + minorRadius * sinMinorAngle * localZVector;

			majorAngle += majorTol;
		}
		minorAngle += minorTol;
	}

	return vertexList;
}

rg_ToroidalSurface3D& rg_ToroidalSurface3D::operator=(const rg_ToroidalSurface3D& sourceObj)
{
	minorRadius = sourceObj.minorRadius;
	majorRadius = sourceObj.majorRadius;
	localOrigin = sourceObj.localOrigin;
	localXAxis = sourceObj.localXAxis;
	localYAxis = sourceObj.localYAxis;
	localZAxis = sourceObj.localZAxis;

	return (*this);
}


