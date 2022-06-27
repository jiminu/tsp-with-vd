#include "rg_CylindericalSurface3D.h"

// Constructor & destructor

rg_CylindericalSurface3D::rg_CylindericalSurface3D()
{
}

rg_CylindericalSurface3D::~rg_CylindericalSurface3D()
{
}

// Member functions

rg_Point3D rg_CylindericalSurface3D::getlocalOrigin() const
{
	return localOrigin;
}

rg_Point3D rg_CylindericalSurface3D::getlocalXAxis() const
{
	return localXAxis;
}

rg_Point3D rg_CylindericalSurface3D::getlocalYAxis() const
{
	return localYAxis;
}

rg_Point3D rg_CylindericalSurface3D::getlocalZAxis() const
{
	return localZAxis;
}

rg_REAL    rg_CylindericalSurface3D::getRadius() const
{
	return radius;
}

rg_REAL    rg_CylindericalSurface3D::getHeight() const
{
	return height;
}

void rg_CylindericalSurface3D::setLocalOrigin(const rg_Point3D& localOrg)
{
	localOrigin = localOrg;
}

void rg_CylindericalSurface3D::setLocalXAxis(const rg_Point3D& localX)
{
	localXAxis = localX;
}

void rg_CylindericalSurface3D::setLocalYAxis(const rg_Point3D& localY)
{
	localYAxis = localY;
}

void rg_CylindericalSurface3D::setLocalZAxis(const rg_Point3D& localZ)
{
	localZAxis = localZ;
}

void rg_CylindericalSurface3D::setRadius(const rg_REAL& r)
{
	radius = r;
}

void rg_CylindericalSurface3D::setHeight(const rg_REAL& h)
{
	height = h;
}

rg_Point3D** rg_CylindericalSurface3D::triangulateSurface(rg_INT& numOfRow, rg_INT& numOfCol, rg_REAL angleTol)
{
	numOfRow = 2; // v-direction
	numOfCol = rg_INT(360 / angleTol); // u-direction

	rg_Point3D** vertexList = new rg_Point3D* [numOfRow];

	for(rg_INT i = 0;i < numOfRow;i++)
		vertexList[ i ] = new rg_Point3D[numOfCol];

	rg_Point3D localXVector = localXAxis;
	rg_Point3D localYVector = localYAxis;
	rg_Point3D localZVector = localZAxis;

	//rg_REAL angle = 0.0; // parameter u

	rg_REAL PI = 3.1415926535;

	for(rg_INT indexOfRow = 0;indexOfRow < numOfRow;indexOfRow++)
	{
		rg_REAL angle = 0.0; // parameter u
		for(rg_INT indexOfCol = 0;indexOfCol < numOfCol;indexOfCol++)
		{
            rg_Point3D pt=localOrigin;
            rg_REAL cosAngle=cos( angle );
            rg_REAL sinAngle=sin(angle);
            pt+=radius * (cosAngle* localXVector + sinAngle * localYVector); 
            pt+= ((- height / 2 ) + (height * indexOfRow)) * localZVector;
			vertexList[indexOfRow][indexOfCol]= pt;
			angle += angleTol;
		}
	}

	return vertexList;
}


