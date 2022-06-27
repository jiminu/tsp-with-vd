#include "rg_BoundingBox3D.h"
#include "rg_RelativeOp.h"

rg_BoundingBox3D::rg_BoundingBox3D()
: minPt(rg_MAX_REAL, rg_MAX_REAL, rg_MAX_REAL),
  maxPt(-rg_MAX_REAL, -rg_MAX_REAL, -rg_MAX_REAL)
{
}

rg_BoundingBox3D::rg_BoundingBox3D(const rg_Point3D& tMinPt, const rg_Point3D& tMaxPt)
: minPt(tMinPt), maxPt(tMaxPt)
{}

rg_BoundingBox3D::rg_BoundingBox3D(const rg_BoundingBox3D& temp)
: minPt(temp.minPt), 
  maxPt(temp.maxPt)
{}

rg_BoundingBox3D::rg_BoundingBox3D( rg_dList<Sphere>& spheres )
{
	constructBoxByAddingSpheres( spheres );
}

rg_BoundingBox3D::~rg_BoundingBox3D()
{}


//  get & set functions.
rg_Point3D rg_BoundingBox3D::getMinPt() const
{
	return minPt;
}

rg_Point3D rg_BoundingBox3D::getMaxPt() const
{
	return maxPt;
}

rg_Point3D rg_BoundingBox3D::getCenterPt() const
{
	rg_Point3D output(0.0,0.0,0.0);

	output= (minPt+maxPt)*0.5;

	return output;
}

rg_REAL    rg_BoundingBox3D::getXLength() const
{	
	rg_REAL output=maxPt.getX()-minPt.getX();

	if ( rg_LT(output,0.0) )
	{
		output=0.0;
	}
	
	return output;
}

rg_REAL    rg_BoundingBox3D::getYLength() const   
{
	rg_REAL output=maxPt.getY()-minPt.getY();

	if ( rg_LT(output,0.0) )
	{
		output=0.0;
	}
	
	return output;
}

rg_REAL    rg_BoundingBox3D::getZLength() const   
{
	rg_REAL output=maxPt.getY()-minPt.getZ();

	if ( rg_LT(output,0.0) )
	{
		output=0.0;
	}
	
	return output;

}

rg_REAL    rg_BoundingBox3D::getLongestLength() const
{
	rg_Point3D longestVector=maxPt-minPt;

	rg_REAL output=0.0;
	if (   rg_GE( longestVector.getX(),0.0)
		&& rg_GE( longestVector.getY(),0.0)
		&& rg_GE( longestVector.getZ(),0.0) )
	{
		output=longestVector.magnitude();
	}

	return output;
}

rg_BoundingBox3D rg_BoundingBox3D::evaluateRelativeOffset(const rg_REAL& ratio) const
{
	if ( isNull() )
	{
		return rg_BoundingBox3D();
	}
//    rg_REAL    dist=getLongestLength()*ratio;
	rg_Point3D vector=maxPt-minPt;
	rg_Point3D tMinPt=minPt-vector*ratio;
	rg_Point3D tMaxPt=maxPt+vector*ratio;

	return rg_BoundingBox3D(tMinPt,tMaxPt);
}

rg_BoundingBox3D rg_BoundingBox3D::evaluateAbsoluteOffset(const rg_REAL& dist) const
{
	if ( isNull() )
	{
		return rg_BoundingBox3D();
	}

	rg_Point3D vector(dist,dist,dist);
	rg_Point3D tMinPt=minPt-vector;
	rg_Point3D tMaxPt=maxPt+vector;

	return rg_BoundingBox3D(tMinPt,tMaxPt);
}


void    rg_BoundingBox3D::setMinPt(const rg_Point3D& tMinPt)
{
	minPt = tMinPt;
}

void    rg_BoundingBox3D::setMaxPt(const rg_Point3D& tMaxPt)
{
	maxPt= tMaxPt;
}

void    rg_BoundingBox3D::setAll(const rg_Point3D& tMinPt, const rg_Point3D& tMaxPt)
{
	minPt = tMinPt;
	maxPt = tMaxPt;
}

void    rg_BoundingBox3D::setAll(const rg_BoundingBox3D& temp)
{
	minPt = temp.minPt;
	maxPt = temp.maxPt;
}

void    rg_BoundingBox3D::reset()
{
	minPt = rg_Point3D(rg_MAX_REAL, rg_MAX_REAL, rg_MAX_REAL);
	maxPt = rg_Point3D(-rg_MAX_REAL, -rg_MAX_REAL, -rg_MAX_REAL);
}

//  operator overloading
rg_BoundingBox3D& rg_BoundingBox3D::operator =(const rg_BoundingBox3D& temp)
{
	if ( this == &temp )
	{
		return *this;
	}
	
	setAll( temp );

	return *this;
}

/*
//  bounding box 3d of curve and surface.
void rg_BoundingBox3D::calculateBoundingBox3DOfBezierSurface(const rg_BzSurface3D& aBzSurface)
{
	rg_REAL top    = -rg_MAX_REAL;	// max of z
	rg_REAL left   = rg_MAX_REAL;		// min of y
	rg_REAL front  = -rg_MAX_REAL;	// max of x
	rg_REAL bottom = rg_MAX_REAL;		// min of z
	rg_REAL right  = -rg_MAX_REAL;	// max of y
	rg_REAL back   = rg_MAX_REAL;		// min of x

	rg_INT uDegree = aBzSurface.getDegreeOfU();
	rg_INT vDegree = aBzSurface.getDegreeOfV();

	for (rg_INDEX i=0; i<uDegree; i++)  
	{
		for (rg_INDEX j=0; j<vDegree; j++)
		{
			rg_Point3D pt = aBzSurface.getCtrlPt(i, j);
			if ( rg_GE(pt.getZ(), top) )	{
				top = pt.getZ();
			}

			if ( rg_LE(pt.getY(), left) )	{
				left = pt.getY();
			}

			if ( rg_GE(pt.getX(), front) )	{
				front = pt.getX();
			}

			if ( rg_LE(pt.getZ(), bottom) )	{
				bottom = pt.getZ();
			}

			if ( rg_GE(pt.getY(), right) )	{
				right = pt.getY();
			}

			if ( rg_LE(pt.getX(), back) )	{
				back = pt.getX();
			}
		}
	}

	setBoundingBox3D( rg_Point3D( front, left, top),
                      rg_Point3D( back, right, bottom) );
}
*/

void rg_BoundingBox3D::contain(const rg_INT& numPts, rg_Point3D* pts)
{
	for (rg_INDEX i=0; i<numPts; i++)
	{
		if ( rg_GT(pts[i].getX(), maxPt.getX()) )	
		{
			maxPt.setX( pts[i].getX() );
		}

		if ( rg_GT(pts[i].getY(), maxPt.getY()) )	
		{
			maxPt.setY( pts[i].getY() );
		}

		if ( rg_GE(pts[i].getZ(), maxPt.getZ()) )	
		{
			maxPt.setZ( pts[i].getZ() );
		}

		if ( rg_LT(pts[i].getX(), minPt.getX()) )	
		{
			minPt.setX( pts[i].getX() );
		}

		if ( rg_LT(pts[i].getY(), minPt.getY()) )	
		{
			minPt.setY( pts[i].getY() );
		}

		if ( rg_LE(pts[i].getZ(), minPt.getZ()) )	
		{
			minPt.setZ( pts[i].getZ() );
		}
	}

}

void rg_BoundingBox3D::contain(const rg_Point3D& pt)
{
		if ( rg_GT(pt.getX(), maxPt.getX()) )	
		{
			maxPt.setX( pt.getX() );
		}

		if ( rg_GT(pt.getY(), maxPt.getY()) )	
		{
			maxPt.setY( pt.getY() );
		}

		if ( rg_GE(pt.getZ(), maxPt.getZ()) )	
		{
			maxPt.setZ( pt.getZ() );
		}

		if ( rg_LT(pt.getX(), minPt.getX()) )	
		{
			minPt.setX( pt.getX() );
		}

		if ( rg_LT(pt.getY(), minPt.getY()) )	
		{
			minPt.setY( pt.getY() );
		}

		if ( rg_LE(pt.getZ(), minPt.getZ()) )	
		{
			minPt.setZ( pt.getZ() );
		}

}


rg_FLAG  rg_BoundingBox3D::doContain(const rg_Point3D& pt) const
{
    if ( isNull() )
    {
        return rg_FALSE;
    }

	rg_FLAG output=    rg_BTOR( minPt.getX(), pt.getX(), maxPt.getX() ) 
                    && rg_BTOR( minPt.getY(), pt.getY(), maxPt.getY() )
                    && rg_BTOR( minPt.getZ(), pt.getZ(), maxPt.getZ() );

    return output;
}


rg_FLAG  rg_BoundingBox3D::isNull() const
{
	rg_FLAG output=    (    rg_GE( minPt.getX() , maxPt.getX() )
		                 && rg_GE( minPt.getY() , maxPt.getY() )
		                 && rg_GE( minPt.getZ() , maxPt.getZ() ) );

    return output;
}

void rg_BoundingBox3D::updateBoxByAddingSphere(const Sphere& sphere)
{
	rg_Point3D center = sphere.getCenter();
	double  radius = sphere.getRadius();

	double xCoordOfPoint = center.getX();
	double yCoordOfPoint = center.getY();
	double zCoordOfPoint = center.getZ();

	// update minimum point
	if(xCoordOfPoint - radius < minPt.getX())
		minPt.setX( xCoordOfPoint - radius );

	if(yCoordOfPoint - radius < minPt.getY())
		minPt.setY( yCoordOfPoint - radius );

	if(zCoordOfPoint - radius < minPt.getZ())
		minPt.setZ( zCoordOfPoint - radius );

	// update maximum point
	if(xCoordOfPoint + radius > maxPt.getX())
		maxPt.setX( xCoordOfPoint + radius );

	if(yCoordOfPoint + radius > maxPt.getY())
		maxPt.setY( yCoordOfPoint + radius );

	if(zCoordOfPoint + radius > maxPt.getZ())
		maxPt.setZ( zCoordOfPoint + radius );
}

void rg_BoundingBox3D::constructBoxByAddingSpheres(Sphere* spheres, const int& numSpheres)
{
	if(spheres == NULL)
		return;

	int i;
	for (i = 0;i < numSpheres;i++)
	{
		updateBoxByAddingSphere(spheres[ i ]);
	}

	delete [] spheres;
}

void rg_BoundingBox3D::constructBoxByAddingSpheres(rg_dList<Sphere>& spheres)
{
	spheres.reset4Loop();
	while (spheres.setNext4Loop())
	{
		Sphere currSphere = spheres.getEntity();
		updateBoxByAddingSphere(currSphere);
	}
}
	

