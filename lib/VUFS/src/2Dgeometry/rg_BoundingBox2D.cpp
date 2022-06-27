#include "rg_BoundingBox2D.h"
#include "Ellipse2D.h"

// construtors and destructor
rg_BoundingBox2D::rg_BoundingBox2D()
{
    min=InfinitPt2D;
    max=-InfinitPt2D;
}

rg_BoundingBox2D::rg_BoundingBox2D(const rg_Point2D& pt)
{
    min=pt;
    max=pt;
}
rg_BoundingBox2D::rg_BoundingBox2D(const rg_Point2D& tMinPt,const rg_Point2D& tMaxPt)
{
    min=tMinPt;
    max=tMaxPt;
}

rg_BoundingBox2D::rg_BoundingBox2D(rg_dList<rg_Circle2D> &circles)
{
	constructBoxByAddingCircles(circles);
}

rg_BoundingBox2D::rg_BoundingBox2D(const rg_BoundingBox2D& temp)
{
    min=temp.min;
    max=temp.max;
}

rg_BoundingBox2D::~rg_BoundingBox2D()
{
}

rg_Point2D rg_BoundingBox2D::getMinPt() const
{
    return min;
}

rg_Point2D rg_BoundingBox2D::getMaxPt() const
{
    return max;
}

void       rg_BoundingBox2D::setMinPt(const rg_Point2D& tMin)
{
	min=tMin;
}

void       rg_BoundingBox2D::setMaxPt(const rg_Point2D& tMax)
{
	max=tMax;
}


rg_Point2D rg_BoundingBox2D::getCenterPt() const
{
    return (min + max)*0.5;
}



rg_REAL    rg_BoundingBox2D::evaluateXLength() const
{
	rg_REAL output=max.getX()-min.getX();

	if ( rg_LT(output,0.0) )
	{
		output=0.0;
	}

	return output;
}

rg_REAL    rg_BoundingBox2D::evaluateLongestLength() const
{
    rg_REAL output=0.0;
    if ( !isNull() )
    {
        output=min.distance(max);
    }

    return output;
}

rg_REAL    rg_BoundingBox2D::evaluateYLength() const
{
	rg_REAL output=max.getY()-min.getY();

	if ( rg_LT(output,0.0) )
	{
		output=0.0;
	}

	return output;
}

rg_REAL    rg_BoundingBox2D::evaluateAspectRatio() const
{
    rg_Point2D vector=max-min;

    return vector.getX()/vector.getY();
}

void rg_BoundingBox2D::contain(const rg_Point2D& pt)
{
    if ( min == InfinitPt2D )
    {
        min=pt;
        max=pt;
    }
    else
    {
        const rg_REAL ptX=pt.getX();
        const rg_REAL ptY=pt.getY();
        const rg_REAL minX=min.getX();
        const rg_REAL minY=min.getY();
        const rg_REAL maxX=max.getX();
        const rg_REAL maxY=max.getY();

        if ( minX > ptX )
        {
            min.setX(ptX);
        }
        else if( maxX < ptX )
        {
            max.setX(ptX);
        }

        if ( minY > ptY )
        {
            min.setY(ptY);
        }
        else if( maxY < ptY )
        {
            max.setY(ptY);
        }
    }
}

void rg_BoundingBox2D::contain(const rg_BoundingBox2D& tBox)
{
    if ( tBox.isNull() )
    {
        return;
    }

    if ( isNull() )
    {
        *this=tBox;
        return;
    }
    
    if ( rg_GT( min.getX(), tBox.min.getX() ) )
    {
        min.setX(tBox.min.getX());
    }
    if ( rg_GT( min.getY(), tBox.min.getY() ) )
    {
        min.setY(tBox.min.getY());
    }

    if ( rg_LT( max.getX(), tBox.max.getX() ) )
    {
        max.setX(tBox.max.getX());
    }
    if ( rg_LT( max.getY(), tBox.max.getY() ) )
    {
        max.setY(tBox.max.getY());
    }
}



rg_BoundingBox2D rg_BoundingBox2D::evaluateRelativeOffset(const rg_REAL& ratio) const
{
	if ( isNull() )
	{
		return rg_BoundingBox2D();
	}

    rg_REAL    dist=(max.distance(min))*ratio;
	rg_Point2D vector(dist,dist);
	rg_Point2D tMinPt=min-vector;
	rg_Point2D tMaxPt=max+vector;

	return rg_BoundingBox2D(tMinPt,tMaxPt);
}

rg_BoundingBox2D rg_BoundingBox2D::evaluateAbsoluteOffset(const rg_REAL& dist) const
{
	if ( isNull() )
	{
		return rg_BoundingBox2D();
	}

	rg_Point2D vector(dist,dist);
	rg_Point2D tMinPt=min-vector;
	rg_Point2D tMaxPt=max+vector;

	return rg_BoundingBox2D(tMinPt,tMaxPt);
}

void   rg_BoundingBox2D::refitAspectRatio(const rg_REAL& ratio)
{
	const rg_REAL currentRatio=evaluateAspectRatio();


	if (   rg_EQ(ratio,currentRatio)
		|| rg_LE(ratio,0.0) ) 
	{
		return;
	}

	rg_Point2D halfExtension;
	rg_REAL    width=evaluateXLength();
	//rg_REAL    height=evaluateYLength();

	rg_REAL halfExtendedSize=0.5*width*(ratio/currentRatio -1.0 );
	halfExtension = rg_Point2D(halfExtendedSize,0.0);

	/*
	if ( rg_GT( ratio, currentRatio ) )
	{
		//  case that y must be exteded.
		rg_REAL halfExtendedSize=0.5*width*(ratio/currentRatio -1.0 );
		halfExtension = rg_Point2D(halfExtendedSize,0.0);
	}
	else 
	{
		rg_REAL halfExtendedSize=0.5*height*(currentRatio/ratio -1.0 );
		halfExtension = rg_Point2D(0.0,halfExtendedSize);
	}
*/
	min=min-halfExtension;
	max=max+halfExtension;
}




        
void rg_BoundingBox2D::contain(const rg_Point2D*  pts,
                          const rg_INT&      size )
{
    if ( pts != rg_NULL )
    {
        for( rg_INT i=0 ; i < size ; i++ )
        {
            contain( pts[i] );
        }
    }

/*  if ( pts != rg_NULL )
    {
        rg_REAL minX=pts[0].getX();
        rg_REAL maxX=min;
        rg_REAL minY=pts[0].getY();
        rg_REAL maxY=minY;

        if ( min == rg_NULL )
        {
            min= new rg_Point2D(pts[0]);
            max= new rg_Point2D(pts[0]);
        }

        for( rg_INT i=1 ; i < size ; i++ )
        {
            rg_REAL thisX=pts[i].getX();
            rg_REAL thisY=pts[i].getY();

            if( minX > thisX )
            {
                minX=thisX;
            }
            else if ( maxX < thisX )
            {
                maxX=thisX;
            }

            if ( minY > ptY )
            {
                minY=thisY;
            }
            else if( maxY < ptY )
            {
                maxY=thisY;
            }
        }
        
        contain(rg_Point2D(minX,minY));
        contain(rg_Point2D(maxX,maxY));
    }
*/
}



void rg_BoundingBox2D::constructBoxByAddingCircles(rg_dList<rg_Circle2D>& circles)
{
	circles.reset4Loop();
	while (circles.setNext4Loop())
	{
		rg_Circle2D currCircle = circles.getEntity();
		updateBoxByAddingCircle(currCircle);
	}
}

void rg_BoundingBox2D::updateBoxByAddingCircle(const rg_Circle2D& circle)
{
	rg_Point2D center = circle.getCenterPt();
	rg_REAL    radius = circle.getRadius();

	rg_REAL xCoordOfPoint = center.getX();
	rg_REAL yCoordOfPoint = center.getY();

	// update minimum point
	if(xCoordOfPoint - radius < min.getX())
		min.setX( xCoordOfPoint - radius );

	if(yCoordOfPoint - radius < min.getY())
		min.setY( yCoordOfPoint - radius );

	// update maximum point
	if(xCoordOfPoint + radius > max.getX())
		max.setX( xCoordOfPoint + radius );

	if(yCoordOfPoint + radius > max.getY())
		max.setY( yCoordOfPoint + radius );
}


void rg_BoundingBox2D::updateBoxByAddingEllipse(const Ellipse2D& ellipse)
{
	rg_Point2D center = ellipse.get_centerPt();
	rg_REAL    majorAxisHalfLength, minorAxisHalfLength;
	ellipse.get_semi_axis_length(majorAxisHalfLength, minorAxisHalfLength);

	rg_REAL xCoordOfPoint = center.getX();
	rg_REAL yCoordOfPoint = center.getY();

	// update minimum point
	if (xCoordOfPoint - majorAxisHalfLength < min.getX())
		min.setX(xCoordOfPoint - majorAxisHalfLength);

	if (yCoordOfPoint - minorAxisHalfLength < min.getY())
		min.setY(yCoordOfPoint - minorAxisHalfLength);

	// update maximum point
	if (xCoordOfPoint + majorAxisHalfLength > max.getX())
		max.setX(xCoordOfPoint + majorAxisHalfLength);

	if (yCoordOfPoint + minorAxisHalfLength > max.getY())
		max.setY(yCoordOfPoint + minorAxisHalfLength);

}


void rg_BoundingBox2D::get_boundary_points_in_CCW(list<rg_Point2D>& boundaryPoints)
{
	boundaryPoints.push_back(min);
	boundaryPoints.push_back(rg_Point2D(max.getX(), min.getY()));
	boundaryPoints.push_back(max);
	boundaryPoints.push_back(rg_Point2D(min.getX(), max.getY()));
}

void rg_BoundingBox2D::get_boundary_points_in_CW(list<rg_Point2D>& boundaryPoints)
{
	boundaryPoints.push_back(min);
	boundaryPoints.push_back(rg_Point2D(min.getX(), max.getY()));
	boundaryPoints.push_back(max);
	boundaryPoints.push_back(rg_Point2D(max.getX(), min.getY()));		
}


void rg_BoundingBox2D::reset()
{
    min= InfinitPt2D; 
	max= -InfinitPt2D;
}


// operator overloading
rg_BoundingBox2D& rg_BoundingBox2D::operator=(const rg_BoundingBox2D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    min=temp.min;
    max=temp.max;

    return *(this);
}

// ETC functions
rg_FLAG rg_BoundingBox2D::isOverlapped(const rg_BoundingBox2D& temp) const
{
	if(   ( rg_GT(min.getX(), temp.max.getX(), resNeg12) )
       || ( rg_LT(max.getX(), temp.min.getX(), resNeg12) )
       || ( rg_GT(min.getY(), temp.max.getY(), resNeg12) )
       || ( rg_LT(max.getY(), temp.min.getY(), resNeg12) ) )

	{
		return rg_FALSE;
	}
	else 
	{
		return rg_TRUE;
	}
}

rg_FLAG rg_BoundingBox2D::isNull() const
{
    if ( min == InfinitPt2D )
    {
        return rg_TRUE;
    }
    else
    {
        return rg_FALSE;
    }
}


bool rg_BoundingBox2D::doesContain(const rg_Point2D& point, const double& res) const
{
    if (   rg_GE(point.getX(), min.getX(), res)
        && rg_LE(point.getX(), max.getX(), res)
        && rg_GE(point.getY(), min.getY(), res)
        && rg_LE(point.getY(), max.getY(), res))
    {
        return true;
    }
    else
    {
        return false;
    }
}


/*

rg_Box3D *rg_IntervalIntersector::splitBox(const rg_Box3D& box, const rg_BzCurve2D& curve)
{
	rg_Box3D *returnBox = new rg_Box3D[2];

	returnBox[0].SP = box.SP;
	returnBox[1].EP = box.EP;

	returnBox[0].ST = box.ST;
	returnBox[1].ET = box.ET;

	rg_REAL  midParam = (box.ST+box.ET)/2.0;
	rg_Point2D midPoint = curve.evaluatePt(midParam);
	returnBox[0].EP  = returnBox[1].SP = midPoint;
	returnBox[0].ET  = returnBox[1].ST = midParam;

	return returnBox;
}

void rg_IntervalIntersector::changingBox(rg_Box3D &box)
{
	rg_Point2D leftBottom;
	rg_Point2D rightTop;
	if( rg_LT(box.SP.getX(), box.EP.getX()) )
	{
		if( rg_LT(box.SP.getY(), box.EP.getY()) )
		{
			return;
		}
		else 
		{
			leftBottom.setPoint(box.SP.getX(), box.EP.getY());
			rightTop.setPoint(  box.EP.getX(), box.SP.getY());
		}
	}

	else
	{
		if( rg_LT(box.SP.getY(), box.EP.getY()) )
		{
			leftBottom.setPoint(box.EP.getX(), box.SP.getY());
			rightTop.setPoint(  box.SP.getX(), box.EP.getY());
		}
		else
		{
			leftBottom.setPoint(box.EP);
			rightTop.setPoint(  box.SP);
		}
	}
	box.SP = leftBottom;
	box.EP = rightTop;
}

BOOL rg_IntervalIntersector::isOverlapping(rg_Box3D box_s, rg_Box3D box_t)
{
	changingBox(box_s);
	changingBox(box_t);

	if(    rg_GT(box_s.SP.getX(), box_t.SP.getX(), resNeg12)
		&& rg_GT(box_s.SP.getX(), box_t.EP.getX(), resNeg12) )
	{
		return rg_FALSE;
	}
	else if( rg_LT(box_s.EP.getX(), box_t.SP.getX(), resNeg12)
		  && rg_LT(box_s.EP.getX(), box_t.EP.getX(), resNeg12) )
	{
		return rg_FALSE;
	}
	else if(rg_GT(box_s.SP.getY(), box_t.SP.getY(), resNeg12)
		 && rg_GT(box_s.SP.getY(), box_t.EP.getY(), resNeg12) )
	{
		return rg_FALSE;
	}
	else if( rg_LT(box_s.EP.getY(), box_t.SP.getY(), resNeg12)
		  && rg_LT(box_s.EP.getY(), box_t.EP.getY(), resNeg12) )
	{
		return rg_FALSE;
	}
	else 
	{
		return rg_TRUE;
	}
}
*/
