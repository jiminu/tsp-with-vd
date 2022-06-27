#include "rg_RevolvedSurface.h"
#include "rg_TMatrix3D.h"


rg_RevolvedSurface::rg_RevolvedSurface()
:rg_NURBSplineSurface3D(),generatrix(),axis(rg_Point3D(0.0,0.0,0.0), rg_Point3D(0.0,0.0,1.0)),angle(0.0)
{
}

rg_RevolvedSurface::rg_RevolvedSurface(const rg_Outline3D& tGeneratix)
:rg_NURBSplineSurface3D(),generatrix(tGeneratix),axis(rg_Point3D(0.0,0.0,0.0), rg_Point3D(0.0,0.0,1.0)),angle(0.0)
{
}

rg_RevolvedSurface::rg_RevolvedSurface(const rg_Outline3D& tGeneratix,
                                 const rg_Line3D&   tAxis,
                                 const rg_REAL&     tAngle)
:rg_NURBSplineSurface3D(),generatrix(tGeneratix),axis(tAxis),angle(tAngle)
{
}

rg_RevolvedSurface::rg_RevolvedSurface(const rg_RevolvedSurface& temp)
:rg_NURBSplineSurface3D(temp), generatrix(temp.generatrix), axis(temp.axis), angle( temp.angle)
{
}

rg_RevolvedSurface::~rg_RevolvedSurface()
{
}
// get function
rg_Outline3D rg_RevolvedSurface::getGeneratrix() const
{   
    return generatrix;
}

rg_Line3D   rg_RevolvedSurface::getAxis() const
{
    return axis;
}

rg_REAL  rg_RevolvedSurface::getAngle() const
{
    return angle;
}

// set function
void   rg_RevolvedSurface::setGeneratrix(const rg_Outline3D& tGeneratrix)
{
    generatrix=tGeneratrix;
}

void   rg_RevolvedSurface::setAxis(const rg_Line3D& tAxis)
{
    axis=tAxis;
}

void   rg_RevolvedSurface::setAngle(const rg_REAL& tAngle)
{
    angle=tAngle;
}

void   rg_RevolvedSurface::makeSurface() 
{
    // delete weight vector
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet(); 
   
    // get number of partition
    const rg_REAL PHI=atan(1.0)*4.0;
    const rg_REAL ratioOfPartition=angle/(PHI*0.5);
    rg_INT        numOfPartition=(rg_INT) ratioOfPartition;
    numOfPartition = ( rg_EQ( ratioOfPartition,numOfPartition ) ) ? numOfPartition:numOfPartition+1;

    row = generatrix.getNumOfCtrlPts();
    col = 2*numOfPartition+1;
    const rg_INT  uOrder= generatrix.getOrder();
    const rg_INT  vOrder= 3;

    
    const rg_REAL stepAngle=angle/numOfPartition;
    const rg_REAL halfStepAngle=stepAngle/2.0;
    const rg_REAL cosOfHalfStepAngle=cos( halfStepAngle);
    rg_REAL*  vKnots=new rg_REAL[col+vOrder];

    const rg_Point3D localZAxis=axis.getUnitDirection();

    rg_Point3D* center = new rg_Point3D[row];
    rg_REAL*  radius=new rg_REAL[row];
    
    rg_INT i = 0;
	for( i=0; i < row; i++ )
    {
        rg_Point3D startPt=generatrix.getCtrlPt(i);
        rg_Point3D centerOfAxis=axis.getSP();
        rg_Point3D roof=(startPt-centerOfAxis);
        rg_REAL height=roof%localZAxis;
        center[i]=centerOfAxis+height*localZAxis;
        radius[i]=(startPt-center[i]).magnitude();
    }

    
    // make knot vector 
    for( i=0 ; i < vOrder; i++)
    {
        vKnots[0]=0.0;
        vKnots[col+i]=1.0;
    }
    for( i=0; i < numOfPartition; i++ ) 
    {
        vKnots[2*i+1] = i/(rg_REAL)numOfPartition;
        vKnots[2*i+2] = i/(rg_REAL)numOfPartition;
    }
    rg_REAL* uKnots= generatrix.getKnotVector();

    rg_BSplineSurface3D::setOrderOfSurface(uOrder,vOrder);
    rg_NURBSplineSurface3D::setControlNetAndWeight(row,col);
 

    for( i=0; i < row; i++ )
    {
        rg_TMatrix3D mat;
        mat.rotateArbitraryAxis(localZAxis,angle);
        const rg_Point3D startPt=generatrix.getCtrlPt(i);
        const rg_Point3D end=mat*(startPt-center[i])+center[i];
        for(rg_INT j=0; j < numOfPartition; j++ )
        {
            rg_TMatrix3D transform1;
            rg_TMatrix3D transform2;
            transform1.rotateArbitraryAxis( localZAxis, j*stepAngle );
            transform2.rotateArbitraryAxis( localZAxis, (j+1)*stepAngle);

            const rg_Point3D front=transform1*(startPt-center[i])+center[i];
            const rg_Point3D rear=transform2*(startPt-center[i])+center[i];
            const rg_Point3D bisectVector=(front+rear)/2.0-center[i];
			rg_Point3D middle;
			if(rg_EQ(bisectVector.magnitude(), 0))
			{
				middle = (front+rear)/2.;
			}
			else
			{
				middle = radius[i]/cosOfHalfStepAngle*bisectVector.getUnitVector()+ center[i]; 
			}


            rg_BSplineSurface3D::setPointOnControlNet(i,2*j,front);
            rg_BSplineSurface3D::setPointOnControlNet(i,2*j+1,middle);

            rg_REAL frontWeight=generatrix.getWeight(i);
            rg_REAL rearWeight=generatrix.getWeight(i)*cosOfHalfStepAngle;
            rg_NURBSplineSurface3D::setWeight(i,2*j,frontWeight);
            rg_NURBSplineSurface3D::setWeight(i,2*j+1,rearWeight);
        }

        rg_BSplineSurface3D::setPointOnControlNet(i,col-1,end);
        rg_NURBSplineSurface3D::setWeight(i,col-1,generatrix.getWeight(i));

    }
    rg_NUBSplineSurface3D::setKnotVectorOfU(row+uOrder,uKnots);
    rg_NUBSplineSurface3D::setKnotVectorOfV(col+vOrder,vKnots);
	delete[] uKnots;
    delete[] vKnots;
    delete[] center;
    delete[] radius;

}
    
// operator overloading
rg_RevolvedSurface& rg_RevolvedSurface::operator=(const rg_RevolvedSurface& temp)
{
    if ( &temp == this )
    {
        return *this;
    }
    rg_NURBSplineSurface3D::operator=(temp);
    generatrix=temp.generatrix;
    axis=temp.axis; //
    angle=temp.angle;

    return *this;
}


