//********************************************************************
//
//	  FILENAME    : BSplineSurface.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_BSplineSurface3D 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************


#include <stdio.h>

#include "rg_RelativeOp.h"

#include "rg_BSplineSurface3D.h"

////	Constructor & Destructor
rg_BSplineSurface3D::rg_BSplineSurface3D()
: rg_Surface(), u_order(4), v_order(4),
  rowOfControlNet(0), columnOfControlNet(0), control_net(rg_NULL)
{
}

rg_BSplineSurface3D::rg_BSplineSurface3D( const unsigned rg_INT &newID, 
								    const rg_Planarity    &newPlanarity )
: rg_Surface(newID, newPlanarity), u_order(4), v_order(4),
  rowOfControlNet(0), columnOfControlNet(0), control_net(rg_NULL)
{
}

rg_BSplineSurface3D::rg_BSplineSurface3D( const unsigned rg_INT &newID, 
								    const rg_Planarity    &newPlanarity, 
						            const rg_INT          &row, 
									const rg_INT          &col )
: rg_Surface(newID, newPlanarity), u_order(4), v_order(4),
  rowOfControlNet(row), columnOfControlNet(col)
{
	control_net = new rg_Point3D* [rowOfControlNet];
	for (rg_INT i=0; i<rowOfControlNet; i++)
		control_net[i] = new rg_Point3D [columnOfControlNet];
}

rg_BSplineSurface3D::rg_BSplineSurface3D( const unsigned rg_INT &newID, 
								    const rg_Planarity    &newPlanarity, 
						            const rg_ORDER        &uOrder, 
									const rg_ORDER        &vOrder )
: rg_Surface(newID, newPlanarity), u_order(uOrder), v_order(vOrder),
  rowOfControlNet(0), columnOfControlNet(0), control_net(rg_NULL)
{
}

rg_BSplineSurface3D::rg_BSplineSurface3D( const rg_INT &row, 
								    const rg_INT &col )
: rg_Surface(), u_order(4), v_order(4),
  rowOfControlNet(row), columnOfControlNet(col) 
{
	control_net = new rg_Point3D* [rowOfControlNet];
	for (rg_INT i=0; i<rowOfControlNet; i++)
		control_net[i] = new rg_Point3D [columnOfControlNet];
}


rg_BSplineSurface3D::rg_BSplineSurface3D( const rg_ORDER &uOrder, 
								    const rg_ORDER &vOrder )
: rg_Surface(), u_order(uOrder), v_order(vOrder),
  rowOfControlNet(0), columnOfControlNet(0), control_net(rg_NULL)
{
}

rg_BSplineSurface3D::rg_BSplineSurface3D( const rg_INT   &row, 
								    const rg_INT   &col,
						            const rg_ORDER &uOrder, 
									const rg_ORDER &vOrder )
: rg_Surface(), u_order(uOrder), v_order(vOrder),
  rowOfControlNet(row), columnOfControlNet(col)
{
	control_net = new rg_Point3D* [rowOfControlNet];
	for (rg_INT i=0; i<rowOfControlNet; i++)
		control_net[i] = new rg_Point3D [columnOfControlNet];
}

rg_BSplineSurface3D::rg_BSplineSurface3D( const unsigned rg_INT &newID, 
								    const rg_Planarity    &newPlanarity, 
						            const rg_INT          &row, 
									const rg_INT          &col, 
						            const rg_ORDER        &uOrder, 
									const rg_ORDER        &vOrder )
: rg_Surface(newID, newPlanarity), u_order(uOrder), v_order(vOrder),
  rowOfControlNet(row), columnOfControlNet(col)
{
	control_net = new rg_Point3D* [rowOfControlNet];
	for (rg_INT i=0; i<rowOfControlNet; i++)
		control_net[i] = new rg_Point3D [columnOfControlNet];
}

////  Constructor      : March 13 1997
rg_BSplineSurface3D::rg_BSplineSurface3D( const unsigned rg_INT &newID, 
								    const rg_Planarity    &newPlanarity, 
						            const rg_INT          &row, 
									const rg_INT          &col, 
						            const rg_ORDER        &uOrder, 
									const rg_ORDER        &vOrder,
                                    rg_Point3D**             ctrlNet)
: rg_Surface(newID, newPlanarity), 
  u_order(uOrder), 
  v_order(vOrder),
  rowOfControlNet(row),
  columnOfControlNet(col)
{
	control_net = new rg_Point3D* [row];
	for (rg_INT i=0; i<row; i++)
    {
		control_net[i] = new rg_Point3D [col];
       for (rg_INT j=0; j<col; j++)
            control_net[i][j] = ctrlNet[i][j];
    }
}

////  Copy Constructor : March 13 1997
rg_BSplineSurface3D::rg_BSplineSurface3D( const rg_BSplineSurface3D &surface )
:  rg_Surface(surface.getID(), surface.getPlanarity() ),
   u_order(surface.u_order), 
   v_order(surface.v_order),
   rowOfControlNet(surface.rowOfControlNet),
   columnOfControlNet(surface.columnOfControlNet)
{
	control_net = new rg_Point3D* [surface.rowOfControlNet];
	for (rg_INT i=0; i<surface.rowOfControlNet; i++)
    {
		control_net[i] = new rg_Point3D [surface.columnOfControlNet];
        for (rg_INT j=0; j<surface.columnOfControlNet; j++)
            control_net[i][j] = surface.control_net[i][j];
    }
}

rg_BSplineSurface3D::~rg_BSplineSurface3D()
{
	if (control_net != rg_NULL)
	{
		for (rg_INT i=0; i<rowOfControlNet; i++)
			delete [] control_net[i];

		delete [] control_net;
	}
}

////	Get Functions.
rg_INT rg_BSplineSurface3D::isValidSurface() const
{
	if (control_net == rg_NULL)
		return rg_FALSE;

	return rg_TRUE;
}

rg_INT rg_BSplineSurface3D::getRowOfControlNet() const
{
	return rowOfControlNet;
}

rg_INT rg_BSplineSurface3D::getColumnOfControlNet() const
{
	return columnOfControlNet;
}

rg_ORDER rg_BSplineSurface3D::getOrderOfU() const
{
	return u_order;
}

rg_ORDER rg_BSplineSurface3D::getOrderOfV() const
{
	return v_order;
}

rg_Point3D rg_BSplineSurface3D::getPointOnControlNet( const rg_INT &row, 
											  const rg_INT &col ) const
{
	return control_net[row][col];
}

rg_Point3D** rg_BSplineSurface3D::getControlNet() const
{
    return control_net;
}

rg_REAL rg_BSplineSurface3D::getKnotValueOfU( const rg_INT &kIndex ) const
{
	rg_INT row    = rowOfControlNet;
	rg_INT uOrder = (rg_INT) u_order;

	if ( kIndex < ( row + uOrder) )
	{
		if ( kIndex < uOrder )
			return  0.;
		else if (kIndex < row )
			return (kIndex - uOrder + 1.) / (row - uOrder + 1.);
		else
			return 1.;
	}
	else
		return -1.;
}

rg_REAL*  rg_BSplineSurface3D::getKnotVectorOfU() const
{
    const rg_INT uOrder= rg_BSplineSurface3D::getOrderOfU();
    const rg_INT rowOfCtrlPts= rg_BSplineSurface3D::getRowOfControlNet();
    const rg_INT numOfKnotsOfU=uOrder+rowOfCtrlPts;
    rg_REAL* output= new rg_REAL[numOfKnotsOfU];

    for( rg_INT i=0; i <numOfKnotsOfU; i++ )
    {
        output[i]=rg_BSplineSurface3D::getKnotValueOfU(i);
    }
    return output;
}

rg_REAL*  rg_BSplineSurface3D::getKnotVecotrOfV() const
{
    const rg_INT vOrder= rg_BSplineSurface3D::getOrderOfV();
    const rg_INT colOfCtrlPts= rg_BSplineSurface3D::getColumnOfControlNet();
    const rg_INT numOfKnotsOfV=vOrder+colOfCtrlPts;
    rg_REAL* output= new rg_REAL[numOfKnotsOfV];

    for( rg_INT i=0; i <numOfKnotsOfV; i++ )
    {
        output[i]=rg_BSplineSurface3D::getKnotValueOfV(i);
    }
    return output;
}


rg_REAL rg_BSplineSurface3D::getKnotValueOfV( const rg_INT &kIndex ) const
{
	rg_INT col    = columnOfControlNet;
	rg_INT vOrder = (rg_INT) v_order;

	if ( kIndex < (col + vOrder) )
	{
		if ( kIndex < vOrder )
			return  0.;
		else if (kIndex < col )
			return (kIndex - vOrder + 1.) / (col + vOrder + 1.);
		else
			return 1.;
	}
	else
		return -1.;
}


////	Set Functions.
void rg_BSplineSurface3D::setControlNet( const rg_INT &row, 
									  const rg_INT &col )
{
	if (control_net != rg_NULL)
	{
		for (rg_INT i=0; i<rowOfControlNet; i++)
			delete [] control_net[i];

		delete [] control_net;
	}

	rowOfControlNet	   = row;
	columnOfControlNet = col;

	control_net = new rg_Point3D* [rowOfControlNet];
	for (rg_INT i=0; i<rowOfControlNet; i++)
		control_net[i] = new rg_Point3D [columnOfControlNet];
}

void rg_BSplineSurface3D::setControlNet( const rg_INT &row, 
									  const rg_INT &col, 
									  rg_Point3D**   controlNet)
{
	if (control_net != rg_NULL)
	{
		for (rg_INT i=0; i<rowOfControlNet; i++)
			delete [] control_net[i];

		delete [] control_net;
	}

	rowOfControlNet		= row;
	columnOfControlNet	= col;

    control_net = new rg_Point3D* [row];
    for (rg_INT i=0; i<row; i++)
    {
        control_net[i] = new rg_Point3D [col];
        for (rg_INT j=0; j<col; j++)
            control_net[i][j] = controlNet[i][j];
    }
}

void rg_BSplineSurface3D::setOrderOfU( const rg_ORDER &uOrder )
{
	u_order = uOrder;
}

void rg_BSplineSurface3D::setOrderOfV( const rg_ORDER &vOrder )
{
	v_order = vOrder;
}

void rg_BSplineSurface3D::setOrderOfSurface( const rg_ORDER &uOrder, 
										  const rg_ORDER &vOrder )
{
	u_order = uOrder;
	v_order = vOrder;
}

void rg_BSplineSurface3D::setPointOnControlNet( const rg_INT   &row,	
											 const rg_INT   &col, 
											 const rg_Point3D &point )
{
	control_net[row][col] = point;
}


////	Operating & Calculating.
rg_REAL rg_BSplineSurface3D::evaluateBasisFuncU( const rg_INDEX     &index, 
                                           const rg_PARAMETER &u,
                                           const rg_ORDER     &uOrder)
{
	rg_INT row    = rowOfControlNet - 1;
//    if (uOrder == -1)
//        uOrder = (rg_INT) u_order;
	
	if (    ( index == 0 && rg_EQ( u, getKnotValueOfU(0) ) )
		 || ( index == row && rg_EQ( u, getKnotValueOfU(row + uOrder) ) ) )
	{
		return 1.0;
	}
	else if (    rg_LT( u, getKnotValueOfU(index) )
		      || rg_GE( u, getKnotValueOfU(index + uOrder) ) )
	{
		return 0.0;
	}
	else
	{}

	rg_REAL Uleft  = 0;
	rg_REAL Uright = 0;
	rg_REAL temp   = 0;
	rg_REAL saved  = 0;

	rg_REAL* triN = new rg_REAL [uOrder];
	
	rg_INT i = 0;
	for (i=0; i<uOrder; i++)
	{
		if (    rg_LT( u, getKnotValueOfU(index + i + 1) ) 
			 && rg_GE( u, getKnotValueOfU(index + i) ) )
			triN[i] = 1.0;
		else 
			triN[i] = 0.0;
	}

	for (i=1; i<uOrder; i++)
	{
		if (triN[0] == 0.0) 
			saved = 0.0;
		else 
		{
			saved = ( (u - getKnotValueOfU(index)) * triN[0] )
			        / ( getKnotValueOfU(index+i) - getKnotValueOfU(index) );
		}
	
		for (rg_INT j=0; j<(uOrder-i); j++)
		{
			Uleft  = getKnotValueOfU(index + j + 1);
			Uright = getKnotValueOfU(index + j + i + 1);
		
			if (triN[j+1] == 0.0)
			{
				triN[j] = saved;
				saved   = 0.0;
			}
			else
			{
				temp    = triN[j+1] / (Uright - Uleft);
				triN[j] = saved + (Uright - u) * temp;
				saved   = (u - Uleft) * temp;
			}
		}
	}
	
	rg_REAL returnValue = triN[0];
	delete [] triN;

	return returnValue;
}

rg_REAL rg_BSplineSurface3D::evaluateBasisFuncV( const rg_INDEX     &index, 
                                           const rg_PARAMETER &v,
                                           const rg_ORDER     &vOrder)
{
	rg_INT col    = columnOfControlNet - 1;
//    if (vOrder == -1)
//        vOrder = (rg_INT) v_order;

	if (    ( index == 0 && rg_EQ( v, getKnotValueOfV(0) ) )
		 || ( index == col && rg_EQ( v, getKnotValueOfV(col + vOrder) ) ) )
	{
		return 1.0;
	}
	else if (    rg_LT( v, getKnotValueOfV(index) )
		      || rg_GE( v, getKnotValueOfV(index + vOrder) ) )
	{
		return 0.0;
	}
	else
	{}

	rg_REAL Uleft  = 0;
	rg_REAL Uright = 0;
	rg_REAL temp   = 0;
	rg_REAL saved  = 0;

	rg_REAL* triN = new rg_REAL [v_order];
	
	rg_INT i = 0;
	for (i=0; i<vOrder; i++)
	{
		if (    rg_LT( v, getKnotValueOfV(index + i + 1) ) 
			 && rg_GE( v, getKnotValueOfV(index + i) ) )
		{
			triN[i] = 1.0;
		}
		else 
			triN[i] = 0.0;
	}

	for (i=1; i<vOrder; i++)
	{
		if (triN[0] == 0.0) 
			saved = 0.0;
		else 
			saved = ( (v - getKnotValueOfV(index)) * triN[0] )
			        / ( getKnotValueOfV(index + i) - getKnotValueOfV(index) );
	
		for (rg_INT j=0; j<(vOrder-i); j++)
		{
			Uleft  = getKnotValueOfV(index + j + 1);
			Uright = getKnotValueOfV(index + j + i + 1);

			if (triN[j+1] == 0.0)
			{
				triN[j] = saved;
				saved = 0.0;
			}
			else
			{
				temp    = triN[j+1] / (Uright - Uleft);
				triN[j] = saved + (Uright - v) * temp;
				saved   = (v - Uleft) * temp;
			}
		}
	}
	
	rg_REAL returnValue = triN[0];
	delete [] triN;

	return returnValue;
}

//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D rg_BSplineSurface3D::evaluatePt( const rg_PARAMETER &u, 
                                    const rg_PARAMETER &v )
{
    rg_Point3D ptOnSurface;

    rg_REAL uBasisValue = 0.;
    rg_REAL vBasisValue = 0.;

    for	(rg_INT i=0; i<rowOfControlNet; i++)
    {
        uBasisValue = evaluateBasisFuncU(i, u, u_order);
        if ( rg_NZERO(uBasisValue) )
        {
            for (rg_INT j=0; j<columnOfControlNet; j++)
    		{
	    		vBasisValue = evaluateBasisFuncV(j, v, v_order);
                if ( rg_NZERO(vBasisValue) )
                    ptOnSurface += control_net[i][j] * uBasisValue * vBasisValue;

            }
        }
	}

	return ptOnSurface;
}

////    Derivative.
rg_BSplineSurface3D rg_BSplineSurface3D::derivativeSurfaceOfU()
{
    rg_BSplineSurface3D derivativeSurface( rowOfControlNet-1,
                                        columnOfControlNet,
                                        u_order-1,
                                        v_order );

    rg_Point3D newCtrlPt;
    for (rg_INT i=0; i<rowOfControlNet-1; i++)
    {
        for (rg_INT j=0; j<columnOfControlNet; j++)
        {
            newCtrlPt 
                = ( u_order-1 )
                  * ( control_net[i+1][j] - control_net[i][j])
                  / ( getKnotValueOfU(i+u_order) - getKnotValueOfU(i+1));

            derivativeSurface.control_net[i][j] = newCtrlPt;
        }
    }     

    return derivativeSurface;
}

rg_BSplineSurface3D rg_BSplineSurface3D::derivativeSurfaceOfV()
{
    rg_BSplineSurface3D derivativeSurface( rowOfControlNet,
                                        columnOfControlNet-1,
                                        u_order,
                                        v_order-1 );

    rg_Point3D newCtrlPt;
    for (rg_INT i=0; i<rowOfControlNet; i++)
    {
        for (rg_INT j=0; j<columnOfControlNet-1; j++)
        {
            newCtrlPt 
                = ( v_order-1 )
                  * ( control_net[i][j+1] - control_net[i][j])
                  / ( getKnotValueOfV(j+v_order) - getKnotValueOfV(j+1));

            derivativeSurface.control_net[i][j] = newCtrlPt;
        }
    }     

    return derivativeSurface;
}

rg_BSplineSurface3D rg_BSplineSurface3D::derivativeSurfaceOfUV()
{
    rg_BSplineSurface3D derivativeSurface( rowOfControlNet-1,
                                        columnOfControlNet-1,
                                        u_order-1,
                                        v_order-1 );

    rg_Point3D newCtrlPt;
    for (rg_INT i=0; i<rowOfControlNet-1; i++)
    {
        for (rg_INT j=0; j<columnOfControlNet-1; j++)
        {
            newCtrlPt 
                = ( u_order-1 )*( v_order-1 )
                  * ( control_net[i+1][j+1] - control_net[i][j+1] 
                      - control_net[i+1][j] + control_net[i][j] )
                  / ( ( getKnotValueOfU(i+u_order)-getKnotValueOfV(i+1))
                      * ( getKnotValueOfV(j+v_order)-getKnotValueOfV(j+1)) );

            derivativeSurface.control_net[i][j] = newCtrlPt;
        }
    }     

    return derivativeSurface;
}

////    Operator Overloading.

//  April 7 1997 : made
rg_BSplineSurface3D& rg_BSplineSurface3D::operator =(const rg_BSplineSurface3D &surface)
{
    if ( this == &surface )
        return *this;

    rg_Surface::setID( surface.getID() );
    rg_Surface::setPlanarity( surface.getPlanarity() );

    //  Set of control net.
    if (control_net != rg_NULL)
	{
		for (rg_INT i=0; i<rowOfControlNet; i++)
			delete [] control_net[i];

		delete [] control_net;
	}

    rowOfControlNet    = surface.rowOfControlNet;
    columnOfControlNet = surface.columnOfControlNet;

    control_net = new rg_Point3D* [rowOfControlNet];
    for (rg_INT i=0; i<rowOfControlNet; i++)
    {
        control_net[i] = new rg_Point3D[columnOfControlNet];
        for (rg_INT j=0; j<columnOfControlNet; j++)
            control_net[i][j] = surface.control_net[i][j];
    }

    //  Set of order.
    u_order = surface.u_order;
    v_order = surface.v_order;

    return *this;
}

///	File-Out Function.

//void rg_BSplineSurface3D::fileOut(const char* fileName)
//{
//	FILE* fout;
//	if ((fout = fopen(fileName, "a+")) == rg_NULL)
//		return;
//
//	fprintf( fout, "u_degree\n" );
//	fprintf( fout, "%d\n", u_order-1 );
//	fprintf( fout, "v_degree\n" );
//	fprintf( fout, "%d\n", v_order-1 );
//
//	for (rg_INT i=0; i<rowOfControlNet; i++)
//	{
//		rg_INT j = 0;
//		for (j=0; j<(columnOfControlNet-1); j++)
//			fprintf( fout, "%12.6f %12.6f %12.6f\n", control_net[i][j].getX(), 
//			               control_net[i][j].getY(), control_net[i][j].getZ() );
//
//		if (i<(rowOfControlNet-1))
//			fprintf( fout, "%12.6f %12.6f %12.6f;\n", control_net[i][j].getX(), 
//			               control_net[i][j].getY(), control_net[i][j].getZ() );
//		else
//			fprintf( fout, "%12.6f %12.6f %12.6f:\n", control_net[i][j].getX(), 
//			               control_net[i][j].getY(), control_net[i][j].getZ() );
//	}
//	
//	fclose( fout );
//}







