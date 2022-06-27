#include <math.h>
//#include <fstream>
#include <time.h>

#include "rg_BzIntersector.h"
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"
#include "rg_IntersectFunc.h"
#include "sortFunc.h"

#include "rg_Point3D.h" // test

const rg_REAL pi = 4.0*atan(1.0);

rg_BzIntersector::rg_BzIntersector()
{

}
 
rg_BzIntersector::~rg_BzIntersector()
{

}

rg_dList<rg_Point2D> rg_BzIntersector::intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, 
														rg_dList<rg_REAL*> &seedParam4TwoCurve)
{

    rg_dList<rg_Point2D> intersectPointList;

	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime, EndTime;
    StartTime = clock();


	//for(rg_INT i=0; i<20; i++)
	//{

        intersectPointList.removeAll();
  
        rg_dList<rg_REAL>    param_s = findCharParam(curve_s);
        rg_dList<rg_REAL>    param_t = findCharParam(curve_t);

        rg_dList<rg_RQBzCurve2D> rqBzCurveList_s    = approximateRQBzCurves(curve_s, param_s);
        rg_dList<rg_RQBzCurve2D> rqBzCurveList_t    = approximateRQBzCurves(curve_t, param_t);

	    seedParam4TwoCurve = makeSeed(rqBzCurveList_s, param_s, 
								      rqBzCurveList_t, param_t);

//    rg_dList<rg_Point2D> intersectPointList;
        if(seedParam4TwoCurve.getSize()==0) return intersectPointList;

        seedParam4TwoCurve.reset();

	    do
	    {
		    rg_REAL *parameter = seedParam4TwoCurve.getEntity();
		    rg_Point2D intersectPt = iterationWithSeed(curve_s, parameter[0], 
			                                        curve_t, parameter[1]);
		    
		    intersectPointList.add(intersectPt);
		    seedParam4TwoCurve.setCurrentNext();
	    } while(!seedParam4TwoCurve.isHead());

	//}

	EndTime = clock();
//	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
//    ofstream fout("timeOut.dat", ios::app);
//    fout << time << endl;

	return intersectPointList;
}








rg_dList<rg_Point2D> rg_BzIntersector::intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, 
														rg_dList<rg_REAL*> &seedParam4TwoCurve, rg_REAL& time)
{

    rg_dList<rg_Point2D> intersectPointList;

	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime, EndTime;
    StartTime = clock();


	//for(rg_INT i=0; i<20; i++)
	//{

        intersectPointList.removeAll();
  
        rg_dList<rg_REAL>    param_s = findCharParam(curve_s);
        rg_dList<rg_REAL>    param_t = findCharParam(curve_t);

        rg_dList<rg_RQBzCurve2D> rqBzCurveList_s    = approximateRQBzCurves(curve_s, param_s);
        rg_dList<rg_RQBzCurve2D> rqBzCurveList_t    = approximateRQBzCurves(curve_t, param_t);

	    seedParam4TwoCurve = makeSeed(rqBzCurveList_s, param_s, 
								      rqBzCurveList_t, param_t);

//    rg_dList<rg_Point2D> intersectPointList;
        if(seedParam4TwoCurve.getSize()==0) return intersectPointList;

        seedParam4TwoCurve.reset();

	    do
	    {
		    rg_REAL *parameter = seedParam4TwoCurve.getEntity();
		    rg_Point2D intersectPt = iterationWithSeed(curve_s, parameter[0], 
			                                        curve_t, parameter[1]);
		    
		    intersectPointList.add(intersectPt);
		    seedParam4TwoCurve.setCurrentNext();
	    } while(!seedParam4TwoCurve.isHead());

	//}

	EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
//    ofstream fout("timeOut.dat", ios::app);
//    fout << time << endl;

	return intersectPointList;
}











/*
rg_dList<rg_Point2D> rg_BzIntersector::intersectBzCurveVsBzCurve(rg_BzCurve2D &curve_s, rg_BzCurve2D &curve_t, 
														rg_dList<rg_REAL*> &seedParam4TwoCurve)
{

    rg_dList<rg_Point2D> intersectPointList;

	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime = clock(); 
	for(rg_INT i=0; i<10; i++)
	{

    rg_dList<rg_REAL>    param_s = findCharParam(curve_s);
	rg_dList<rg_REAL>    param_t = findCharParam(curve_t);
	rg_dList<rg_RQBzCurve2D> rqBzCurveList_s    = approximateRQBzCurves(curve_s, param_s);
	rg_dList<rg_RQBzCurve2D> rqBzCurveList_t    = approximateRQBzCurves(curve_t, param_t);
	seedParam4TwoCurve = makeSeed(rqBzCurveList_s, param_s, 
								  rqBzCurveList_t, param_t);

//    rg_dList<rg_Point2D> intersectPointList;
    if(seedParam4TwoCurve.getSize()==0) return intersectPointList;

	seedParam4TwoCurve.reset();

	do
	{
		rg_REAL *parameter = seedParam4TwoCurve.getEntity();
		rg_Point2D intersectPt = iterationWithSeed(curve_s, parameter[0], 
			                                    curve_t, parameter[1]);
		
		intersectPointList.add(intersectPt);
		seedParam4TwoCurve.setCurrentNext();
	} while(!seedParam4TwoCurve.isHead());

	}

	clock_t EndTime = clock();
	rg_REAL ComputingTime = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	ofstream fout("out.dat", ios::app);
	fout << " ---- Our Alg. ----" << endl;
    fout << " degree of C1 : " << curve_s.getDegree();
    fout << ",  degree of C2 : " << curve_t.getDegree() << endl;
	fout << " # of intersection points  : " << intersectPointList.getSize() << endl;
	fout << " computation time/10.0  : " << ComputingTime/10 << endl ;

	return intersectPointList;
}
*/

rg_dList<rg_RQBzCurve2D> rg_BzIntersector::approximateRQBzCurves(const rg_BzCurve2D &curve, rg_dList<rg_REAL> &param)
{
	rg_dList<rg_RQBzCurve2D> rqBzCurveList;

	param.reset();
	rg_REAL t0 = param.getEntity();
	param.setCurrentNext();
	
	do
	{
		rg_REAL t1 = param.getEntity();
		rqBzCurveList.add(makeOneRQBzCurve(curve, t0, t1));
		t0 = t1;
		param.setCurrentNext();
	} while(!param.isHead());

	return rqBzCurveList;
}

rg_RQBzCurve2D rg_BzIntersector::makeOneRQBzCurve(const rg_BzCurve2D &curve, 
										  const rg_REAL &t0,
										  const rg_REAL &t1)
{
	rg_BzCurve2D  curveDev = curve.makeDerivative();
	rg_Point2D pt[2];
	rg_Point2D pt_prime[2];

	pt[0]       = curve.evaluatePt(t0);
	pt_prime[0] =  curveDev.evaluatePt(t0);

	pt[1]       = curve.evaluatePt(t1);
	pt_prime[1] =  curveDev.evaluatePt(t1);

	rg_RQBzCurve2D rqcurve;
	rqcurve.makeRQBezier(pt[0], pt_prime[0], 
		                 pt[1], pt_prime[1],
						 rg_Point2D(curve.evaluatePt((t0+t1)/2.0)));

	return rqcurve;
}

rg_dList<rg_REAL*> rg_BzIntersector::makeSeed(rg_dList<rg_RQBzCurve2D> &rqcurve_s, rg_dList<rg_REAL> subParam_s,
	                                        rg_dList<rg_RQBzCurve2D> &rqcurve_t, rg_dList<rg_REAL> subParam_t)
{
	rg_dList <rg_REAL*> seedParam4TwoCurve;
	rg_dList <rg_Point2D> intersectPointList;

	rqcurve_s.reset();
	subParam_s.reset();
	rg_REAL s0 = subParam_s.getEntity();
	subParam_s.setCurrentNext();

	do
	{
		rg_REAL s1 = subParam_s.getEntity();

		rqcurve_t.reset();
		subParam_t.reset();
		rg_REAL t0 = subParam_t.getEntity();
		subParam_t.setCurrentNext();
		
		do
		{
			rg_REAL t1 = subParam_t.getEntity();
	
			// Using the bounding box, 
            // it can be omitted to  compute the intersetions between two curves.
            rg_RQBzCurve2D   rqbz_s = rqcurve_s.getEntity();
            rg_RQBzCurve2D   rqbz_t = rqcurve_t.getEntity();
            rg_BoundingBox2D box_s  = rqbz_s.makeBoundingBox();
            rg_BoundingBox2D box_t  = rqbz_t.makeBoundingBox();
			if( box_s.isOverlapped(box_t) )
			{
				rg_IntersectFunc::intersectRQBzCurveVsRQBzCurve(rqbz_s, s0, s1, 
											  rqbz_t, t0, t1, 
											  intersectPointList, 
											  seedParam4TwoCurve);
			}
			t0 = t1;
			rqcurve_t.setCurrentNext();
			subParam_t.setCurrentNext();
		} while(!rqcurve_t.isHead());

		s0 = s1;
		rqcurve_s.setCurrentNext();
		subParam_s.setCurrentNext();
	} while(!rqcurve_s.isHead());

	return seedParam4TwoCurve;
}

rg_Point2D rg_BzIntersector::iterationWithSeed(const rg_BzCurve2D &curve_s, const rg_REAL &param_s,
	                                     const rg_BzCurve2D &curve_t, const rg_REAL &param_t)
{
/*
	rg_REAL newParam_s = param_s;
	rg_REAL newParam_t = param_t;

	rg_Point2D initialPtOnCurve_s = curve_s.evaluatePt(newParam_s);
	rg_Point2D initialPtOnCurve_t = curve_t.evaluatePt(newParam_t);

	rg_Point2D ptOnCurve_s = initialPtOnCurve_s;
	rg_Point2D ptOnCurve_t = initialPtOnCurve_t;

	rg_Point2D distance = ptOnCurve_s - ptOnCurve_t;

	rg_BzCurve2D devCurve_s = curve_s.makeDerivative();
	rg_BzCurve2D devCurve_t = curve_t.makeDerivative();

	rg_INT count = 0;
	while( rg_GT(distance.magnitudeSquare(), 1.0e-12) && (count < 8) )
	{
		rg_Point2D tangentVector_s = devCurve_s.evaluatePt(newParam_s);
		rg_Point2D tangentVector_t = devCurve_t.evaluatePt(newParam_t);

		rg_REAL temp = tangentVector_s*tangentVector_t;
		newParam_s += (ptOnCurve_t-ptOnCurve_s)*tangentVector_t/temp;
		newParam_t += (ptOnCurve_s-ptOnCurve_t)*tangentVector_s/(-temp);
		
		ptOnCurve_s.setPoint( curve_s.evaluatePt(newParam_s) );
		ptOnCurve_t.setPoint( curve_t.evaluatePt(newParam_t) );

		distance = ptOnCurve_s - ptOnCurve_t;
		count++;
	}

	return (ptOnCurve_s + ptOnCurve_t)/2;
*/

	rg_REAL newParam_s = param_s;
	rg_REAL newParam_t = param_t;

	rg_Point2D initialPtOnCurve_s = curve_s.evaluatePt(newParam_s);
	rg_Point2D initialPtOnCurve_t = curve_t.evaluatePt(newParam_t);

	rg_Point2D ptOnCurve_s = initialPtOnCurve_s;
	rg_Point2D ptOnCurve_t = initialPtOnCurve_t;

	rg_Point2D distance = ptOnCurve_s - ptOnCurve_t;

	rg_BzCurve2D devCurve_s = curve_s.makeDerivative();
	rg_BzCurve2D devCurve_t = curve_t.makeDerivative();

	rg_INT count = 0;
	while( rg_GT(distance.magnitudeSquare(), 1.0e-12) && (count < 8) )
	{
		rg_Point2D tangentVector_s = devCurve_s.evaluatePt(newParam_s);
		rg_Point2D tangentVector_t = devCurve_t.evaluatePt(newParam_t);

		rg_REAL temp = tangentVector_s*tangentVector_t;
		newParam_s += (ptOnCurve_t-ptOnCurve_s)*tangentVector_t/temp;
		newParam_t += (ptOnCurve_s-ptOnCurve_t)*tangentVector_s/(-temp);
		
		ptOnCurve_s.setPoint( curve_s.evaluatePt(newParam_s) );
		ptOnCurve_t.setPoint( curve_t.evaluatePt(newParam_t) );

		distance = ptOnCurve_s - ptOnCurve_t;
		count++;
	}

	return (ptOnCurve_s + ptOnCurve_t)/2;

}

rg_dList<rg_REAL> rg_BzIntersector::findCharParam(const rg_BzCurve2D &curve)
{
	rg_BzCurve2D ht     = curve.makeDerivative();
	rg_BzCurve2D ht_dev =    ht.makeDerivative();

	rg_Polynomial *poly_ht     =     ht.convertBzCurve2Polynomial();
	rg_Polynomial *poly_ht_dev = ht_dev.convertBzCurve2Polynomial();
	rg_Polynomial  polynomial  = poly_ht[0]*poly_ht_dev[1] - poly_ht[1]*poly_ht_dev[0];

    delete[] poly_ht;
    delete[] poly_ht_dev;

	rg_ComplexNumber *root = polynomial.solve();

    rg_INT     numOfTotalRoot = polynomial.getDegree();
	rg_REAL *r = new rg_REAL[numOfTotalRoot];
	rg_INT i = 0;
	for(i = 0; i < numOfTotalRoot; i++) r[i] = 0.0;

	rg_INT j = 0;
	for(i = 0; i < numOfTotalRoot; i++)
	{
		if(root[i].isPureRealNumber() && rg_BTORexclusive(0.0, root[i].getRealNumber(), 1.0))
		{
   			r[j++] = root[i].getRealNumber();
        }
	}

    delete[] root;
	QuickSort(r, 0, numOfTotalRoot-1);

    rg_dList<rg_REAL> inflectParam;
	for(i = 0; i < numOfTotalRoot; i++)
	{
		if( rg_NZERO(r[i]) )
			inflectParam.add(r[i]);
	}
    delete[] r;

    inflectParam.addHead(0.0);
	inflectParam.add(1.0);
    inflectParam.reset();

    rg_dList<rg_REAL> charParam;
   
    rg_REAL t0 = inflectParam.getEntity();
    charParam.add(t0);
    inflectParam.setCurrentNext();

    do
    {
        rg_REAL t1 = inflectParam.getEntity();
        
        makeSimpleParam(ht, t0, t1, charParam);
        charParam.add(t1);
        t0 = t1;
        inflectParam.setCurrentNext();
    } while(!inflectParam.isHead());

	return charParam;
}

void rg_BzIntersector::makeSimpleParam(const rg_BzCurve2D& curve, const rg_REAL& t0, const rg_REAL& t1,
                                    rg_dList<rg_REAL>& charParam)
{
    rg_Point2D v1 = curve.evaluatePt(t0);
    rg_Point2D v2 = curve.evaluatePt(t1);
    rg_Point2D p  = curve.evaluatePt((t0+t1)/2.0);

    rg_REAL  D  = v1.getX()*v2.getY() - v2.getX()*v1.getY();

    rg_REAL  a  = (v2.getY()*p.getX() - v2.getX()*p.getY())/D;
    rg_REAL  b  = (v1.getX()*p.getY() - v1.getY()*p.getX())/D;

    rg_REAL  theta = acos(v1%v2/v1.magnitude()/v2.magnitude());

    if(rg_POS(a) && rg_POS(b))
    {
        return;
    }

    else
    {
        theta = 2*pi - theta;
        
        rg_REAL *t = getConjugateTangentParam(curve, t0, t1);

        charParam.add(t[0]);
        charParam.add(t[1]);
            
        delete[] t;
    }
}

rg_REAL* rg_BzIntersector::getConjugateTangentParam(const rg_BzCurve2D& curve, 
                                                const rg_REAL& t0,
                                                const rg_REAL& t1)
{
    rg_Point2D v0 = curve.evaluatePt(t0);
    rg_Point2D v1 = curve.evaluatePt(t1);
 
    rg_REAL tangent0 = v0.getY()/v0.getX();
    rg_REAL tangent1 = v1.getY()/v1.getX();

    rg_Polynomial *curvePolynomial = curve.convertBzCurve2Polynomial();
    rg_Polynomial p0 = curvePolynomial[1] - tangent0*curvePolynomial[0];
    rg_Polynomial p1 = curvePolynomial[1] - tangent1*curvePolynomial[0];

    delete[] curvePolynomial;

    rg_ComplexNumber *root0 = p0.solve();
    rg_ComplexNumber *root1 = p1.solve();

    rg_REAL *result = new rg_REAL[2];
    result[0] = t0;
    result[1] = t1;
    
    const rg_REAL errBound = 0.1;
        
    for(rg_INT i=0; i<p0.getDegree(); i++)
    {
        if( root0[i].isPureRealNumber() )
        {
            rg_REAL r = root0[i].getRealNumber();
            if(rg_BTORexclusive(result[0], r, t1) ) 
            {
                result[1] = r;
            }
        }

        if( root1[i].isPureRealNumber() )
        {
            rg_REAL r = root1[i].getRealNumber();
            if(rg_BTORexclusive(t0, r, result[1]) )
            {
                result[0] = r;
            }
        }
    }

    delete[] root0;
    delete[] root1;
    return result;
}


