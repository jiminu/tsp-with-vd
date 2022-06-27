#include "rg_IntervalIntersector.h"
#include "rg_Point2D.h"
#include "rg_ComplexNumber.h"
#include "sortFunc.h"
#include "rg_RelativeOp.h"

#include <time.h>
//#include <fstream>

rg_IntervalIntersector::rg_IntervalIntersector()
{

}

rg_IntervalIntersector::rg_IntervalIntersector(const rg_BzCurve2D &curve1, const rg_BzCurve2D &curve2, const rg_dList<rg_Point2D> &ptlist)
{
	curve_s = curve1;
	curve_t = curve2;
	intersectionPointList = ptlist;
}

rg_IntervalIntersector::~rg_IntervalIntersector()
{
}

void rg_IntervalIntersector::intersectBzCurveVsBzCurve(const rg_BzCurve2D& curve1, const rg_BzCurve2D& curve2, rg_REAL &time)
{
	curve_s = curve1;
	curve_t = curve2;

	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime = clock(); 
	//for(rg_INT i=0; i<20; i++)
	//{

		intersectionPointList.removeAll();
		rg_dList<rg_Box2DWithParameter> box1 = predetermineBox(curve_s);
		rg_dList<rg_Box2DWithParameter> box2 = predetermineBox(curve_t);

		box1.reset();
		do
		{
			box2.reset();

			do
			{
				boxingLoop(box1.getEntity(), box2.getEntity());
				box2.setCurrentNext();
			} while(!box2.isHead());

			box1.setCurrentNext();
		} while(!box1.isHead());

	//}

	clock_t EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
}

rg_dList<rg_Box2DWithParameter> rg_IntervalIntersector::predetermineBox(const rg_BzCurve2D& curve)
{
	rg_dList<rg_Box2DWithParameter>    box;
	rg_dList<rg_REAL> param = getLocalMaxMinPointsOntheXY(curve);

    if(param.getSize()==0)
    {
        rg_Box2DWithParameter tm_box(curve.evaluatePt(0.0),
                   curve.evaluatePt(1.0),
                   0.0,
                   1.0);
  //      changingBox(tm_box);
        box.add(tm_box);

        return box;
    }

	rg_REAL tm_param = 0.0;
	param.reset();
	
	do
	{
		rg_Box2DWithParameter tm_box(curve.evaluatePt(tm_param), 
			    curve.evaluatePt(param.getEntity()),
				tm_param, 
				param.getEntity());
//      changingBox(tm_box);
		box.add(tm_box);
		tm_param = param.getEntity();
		param.setCurrentNext();
	} while(!param.isHead());

	rg_Box2DWithParameter tm_box(curve.evaluatePt(tm_param), 
			    curve.evaluatePt(1.0),
				tm_param, 
				1.0);
//    changingBox(tm_box);
	box.add(tm_box);
	
	return box;
}

void rg_IntervalIntersector::boxingLoop(rg_Box2DWithParameter box_s, rg_Box2DWithParameter box_t)
{
	if( isOverlapping(box_s, box_t) )
	{
		rg_Point2D point;
		if(isTerminalCondition(box_s, box_t, point))
		{
			intersectionPointList.add(point);
		}

		else
		{
			rg_Box2DWithParameter *box_s_sub = splitBox(box_s, curve_s);
			rg_Box2DWithParameter *box_t_sub = splitBox(box_t, curve_t);
			
            boxingLoop(box_s_sub[0], box_t_sub[0]);
			boxingLoop(box_s_sub[0], box_t_sub[1]);
			boxingLoop(box_s_sub[1], box_t_sub[0]);
			boxingLoop(box_s_sub[1], box_t_sub[1]);

            delete[] box_s_sub;
            delete[] box_t_sub;
		}
	}

	else 
	{
		return;
	}
}

// 오버랩이 되는 지를 확인할 때 주의해야 할 점은 
// 시작점과 끝점이 어떤 관계를 가지고 있는 지를 확인해야 한다는 것이다. 
// 따라서 다음의 함수의 argument는 local로 두도록 한다. 
rg_INT rg_IntervalIntersector::isOverlapping(rg_Box2DWithParameter box_s, rg_Box2DWithParameter box_t)
{
	changingBox(box_s);
	changingBox(box_t);
/*
	if(   ( rg_GT(box_s.SP.getX(), box_t.EP.getX(), resNeg12) )
       || ( rg_LT(box_s.EP.getX(), box_t.SP.getX(), resNeg12) )
       || ( rg_GT(box_s.SP.getY(), box_t.EP.getY(), resNeg12) )
       || ( rg_LT(box_s.EP.getY(), box_t.SP.getY(), resNeg12) ) )

	{
		return rg_FALSE;
	}
	else 
	{
		return rg_TRUE;
	}
*/
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

void rg_IntervalIntersector::changingBox(rg_Box2DWithParameter &box)
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

rg_INT rg_IntervalIntersector::isTerminalCondition(const rg_Box2DWithParameter &box_s, 
									     const rg_Box2DWithParameter &box_t, 
										 rg_Point2D &point)
{
	rg_Point2D ptOnCurve_s = curve_s.evaluatePt((box_s.ST+box_s.ET)/2.0);
	rg_Point2D ptOnCurve_t = curve_t.evaluatePt((box_t.ST+box_t.ET)/2.0);

	rg_Point2D distance(ptOnCurve_s - ptOnCurve_t);
	rg_REAL test = distance.magnitude();
 	if(	rg_GT(distance.magnitude(), 1.0e-8) ) 
	{
		return rg_FALSE;
	}

	else
	{
		point = ptOnCurve_s;
		return rg_TRUE;
	}
}

rg_Box2DWithParameter *rg_IntervalIntersector::splitBox(const rg_Box2DWithParameter& box, const rg_BzCurve2D& curve)
{
	rg_Box2DWithParameter *returnBox = new rg_Box2DWithParameter[2];

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

rg_dList<rg_Point2D> rg_IntervalIntersector::getIntersectionPointList()
{
	return intersectionPointList;
}

rg_dList<rg_REAL> rg_IntervalIntersector::getLocalMaxMinPointsOntheXY(const rg_BzCurve2D& curve)
{
	rg_INT      n          = curve.getDegree();
	rg_BzCurve2D  ht         = curve.makeDerivative();

	rg_Polynomial *poly = ht.convertBzCurve2Polynomial();

	rg_dList<rg_REAL> result;
	rg_ComplexNumber *x = poly[0].solve();
	rg_ComplexNumber *y = poly[1].solve();

	rg_INT     numOfTotalRoot = poly[0].getDegree();
	delete[] poly;

	rg_REAL *root = new rg_REAL[2*numOfTotalRoot];

	rg_INT i = 0;
	for(i = 0; i < 2*numOfTotalRoot; i++)
	{
		root[i] = 0.0;
	}

	rg_INT j = 0;
	for(i = 0; i < numOfTotalRoot; i++)
	{
		if(x[i].isPureRealNumber() && rg_BTORexclusive(0.0, x[i].getRealNumber(), 1.0))
		{
			root[j++] = x[i].getRealNumber();
		}

		if(y[i].isPureRealNumber() && rg_BTORexclusive(0.0, y[i].getRealNumber(), 1.0))
		{
			root[j++] = y[i].getRealNumber();
		}
	}

	QuickSort(root, 0, 2*numOfTotalRoot-1);

	for(i = 0; i < 2*numOfTotalRoot; i++)
	{
		if( rg_NZERO(root[i]) )
			result.add(root[i]);
	}
    delete[] root;

	return result;
}


