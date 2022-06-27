#include "EnclosingCircleOfCircles.h"


EnclosingCircleOfCircles::EnclosingCircleOfCircles()
{
	inputCircles = NULL;
	squaredNorm = NULL;
	numOfCircles = 0;
	epsilon = 0.001;
	bComputed = false;
}

EnclosingCircleOfCircles::EnclosingCircleOfCircles(const list<rg_Circle2D*>& circles)
{
	list<rg_Circle2D*>::const_iterator iter = circles.begin();
	numOfCircles = circles.size();
	inputCircles = new rg_Circle2D* [numOfCircles];
	squaredNorm = new rg_REAL [numOfCircles];

	for(int i=0; i<numOfCircles; i++)
	{
		inputCircles[i] = *iter;
		squaredNorm[i] = inputCircles[i]->getCenterPt().magnitudeSquare();

		++iter;
	}

	epsilon = 0.001;
	bComputed = false;
}

EnclosingCircleOfCircles::EnclosingCircleOfCircles(list<rg_Circle2D>& circles)
{
	list<rg_Circle2D>::iterator iter = circles.begin();
	numOfCircles = circles.size();
	inputCircles = new rg_Circle2D* [numOfCircles];
	squaredNorm = new rg_REAL [numOfCircles];

	for(int i=0; i<numOfCircles; i++)
	{
		inputCircles[i] = &(*iter);
		squaredNorm[i] = inputCircles[i]->getCenterPt().magnitudeSquare();
		++iter;
	}

	epsilon = 0.001;
	bComputed = false;
}


bool	EnclosingCircleOfCircles::isComputed() const
{
	return bComputed;
}

rg_INT	EnclosingCircleOfCircles::getNumOfCircles() const
{
	return numOfCircles;
}


void EnclosingCircleOfCircles::setEpsilon(const rg_REAL& e)
{
	epsilon = e;
}

void EnclosingCircleOfCircles::setCircles(rg_Circle2D** circles, const rg_INT& numCircles)
{
	inputCircles = circles;
	numOfCircles = numCircles;
}

void EnclosingCircleOfCircles::setCircles(rg_Circle2D* circles, const rg_INT& numCircles)
{
	numOfCircles = numCircles;
	if(inputCircles != NULL)
		delete [] inputCircles;
	inputCircles = new rg_Circle2D* [numOfCircles];

	if(squaredNorm != NULL)
		delete [] squaredNorm;

	for(int i=0; i<numOfCircles; i++)
	{
		inputCircles[i] = &circles[i];
		squaredNorm[i] = inputCircles[i]->getCenterPt().magnitudeSquare();
	}
}

void EnclosingCircleOfCircles::setCircles(const list<rg_Circle2D*>& circles)
{
	list<rg_Circle2D*>::const_iterator iter = circles.begin();
	numOfCircles = circles.size();
	if(inputCircles != NULL)
		delete [] inputCircles;
	inputCircles = new rg_Circle2D* [numOfCircles];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfCircles];

	for(int i=0; i<numOfCircles; i++)
	{
		inputCircles[i] = *iter;
		squaredNorm[i] = inputCircles[i]->getCenterPt().magnitudeSquare();

		++iter;
	}
}

void EnclosingCircleOfCircles::setCircles(list<rg_Circle2D>& circles)
{
	list<rg_Circle2D>::iterator iter = circles.begin();
	numOfCircles = circles.size();
	if(inputCircles != NULL)
		delete [] inputCircles;
	inputCircles = new rg_Circle2D* [numOfCircles];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfCircles];

	for(int i=0; i<numOfCircles; i++)
	{
		inputCircles[i] = &(*iter);
		squaredNorm[i] = inputCircles[i]->getCenterPt().magnitudeSquare();

		++iter;
	}
}


bool EnclosingCircleOfCircles::computeEnclosingCircle()
{
	if(bComputed)
		return true;
	if( numOfCircles < 1)
		return false;
	if( numOfCircles == 1)
	{
		this->setCircle(inputCircles[0]->getCenterPt(), inputCircles[0]->getRadius());
		return true;
	}

	rg_Point2D currCenter = inputCircles[0]->getCenterPt();

	rg_REAL r = getIndexNMaxDistWithFiltering(currCenter).second/2.;
	rg_REAL delta = r;
	int numIter = 0;
	while(delta > epsilon)
	{
		numIter = (int)(2); // O(1/delta) iterations
		for(int i = 0; i<numIter; i++)
		{
			pair<int,rg_REAL> currMaxPair = getIndexNMaxDistWithFiltering(currCenter);

			int m = currMaxPair.first;
			rg_REAL maxDist = currMaxPair.second;

			//move current circle until it touches farthest point pts[m]
			currCenter = currCenter+(inputCircles[m]->getCenterPt()-currCenter).getUnitVector()*(maxDist-r);
		}
		rg_REAL s = getIndexNMaxDistWithFiltering(currCenter).second - r;

		rg_REAL newDelta = delta * 3./4.;
		if(s <= newDelta)
		{
			delta = newDelta;
		}
		else
		{
			r = r + delta/4.;
			delta = newDelta;
		}
	}

	this->setCircle(currCenter,r+delta);
	bComputed = true;

	return true;
}

bool EnclosingCircleOfCircles::computeEnclosingCircle(const rg_REAL& eps)
{
	bComputed = false;
	epsilon = eps;
	return computeEnclosingCircle();
}

pair<int,rg_REAL> EnclosingCircleOfCircles::getIndexNMaxDistWithFiltering(const rg_Point2D& query)
{
	int maxIndex = -1;
	rg_REAL maxDist = -1;

	rg_REAL currSqrdNorm = query.magnitudeSquare();
	rg_REAL upperBound = 0;

	rg_REAL currDist = 0;

	for(int i=0; i<numOfCircles; i++)
	{
		upperBound = currSqrdNorm + squaredNorm[i]+2*sqrt(currSqrdNorm*squaredNorm[i]) + inputCircles[i]->getRadius();
		if(maxDist > upperBound)
			continue;

		currDist = query.distance(inputCircles[i]->getCenterPt())+inputCircles[i]->getRadius();

		if( currDist > maxDist )
		{
			maxIndex = i;
			maxDist = currDist;
		}
	}

	return pair<int,rg_REAL>(maxIndex, maxDist);

}
