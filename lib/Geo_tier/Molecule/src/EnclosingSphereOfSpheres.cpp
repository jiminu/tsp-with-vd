#include "EnclosingSphereOfSpheres.h"
using namespace V::GeometryTier;


EnclosingSphereOfSpheres::EnclosingSphereOfSpheres()
{
	inputSpheres = NULL;
	squaredNorm = NULL;
	numOfSpheres = 0;
	epsilon = 0.001;
	bComputed = false;
}

EnclosingSphereOfSpheres::EnclosingSphereOfSpheres(const EnclosingSphereOfSpheres& enclosingSphereOfSpheres)
{

}

EnclosingSphereOfSpheres::EnclosingSphereOfSpheres(const list<Sphere*>& spheres)
{
	list<Sphere*>::const_iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	inputSpheres = new Sphere* [numOfSpheres];
	squaredNorm = new rg_REAL [numOfSpheres];

	for(int i=0; i<numOfSpheres; i++)
	{
		inputSpheres[i] = *iter;
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();

		++iter;
	}

	epsilon = 0.001;
	bComputed = false;
}

EnclosingSphereOfSpheres::EnclosingSphereOfSpheres(list<Sphere>& spheres)
{
	list<Sphere>::iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	inputSpheres = new Sphere* [numOfSpheres];
	squaredNorm = new rg_REAL [numOfSpheres];

	for(int i=0; i<numOfSpheres; i++)
	{
		inputSpheres[i] = &(*iter);
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();
		++iter;
	}

	epsilon = 0.001;
	bComputed = false;
}

EnclosingSphereOfSpheres::EnclosingSphereOfSpheres(const list<Atom*>& spheres)
{
	list<Atom*>::const_iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	inputSpheres = new Sphere* [numOfSpheres];
	squaredNorm = new rg_REAL [numOfSpheres];

	Atom* currAtom = NULL;
	for(int i=0; i<numOfSpheres; i++)
	{
		currAtom = *iter;
		inputSpheres[i] = new Sphere(currAtom->getAtomBall());
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();
		++iter;
	}

	epsilon = 0.001;
	bComputed = false;
}

EnclosingSphereOfSpheres::EnclosingSphereOfSpheres(list<Atom>& spheres)
{
	list<Atom>::iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	inputSpheres = new Sphere* [numOfSpheres];
	squaredNorm = new rg_REAL [numOfSpheres];

	Atom* currAtom = NULL;
	for(int i=0; i<numOfSpheres; i++)
	{
		currAtom = &(*iter);
		inputSpheres[i] = new Sphere(currAtom->getAtomBall());
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();

		++iter;
	}

	epsilon = 0.001;
	bComputed = false;
}

EnclosingSphereOfSpheres::~EnclosingSphereOfSpheres()
{
	if(inputSpheres != NULL)
	{
		delete [] inputSpheres;
	}
	if(squaredNorm != NULL)
	{
		delete [] squaredNorm;
	}

	inputSpheres = NULL;
	squaredNorm = NULL;
	numOfSpheres = 0;
	epsilon = 0.001;
	bComputed = false;
}

bool	EnclosingSphereOfSpheres::isComputed() const
{
	return bComputed;
}

rg_INT	EnclosingSphereOfSpheres::getNumOfSpheres() const
{
	return numOfSpheres;
}


void EnclosingSphereOfSpheres::setEpsilon(const rg_REAL& e)
{
	epsilon = e;
}

void EnclosingSphereOfSpheres::setSpheres(Sphere** spheres, const rg_INT& numSpheres)
{
	inputSpheres = spheres;
	numOfSpheres = numSpheres;
}

void EnclosingSphereOfSpheres::setSpheres(Sphere* spheres, const rg_INT& numSpheres)
{
	numOfSpheres = numSpheres;
	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];

	if(squaredNorm != NULL)
		delete [] squaredNorm;

	for(int i=0; i<numOfSpheres; i++)
	{
		inputSpheres[i] = &spheres[i];
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();
	}
}

void EnclosingSphereOfSpheres::setSpheres(const list<Sphere*>& spheres)
{
	list<Sphere*>::const_iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfSpheres];

	for(int i=0; i<numOfSpheres; i++)
	{
		inputSpheres[i] = *iter;
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();

		++iter;
	}
}

void EnclosingSphereOfSpheres::setSpheres(list<Sphere>& spheres)
{
	list<Sphere>::iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfSpheres];

	for(int i=0; i<numOfSpheres; i++)
	{
		inputSpheres[i] = &(*iter);
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();

		++iter;
	}
}

void EnclosingSphereOfSpheres::setSpheres(rg_dList<Sphere>& spheres)
{
	numOfSpheres = spheres.getSize();
	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];

	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfSpheres];

	rg_INT i = 0;
	spheres.reset4Loop();
	while(spheres.setNext4Loop())
	{
		inputSpheres[i] = spheres.getpEntity();
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();
		i++;
	}
}

void EnclosingSphereOfSpheres::setSpheres(const list<Atom*>& spheres)
{
	list<Atom*>::const_iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfSpheres];

	Atom* currAtom = NULL;
	for(int i=0; i<numOfSpheres; i++)
	{
		currAtom = *iter;
		inputSpheres[i] = new Sphere(currAtom->getAtomBall());
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();

		++iter;
	}
}

void EnclosingSphereOfSpheres::setSpheres(list<Atom>& spheres)
{
	list<Atom>::iterator iter = spheres.begin();
	numOfSpheres = spheres.size();
	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfSpheres];

	Atom* currAtom = NULL;
	for(int i=0; i<numOfSpheres; i++)
	{
		currAtom = &(*iter);
		inputSpheres[i] = new Sphere(currAtom->getAtomBall());
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();

		++iter;
	}
}

void EnclosingSphereOfSpheres::setSpheres(const rg_dList<Atom*>& spheres)
{
	rg_INT size = spheres.getSize();
	if(size == 0)
		return;

	numOfSpheres = size;

	if(inputSpheres != NULL)
		delete [] inputSpheres;
	inputSpheres = new Sphere* [numOfSpheres];
	
	if(squaredNorm != NULL)
		delete [] squaredNorm;
	squaredNorm = new rg_REAL [numOfSpheres];

	rg_INDEX i = 0;
	spheres.reset4Loop();
	while(spheres.setNext4Loop())
	{
		Atom* currAtom = spheres.getEntity();
		inputSpheres[i] = currAtom->getpAtomBall();
		squaredNorm[i] = inputSpheres[i]->getCenter().squaredMagnitude();
		i++;
	}	
}

bool EnclosingSphereOfSpheres::computeEnclosingSphere()
{
	if(bComputed)
		return true;
	if( numOfSpheres < 1)
		return false;
	if( numOfSpheres == 1)
	{
		this->setSphere(inputSpheres[0]->getCenter(), inputSpheres[0]->getRadius());
		return true;
	}

	rg_Point3D currCenter = inputSpheres[0]->getCenter();

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

			//move current sphere until it touches farthest point pts[m]
			currCenter = currCenter+(inputSpheres[m]->getCenter()-currCenter).getUnitVector()*(maxDist-r);
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

	this->setSphere(currCenter,r+delta);
	bComputed = true;

	return true;
}

bool EnclosingSphereOfSpheres::computeEnclosingSphere(const rg_REAL& eps)
{
	bComputed = false;
	epsilon = eps;
	return computeEnclosingSphere();
}

pair<int,rg_REAL> EnclosingSphereOfSpheres::getIndexNMaxDistWithFiltering(const rg_Point3D& query)
{
	int maxIndex = -1;
	rg_REAL maxDist = -1;

	rg_REAL currSqrdNorm = query.squaredMagnitude();
	rg_REAL upperBound = 0;

	rg_REAL currDist = 0;

	for(int i=0; i<numOfSpheres; i++)
	{
		upperBound = currSqrdNorm + squaredNorm[i]+2*sqrt(currSqrdNorm*squaredNorm[i]) + inputSpheres[i]->getRadius();
		if(maxDist > upperBound)
			continue;

		currDist = query.distance(inputSpheres[i]->getCenter())+inputSpheres[i]->getRadius();

		if( currDist > maxDist )
		{
			maxIndex = i;
			maxDist = currDist;
		}
	}

	return pair<int,rg_REAL>(maxIndex, maxDist);

}
