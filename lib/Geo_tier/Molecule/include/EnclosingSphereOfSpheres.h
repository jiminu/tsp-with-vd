#ifndef _ENCLOSING_SPHERE_H_
#define _ENCLOSING_SPHERE_H_

#include "Sphere.h"
#include "rg_Atom.h"
#include <list>

using namespace std;


namespace V {
namespace GeometryTier {



class EnclosingSphereOfSpheres : public Sphere
{
private:
	Sphere**	inputSpheres;
	rg_INT		numOfSpheres;
	rg_REAL		epsilon;  // (1+epsilon)-approximation : optimal solution 과의 오차가 epsilon 안에 있다는 의미.
	bool		bComputed;

	rg_REAL*	squaredNorm; //used for filtering when computing farthest sphere from a point

public:
	EnclosingSphereOfSpheres();
	EnclosingSphereOfSpheres(const EnclosingSphereOfSpheres& enclosingSphereOfSpheres);
	EnclosingSphereOfSpheres(const list<Sphere*>& spheres);
	EnclosingSphereOfSpheres(      list<Sphere>& spheres);
	EnclosingSphereOfSpheres(const list<Atom*>& spheres);
	EnclosingSphereOfSpheres(      list<Atom>& spheres);

	virtual ~EnclosingSphereOfSpheres();

	bool	isComputed() const;
	rg_INT	getNumOfSpheres() const;

	void setEpsilon(const rg_REAL& e);
	void setSpheres(Sphere** spheres, const rg_INT& numSpheres);
	void setSpheres(Sphere* spheres, const rg_INT& numSpheres);
	void setSpheres(const list<Sphere*>& spheres);
	void setSpheres(      list<Sphere>& spheres);
	void setSpheres(      rg_dList<Sphere>& spheres);
	void setSpheres(const list<Atom*>& spheres);
	void setSpheres(      list<Atom>& spheres);
	void setSpheres(const rg_dList<Atom*>& spheres);

	bool computeEnclosingSphere();
	bool computeEnclosingSphere(const rg_REAL& eps);
	
	pair<int,rg_REAL> getIndexNMaxDistWithFiltering(const rg_Point3D& query);

};



} // namespace GeometryTier
} // namespace V


#endif

