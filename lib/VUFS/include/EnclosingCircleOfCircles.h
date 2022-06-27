#include "rg_Circle2D.h"
#include <list>

using namespace std;

class EnclosingCircleOfCircles : public rg_Circle2D
{
private:
	rg_Circle2D**	inputCircles;
	rg_INT		numOfCircles;
	rg_REAL		epsilon;  // (1+epsilon)-approximation : optimal solution ���� ������ epsilon �ȿ� �ִٴ� �ǹ�.
	bool		bComputed;

	rg_REAL*	squaredNorm; //used for filtering when computing farthest circle from a point

public:
	EnclosingCircleOfCircles();
	EnclosingCircleOfCircles(const list<rg_Circle2D*>& circles);
	EnclosingCircleOfCircles(      list<rg_Circle2D>& circles);

	bool	isComputed() const;
	rg_INT	getNumOfCircles() const;

	void setEpsilon(const rg_REAL& e);
	void setCircles(rg_Circle2D** circles, const rg_INT& numCircles);
	void setCircles(rg_Circle2D* circles, const rg_INT& numCircles);
	void setCircles(const list<rg_Circle2D*>& circles);
	void setCircles(      list<rg_Circle2D>& circles);

	bool computeEnclosingCircle();
	bool computeEnclosingCircle(const rg_REAL& eps);
	
	pair<int,rg_REAL> getIndexNMaxDistWithFiltering(const rg_Point2D& query);

};
