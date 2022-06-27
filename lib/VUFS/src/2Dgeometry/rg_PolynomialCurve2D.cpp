#include "rg_PolynomialCurve2D.h"


rg_PolynomialCurve2D::rg_PolynomialCurve2D()
{
    degree=-1;
    coefficients=rg_NULL;
}

rg_PolynomialCurve2D::rg_PolynomialCurve2D(const rg_DEGREE &tDegree)
{
    degree=tDegree;
    coefficients=new rg_Point2D[degree+1];
}

rg_PolynomialCurve2D::rg_PolynomialCurve2D(const rg_DEGREE &tDegree, const rg_Point2D* tCoefficients)
{
    degree=tDegree;
    coefficients=new rg_Point2D[degree+1];
    for( rg_INT i=0; i < degree+1; i++ )
    {
        coefficients[i]=tCoefficients[i];
    }
}

rg_PolynomialCurve2D::rg_PolynomialCurve2D(const rg_PolynomialCurve2D &temp)
{
    degree=temp.degree;
    coefficients=new rg_Point2D[degree+1];
    for( rg_INT i=0; i < degree+1; i++ )
    {
        coefficients[i]=temp.coefficients[i];
    }

}

rg_PolynomialCurve2D::~rg_PolynomialCurve2D()
{
    removeAll();
}

void rg_PolynomialCurve2D::removeAll()
{
    if ( coefficients != rg_NULL )
    {
        delete[] coefficients;
        coefficients=rg_NULL;
    }
    degree=-1;
}

//Operations
rg_Point2D    rg_PolynomialCurve2D::evaluatePt(const rg_PARAMETER& t) const
{
	rg_Point2D result(0.0,0.0);

    for( rg_INDEX i=0; i < degree+1; i++ )
    {
        result = result* t + coefficients[degree-i];
    }

	return result;
}

//Access elements
rg_DEGREE   rg_PolynomialCurve2D::getDegree() const
{
    return degree;
}

rg_Point2D  rg_PolynomialCurve2D::getCoefficient(const rg_DEGREE& i) const
{
    return coefficients[i];
}

rg_Point2D* rg_PolynomialCurve2D::getCoefficients() const
{
    if ( coefficients == rg_NULL )
    {
        return rg_NULL;
    }

    rg_Point2D* output=new rg_Point2D[degree+1];
    for( rg_INT i=0; i < degree+1; i++ )
    {
        output[i]=coefficients[i];
    }
    return output;
}


void     rg_PolynomialCurve2D::setDegree(const rg_DEGREE& tDegree)
{
    if ( degree == tDegree )
    {
        return;
    }

    if ( coefficients != rg_NULL )
    {
        removeAll();
    }
    degree=tDegree;
    coefficients=new rg_Point2D[degree+1];
}

void     rg_PolynomialCurve2D::setCoefficient(const rg_INDEX& i, const rg_Point2D& coefficient)
{
    coefficients[i]=coefficient;
}

void     rg_PolynomialCurve2D::setCoefficients(const rg_Point2D* tCoefficients)
{
    for( rg_INT i=0; i < degree+1; i++ )
    {
        coefficients[i]=tCoefficients[i];
    }
}
void     rg_PolynomialCurve2D::setCurve(const rg_PolynomialCurve2D& temp)
{
    if ( temp.degree != degree )
    {
        if ( coefficients != rg_NULL )
        {
            removeAll();
        }
        degree=temp.degree;
        coefficients=new rg_Point2D[degree+1];
    }
    for( rg_INT i=0; i < degree+1; i++ )
    {
        coefficients[i]=temp.coefficients[i];
    }
}

rg_PolynomialCurve2D& rg_PolynomialCurve2D::operator=(const rg_PolynomialCurve2D& temp)
{
    if ( this != &temp )
    {
        setCurve(temp);
    }

    return *this;
}   


