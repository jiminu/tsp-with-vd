#include "rg_RPolynomialCurve2D.h"

void rg_RPolynomialCurve2D::removeAll()
{
    if ( coefficientsOfDenominator != rg_NULL )
    {
        delete[] coefficientsOfDenominator;
    }

    rg_PolynomialCurve2D::removeAll();


}


rg_RPolynomialCurve2D::rg_RPolynomialCurve2D()
:rg_PolynomialCurve2D(), coefficientsOfDenominator(rg_NULL)
{
}

rg_RPolynomialCurve2D::rg_RPolynomialCurve2D(const rg_DEGREE &tDegree)
:rg_PolynomialCurve2D(tDegree)
{
    coefficientsOfDenominator=new rg_REAL[tDegree+1];
    coefficientsOfDenominator[0]=1.0;
    for( rg_INT i=1 ; i < tDegree+1; i++)
    {
        coefficientsOfDenominator[i]=0.0;

    }
}

rg_RPolynomialCurve2D::rg_RPolynomialCurve2D(const rg_PolynomialCurve2D &tCurve)
:rg_PolynomialCurve2D(tCurve)
{
    const rg_INT degree=rg_PolynomialCurve2D::getDegree();
    coefficientsOfDenominator=new rg_REAL[degree+1];
    coefficientsOfDenominator[0]=1.0;
    for( rg_INT i=1 ; i < degree+1; i++)
    {
        coefficientsOfDenominator[i]=0.0;

    }

}

rg_RPolynomialCurve2D::rg_RPolynomialCurve2D(const rg_RPolynomialCurve2D &temp)
:rg_PolynomialCurve2D(temp)
{
    const rg_INT degree=rg_PolynomialCurve2D::getDegree();
    coefficientsOfDenominator=new rg_REAL[degree+1];
    for( rg_INT i=0 ; i < degree+1; i++)
    {
        coefficientsOfDenominator[i]=temp.coefficientsOfDenominator[i];

    }
}

rg_RPolynomialCurve2D::~rg_RPolynomialCurve2D()
{
    removeAll();
}


//Operations
rg_Point2D    rg_RPolynomialCurve2D::evaluatePt(const rg_PARAMETER& t) const
{
    rg_Point2D output=rg_PolynomialCurve2D::evaluatePt(t);
	rg_REAL denominator=0.0;
    const rg_INT degree=rg_PolynomialCurve2D::getDegree();
    for( rg_INDEX i=0; i < degree+1; i++ )
    {
        denominator = denominator* t + coefficientsOfDenominator[degree-i];
    }                              
    
    return output/denominator;
}


rg_RPolynomialCurve2D rg_RPolynomialCurve2D::getDerivative() const
{
    const rg_INT degree=rg_PolynomialCurve2D::getDegree();
    rg_Point2D* dNumerator=new rg_Point2D [degree];
    rg_REAL*    dDenominator=new rg_REAL[degree];
    
    rg_INT i = 0;
    for( i=0; i < degree; i++ )
    {
        dNumerator[i]=((rg_REAL)(i+1))*rg_PolynomialCurve2D::getCoefficient(i+1);
        dDenominator[i]=((rg_REAL)(i+1))*getCoefficientOfDenominator(i+1);
    }

    rg_Point2D* NumeratorOfDerivative=new rg_Point2D[2*degree+1];
    rg_REAL*    DenominatorOfDerivative=new rg_REAL[2*degree+1];
    for ( i=0; i < 2*degree+1; i++ )
    {
        NumeratorOfDerivative[i].setPoint(0.0,0.0);
        DenominatorOfDerivative[i]=0.0;
    }

    for ( i=0; i < degree+1; i++ )
    {
	rg_INT j = 0;
        for( j=0; j < degree; j++)
        {
            NumeratorOfDerivative[i+j] += (  dNumerator[j]*getCoefficientOfDenominator(i)
                                           - dDenominator[j]*getCoefficient(i) );
            DenominatorOfDerivative[i+j] +=getCoefficientOfDenominator(i)*getCoefficientOfDenominator(j);
        }
        DenominatorOfDerivative[i+j]+=getCoefficientOfDenominator(i)*getCoefficientOfDenominator(j);
    }
    rg_RPolynomialCurve2D output(2*degree);
    output.setCoefficients(NumeratorOfDerivative);
    output.setCoefficientsOfDenominator(DenominatorOfDerivative);

    delete[] dNumerator;
    delete[] dDenominator;
    delete[] NumeratorOfDerivative;
    delete[] DenominatorOfDerivative;

    return output;


}
//Access elements
rg_REAL     rg_RPolynomialCurve2D::getCoefficientOfDenominator(const rg_INT& index) const 
{
    return coefficientsOfDenominator[index];
}

rg_REAL*    rg_RPolynomialCurve2D::getCoefficientsOfDenominator()
{
    if ( coefficientsOfDenominator == rg_NULL )
    {
        return rg_NULL;
    }

    const rg_INT degree=rg_PolynomialCurve2D::getDegree();
    rg_REAL* output=new rg_REAL[degree+1];

    for ( rg_INT i=0; i < degree+1; i++ )
    {
        output[i]=coefficientsOfDenominator[i];
    }

    return output;
}

void     rg_RPolynomialCurve2D::setCoefficientOfDenominator(const rg_INDEX& i, const rg_REAL& tCoefficient)
{
    coefficientsOfDenominator[i]=tCoefficient;
}

void     rg_RPolynomialCurve2D::setCoefficientsOfDenominator(const rg_REAL* tCoefficients)
{
    const rg_INT degree=rg_PolynomialCurve2D::getDegree();
    for ( rg_INT i=0; i < degree+1; i++ )
    {
        coefficientsOfDenominator[i]=tCoefficients[i];
    }
}

void     rg_RPolynomialCurve2D::setCurve(const rg_RPolynomialCurve2D& temp)
{
    removeAll();
    rg_PolynomialCurve2D::setCurve(temp);
    rg_INT degree =rg_PolynomialCurve2D::getDegree();

    coefficientsOfDenominator=new rg_REAL[degree+1];
    setCoefficientsOfDenominator(temp.coefficientsOfDenominator);
}

rg_RPolynomialCurve2D& rg_RPolynomialCurve2D::operator=(const rg_RPolynomialCurve2D& temp)
{
    if ( this != &temp )
    {
        setCurve(temp);
    }
    return *this;
}


