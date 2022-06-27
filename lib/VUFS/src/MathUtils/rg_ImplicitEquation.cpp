/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_ImplicitEquation.cpp
//	  
//    DESCRIPTION : 
//           This is the implementaion of the class rg_ImplicitEquation
//           which defines the inplicit form of curve
//                          
//	  CLASS NAME  : rg_ImplicitEquation
//
//    BASE CLASS  : None
//      
//
//    AUTHOR      : Tae-Bum Jang, Soon-Woong Lee
//
//    START DATE  : 8.18. 1997    
//
//    HISTORY     :
//
//          By Lee, Soon-Woong, 18 Aug. 1997
//				make constructor : rg_ImplicitEquation( const rg_DEGREE& tDegree,
//                                                   const rg_REAL* coeff);
//          By Lee, Soon-Woong, 5 Sep. 1997
//              make function : rg_REAL substitute(const rg_REAL& x
//                                              const rg_REAL& y) const;
//          By Ryu, Jung-Hyun,  Feb. 17. 1998
//              make operator : rg_REAL* operator[](const rg_INDEX& index);
//
//          By Ryu, Jung-Hyun,  Feb. 19. 1998
//              make operator : rg_ImplicitEquation operator*(const rg_ImplicitEquation& operand);
//              update method: rg_REAL rg_ImplicitEquation::getCoeff( const rg_INT& powerX, const rg_INT& powerY ) const
//              make operator : rg_ImplicitEquation operator+=(const rg_ImplicitEquation& temp) const
//
//          By Ryu, Jung-Hyun,  Feb. 21. 1998
//              update method : rg_ImplicitEquation  operator+( const rg_ImplicitEquation&   temp ) const;
//                              rg_ImplicitEquation  operator-( const rg_ImplicitEquation&   temp ) const;
//	                            rg_ImplicitEquation  operator*( const rg_ImplicitEquation& operand) const;
//
//          By Ryu, Jung-Hyun,  Mar. 6. 1998
//	            rg_ImplicitEquation  operator^( const rg_DEGREE index) const;
//              
/////////////////////////////////////////////////////////////////////

#include <math.h>
#include "rg_ImplicitEquation.h"

//// private functions 
void rg_ImplicitEquation::removeAll()
{
    if ( coefficient == rg_NULL )
    {
        degree=-1;
        return;
    }

    for( rg_INDEX i=0; i < degree+1; i++ )
    {
        delete []( coefficient[i] );
    }
    delete []coefficient;

    degree=-1;
    coefficient=rg_NULL;
}

//// constructors and destructor         //////////////////////                  
rg_ImplicitEquation::rg_ImplicitEquation()
{
    coefficient=rg_NULL;
    degree=-1;
}

rg_ImplicitEquation::rg_ImplicitEquation( const rg_DEGREE& tDegree)
{
    degree=tDegree;

    coefficient = new rg_REAL*[degree+1];
    for( rg_INDEX i=0 ; i < degree+1; i++ )
    {
        coefficient[i] = new rg_REAL[degree+1-i];
        for( rg_INDEX j=0; j < degree+1-i; j++ )
        {
            coefficient[i][j]=0.0;
        }
    }
}

rg_ImplicitEquation::rg_ImplicitEquation( const rg_DEGREE& tDegree, const rg_REAL* coeff)
{
    degree=tDegree;

    coefficient = new rg_REAL*[degree+1];
    rg_INT i = 0;
	for( i=0 ; i < degree+1; i++ )
    {
        coefficient[i] = new rg_REAL[degree+1-i];
    }

    rg_INDEX k=0;
    for( i=0; i < degree+1; i++)
    {
        for( rg_INDEX j=0; j < i+1; j++)
        {
            coefficient[i-j][j] = coeff[k++];
        }
    }
}

rg_ImplicitEquation::rg_ImplicitEquation( const rg_ImplicitEquation& temp )
{
    degree=temp.degree;

    coefficient = new rg_REAL*[degree+1];
    for( rg_INDEX i=0 ; i < degree+1; i++ )
    {
        coefficient[i] = new rg_REAL[degree+1-i];
        for( rg_INDEX j=0; j < degree+1-i; j++ )
        {
            coefficient[i][j]=temp.coefficient[i][j];
        }
    }
}

rg_ImplicitEquation::~rg_ImplicitEquation()
{
    removeAll();
}

//// get functions                       //////////////////////
rg_INT rg_ImplicitEquation::getSizeOfCoeff() const
{
    return (degree+1)*(degree+2)/2; 
}

rg_REAL* rg_ImplicitEquation::getAllCoeff() const
{
    rg_INT size=getSizeOfCoeff();

    rg_REAL* output=rg_NULL;
    output = new rg_REAL[size];
    rg_INT i = 0;
	for(i=0 ; i < size; i++ )
    {
        output[i]=0.0;
    }

    rg_INDEX k=0;
    for( i=0; i < degree+1; i++)
    {
        for( rg_INDEX j=0; j < i+1; j++)
        {
            output[k]=coefficient[i-j][j];
            k++;
        }
    }

    return output;
}
            
rg_REAL rg_ImplicitEquation::getCoeff( const rg_INT& powerX,
                                 const rg_INT& powerY ) const
{
	if((powerX + powerY > degree) || 
	   (powerX < 0) || (powerY < 0))
	{
		return 0.0;
	}
	else
	{
		return coefficient[powerX][powerY];
	}
}

rg_DEGREE rg_ImplicitEquation::getDegree() const
{
    return degree;
}

rg_REAL rg_ImplicitEquation::evaluateImpEquation(const rg_REAL& valueForX, const rg_REAL& valueForY) const
{
	rg_REAL result = 0.0;
    for( rg_INDEX i=0 ; i < degree+1; i++ )
    {
        for( rg_INDEX j=0; j < degree+1-i; j++ )
        {
            result += coefficient[i][j] * pow(valueForX, i) * pow(valueForY, j);
        }
    }

	return result;
}

//// set functions                       //////////////////////
void rg_ImplicitEquation::setAllCoeff( const rg_REAL*  tCoeff,
                                    const rg_INT      size )
{
    if ( size != getSizeOfCoeff() )
    {
        return;
    }

    rg_INDEX k=0;
    for( rg_INDEX i=0; i < degree+1; i++)
    {
        for( rg_INDEX j=0; j < i+1; j++)
        {
            coefficient[i-j][j]=tCoeff[k];
            k++;
        }
    }
}

void rg_ImplicitEquation::setCoeff( const rg_INT&  powerX,
                                 const rg_INT&  powerY,
                                 const rg_REAL& tCoeff )
{
    if (   powerX+powerY > degree
        || powerX        < 0
        || powerY        < 0      )
    {
        return;
    }

    coefficient[powerX][powerY]=tCoeff;
}
void rg_ImplicitEquation::setDegree( const rg_DEGREE& tDegree)
{
    removeAll();
    degree=tDegree;

    coefficient = new rg_REAL*[degree+1];
    for( rg_INT i=0 ; i < degree+1; i++ )
    {
        coefficient[i] = new rg_REAL[degree+1-i];
        for( rg_INT j=0; j < degree+1-i; j++ )
        {
            coefficient[i][j]=0.0;
        }
    }
}

// reduce the storage by deleting the all-zero-coefficient terms in some degree
// (note that this case may happen during the operation among the implicit equations!!)
void rg_ImplicitEquation::reconstruct()
{
	rg_DEGREE newDegree = degree;

	rg_FLAG isConsecutive = rg_TRUE;
	rg_INDEX k = 0;
	for(rg_INDEX j = degree;j >= 0;j--)
	{
		rg_FLAG isAllZero = rg_TRUE;
		for(rg_INDEX i = 0;i <= j;i++)
		{
			if(coefficient[ i ][degree - k - i] != 0)
				isAllZero = rg_FALSE;
		}
		k++;
		if(isAllZero)
			newDegree--;
		else
			isConsecutive = rg_FALSE;
		
		if(!isConsecutive)
			break;
	}

	rg_ImplicitEquation temp(newDegree);

    for( rg_INDEX i=0 ; i < newDegree+1; i++ )
    {
        for( rg_INDEX j=0; j < newDegree+1-i; j++ )
        {
            temp[i][j] = coefficient[i][j];
        }
    }
	removeAll();
	(*this) = temp;
}

//// function for substitutions          //////////////////////
rg_Polynomial rg_ImplicitEquation::substitute( const rg_Polynomial& polynomialX,
                                         const rg_Polynomial& polynomialY ) const
{
    rg_Polynomial output(degree*polynomialX.getDegree());

    for( rg_INDEX i=0 ; i < degree+1 ; i++ )
    {
        for( rg_INDEX j=0 ; j < degree+1-i; j++)
        {
            output =  output
                    + coefficient[i][j]*polynomialX.power(i)*polynomialY.power(j);
        }
    }
    return output;
}

rg_REAL rg_ImplicitEquation::substitute( const rg_REAL& x,
                                   const rg_REAL& y) const
{
    rg_REAL output = 0;

    for( rg_INDEX i=0 ; i < degree+1 ; i++ )
    {
        for( rg_INDEX j=0 ; j < degree+1-i; j++)
        {
            output =  output
                    + coefficient[i][j]*pow(x, i)*pow(y, j);
        }
    }
    return output;
}

//// overloaded operator
rg_ImplicitEquation  rg_ImplicitEquation::operator+( const rg_ImplicitEquation&   temp ) const
{
	rg_DEGREE maxDegree;
	if(degree > temp.degree)
	{
		maxDegree = degree;
	}
	else
	{
		maxDegree = temp.degree;
	}

    rg_ImplicitEquation output(maxDegree);

	for(rg_INT i = 0;i < maxDegree + 1;i++)
	{
		for(rg_INT j = 0;j < maxDegree + 1 - i;j++)
		{
			if(j >= degree + 1 - i)
			{
				output.coefficient[ i ][ j ] = temp.coefficient[ i ][ j ];
			}
			else if(j >= temp.degree + 1 - i)
			{
				output.coefficient[ i ][ j ] = coefficient[ i ][ j ];
			}
			else
			{
				output.coefficient[ i ][ j ] = coefficient[ i ][ j ] + temp.coefficient[ i ][ j ];
			}
		}
	}
/*
    if ( degree == temp.degree )
    {
        output.setDegree(degree);
        for( rg_INDEX i=0 ; i < degree+1; i++ )
        {
            for( rg_INDEX j=0; j < degree+1-i; j++)
            {
                output.coefficient[i][j] +=   coefficient[i][j]
                                            + temp.coefficient[i][j];
            }
        }
    }
*/

	output.reconstruct();

    return output;
}


rg_ImplicitEquation&  rg_ImplicitEquation::operator+=( const rg_ImplicitEquation&   temp )
{
    //removeAll();
    //setDegree(temp.degree);

	/*
    for( rg_INDEX i=0 ; i < degree+1; i++ )
    {
        for( rg_INDEX j=0; j < degree+1-i; j++ )
        {
            coefficient[i][j]=coefficient[i][j]+temp.coefficient[i][j];
        }
    }*/
	
	*this = *this + temp;


    return *this;

	//return ((*this) + temp);
}


rg_ImplicitEquation  rg_ImplicitEquation::operator-( const rg_ImplicitEquation&   temp ) const
{
	rg_DEGREE maxDegree;
	if(degree > temp.degree)
	{
		maxDegree = degree;
	}
	else
	{
		maxDegree = temp.degree;
	}

    rg_ImplicitEquation output(maxDegree);

	for(rg_INT i = 0;i < maxDegree + 1;i++)
	{
		for(rg_INT j = 0;j < maxDegree + 1 - i;j++)
		{
			if(j >= degree + 1 - i)
			{
				output.coefficient[ i ][ j ] = - temp.coefficient[ i ][ j ];
			}
			else if(j >= temp.degree + 1 - i)
			{
				output.coefficient[ i ][ j ] = coefficient[ i ][ j ];
			}
			else
			{
				output.coefficient[ i ][ j ] = coefficient[ i ][ j ] - temp.coefficient[ i ][ j ];
			}
		}
	}

/*
    if ( degree == temp.degree )
    {
        output.setDegree(degree);
        for( rg_INDEX i=0 ; i < degree+1; i++ )
        {
            for( rg_INDEX j=0; j < degree+1-i; j++)
            {
                output.coefficient[i][j] =   coefficient[i][j]
                                           - temp.coefficient[i][j];
            }
        }
    }
*/
	output.reconstruct();

    return output;
}

rg_ImplicitEquation  rg_ImplicitEquation::operator*( const rg_REAL& scalar ) const
{
    rg_ImplicitEquation output(degree);

    for( rg_INDEX i=0 ; i < degree+1; i++ )
    {
        for( rg_INDEX j=0; j < degree+1-i; j++ )
        {
            output.coefficient[i][j] =
                          scalar*coefficient[i][j];
        }
    }

	output.reconstruct();

    return output;
}


// The following is the relationship 
// between degree and storage
// storage = 1(constant) + 2 dimensional * (degree) + Summation{i = degree}(i - 1)
// for example degree = 4
// storage = 1 + 2 * 4 + {(1-1) + (2-1) + (3-1) + (4-1)} = 15
// Of course, this figure is equal to
// the total number of cell for upper triangular of 5 x 5 matrix.

rg_ImplicitEquation  rg_ImplicitEquation::operator*( const rg_ImplicitEquation& operand) const
{
	rg_DEGREE newDegree = degree + operand.degree;
	rg_ImplicitEquation output(newDegree);
	//rg_INT k = degree, l = operand.degree;

	for(rg_INT i = 0;i < degree + 1;i++)
	{
		for(rg_INT j = 0;j < degree + 1 - i;j++)
		{
			for(rg_INT k = 0;k < operand.degree + 1;k++)
			{
				for(rg_INT l = 0;l < operand.degree + 1 - k;l++)
				{
					output[i + k][j + l] += getCoeff(i, j)*operand.getCoeff(k, l);
				}
			}
		}
	}

	output.reconstruct();

	return output;
}

rg_ImplicitEquation rg_ImplicitEquation::operator^( const rg_DEGREE index) const
{
	rg_ImplicitEquation temp(0);
	temp[ 0 ][ 0 ] = 1.0;

	for(rg_INT i = 0;i < index;i++)
		temp = temp * (*this);
		
	return temp;
}

rg_ImplicitEquation& rg_ImplicitEquation::operator=( const rg_ImplicitEquation&   temp ) 
{
    removeAll();
    setDegree(temp.degree);

    for( rg_INDEX i=0 ; i < degree+1; i++ )
    {
        for( rg_INDEX j=0; j < degree+1-i; j++ )
        {
            coefficient[i][j]=temp.coefficient[i][j];
        }
    }

    return *this;
}

rg_REAL* rg_ImplicitEquation::operator[](const rg_INDEX& index)
{
	return coefficient[ index ];
}

//// auxilary fuctions to test          //////////////////////
void rg_ImplicitEquation::show() const
{
	/*
    if ( coefficient == rg_NULL )
    {
        return;
    }

    for( rg_INDEX i=0; i < degree+1; i++)
    {
        for( rg_INDEX j=0; j < degree+1-i; j++)
        {
            cout.width(5);
            cout<<coefficient[i][j]<<"x^"<<i<<"y^"<<j<< "  ";
        }
        cout<<endl;
    }
	*/
}
//// friend functions and classes       //////////////////////
rg_ImplicitEquation  operator*( const rg_REAL&             scalar,
                             const rg_ImplicitEquation&   temp )
{
    rg_ImplicitEquation output(temp.degree);

    for( rg_INDEX i=0 ; i < temp.degree+1; i++ )
    {
        for( rg_INDEX j=0; j < temp.degree+1-i; j++ )
        {
            output.coefficient[i][j] =
                     scalar*temp.coefficient[i][j];
        }
    }

	output.reconstruct();

    return output;
}



