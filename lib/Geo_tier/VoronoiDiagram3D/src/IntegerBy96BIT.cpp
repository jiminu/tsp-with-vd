#include "IntegerBy96BIT.h"
using namespace V::GeometryTier;

#include <math.h>

IntegerBy96BIT::IntegerBy96BIT()
: m_sign(1)
{
	for (rg_INT i=0; i<3; i++ )
		m_significantFigure[i] = 0;
}



IntegerBy96BIT::IntegerBy96BIT(const rg_INT& intValue)
{
	m_significantFigure[0] = 0;
	m_significantFigure[1] = 0;
    if ( intValue >= 0)  {
        m_sign = 1;
    	m_significantFigure[2] = intValue;
    }
	else  {
		m_sign = -1;
		m_significantFigure[2] = -intValue;
	}
}



IntegerBy96BIT::IntegerBy96BIT(const rg_REAL& realValue, const rg_REAL& validMantissa)
{
	rg_REAL scaledReal = realValue*validMantissa;
	int     intValue   = (int) scaledReal;

	m_significantFigure[0] = 0;
	m_significantFigure[1] = 0;
	if ( intValue >= 0 )  {
		m_sign = 1;
		m_significantFigure[2] = intValue;
	}
	else  {
		m_sign = -1;
		m_significantFigure[2] = -intValue;
	}
}



IntegerBy96BIT::IntegerBy96BIT(const IntegerBy96BIT& lInteger)
{
	m_sign = lInteger.m_sign;
	for (rg_INT i=0; i<3; i++)
		m_significantFigure[i] = lInteger.m_significantFigure[i];
}



IntegerBy96BIT::~IntegerBy96BIT()
{
}




rg_INT  IntegerBy96BIT::getSign() const
{
	return m_sign;
}



unsigned int*   IntegerBy96BIT::getSignificantFigure() 
{
	return m_significantFigure;
}




void IntegerBy96BIT::setSign(const rg_INT& sign)
{
	m_sign = sign;
}



void IntegerBy96BIT::setSignificantFigure(unsigned int* significantFigure)
{
	for (rg_INT i=0; i<3; i++)
		m_significantFigure[i] = significantFigure[i];
}



void IntegerBy96BIT::setIntegerBy96Bit(const rg_INT& intValue)
{
	m_significantFigure[0] = 0;
	m_significantFigure[1] = 0;
    if ( intValue >= 0)  {
        m_sign = 1;
    	m_significantFigure[2] = intValue;
    }
	else  {
		m_sign = -1;
		m_significantFigure[2] = -intValue;
	}
}


IntegerBy96BIT IntegerBy96BIT::ABS(const IntegerBy96BIT& lInteger) const
{
    IntegerBy96BIT result;
    result.m_sign = 1;
	for (rg_INT i=0; i<3; i++)
		result.m_significantFigure[i] = lInteger.m_significantFigure[i];

    return result;
}



rg_FLAG IntegerBy96BIT::isZERO()
{
    if ( m_significantFigure[2] == 0 && m_significantFigure[1] == 0 && m_significantFigure[0] == 0)
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG IntegerBy96BIT::isNONZERO()
{
    if ( m_significantFigure[2] == 0 && m_significantFigure[1] == 0 && m_significantFigure[0] == 0)
        return rg_FALSE;
    else
        return rg_TRUE;
}



rg_FLAG IntegerBy96BIT::isGTZERO()
{
    if ( isNONZERO() && m_sign == 1 )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_REAL IntegerBy96BIT::convertToREAL() const
{
    rg_REAL realValue = m_sign*(m_significantFigure[2] + m_significantFigure[1]*pow(2.,32) + m_significantFigure[0]*pow(2.,64));

    return realValue;
}





//IntegerBy96BIT IntegerBy96BIT::multiple(const rg_REAL& realValue1, const rg_REAL& realValue2, 
//                                        const rg_REAL& realValue3, const rg_REAL& validMantissa)
//{
//}



IntegerBy96BIT& IntegerBy96BIT::operator =(const IntegerBy96BIT& lInteger)
{
	if ( this == &lInteger )
		return *this;

	m_sign = lInteger.m_sign;
	for (rg_INT i=0; i<3; i++)
		m_significantFigure[i] = lInteger.m_significantFigure[i];

	return *this;
}




IntegerBy96BIT IntegerBy96BIT::operator +(const IntegerBy96BIT& lInteger) const
{
    rg_INT       Shift16Bit    = 16;
    unsigned int maskFor16bit  = 65535;
    unsigned int maxIntOf32bit = (int)(pow(2.,32) - 1);
    unsigned int maskForOver16bit  = (int)((pow(2.,32) - 1) - (pow(2.,16) -1));

	IntegerBy96BIT result;
	if ( m_sign == lInteger.m_sign )  {
		result.m_sign = m_sign;

        unsigned int sumToJumpNext[4] = {0, 0, 0, 0};
        for (rg_INT i=2; i>=0; i-- )  {
            unsigned int dividedIntValue[2][2];
            dividedIntValue[0][0] = (m_significantFigure[i]>>Shift16Bit)&maskFor16bit;
            dividedIntValue[1][0] = (lInteger.m_significantFigure[i]>>Shift16Bit)&maskFor16bit;
            dividedIntValue[0][1] = m_significantFigure[i]&maskFor16bit;
            dividedIntValue[1][1] = lInteger.m_significantFigure[i]&maskFor16bit;

            unsigned int sumUnder16bit      = dividedIntValue[0][1] + dividedIntValue[1][1]
                                              + sumToJumpNext[i+1];
            unsigned int sumToJumpOver16bit = (sumUnder16bit>>Shift16Bit)&maskFor16bit;
            unsigned int sumOver16bit       = dividedIntValue[0][0] + dividedIntValue[1][0]
                                              + sumToJumpOver16bit;
            sumToJumpNext[i] = (sumOver16bit>>Shift16Bit)&maskFor16bit;
            result.m_significantFigure[i] = (sumUnder16bit&maskFor16bit) + ((sumOver16bit<<Shift16Bit)&maskForOver16bit);
        }
	}
	else  {
        IntegerBy96BIT absThis = ABS(*this);
        IntegerBy96BIT absThat = ABS(lInteger);

        if ( absThis == absThat ) {
        }
        else if ( absThis < absThat )  {
            result.m_sign = lInteger.m_sign;
            for (rg_INT i=2; i>=1; i--)  {
                if ( absThat.m_significantFigure[i] < absThis.m_significantFigure[i] )  {
                    absThat.m_significantFigure[i-1] = absThat.m_significantFigure[i-1]-1;

                    result.m_significantFigure[i] = maxIntOf32bit - absThis.m_significantFigure[i]
                                                    + 1           + absThat.m_significantFigure[i];
                }
                else {
                    result.m_significantFigure[i] = absThat.m_significantFigure[i] - absThis.m_significantFigure[i];
                }
            }
            result.m_significantFigure[0] = absThat.m_significantFigure[0] - absThis.m_significantFigure[0];
        }
        else  {
            result.m_sign = m_sign;
            for (rg_INT i=2; i>=1; i--)  {
                if ( absThis.m_significantFigure[i] < absThat.m_significantFigure[i] )  {
                    absThis.m_significantFigure[i-1] = absThis.m_significantFigure[i-1]-1;

                    result.m_significantFigure[i] = maxIntOf32bit - absThat.m_significantFigure[i]
                                                    + 1           + absThis.m_significantFigure[i];
                }
                else {
                    result.m_significantFigure[i] = absThis.m_significantFigure[i] - absThat.m_significantFigure[i];
                }
            }
            result.m_significantFigure[0] = absThis.m_significantFigure[0] - absThat.m_significantFigure[0];
        }
	}

    return result;
}



IntegerBy96BIT IntegerBy96BIT::operator -(const IntegerBy96BIT& lInteger) const
{
    IntegerBy96BIT copylInteger(lInteger);
    copylInteger.m_sign = copylInteger.m_sign*(-1);

    return ((*this)+copylInteger);
}



IntegerBy96BIT IntegerBy96BIT::operator *(const IntegerBy96BIT& lInteger) const
{
    rg_INT       Shift16Bit    = 16;
    unsigned int maskFor16bit  = 65535;
    //unsigned int maxIntOf32bit = (int)(pow(2.,32) - 1);

    
    IntegerBy96BIT result;

    rg_INT numNonZeroGroup = 0;
	rg_INT i=0;
	for ( i=0; i<3; i++)  {
        if ( m_significantFigure[i] != 0 )
            numNonZeroGroup++;
        if ( lInteger.m_significantFigure[i] != 0 )
            numNonZeroGroup++;
    }

    if ( numNonZeroGroup > 3 )  {
        //  can't solve this * operator;
        result.m_sign = 0;
        return result;
    }

    result.m_sign = m_sign*lInteger.m_sign;

    unsigned int dividedIntValue[2][4];
    dividedIntValue[0][3] = (m_significantFigure[1]>>Shift16Bit)&maskFor16bit;
    dividedIntValue[0][2] = m_significantFigure[1]&maskFor16bit;
    dividedIntValue[0][1] = (m_significantFigure[2]>>Shift16Bit)&maskFor16bit;
    dividedIntValue[0][0] = m_significantFigure[2]&maskFor16bit;
    dividedIntValue[1][3] = (lInteger.m_significantFigure[1]>>Shift16Bit)&maskFor16bit;
    dividedIntValue[1][2] = lInteger.m_significantFigure[1]&maskFor16bit;
    dividedIntValue[1][1] = (lInteger.m_significantFigure[2]>>Shift16Bit)&maskFor16bit;
    dividedIntValue[1][0] = lInteger.m_significantFigure[2]&maskFor16bit;

    unsigned int multipleEachDIV[4][7];
    for ( i=0; i<4; i++)  {
        for ( rg_INT j=0; j<7; j++)
            multipleEachDIV[i][j] = 0;
    }

    for ( i=0; i<4; i++)  {
        for ( rg_INT j=0; j<4; j++)  {
            multipleEachDIV[i][i+j] = dividedIntValue[0][j]*dividedIntValue[1][i];
        }
    }

    unsigned int under16Bit[7];
    unsigned int over16Bit[7];
    for ( i=0; i<7; i++)  {
        unsigned int sumOfMultiple = multipleEachDIV[0][i] + multipleEachDIV[1][i] 
                                     + multipleEachDIV[2][i] + multipleEachDIV[3][i];
        under16Bit[i] = sumOfMultiple&maskFor16bit;
        over16Bit[i]  = (sumOfMultiple>>Shift16Bit)&maskFor16bit;
    }

    unsigned int sumUnderToJumpNext   = 0;
    unsigned int sumUnderToRemainHere = 0;
    unsigned int sumOverToJumpNext    = 0;
    unsigned int sumOverToRemainHere  = 0;
    
    unsigned int sumUnder16Bit = under16Bit[0];
    sumUnderToJumpNext   = (sumUnder16Bit>>Shift16Bit)&maskFor16bit;
    sumUnderToRemainHere = sumUnder16Bit&maskFor16bit;
    
    unsigned int sumOver16Bit  = under16Bit[1] + over16Bit[0] + sumUnderToJumpNext;
    sumOverToJumpNext    = (sumOver16Bit>>Shift16Bit)&maskFor16bit;
    sumOverToRemainHere  = sumOver16Bit&maskFor16bit;

    result.m_significantFigure[2] = (sumOverToRemainHere<<Shift16Bit) + sumUnderToRemainHere;
    

    sumUnder16Bit = under16Bit[2] + over16Bit[1] + sumOverToJumpNext;
    sumUnderToJumpNext   = (sumUnder16Bit>>Shift16Bit)&maskFor16bit;
    sumUnderToRemainHere = sumUnder16Bit&maskFor16bit;
    
    sumOver16Bit  = under16Bit[3] + over16Bit[2] + sumUnderToJumpNext;
    sumOverToJumpNext    = (sumOver16Bit>>Shift16Bit)&maskFor16bit;
    sumOverToRemainHere  = sumOver16Bit&maskFor16bit;

    result.m_significantFigure[1] = (sumOverToRemainHere<<Shift16Bit) + sumUnderToRemainHere;

    sumUnder16Bit = under16Bit[4] + over16Bit[3] + sumOverToJumpNext;
    sumUnderToJumpNext   = (sumUnder16Bit>>Shift16Bit)&maskFor16bit;
    sumUnderToRemainHere = sumUnder16Bit&maskFor16bit;
    
    sumOver16Bit  = under16Bit[5] + over16Bit[4] + sumUnderToJumpNext;
    sumOverToJumpNext    = (sumOver16Bit>>Shift16Bit)&maskFor16bit;
    sumOverToRemainHere  = sumOver16Bit&maskFor16bit;
    result.m_significantFigure[0] = (sumOverToRemainHere<<Shift16Bit) + sumUnderToRemainHere;

    return result;
}



//IntegerBy96BIT IntegerBy96BIT::operator /(const IntegerBy96BIT& lInteger)
//{
//}



rg_FLAG IntegerBy96BIT::operator <(const IntegerBy96BIT& lInteger)
{
    if ( m_sign < lInteger.m_sign )
        return rg_TRUE;
    else if ( m_sign > lInteger.m_sign )
        return rg_FALSE;
    else  {
        if ( m_sign > 0 )  {
            for ( rg_INT i=0; i<3; i++ )  {
                if ( m_significantFigure[i] < lInteger.m_significantFigure[i] )
                    return rg_TRUE;
                else if (m_significantFigure[i] > lInteger.m_significantFigure[i])
                    return rg_FALSE;
                else {}
            }

            return rg_FALSE;
        }
        else  {
            for ( rg_INT i=0; i<3; i++ )  {
                if ( m_significantFigure[i] < lInteger.m_significantFigure[i] )
                    return rg_FALSE;
                else if (m_significantFigure[i] > lInteger.m_significantFigure[i])
                    return rg_TRUE;
                else {}
            }

            return rg_FALSE;
        }
    }
}




rg_FLAG IntegerBy96BIT::operator<=(const IntegerBy96BIT& lInteger)
{
    if ( m_sign < lInteger.m_sign )
        return rg_TRUE;
    else if ( m_sign > lInteger.m_sign )
        return rg_FALSE;
    else  {
        if ( m_sign > 0 )  {
            for ( rg_INT i=0; i<3; i++ )  {
                if ( m_significantFigure[i] < lInteger.m_significantFigure[i] )
                    return rg_TRUE;
                else if (m_significantFigure[i] > lInteger.m_significantFigure[i])
                    return rg_FALSE;
                else {}
            }

            return rg_TRUE;
        }
        else  {
            for ( rg_INT i=0; i<3; i++ )  {
                if ( m_significantFigure[i] < lInteger.m_significantFigure[i] )
                    return rg_FALSE;
                else if (m_significantFigure[i] > lInteger.m_significantFigure[i])
                    return rg_TRUE;
                else {}
            }

            return rg_TRUE;
        }
    }
}



rg_FLAG IntegerBy96BIT::operator==(const IntegerBy96BIT& lInteger)
{
    if ( m_sign != lInteger.m_sign )
        return rg_FALSE;
    else {
        for ( rg_INT i=0; i<3; i++ )  {
            if ( m_significantFigure[i] != lInteger.m_significantFigure[i] )
                return rg_FALSE;
        }

        return rg_TRUE;
    }
}



/*
rg_FLAG testLastPointForBeingOnBoundaryOfConvexHullBySymbolicPerturbation(
                              const rg_Point3D& point1, const rg_Point3D& point2, 
                              const rg_Point3D& point3, const rg_Point3D& lastPoint, 
                              const rg_REAL& validMantissa)
{
    
    // | 1  point1.x     point1.y    point1.z    |
    // | 1  point2.x     point2.y    point2.z    |
    // | 1  point3.x     point3.y    point3.z    |
    // | 1  lastPoint.x  lastPoint.y lastPoint.z |

    IntegerBy96BIT row1[3] = { IntegerBy96BIT(point1.getX(), validMantissa), IntegerBy96BIT(point1.getY(), validMantissa), IntegerBy96BIT(point1.getZ(), validMantissa) };
    IntegerBy96BIT row2[3] = { IntegerBy96BIT(point2.getX(), validMantissa), IntegerBy96BIT(point2.getY(), validMantissa), IntegerBy96BIT(point2.getZ(), validMantissa) };
    IntegerBy96BIT row3[3] = { IntegerBy96BIT(point3.getX(), validMantissa), IntegerBy96BIT(point3.getY(), validMantissa), IntegerBy96BIT(point3.getZ(), validMantissa) };
    IntegerBy96BIT row4[3] = { IntegerBy96BIT(lastPoint.getX(), validMantissa), IntegerBy96BIT(lastPoint.getY(), validMantissa), IntegerBy96BIT(lastPoint.getZ(), validMantissa) };

    
    IntegerBy96BIT result = computeDeterminantOfMatrixForConvexHullTest(row1, row2, row3, row4);




    if ( result.isNONZERO() )  {
        if ( result.isGTZERO() )
            return rg_TRUE;
        else
            return rg_FALSE;
    }
    else  {
        for ( rg_INT i=0; i<14; i++)  {
            result = computeDeterminantAccordingToDegreeOfPerturbingError(row1, row2, row3, row4, i);

            if ( result.isNONZERO() )  {
                if ( result.isGTZERO() )
                    return rg_TRUE;
                else
                    return rg_FALSE;
            }
        }        
    }
}
*/



IntegerBy96BIT V::GeometryTier::computeDeterminantOfMatrixForConvexHullTest(IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                                           IntegerBy96BIT* row3, IntegerBy96BIT* row4)
{
    // | 1  row1[0]  row1[1]  row1[2]  |
    // | 1  row2[0]  row2[1]  row2[2]  |
    // | 1  row3[0]  row3[1]  row3[2]  |
    // | 1  row4[0]  row4[1]  row4[2]  |
    //
    //   | row2[0]  row2[1]  row2[2] |   | row1[0]  row1[1]  row1[2] |
    // = | row3[0]  row3[1]  row3[2] | - | row3[0]  row3[1]  row3[2] |
    //   | row4[0]  row4[1]  row4[2] |   | row4[0]  row4[1]  row4[2] |
    //
    //     | row1[0]  row1[1]  row1[2] |   | row1[0]  row1[1]  row1[2] |
    //   + | row2[0]  row2[1]  row2[2] | - | row2[0]  row2[1]  row2[2] |
    //     | row4[0]  row4[1]  row4[2] |   | row3[0]  row3[1]  row3[2] |

    IntegerBy96BIT result;
    IntegerBy96BIT subResult1 = computeDeterminantOf3By3Matrix(row2, row3, row4);
    IntegerBy96BIT subResult2 = computeDeterminantOf3By3Matrix(row1, row3, row4);
    IntegerBy96BIT subResult3 = computeDeterminantOf3By3Matrix(row1, row2, row4);
    IntegerBy96BIT subResult4 = computeDeterminantOf3By3Matrix(row1, row2, row3);

    result = subResult1 - subResult2 + subResult3 - subResult4;

    return result;
}



IntegerBy96BIT V::GeometryTier::computeDeterminantOf3By3Matrix(IntegerBy96BIT* row1, IntegerBy96BIT* row2, IntegerBy96BIT* row3)
{
    // | row1[0]  row1[1]  row1[2] |
    // | row2[0]  row2[1]  row2[2] |
    // | row3[0]  row3[1]  row3[2] |
    //
    // = row1[0]| row2[1]  row2[2] | - row2[0]| row1[1]  row1[2] | + row3[0]| row1[1]  row1[2] | 
    //          | row3[1]  row3[2] |          | row3[1]  row3[2] |          | row2[1]  row2[2] |

    IntegerBy96BIT result;
    IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row2[1], row2[2], row3[1], row3[2]);
    IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row1[1], row1[2], row3[1], row3[2]);
    IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row1[1], row1[2], row2[1], row2[2]);


    subResult1 = row1[0]*subResult1;
    subResult2 = row2[0]*subResult2;
    subResult3 = row3[0]*subResult3;

    result = subResult1 - subResult2 + subResult3;

    return result;

}



IntegerBy96BIT V::GeometryTier::computeDeterminantOf2By2Matrix(const IntegerBy96BIT& ele1, const IntegerBy96BIT& ele2, 
                                              const IntegerBy96BIT& ele3, const IntegerBy96BIT& ele4)
{
    // | ele1  ele2 | = ele1*ele4 - ele2*ele3
    // | ele3  ele4 |

    IntegerBy96BIT result;
    IntegerBy96BIT subResult1 = ele1*ele4;
    IntegerBy96BIT subResult2 = ele2*ele3;

    result = subResult1 - subResult2;

    return result;
}



IntegerBy96BIT V::GeometryTier::computeDeterminantForConvexHullTestBySymbolicPerturbation(
                                     IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                     IntegerBy96BIT* row3, IntegerBy96BIT* row4)
{
    IntegerBy96BIT result;

    for ( rg_INT degree=0; degree<15; degree++)  {
        result = computeDeterminantAccordingToDegreeOfPerturbingError(row1, row2, row3, row4, degree);

        if ( result.isNONZERO() )  
            break;
    }        

    return result;
}


IntegerBy96BIT V::GeometryTier::computeDeterminantAccordingToDegreeOfPerturbingError(
                                       IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                       IntegerBy96BIT* row3, IntegerBy96BIT* row4, 
                                       const rg_INT& degreeOrderOfSymbol)
{
    //          | x_i  y_i  z_i |   | 1  row1[0]  row1[1]  row1[2]  |
    //          | x_j  y_j  z_j |   | 1  row2[0]  row2[1]  row2[2]  |
    // matrix = | x_k  y_k  z_k | = | 1  row3[0]  row3[1]  row3[2]  |
    //          | x_l  y_l  z_l |   | 1  row4[0]  row4[1]  row4[2]  |
    IntegerBy96BIT result;

    switch (degreeOrderOfSymbol)  {
        case 0:  
            result = computeDeterminantOfMatrixForConvexHullTest(row1, row2, row3, row4);
            break;

        case 1:  
            {
                // | 1  y_j  z_j |
                // | 1  y_k  z_k | 
                // | 1  y_l  z_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row3[1], row3[2], row4[1], row4[2]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row2[1], row2[2], row4[1], row4[2]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row2[1], row2[2], row3[1], row3[2]);

                result = subResult1 - subResult2 + subResult3;
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 2:  
            {
                // | 1  y_i  z_i |
                // | 1  y_k  z_k | 
                // | 1  y_l  z_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row3[1], row3[2], row4[1], row4[2]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row1[1], row1[2], row4[1], row4[2]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row1[1], row1[2], row3[1], row3[2]);

                result = subResult1 - subResult2 + subResult3;
            }
            break;

        case 3:  
            {
                // | 1  y_i  z_i |
                // | 1  y_j  z_j | 
                // | 1  y_l  z_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row2[1], row2[2], row4[1], row4[2]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row1[1], row1[2], row4[1], row4[2]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row1[1], row1[2], row2[1], row2[2]);

                result = subResult1 - subResult2 + subResult3;
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 4:  
            {
                // | 1  x_j  z_j |
                // | 1  x_k  z_k | 
                // | 1  x_l  z_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row3[0], row3[2], row4[0], row4[2]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row2[0], row2[2], row4[0], row4[2]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row2[0], row2[2], row3[0], row3[2]);

                result = subResult1 - subResult2 + subResult3;
            }
            break;

        case 5:  
            {
                // | 1  z_k | 
                // | 1  z_l |
                result = row4[2] - row3[2];
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 6:  
            {
                // | 1  z_j | 
                // | 1  z_l |
                result = row4[2] - row2[2];
            }
            break;

        case 7:  
            {
                // | 1  x_i  z_i |
                // | 1  x_k  z_k | 
                // | 1  x_l  z_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row3[0], row3[2], row4[0], row4[2]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row1[0], row1[2], row4[0], row4[2]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row1[0], row1[2], row3[0], row3[2]);

                result = subResult1 - subResult2 + subResult3;
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 8:  
            {
                // | 1  z_i | 
                // | 1  z_l |
                result = row4[2] - row1[2];
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 9:  
            {
                // | 1  x_i  z_i |
                // | 1  x_j  z_j | 
                // | 1  x_l  z_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row2[0], row2[2], row4[0], row4[2]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row1[0], row1[2], row4[0], row4[2]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row1[0], row1[2], row2[0], row2[2]);

                result = subResult1 - subResult2 + subResult3;
            }
            break;

        case 10:  
            {
                // | 1  x_j  y_j |
                // | 1  x_k  y_k | 
                // | 1  x_l  y_l | 
                IntegerBy96BIT subResult1 = computeDeterminantOf2By2Matrix(row3[0], row3[1], row4[0], row4[1]);
                IntegerBy96BIT subResult2 = computeDeterminantOf2By2Matrix(row2[0], row2[1], row4[0], row4[1]);
                IntegerBy96BIT subResult3 = computeDeterminantOf2By2Matrix(row2[0], row2[1], row3[0], row3[1]);

                result = subResult1 - subResult2 + subResult3;
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 11:  
            {
                // | 1  y_k | 
                // | 1  y_l |
                result = row4[1] - row3[1];
            }
            break;

        case 12:  
            {
                // | 1  y_j | 
                // | 1  y_l |
                result = row4[1] - row2[1];
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 13:  
            {
                // | 1  x_k | 
                // | 1  x_l |
                result = row4[0] - row3[0];
                result.m_sign = result.m_sign*-1;
            }
            break;

        case 14:  
            {
                result = IntegerBy96BIT(1);
            }
            break;

        default:
            break;
    }

    return result;
}



double V::GeometryTier::computeDeterminantReal33Matrix(const rg_Point3D& point1, const rg_Point3D& point2,
                                                     const rg_Point3D& point3)
{
    // | row1[0]  row1[1]  row1[2] |
    // | row2[0]  row2[1]  row2[2] |
    // | row3[0]  row3[1]  row3[2] |
    //
    // = row1[0]| row2[1]  row2[2] | - row2[0]| row1[1]  row1[2] | + row3[0]| row1[1]  row1[2] |
    //          | row3[1]  row3[2] |          | row3[1]  row3[2] |          | row2[1]  row2[2] |

    double sub1_2by2 = (point2.getY()*point3.getZ()) - (point2.getZ()*point3.getY());
    double sub2_2by2 = (point1.getY()*point3.getZ()) - (point1.getZ()*point3.getY());
    double sub3_2by2 = (point1.getY()*point2.getZ()) - (point1.getZ()*point2.getY());

    double subresult1 = point1.getX()*sub1_2by2;
    double subresult2 = point2.getX()*sub2_2by2;
    double subresult3 = point3.getX()*sub3_2by2;

    double result = subresult1 - subresult2 + subresult3;

    return result;
}



double V::GeometryTier::computeDeterminant( const rg_Point3D& point1, const rg_Point3D& point2, 
                           const rg_Point3D& point3, const rg_Point3D& point4, double* subResult)
{
    // | 1  row1[0]  row1[1]  row1[2]  |
    // | 1  row2[0]  row2[1]  row2[2]  |
    // | 1  row3[0]  row3[1]  row3[2]  |
    // | 1  row4[0]  row4[1]  row4[2]  |
    //
    //   | row2[0]  row2[1]  row2[2] |   | row1[0]  row1[1]  row1[2] |
    // = | row3[0]  row3[1]  row3[2] | - | row3[0]  row3[1]  row3[2] |
    //   | row4[0]  row4[1]  row4[2] |   | row4[0]  row4[1]  row4[2] |
    //
    //     | row1[0]  row1[1]  row1[2] |   | row1[0]  row1[1]  row1[2] |
    //   + | row2[0]  row2[1]  row2[2] | - | row2[0]  row2[1]  row2[2] |
    //     | row4[0]  row4[1]  row4[2] |   | row3[0]  row3[1]  row3[2] |

    subResult[0] = computeDeterminantReal33Matrix(point2, point3, point4);
    subResult[1] = computeDeterminantReal33Matrix(point1, point3, point4);
    subResult[2] = computeDeterminantReal33Matrix(point1, point2, point4);
    subResult[3] = computeDeterminantReal33Matrix(point1, point2, point3);

    double result = subResult[0] - subResult[1] + subResult[2] - subResult[3];

    return result;
}




