#ifndef _INTEGERBY96BIT_H
#define _INTEGERBY96BIT_H

#include "rg_Const.h"
#include "rg_Point3D.h"


namespace V {

namespace GeometryTier {


class IntegerBy96BIT
{
private:
	rg_INT  m_sign;
	unsigned int m_significantFigure[3];

public:
	IntegerBy96BIT();
	IntegerBy96BIT(const rg_INT& intValue);
	IntegerBy96BIT(const rg_REAL& realValue, const rg_REAL& validMantissa);
	IntegerBy96BIT(const IntegerBy96BIT& lInteger);
	~IntegerBy96BIT();

	rg_INT  getSign() const;
	unsigned int* getSignificantFigure();

	void setSign(const rg_INT& sign);
	void setSignificantFigure(unsigned int* significantFigure);
    void setIntegerBy96Bit(const rg_INT& intValue);

    IntegerBy96BIT ABS(const IntegerBy96BIT& lInteger) const;
    rg_FLAG        isZERO();
    rg_FLAG        isNONZERO();
    rg_FLAG        isGTZERO();

    rg_REAL convertToREAL() const;
    
    IntegerBy96BIT multiple(const rg_REAL& realValue1, const rg_REAL& realValue2, 
                            const rg_REAL& realValue3, const rg_REAL& validMantissa);

	IntegerBy96BIT& operator =(const IntegerBy96BIT& lInteger);
	IntegerBy96BIT  operator +(const IntegerBy96BIT& lInteger) const;
	IntegerBy96BIT  operator -(const IntegerBy96BIT& lInteger) const;
	IntegerBy96BIT  operator *(const IntegerBy96BIT& lInteger) const;
	//IntegerBy96BIT  operator /(const IntegerBy96BIT& lInteger);

    rg_FLAG          operator <(const IntegerBy96BIT& lInteger);
    rg_FLAG          operator<=(const IntegerBy96BIT& lInteger);
    rg_FLAG          operator==(const IntegerBy96BIT& lInteger);

    //friend rg_FLAG testLastPointForBeingOnBoundaryOfConvexHullBySymbolicPerturbation(
    //                          const rg_Point3D& point1, const rg_Point3D& point2, 
    //                          const rg_Point3D& point3, const rg_Point3D& lastPoint, 
    //                          const rg_REAL& validMantissa);
    friend IntegerBy96BIT computeDeterminantForConvexHullTestBySymbolicPerturbation(
                                     IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                     IntegerBy96BIT* row3, IntegerBy96BIT* row4);
        friend IntegerBy96BIT computeDeterminantAccordingToDegreeOfPerturbingError(
                                     IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                     IntegerBy96BIT* row3, IntegerBy96BIT* row4, 
                                     const rg_INT& degreeOrderOfSymbol);

        friend IntegerBy96BIT computeDeterminantOfMatrixForConvexHullTest(IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                                                          IntegerBy96BIT* row3, IntegerBy96BIT* row4);
        friend IntegerBy96BIT computeDeterminantOf3By3Matrix(IntegerBy96BIT* row1, IntegerBy96BIT* row2, IntegerBy96BIT* row3);
        friend IntegerBy96BIT computeDeterminantOf2By2Matrix(const IntegerBy96BIT& ele1, const IntegerBy96BIT& ele2, 
                                                             const IntegerBy96BIT& ele3, const IntegerBy96BIT& ele4);

        friend double computeDeterminant( const rg_Point3D& point1, const rg_Point3D& point2, 
                                          const rg_Point3D& point3, const rg_Point3D& lastPoint, double* subResult);
        friend double computeDeterminantReal33Matrix(const rg_Point3D& point1, const rg_Point3D& point2, 
                                                     const rg_Point3D& point3);
        
};


IntegerBy96BIT computeDeterminantForConvexHullTestBySymbolicPerturbation(
                         IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                         IntegerBy96BIT* row3, IntegerBy96BIT* row4);
IntegerBy96BIT computeDeterminantAccordingToDegreeOfPerturbingError(
                         IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                         IntegerBy96BIT* row3, IntegerBy96BIT* row4, 
                         const rg_INT& degreeOrderOfSymbol);

IntegerBy96BIT computeDeterminantOfMatrixForConvexHullTest(IntegerBy96BIT* row1, IntegerBy96BIT* row2, 
                                                             IntegerBy96BIT* row3, IntegerBy96BIT* row4);
IntegerBy96BIT computeDeterminantOf3By3Matrix(IntegerBy96BIT* row1, IntegerBy96BIT* row2, IntegerBy96BIT* row3);
IntegerBy96BIT computeDeterminantOf2By2Matrix(const IntegerBy96BIT& ele1, const IntegerBy96BIT& ele2, 
                                                 const IntegerBy96BIT& ele3, const IntegerBy96BIT& ele4);

double computeDeterminant( const rg_Point3D& point1, const rg_Point3D& point2, 
                             const rg_Point3D& point3, const rg_Point3D& lastPoint, double* subResult);
double computeDeterminantReal33Matrix(const rg_Point3D& point1, const rg_Point3D& point2, 
                                         const rg_Point3D& point3);
        

} // namespace GeometryTier

} // namespace V


#endif

