#include "EdgeFlipTimeCalculator.h"
using namespace V::GeometryTier;

#include "JenkinsTraub.h"
#include <fstream>
using namespace std;

rg_Polynomial EdgeFlipTimeCalculator::make_polynomial_equation_for_vedge_flipping(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder)
{
     rg_Polynomial result(0);

    if (all_four_disks_have_same_radius(fourDiskGeneratorsInNonIncreasingOrder))
    {
        rg_Polynomial XYMat[3][3];

        fill_elements_of_all_polynomial_XY_matrix(FOR_FINDING_TIME, fourDiskGeneratorsInNonIncreasingOrder, XYMat);

        rg_Polynomial XYDet = make_determinant_polynomial(XYMat);

        //make final equation
        result = XYDet;
    }
    else
    {
        rg_Polynomial XMat[3][3], YMat[3][3], XYMat[3][3];

        fill_elements_of_all_polynomial_matrices(FOR_FINDING_TIME, fourDiskGeneratorsInNonIncreasingOrder, XMat, YMat, XYMat);

        rg_Polynomial XDet  = make_determinant_polynomial(XMat);
        rg_Polynomial YDet  = make_determinant_polynomial(YMat);
        rg_Polynomial XYDet = make_determinant_polynomial(XYMat);

        //make final equation
        result = (XDet ^ 2) + (YDet ^ 2) - (XYDet ^ 2);
    }

    return result;
}


std::pair<bool, double> EdgeFlipTimeCalculator::solve_polynomial_equation_for_vedge_flipping(
    const PolynomialSolverType& solverType, 
    const double& mostRecentlyEventOccurringTime, 
    const rg_Polynomial& polynomialEq,
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    const bool& isThisVEdgeFlipping)
{
#ifdef WRITING_POLYNOMIAL_RELATED_INFORMATION
    ___store_current_time(mostRecentlyEventOccurringTime);
#endif

    pair<bool, double> nextEdgeFlipTimeRelativeToCurrTime (false, 0.0);

    switch (solverType)
    {
        case LABSRC:
        {
	        nextEdgeFlipTimeRelativeToCurrTime = solve_polynomial_equation_for_vedge_flipping_by_LABSRC(polynomialEq,
                                                                                          fourDiskGeneratorsInNonIncreasingOrder,
                                                                                          isThisVEdgeFlipping); 
            break;
        }

        case JENKINS_TRAUB:
        {
            nextEdgeFlipTimeRelativeToCurrTime = solve_polynomial_equation_for_vedge_flipping_by_JENKINS_TRAUB(polynomialEq,
                                                                                                 fourDiskGeneratorsInNonIncreasingOrder,
                                                                                                 isThisVEdgeFlipping);
            break;
        }

        case STURMS_SEQ :
        {
            nextEdgeFlipTimeRelativeToCurrTime = solve_polynomial_equation_for_vedge_flipping_by_STURMS_SEQUENCE(polynomialEq,
                                                                                                   fourDiskGeneratorsInNonIncreasingOrder,
                                                                                                   isThisVEdgeFlipping);         
            break;
        }
    }

    if (nextEdgeFlipTimeRelativeToCurrTime.first)
    {
        nextEdgeFlipTimeRelativeToCurrTime.second += mostRecentlyEventOccurringTime;
    }


    return nextEdgeFlipTimeRelativeToCurrTime;
}


bool EdgeFlipTimeCalculator::all_four_disks_have_same_radius(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder)
{
    //check all circle's radius are same
    double largestRadius = fourDiskGeneratorsInNonIncreasingOrder.front()->getCircle().getRadius();

    int sameRadiusNum = 0;

    for (int i = 1; i <= 3; i++)
    {
        double currRadius = fourDiskGeneratorsInNonIncreasingOrder.at(i)->getCircle().getRadius();

        if (rg_EQ(currRadius, largestRadius, TOLERANCE_OF_EQUAL_RADIUS))
        {
            sameRadiusNum++;
        }
    }

    if (sameRadiusNum == 3)
    {
        return true;
    }
    else {
        return false;
    }
}


void EdgeFlipTimeCalculator::fill_elements_of_all_polynomial_XY_matrix(
    EdgeFlippingFormulaGeneratingType mode, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    rg_Polynomial(*XYMat)[3])
{
    //element equation
    rg_Polynomial x1_x4(1);  //x1_x4 this implies x1 - x4
    rg_Polynomial x2_x4(1);
    rg_Polynomial x3_x4(1);
    rg_Polynomial y1_y4(1);
    rg_Polynomial y2_y4(1);
    rg_Polynomial y3_y4(1);
    rg_Polynomial r1_r4(0);
    rg_Polynomial r2_r4(0);
    rg_Polynomial r3_r4(0);
    rg_Polynomial p1(2);
    rg_Polynomial p2(2);
    rg_Polynomial p3(2);


    fill_elements_for_polynomial_matrix(fourDiskGeneratorsInNonIncreasingOrder,
                                        x1_x4, x2_x4, x3_x4,
                                        y1_y4, y2_y4, y3_y4,
                                        r1_r4, r2_r4, r3_r4,
                                        p1,    p2,    p3);
	//make equation 3
	switch (mode)
	{
	case FOR_FINDING_TIME:
		XYMat[0][0] = x1_x4;
		XYMat[0][1] = y1_y4;
		XYMat[0][2] = p1;
		XYMat[1][0] = x2_x4;
		XYMat[1][1] = y2_y4;
		XYMat[1][2] = p2;
		XYMat[2][0] = x3_x4;
		XYMat[2][1] = y3_y4;
		XYMat[2][2] = p3;
		break;

	case FOR_VERIFING_TIME:
		XYMat[0][0] = x1_x4;
		XYMat[0][1] = y1_y4;
		XYMat[0][2] = r1_r4;
		XYMat[1][0] = x2_x4;
		XYMat[1][1] = y2_y4;
		XYMat[1][2] = r2_r4;
		XYMat[2][0] = x3_x4;
		XYMat[2][1] = y3_y4;
		XYMat[2][2] = r3_r4;
		break;

	default:
		break;
	}
}


void EdgeFlipTimeCalculator::fill_elements_of_all_polynomial_matrices(
    EdgeFlippingFormulaGeneratingType mode, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    rg_Polynomial(*XMat)[3], 
    rg_Polynomial(*YMat)[3], 
    rg_Polynomial(*XYMat)[3])
{
    //element equation
    rg_Polynomial x1_x4(1);  //x1_x4 this implies x1 - x4
    rg_Polynomial x2_x4(1);
    rg_Polynomial x3_x4(1);
    rg_Polynomial y1_y4(1);
    rg_Polynomial y2_y4(1);
    rg_Polynomial y3_y4(1);
    rg_Polynomial r1_r4(0);
    rg_Polynomial r2_r4(0);
    rg_Polynomial r3_r4(0);
    rg_Polynomial p1(2);
    rg_Polynomial p2(2);
    rg_Polynomial p3(2);


    fill_elements_for_polynomial_matrix(fourDiskGeneratorsInNonIncreasingOrder,
                                        x1_x4, x2_x4, x3_x4,
                                        y1_y4, y2_y4, y3_y4,
                                        r1_r4, r2_r4, r3_r4,
                                        p1,    p2,    p3);

	//make equation 1
	XMat[0][0] = x1_x4;
	XMat[0][1] = r1_r4;
	XMat[0][2] = p1;
	XMat[1][0] = x2_x4;
	XMat[1][1] = r2_r4;
	XMat[1][2] = p2;
	XMat[2][0] = x3_x4;
	XMat[2][1] = r3_r4;
	XMat[2][2] = p3;

	//make equation 2
	YMat[0][0] = y1_y4;
	YMat[0][1] = r1_r4;
	YMat[0][2] = p1;
	YMat[1][0] = y2_y4;
	YMat[1][1] = r2_r4;
	YMat[1][2] = p2;
	YMat[2][0] = y3_y4;
	YMat[2][1] = r3_r4;
	YMat[2][2] = p3;

	//make equation 3
	switch (mode)
	{
	case FOR_FINDING_TIME:
		XYMat[0][0] = x1_x4;
		XYMat[0][1] = y1_y4;
		XYMat[0][2] = p1;
		XYMat[1][0] = x2_x4;
		XYMat[1][1] = y2_y4;
		XYMat[1][2] = p2;
		XYMat[2][0] = x3_x4;
		XYMat[2][1] = y3_y4;
		XYMat[2][2] = p3;
		break;

	case FOR_VERIFING_TIME:
		XYMat[0][0] = x1_x4;
		XYMat[0][1] = y1_y4;
		XYMat[0][2] = r1_r4;
		XYMat[1][0] = x2_x4;
		XYMat[1][1] = y2_y4;
		XYMat[1][2] = r2_r4;
		XYMat[2][0] = x3_x4;
		XYMat[2][1] = y3_y4;
		XYMat[2][2] = r3_r4;
		break;

	default:
		break;
	}
}



void EdgeFlipTimeCalculator::fill_elements_for_polynomial_matrix(const vector<DynamicDisk*>& movingDiskSet,
                                                                 rg_Polynomial& x1_x4, rg_Polynomial& x2_x4, rg_Polynomial& x3_x4, 
                                                                 rg_Polynomial& y1_y4, rg_Polynomial& y2_y4, rg_Polynomial& y3_y4, 
                                                                 rg_Polynomial& r1_r4, rg_Polynomial& r2_r4, rg_Polynomial& r3_r4, 
                                                                 rg_Polynomial& p1,    rg_Polynomial& p2,    rg_Polynomial& p3)
{
    x1_x4.setCoefficient(0, movingDiskSet[0]->getCircle().getCenterPt().getX() - movingDiskSet[3]->getCircle().getCenterPt().getX());
    x1_x4.setCoefficient(1, movingDiskSet[0]->getVelocityVectorX() - movingDiskSet[3]->getVelocityVectorX());

    x2_x4.setCoefficient(0, movingDiskSet[1]->getCircle().getCenterPt().getX() - movingDiskSet[3]->getCircle().getCenterPt().getX());
    x2_x4.setCoefficient(1, movingDiskSet[1]->getVelocityVectorX() - movingDiskSet[3]->getVelocityVectorX());

    x3_x4.setCoefficient(0, movingDiskSet[2]->getCircle().getCenterPt().getX() - movingDiskSet[3]->getCircle().getCenterPt().getX());
    x3_x4.setCoefficient(1, movingDiskSet[2]->getVelocityVectorX() - movingDiskSet[3]->getVelocityVectorX());


    y1_y4.setCoefficient(0, movingDiskSet[0]->getCircle().getCenterPt().getY() - movingDiskSet[3]->getCircle().getCenterPt().getY());
    y1_y4.setCoefficient(1, movingDiskSet[0]->getVelocityVectorY() - movingDiskSet[3]->getVelocityVectorY());

    y2_y4.setCoefficient(0, movingDiskSet[1]->getCircle().getCenterPt().getY() - movingDiskSet[3]->getCircle().getCenterPt().getY());
    y2_y4.setCoefficient(1, movingDiskSet[1]->getVelocityVectorY() - movingDiskSet[3]->getVelocityVectorY());

    y3_y4.setCoefficient(0, movingDiskSet[2]->getCircle().getCenterPt().getY() - movingDiskSet[3]->getCircle().getCenterPt().getY());
    y3_y4.setCoefficient(1, movingDiskSet[2]->getVelocityVectorY() - movingDiskSet[3]->getVelocityVectorY());

    if (movingDiskSet[0]->this_disk_is_container() && movingDiskSet[1]->getCircle().isIncludedIn(movingDiskSet[0]->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        r1_r4.setCoefficient(0, -movingDiskSet[0]->getCircle().getRadius() - movingDiskSet[3]->getCircle().getRadius());
    }
    else
    {
        r1_r4.setCoefficient(0, movingDiskSet[0]->getCircle().getRadius() - movingDiskSet[3]->getCircle().getRadius());
    }

    r2_r4.setCoefficient(0, movingDiskSet[1]->getCircle().getRadius() - movingDiskSet[3]->getCircle().getRadius());
    r3_r4.setCoefficient(0, movingDiskSet[2]->getCircle().getRadius() - movingDiskSet[3]->getCircle().getRadius());



    p1 = (x1_x4 ^ 2) + (y1_y4 ^ 2) - (r1_r4 ^ 2);
	p2 = (x2_x4 ^ 2) + (y2_y4 ^ 2) - (r2_r4 ^ 2);
	p3 = (x3_x4 ^ 2) + (y3_y4 ^ 2) - (r3_r4 ^ 2);
}



rg_Polynomial EdgeFlipTimeCalculator::make_determinant_polynomial(rg_Polynomial(*matrix)[3])
{
	return (matrix[0][0])*(matrix[1][1])*(matrix[2][2])
		 + (matrix[0][1])*(matrix[1][2])*(matrix[2][0])
		 + (matrix[0][2])*(matrix[1][0])*(matrix[2][1])
		 - (matrix[0][2])*(matrix[1][1])*(matrix[2][0])
		 - (matrix[0][1])*(matrix[1][0])*(matrix[2][2])
		 - (matrix[0][0])*(matrix[1][2])*(matrix[2][1]);
}



pair<bool, double> EdgeFlipTimeCalculator::solve_polynomial_equation_for_vedge_flipping_by_LABSRC(
    const rg_Polynomial & polynomialEq,  
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    const bool& isThisVEdgeFlipping)
{
    pair<bool, double> output(false, DBL_MAX);
  
    rg_ComplexNumber* roots = NULL;

    const int polynomialDegree = polynomialEq.getDegree();

    rg_Polynomial temp_polynomialEq(polynomialEq);
    roots = temp_polynomialEq.solve();

    //3. find valid real roots from all complex roots 
    list<double> positiveRealRoots;
    find_positive_real_roots(roots, polynomialDegree, positiveRealRoots);

    if (isThisVEdgeFlipping)
    {
        double minAbsoluteRealNumber = find_min_absolute_real_root(roots, polynomialDegree);
        erase_numerically_zero_root_from_positive_real_roots(positiveRealRoots, minAbsoluteRealNumber);
    }

    //4. check each real roots by a gavrilova condition
    double solution = DBL_MAX;
    for (list<double>::const_iterator it_RealNumber = positiveRealRoots.begin(); it_RealNumber != positiveRealRoots.end(); it_RealNumber++)
    {
        double realRoot = *it_RealNumber;

        if (four_disks_are_circumscribing(realRoot, fourDiskGeneratorsInNonIncreasingOrder))
        {
            if (realRoot < solution)
            {
                solution = realRoot;
            }
        }
    }

    if (solution != DBL_MAX)
    {
        output.first = true;
        output.second = solution;
    }


#ifdef WRITING_POLYNOMIAL_RELATED_INFORMATION
    ___check_roots_of_polynomial_for_statistics(polynomialDegree, roots, fourDiskGeneratorsInNonIncreasingOrder);

    char type = ___change_bool_type_of_flipping_status_to_char(isThisVEdgeFlipping);
    ___store_this_polynomial_related_info(fourDiskGeneratorsInNonIncreasingOrder, polynomialEq, roots, type);
#endif

    delete[] roots;

    return output;
}



pair<bool, double> EdgeFlipTimeCalculator::solve_polynomial_equation_for_vedge_flipping_by_JENKINS_TRAUB(
    const rg_Polynomial& polynomialEq, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    const bool& isThisVEdgeFlipping)
{
    pair<bool, double> output(false, DBL_MAX);

    //1. find roots of polynomial
    const int polynomialDegree = polynomialEq.getDegree();

    double* zeror = new double[polynomialDegree];
    double* zeroi = new double[polynomialDegree];
    int* info     = new int[polynomialDegree + 1];

    double* coeffs = new double[polynomialDegree + 1];
    for (int i = 0; i <= polynomialDegree; i++)
    {
        coeffs[i] = polynomialEq.getCoefficient(polynomialDegree - i);
    }

    int numRoot = rpoly(coeffs, polynomialDegree, zeror, zeroi, info);

    //2. convert root's format as rg_ComplexNumber
    rg_ComplexNumber* roots = new rg_ComplexNumber[polynomialDegree];

    for (int i = 0; i < polynomialDegree; i++)
    {
        roots[i].setRealNumber(zeror[i]);
        roots[i].setImaginaryNumber(zeroi[i]);
    }

    //3. find valid real roots from all complex roots 
    list<double> positiveRealRoots;
    find_positive_real_roots(roots, polynomialDegree, positiveRealRoots);

    if (isThisVEdgeFlipping)
    {
        double minAbsoluteRealNumber = find_min_absolute_real_root(roots, polynomialDegree);
        erase_numerically_zero_root_from_positive_real_roots(positiveRealRoots, minAbsoluteRealNumber);
    }

    //4. check each real roots by a gavrilova condition
    double solution = DBL_MAX;
    for (list<double>::const_iterator it_RealNumber = positiveRealRoots.begin(); it_RealNumber != positiveRealRoots.end(); it_RealNumber++)
    {
        double realRoot = *it_RealNumber;

        if (four_disks_are_circumscribing(realRoot, fourDiskGeneratorsInNonIncreasingOrder))
        {
            if (realRoot < solution)
            {
                solution = realRoot;
            }
        }
    }

    if (solution != DBL_MAX)
    {
        output.first = true;
        output.second = solution;
    }

#ifdef WRITING_POLYNOMIAL_RELATED_INFORMATION
    ___check_roots_of_polynomial_for_statistics(polynomialDegree, roots, fourDiskGeneratorsInNonIncreasingOrder);

    char type = ___change_bool_type_of_flipping_status_to_char(isThisVEdgeFlipping);
    ___store_this_polynomial_related_info(fourDiskGeneratorsInNonIncreasingOrder, polynomialEq, roots, type);
#endif

    delete[] zeror;
    delete[] zeroi;
    delete[] info;
    delete[] coeffs;
    delete[] roots;



    return output;
}


std::pair<bool, double> EdgeFlipTimeCalculator::solve_polynomial_equation_for_vedge_flipping_by_STURMS_SEQUENCE(
    const rg_Polynomial& polynomialEq, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    const bool& isThisVEdgeFlipping, 
    const double& upperBoundOfTimeInterval /*= DBL_MAX*/)
{
    pair<bool, double> output(false, DBL_MAX);

    //. find flipping time by using sturm's theory
    bool   solutionIsComputed = false;
    double solution           = DBL_MAX;

    list<rg_Polynomial> sturmsChain = make_sturms_chain(polynomialEq);

    pair<double, double> currInterval;

    if (isThisVEdgeFlipping)
    {
        double lowerBound = find_lower_bound_for_flipping_VEdge(polynomialEq, sturmsChain, -0.0382, 0.0618);
        currInterval      = make_pair(lowerBound, upperBoundOfTimeInterval);
    }
    else
    {
        currInterval = make_pair(0.0, upperBoundOfTimeInterval);
    }


    list<pair<double, double>> searchIntervalStack;
    searchIntervalStack.push_back(currInterval);

    while (!searchIntervalStack.empty())
    {
        currInterval = searchIntervalStack.back();
        searchIntervalStack.pop_back();

        int numOfRoot = find_number_of_roots_in_interval_using_sturms_chain(sturmsChain, currInterval.first, currInterval.second);

        switch (numOfRoot)
        {
        case 0:
        {
            break;
        }

        case 1:
        {
            double root = do_false_position_method(polynomialEq, 
                                                   currInterval.first, currInterval.second, 
                                                   MAX_ITERATION_OF_FALSE_POSITION, TOLERANCE_OF_FALSE_POSITION);

            if (four_disks_are_circumscribing(root, fourDiskGeneratorsInNonIncreasingOrder))
            {
                solution = root;
                solutionIsComputed = true;
            }

            break;
        }

        default:
        {
            double halfPointOfInterval = (currInterval.first + currInterval.second)/2.0;

            if (halfPointOfInterval < currInterval.first)
            {
                halfPointOfInterval = currInterval.first;
            }
            else if (halfPointOfInterval > currInterval.second)
            {
                halfPointOfInterval = currInterval.second;
            }

            pair<double, double> lowerHalfInterval = make_pair(currInterval.first, halfPointOfInterval);
            pair<double, double> upperHalfInterval = make_pair(halfPointOfInterval, currInterval.second);

            searchIntervalStack.push_back(upperHalfInterval);
            searchIntervalStack.push_back(lowerHalfInterval);

            break;
        }
        }


        if (solutionIsComputed)
        {
            break;
        }
    }
 
    if (solution != DBL_MAX)
    {
        output.first  = true;
        output.second = solution;
    }

    return output;
}


void EdgeFlipTimeCalculator::find_positive_real_roots(rg_ComplexNumber* roots,
                                                       const int & numOfRoots,
                                                       list<double>& positiveRealRoots)
{
    list<double> tempPositiveRealRoots;

    for (int i = 0; i < numOfRoots; i++)
    {
        if (roots[i].isPureRealNumber())
        {
            double realPart = roots[i].getRealNumber();
            tempPositiveRealRoots.push_back(realPart);
        }
    }

    for (list<double>::const_iterator it_RealRoot = tempPositiveRealRoots.begin(); it_RealRoot != tempPositiveRealRoots.end(); it_RealRoot++)
    {
        double realRoot = *it_RealRoot;

        if (realRoot > 0.0)
        {
            positiveRealRoots.push_back(realRoot);
        }
    }
}


double EdgeFlipTimeCalculator::find_min_absolute_real_root(rg_ComplexNumber * roots, const int & numOfRoots)
{
    double minAbsoluteRealNumber = DBL_MAX;

    for (int i = 0; i < numOfRoots; i++)
    {
        if (roots[i].isPureRealNumber())
        {
            double realPart = roots[i].getRealNumber();
            
            if (abs(realPart) < minAbsoluteRealNumber)
            {
                minAbsoluteRealNumber = abs(realPart);
            }
        }
    }

    return minAbsoluteRealNumber;
}



void EdgeFlipTimeCalculator::erase_numerically_zero_root_from_positive_real_roots(
    list<double>& positiveRealRoots, 
    const double& minAbsoluteRealNumber) 
{
    if (minAbsoluteRealNumber < TOLERANCE_OF_EDGE_FLIP_ROOT)
    {
        for (list<double>::const_iterator it_RealRoot = positiveRealRoots.begin(); it_RealRoot != positiveRealRoots.end(); it_RealRoot++)
        {
            double realRoot = *it_RealRoot;

            if (realRoot == minAbsoluteRealNumber)
            {
                positiveRealRoots.erase(it_RealRoot);
                break;
            }
        }
    }

}


bool EdgeFlipTimeCalculator::four_disks_are_circumscribing(
    const double& flippingTime, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder)
{
	if (fourDiskGeneratorsInNonIncreasingOrder.size() != 4)
	{
		return false;
	}

	
	if (all_four_disks_have_same_radius(fourDiskGeneratorsInNonIncreasingOrder))
    {
		return true;
	}


    //1999. Gavrilova. Swap Conditions ... method 

	rg_Polynomial XMat[3][3], YMat[3][3], XYMat[3][3];

	fill_elements_of_all_polynomial_matrices(FOR_VERIFING_TIME, fourDiskGeneratorsInNonIncreasingOrder, XMat, YMat, XYMat);

	rg_Polynomial XDet = make_determinant_polynomial(XMat);
	rg_Polynomial YDet = make_determinant_polynomial(YMat);
	rg_Polynomial XYDet = make_determinant_polynomial(XYMat);

	rg_Polynomial xi_x4[3];
    rg_Polynomial yi_y4[3];
    rg_Polynomial pi[3];

	xi_x4[0] = XMat[0][0];
	xi_x4[1] = XMat[1][0];
	xi_x4[2] = XMat[2][0];

    yi_y4[0] = YMat[0][0];
	yi_y4[1] = YMat[1][0];
	yi_y4[2] = YMat[2][0];

	pi[0]    = XMat[0][2];
	pi[1]    = XMat[1][2];
	pi[2]    = XMat[2][2];

    
    bool isInequalitySatisfied       = false;
    int  conditionOfOutterSideDisks  = 3;

    //DynamicDisk* largestMovingDisk       = fourDiskGeneratorsInNonIncreasingOrder.front();
    DynamicDisk* largestMovingDisk       = fourDiskGeneratorsInNonIncreasingOrder[0];
    DynamicDisk* secondLargestMovingDisk = fourDiskGeneratorsInNonIncreasingOrder[1];

    if (largestMovingDisk->this_disk_is_container() && secondLargestMovingDisk->getCircle().isIncludedIn(largestMovingDisk->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        conditionOfOutterSideDisks = 2;
    }

    int countOfOutterSideDisk = 0;

    for (int i = 0; i < 3; i++)
    {
        double pi_Val = pi[i].evaluatePolynomial(flippingTime);

        double result = ((xi_x4[i])*YDet - (yi_y4[i])*XDet).evaluatePolynomial(flippingTime) / XYDet.evaluatePolynomial(flippingTime);

        if ((result - TOLERANCE_OF_GAVRILOVA_INEQUALITY < pi_Val))
        {
            countOfOutterSideDisk++;
        }
    }

    if (countOfOutterSideDisk == conditionOfOutterSideDisks)
    {
        isInequalitySatisfied = true;
    }
    
    return isInequalitySatisfied;
}


double EdgeFlipTimeCalculator::do_false_position_method(const rg_Polynomial & polynomial, 
                                                          const double& lowerBound, 
                                                          const double& upperBound, 
                                                          const int& maxIteration, 
                                                          const double& tolerance)
{
    double a     = lowerBound;
    double b     = upperBound;
    double f_a   = polynomial.evaluatePolynomial(a);
    double f_b   = polynomial.evaluatePolynomial(b);

    double c     = 0.0;
    double f_c   = 0.0;
    double c_old = 0.0;

    double currRelativeErr = 0.0;

    int a_counter = 0;
    int b_counter = 0;

    if (f_a * f_b > 0)
    {
        return DBL_MAX;
    }

    int prevOption = -1;
    int currOption = -1;

    for (int i = 0; i < maxIteration; i++)
    {
        c_old = c;
        c     = (a*f_b - b*f_a) / (f_b - f_a);

        if ((c > b) || (c < a))
        {
            c = (a + b) / 2;
        }

        f_c = polynomial.evaluatePolynomial(c);


        if (c != 0)
        {
            currRelativeErr = abs((c - c_old) / c) * 100;
        }


        double sign = f_a * f_c;

        if (sign > 0) // b is fixed
        {
            a   = c;
            f_a = f_c;

            a_counter  = 0;
            b_counter += 1;

            if (b_counter >= 2)
            {
                f_b = f_b / 2;
            }
        }
        else if(sign < 0) // a is fixed
        {
            b   = c;
            f_b = f_c;

            b_counter  = 0;
            a_counter += 1;

            if (a_counter >= 2)
            {
                f_a = f_a / 2;
            }
        }
        else
        {
            currRelativeErr = 0.0;
        }

        if ((currRelativeErr < tolerance) || (abs(f_c) < tolerance))
        {
            break;
        }
    }

    return c;
}


int EdgeFlipTimeCalculator::find_number_of_roots_in_interval_using_sturms_chain(
    const list<rg_Polynomial>& sturmsChain, 
    const double& intervalBegin, 
    const double& intervalEnd)
{
    if (intervalBegin >= intervalEnd)
    {
        return 0;
    }

    int sturmTestAtIntervalBegin = count_sign_transition_in_sturms_chain(sturmsChain, intervalBegin);
    int sturmTestAtIntervalEnd   = count_sign_transition_in_sturms_chain(sturmsChain, intervalEnd);

    int numberOfRootsInInterval = sturmTestAtIntervalBegin - sturmTestAtIntervalEnd;
    return numberOfRootsInInterval;
}



list<rg_Polynomial> EdgeFlipTimeCalculator::make_sturms_chain(const rg_Polynomial& polynomial)
{
    list<rg_Polynomial> sturmsChain;
    sturmsChain.push_back(polynomial);
    sturmsChain.push_back(polynomial.makeDerivative());

    list<rg_Polynomial>::const_iterator itForDividend = sturmsChain.begin();
    list<rg_Polynomial>::const_iterator itForDivider = ++sturmsChain.begin();

    while ((*itForDivider).getDegree() > 0)
    {
        const rg_Polynomial& dividend = *itForDividend;
        const rg_Polynomial& divider  = *itForDivider;

        pair<rg_Polynomial, rg_Polynomial> polynomialDivisionResult = divide_polynomial(dividend, divider);
        rg_Polynomial& quotient  = polynomialDivisionResult.first;
        rg_Polynomial& remainder = polynomialDivisionResult.second;

        sturmsChain.push_back(rg_Polynomial(remainder.getDegree()) - remainder);
        itForDividend++;
        itForDivider++;
    }

    //cysong
    for (list<rg_Polynomial>::iterator it_Polynomial = sturmsChain.begin(); it_Polynomial != sturmsChain.end(); it_Polynomial++)
    {
        rg_Polynomial& currPolynomial = *it_Polynomial;

        int    currPolynomialDegree = currPolynomial.getDegree();
        double coeffOfMaxOrder      = currPolynomial.getCoefficient(currPolynomialDegree);

        if (rg_NE(currPolynomial.getCoefficient(coeffOfMaxOrder), 1.0))
        {
            currPolynomial.setCoefficient(currPolynomialDegree, 1.0);

            for (int i = 0; i < currPolynomialDegree; i++)
            {
                double newCoefficient = currPolynomial.getCoefficient(i) / coeffOfMaxOrder;
                currPolynomial.setCoefficient(i, newCoefficient);
            }
        }
    }

    return sturmsChain;
}


pair<rg_Polynomial, rg_Polynomial> EdgeFlipTimeCalculator::divide_polynomial(
    const rg_Polynomial& dividend, 
    const rg_Polynomial& divider)
{
    //returning polynomial pair means <quotient, remainder>
    rg_Polynomial quotient(dividend.getDegree() - divider.getDegree());
    double dividerLeadingCoefficient = divider.getCoefficient(divider.getDegree());

    rg_Polynomial remainder = dividend;
    while (remainder.getDegree() >= divider.getDegree())
    {
        double remainderLeadingCoefficient = remainder.getCoefficient(remainder.getDegree());
        double leadingCoeffsRatio = remainderLeadingCoefficient / dividerLeadingCoefficient;
        rg_Polynomial tempDividend = remainder;
        int quotientTermDegree = tempDividend.getDegree() - divider.getDegree();
        quotient.setCoefficient(quotientTermDegree, leadingCoeffsRatio);
        //cout << "Quotient at degree[" << quotientTermDegree << "]: " << leadingCoeffsRatio << endl;
        //cout << "Remainder degree: " << remainder.getDegree() - 1 << endl;
        remainder.setDegree(remainder.getDegree() - 1);

        for (int i = quotientTermDegree; i < tempDividend.getDegree(); i++)
        {
            //cout << "Remainder term [" << i << "]: " << tempDividend.getCoefficient(i) << " - " << quotient.getCoefficient(quotientTermDegree) << " * " << divider.getCoefficient(i - quotientTermDegree) <<" = "<< tempDividend.getCoefficient(i) - quotient.getCoefficient(quotientTermDegree)*divider.getCoefficient(i - quotientTermDegree)<<endl;
            remainder.setCoefficient(i, tempDividend.getCoefficient(i) - quotient.getCoefficient(quotientTermDegree)*divider.getCoefficient(i - quotientTermDegree));
        }
        for (int i = 0; i < quotientTermDegree; i++)
        {
            //cout << "Remainder term [" << i << "]: " << tempDividend.getCoefficient(i) << endl;
            remainder.setCoefficient(i, tempDividend.getCoefficient(i));
        }
        remainder.reconstruct();
    }

    return make_pair(quotient, remainder);
}


int EdgeFlipTimeCalculator::count_sign_transition_in_sturms_chain(
    const list<rg_Polynomial>& sturmsChain, 
    const double& valueForEvaluation)
{
    //cout << "Count sign transition for " << valueForEvaluation << endl;
    int signTransitionCounter = 0;
    double lastEvaluation = 0;
    int chainCounter = 0;
    for (list<rg_Polynomial>::const_iterator itForChain = sturmsChain.begin(); itForChain != sturmsChain.end(); itForChain++)
    {
        double currEvaluation = (*itForChain).evaluatePolynomial(valueForEvaluation);
        //cout << "p[" << chainCounter << "]: " << currEvaluation << ", ";
        chainCounter++;

        if (lastEvaluation * currEvaluation < 0)
        {
            signTransitionCounter++;
        }
        lastEvaluation = currEvaluation;
    }
    //cout << endl;

    return signTransitionCounter;
}


double EdgeFlipTimeCalculator::find_lower_bound_for_flipping_VEdge(
    const rg_Polynomial & polynomial, 
    const list<rg_Polynomial>& sturmsChain, 
    const double & lowerInterval, 
    const double & upperInterval)
{
    bool   isLowerBoundFound = false;
    double lowerBound = upperInterval;

    pair<double, double> currInterval = make_pair(lowerInterval, upperInterval);

    list<pair<double, double>> searchIntervalStack;
    searchIntervalStack.push_back(currInterval);

    while (!searchIntervalStack.empty())
    {
        currInterval = searchIntervalStack.back();
        searchIntervalStack.pop_back();

        int numOfRoot = find_number_of_roots_in_interval_using_sturms_chain(sturmsChain, currInterval.first, currInterval.second);

        switch (numOfRoot)
        {
            case 0:
            {
                break;
            }

            case 1:
            {
                if ((currInterval.first <= 0.0) && (currInterval.second >= 0.0))
                {
                    lowerBound = currInterval.second;
                    isLowerBoundFound = true;
                }

                break;
            }

            default:
            {
                double halfPointOfInterval = (currInterval.first + currInterval.second) / 2.0;

                if (halfPointOfInterval < currInterval.first)
                {
                    halfPointOfInterval = currInterval.first;
                }
                else if (halfPointOfInterval > currInterval.second)
                {
                    halfPointOfInterval = currInterval.second;
                }

                pair<double, double> lowerHalfInterval = make_pair(currInterval.first, halfPointOfInterval);
                pair<double, double> upperHalfInterval = make_pair(halfPointOfInterval, currInterval.second);

                searchIntervalStack.push_back(upperHalfInterval);
                searchIntervalStack.push_back(lowerHalfInterval);

                break;
            }
        }


        if (isLowerBoundFound)
        {
            break;
        }
    }

    return lowerBound;
}



#ifdef WRITING_POLYNOMIAL_RELATED_INFORMATION


int    EdgeFlipTimeCalculator::numOfWriting     = 0;
string EdgeFlipTimeCalculator::outputFolderPath = ".\\";

const int    EdgeFlipTimeCalculator::LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY         = 100000;
const string EdgeFlipTimeCalculator::polynomialRootsFileNameWithExtension       = "polynomial_root_counts.txt";
const string EdgeFlipTimeCalculator::polynomialCoefficientFileNameWithExtension = "polynomial_coefficients.txt";
const string EdgeFlipTimeCalculator::movingDiskFileNameWithExtension            = "dynamic_disks.txt";
const string EdgeFlipTimeCalculator::realRootFileNameWithExtension              = "real_roots.txt";
const string EdgeFlipTimeCalculator::realRootEvaluationFileNameWithExtension    = "real_roots_evalution.txt";
const string EdgeFlipTimeCalculator::VEdgeFlippingStatusFileNameWithExtension   = "flipping_status.txt";


void EdgeFlipTimeCalculator::___count_polynomial_root_type(const RootTypeOfPolynomial& rootType)
{
    if (m_CountStatisticsOfRoots.size() == 0)
    {
        m_CountStatisticsOfRoots.reset(RootTypeOfPolynomial::POLY_ROOT_COUNT_SIZE);
    }

    m_CountStatisticsOfRoots.addCount(rootType);
}


void EdgeFlipTimeCalculator::___check_roots_of_polynomial_for_statistics(
    const int& numOfRoots, 
    rg_ComplexNumber* roots, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder)
{
    int numOfRealRoots                             = 0;
    int numOfPositiveRealRoots                     = 0;
    int numOfPositiveRealRootsSatisfyingInequality = 0;

    for (int i = 0; i < numOfRoots; i++)
    {
        if (roots[i].isPureRealNumber())
        {
            numOfRealRoots++;

            if ((roots[i].getRealNumber() > TOLERANCE_OF_EDGE_FLIP_ROOT))
            {
                numOfPositiveRealRoots++;

                if (EdgeFlipTimeCalculator::four_disks_are_circumscribing(roots[i].getRealNumber(), fourDiskGeneratorsInNonIncreasingOrder))
                {
                    numOfPositiveRealRootsSatisfyingInequality++;
                }
            }
        }
    }


    switch (numOfRealRoots)
    {
    case 0:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_0);
        break;
    case 1:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_1);
        break;
    case 2:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_2);
        break;
    case 3:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_3);
        break;
    case 4:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_4);
        break;
    case 5:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_5);
        break;
    case 6:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_6);
        break;
    case 7:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_7);
        break;
    case 8:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_8);
        break;
    default:
        break;
    }


    switch (numOfPositiveRealRoots)
    {
    case 0:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_0);
        break;
    case 1:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_1);
        break;
    case 2:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_2);
        break;
    case 3:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_3);
        break;
    case 4:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_4);
        break;
    case 5:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_5);
        break;
    case 6:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_6);
        break;
    case 7:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_7);
        break;
    case 8:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_8);
        break;
    default:
        break;
    }


    switch (numOfPositiveRealRootsSatisfyingInequality)
    {
    case 0:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_0);
        break;
    case 1:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_1);
        break;
    case 2:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_2);
        break;
    case 3:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_3);
        break;
    case 4:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_4);
        break;
    case 5:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_5);
        break;
    case 6:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_6);
        break;
    case 7:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_7);
        break;
    case 8:
        ___count_polynomial_root_type(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_8);
        break;
    default:
        break;
    }

}


void EdgeFlipTimeCalculator::___store_current_time(const double& currentTime)
{
    m_SimulationTimesAtMakingPolynomial.push_back(currentTime);
}


void EdgeFlipTimeCalculator::___store_this_polynomial_related_info(
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder,
    const rg_Polynomial& polynomial,
    const rg_ComplexNumber* roots,
    const char& isThisFlippingVEdge)
{
    if (m_Polynomials.size() >= LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY)
    {
        ___write_polynomial_related_info_in_temporal_storage();
    }

    m_Polynomials.push_back(polynomial);
    m_FlippingStatusForVEdge.push_back(isThisFlippingVEdge);

    vector<DynamicDisk> quadrupletOfMovingDisks;

    for (int i = 0; i < fourDiskGeneratorsInNonIncreasingOrder.size(); i++)
    {
        DynamicDisk* currMovingDisk = fourDiskGeneratorsInNonIncreasingOrder.at(i);

        quadrupletOfMovingDisks.push_back(*currMovingDisk);
    }

    m_QuadrupletsOfMovingDisks.push_back(quadrupletOfMovingDisks);


    vector<double> realRoots;

    for (int i = 0; i < polynomial.getDegree(); i++)
    {
        if (roots[i].isPureRealNumber())
        {
            realRoots.push_back(roots[i].getRealNumber());
        }
    }

    m_RealRoots.push_back(realRoots);
}


void EdgeFlipTimeCalculator::___clear_all_polynomial_related_info()
{
    //m_CountStatisticsOfRoots.reset(0);
    m_SimulationTimesAtMakingPolynomial.clear();
    m_QuadrupletsOfMovingDisks.clear();
    m_Polynomials.clear();
    m_RealRoots.clear();
    m_FlippingStatusForVEdge.clear();
}


void EdgeFlipTimeCalculator::___write_polynomial_related_info_in_temporal_storage()
{
    if (numOfWriting == 0)
    {
        ___write_headers_for_each();
    }

    ___write_dynamic_disks_info();
    ___write_polynomial_coefficient_info();
    ___write_polynomial_roots_info();
    ___write_real_roots_info();
    ___write_roots_evaluation_info();
    ___write_VEdge_flipping_status();
    ___clear_all_polynomial_related_info();

    ++numOfWriting;
}


void EdgeFlipTimeCalculator::___write_headers_for_each()
{
    string fileName1(outputFolderPath + movingDiskFileNameWithExtension);
    ofstream fout1(fileName1.c_str(), ios_base::out);

    fout1 << "Number"  << "\t" << "Time"    << "\t";
    fout1 << "Disk1-x" << "\t" << "Disk1-y" << "\t" << "Disk1-radius" << "\t" << "Velocity Vector1-x" << "\t" << "Velocity Vector1-y" << "\t";
    fout1 << "Disk2-x" << "\t" << "Disk2-y" << "\t" << "Disk2-radius" << "\t" << "Velocity Vector2-x" << "\t" << "Velocity Vector2-y" << "\t";
    fout1 << "Disk3-x" << "\t" << "Disk3-y" << "\t" << "Disk3-radius" << "\t" << "Velocity Vector3-x" << "\t" << "Velocity Vector3-y" << "\t";
    fout1 << "Disk4-x" << "\t" << "Disk4-y" << "\t" << "Disk4-radius" << "\t" << "Velocity Vector4-x" << "\t" << "Velocity Vector4-y" << endl;

    fout1.close();


    string fileName2(outputFolderPath + polynomialCoefficientFileNameWithExtension);
    ofstream fout2(fileName2.c_str(), ios_base::out);

    fout2 << "Number" << "\t" << "Time" << "\t";
    fout2 << "Degree" << "\t" << "Coef-0" << "\t" << "Coef-1" << "\t" << "Coef-2" << "\t" << "Coef-3" << "\t";
    fout2 << "Coef-4" << "\t" << "Coef-5" << "\t" << "Coef-6" << "\t" << "Coef-7" << "\t" << "Coef-8" << endl;

    fout2.close();


    string fileName3(outputFolderPath + realRootFileNameWithExtension);
    ofstream fout3(fileName3.c_str(), ios_base::out);

    fout3 << "Number" << "\t" << "Time"              << "\t";
    fout3 << "Degree" << "\t" << "Num of Real Roots" << "\t" << "Root-1" << "\t" << "Root-2" << "\t" << "Root-3" << "\t" ;
    fout3 << "Root-4" << "\t" << "Root-5"            << "\t" << "Root-6" << "\t" << "Root-7" << "\t" << "Root-8" << endl;

    fout3.close();


    string fileName4(outputFolderPath + realRootEvaluationFileNameWithExtension);
    ofstream fout4(fileName4.c_str(), ios_base::out);

    fout4 << "##########################################################################################################"   << endl;
    fout4 << "# *. #POSITIVE    : NUMBER OF POSITIVE REAL ROOTS (TOLERANCE_OF_ROOT_IN_DYNAMIC_VD)"                          << endl;
    fout4 << "# *. #INEQULITY   : NUMBER OF POSITIVE REAL ROOTS SATISFYING INEQAULITY (TOLERANCE_OF_GAVRILOVA_INEQUALITY_IN_DYNAMIC_VD)"  << endl;
    fout4 << "# *. #NEAR 0      : NUMBER OF ROOTS NEAR O (TOLERANCE_OF_ROOT_IN_DYNAMIC_VD)"                                 << endl;
    fout4 << "# *. #DOUBLE ROOTS: NUMBER OF DOUBLE ROOTS (TOLERANCE_OF_ROOT_IN_DYNAMIC_VD)"                                 << endl;
    fout4 << "# *. ROOT         : ROOT OF POLYNOMIAL (Root1, 2, ..., and 8)"                                                << endl;
    fout4 << "# *. EVAL         : POLYNOMIAL EVAULATION BY INSERTING (Root1, 2, ..., and 8)"                                << endl;
    fout4 << "# *. INEQULITY    : IF SATISFY THE INEQUALITY, CHECKED (O)                      "                             << endl;
    fout4 << "##########################################################################################################"   << endl << endl;


    fout4 << "Number"     << "\t" << "Time"       << "\t" << "Degree"       << "\t" << "Number of Real Roots" << "\t";
    fout4 << "#Positive"  << "\t" << "#Inequality"<< "\t" << "#Near 0"      << "\t" << "#Double Roots" << "\t";
    fout4 << "Root1"      << "\t" << "Root2"      << "\t" << "Root3"        << "\t" << "Root4"         << "\t";
    fout4 << "Root5"      << "\t" << "Root6"      << "\t" << "Root7"        << "\t" << "Root8"         << "\t";
    fout4 << "Eval1"      << "\t" << "Eval2"      << "\t" << "Eval3"        << "\t" << "Eval4"         << "\t";
    fout4 << "Eval5"      << "\t" << "Eval6"      << "\t" << "Eval7"        << "\t" << "Eval8"         << "\t";
    fout4 <<"Ineqaulity1" << "\t" << "Ineqaulity2"<< "\t" << "Ineqaulity3"  << "\t" << "Ineqaulity4"   << "\t";
    fout4 <<"Ineqaulity5" << "\t" << "Ineqaulity6"<< "\t" << "Ineqaulity7"  << "\t" << "Ineqaulity8"   << endl;

    fout4.close();
}


void EdgeFlipTimeCalculator::___write_polynomial_roots_info()
{
    string fileName(outputFolderPath + movingDiskFileNameWithExtension);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);

    fout << "[ 1. REAL ROOTS COUNT ]" << endl << endl;

    fout << "1. Number 0 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_0) << endl;
    fout << "2. Number 1 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_1) << endl;
    fout << "3. Number 2 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_2) << endl;
    fout << "4. Number 3 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_3) << endl;
    fout << "5. Number 4 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_4) << endl;
    fout << "6. Number 5 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_5) << endl;
    fout << "7. Number 6 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_6) << endl;
    fout << "8. Number 7 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_7) << endl;
    fout << "9. Number 8 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_REAL_ROOT_8) << endl;

    fout << endl << endl;


    fout << "[ 2. POSITIVE REAL ROOTS COUNT]" << endl << endl;

    fout << "1. Number 0 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_0) << endl;
    fout << "2. Number 1 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_1) << endl;
    fout << "3. Number 2 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_2) << endl;
    fout << "4. Number 3 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_3) << endl;
    fout << "5. Number 4 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_4) << endl;
    fout << "6. Number 5 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_5) << endl;
    fout << "7. Number 6 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_6) << endl;
    fout << "8. Number 7 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_7) << endl;
    fout << "9. Number 8 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_8) << endl;

    fout << endl << endl;


    fout << "[ 3. POSITIVE REAL ROOTS SATISFYING INEQUALITY COUNT]" << endl << endl;

    fout << "1. Number 0 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_0) << endl;
    fout << "2. Number 1 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_1) << endl;
    fout << "3. Number 2 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_2) << endl;
    fout << "4. Number 3 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_3) << endl;
    fout << "5. Number 4 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_4) << endl;
    fout << "6. Number 5 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_5) << endl;
    fout << "7. Number 6 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_6) << endl;
    fout << "8. Number 7 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_7) << endl;
    fout << "9. Number 8 :" << "\t" << m_CountStatisticsOfRoots.count(POLY_ROOT_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_8) << endl;

    fout << endl << endl;

    fout.close();
}


void EdgeFlipTimeCalculator::___write_dynamic_disks_info()
{
    string fileName(outputFolderPath + movingDiskFileNameWithExtension);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);

    //ofstream fout(movingDiskFileNameWithPath.c_str(), ios_base::out | ios_base::app);

    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(6);


    for (int i = 0; i < m_QuadrupletsOfMovingDisks.size(); i++)
    {
        const vector<DynamicDisk>& currQuadruplet = m_QuadrupletsOfMovingDisks.at(i);

        fout << LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY * numOfWriting + (i + 1) << "\t" << m_SimulationTimesAtMakingPolynomial.at(i) << "\t";

        for (int j = 0; j < currQuadruplet.size(); j++)
        {
            const DynamicDisk currDisk = currQuadruplet.at(j);

            fout << currDisk.getCircle().getCenterPt().getX() << "\t" << currDisk.getCircle().getCenterPt().getY() << "\t" << currDisk.getCircle().getRadius() << "\t";
            fout << currDisk.getVelocityVectorX() << "\t" << currDisk.getVelocityVectorY() << "\t";
        }

        fout << endl;
    }

    fout.close();
}


void EdgeFlipTimeCalculator::___write_polynomial_coefficient_info()
{
    string fileName(outputFolderPath + polynomialCoefficientFileNameWithExtension);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);

    //ofstream fout(polynomialCoefficientFileNameWithPath.c_str(), ios_base::out | ios_base::app);

    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(12);

    for (int i = 0; i < m_Polynomials.size(); i++)
    {
        const rg_Polynomial& currPolynomial = m_Polynomials.at(i);

        const int degree = currPolynomial.getDegree();

        fout << LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY * numOfWriting + (i + 1) << "\t" << m_SimulationTimesAtMakingPolynomial.at(i) << "\t" << degree << "\t";

        for (int j = 0; j <= degree; j++)
        {
            fout << currPolynomial.getCoefficient(j) << "\t";
        }

        fout << endl;
    }

    fout.close();
}



void EdgeFlipTimeCalculator::___write_real_roots_info()
{
    string fileName(outputFolderPath + realRootFileNameWithExtension);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);

    //ofstream fout(realRootFileNameWithPath.c_str(), ios_base::out | ios_base::app);

    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(6);

    for (int i = 0; i < m_RealRoots.size(); i++)
    {
        const vector<double>& currRealRoots = m_RealRoots.at(i);
        const rg_Polynomial& currPolynomial = m_Polynomials.at(i);

        const int degree = currPolynomial.getDegree();
        const int numOfRealRoots = currRealRoots.size();

        fout << LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY * numOfWriting + (i + 1) << "\t" << m_SimulationTimesAtMakingPolynomial.at(i) << "\t" << degree << "\t" << numOfRealRoots << "\t";

        for (int j = 0; j < numOfRealRoots; j++)
        {
            fout << currRealRoots.at(j) << "\t";
        }

        fout << endl;
    }

    fout.close();
}


void EdgeFlipTimeCalculator::___write_roots_evaluation_info()
{
    string fileName(outputFolderPath + realRootEvaluationFileNameWithExtension);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);

    //ofstream fout(realRootEvaluationFileNameWithPath.c_str(), ios_base::out | ios_base::app);

    for (int i = 0; i < m_RealRoots.size(); i++)
    {
        const vector<double>& currRealRoots = m_RealRoots.at(i);
        const rg_Polynomial& currPolynomial = m_Polynomials.at(i);
        const vector<DynamicDisk>& currQuadruplet = m_QuadrupletsOfMovingDisks.at(i);

        vector<DynamicDisk*> currQuadrupletOfMovingDiskPointers;

        for (int j = 0; j < currQuadruplet.size(); j++)
        {
            DynamicDisk* currMovingDiskPointer = const_cast<DynamicDisk*>(&currQuadruplet.at(j));
            currQuadrupletOfMovingDiskPointers.push_back(currMovingDiskPointer);
        }

        vector<bool> isSatisfyingInequality;

        const double currTime = m_SimulationTimesAtMakingPolynomial.at(i);
        const int    degree = currPolynomial.getDegree();
        const int    numOfRealRoots = currRealRoots.size();

        int    numOfPostiveRealRoots = 0;
        int    numOfPostivieRealRootsSatisfyingInequality = 0;
        int    numOfRealRootsNearZero = 0;
        int    numOfDoubleRealRoots = 0;

        for (int j = 0; j < numOfRealRoots; j++)
        {
            const double currRealRoot = currRealRoots.at(j);

            if (currRealRoot > TOLERANCE_OF_EDGE_FLIP_ROOT)
            {
                numOfPostiveRealRoots++;

                if (four_disks_are_circumscribing(currRealRoot, currQuadrupletOfMovingDiskPointers))
                {
                    numOfPostivieRealRootsSatisfyingInequality++;

                    isSatisfyingInequality.push_back(true);
                }
                else
                {
                    isSatisfyingInequality.push_back(false);
                }
            }
            else
            {
                isSatisfyingInequality.push_back(false);
            }

            //near zero check
            if (rg_EQ(currRealRoot, 0.0, TOLERANCE_OF_EDGE_FLIP_ROOT))
            {
                numOfRealRootsNearZero++;
            }
        }


        //double roots check

        for (int j = 0; j < numOfRealRoots; j++)
        {
            for (int k = j + 1; k < numOfRealRoots; k++)
            {
                if (rg_EQ(currRealRoots.at(k), currRealRoots.at(j), TOLERANCE_OF_EDGE_FLIP_ROOT))
                {
                    numOfDoubleRealRoots++;
                }
            }
        }

        fout << LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY * numOfWriting + (i + 1) << "\t" << currTime << "\t" << degree << "\t";
        fout << numOfRealRoots << "\t" << numOfPostiveRealRoots << "\t" << numOfPostivieRealRootsSatisfyingInequality << "\t";
        fout << numOfRealRootsNearZero << "\t" << numOfDoubleRealRoots << "\t";

        fout.setf(std::ios::fixed, std::ios::floatfield);
        fout.precision(6);

        for (int j = 0; j < 8; j++)
        {
            if (j < numOfRealRoots)
            {
                fout << currRealRoots.at(j) << "\t";
            }
            else
            {
                fout << "\t";
            }
        }

        for (int j = 0; j < 8; j++)
        {
            if (j < numOfRealRoots)
            {
                fout << currPolynomial.evaluatePolynomial(currRealRoots.at(j)) << "\t";
            }
            else
            {
                fout << "\t";
            }
        }

        for (int j = 0; j < 8; j++)
        {
            if (j < numOfRealRoots)
            {
                switch (isSatisfyingInequality.at(j))
                {
                case true:
                    fout << "O" << "\t";
                    break;

                case false:
                    fout << "X" << "\t";
                    break;
                }
            }
            else
            {
                fout << "\t";
            }
        }

        fout << endl;
    }

    fout.close();
}



void EdgeFlipTimeCalculator::___write_VEdge_flipping_status()
{
    string fileName(outputFolderPath + VEdgeFlippingStatusFileNameWithExtension);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);

    //ofstream fout(realRootFileNameWithPath.c_str(), ios_base::out | ios_base::app);

    for (int i = 0; i < m_FlippingStatusForVEdge.size(); i++)
    {
        fout << LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY * numOfWriting + (i + 1) << "\t" << m_FlippingStatusForVEdge.at(i) << endl;
    }

    fout.close();
}


char EdgeFlipTimeCalculator::___change_bool_type_of_flipping_status_to_char(const bool& type)
{
    if (type == true)
    {
        return 'T';
    }
    else
    {
        return 'F';
    }
}


#endif //WRITING_POLYNOMIAL_RELATED_INFORMATION
