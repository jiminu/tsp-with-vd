#ifndef _EDGE_FLIP_TIME_CALCULATOR_
#define _EDGE_FLIP_TIME_CALCULATOR_

//#define WRITING_POLYNOMIAL_RELATED_INFORMATION

#include "SimpleTypesForDVD2.h"
#include "DefaultValuesForDVD2.h"
#include "rg_Polynomial.h"
#include "rg_ComplexNumber.h"
#include "DynamicDisk.h"
#include "CountStatistics.h"
#include <vector>
#include <list>
using namespace std;


namespace V {
namespace GeometryTier {



class EdgeFlipTimeCalculator
{
public:
    static rg_Polynomial      make_polynomial_equation_for_vedge_flipping(
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);

    static pair<bool, double> solve_polynomial_equation_for_vedge_flipping( 
        const PolynomialSolverType& solverType,     
        const double& mostRecentlyEventOccurringTime, 
        const rg_Polynomial& polynomialEq, 
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
        const bool& isThisVEdgeFlipping);

    static pair<bool, double> solve_polynomial_equation_for_vedge_flipping_by_LABSRC(
        const rg_Polynomial& polynomialEq, 
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
        const bool& isThisVEdgeFlipping);
   
    static pair<bool, double> solve_polynomial_equation_for_vedge_flipping_by_JENKINS_TRAUB(
        const rg_Polynomial& polynomialEq, 
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
        const bool& isThisVEdgeFlipping);
   
    static pair<bool, double> solve_polynomial_equation_for_vedge_flipping_by_STURMS_SEQUENCE( 
        const rg_Polynomial& polynomialEq, 
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
        const bool& isThisVEdgeFlipping, 
        const double& upperBoundOfTimeInterval = DBL_MAX);

    static void fill_elements_of_all_polynomial_XY_matrix(
        EdgeFlippingFormulaGeneratingType mode, 
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
        rg_Polynomial(*XYMat)[3]);
   
    static void fill_elements_of_all_polynomial_matrices(
        EdgeFlippingFormulaGeneratingType mode, 
        const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
        rg_Polynomial(*XMat)[3], 
        rg_Polynomial(*YMat)[3], 
        rg_Polynomial(*XYMat)[3]);
    
    static void fill_elements_for_polynomial_matrix(
                    const vector<DynamicDisk*>& movingDiskSet,
                    rg_Polynomial& x1_x4, rg_Polynomial& x2_x4, rg_Polynomial& x3_x4,
                    rg_Polynomial& y1_y4, rg_Polynomial& y2_y4, rg_Polynomial& y3_y4,
                    rg_Polynomial& r1_r4, rg_Polynomial& r2_r4, rg_Polynomial& r3_r4,
                    rg_Polynomial& p1, rg_Polynomial& p2, rg_Polynomial& p3);

    static rg_Polynomial make_determinant_polynomial(rg_Polynomial(*matrix)[3]);
    
    static void   find_positive_real_roots(rg_ComplexNumber* roots, const int& numOfRoots, list<double>& realRoots);
    
    static double find_min_absolute_real_root(rg_ComplexNumber* roots, const int& numOfRoots);
    
    static void   erase_numerically_zero_root_from_positive_real_roots(list<double>& realRoots, const double& minAbsoluteRealNumber);
    
    static bool   four_disks_are_circumscribing(const double& flippingTime, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);
    
    static bool   all_four_disks_have_same_radius(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);

    //For STURMS_SEQUENCE solver which is not good 
    static double do_false_position_method(const rg_Polynomial& polynomial, 
                                           const double& lowerBound,
                                           const double& upperBound,
                                           const int& maxIteration,
                                           const double& tolerance);

    static int find_number_of_roots_in_interval_using_sturms_chain(
        const list<rg_Polynomial>& sturmsChain, 
        const double& intervalBegin, 
        const double& intervalEnd);
   
    static list<rg_Polynomial> make_sturms_chain(const rg_Polynomial& polynomial);
    
    static pair<rg_Polynomial, rg_Polynomial> divide_polynomial(
        const rg_Polynomial& dividend, 
        const rg_Polynomial& divider);

    static int count_sign_transition_in_sturms_chain(
        const list<rg_Polynomial>& sturmsChain, 
        const double& valueForEvaluation);
   
    static double find_lower_bound_for_flipping_VEdge(
        const rg_Polynomial& polynomial, 
        const list<rg_Polynomial>& sturmsChain, 
        const double& lowerInterval, 
        const double& upperInterval);



#ifdef WRITING_POLYNOMIAL_RELATED_INFORMATION
    static int    numOfWriting;
    static string outputFolderPath;

    static const int    LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY   ;
    static const string polynomialRootsFileNameWithExtension      ;
    static const string polynomialCoefficientFileNameWithExtension;
    static const string movingDiskFileNameWithExtension           ;
    static const string realRootFileNameWithExtension             ;
    static const string realRootEvaluationFileNameWithExtension   ;
    static const string VEdgeFlippingStatusFileNameWithExtension                  ;

    //For getting disk and polynomial information during simulation by writing file 
    static CountStatistics                     m_CountStatisticsOfRoots;
    static vector<double>                      m_SimulationTimesAtMakingPolynomial;
    static vector<vector<DynamicDisk>>         m_QuadrupletsOfMovingDisks;
    static vector<rg_Polynomial>               m_Polynomials;
    static vector<char>                        m_FlippingStatusForVEdge;
    static vector<vector<double>>              m_RealRoots;
    
    static void  ___count_polynomial_root_type(const RootTypeOfPolynomial& rootType);

    static void  ___check_roots_of_polynomial_for_statistics(
                    const int& numOfRoots, 
                    rg_ComplexNumber* roots, 
                    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);

    static void ___store_current_time(const double& currentTime);

    static void ___store_this_polynomial_related_info(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
                                               const rg_Polynomial& polynomial,
                                               const rg_ComplexNumber* roots,
                                               const char& isThisFlippingVEdge);

    static void ___write_polynomial_related_info_in_temporal_storage();

private:
    static void ___write_headers_for_each();
    static void ___write_polynomial_roots_info();
    static void ___write_dynamic_disks_info();
    static void ___write_polynomial_coefficient_info();
    static void ___write_real_roots_info();
    static void ___write_roots_evaluation_info();
    static void ___write_VEdge_flipping_status();
    static char ___change_bool_type_of_flipping_status_to_char(const bool& type);
    static void ___clear_all_polynomial_related_info();

#endif //WRITING_POLYNOMIAL_RELATED_INFORMATION
};



}
}

#endif
