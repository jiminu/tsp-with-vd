#ifndef _RANDOM_DISK_GENERATOR_
#define _RANDOM_DISK_GENERATOR_

#include "Disk.h"
#include "DynamicDisk.h"
#include "BucketForDisk.h"
#include <list>
using namespace std;

class RandomDiskGenerator
{
public: 
    static int DEFAULT_CONTAINER_ID;

    //Function1. Disk generation
    static void generate_random_disks_in_circle_container(const int& numDisks,
                                                          const double& minRadius,
                                                          const double& maxRadius,
                                                          const double& ratioOfSumOfDisksAreaVsContainerArea,
                                                          const bool& diskCanBeIntersected,
                                                          list<Disk>& generatedDisks,
                                                          const rg_Point2D& containerCenterPt = rg_Point2D(0.0,0.0),
                                                          const bool& containCircleContainerInGeneratedDiskList = false,
                                                          const int& seedNumber = 0);

    static void generate_random_disks_in_circle_container(const int& numDisks,
                                                          const double& minRadius,
                                                          const double& maxRadius,
                                                          const bool& diskCanBeIntersected,
                                                          const rg_Circle2D& container,
                                                          list<Disk>& generatedDisks,
                                                          const bool& containCircleContainerInGeneratedDiskList = false,
                                                          const int& seedNumber = 0);


    static void generate_random_moving_disks_in_circle_container(const int& numDisks,
                                                                 const double& minRadius,
                                                                 const double& maxRadius,
                                                                 const double& minMagnitudeOfVelocityVector,
                                                                 const double& maxMagnitudeOfVelocityVector,
                                                                 const double& ratioOfSumOfDisksAreaVsContainerArea,
                                                                 const bool& diskCanBeIntersected,
                                                                 list<DynamicDisk>& generatedMovingDisks,
                                                                 const rg_Point2D& containerCenterPt = rg_Point2D(0.0,0.0),
                                                                 const bool& containCircleContainerInGeneratedDiskList = false,
                                                                 const int& seedNumber = 0);

    static void generate_random_moving_disks_in_circle_container(const int& numDisks,
                                                                 const double& minRadius,
                                                                 const double& maxRadius,
                                                                 const double& minMagnitudeOfVelocityVector,
                                                                 const double& maxMagnitudeOfVelocityVector,
                                                                 const bool& diskCanBeIntersected,
                                                                 const rg_Circle2D& container,
                                                                 list<DynamicDisk>& generatedMovingDisks,
                                                                 const bool& containCircleContainerInGeneratedDiskList = false,
                                                                 const int& seedNumber = 0);


    static void generate_random_disks_in_bounding_box(const int& numDisks,
                                                      const double& minRadius,
                                                      const double& maxRadius,
                                                      const double& ratioOfSumOfDisksAreaVsContainerArea,
                                                      const bool& diskCanBeIntersected,
                                                      list<Disk>& generatedDisks,
                                                      const int& seedNumber = 0);

    static void generate_random_moving_disks_in_bounding_box(const int& numDisks,
                                                             const double& minRadius,
                                                             const double& maxRadius,
                                                             const double& minMagnitudeOfVelocityVector,
                                                             const double& maxMagnitudeOfVelocityVector,
                                                             const double& ratioOfSumOfDisksAreaVsContainerArea,
                                                             const bool& diskCanBeIntersected,
                                                             list<DynamicDisk>& generatedMovingDisks,
                                                             const int& seedNumber = 0);

    //Function2. File writing

private:
    static void generate_random_disks_in_circle_container_centered_at_origin(const int& numDisks,
                                                                             const double& minRadius,
                                                                             const double& maxRadius,
                                                                             const double& containerRadius,
                                                                             const bool& diskCanBeIntersected,
                                                                             BucketForDisk& bucketForDisk,
                                                                             list<Disk>& generatedDisks);


    static double compute_radius_of_bounding_circle(const int& numDisks,
                                                    const double& minRadius,
                                                    const double& maxRadius,
                                                    const double& ratioOfSumOfDisksAreaVsContainerArea);
    
    static double compute_length_of_one_side_of_square_bounding_box(const int& numOfDisks, 
                                                                    const double& minRadius,
                                                                    const double& maxRadius,
                                                                    const double& ratioOfSumOfDisksAreaVsContainerArea);

    static bool this_disk_is_in_circle_container(const Disk& disk, const rg_Circle2D& circleContainer);
    
    static bool this_disk_is_in_bounding_box(const Disk& disk, const rg_Point2D& minPt, const rg_Point2D& maxPt);

    static bool this_is_included_in_other_disks_vice_versa(const Disk& disk, const BucketForDisk& bucketForDisk);

    static bool this_disk_is_intersected_with_other_disks(const Disk& disk, const BucketForDisk& bucketForDisk);

    static bool these_two_disks_are_intersected(const rg_Circle2D& circle1, const rg_Circle2D& circle2, const double& res);

    static void allocate_random_velocity_vector_to_disks(
        const double& minMagnitudeOfVelocityVector,const double& maxMagnitudeOfVelocityVector, 
        const list<Disk>& randomDisks, list<DynamicDisk>& randomMovingDisks);

    static void translate_disks(list<Disk>& generatedDisks, const rg_Point2D& translationIncrement);
};



#endif // 
