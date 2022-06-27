#ifndef _DISK_COLLISION_TIME_CALCULATOR_
#define _DISK_COLLISION_TIME_CALCULATOR_

#include "DynamicDisk.h"

namespace V {
namespace GeometryTier {


class DiskCollisionTimeCalculator
{
public:
    static pair<bool, double> find_possible_collision_time_of_disks_after(
        const double& mostRecentlyEventOccurringTime, const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);

    static pair<bool, double> find_possible_collision_time_of_disks_DISK_N_DISK( const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    
    static pair<bool, double> find_possible_collision_time_of_disks_DISK_N_CONTAINER_INCLUDING_THE_DISK( const DynamicDisk* movingDisk, const DynamicDisk* container);

    //For checking whether two disks will be collided. 
    static inline bool two_disks_are_moving_to_opposite_side_relatively(
        const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    
    static inline bool two_disks_will_be_collided(
        const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
};


}
}
#endif