#include "DiskCollisionTimeCalculator.h"
#include "SimpleTypesForDVD2.h"
#include "DefaultValuesForDVD2.h"
using namespace V::GeometryTier;

std::pair<bool, double> DiskCollisionTimeCalculator::find_possible_collision_time_of_disks_after(
    const double& mostRecentlyEventOccurringTime, const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2)
{
    pair<bool, double> nextCollisionTime(false, DBL_MAX);

    DVD_StatusOfTwoDisks statusOfTwoDisks = find_state_of_two_disks(const_cast<DynamicDisk*>(movingDisk1), 
                                                                    const_cast<DynamicDisk*>(movingDisk2));
    switch (statusOfTwoDisks)
    {
        case DISK_N_DISK:
        {
            nextCollisionTime = find_possible_collision_time_of_disks_DISK_N_DISK(
                movingDisk1, movingDisk2);

            break;
        }

        case DISK_N_CONTAINER:
        {
            if (movingDisk1->getCircle().isIncludedIn(movingDisk2->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
            {
                 nextCollisionTime = find_possible_collision_time_of_disks_DISK_N_CONTAINER_INCLUDING_THE_DISK(
                     movingDisk1, movingDisk2);
            }
            else
            {
                nextCollisionTime = find_possible_collision_time_of_disks_DISK_N_DISK(
                    movingDisk1, movingDisk2);
            }

            break;
        }
        
        
        case CONTAINER_N_DISK:
        {
            if (movingDisk2->getCircle().isIncludedIn(movingDisk1->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
            {
                nextCollisionTime = find_possible_collision_time_of_disks_DISK_N_CONTAINER_INCLUDING_THE_DISK(
                    movingDisk2, movingDisk1);
            }
            else
            {
                nextCollisionTime = find_possible_collision_time_of_disks_DISK_N_DISK(
                    movingDisk1, movingDisk2);
            }

            break;
        }

        case CONTAINER_N_CONTAINER:
        {
            nextCollisionTime = pair<bool, double>(false, DBL_MAX);
            break;
        }
    }

    if (nextCollisionTime.first)
    {
        nextCollisionTime.second += mostRecentlyEventOccurringTime;
    }

    return nextCollisionTime;
}


std::pair<bool, double> DiskCollisionTimeCalculator::find_possible_collision_time_of_disks_DISK_N_DISK(
    const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2)
{
    rg_Point2D c1_c2 = movingDisk1->getCircle().getCenterPt() - movingDisk2->getCircle().getCenterPt();
    rg_Point2D v1_v2 = movingDisk1->getVelocityVector()       - movingDisk2->getVelocityVector();

    double r1 = movingDisk1->getCircle().getRadius() + movingDisk1->getOffset();
    double r2 = movingDisk2->getCircle().getRadius() + movingDisk2->getOffset();

    double A = 0.0;
    double B = 0.0;
    double collisionTime = 0.0;

    pair<bool, double> nextCollisionTime(false, DBL_MAX);

    if (rg_ZERO(r1) && rg_ZERO(r2)) // if there are two points ( no radius size )
    {
        nextCollisionTime = pair<bool, double>(false, DBL_MAX);
    }
    else if (v1_v2.magnitude() < rg_MATH_RES) //2017.10.03 cysong
    {
        nextCollisionTime = pair<bool, double>(false, DBL_MAX);
    }
    else if (there_is_intersection_btw(rg_Circle2D(movingDisk1->getCircle().getCenterPt(), movingDisk1->getCircle().getRadius() + movingDisk1->getOffset()),
                                       rg_Circle2D(movingDisk2->getCircle().getCenterPt(), movingDisk2->getCircle().getRadius() + movingDisk2->getOffset()),
                                       TOLERANCE_OF_DISK_INTERSECTION))
    {
        //return pair<bool, double>(true, mostRecentlyEventOccurringTime);
        nextCollisionTime = pair<bool, double>(false, DBL_MAX);
    }
    else
    {
        nextCollisionTime = pair<bool, double>(false, DBL_MAX);

        //check1. necessary condition for the collision btw movingDisk1 and 2.
        if (two_disks_are_moving_to_opposite_side_relatively(movingDisk1, movingDisk2)) //are_two_particles_moving_to_opposite_side()
        {
            //check2. sufficient condition for the collision btw movingDisk1 and 2.
            if (two_disks_will_be_collided(movingDisk1, movingDisk2)) //will_two_particles_be_collided()
            {
                A = (v1_v2 % c1_c2);
                B = sqrt(pow(v1_v2 % c1_c2, 2) - v1_v2.magnitudeSquare()*(c1_c2.magnitudeSquare() - pow(r1 + r2, 2)));

                //[Start] ADDED at Oct 04, 2019 for numericall stability.
                double collisionTime = (-A - B) / v1_v2.magnitudeSquare();
                if (collisionTime < 0.0)
                {
                    collisionTime = 0.0;
                }
                //[End] ADDED at Oct 04, 2019 for numericall stability.
                nextCollisionTime = make_pair(true, collisionTime);
            }
        }
    }

    return nextCollisionTime;
}


std::pair<bool, double> DiskCollisionTimeCalculator::find_possible_collision_time_of_disks_DISK_N_CONTAINER_INCLUDING_THE_DISK(
    const DynamicDisk* movingDisk, const DynamicDisk* container)
{
    rg_Point2D c1_c2 = movingDisk->getCircle().getCenterPt() - container->getCircle().getCenterPt();
    rg_Point2D v1_v2 = movingDisk->getVelocityVector()       - container->getVelocityVector();

    double r1 = movingDisk->getCircle().getRadius() + movingDisk->getOffset();
    double r2 = container->getCircle().getRadius()  + container->getOffset();

    double A = 0.0;
    double B = 0.0;
    double collisionTime = 0.0;

    if (v1_v2.magnitude() < rg_MATH_RES)
    {
        return pair<bool, double>(false, DBL_MAX);
    }
    else
    {
        A = (v1_v2 % c1_c2);
        B = sqrt(pow(v1_v2 % c1_c2, 2) - v1_v2.magnitudeSquare()*(c1_c2.magnitudeSquare() - pow(r1 - r2, 2)));
        collisionTime = (-A + B) / v1_v2.magnitudeSquare();

        return pair<bool, double>(true, collisionTime);
    }
}



bool DiskCollisionTimeCalculator::two_disks_are_moving_to_opposite_side_relatively(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2)
{
    rg_Point2D c1_c2 = movingDisk1->getCircle().getCenterPt() - movingDisk2->getCircle().getCenterPt();
    rg_Point2D v1_v2 = movingDisk1->getVelocityVector() - movingDisk2->getVelocityVector();

    //check1. necessary condition for the collision btw movingDisk1 and 2.
    if (v1_v2 % c1_c2 < 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}



bool DiskCollisionTimeCalculator::two_disks_will_be_collided(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2)
{
    rg_Point2D c1_c2 = movingDisk1->getCircle().getCenterPt() - movingDisk2->getCircle().getCenterPt();
    rg_Point2D v1_v2 = movingDisk1->getVelocityVector() - movingDisk2->getVelocityVector();

    double r1 = movingDisk1->getCircle().getRadius() + movingDisk1->getOffset() * 0.5;
    double r2 = movingDisk2->getCircle().getRadius() + movingDisk2->getOffset() * 0.5;

    //double r1 = movingDisk1->getCircle().getRadius() + movingDisk1->getOffset();
    //double r2 = movingDisk2->getCircle().getRadius() + movingDisk2->getOffset();

    //check2. sufficient condition for the collision btw movingDisk1 and 2.
    if (pow(v1_v2 % c1_c2, 2) - (v1_v2 % v1_v2) * ((c1_c2 % c1_c2) - pow(r1 + r2, 2)) > TOLERANCE_OF_COLLISION_PREDICTION) //will_two_particles_be_collided()
    {
        return true;
    }
    else
    {
        return false;
    }
}

