#include "RandomDiskGenerator.h"
#include <cstdlib>
using namespace std;

int RandomDiskGenerator::DEFAULT_CONTAINER_ID = -1;

void RandomDiskGenerator::generate_random_disks_in_circle_container(
    const int& numDisks, const double& minRadius, const double& maxRadius, 
    const double& ratioOfSumOfDisksAreaVsContainerArea, const bool& diskCanBeIntersected, 
    list<Disk>& generatedDisks, const rg_Point2D& containerCenterPt /*= rg_Point2D(0.0,0.0)*/, const bool& containCircleContainerInGeneratedDiskList /* = false */, const int& seedNumber /*= 0*/)
{
    srand(seedNumber);

    double containerRadius = compute_radius_of_bounding_circle(
                    numDisks, minRadius, maxRadius, ratioOfSumOfDisksAreaVsContainerArea);

    BucketForDisk bucketForDisk;
    rg_Point2D maxPtOfBoundingBox(containerRadius, containerRadius);
    rg_Point2D minPtOfBoundingBox = -maxPtOfBoundingBox;
    bucketForDisk.setBoundingBox(maxPtOfBoundingBox, minPtOfBoundingBox);
    bucketForDisk.setBucketCellSize(maxRadius*2.0 + 0.5);

    generate_random_disks_in_circle_container_centered_at_origin(
        numDisks, minRadius, maxRadius, containerRadius, 
        diskCanBeIntersected, bucketForDisk, generatedDisks);

    translate_disks(generatedDisks, containerCenterPt);

    if (containCircleContainerInGeneratedDiskList)
    {
        Disk circleContainer(DEFAULT_CONTAINER_ID, containerCenterPt, containerRadius);
        generatedDisks.push_front(circleContainer);
    }
}



void RandomDiskGenerator::generate_random_disks_in_circle_container_centered_at_origin(const int& numDisks, const double& minRadius, const double& maxRadius, const double& containerRadius, const bool& diskCanBeIntersected, BucketForDisk& bucketForDisk, list<Disk>& generatedDisks)
{
    rg_Circle2D circleContainer(0.0, 0.0, containerRadius);

    double intervalOfRadius            = maxRadius - minRadius;
    double currDiskRadius              = 0.0;
    bool   isRadiusAlreadySet          = false;
    int    numOfGeneratedDisksUntilNow = 0;

    while (numOfGeneratedDisksUntilNow < numDisks)
    {
        double x = ((double)rand() / (double)RAND_MAX - 0.5) * containerRadius * 2.0;
        double y = ((double)rand() / (double)RAND_MAX - 0.5) * containerRadius * 2.0;

        if (!isRadiusAlreadySet)
        {
            currDiskRadius = minRadius + intervalOfRadius * ((double)rand() / (double)RAND_MAX);
            isRadiusAlreadySet = true;
        }

        Disk currDisk(numOfGeneratedDisksUntilNow + 1, x, y, currDiskRadius);


        if (this_disk_is_in_circle_container(currDisk, circleContainer))
        {
            if (diskCanBeIntersected)
            {
                if (!this_is_included_in_other_disks_vice_versa(currDisk, bucketForDisk))
                {
                    generatedDisks.push_back(currDisk);
                    bucketForDisk.addDiskIntoBucket(&(generatedDisks.back()));
                    isRadiusAlreadySet = false;
                    ++numOfGeneratedDisksUntilNow;
                }
            }
            else
            {
                if (!this_disk_is_intersected_with_other_disks(currDisk, bucketForDisk))
                {
                    generatedDisks.push_back(currDisk);
                    bucketForDisk.addDiskIntoBucket(&(generatedDisks.back()));
                    isRadiusAlreadySet = false;
                    ++numOfGeneratedDisksUntilNow;
                }
            }
        }
    }
}




void RandomDiskGenerator::generate_random_disks_in_circle_container(
    const int& numDisks, const double& minRadius, const double& maxRadius, 
    const bool& diskCanBeIntersected, const rg_Circle2D& container, list<Disk>& generatedDisks, 
    const bool& containCircleContainerInGeneratedDiskList /*= false*/, const int& seedNumber /*= 0*/)
{
    srand(seedNumber);

    double containerRadius = container.getRadius();

    BucketForDisk bucketForDisk;
    rg_Point2D maxPtOfBoundingBox(containerRadius, containerRadius);
    rg_Point2D minPtOfBoundingBox = -maxPtOfBoundingBox;
    bucketForDisk.setBoundingBox(maxPtOfBoundingBox, minPtOfBoundingBox);
    bucketForDisk.setBucketCellSize(maxRadius * 2.0 + 0.5);

    generate_random_disks_in_circle_container_centered_at_origin(
        numDisks, minRadius, maxRadius, containerRadius,
        diskCanBeIntersected, bucketForDisk, generatedDisks);

    translate_disks(generatedDisks, container.getCenterPt());

    if (containCircleContainerInGeneratedDiskList)
    {
        Disk circleContainer(DEFAULT_CONTAINER_ID, container);
        generatedDisks.push_front(circleContainer);
    }
}


void RandomDiskGenerator::generate_random_moving_disks_in_circle_container(
    const int& numDisks, const double& minRadius, const double& maxRadius, 
    const double& minMagnitudeOfVelocityVector, const double& maxMagnitudeOfVelocityVector, 
    const double& ratioOfSumOfDisksAreaVsContainerArea, const bool& diskCanBeIntersected, 
    list<DynamicDisk>& generatedMovingDisks, 
    const rg_Point2D& containerCenterPt /*= rg_Point2D(0.0,0.0)*/, 
    const bool& containCircleContainerInGeneratedDiskList /* = false */, 
    const int& seedNumber /*= 0*/)
{
    list<Disk> randomDisks;
    generate_random_disks_in_circle_container(
        numDisks, minRadius, maxRadius, ratioOfSumOfDisksAreaVsContainerArea, 
        diskCanBeIntersected, randomDisks, containerCenterPt, containCircleContainerInGeneratedDiskList, seedNumber);

    allocate_random_velocity_vector_to_disks(minMagnitudeOfVelocityVector, maxMagnitudeOfVelocityVector,
        randomDisks, generatedMovingDisks);
}


void RandomDiskGenerator::generate_random_moving_disks_in_circle_container(
    const int& numDisks, const double& minRadius, const double& maxRadius, 
    const double& minMagnitudeOfVelocityVector, const double& maxMagnitudeOfVelocityVector, 
    const bool& diskCanBeIntersected, const rg_Circle2D& container, 
    list<DynamicDisk>& generatedMovingDisks, 
    const bool& containCircleContainerInGeneratedDiskList /*= false*/, const int& seedNumber /*= 0*/)
{
    list<Disk> randomDisks;
    generate_random_disks_in_circle_container(
        numDisks, minRadius, maxRadius, diskCanBeIntersected, container, randomDisks, 
        containCircleContainerInGeneratedDiskList, seedNumber);

    allocate_random_velocity_vector_to_disks(minMagnitudeOfVelocityVector, maxMagnitudeOfVelocityVector,
        randomDisks, generatedMovingDisks);
}


void RandomDiskGenerator::generate_random_disks_in_bounding_box(
    const int& numDisks, const double& minRadius, const double& maxRadius, 
    const double& ratioOfSumOfDisksAreaVsContainerArea, const bool& diskCanBeIntersected, 
    list<Disk>& generatedDisks, const int& seedNumber /*= 0*/)
{
    srand(seedNumber);
 
    double lengthOfOneSideOfBoundingBox = compute_length_of_one_side_of_square_bounding_box(
        numDisks, minRadius, maxRadius, ratioOfSumOfDisksAreaVsContainerArea);

    BucketForDisk bucketForDisk;
    rg_Point2D maxPtOfBoundingBox(lengthOfOneSideOfBoundingBox / 2.0, lengthOfOneSideOfBoundingBox / 2.0);
    rg_Point2D minPtOfBoundingBox = -maxPtOfBoundingBox;
    bucketForDisk.setBoundingBox(maxPtOfBoundingBox, minPtOfBoundingBox);
    bucketForDisk.setBucketCellSize(maxRadius*2.0 + 0.5);


    double intervalOfRadius   = maxRadius - minRadius;
    double currDiskRadius     = 0.0;
    bool   isRadiusAlreadySet = false;
    int    numOfDisksInVector = 0;

    while (numOfDisksInVector < numDisks)
    {
        double x = ((double)rand() / (double)RAND_MAX  - 0.5) * lengthOfOneSideOfBoundingBox;
        double y = ((double)rand() / (double)RAND_MAX  - 0.5) * lengthOfOneSideOfBoundingBox;

        if (!isRadiusAlreadySet)
        {
            currDiskRadius     = minRadius + intervalOfRadius * ((double)rand() / (double)RAND_MAX);
            isRadiusAlreadySet = true;
        }

        Disk currDisk(numOfDisksInVector + 1, x, y, currDiskRadius);

        if (this_disk_is_in_bounding_box(currDisk, minPtOfBoundingBox, maxPtOfBoundingBox))
        {
            if (diskCanBeIntersected)
            {
                if (!this_is_included_in_other_disks_vice_versa(currDisk, bucketForDisk))
                {
                    generatedDisks.push_back(currDisk);
                    bucketForDisk.addDiskIntoBucket(&(generatedDisks.back()));
                    isRadiusAlreadySet                 = false;
                    ++numOfDisksInVector;
                }
            }
            else
            {
                if (!this_disk_is_intersected_with_other_disks(currDisk, bucketForDisk))
                {
                    generatedDisks.push_back(currDisk);
                    bucketForDisk.addDiskIntoBucket(&(generatedDisks.back()));
                    isRadiusAlreadySet                 = false;
                    ++numOfDisksInVector;
                }
            }
        }
    }

}


void RandomDiskGenerator::generate_random_moving_disks_in_bounding_box(
    const int& numDisks, const double& minRadius, const double& maxRadius, 
    const double& minMagnitudeOfVelocityVector, const double& maxMagnitudeOfVelocityVector, 
    const double& ratioOfSumOfDisksAreaVsContainerArea, const bool& diskCanBeIntersected, 
    list<DynamicDisk>& generatedMovingDisks, const int& seedNumber /*= 0*/)
{
    list<Disk> randomDisks;
    generate_random_disks_in_bounding_box(
        numDisks, minRadius, maxRadius, ratioOfSumOfDisksAreaVsContainerArea,
        diskCanBeIntersected, randomDisks, seedNumber);

    allocate_random_velocity_vector_to_disks(
        minMagnitudeOfVelocityVector, maxMagnitudeOfVelocityVector,
        randomDisks, generatedMovingDisks);
}



double RandomDiskGenerator::compute_radius_of_bounding_circle(const int& numDisks, const double& minRadius, const double& maxRadius, const double& ratioOfSumOfDisksAreaVsContainerArea)
{
    double sumOfExpectedAreaOfDiskDevidedPI = 0.0;
    double squareOfRadiusOfBoundungCircle   = 0.0;
    double radiusOfOuterBoundingCircle      = 0.0;
    double intervalOfRadius = maxRadius - minRadius;

    for (int i = 1; i <= numDisks; i++)
    {
        double currRadius = minRadius + intervalOfRadius * ((double)i / (double)(numDisks + 1));
        double currAreaDevidedPI = currRadius * currRadius;

        sumOfExpectedAreaOfDiskDevidedPI += currAreaDevidedPI;
    }


    squareOfRadiusOfBoundungCircle = sumOfExpectedAreaOfDiskDevidedPI / ratioOfSumOfDisksAreaVsContainerArea;
    radiusOfOuterBoundingCircle    = sqrt(squareOfRadiusOfBoundungCircle);

    return radiusOfOuterBoundingCircle;
}


double RandomDiskGenerator::compute_length_of_one_side_of_square_bounding_box(
    const int& numOfDisks, const double& minRadius, const double& maxRadius, 
    const double& ratioOfSumOfDisksAreaVsContainerArea)
{
    double expectedTotalAreaOfDisksDevidedPI = 0.0;
    double intervalOfRadius                  = maxRadius - minRadius;

    for (int i = 1; i <= numOfDisks; i++)
    {
        double currRadius            = minRadius + intervalOfRadius * ((double)i / (double)(numOfDisks + 1));
        double currDiskAreaDevidedPI = currRadius * currRadius;

        expectedTotalAreaOfDisksDevidedPI += currDiskAreaDevidedPI;
    }

    double expectedTotalAreaOfDisks = expectedTotalAreaOfDisksDevidedPI * rg_PI;
    double boundingBoxArea          = expectedTotalAreaOfDisks / ratioOfSumOfDisksAreaVsContainerArea;
    double lengthOfOneSideOfBoundingBox = sqrt(boundingBoxArea);

    return lengthOfOneSideOfBoundingBox;
}


bool RandomDiskGenerator::this_disk_is_in_circle_container(const Disk& disk, const rg_Circle2D& circleContainer)
{
    if (disk.isIncludedIn(circleContainer))
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool RandomDiskGenerator::this_disk_is_in_bounding_box(const Disk& disk, const rg_Point2D& minPt, const rg_Point2D& maxPt)
{
    if(  disk.getCenterPt().getX() - disk.getRadius() > minPt.getX() 
      && disk.getCenterPt().getX() + disk.getRadius() < maxPt.getX()
      && disk.getCenterPt().getY() - disk.getRadius() > minPt.getY()
      && disk.getCenterPt().getY() + disk.getRadius() < maxPt.getY())
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool RandomDiskGenerator::this_is_included_in_other_disks_vice_versa(const Disk& disk, const BucketForDisk& bucketForDisk)
{
    bool thisDiskIsIncludedInObstacleDisksViceVersa = false;

    BucketCellForDisks* bucketCellIncludingPoint = bucketForDisk.getBucketCellIncludingPoint(disk.getCenterPt());

    list<BucketCellForDisks*> maxNineCellsAroundTheCellIncludingPoint;
    bucketForDisk.getMax9BucketCellsIncludingNAroundThisCell(bucketCellIncludingPoint, maxNineCellsAroundTheCellIncludingPoint);

    for(list<BucketCellForDisks*>::const_iterator it_BucketCell = maxNineCellsAroundTheCellIncludingPoint.begin();
        it_BucketCell != maxNineCellsAroundTheCellIncludingPoint.end();
        ++it_BucketCell)
    {
        BucketCellForDisks* currBucketCell                           = *it_BucketCell;
        const list<rg_Circle2D*>& containedObstacleDisksInBucketCell = currBucketCell->get_contained_disks();

        for (list<rg_Circle2D*>::const_iterator it_Obstacle = containedObstacleDisksInBucketCell.begin();
            it_Obstacle != containedObstacleDisksInBucketCell.end();
            it_Obstacle++)
        {
            const rg_Circle2D* currObstacle = *it_Obstacle;

            if (currObstacle->isIncludedIn(disk) || disk.isIncludedIn(*currObstacle))
            {
                thisDiskIsIncludedInObstacleDisksViceVersa = true;
                break;
            }
        }

        if (thisDiskIsIncludedInObstacleDisksViceVersa)
        {
            break;
        }
    }

    return thisDiskIsIncludedInObstacleDisksViceVersa;
}


bool RandomDiskGenerator::this_disk_is_intersected_with_other_disks(const Disk& disk, const BucketForDisk& bucketForDisk)
{
    bool thisDiskIsIntersectedWithOtherDisks = false;

    BucketCellForDisks* bucketCellIncludingPoint = bucketForDisk.getBucketCellIncludingPoint(disk.getCenterPt());

    list<BucketCellForDisks*> maxNineCellsAroundTheCellIncludingPoint;
    bucketForDisk.getMax9BucketCellsIncludingNAroundThisCell(bucketCellIncludingPoint, maxNineCellsAroundTheCellIncludingPoint);

    for(list<BucketCellForDisks*>::const_iterator it_BucketCell = maxNineCellsAroundTheCellIncludingPoint.begin();
        it_BucketCell != maxNineCellsAroundTheCellIncludingPoint.end();
        ++it_BucketCell)
    {
        BucketCellForDisks* currBucketCell                           = *it_BucketCell;
        const list<rg_Circle2D*>& containedObstacleDisksInBucketCell = currBucketCell->get_contained_disks();

        for (list<rg_Circle2D*>::const_iterator it_Obstacle = containedObstacleDisksInBucketCell.begin();
            it_Obstacle != containedObstacleDisksInBucketCell.end();
            it_Obstacle++)
        {
            const rg_Circle2D* currObstacle = *it_Obstacle;

            if (these_two_disks_are_intersected(disk, *currObstacle, resNeg3))
            {
                thisDiskIsIntersectedWithOtherDisks = true;
                break;
            }
        }

        if (thisDiskIsIntersectedWithOtherDisks)
        {
            break;
        }
    }

    return thisDiskIsIntersectedWithOtherDisks;
}


bool RandomDiskGenerator::these_two_disks_are_intersected(const rg_Circle2D& circle1, const rg_Circle2D& circle2, const double& res)
{
    return circle1.isIntersectWith(circle2, res);
}


void RandomDiskGenerator::allocate_random_velocity_vector_to_disks(
    const double& minMagnitudeOfVelocityVector, 
    const double& maxMagnitudeOfVelocityVector, 
    const list<Disk>& randomDisks, 
    list<DynamicDisk>& randomMovingDisks)
{
    double intervalOfMagnitude = maxMagnitudeOfVelocityVector - minMagnitudeOfVelocityVector;

    for(list<Disk>::const_iterator it_Disk = randomDisks.begin();
        it_Disk != randomDisks.end();
        ++it_Disk)
    {
        const Disk& currDisk = *it_Disk;

        if (currDisk.getID() == DEFAULT_CONTAINER_ID)
        {
            randomMovingDisks.push_back(DynamicDisk(currDisk.getID(), currDisk, rg_Point2D(0.0,0.0), 0.0));
        }
        else
        {
            double theta         = 2 * rg_PI * ((double)rand() / (double)RAND_MAX);
            double magnitude     = minMagnitudeOfVelocityVector + intervalOfMagnitude * ((double)rand() / (double)RAND_MAX);

            rg_Point2D velocityVector(cos(theta)* magnitude, sin(theta)*magnitude);
            randomMovingDisks.push_back(DynamicDisk(currDisk.getID(), currDisk, velocityVector, 0.0));
        }
    }
}


void RandomDiskGenerator::translate_disks(list<Disk>& generatedDisks, const rg_Point2D& translationIncrement)
{
    for (list<Disk>::iterator it_Disk = generatedDisks.begin();
        it_Disk != generatedDisks.end();
        ++it_Disk)
    {
        Disk& currDisk = *it_Disk;
        rg_Point2D currCenterPt = currDisk.getCenterPt();
        currDisk.setCenterPt(currCenterPt + translationIncrement);
    }
}