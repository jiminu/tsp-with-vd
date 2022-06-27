#ifndef _COLLISION_HANDLING_TYPE_
#define _COLLISION_HANDLING_TYPE_  

#include "DefaultValuesForDVD2.h"
#include "DynamicDisk.h"
#include <vector>
#include <random>
#include <ctime>
using namespace std;


namespace V {
namespace GeometryTier {


class DiskCollisionHandler
{
public:
	enum CollisionType { PHYSICAL_COLLISION, 
					     INITIALLY_MOVED_OR_FIXED_COLLISION, 
						 AVOIDANCE /*do not guarantee correct sol.*/, 
					     PHYSICAL_COLLISION_N_RANDOM_SPEED,
					     STEADY};



private:
	CollisionType  m_CollisionType;
	double         m_CoefficientOfRestitution;

public:
	DiskCollisionHandler();
	DiskCollisionHandler(const CollisionType& collisionType, const double& COR);
	DiskCollisionHandler(const DiskCollisionHandler& collisionHandlingType);
	DiskCollisionHandler& operator= (const DiskCollisionHandler& collisionHandlingType);

	/************************************************************************/
	/* Getter, Setter                                                      */
	/************************************************************************/
    inline double         get_coefficient_of_restitution()                       const;
    inline CollisionType  get_collision_type()                                   const;

	inline void  set_coefficient_of_restitution(const double& COR);
	inline void  set_collision_type(const CollisionType& collisionType);


    /************************************************************************/
    /* Main                                                                 */
    /************************************************************************/
	void compute_velocity_vectors_of_disks_after_collided(
				const DynamicDisk* movingDisk1,
				const DynamicDisk* movingDisk2,
				rg_Point2D& outputVector1,
				rg_Point2D& outputVector2) const;

    void compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION(
		const DynamicDisk* movingDisk1, 
		const DynamicDisk* movingDisk2, 
		rg_Point2D& outputVector1, 
		rg_Point2D& outputVector2) const;
    
	void compute_velocity_vectors_of_disks_after_collided_for_INTIALLY_MOVED_OR_FIXED_COLLISION(
		const DynamicDisk* movingDisk1, 
		const DynamicDisk* movingDisk2, 
		rg_Point2D& outputVector1, 
		rg_Point2D& outputVector2) const;
    
	void compute_velocity_vectors_of_disks_after_collided_for_AVOIDANCE(
		const DynamicDisk* movingDisk1, 
		const DynamicDisk* movingDisk2, 
		rg_Point2D& outputVector1, 
		rg_Point2D& outputVector2) const;
    
	void compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION_N_RANDOM_SPEED(
		const DynamicDisk* movingDisk1, 
		const DynamicDisk* movingDisk2, 
		rg_Point2D& outputVector1, 
		rg_Point2D& outputVector2) const;

    inline void compute_velocity_vectors_of_disks_after_collided_STEADY(
        const DynamicDisk* movingDisk1,
        const DynamicDisk* movingDisk2,
        rg_Point2D& outputVector1,
        rg_Point2D& outputVector2) const;

private:
    void find_optimal_velocity_vectors_for_PHYSICAL_COLLISION_N_RANDOM_SPEED(
		const DynamicDisk* movingDisk1, 
		const DynamicDisk* movingDisk2, 
		rg_Point2D& outputVector1, 
		rg_Point2D& outputVector2) const;

    inline void set_random_seed_number_for_simulation() const;

    inline double get_randomly_generated_zero_to_one() const;

	void copy_from(const DiskCollisionHandler& collisionHandlingType);
};


inline void DiskCollisionHandler::set_random_seed_number_for_simulation() const
{
    srand(time(NULL));
}

inline double DiskCollisionHandler::get_randomly_generated_zero_to_one() const
{
    return (double)(rand()) / (double)(RAND_MAX);
}

inline double DiskCollisionHandler::get_coefficient_of_restitution() const
{
	return m_CoefficientOfRestitution;
}

inline DiskCollisionHandler::CollisionType DiskCollisionHandler::get_collision_type() const
{
	return m_CollisionType;
}


inline void DiskCollisionHandler::set_coefficient_of_restitution(const double& COR)
{
	m_CoefficientOfRestitution = COR;
}


inline void DiskCollisionHandler::set_collision_type(const CollisionType& collisionType)
{
	m_CollisionType = collisionType;
}


inline void DiskCollisionHandler::compute_velocity_vectors_of_disks_after_collided_STEADY(
    const DynamicDisk* movingDisk1,
    const DynamicDisk* movingDisk2,
    rg_Point2D& outputVector1,
    rg_Point2D& outputVector2) const
{
    outputVector1 = movingDisk1->getVelocityVector();
    outputVector2 = movingDisk2->getVelocityVector();
}



}
}


#endif 