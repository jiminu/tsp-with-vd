#include "DiskCollisionHandler.h"
using namespace V::GeometryTier;

DiskCollisionHandler::DiskCollisionHandler()
{
	m_CollisionType            = PHYSICAL_COLLISION;
	m_CoefficientOfRestitution = 1.0;

    set_random_seed_number_for_simulation();
}


DiskCollisionHandler::DiskCollisionHandler(const DiskCollisionHandler& collisionHandlingType)
{
	copy_from(collisionHandlingType);
}


DiskCollisionHandler::DiskCollisionHandler(const CollisionType& collisionType, const double& COR)
{
	m_CollisionType            = collisionType;
	m_CoefficientOfRestitution = COR;
}


DiskCollisionHandler& DiskCollisionHandler::operator=(const DiskCollisionHandler& collisionHandlingType)
{
	if (this != &collisionHandlingType)
	{
		copy_from(collisionHandlingType);
	}

	return *this;
}


void DiskCollisionHandler::copy_from(const DiskCollisionHandler& collisionHandlingType)
{
	m_CollisionType            = collisionHandlingType.m_CollisionType;
	m_CoefficientOfRestitution = collisionHandlingType.m_CoefficientOfRestitution;
}


void DiskCollisionHandler::compute_velocity_vectors_of_disks_after_collided(
    const DynamicDisk* movingDisk1, 
    const DynamicDisk* movingDisk2, 
    rg_Point2D& outputVector1, 
    rg_Point2D& outputVector2) const
{
	switch (m_CollisionType)
    {
        case DiskCollisionHandler::PHYSICAL_COLLISION:
        {
            compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION(movingDisk1, movingDisk2, outputVector1, outputVector2);
            break;
        }

        case DiskCollisionHandler::INITIALLY_MOVED_OR_FIXED_COLLISION:
        {
            compute_velocity_vectors_of_disks_after_collided_for_INTIALLY_MOVED_OR_FIXED_COLLISION(movingDisk1, movingDisk2, outputVector1, outputVector2);
            break;
        }

        case DiskCollisionHandler::AVOIDANCE:
        {
            compute_velocity_vectors_of_disks_after_collided_for_AVOIDANCE(movingDisk1, movingDisk2, outputVector1, outputVector2);
            break;
        }

        case DiskCollisionHandler::PHYSICAL_COLLISION_N_RANDOM_SPEED:
        {
            compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION_N_RANDOM_SPEED(movingDisk1, movingDisk2, outputVector1, outputVector2);
            break;
        }

        case DiskCollisionHandler::STEADY:
        {
            compute_velocity_vectors_of_disks_after_collided_STEADY(movingDisk1, movingDisk2, outputVector1, outputVector2);
            break;
        }
    }
}


void DiskCollisionHandler::compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION(
    const DynamicDisk* movingDisk1, 
    const DynamicDisk* movingDisk2, 
    rg_Point2D& outputVector1, 
    rg_Point2D& outputVector2) const
{
    double e = m_CoefficientOfRestitution;

    rg_Point2D c1 = movingDisk1->getCircle().getCenterPt();
    rg_Point2D c2 = movingDisk2->getCircle().getCenterPt();

    rg_Point2D v1 = rg_Point2D(movingDisk1->getVelocityVectorX(), movingDisk1->getVelocityVectorY());
    rg_Point2D v2 = rg_Point2D(movingDisk2->getVelocityVectorX(), movingDisk2->getVelocityVectorY());

    double m1     = pow(movingDisk1->getCircle().getRadius(), 2);
    double m2     = pow(movingDisk2->getCircle().getRadius(), 2);

    if (movingDisk2->getCircle().isIncludedIn(movingDisk1->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        rg_Point2D basisVector((c2 - c1).getUnitVector());
        rg_Point2D newVector(v2 - 2 * (v2 % basisVector) * basisVector);

        outputVector1.setPoint(0.0, 0.0);
        outputVector2.setPoint(e * newVector);
    }
    else if (movingDisk1->getCircle().isIncludedIn(movingDisk2->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        rg_Point2D basisVector((c1 - c2).getUnitVector());
        rg_Point2D newVector(v1 - 2 * (v1 % basisVector) * basisVector);

        outputVector1.setPoint(e * newVector);
        outputVector2.setPoint(0.0, 0.0);
    }
    else
    {
        /*

            |A|*cost*B / |B| = (A dot B) / |B| * B


            solve equations

            known: m1,m2,v1,v2, unknown: u1, u2

            1. m1v1 + m2v2 = m1u1 + m2u2
            2. e = (u1 - u2) / (v2 - v1)

            u1 = v1 + (v2 - v1)*(e + 1) * m2 / (m1 + m2)
            u2 = v2 + (v1 - v2)*(e + 1) * m1 / (m1 + m2)

            */

        /*
        const double v1_X = v1.getX();
        const double v2_X = v2.getX();
        const double v1_Y = v1.getY();
        const double v2_Y = v2.getY();

        const rg_Point2D basisVector = (c2 - c1).getUnitVector();

        const double cosTheta = basisVector.getX();
        const double sinTheta = basisVector.getY();

        double v1_X_on_the_basis_vector_after_collision;
        double v2_X_on_the_basis_vector_after_collision;
        double v1_Y_on_the_basis_vector_after_collision;
        double v2_Y_on_the_basis_vector_after_collision;
        
        double v1_X_after_collision;
        double v2_X_after_collision;
        double v1_Y_after_collision;
        double v2_Y_after_collision;


        v1_X_on_the_basis_vector_after_collision = (m1 - e * m2)/(m1 + m2) * (v1_X * cosTheta + v1_Y * sinTheta) 
                                                    + (m2 + e * m2)/(m1 + m2) * (v2_X * cosTheta + v2_Y * sinTheta);
        v2_X_on_the_basis_vector_after_collision = (m1 + e * m1)/(m1 + m2) * (v1_X * cosTheta + v1_Y * sinTheta) 
                                                    + (m2 - e * m1)/(m1 + m2) * (v2_X * cosTheta + v2_Y * sinTheta);

        v1_Y_on_the_basis_vector_after_collision = v1_Y * cosTheta - v1_X * sinTheta;

        v2_Y_on_the_basis_vector_after_collision = v2_Y * cosTheta - v2_X * sinTheta;
         
        v1_X_after_collision = v1_X_on_the_basis_vector_after_collision * cosTheta 
                                - v1_Y_on_the_basis_vector_after_collision * sinTheta;
        
        v1_Y_after_collision = v1_X_on_the_basis_vector_after_collision * sinTheta 
                                + v1_Y_on_the_basis_vector_after_collision * cosTheta;
        
        v2_X_after_collision = v2_X_on_the_basis_vector_after_collision * cosTheta
                                - v2_Y_on_the_basis_vector_after_collision * sinTheta;
        
        v2_Y_after_collision = v2_X_on_the_basis_vector_after_collision * sinTheta 
                                + v2_Y_on_the_basis_vector_after_collision * cosTheta;
        
        outputVector1.setPoint(v1_X_after_collision, v1_Y_after_collision);
        outputVector2.setPoint(v2_X_after_collision, v2_Y_after_collision);
        */

        
        double m1 = pow(movingDisk1->getCircle().getRadius(), 2);
        double m2 = pow(movingDisk2->getCircle().getRadius(), 2);

        const rg_Point2D basisVector = (c2 - c1).getUnitVector();

        double vb1_x = (v1 % basisVector);
        double vb2_x = (v2 % basisVector);

        double vb1_x_afterCollision = vb1_x + (vb2_x - vb1_x)*(e + 1) * m2 / (m1 + m2);
        double vb2_x_afterCollision = vb2_x + (vb1_x - vb2_x)*(e + 1) * m1 / (m1 + m2);

        double cosTheta = basisVector.getX();
        double sinTheta = basisVector.getY();

        double vb1_y_afterCollision = v1.getY()*cosTheta - v1.getX()*sinTheta;
        double vb2_y_afterCollision = v2.getY()*cosTheta - v2.getX()*sinTheta;

        double v1_x_afterCollision = vb1_x_afterCollision * cosTheta - vb1_y_afterCollision * sinTheta;
        double v1_y_afterCollision = vb1_x_afterCollision * sinTheta + vb1_y_afterCollision * cosTheta;

        double v2_x_afterCollision = vb2_x_afterCollision * cosTheta - vb2_y_afterCollision * sinTheta;
        double v2_y_afterCollision = vb2_x_afterCollision * sinTheta + vb2_y_afterCollision * cosTheta;

        if (movingDisk1->this_disk_is_container())
        {
            outputVector1.setPoint(0.0, 0.0);
            outputVector2.setPoint(v2_x_afterCollision, v2_y_afterCollision);
            outputVector2 = outputVector2.getUnitVector() * v2.magnitude();
        }
        else if (movingDisk2->this_disk_is_container())
        {
            outputVector1.setPoint(v1_x_afterCollision, v1_y_afterCollision);
            outputVector1 = outputVector1.getUnitVector() * v1.magnitude();
            outputVector2.setPoint(0.0, 0.0);
        }
        else
        {
            outputVector1.setPoint(v1_x_afterCollision, v1_y_afterCollision);
            outputVector2.setPoint(v2_x_afterCollision, v2_y_afterCollision);
        }
    }
}


void DiskCollisionHandler::compute_velocity_vectors_of_disks_after_collided_for_INTIALLY_MOVED_OR_FIXED_COLLISION(
    const DynamicDisk* movingDisk1, 
    const DynamicDisk* movingDisk2, 
    rg_Point2D& outputVector1, 
    rg_Point2D& outputVector2) const
{
     double e = m_CoefficientOfRestitution;

    rg_Point2D c1 = movingDisk1->getCircle().getCenterPt();
    rg_Point2D c2 = movingDisk2->getCircle().getCenterPt();

    rg_Point2D v1 = rg_Point2D(movingDisk1->getVelocityVectorX(), movingDisk1->getVelocityVectorY());
    rg_Point2D v2 = rg_Point2D(movingDisk2->getVelocityVectorX(), movingDisk2->getVelocityVectorY());

    double m1 = pow(movingDisk1->getCircle().getRadius(), 2);
    double m2 = pow(movingDisk2->getCircle().getRadius(), 2);


    if (movingDisk1->getVelocityVector().magnitude() < resNeg3 )
    {
        rg_Point2D basisVector((c2 - c1).getUnitVector());
        rg_Point2D newVector(v2 - 2 * (v2 % basisVector) * basisVector);

        outputVector1.setPoint(0.0, 0.0);
        outputVector2.setPoint(e * newVector);
    }
    else if (movingDisk2->getVelocityVector().magnitude() < resNeg3)
    {
        rg_Point2D basisVector((c1 - c2).getUnitVector());
        rg_Point2D newVector(v1 - 2 * (v1 % basisVector) * basisVector);

        outputVector1.setPoint(e * newVector);
        outputVector2.setPoint(0.0, 0.0);
    }
    else if ((movingDisk1->getVelocityVector().magnitude() >= resNeg3) && (movingDisk2->getVelocityVector().magnitude() >= resNeg3))
    {
        
        const rg_Point2D basisVector = (c2 - c1).getUnitVector();

        double vb1_x = (v1 % basisVector);
        double vb2_x = (v2 % basisVector);

        double vb1_x_afterCollision = vb1_x + (vb2_x - vb1_x)*(e + 1) * m2 / (m1 + m2);
        double vb2_x_afterCollision = vb2_x + (vb1_x - vb2_x)*(e + 1) * m1 / (m1 + m2);

        double cosTheta = basisVector.getX();
        double sinTheta = basisVector.getY();

        double vb1_y_afterCollision = v1.getY()*cosTheta - v1.getX()*sinTheta;
        double vb2_y_afterCollision = v2.getY()*cosTheta - v2.getX()*sinTheta;

        double v1_x_afterCollision = vb1_x_afterCollision*cosTheta - vb1_y_afterCollision*sinTheta;
        double v1_y_afterCollision = vb1_x_afterCollision*sinTheta + vb1_y_afterCollision*cosTheta;

        double v2_x_afterCollision = vb2_x_afterCollision*cosTheta - vb2_y_afterCollision*sinTheta;
        double v2_y_afterCollision = vb2_x_afterCollision*sinTheta + vb2_y_afterCollision*cosTheta;


        outputVector1.setPoint(v1_x_afterCollision, v1_y_afterCollision);
        outputVector2.setPoint(v2_x_afterCollision, v2_y_afterCollision);
    }
}


void DiskCollisionHandler::compute_velocity_vectors_of_disks_after_collided_for_AVOIDANCE(
    const DynamicDisk* movingDisk1, 
    const DynamicDisk* movingDisk2, 
    rg_Point2D& outputVector1, 
    rg_Point2D& outputVector2) const
{
     double e = m_CoefficientOfRestitution;

    rg_Point2D c1 = movingDisk1->getCircle().getCenterPt();
    rg_Point2D c2 = movingDisk2->getCircle().getCenterPt();

    rg_Point2D v1 = rg_Point2D(movingDisk1->getVelocityVectorX(), movingDisk1->getVelocityVectorY());
    rg_Point2D v2 = rg_Point2D(movingDisk2->getVelocityVectorX(), movingDisk2->getVelocityVectorY());

    double m1 = pow(movingDisk1->getCircle().getRadius(), 2);
    double m2 = pow(movingDisk2->getCircle().getRadius(), 2);


    if (movingDisk1->getID() == -1 || movingDisk2->getID() == -1)
    {
        rg_Point2D basisVector((c2 - c1).getUnitVector());
        rg_Point2D newVector(v2 - 2 * (v2 % basisVector) * basisVector);

        outputVector1.setPoint(0.0, 0.0);
        outputVector2.setPoint(0.0, 0.0);
    }
    else if (movingDisk1->getID() == -2)
    {
        rg_Point2D vecFromC1ToC2 = c2 - c1;
        double angleOfC2         = angleFromVec1toVec2(v1, vecFromC1ToC2);

        rg_Point2D v2_adjusted;
        double adjustingAngle;

        if (angleOfC2 >= 0 && angleOfC2 < rg_PI)
        {
            adjustingAngle = rg_PI / 2.0;
        }
        else 
        {
            adjustingAngle = -rg_PI / 2.0;
        }

        v2_adjusted.setX(v1.getX()* cos(adjustingAngle) - v1.getY() * sin(adjustingAngle));
        v2_adjusted.setY(v1.getX()* sin(adjustingAngle) + v1.getY() * cos(adjustingAngle));

        outputVector1 = v1;
        outputVector2 = v2_adjusted * 2.0;

    }
    else if (movingDisk2->getID() == -2)
    {
        rg_Point2D vecFromC2ToC1 = c1 - c2;
        double angleOfC1         = angleFromVec1toVec2(v2, vecFromC2ToC1);

        rg_Point2D v1_adjusted;
        double adjustingAngle;

        if (angleOfC1 >= 0 && angleOfC1 < rg_PI)
        {
            adjustingAngle = rg_PI / 2.0;
        }
        else
        {
            adjustingAngle = -rg_PI / 2.0;
        }
        v1_adjusted.setX(v2.getX()* cos(adjustingAngle) - v2.getY() * sin(adjustingAngle));
        v1_adjusted.setY(v2.getX()* sin(adjustingAngle) + v2.getY() * cos(adjustingAngle));

        outputVector1 = v1_adjusted * 2.0;
        outputVector2 = v2;
    }
    else
    {

        double m1 = pow(movingDisk1->getCircle().getRadius(), 2);
        double m2 = pow(movingDisk2->getCircle().getRadius(), 2);

        const rg_Point2D basisVector = (c2 - c1).getUnitVector();

        double vb1_x = (v1 % basisVector);
        double vb2_x = (v2 % basisVector);

        double vb1_x_afterCollision = vb1_x + (vb2_x - vb1_x)*(e + 1) * m2 / (m1 + m2);
        double vb2_x_afterCollision = vb2_x + (vb1_x - vb2_x)*(e + 1) * m1 / (m1 + m2);

        double cosTheta = basisVector.getX();
        double sinTheta = basisVector.getY();

        double vb1_y_afterCollision = v1.getY()*cosTheta - v1.getX()*sinTheta;
        double vb2_y_afterCollision = v2.getY()*cosTheta - v2.getX()*sinTheta;

        double v1_x_afterCollision = vb1_x_afterCollision * cosTheta - vb1_y_afterCollision * sinTheta;
        double v1_y_afterCollision = vb1_x_afterCollision * sinTheta + vb1_y_afterCollision * cosTheta;

        double v2_x_afterCollision = vb2_x_afterCollision * cosTheta - vb2_y_afterCollision * sinTheta;
        double v2_y_afterCollision = vb2_x_afterCollision * sinTheta + vb2_y_afterCollision * cosTheta;

        outputVector1.setPoint(v1_x_afterCollision, v1_y_afterCollision);
        outputVector2.setPoint(v2_x_afterCollision, v2_y_afterCollision);  
    }
}


void DiskCollisionHandler::compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION_N_RANDOM_SPEED(
    const DynamicDisk* movingDisk1, 
    const DynamicDisk* movingDisk2, 
    rg_Point2D& outputVector1, 
    rg_Point2D& outputVector2) const
{
    double e = m_CoefficientOfRestitution;

    rg_Point2D c1 = movingDisk1->getCircle().getCenterPt();
    rg_Point2D c2 = movingDisk2->getCircle().getCenterPt();

    rg_Point2D v1 = rg_Point2D(movingDisk1->getVelocityVectorX(), movingDisk1->getVelocityVectorY());
    rg_Point2D v2 = rg_Point2D(movingDisk2->getVelocityVectorX(), movingDisk2->getVelocityVectorY());

    double m1     = pow(movingDisk1->getCircle().getRadius(), 2);
    double m2     = pow(movingDisk2->getCircle().getRadius(), 2);

    if (movingDisk2->getCircle().isIncludedIn(movingDisk1->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        rg_Point2D basisVector((c2 - c1).getUnitVector());
        rg_Point2D newVector(v2 - 2 * (v2 % basisVector) * basisVector);

        outputVector1.setPoint(0.0, 0.0);
        outputVector2.setPoint(e * newVector);
    }
    else if (movingDisk1->getCircle().isIncludedIn(movingDisk2->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        rg_Point2D basisVector((c1 - c2).getUnitVector());
        rg_Point2D newVector(v1 - 2 * (v1 % basisVector) * basisVector);

        outputVector1.setPoint(e * newVector);
        outputVector2.setPoint(0.0, 0.0);
    }
    else
    {
        /*

            |A|*cost*B / |B| = (A dot B) / |B| * B


            solve equations

            known: m1,m2,v1,v2, unknown: u1, u2

            1. m1v1 + m2v2 = m1u1 + m2u2
            2. e = (u1 - u2) / (v2 - v1)

            u1 = v1 + (v2 - v1)*(e + 1) * m2 / (m1 + m2)
            u2 = v2 + (v1 - v2)*(e + 1) * m1 / (m1 + m2)

            */

        double m1 = pow(movingDisk1->getCircle().getRadius(), 2);
        double m2 = pow(movingDisk2->getCircle().getRadius(), 2);

        const rg_Point2D basisVector = (c2 - c1).getUnitVector();

        double vb1_x = (v1 % basisVector);
        double vb2_x = (v2 % basisVector);

        double vb1_x_afterCollision = vb1_x + (vb2_x - vb1_x)*(e + 1) * m2 / (m1 + m2);
        double vb2_x_afterCollision = vb2_x + (vb1_x - vb2_x)*(e + 1) * m1 / (m1 + m2);

        double cosTheta = basisVector.getX();
        double sinTheta = basisVector.getY();

        double vb1_y_afterCollision = v1.getY()*cosTheta - v1.getX()*sinTheta;
        double vb2_y_afterCollision = v2.getY()*cosTheta - v2.getX()*sinTheta;

        double v1_x_afterCollision = vb1_x_afterCollision * cosTheta - vb1_y_afterCollision * sinTheta;
        double v1_y_afterCollision = vb1_x_afterCollision * sinTheta + vb1_y_afterCollision * cosTheta;

        double v2_x_afterCollision = vb2_x_afterCollision * cosTheta - vb2_y_afterCollision * sinTheta;
        double v2_y_afterCollision = vb2_x_afterCollision * sinTheta + vb2_y_afterCollision * cosTheta;

        if (movingDisk1->this_disk_is_container())
        {
            outputVector1.setPoint(0.0, 0.0);
            outputVector2.setPoint(v2_x_afterCollision, v2_y_afterCollision);
            outputVector2 = outputVector2.getUnitVector() * get_randomly_generated_zero_to_one();
        }
        else if (movingDisk2->this_disk_is_container())
        {
            outputVector1.setPoint(v1_x_afterCollision, v1_y_afterCollision);
            outputVector1 = outputVector1.getUnitVector() * get_randomly_generated_zero_to_one();
            outputVector2.setPoint(0.0, 0.0);
        }
        else
        {
            outputVector1.setPoint(v1_x_afterCollision, v1_y_afterCollision);
            outputVector2.setPoint(v2_x_afterCollision, v2_y_afterCollision);

            find_optimal_velocity_vectors_for_PHYSICAL_COLLISION_N_RANDOM_SPEED(movingDisk1, movingDisk2, outputVector1, outputVector2);
        }
    }
}


void DiskCollisionHandler::find_optimal_velocity_vectors_for_PHYSICAL_COLLISION_N_RANDOM_SPEED(
    const DynamicDisk* movingDisk1, 
    const DynamicDisk* movingDisk2, 
    rg_Point2D& outputVector1, 
    rg_Point2D& outputVector2) const
{
    bool optimalVelocityVectorsAreFound = false;

    while (!optimalVelocityVectorsAreFound)
    {
        rg_Point2D candidateVector1 = outputVector1.getUnitVector() * get_randomly_generated_zero_to_one();
		rg_Point2D candidateVector2 = outputVector2.getUnitVector() * get_randomly_generated_zero_to_one();

		rg_Point2D c1_c2 = movingDisk1->getCircle().getCenterPt() - movingDisk2->getCircle().getCenterPt();
		rg_Point2D v1_v2 = candidateVector1 - candidateVector2;

		double r1 = movingDisk1->getCircle().getRadius() + movingDisk1->getOffset();
		double r2 = movingDisk2->getCircle().getRadius() + movingDisk2->getOffset();

		double A = 0.0;
		double B = 0.0;
		double collisionTime = 0.0;

		if (v1_v2.magnitude() >= rg_MATH_RES) //2017.10.03 cysong
		{
			if (pow(v1_v2 % c1_c2, 2) - (v1_v2 % v1_v2) * ((c1_c2 % c1_c2) - pow(r1 + r2, 2)) > TOLERANCE_OF_COLLISION_PREDICTION) //will_two_particles_be_collided()
			{
				A = (v1_v2 % c1_c2);
				B = sqrt(pow(v1_v2 % c1_c2, 2) - v1_v2.magnitudeSquare() * (c1_c2.magnitudeSquare() - pow(r1 + r2, 2)));

				if (rg_EQ(-A + B, 0.0, resNeg3))
				{
                    outputVector1 = candidateVector1;
                    outputVector2 = candidateVector2;
                    optimalVelocityVectorsAreFound = true;
				}
			}
		}
    }
}

