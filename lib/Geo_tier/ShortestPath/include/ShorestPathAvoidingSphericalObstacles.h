#ifndef _SHORTEST_PATH_AVOIDING_SPHERICAL_OBSTACLES_H
#define _SHORTEST_PATH_AVOIDING_SPHERICAL_OBSTACLES_H


#include "rg_Point3D.h"
#include "Sphere.h"
#include "AxisAlignedBox.h"
#include "rg_NURBSplineCurve3D.h"
#include "Plane.h"

#include <list>
using namespace std;


const   string  FF_STARTINGPOINT   = "STARTING_POINT";
const   string  FF_DESTINATION     = "DESTINATION";
const   string  FF_VEHICLE         = "VEHICLE";
const   string  FF_OBSTACLE        = "OBSTACLE";

const double  PHI = 3.141592653589793238462643383279; // ø¯¡÷¿≤

class ShorestPathAvoidingSphericalObstacles
{
private:

private:
    rg_Point3D              m_starting_point;
    rg_Point3D              m_destination;
    Sphere                  m_vehicle;

    list<Sphere>            m_obstacles;

    rg_NURBSplineCurve3D    m_shortest_path;


    list<Sphere>            m_obstacles_intersecting_straight_line_path;
    list<Circle3D>          m_obstacles_in_plane;
    vector<rg_Point3D>      m_points_on_the_shortest_path;

public:
    ShorestPathAvoidingSphericalObstacles();
    ShorestPathAvoidingSphericalObstacles(const ShorestPathAvoidingSphericalObstacles& path);
    ~ShorestPathAvoidingSphericalObstacles();

    void    clear();

    ShorestPathAvoidingSphericalObstacles& operator =(const ShorestPathAvoidingSphericalObstacles& path);

    rg_Point3D              starting_point() const;
    rg_Point3D              destination() const;
    Sphere                  vehicle() const;
    const list<Sphere>&     obstacles() const;
    const list<Sphere>&     obstacles_intersecting_straight_line_path() const;
    const list<Circle3D>&   obstacles_in_plane() const;

    const vector<rg_Point3D>&    points_on_the_shortest_path() const;

    void                set_starting_point(const rg_Point3D& startingPoint);
    void                set_destination(const rg_Point3D& destination);
    void                set_vehicle(const Sphere& vehicle);
    void                set_obstacles(const list<Sphere>& obstacles);
    void                set_problem(const rg_Point3D& startingPoint, const rg_Point3D& destination, const Sphere& vehicle, const list<Sphere>& obstacles);

    double              approximated_path_length() const;
    void                points_on_the_shortest_path(const int& numPts, vector<rg_Point3D>& pointsOnPath) const;

    AxisAlignedBox      compute_bounding_space() const;

    bool                find_the_shortest_path();


    void                read(const string& filename);
    void                write(const string& filename) const;

    void                write_way_points_in_local_coordinates(const string& filename) const;
    void                write_way_points(const string& filename,
                            const double& fence_length, const double& fence_width, const double& fence_height, 
                            const double& home_latitude, const double& home_longitude, const double& home_altitude) const;
    void                write_fence(const string& filename,
                            const double& fence_length, const double& fence_width, const double& fence_height,
                            const double& home_latitude, const double& home_longitude, const double& home_altitude) const;

private:
    unsigned int        find_obstacles_intersecting_straight_line_path(list<Sphere>& intersectingObstacles);
    Plane               define_good_plane_through_straight_line_path(const list<Sphere>& intersectingObstacles);
    void                find_intersections_between_plane_and_obstacles(const Plane& plane4Path, list<Circle3D>& intersections);
    void                make_matrix_to_transform_plane_into_XYplane(const Plane& plane4Path, rg_TMatrix3D& transformMatrix, rg_TMatrix3D& inverseMatrix);

    void                calculate_weights_of_obstacles_intersecting_line_path_for_determining_plane(const vector<double>& distances, vector<double>& weights, const double& lambda = 5.0);
    //void                find_shortest_path_on_plane(const Plane& plane4Path);

    double              degreeToRadian(double a) const { return (PHI * a / 180.0); }


};


inline  rg_Point3D          ShorestPathAvoidingSphericalObstacles::starting_point() const { return m_starting_point; }
inline  rg_Point3D          ShorestPathAvoidingSphericalObstacles::destination() const { return m_destination; }
inline  Sphere              ShorestPathAvoidingSphericalObstacles::vehicle() const { return m_vehicle; }
inline  const list<Sphere>& ShorestPathAvoidingSphericalObstacles::obstacles() const { return m_obstacles; }
inline  const list<Sphere>& ShorestPathAvoidingSphericalObstacles::obstacles_intersecting_straight_line_path() const { return m_obstacles_intersecting_straight_line_path; }
inline  const list<Circle3D>&     ShorestPathAvoidingSphericalObstacles::obstacles_in_plane() const { return m_obstacles_in_plane; }
inline  const vector<rg_Point3D>& ShorestPathAvoidingSphericalObstacles::points_on_the_shortest_path() const { return m_points_on_the_shortest_path; }


inline  void                ShorestPathAvoidingSphericalObstacles::set_starting_point(const rg_Point3D& startingPoint) { m_starting_point = startingPoint; }
inline  void                ShorestPathAvoidingSphericalObstacles::set_destination(const rg_Point3D& destination) { m_destination = destination; }
inline  void                ShorestPathAvoidingSphericalObstacles::set_vehicle(const Sphere& vehicle) { m_vehicle = vehicle; }
inline  void                ShorestPathAvoidingSphericalObstacles::set_obstacles(const list<Sphere>& obstacles) { m_obstacles = obstacles; }

inline  double              ShorestPathAvoidingSphericalObstacles::approximated_path_length() const { return m_shortest_path.evaluateCurveLengthByCtrlPts(); }

#endif

