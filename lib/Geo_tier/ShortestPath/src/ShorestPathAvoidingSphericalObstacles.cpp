
#include "ShorestPathAvoidingSphericalObstacles.h"
#include "LineSegment3D.h"
#include "ShortestPathFinder2D.h"

#include <fstream>
#include <string>
#include <cmath>
using namespace std;



ShorestPathAvoidingSphericalObstacles::ShorestPathAvoidingSphericalObstacles()
{
}



ShorestPathAvoidingSphericalObstacles::ShorestPathAvoidingSphericalObstacles(const ShorestPathAvoidingSphericalObstacles& path)
: m_starting_point(path.m_starting_point)
, m_destination(path.m_destination)
, m_vehicle(path.m_vehicle)
, m_obstacles(path.m_obstacles)
, m_shortest_path(path.m_shortest_path)
{
}



ShorestPathAvoidingSphericalObstacles::~ShorestPathAvoidingSphericalObstacles()
{
}



void    ShorestPathAvoidingSphericalObstacles::clear()
{
    m_obstacles.clear();
    m_shortest_path.removeAll();

    m_obstacles_intersecting_straight_line_path.clear();
    m_obstacles_in_plane.clear();
    m_points_on_the_shortest_path.clear();
}



ShorestPathAvoidingSphericalObstacles& ShorestPathAvoidingSphericalObstacles::operator =(const ShorestPathAvoidingSphericalObstacles& path)
{
    if (this != &path) {
        m_starting_point = path.m_starting_point;
        m_destination    = path.m_destination;
        m_vehicle        = path.m_vehicle;
        m_obstacles      = path.m_obstacles;

        m_shortest_path  = path.m_shortest_path;
    }

    return *this;
}



void ShorestPathAvoidingSphericalObstacles::set_problem(const rg_Point3D& startingPoint, const rg_Point3D& destination, const Sphere& vehicle, const list<Sphere>& obstacles) 
{
    clear();

    m_starting_point = startingPoint;
    m_destination    = destination;
    m_vehicle        = vehicle;
    
    for (list<Sphere>::const_iterator i_obstacle = obstacles.begin(); i_obstacle != obstacles.end(); ++i_obstacle) {
        Sphere currObstacle = (*i_obstacle);
        m_obstacles.push_back( Sphere(currObstacle.getCenter(), currObstacle.getRadius() + vehicle.getRadius()) );
    }
}



void                ShorestPathAvoidingSphericalObstacles::points_on_the_shortest_path(const int& numPts, vector<rg_Point3D>& pointsOnPath) const
{
    rg_Point3D* ptsOnCurve = m_shortest_path.evaluatePtsInEvenParameter(numPts);

    pointsOnPath.resize(numPts);
    for (int i = 0; i < numPts; ++i) {
        pointsOnPath[i] = ptsOnCurve[i];
    }
    delete[] ptsOnCurve;
}



AxisAlignedBox      ShorestPathAvoidingSphericalObstacles::compute_bounding_space() const
{
    AxisAlignedBox boundingBox;

    for (list<Sphere>::const_iterator i_obstacle = m_obstacles.begin(); i_obstacle != m_obstacles.end(); ++i_obstacle) {
        boundingBox.update( *i_obstacle );
    }

    boundingBox.update(m_vehicle);
    boundingBox.update(m_starting_point);
    boundingBox.update(m_destination);


    return boundingBox;
}



bool ShorestPathAvoidingSphericalObstacles::find_the_shortest_path()
{
    bool    isThereAPath = true;
    m_obstacles_intersecting_straight_line_path.clear();
    m_obstacles_in_plane.clear();
    m_points_on_the_shortest_path.clear();



    unsigned int numObstaclesIntersectingStraightLinePath = find_obstacles_intersecting_straight_line_path(m_obstacles_intersecting_straight_line_path);

    if (numObstaclesIntersectingStraightLinePath == 0) {
        m_shortest_path.formLine(m_starting_point, m_destination);

        m_points_on_the_shortest_path.push_back(m_starting_point);
        m_points_on_the_shortest_path.push_back(m_destination);

        isThereAPath = true;
    }
    else {
        Plane plane4Path = define_good_plane_through_straight_line_path(m_obstacles_intersecting_straight_line_path);

        //  we will use a cylindrical filter.
        //list<Circle3D> intersections;
        find_intersections_between_plane_and_obstacles(plane4Path, m_obstacles_in_plane);



        //  compute this problem in 2D
        rg_TMatrix3D transformMatrix;
        rg_TMatrix3D inverseMatrix;
        make_matrix_to_transform_plane_into_XYplane(plane4Path, transformMatrix, inverseMatrix);


        rg_Circle2D startPtIn2D( (transformMatrix*m_starting_point).evaluatePt2D(), 0.0);
        rg_Circle2D endPtIn2D(   (transformMatrix*m_destination).evaluatePt2D(), 0.0);
        rg_Circle2D vehicleIn2D( (transformMatrix*m_vehicle.getCenter()).evaluatePt2D(), 0.0);
        list<rg_Circle2D> obstaclesIn2D;
        for (list<Circle3D>::iterator i_intersect = m_obstacles_in_plane.begin(); i_intersect != m_obstacles_in_plane.end(); ++i_intersect) {
            obstaclesIn2D.push_back(rg_Circle2D((transformMatrix*i_intersect->getCenter()).evaluatePt2D(), i_intersect->getRadius()));
        }


        ShortestPathFinder2D solverIn2D;
        isThereAPath = solverIn2D.compute_geodesic_path_based_on_VD(startPtIn2D, endPtIn2D, vehicleIn2D, obstaclesIn2D);

        if (isThereAPath) {

            list<EdgeForSolutionPoolGraph*> pathIn2D;
            solverIn2D.get_geodesic_path_by_BFS(pathIn2D);


            VertexForSolutionPoolGraph* prevVtx = pathIn2D.front()->get_start_vtx();
            rg_Point2D startingPtIn2D = prevVtx->get_tangent_point_to_disk();
            if (!(startingPtIn2D == startPtIn2D.getCenterPt())) {
                prevVtx = pathIn2D.front()->get_end_vtx();
                startingPtIn2D = prevVtx->get_tangent_point_to_disk();
            }
            rg_Point3D startingPtIn3D(startingPtIn2D.getX(), startingPtIn2D.getY(), 0.0);



            m_points_on_the_shortest_path.push_back(inverseMatrix*startingPtIn3D);
            for (list<EdgeForSolutionPoolGraph*>::iterator i_edge = pathIn2D.begin(); i_edge != pathIn2D.end(); ++i_edge) {
                EdgeForSolutionPoolGraph*   currEdge = *i_edge;
                VertexForSolutionPoolGraph* nextVtx = currEdge->get_opposite_vtx(prevVtx);

                rg_Point2D endPtIn2D = nextVtx->get_tangent_point_to_disk();
                rg_Point3D endPtIn3D(endPtIn2D.getX(), endPtIn2D.getY(), 0.0);

                m_points_on_the_shortest_path.push_back(inverseMatrix*endPtIn3D);

                prevVtx = nextVtx;
            }




            int numPts = m_points_on_the_shortest_path.size();
            rg_Point3D* pts = new rg_Point3D[numPts];
            for (int i = 0; i < numPts; ++i) {
                pts[i] = m_points_on_the_shortest_path[i];
            }

            m_shortest_path.curveInterpolation(numPts, pts);
            delete[] pts;
        }
    }

    if (isThereAPath) {

        double      curveLength = m_shortest_path.evaluateCurveLengthByCtrlPts();
        int numPtsOnCurve = (curveLength > 100.0) ? (int)curveLength : 100;
        rg_Point3D* ptsOnCurve = m_shortest_path.evaluatePtsInEvenParameter(numPtsOnCurve);

        m_points_on_the_shortest_path.clear();
        m_points_on_the_shortest_path.resize(numPtsOnCurve);
        for (int i = 0; i < numPtsOnCurve; ++i) {
            m_points_on_the_shortest_path[i] = ptsOnCurve[i];
        }
        delete[] ptsOnCurve;

    }

    return isThereAPath;

}



unsigned int     ShorestPathAvoidingSphericalObstacles::find_obstacles_intersecting_straight_line_path(list<Sphere>& intersectingObstacles)
{
    LineSegment3D straight_path(m_starting_point, m_destination);

    for (list<Sphere>::iterator i_obstacle = m_obstacles.begin(); i_obstacle != m_obstacles.end(); ++i_obstacle) {
        Sphere currObstacle = (*i_obstacle);

        if (currObstacle.hasIntersectionWith(straight_path)) {
            intersectingObstacles.push_back(currObstacle);
        }
    }

    return intersectingObstacles.size();
}



Plane   ShorestPathAvoidingSphericalObstacles::define_good_plane_through_straight_line_path(const list<Sphere>& intersectingObstacles)
{
    rg_Point3D  normal;


    rg_Point3D  vectorFromS2T = m_destination - m_starting_point;
    rg_Point3D  vectorOfZ(1.0, 0.0, 0.0);

    normal = vectorFromS2T.crossProduct(vectorOfZ);
    normal.normalize();

    if (normal.innerProduct(rg_Point3D(0.0, 0.0, 1.0)) < 0.0) {
        normal = -1.0*normal;
    }

    /*
    LineSegment3D straight_path(m_starting_point, m_destination);

    unsigned int  numObstacles = intersectingObstacles.size();
    vector<rg_Point3D> vectorFromFootprintToCenter(numObstacles);
    vector<double>     distanceToFootprint(numObstacles);
    vector<double>     weightOfObstacle;
    int i = 0;
    for (list<Sphere>::const_iterator i_obstacle = intersectingObstacles.begin(); i_obstacle != intersectingObstacles.end(); ++i_obstacle, ++i) {
        rg_REAL param = 0.0;
        rg_Point3D footprint = straight_path.computeProjectedPointForPoint3D(i_obstacle->getCenter(), param);
        vectorFromFootprintToCenter[i] = i_obstacle->getCenter() - footprint;
        distanceToFootprint[i] = m_starting_point.distance(footprint);
    }

    calculate_weights_of_obstacles_intersecting_line_path_for_determining_plane(distanceToFootprint, weightOfObstacle, 5.);

    for (i = 0; i < numObstacles; ++i) {
        normal += (weightOfObstacle[i] * vectorFromFootprintToCenter[i]);
    }


    if (normal == rg_Point3D()) {
        //  make a plane perpendicular to the xy-plane. 
        rg_Point3D  vectorFromS2T = m_destination - m_starting_point;
        rg_Point3D  vectorOfZ(0.0, 0.0, 1.0);

        normal = vectorFromS2T.crossProduct(vectorOfZ);
    }
    */

    return Plane(normal, m_starting_point);
    
}



void    ShorestPathAvoidingSphericalObstacles::calculate_weights_of_obstacles_intersecting_line_path_for_determining_plane(
            const vector<double>& distances, vector<double>& weights, const double& lambda)
{
    double sumOfWeights = 0.0;

    for (vector<double>::const_iterator itForDistance = distances.begin(); itForDistance != distances.end(); ++itForDistance)
    {
        double weight = lambda * exp(-lambda * (*itForDistance));
        weights.push_back(weight);
        sumOfWeights += weight;
    }

    for (vector<double>::iterator itForWeight = weights.begin(); itForWeight != weights.end(); ++itForWeight)
    {
        (*itForWeight) = (*itForWeight) / sumOfWeights;
    }

}


void    ShorestPathAvoidingSphericalObstacles::find_intersections_between_plane_and_obstacles(const Plane& plane4Path, list<Circle3D>& intersections)
{
    for (list<Sphere>::iterator i_obstacle = m_obstacles.begin(); i_obstacle != m_obstacles.end(); ++i_obstacle) {
        Sphere currObstacle = (*i_obstacle);

        Circle3D intersectingCircle;
        if (currObstacle.intersect(plane4Path, intersectingCircle) == rg_TRUE) {
            intersections.push_back(intersectingCircle);
        }
    }
}



void    ShorestPathAvoidingSphericalObstacles::make_matrix_to_transform_plane_into_XYplane(const Plane& plane4Path, rg_TMatrix3D& transformMatrix, rg_TMatrix3D& inverseMatrix)
{
    rg_TMatrix3D forward[3], backward[3];
    forward[0].translate(-m_starting_point);
    backward[0].translate(m_starting_point);

    rg_Point3D   axis = (m_destination - m_starting_point).getUnitVector();

    rg_FLAG bReverse = rg_EQ(plane4Path.getNormal().innerProduct(rg_Point3D(0.0, 0.0, 1.0)), -1);
    if (!bReverse) {
        forward[1].rotate(plane4Path.getNormal(), rg_Point3D(0.0, 0.0, 1.0));
    }
    else {
        forward[1].rotateY(rg_PI);
    }
    forward[2].rotate(forward[1] * axis, rg_Point3D(1.0, 0.0, 0.0));

    if (!bReverse) {
        backward[1].rotate(rg_Point3D(0.0, 0.0, 1.0), plane4Path.getNormal());
    }
    else {
        backward[1].rotateY(-rg_PI);
    }
    backward[2].rotate(rg_Point3D(1.0, 0.0, 0.0), forward[1] * axis);

    transformMatrix = forward[2] * forward[1] * forward[0];
    inverseMatrix = backward[0] * backward[1] * backward[2];
}



void    ShorestPathAvoidingSphericalObstacles::read(const string& filename)
{
    clear();


    const int   BUFFER_SIZE = 300;
    char	    lineData[BUFFER_SIZE];
    char*	    seps = " :\t\n";
    char*	    token = NULL;
    char*       context = NULL;


    ifstream fin;
    fin.open(filename.c_str());

    while (!fin.eof()) {
        fin.getline(lineData, BUFFER_SIZE);

        token = strtok(lineData, seps);
        if (token == rg_NULL) {
            continue;
        }

        if (FF_STARTINGPOINT.compare(token) == 0) {
            double  x = atof(strtok(NULL, seps));
            double  y = atof(strtok(NULL, seps));
            double  z = atof(strtok(NULL, seps));
            m_starting_point.setPoint(x, y, z);
        }        
        else if (FF_DESTINATION.compare(token) == 0) {
            double  x = atof(strtok(NULL, seps));
            double  y = atof(strtok(NULL, seps));
            double  z = atof(strtok(NULL, seps));
            m_destination.setPoint(x, y, z);
        }
        else if (FF_VEHICLE.compare(token) == 0) {
            double  x = atof(strtok(NULL, seps));
            double  y = atof(strtok(NULL, seps));
            double  z = atof(strtok(NULL, seps));
            double  r = atof(strtok(NULL, seps));
            m_vehicle.setSphere(x, y, z, r);
        }
        else if (FF_OBSTACLE.compare(token) == 0) {
            double  x = atof(strtok(NULL, seps));
            double  y = atof(strtok(NULL, seps));
            double  z = atof(strtok(NULL, seps));
            double  r = atof(strtok(NULL, seps));
            m_obstacles.push_back( Sphere(x, y, z, r + m_vehicle.getRadius() ) );
        }
    }

    fin.close();
}



void    ShorestPathAvoidingSphericalObstacles::write(const string& filename) const
{
    //const   string  FF_STARTINGPOINT = "STARTING_POINT";
    //const   string  FF_DESTINATION = "DESTINATION";
    //const   string  FF_VEHICLE = "VEHICLE";
    //const   string  FF_OBSTACLE = "OBSTACLE";

    ofstream fout;
    fout.open(filename.c_str());

    fout << FF_STARTINGPOINT.c_str() << "\t" << m_starting_point.getX() << "\t" << m_starting_point.getY() << "\t" << m_starting_point.getZ() << endl;
    fout << FF_DESTINATION.c_str()   << "\t" << m_destination.getX() << "\t" << m_destination.getY() << "\t" << m_destination.getZ() << endl;
    fout << FF_VEHICLE.c_str()       << "\t" << m_vehicle.getCenter().getX() << "\t" << m_vehicle.getCenter().getY() << "\t" << m_vehicle.getCenter().getZ() << "\t" << m_vehicle.getRadius() << endl;

    double vehicleRadius = m_vehicle.getRadius();
    for (list<Sphere>::const_iterator i_obstacle = m_obstacles.begin(); i_obstacle != m_obstacles.end(); ++i_obstacle) {
        fout << FF_OBSTACLE.c_str() << "\t" 
             << i_obstacle->getCenter().getX() << "\t" << i_obstacle->getCenter().getY() << "\t" << i_obstacle->getCenter().getZ() 
             << "\t" << (i_obstacle->getRadius() - vehicleRadius) << endl;
    }

    fout.close();
}


void    ShorestPathAvoidingSphericalObstacles::write_way_points_in_local_coordinates(const string& filename) const
{
    ofstream fout;
    fout.open(filename.c_str());

    int numPoints = m_points_on_the_shortest_path.size();
    for (int i = 0; i < numPoints; ++i) {
        rg_Point3D pt = m_points_on_the_shortest_path[i];

        fout << i + 1 << "\t";
        fout << pt.getX() << "\t" << pt.getY() << "\t" << pt.getZ() << endl;
    }

    fout.close();
}


void    ShorestPathAvoidingSphericalObstacles::write_way_points(const string& filename,
            const double& fence_length, const double& fence_width, const double& fence_height,
            const double& home_latitude, const double& home_longitude, const double& home_altitude) const
{
    double  RADIUS_EARTH_EQUATOR = 6378245.0;    // 지구 적도 반경 (단위:m)
    double  RADIUS_EARTH_POLAR = 6356863.0;    // 지구 극 반경 (단위:m)
    double  RATIO_MILE_TO_METER = 1609.344;    // 1 mile을 meter로 변환
    double  NORTH_TILT_MAGNETIC_TO_MAP = 6.5;
    int     RADAR_ACP_COUNT = 4096;     // ACP 분해각 개수
                                                 // 전역 변수 정의
    double  Base_Dist_Latitude;  // 현 위도에서의 위도 1도에 대한 거리 (단위:m)
    double  Base_Dist_Longuitude;  // 현 위도에서의 경도 1도에 대한 거리 (단위:m)
    double  Radius_Earth;   // 현 위도에서의 지구 반경

    //double  RADAR_LATITUDE = 35.181692;    // RADAR 위치 : 위도
    //double  RADAR_LONGUITUDE = 128.950575;    // RADAR 위치 : 경도


    Radius_Earth         = RADIUS_EARTH_POLAR + (RADIUS_EARTH_EQUATOR - RADIUS_EARTH_POLAR) * sin( degreeToRadian(home_latitude));
    Base_Dist_Longuitude = 2.0 * PHI * Radius_Earth / 360.0 * cos(degreeToRadian(home_latitude));
    Base_Dist_Latitude   = 2.0 * PHI * Radius_Earth / 360.0;




    ofstream fout;
    fout.open(filename.c_str());

    //fence
    //37.555447	127.047615
    //37.555763	127.047142
    //37.555305	127.047028
    //37.555130	127.048119
    //37.555565	127.048195
    //37.555763	127.047142


    // way points
    fout.precision(7);
    fout << "QGC WPL 110" << endl;
    fout << "0	0	0	16	0.000000	0.000000	0.000000	0.000000" << "\t" 
         << fixed << home_latitude << "\t" << fixed << home_longitude << "\t" << "0.000000	1" << endl;
    
    rg_Point3D  origin = m_points_on_the_shortest_path[0];
    int numPoints = m_points_on_the_shortest_path.size();
    for (int i = 0; i < numPoints; ++i) {
        rg_Point3D pt = m_points_on_the_shortest_path[i] - origin;

        double pt_latitude  = home_latitude  + (pt.getY() / Base_Dist_Latitude);
        double pt_longitude = home_longitude + (pt.getX() / Base_Dist_Longuitude);
        double pt_altitude  = home_altitude  + pt.getZ();
        fout << i+1 << "\t";
        fout << "0	3	16	0.000000	0.000000	0.000000	0.000000" << "\t";
        fout << fixed << pt_latitude << "\t";
        fout << fixed << pt_longitude << "\t";
        fout << fixed << pt_altitude << "\t";
        fout << "1" << endl;
    }

    fout.close();
}



void    ShorestPathAvoidingSphericalObstacles::write_fence(const string& filename,
        const double& fence_length, const double& fence_width, const double& fence_height,
        const double& home_latitude, const double& home_longitude, const double& home_altitude) const
{
    double  RADIUS_EARTH_EQUATOR = 6378245.0;    // 지구 적도 반경 (단위:m)
    double  RADIUS_EARTH_POLAR = 6356863.0;    // 지구 극 반경 (단위:m)
    double  RATIO_MILE_TO_METER = 1609.344;    // 1 mile을 meter로 변환
    double  NORTH_TILT_MAGNETIC_TO_MAP = 6.5;
    int     RADAR_ACP_COUNT = 4096;     // ACP 분해각 개수
                                        // 전역 변수 정의
    double  Base_Dist_Latitude;  // 현 위도에서의 위도 1도에 대한 거리 (단위:m)
    double  Base_Dist_Longuitude;  // 현 위도에서의 경도 1도에 대한 거리 (단위:m)
    double  Radius_Earth;   // 현 위도에서의 지구 반경

                            //double  RADAR_LATITUDE = 35.181692;    // RADAR 위치 : 위도
                            //double  RADAR_LONGUITUDE = 128.950575;    // RADAR 위치 : 경도


    Radius_Earth = RADIUS_EARTH_POLAR + (RADIUS_EARTH_EQUATOR - RADIUS_EARTH_POLAR) * sin(degreeToRadian(home_latitude));
    Base_Dist_Longuitude = 2.0 * PHI * Radius_Earth / 360.0 * cos(degreeToRadian(home_latitude));
    Base_Dist_Latitude = 2.0 * PHI * Radius_Earth / 360.0;



    rg_Point2D  pt[5] = {   rg_Point2D(fence_length / 2, fence_width / 2),
                            rg_Point2D(-fence_length / 2, fence_width / 2),
                            rg_Point2D(-fence_length / 2, -fence_width / 2),
                            rg_Point2D(fence_length / 2, -fence_width / 2), 
                            rg_Point2D(fence_length / 2, fence_width / 2)
                        };

    ofstream fout;
    fout.open(filename.c_str());

    //fence
    //37.555447	127.047615
    //37.555763	127.047142
    //37.555305	127.047028
    //37.555130	127.048119
    //37.555565	127.048195
    //37.555763	127.047142

    fout.precision(7);
    fout << fixed << home_latitude << "\t" << fixed << home_longitude << endl;

    for (int i = 0; i < 5; ++i) {

        double pt_latitude  = home_latitude + (pt[i].getY() / Base_Dist_Latitude);
        double pt_longitude = home_longitude + (pt[i].getX() / Base_Dist_Longuitude);
        fout << fixed << pt_latitude << "\t";
        fout << fixed << pt_longitude << endl;
    }

    fout.close();
}