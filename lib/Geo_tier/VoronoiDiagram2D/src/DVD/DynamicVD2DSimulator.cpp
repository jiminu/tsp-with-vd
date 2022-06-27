/*
* NOTE. EventHistory Read/Write
 - 1) Write: based on t = 0
 - 2) Read: assume that EVENTSEQ data is generated based on t = 0
*/

#include "DynamicVD2DSimulator.h"
using namespace V::GeometryTier;

#include "IOFunctionsForDisks.h"
#include "StringFunctions.h"
#include <string.h>
#include <fstream>
#include <sstream>

using namespace std;

#if defined(_WIN32_WINDOWS) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
#include <windows.h> 
#define localtime_r(_time, _result) localtime_s(_result, _time)
#endif


DynamicVD2DSimulator::DynamicVD2DSimulator()
{
    m_DynamicDisksAreGeneratedInsideOfThisObject  = false;

    m_CompTimeStatistics.reset(DVDS_TIME_SIZE);
}


DynamicVD2DSimulator::DynamicVD2DSimulator(const DynamicVD2DSimulator& DVDSimulator)
{
    copy_from(DVDSimulator);
}


DynamicVD2DSimulator::~DynamicVD2DSimulator()
{
    clear();
}


DynamicVD2DSimulator& DynamicVD2DSimulator::operator=(const DynamicVD2DSimulator& dvdSimulator)
{
    if (this == &dvdSimulator)
    {
        return *this;
    }

    clear();
    copy_from(dvdSimulator);

    return *this;
}




void DynamicVD2DSimulator::clear()
{
    if (m_DynamicDisksAreGeneratedInsideOfThisObject)
    {
        delete_dynamic_disks();
        delete_containers();
    }
    else
    {
        m_DynamicDisks.clear();
        m_AllContainers.clear();
    }

    m_DynamicVD.clear();
    m_CompTimeStatistics.reset(0);
}


void DynamicVD2DSimulator::clear_EVENTSEQ_related_objects()
{
    m_DynamicVD.clear_EVENTSEQ_related_objects();
    m_DynamicVD.prepare_statistics_in_generating_EVENTSEQ();
}


void DynamicVD2DSimulator::copy_from(const DynamicVD2DSimulator& dvdSimulator)
{
    /************************************************************************/
    /* Dynamic objects are created in this function.                        */
    /************************************************************************/
    m_DynamicDisksAreGeneratedInsideOfThisObject = true;

    m_DiskFileNameWthPath = dvdSimulator.m_DiskFileNameWthPath;
    m_CompTimeStatistics  = dvdSimulator.m_CompTimeStatistics;

    m_DynamicDisks        = dvdSimulator.create_dynamic_objects_from_current_dynamic_objects();
    m_AllContainers       = dvdSimulator.create_containers_from_current_containers();

    m_DynamicVD           = dvdSimulator.m_DynamicVD;

    change_user_data_of_generators_in_DVD_to_point_current_dynamic_objects_or_containers();
    make_map_from_new_disks_to_copied_disks_at_time_ZERO_in_DVD();
}


std::list<DynamicDisk*> DynamicVD2DSimulator::create_dynamic_objects_from_current_dynamic_objects() const
{
    list<DynamicDisk*> outputDynamicObjects = create_from(m_DynamicDisks);
    return outputDynamicObjects;
}


std::list<DynamicDisk*> DynamicVD2DSimulator::create_containers_from_current_containers() const
{
    list<DynamicDisk*> outputContainers = create_from(m_AllContainers);
    return outputContainers;
}


std::list<DynamicDisk*> DynamicVD2DSimulator::create_from(const list<DynamicDisk*>& input) const
{
    list<DynamicDisk*> outputDynamicObjects;

    for (list<DynamicDisk*>::const_iterator it_DynamicObject = input.begin();
        it_DynamicObject != input.end();
        ++it_DynamicObject)
    {
        DynamicDisk* currDynamicObject = *it_DynamicObject;
        DynamicDisk* newDyanmicObject  = create_dynamic_object(currDynamicObject);
        outputDynamicObjects.push_back(newDyanmicObject);
    }

    return outputDynamicObjects;
}


void DynamicVD2DSimulator::change_user_data_of_generators_in_DVD_to_point_current_dynamic_objects_or_containers()
{
    change_user_data_of_generators_in_DVD_to_point_current_one(m_DynamicDisks);
    change_user_data_of_generators_in_DVD_to_point_current_one(m_AllContainers);
}


void DynamicVD2DSimulator::change_user_data_of_generators_in_DVD_to_point_current_one(const list<DynamicDisk*>& input)
{
    for (list<DynamicDisk*>::const_iterator it_Object = input.begin();
         it_Object != input.end();
         ++it_Object)
    {
        DynamicDisk* currObject    = *it_Object;
        Generator2D* currGenerator =  m_DynamicVD.get_generator_by_its_input_disk_id(currObject->getID());
        currGenerator->setUserData(currObject);
    }
}


void DynamicVD2DSimulator::make_map_from_new_disks_to_copied_disks_at_time_ZERO_in_DVD()
{
    m_DynamicVD.make_map_from_original_disks_to_copied_disks(m_AllContainers, m_DynamicDisks);
}


void DynamicVD2DSimulator::construct_initial_CIC_VD()
{
    m_DynamicVD.construct_initial_VD(m_AllContainers, m_DynamicDisks);
}


void DynamicVD2DSimulator::generate_EVENTSEQ(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType, const bool& reinstateInitialVD /*= true*/)
{
    m_DynamicVD.generate_EVENTSEQ(lengthOfVoronoiClock, collisionHandlingType, reinstateInitialVD);
}


void DynamicVD2DSimulator::generate_EVENTSEQ(const int& numOfCollisions, const bool& reinstateInitialVD /*= true*/)
{
    m_DynamicVD.generate_EVENTSEQ(numOfCollisions, reinstateInitialVD);
}


void DynamicVD2DSimulator::generate_EVENTSEQ(const DiskCollisionWatcher& collisionWatcher, const bool& reinstateInitialVD /*=true*/)
{
    m_DynamicVD.generate_EVENTSEQ(collisionWatcher,
                                  reinstateInitialVD);
}


DiskCollisionEvent* DynamicVD2DSimulator::generate_EVENTSEQ_by_next_collision(const double& endVoronoiClock, const bool& reinstateInitialVD /*= true*/)
{
    return m_DynamicVD.generate_EVENTSEQ_by_next_collision(endVoronoiClock, reinstateInitialVD);
}


list<int> DynamicVD2DSimulator::find_neighbor_disk_ids(const int& anchorDiskId, const double& startVoronoiClock, const double& endVoronoiClock, const bool& thisIsContainer /*= false*/)
{
    Generator2D* anchorGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(anchorDiskId);
    list<Generator2D*> neighborGenerators = m_DynamicVD.find_neighbor_generators(anchorGenerator, startVoronoiClock, endVoronoiClock, thisIsContainer);

    list<int> neighborDiskIds;
    for (auto it_Generator = neighborGenerators.begin(); it_Generator != neighborGenerators.end(); ++it_Generator)
    {
        Generator2D* currGenerator = *it_Generator;
        DynamicDisk* currDisk      = static_cast<DynamicDisk*>(currGenerator->getUserData());
        neighborDiskIds.push_back(currDisk->getID());
    }

    return neighborDiskIds;
}


list<pair<int, double>> DynamicVD2DSimulator::find_future_collisions_if_target_disk_moves_differently(
    const int& targetDiskId, const rg_Point2D& velocityVector, const list<int>& otherDiskIds)
{
    const double startVoronoiClock = get_current_Voronoi_clock();

    list<pair<int, double>> output;

    Generator2D* targetGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(targetDiskId);
    DynamicDisk* targetDisk      = static_cast<DynamicDisk*>(targetGenerator->getUserData());

    const rg_Point2D originalVelocityVector = targetDisk->getVelocityVector();

    targetDisk->setVelocityVector(velocityVector);

    for (auto it_DiskId = otherDiskIds.begin(); it_DiskId != otherDiskIds.end(); ++it_DiskId)
    {
        const int    currDiskId = *it_DiskId;
        Generator2D* currGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(currDiskId);
        DynamicDisk* currDisk = static_cast<DynamicDisk*>(currGenerator->getUserData());

        pair<bool, double> collisionClock = m_DynamicVD.find_possible_collision_time_of_disks_after(startVoronoiClock, targetDisk, currDisk);

        if (collisionClock.first)
        {
            output.push_back(make_pair(currDiskId, collisionClock.second));

#ifdef DEBUG_DVD_SIM
            cout << targetDiskId << "\t" << currDiskId << "\t" << collisionClock.second << "\n";
#endif
        }
    }

    targetDisk->setVelocityVector(originalVelocityVector);

    output.sort(DynamicVD2DSimulator::compare_collision_clock);

    return output;
}


list<pair<int, double>> DynamicVD2DSimulator::find_collisions_with_neighbors_if_target_disk_moves_differently_at_current_clock(
    const int& targetDiskId, const rg_Point2D& velocityVector, const double& endVoronoiClock)
{
    list<pair<int, double>> output;

    Generator2D* targetGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(targetDiskId);
    DynamicDisk* targetDisk      = static_cast<DynamicDisk*>(targetGenerator->getUserData());

    const rg_Point2D originalVelocityVector = targetDisk->getVelocityVector();

    const double startVoronoiClock = get_current_Voronoi_clock();

    list<int> neighborDiskIds = find_neighbor_disk_ids(targetDiskId, startVoronoiClock, endVoronoiClock, targetDisk->this_disk_is_container());
    
    targetDisk->setVelocityVector(velocityVector);

    for (auto it_DiskId = neighborDiskIds.begin(); it_DiskId != neighborDiskIds.end(); ++it_DiskId)
    {
        const int    currDiskId    = *it_DiskId;
        Generator2D* currGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(currDiskId);
        DynamicDisk* currDisk      = static_cast<DynamicDisk*>(currGenerator->getUserData());

        pair<bool, double> collisionClock = m_DynamicVD.find_possible_collision_time_of_disks_after(startVoronoiClock, targetDisk, currDisk);

        if (collisionClock.first)
        {
            output.push_back(make_pair(currDiskId, collisionClock.second));

#ifdef DEBUG_DVD_SIM
            cout << targetDiskId << "\t" << currDiskId << "\t" << collisionClock.second << "\n";
#endif
        }
    }

    targetDisk->setVelocityVector(originalVelocityVector);

    output.sort(DynamicVD2DSimulator::compare_collision_clock);

    return output;
}




void DynamicVD2DSimulator::generate_EVENTSEQ(const double& lengthOfVoronoiClock, const bool& reinstateInitialVD /*= true*/)
{
    m_DynamicVD.generate_EVENTSEQ(lengthOfVoronoiClock, reinstateInitialVD);
}


void DynamicVD2DSimulator::generate_EVENTSEQ_collision_avoidance(const double& lengthOfVoronoiClock, const double& timeRelativeToCollisionTime, const bool& reinstateInitialVD /*= true*/)
{
	m_DynamicVD.generate_EVENTSEQ_collision_avoidance(lengthOfVoronoiClock, timeRelativeToCollisionTime, reinstateInitialVD);
}


void DynamicVD2DSimulator::increase_Voronoi_clock_either_to_FUTURE_or_to_PAST(const double& VoronoiClockIncrement /*= 1.0*/, const bool& VEdgeGeometryIsComputed /*= false*/)
{
	double targetVoronoiClock = get_current_Voronoi_clock() + VoronoiClockIncrement;
	go_to_target_clock_either_to_FUTURE_or_to_PAST(targetVoronoiClock);
}


void DynamicVD2DSimulator::go_to_target_clock_either_to_FUTURE_or_to_PAST(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed /*= false*/)
{
	double targetClock = targetVoronoiClock;

	if (targetVoronoiClock == get_current_Voronoi_clock())
	{
		return;
	}

	if (targetVoronoiClock > get_length_of_Voronoi_clock())
	{
        targetClock = get_length_of_Voronoi_clock();
	}
	else if (targetVoronoiClock < 0.0)
	{
        targetClock = 0.0;
	}

	___start_to_track_statistics_in_playing_simulation(targetClock);

    construct_VD_at_target_Voronoi_clock(targetClock, VEdgeGeometryIsComputed);

    ___stop_to_track_statistics_in_playing_simulation();
}


void DynamicVD2DSimulator::go_to_next_event_Voronoi_clock_either_to_FUTURE_or_to_PAST(const int& eventNumIncrement /*= 1*/, const bool& VEdgeGeometryIsComputed /*= false*/)
{
	if (eventNumIncrement == 0)
	{
		return;
	}
	else if (m_DynamicVD.get_most_recent_event_idx() == -1 && eventNumIncrement < 0)
	{
		return;
	}
	else if (m_DynamicVD.get_most_recent_event_idx() == m_DynamicVD.get_total_event_num() && eventNumIncrement > 0)
	{
		return;
	}

	if (eventNumIncrement > 0)
	{
		progress_DVD_using_EVENTNUM_to_FUTURE(eventNumIncrement);
	}
	else // (eventNumIncrement < 0)
	{
		progress_DVD_using_EVENTNUM_to_PAST(eventNumIncrement);
	}

	compute_all_VVertex_coordinates_of_VD();

	if (VEdgeGeometryIsComputed)
	{
		compute_VEdges_geometry_of_VD();
	}
}


void DynamicVD2DSimulator::reinstate_the_initial_Voronoi_diagram()
{
    m_DynamicVD.reinstate_the_initial_Voronoi_diagram();
}


void DynamicVD2DSimulator::compute_VEdges_geometry_of_VD()
{
    m_DynamicVD.compute_all_VEdges_equation_of_VD();
}



ReservationNumber DynamicVD2DSimulator::insert_disk_velocity_change_event_into_priorityQ(const double& eventTime, const int& inputDiskId, const rg_Point2D& velocityVector)
{
    Generator2D* generator = m_DynamicVD.get_generator_by_its_input_disk_id(inputDiskId);

    return m_DynamicVD.insert_disk_velocity_change_event_into_priorityQ(eventTime, generator, velocityVector);
}


void DynamicVD2DSimulator::remove_disk_velocity_change_event_from_priorityQ(const ReservationNumber& reservationNum)
{
    m_DynamicVD.remove_disk_velocity_change_event_from_priorityQ(reservationNum);
}


void DynamicVD2DSimulator::create_disks(const list<DynamicDisk>& input, list<DynamicDisk*>& output)
{
    for (list<DynamicDisk>::const_iterator it_Disk = input.begin();
         it_Disk != input.end();
         ++it_Disk)
    {
        const DynamicDisk& currDisk = *it_Disk;
        DynamicDisk* newDisk = create_disk(currDisk);
        output.push_back(newDisk);
    }
}


void DynamicVD2DSimulator::create_dynamic_disks(const list<DynamicDisk>& dynamicDisks)
{
    create_disks(dynamicDisks, m_DynamicDisks);
}


void DynamicVD2DSimulator::create_containers(const list<DynamicDisk>& containers)
{
    create_disks(containers, m_AllContainers);
}


void DynamicVD2DSimulator::delete_disks(list<DynamicDisk*>& disks)
{
    for (list<DynamicDisk*>::iterator it_Disk = disks.begin();
        it_Disk != disks.end();
        ++it_Disk)
    {
        delete_disk(*it_Disk);
    }

    disks.clear();
}


void DynamicVD2DSimulator::delete_dynamic_disks()
{
    delete_disks(m_DynamicDisks);
}


void DynamicVD2DSimulator::delete_containers()
{
    delete_disks(m_AllContainers);
}

void DynamicVD2DSimulator::load_disks(const string& fileNameWithFullPath)
{
    ___start_clock_in_DVD_simulator(DVDS_TIME_LOAD_DISKS);

    m_DynamicDisksAreGeneratedInsideOfThisObject = true;

    set_disk_file_name(fileNameWithFullPath);

    parse_input_file_and_create_dynamic_objects(fileNameWithFullPath);

    initialize_dynamic_disks_for_simulation();

    construct_initial_CIC_VD();

    ___end_clock_in_DVD_simulator(DVDS_TIME_LOAD_DISKS);
}


void DynamicVD2DSimulator::load_disks(const list<DynamicDisk*>& containers, 
                                      const list<DynamicDisk*>& dynamicDisks)
{
    ___start_clock_in_DVD_simulator(DVDS_TIME_LOAD_DISKS);

    set_containers(containers);

    set_dynamic_disks(dynamicDisks);

    initialize_dynamic_disks_for_simulation();

    construct_initial_CIC_VD();

    ___end_clock_in_DVD_simulator(DVDS_TIME_LOAD_DISKS);
}


void DynamicVD2DSimulator::load_disks(const list<DynamicDisk>& containers, const list<DynamicDisk>& dynamicDisks)
{
    ___start_clock_in_DVD_simulator(DVDS_TIME_LOAD_DISKS);

    m_DynamicDisksAreGeneratedInsideOfThisObject = true;

    create_containers(containers);

    create_dynamic_disks(dynamicDisks);

    initialize_dynamic_disks_for_simulation();

    construct_initial_CIC_VD();

    ___end_clock_in_DVD_simulator(DVDS_TIME_LOAD_DISKS);
}


void DynamicVD2DSimulator::initialize_dynamic_disks_for_simulation()
{
    for (list<DynamicDisk*>::const_iterator it_DynamicDisk = m_DynamicDisks.begin();
        it_DynamicDisk != m_DynamicDisks.end();
        it_DynamicDisk++)
    {
        DynamicDisk* dynamicDisk = *it_DynamicDisk;
        dynamicDisk->setMostRecentLocationUpdateTime(0.0);
    }
}

void DynamicVD2DSimulator::set_containers(const list<DynamicDisk*>& containers)
{
    m_AllContainers = containers;
}


void DynamicVD2DSimulator::set_dynamic_disks(const list<DynamicDisk*>& disks)
{
    m_DynamicDisks       = disks;
}


void DynamicVD2DSimulator::remove_disk(DynamicDisk* disk)
{
    m_DynamicVD.remove_disk(disk);
    m_DynamicDisks.remove(disk);
}


void DynamicVD2DSimulator::construct_initial_priority_queues()
{
    m_DynamicVD.construct_initial_priority_queues();
}


void DynamicVD2DSimulator::clear_priority_queues_for_collision_and_flip()
{
    m_DynamicVD.clear_priority_queues_for_collision_and_flip();
}


std::tuple<bool, double, DynamicDisk*, DynamicDisk*> DynamicVD2DSimulator::propagate_DVD_until_next_collision()
{
    return m_DynamicVD.propagate_DVD_until_next_collision();
}


void DynamicVD2DSimulator::go_back_to_target_Voronoi_clock_and_remove_events_after_the_target_Voronoi_clock(const double& targetVoronoiClock)
{
    m_DynamicVD.go_back_to_target_Voronoi_clock_and_remove_events_after_the_target_Voronoi_clock(targetVoronoiClock);
}

void DynamicVD2DSimulator::parse_input_file_and_create_dynamic_objects(const string& fileNameWithFullPath)
{
    list<DynamicDisk> containers;
    list<DynamicDisk> dynamicDisks;
    IOFunctionsForDisks::read(fileNameWithFullPath, containers, dynamicDisks);

    create_disks(containers, m_AllContainers);
    create_disks(dynamicDisks, m_DynamicDisks);
}


void DynamicVD2DSimulator::construct_VD_at_target_Voronoi_clock(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed)
{
    if (this_clock_is_future(targetVoronoiClock))
    {
        progress_DVD_to_FUTURE(targetVoronoiClock);
    }
    else if (this_clock_is_past(targetVoronoiClock))
    {
        progress_DVD_to_PAST(targetVoronoiClock);
    }

    compute_all_VVertex_coordinates_of_VD();

	if (VEdgeGeometryIsComputed)
	{
        compute_VEdges_geometry_of_VD();
	}
}


void DynamicVD2DSimulator::compute_all_VVertex_coordinates_of_VD()
{
    m_DynamicVD.compute_all_VVertex_coordinates_of_VD();
}


void DynamicVD2DSimulator::synchronize_disk_locations_and_current_clock_to_VD(const double& targetVoronoiClock)
{
    m_DynamicVD.synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVD2DSimulator::load_EVENTSEQ(const string& fileNameWithFullPath, const bool& reinstateInitialVD/*= true*/)
{
    ___start_clock_in_DVD_simulator(DVDS_TIME_LOAD_EVENTSEQ);

    vector<EventOfDynamicVD2D*>   eventSequence;
    double                        lengthOfVoronoiClock;
    double                        eventSeqFileVers;
    DiskCollisionHandler          collisionHandlingType;

    extract_DVD_data_from_EVENTSEQ_file(fileNameWithFullPath, eventSequence, lengthOfVoronoiClock, eventSeqFileVers, collisionHandlingType);

    set_input_parameters_to_DVD(lengthOfVoronoiClock, collisionHandlingType);

    set_EVENTSEQ_to_DVD(eventSequence, eventSeqFileVers, reinstateInitialVD);

    ___end_clock_in_DVD_simulator(DVDS_TIME_LOAD_EVENTSEQ);
}


void DynamicVD2DSimulator::load_EVENTSEQ(
    const string& fileNameWithFullPath, 
    const DiskCollisionWatcher& diskCollisionWatcher, 
    const bool& reinstateInitialVD/*= true*/ )
{
    ___start_clock_in_DVD_simulator(DVDS_TIME_LOAD_EVENTSEQ);

    vector<EventOfDynamicVD2D*>   eventSequence;
    double                        lengthOfVoronoiClock;
    double                        eventSeqFileVers;
    DiskCollisionHandler          collisionHandlingType;

    extract_DVD_data_from_EVENTSEQ_file(fileNameWithFullPath, eventSequence, lengthOfVoronoiClock, eventSeqFileVers, collisionHandlingType);

    set_input_parameters_to_DVD(lengthOfVoronoiClock, collisionHandlingType, diskCollisionWatcher);
    
    set_EVENTSEQ_to_DVD(eventSequence, eventSeqFileVers, reinstateInitialVD);

    ___end_clock_in_DVD_simulator(DVDS_TIME_LOAD_EVENTSEQ);
}


void DynamicVD2DSimulator::rewind_EVENTSEQ(const double& targetVoronoiClock)
{
    m_DynamicVD.rewind_EVENTSEQ(targetVoronoiClock);
}


void DynamicVD2DSimulator::extract_DVD_data_from_EVENTSEQ_file(
    const string& fileNameWithFullPath, 
    vector<EventOfDynamicVD2D*>& eventSequence, 
    double& lengthOfVoronoiClock, 
    double& eventSeqFileVers, 
    DiskCollisionHandler& collisionHandlingType) const
{
    ifstream fin(fileNameWithFullPath.c_str(), ios::in);

    if (!fin.is_open())
    {
        return;
    }

    parse_comment_parts_of_EVENTSEQ_file(fin, lengthOfVoronoiClock, eventSeqFileVers, collisionHandlingType);

    parse_event_definitions_of_of_EVENTSEQ_file(fin, eventSequence);

    fin.close();
}


void DynamicVD2DSimulator::set_input_parameters_to_DVD(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType)
{
    m_DynamicVD.set_parameters(lengthOfVoronoiClock, collisionHandlingType);
}


void DynamicVD2DSimulator::set_input_parameters_to_DVD(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType, const DiskCollisionWatcher& collisionWatcher)
{
    m_DynamicVD.set_parameters(lengthOfVoronoiClock, collisionWatcher, collisionHandlingType);
}

void V::GeometryTier::DynamicVD2DSimulator::set_EVENTSEQ_to_DVD(const vector<EventOfDynamicVD2D*>& EVENTSEQ, const double eventSeqFileVers, const bool& reinstateInitialVD/*= true*/)
{
    if (rg_EQ(eventSeqFileVers, 1.0, resNeg4))
    {
        m_DynamicVD.m_Vers1p0EVENTSEQFileUsed = true;
    }
    else
    {
        m_DynamicVD.m_Vers1p0EVENTSEQFileUsed = false;
    }

    m_DynamicVD.m_DuringEVENTSEQGeneration = false;

    m_DynamicVD.set_EVENTSEQ(EVENTSEQ, reinstateInitialVD);
}


void V::GeometryTier::DynamicVD2DSimulator::parse_comment_parts_of_EVENTSEQ_file(ifstream& fin, double& lengthOfVoronoiClock, double& eventSeqFileVers, DiskCollisionHandler& collisionHandlingType) const
{
	const char   CR('\r');
    const char   LF('\n');
    const char   TAB('\t');
    const char   BLANK(' ');

    string delimeters = string(1,CR) + string(1, LF) + string(1, TAB) + string(1, BLANK);

	string currentLine;
	for (int i = 0; i < 2; ++i)
	{
		std::getline(fin, currentLine, LF);
	}
  

    vector<string> tokens;

    //line 3
    std::getline(fin, currentLine, LF);
    StringFunctions::tokens_from_string(currentLine, delimeters, tokens);
    eventSeqFileVers = std::stod(tokens[1]);
    tokens.clear();

    //line 4
    std::getline(fin, currentLine, LF);

    //line 5
    std::getline(fin, currentLine, LF);
    StringFunctions::tokens_from_string(currentLine, delimeters, tokens);
	lengthOfVoronoiClock = std::stod(tokens[1]);
    tokens.clear();

    //line 6
    std::getline(fin, currentLine, LF);
    StringFunctions::tokens_from_string(currentLine, delimeters, tokens);
	collisionHandlingType.set_coefficient_of_restitution(std::stod(tokens[1]));
}


void DynamicVD2DSimulator::parse_event_definitions_of_of_EVENTSEQ_file(ifstream& fin, vector<EventOfDynamicVD2D*>& eventSequence) const
{
	const char   CR('\r');
    const char   LF('\n');
    const char   TAB('\t');
    const char   BLANK(' ');

    string delimeters = string(1,CR) + string(1, LF) + string(1, TAB) + string(1, BLANK);

    string currentLine;
    std::getline(fin, currentLine, LF);

    vector<string> tokens;
	StringFunctions::tokens_from_string(currentLine, delimeters, tokens);

    int numOfEvents = std::stoi(tokens[0]);
    eventSequence.resize(numOfEvents);

    for (int i = 0; i < numOfEvents; i++)
    {
		std::getline(fin, currentLine, LF);
        tokens.clear();
		StringFunctions::tokens_from_string(currentLine, delimeters, tokens);

        EventOfDynamicVD2D* currEvent  = parse_event_type_and_create_event(tokens[0].c_str());
        double               eventTime = parse_to_double(tokens[1].c_str());

        currEvent->set_occurring_time(eventTime);

        switch (currEvent->get_event_type())
        {
            case EventOfDynamicVD2D::EDGE_FLIP:
            {
                parse_and_set_edge_flip_event(dynamic_cast<EdgeFlipEvent*>(currEvent), tokens);
                break;
            }

			case EventOfDynamicVD2D::DISK_COLLISION:
			{
                parse_and_set_disk_collision_event(dynamic_cast<DiskCollisionEvent*>(currEvent), tokens);
				break;
			}

			case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
			{
                parse_and_set_disk_velocity_change_event(dynamic_cast<DiskVelocityChangeEvent*>(currEvent), tokens);
				break;
			}

			case EventOfDynamicVD2D::DISK_HOPPING:
			{
                parse_and_set_disk_hopping_event(dynamic_cast<DiskHoppingEvent*>(currEvent), tokens);
				break;
			}
        }

        eventSequence[i] = currEvent;
    }
}


EventOfDynamicVD2D* DynamicVD2DSimulator::parse_event_type_and_create_event(const char* token) const
{
    switch (token[0])
    {
        case EDGE_FLIP:
        {
            EdgeFlipEvent* edgeFlipEvent = new EdgeFlipEvent;
            edgeFlipEvent->set_whether_this_vedge_to_be_flipped(true);
            return edgeFlipEvent;
        }

        case SHADOW_EDGE_FLIP:
        {
			EdgeFlipEvent* edgeFlipEvent = new EdgeFlipEvent;
			edgeFlipEvent->set_whether_this_vedge_to_be_flipped(false);
			return edgeFlipEvent;
        }

        case COLLISION:
        {
			DiskCollisionEvent* diskCollisionEvent = new DiskCollisionEvent;
			return diskCollisionEvent;
        }

		case VELOCITY_CHANGE:
		{
			DiskVelocityChangeEvent* velocityChangeEvent = new DiskVelocityChangeEvent;
			return velocityChangeEvent;
		}

        case TELEPORT_BTW_CONTAINERS:
        {
			DiskHoppingEvent* diskHoppingEvent = new DiskHoppingEvent;
            diskHoppingEvent->set_does_this_event_happen(true);
			return diskHoppingEvent;
        }

        case BOUNCE_SHADOW_TELPORT:
        {
			DiskHoppingEvent* diskHoppingEvent = new DiskHoppingEvent;
			diskHoppingEvent->set_does_this_event_happen(false);
			return diskHoppingEvent;
        }

        default:
        {
            return NULL;
        }
    }
}




double DynamicVD2DSimulator::parse_to_double(const char* token) const
{
#ifdef WRITE_HEXA_VALUES_FOR_DOUBLE_TYPE
	double eventTime = get_double_type_value_from_hexa_characters_in_order_of_little_endian(const_cast<char*>(token));
#else
	double eventTime = atof(const_cast<char*>(token));
#endif

    return eventTime;
}


void DynamicVD2DSimulator::parse_and_set_edge_flip_event(EdgeFlipEvent* edgeFlipEvent, const vector<string>& tokens) const
{
    int IDsOfFourInputDisksDefiningVEdge[4];

    IDsOfFourInputDisksDefiningVEdge[0] = stoi(tokens[2]);
	IDsOfFourInputDisksDefiningVEdge[1] = stoi(tokens[3]);
	IDsOfFourInputDisksDefiningVEdge[2] = stoi(tokens[4]);
	IDsOfFourInputDisksDefiningVEdge[3] = stoi(tokens[5]);

    EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = m_DynamicVD.get_generator_quadruplet_by_its_id(IDsOfFourInputDisksDefiningVEdge);
    edgeFlipEvent->set_generator_quadruplet_before_flip(generatorPtrQuadruplet);
}


void DynamicVD2DSimulator::parse_and_set_disk_collision_event(DiskCollisionEvent* diskCollisionEvent, const vector<string>& tokens) const
{
	int IDsOfCollidingGenerator[2];
    IDsOfCollidingGenerator[0] = stoi(tokens[2]);
    IDsOfCollidingGenerator[1] = stoi(tokens[3]);

	Generator2D* generator1OfCollisionEvent = m_DynamicVD.get_generator_by_its_input_disk_id(IDsOfCollidingGenerator[0]);
	Generator2D* generator2OfCollisionEvent = m_DynamicVD.get_generator_by_its_input_disk_id(IDsOfCollidingGenerator[1]);

    diskCollisionEvent->set_two_collided_generators(generator1OfCollisionEvent, generator2OfCollisionEvent);

    if (tokens.size() == 8)
    {
        double disk1_X   = parse_to_double(tokens[4].c_str());
        double disk1_Y   = parse_to_double(tokens[5].c_str());
        double disk2_X   = parse_to_double(tokens[6].c_str());
        double disk2_Y   = parse_to_double(tokens[7].c_str());

        diskCollisionEvent->set_disk_centers_at_collision_time(rg_Point2D(disk1_X, disk1_Y), rg_Point2D(disk2_X, disk2_Y));
    }
}


void DynamicVD2DSimulator::parse_and_set_disk_velocity_change_event(DiskVelocityChangeEvent* diskVelocityChangeEvent, const vector<string>& tokens) const
{
    int IDOfTargetDisk      = stoi(tokens[2]);
	double x_velocityVector = parse_to_double(tokens[3].c_str());
	double y_velocityVector = parse_to_double(tokens[4].c_str());

    Generator2D* targetGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(IDOfTargetDisk);
    diskVelocityChangeEvent->set_velocity_changed_generator(targetGenerator);
    diskVelocityChangeEvent->set_velocity_vector_after_change(rg_Point2D(x_velocityVector, y_velocityVector));

    if (tokens.size() == 7)
    {
        double x_pos = parse_to_double(tokens[5].c_str());
        double y_pos = parse_to_double(tokens[6].c_str());

        diskVelocityChangeEvent->set_disk_center_at_velocity_change_time(rg_Point2D(x_pos, y_pos));
    }
}


void DynamicVD2DSimulator::parse_and_set_disk_hopping_event(DiskHoppingEvent* diskHoppingEvent, const vector<string>& tokens) const
{
	int IDOfTargetDisk    = stoi(tokens[2]);
    int IDOfFromContainer = stoi(tokens[3]);
	int IDOfToContainer   = stoi(tokens[4]);
    
    double x_beforeHopping = parse_to_double(tokens[5].c_str());
	double y_beforeHopping = parse_to_double(tokens[6].c_str());
	double x_afterHopping  = parse_to_double(tokens[7].c_str());
	double y_afterHopping  = parse_to_double(tokens[8].c_str());

	Generator2D* targetDiskGenerator    = m_DynamicVD.get_generator_by_its_input_disk_id(IDOfTargetDisk);
	Generator2D* fromContainerGenerator = m_DynamicVD.get_generator_by_its_input_disk_id(IDOfFromContainer);
	Generator2D* toContainerGenerator   = m_DynamicVD.get_generator_by_its_input_disk_id(IDOfToContainer);

	diskHoppingEvent->set_target_generator(targetDiskGenerator);
	diskHoppingEvent->set_in_what_container_this_generator_is_deleted_n_inserted(fromContainerGenerator, toContainerGenerator);
	diskHoppingEvent->set_deleted_n_inserted_position(rg_Point2D(x_beforeHopping, y_beforeHopping), rg_Point2D(x_afterHopping, y_afterHopping));
	diskHoppingEvent->set_velocity_vector_setting(false);
}


string DynamicVD2DSimulator::create_EVENTSEQ_file_name() const
{
    string filePath = StringFunctions::getPathFromFileName(m_DiskFileNameWthPath);
    string fileName = StringFunctions::getFileNameWithoutPathAndExtension(m_DiskFileNameWthPath);

	string historyFileName = filePath + fileName
		                   + string("_T")  + to_string((int)get_length_of_Voronoi_clock())
		                   + string("_RC") + StringFunctions::convertDoubleToString(get_coefficient_of_restitution(), 2)
		                   + string(".evtsq");

	return historyFileName;
}


std::string DynamicVD2DSimulator::create_EVNETSEQ_generation_statistics_file_name() const
{
    string filePath = StringFunctions::getPathFromFileName(m_DiskFileNameWthPath);
    string fileName = StringFunctions::getFileNameWithoutPathAndExtension(m_DiskFileNameWthPath);

	string historyFileName = filePath + fileName
		                   + string("_T")  + to_string((int)get_length_of_Voronoi_clock())
		                   + string("_RC") + StringFunctions::convertDoubleToString(get_coefficient_of_restitution(), 2)
                           + string("_EVENTSEQ_STATISTICS.txt");

	return historyFileName;
}

void DynamicVD2DSimulator::write_EVENTSEQ(const string& fileNameWithFullPath)
{
    //Note: VD should be at t = 0.
    if (m_DynamicVD.get_total_event_num() == 0)
    {
        return;
    }

    ___start_clock_in_DVD_simulator(DVDS_TIME_WRITE_EVENTSEQ);

    ofstream historyFile(fileNameWithFullPath);
    historyFile.setf(std::ios::fixed, std::ios::floatfield);
    historyFile.precision(2);
 
    char       currTime[80];
    time_t     now = time(0);
    struct tm  tstruct;
  
    localtime_r(&now, &tstruct);

    // for more information about date/time format
    strftime(currTime, sizeof(currTime), "%c", &tstruct);

    historyFile << "#SOURCE:"                             << "\t"  << m_DiskFileNameWthPath                         << endl;
    historyFile << "#DVD_GENERATOR:"                      << "\t"  << DVD_2D_VERSION                                << endl;
    historyFile << "#EVENT_SEQUENCE_FILE_FORMAT_VERSION:"  << "\t" << EVENT_SEQUENCE_FILE_FORMAT_VERSION             << endl;
    historyFile << "#TIME:"                               << "\t"  << currTime                                      <<endl;
    historyFile << "#SIMULATION_TIME:"                    << "\t"  << m_DynamicVD.get_length_of_Voronoi_clock() << endl;
    historyFile << "#RESTITUTION_COEFFICIENT:"            << "\t"  << m_DynamicVD.get_coefficient_of_restitution()  << endl;
    
    historyFile.precision(15);
	historyFile << m_DynamicVD.get_total_event_num() << "\n";

    // Store current status of time, event index
    const double currentTime            = get_current_Voronoi_clock();

    // Make VD at the time '0'
    if (currentTime != 0.0)
    {
        reinstate_the_initial_Voronoi_diagram();
        //go_to_target_clock_either_to_FUTURE_or_to_PAST(0.0);
    }
    
    // Scan event and store it
    vector<EventOfDynamicVD2D*> eventSequence;
    m_DynamicVD.get_EVENTSEQ(eventSequence);

    for (vector<EventOfDynamicVD2D*>::const_iterator it_Event = eventSequence.begin();
         it_Event != eventSequence.end();
         it_Event++)
    {
        EventOfDynamicVD2D* currEvent = const_cast<EventOfDynamicVD2D*>(*it_Event);
        switch (currEvent->get_event_type())
        {
            case EventOfDynamicVD2D::EDGE_FLIP:
            {
                write_edge_flipping_event(currEvent, historyFile);
                break;
            }

            case EventOfDynamicVD2D::DISK_COLLISION:
            {
                write_disk_collision_event(currEvent, historyFile);
                break;
            }

			case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
			{
                write_disk_velocity_change_event(currEvent, historyFile);
				break;
			}

            case EventOfDynamicVD2D::DISK_HOPPING:
            {
                write_disk_hopping_event(currEvent, historyFile);
                break;
            }
        }
    }
    
	historyFile.close();

    // Make VD at the stored time 
    if (currentTime != 0.0)
    {
        go_to_target_clock_either_to_FUTURE_or_to_PAST(currentTime);
    }

    ___end_clock_in_DVD_simulator(DVDS_TIME_WRITE_EVENTSEQ);
}


void DynamicVD2DSimulator::write_EVENTSEQ()
{
    string eventSequenceFileNameWithPath 
        = create_EVENTSEQ_file_name();

    write_EVENTSEQ(eventSequenceFileNameWithPath);
}


void DynamicVD2DSimulator::write_statistics_for_generating_EVENTSEQ(const string& fileNameWithPath) const
{
    m_DynamicVD.write_statistics_for_generating_EVENTSEQ(fileNameWithPath, false);
}



void DynamicVD2DSimulator::write_statistics_for_generating_EVENTSEQ() const
{
    string eventSequenceStatisticsFileNameWithPath 
        = create_EVNETSEQ_generation_statistics_file_name();

    write_statistics_for_generating_EVENTSEQ(
        eventSequenceStatisticsFileNameWithPath);
}


void DynamicVD2DSimulator::write_data_loading_time_statistics(const string& fileNameWithPath) const
{
    ofstream fout;
    fout.open(fileNameWithPath.c_str());

    fout << "Load Disks: " << "\t" << get_DVDSIM_time_statistics().timeInSec(DVDS_TIME_LOAD_DISKS)       << endl;
    fout << "Load EVENTSEQ: " << "\t" << get_DVDSIM_time_statistics().timeInSec(DVDS_TIME_LOAD_EVENTSEQ) << endl;

    fout.close();
}


void DynamicVD2DSimulator::write_disk_contact_numbers_by_time_increment(const string& fileNameWithPath) const
{
    ofstream fout;
    fout.open(fileNameWithPath.c_str());

    fout << "Order" << "\t" <<  "Start" << "\t" <<  "End" << "\t" << "counts" << endl;
 
    const CountStatistics& agentContactNumberStatistics 
        = m_DynamicVD.get_count_statistics(PHASE2_COLLISION_BTW_DISKS_COUNT);

    const TimeStatistics<milli>& VoronoiClockIncrement
        = m_DynamicVD.get_time_statistics(PHASE2_SIMULATION_TIME_INCREMNT);


    float currStartClock = 0.0;
    for (int j = 0; j < agentContactNumberStatistics.size(); ++j)
    {
        const int    contactNumbers       = agentContactNumberStatistics.count(j);
        const float  VoronoiClockInterval = VoronoiClockIncrement.time(j);
        const float  currEndClock         = currStartClock + VoronoiClockInterval;

        fout << j + 1          << "\t" 
             << currStartClock << "\t" 
             << currEndClock   << "\t"  
             << contactNumbers << endl;

        currStartClock = currEndClock;
    }

    fout << endl;
    fout.close();
}


void DynamicVD2DSimulator::load_disk_hopping_probability_matrix(const string& fileNameWithFullPath)
{
    vector<float> diskHoppingProbabilityMatrix;

    ifstream  DHPMFileStream(fileNameWithFullPath.c_str(), ios::in);

    if (!DHPMFileStream.is_open())
    {
        return;
    }

    char	lineOfData[300];
    char    seps[4] = { ' ', '\t', '\n', '\0' };
    char*   token   = NULL;
    char*   context = NULL;


    DHPMFileStream.getline(lineOfData, 100);
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    token = strtok_s(lineOfData, seps, &context);
#else
    token = strtok_r(lineOfData, seps, &context);
#endif

    int numOfContainer = atoi(token);

    for (int i = 0; i < numOfContainer; ++i)
    {
        DHPMFileStream.getline(lineOfData, 100);

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
        token = strtok_s(lineOfData, seps, &context);
#else
        token = strtok_r(lineOfData, seps, &context);
#endif

        float probability = atof(token);
        diskHoppingProbabilityMatrix.push_back(probability);

        for (int j = 0; j < numOfContainer - 1; ++j)
        {
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
            token = strtok_s(NULL, seps, &context);
#else
            token = strtok_r(NULL, seps, &context);
#endif
            probability = atof(token);
            diskHoppingProbabilityMatrix.push_back(probability);
        }
    }

    m_DynamicVD.set_disk_hopping_probability_matrix(diskHoppingProbabilityMatrix);
}


//for reading and writing history file
void DynamicVD2DSimulator::write_edge_flipping_event(EventOfDynamicVD2D * anEvent, ofstream & fileStream) const
{
    EdgeFlipEvent* edgeFlippingEvent = static_cast<EdgeFlipEvent*>(anEvent);
    EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = edgeFlippingEvent->get_generator_quadruplet_before_flip();

    DynamicDisk* leftDisk  = (DynamicDisk*)(generatorPtrQuadruplet[0]->getUserData());
    DynamicDisk* rightDisk = (DynamicDisk*)(generatorPtrQuadruplet[1]->getUserData());
    DynamicDisk* headDisk  = (DynamicDisk*)(generatorPtrQuadruplet[2]->getUserData());
    DynamicDisk* tailDisk  = (DynamicDisk*)(generatorPtrQuadruplet[3]->getUserData());

    int IDOfLeftDisk  = leftDisk->getID();
    int IDOfRightDisk = rightDisk->getID();
    int IDOfHeadDisk  = headDisk->getID();
    int IDOfTailDisk  = tailDisk->getID();

    char eventSymbol;
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        eventSymbol = EDGE_FLIP;
    }
    else
    {
        eventSymbol = SHADOW_EDGE_FLIP;
    }

#ifdef WRITE_HEXA_VALUES_FOR_DOUBLE_TYPE
    const int sizeOfDouble = sizeof(double);
    unsigned char charOfEventTime[sizeOfDouble];
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(edgeFlippingEvent->get_occurring_clock(), charOfEventTime);

    fileStream << eventSymbol << "\t";

    write_hexa_characters(fileStream, sizeOfDouble, charOfEventTime);

    fileStream << dec           << "\t"
               << IDOfLeftDisk  << "\t"
               << IDOfRightDisk << "\t"
               << IDOfHeadDisk  << "\t"
               << IDOfTailDisk  << "\n";
#else

    fileStream << eventSymbol                             << "\t"
               << edgeFlippingEvent->get_occurring_clock() << "\t"
               << IDOfLeftDisk                            << "\t"
               << IDOfRightDisk                           << "\t"
               << IDOfHeadDisk                            << "\t"
               << IDOfTailDisk                            << "\n";
#endif
}


void DynamicVD2DSimulator::write_disk_collision_event(EventOfDynamicVD2D * anEvent, ofstream & fileStream) const
{
    DiskCollisionEvent* diskCollisionEvent = static_cast<DiskCollisionEvent*>(anEvent);

    int IDOfDisk1 = ((DynamicDisk*)(diskCollisionEvent->get_disk_generator1()->getUserData()))->getID();
    int IDOfDisk2 = ((DynamicDisk*)(diskCollisionEvent->get_disk_generator2()->getUserData()))->getID();

    char eventSymbol = COLLISION;

#ifdef WRITE_HEXA_VALUES_FOR_DOUBLE_TYPE

    const int sizeOfDouble = sizeof(double);
    unsigned char charOfEventTime[sizeOfDouble];
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskCollisionEvent->get_occurring_clock(), charOfEventTime);

    unsigned char locationPointElement[4][sizeOfDouble];
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskCollisionEvent->get_center_of_disk1().getX(), locationPointElement[0]);
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskCollisionEvent->get_center_of_disk1().getY(), locationPointElement[1]);
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskCollisionEvent->get_center_of_disk2().getX(), locationPointElement[2]);
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskCollisionEvent->get_center_of_disk2().getY(), locationPointElement[3]);

	fileStream << eventSymbol << "\t";

	write_hexa_characters(fileStream, sizeOfDouble, charOfEventTime);

    fileStream << dec       << "\t"
               << IDOfDisk1 << "\t"
               << IDOfDisk2 << "\t";

    for (int i = 0; i < 4; ++i)
    {
        write_hexa_characters(fileStream, sizeOfDouble, locationPointElement[i]);
        fileStream << "\t";
    }
    fileStream << "\n";

#else
    fileStream << eventSymbol                               << "\t" 
               << diskCollisionEvent->get_occurring_clock() << "\t"
               << IDOfDisk1                                 << "\t"
               << IDOfDisk2                                 << "\t"
               << diskCollisionEvent->get_center_of_disk1().getX() << "\t"
               << diskCollisionEvent->get_center_of_disk1().getY() << "\t"
               << diskCollisionEvent->get_center_of_disk2().getX() << "\t"
               << diskCollisionEvent->get_center_of_disk2().getY() << "\n"
#endif
   
}


void DynamicVD2DSimulator::write_disk_velocity_change_event(EventOfDynamicVD2D* anEvent, ofstream& fileStream) const
{
	DiskVelocityChangeEvent* velocityChangeEvent = static_cast<DiskVelocityChangeEvent*>(anEvent);

	int IDOfTargetDisk     = (static_cast<DynamicDisk*>(velocityChangeEvent->get_disk_generator()->getUserData()))->getID();

	char eventSymbol = VELOCITY_CHANGE;

#ifdef WRITE_HEXA_VALUES_FOR_DOUBLE_TYPE
	const int sizeOfDouble = sizeof(double);

	unsigned char charOfEventTime[sizeOfDouble];
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(velocityChangeEvent->get_occurring_clock(), charOfEventTime);

	unsigned char velocityVector[2][sizeOfDouble];
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(velocityChangeEvent->get_future_velocity_vector_of_disk().getX(), velocityVector[0]);
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(velocityChangeEvent->get_future_velocity_vector_of_disk().getY(), velocityVector[1]);

    unsigned char locationPointElement[2][sizeOfDouble];
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(velocityChangeEvent->get_center_of_disk().getX(), locationPointElement[0]);
    get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(velocityChangeEvent->get_center_of_disk().getY(), locationPointElement[1]);

    fileStream << eventSymbol << "\t";
	
    write_hexa_characters(fileStream, sizeOfDouble, charOfEventTime);
	
    fileStream << dec << "\t" << IDOfTargetDisk << "\t";
	
	for (int i = 0; i < 2; ++i)
	{
		write_hexa_characters(fileStream, sizeOfDouble, velocityVector[i]);

		fileStream << "\t";
	}

    for (int i = 0; i < 2; ++i)
    {
        write_hexa_characters(fileStream, sizeOfDouble, locationPointElement[i]);

        fileStream << "\t";
    }

    fileStream << "\n";
#else
    fileStream << eventSymbol                                                      << "\t" 
               << velocityChangeEvent->get_occurring_clock()                       << "\t" 
               << IDOfTargetDisk                                                   << "\t"
		       << velocityChangeEvent->get_future_velocity_vector_of_disk().getX() << "\t"
               << velocityChangeEvent->get_future_velocity_vector_of_disk().getY() << "\t"
               << velocityChangeEvent->get_center_of_disk().getX()                 << "\t"
               << velocityChangeEvent->get_center_of_disk().getY()                 << "\n"
#endif
}


void DynamicVD2DSimulator::write_disk_hopping_event(EventOfDynamicVD2D* anEvent, ofstream& fileStream) const
{
    DiskHoppingEvent* diskHoppingEvent = static_cast<DiskHoppingEvent*>(anEvent);

    int IDOfTargetDisk    = (static_cast<DynamicDisk*>(diskHoppingEvent->get_target_generator()->getUserData()))->getID();
    int IDOfFromContainer = (static_cast<DynamicDisk*>(diskHoppingEvent->get_in_what_container_this_generator_is_deleted()->getUserData()))->getID();
    int IDOfToContainer   = (static_cast<DynamicDisk*>(diskHoppingEvent->get_in_what_container_this_generator_is_inserted()->getUserData()))->getID();

    char eventSymbol;
    if (diskHoppingEvent->get_does_this_event_happen())
    {
        eventSymbol = TELEPORT_BTW_CONTAINERS;
    }
    else
    {
        eventSymbol = BOUNCE_SHADOW_TELPORT;
    }

#ifdef WRITE_HEXA_VALUES_FOR_DOUBLE_TYPE

	const int sizeOfDouble = sizeof(double);
	unsigned char charOfEventTime[sizeOfDouble];
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskHoppingEvent->get_occurring_clock(), charOfEventTime);

	unsigned char locationPointElement[4][sizeOfDouble];
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskHoppingEvent->get_deleted_position().getX(), locationPointElement[0]);
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskHoppingEvent->get_deleted_position().getY(), locationPointElement[1]);
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskHoppingEvent->get_inserted_position().getX(), locationPointElement[2]);
	get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(diskHoppingEvent->get_inserted_position().getY(), locationPointElement[3]);


    
    fileStream << eventSymbol << "\t";

    write_hexa_characters(fileStream, sizeOfDouble, charOfEventTime);

	fileStream << dec << "\t" << IDOfTargetDisk
		              << "\t" << IDOfFromContainer
		              << "\t" << IDOfToContainer
		              << "\t";


	for (int i = 0; i < 4; ++i)
	{
		write_hexa_characters(fileStream, sizeOfDouble, locationPointElement[i]);

		fileStream << "\t";
	}

    fileStream << "\n";
#else
    fileStream << eventSymbol << "\t"
               << diskHoppingEvent->get_occurring_clock()           << "\t"
		       << IDOfTargetDisk                                     << "\t"
		       << IDOfFromContainer                                  << "\t"
		       << IDOfToContainer                                    << "\t"
               << diskHoppingEvent->get_deleted_position().getX()  << "\t"
		       << diskHoppingEvent->get_deleted_position().getY()  << "\t"
		       << diskHoppingEvent->get_inserted_position().getX() << "\t"
               << diskHoppingEvent->get_inserted_position().getY() << "\n"
#endif

}


void DynamicVD2DSimulator::write_average_distances_for_one_disks_to_first_neighbors(const string& fileNameWithPath) const
{
    ofstream fout(fileNameWithPath.c_str());
    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(3);

    list<Generator2D*> allGenerators;
    m_DynamicVD.getGenerators(allGenerators);

    for (list<Generator2D*>::const_iterator it_Generator = allGenerators.begin();
        it_Generator != allGenerators.end();
        it_Generator++)
    {
        Generator2D* generator = *it_Generator;

        //container have ID, -1.
        if (generator->getID() == -1)
        {
            continue;
        }

        double avgDistanceFromCurrDiskToNeighbors = 0.0;
        int    numOfFirstNeighbors = 0;

        list<Generator2D*> neighborGenerators;
        generator->getNeighborGenerators(neighborGenerators, false);
        rg_Circle2D currDisk = generator->getDisk();

        for (list<Generator2D*>::const_iterator it_NeighborGenerator = neighborGenerators.begin();
            it_NeighborGenerator != neighborGenerators.end();
            it_NeighborGenerator++)
        {
            Generator2D* neighborGenerator = *it_NeighborGenerator;

            //container have ID, -1.
            if (neighborGenerator->getID() == -1)
            {
                continue;
            }

            rg_Circle2D neighborDisk = neighborGenerator->getDisk();

            if (!currDisk.isIntersectWith(neighborDisk))
            {
                avgDistanceFromCurrDiskToNeighbors += currDisk.distance(neighborDisk);
                numOfFirstNeighbors++;
            }
        }

        avgDistanceFromCurrDiskToNeighbors = avgDistanceFromCurrDiskToNeighbors / numOfFirstNeighbors;

        fout << avgDistanceFromCurrDiskToNeighbors << endl;
    }

    fout.close();
}



void DynamicVD2DSimulator::write_speed_of_all_disks(const string& fileNameWithPath) const
{
    ofstream fout(fileNameWithPath.c_str());
    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(3);

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    for (list<DynamicDisk*>::const_iterator it_DynamicDisk = dynamicDisks.begin();
        it_DynamicDisk != dynamicDisks.end();
        it_DynamicDisk++)
    {
        const DynamicDisk* disk = (*it_DynamicDisk);

        double speed = disk->getVelocityVector().magnitude();

        fout << speed << endl;
    }

    fout.close();
}



double DynamicVD2DSimulator::compute_density() const
{
    double containerArea  =  pow(get_outer_most_container()->getCircle().getRadius(), 2);
    double sumOfDisksArea = 0.0;

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    for (list<DynamicDisk*>::const_iterator it_MovingDisk = dynamicDisks.begin(); it_MovingDisk != dynamicDisks.end(); it_MovingDisk++)
    {
        double currRadius = (*it_MovingDisk)->getCircle().getRadius();

        sumOfDisksArea += pow(currRadius, 2);
    }

    return sumOfDisksArea / containerArea;
}



rg_Point2D DynamicVD2DSimulator::compute_current_momentum() const
{
    rg_Point2D sumOfMomentum(0.0,0.0);

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    for (list<DynamicDisk*>::const_iterator it_MovingDisk = dynamicDisks.begin();
        it_MovingDisk != dynamicDisks.end();
        it_MovingDisk++)
    {
        const DynamicDisk* currMovingDisk = *it_MovingDisk;

        const double      mass           = pow(currMovingDisk->getCircle().getRadius(), 2) * rg_PI;
        const rg_Point2D& velocityVector = currMovingDisk->getVelocityVector();

        sumOfMomentum = sumOfMomentum + mass * velocityVector;
    }

    return sumOfMomentum;
}


double DynamicVD2DSimulator::compute_current_average_velocity_magnitude() const
{
    double sumOfVelocityMagnitude = 0.0;

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    for (list<DynamicDisk*>::const_iterator it_MovingDisk = dynamicDisks.begin();
        it_MovingDisk != dynamicDisks.end();
        it_MovingDisk++)
    {
        const DynamicDisk* currMovingDisk = *it_MovingDisk;
        
        sumOfVelocityMagnitude += currMovingDisk->getVelocityVector().magnitude();
    }

    return sumOfVelocityMagnitude / dynamicDisks.size();
}


double DynamicVD2DSimulator::compute_current_energy() const
{
    double sumOfEnergy = 0.0;

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    for (list<DynamicDisk*>::const_iterator it_MovingDisk = dynamicDisks.begin();
        it_MovingDisk != dynamicDisks.end();
        it_MovingDisk++)
    {
        const DynamicDisk* currMovingDisk = *it_MovingDisk;

        const double      mass              = pow(currMovingDisk->getCircle().getRadius(), 2) * rg_PI;
        const double      velocityMagSquare = currMovingDisk->getVelocityVector().magnitudeSquare();

        sumOfEnergy = sumOfEnergy + 1.0/2.0 * mass * velocityMagSquare;
    }

    return sumOfEnergy;
}


bool DynamicVD2DSimulator::data_set_is_valid()
{
    if (there_is_a_pair_of_disks_intersected_by_VD())
    {
        cout << "There is at least a pair of moving disks intersected" << endl;
        return false;
    }

    if (all_disks_are_stopped())
    {
        cout << "There is no moving Disk" << endl;
        return false;
    }

    return true;
}


bool DynamicVD2DSimulator::there_is_a_pair_of_disks_intersected_by_VD(const double& offset /*=0.0*/, const double& tolerance /*= rg_MATH_RES*/)
{
    //This works correctly when only one container exist!!!


    ___start_clock_in_DVD_simulator(DVDS_TIME_CHECK_WHETHER_DISKS_ARE_INTERSECTED);

#ifdef USING_BUCKET_FOR_CHECKING_ALL_PAIR_WISE_DISK_INTERSECTION
    
    list<DynamicDisk*> inputDynamicDisks;
    get_dynamic_disks(inputDynamicDisks);
    
    list<rg_Circle2D> inputDisksReflectingOffset;
    for(list<DynamicDisk*>::const_iterator it_DynamicDisk = inputDynamicDisks.begin();
        it_DynamicDisk != inputDynamicDisks.end();
        ++it_DynamicDisk)
    {
        DynamicDisk* dynamicDisk = (*it_DynamicDisk);
        if (dynamicDisk != get_most_outer_container())
        {
            rg_Circle2D newDisk(dynamicDisk->getCircle().getCenterPt(), dynamicDisk->getCircle().getRadius() + offset);
            inputDisksReflectingOffset.push_back(newDisk);
        }
    }

    BucketForDisk bucketForDisks(inputDisksReflectingOffset);
    bool thereIsIntersection = bucketForDisks.there_is_a_pair_of_disks_intersected(TOLERANCE_OF_DISK_INTERSECTION);

    return thereIsIntersection;
#else

    bool thereIsIntersection = false;

    list<VEdge2D*> allVEdges;
    if (m_AllContainers.size() == 1)
    {
        m_DynamicVD.getVoronoiEdges(allVEdges);
    }
    else
    {
        m_DynamicVD.getAllEdges_hierarchy(allVEdges);
    }

    for (list<VEdge2D*>::const_iterator it_VEdge = allVEdges.begin(); it_VEdge != allVEdges.end(); it_VEdge++)
    {
        if ((*it_VEdge)->getStatus() == RED_E)
        {
            continue;
        }

        DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>((*it_VEdge)->getRightFace()->getGenerator()->getUserData());
        DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>((*it_VEdge)->getLeftFace()->getGenerator()->getUserData());

        const rg_Circle2D& disk1 = movingDisk1->getCircle();
        const rg_Circle2D& disk2 = movingDisk2->getCircle();

        if (movingDisk1->this_disk_is_container())
        {
            thereIsIntersection = !disk2.isIncludedIn(disk1, TOLERANCE_OF_DISK_INTERSECTION);
        }
        else if (movingDisk2->this_disk_is_container())
        {
            thereIsIntersection = !disk1.isIncludedIn(disk2, TOLERANCE_OF_DISK_INTERSECTION);
        }
        else
        {
            thereIsIntersection = disk1.isIntersectWith(disk2, TOLERANCE_OF_DISK_INTERSECTION);
        }

        if (thereIsIntersection)
        {
            break;
        }
    }

    return thereIsIntersection;
#endif

    ___end_clock_in_DVD_simulator(DVDS_TIME_CHECK_WHETHER_DISKS_ARE_INTERSECTED);
}



bool DynamicVD2DSimulator::all_disks_are_stopped() const
{
    int numOfNotMovingDisks = 0;

    for (list<DynamicDisk*>::const_iterator it_DynamicDisk = m_DynamicDisks.begin();
         it_DynamicDisk != m_DynamicDisks.end();
         it_DynamicDisk++)
    {
        const DynamicDisk* currMovingDisk = *it_DynamicDisk;

        if (rg_EQ(currMovingDisk->getVelocityVectorX(), 0.0f, TOLERANCE_OF_ZERO_VELOCITY) 
         && rg_EQ(currMovingDisk->getVelocityVectorY(), 0.0f, TOLERANCE_OF_ZERO_VELOCITY))
        {
            numOfNotMovingDisks++;
        }
    }

    if (numOfNotMovingDisks == m_DynamicDisks.size())
    {
        return true;
    }
    else
    {
        return false;
    }
}



void DynamicVD2DSimulator::progress_DVD_to_FUTURE(const double& targetVoronoiClock)
{
    while (next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(targetVoronoiClock))
    {
        handle_next_DVD_event_to_FUTURE();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVD2DSimulator::progress_DVD_to_PAST(const double& targetVoronoiClock)
{
    while (  next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(targetVoronoiClock)
          || next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock())
    {
        reversely_handle_next_DVD_event_to_PAST();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVD2DSimulator::progress_DVD_using_EVENTNUM_to_FUTURE(const int& eventNumIncrement)
{
    bool targetClockIsEndClock = false;

    for (int i = 0; i < eventNumIncrement; ++i)
    {
        if (!next_DVD_event_to_FUTURE_exists())
        {
            targetClockIsEndClock = true;
            break;
        }

        handle_next_DVD_event_to_FUTURE();
    }

    double targetVoronoiClock = DBL_MAX;
    if (targetClockIsEndClock)
    {
        targetVoronoiClock = get_length_of_Voronoi_clock();
    }
    else
    {
        targetVoronoiClock = get_most_recent_event_of_DVD()->get_occurring_clock();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVD2DSimulator::progress_DVD_using_EVENTNUM_to_PAST(const int& eventNumIncrement)
{
    int absEventNumIncrement = abs(eventNumIncrement);

    const double mostRecentEventTime = get_most_recent_event_of_DVD()->get_occurring_clock();
    if ((get_current_Voronoi_clock() != mostRecentEventTime))
    {
        synchronize_disk_locations_and_current_clock_to_VD(mostRecentEventTime);
        --absEventNumIncrement;
    }


    bool targetClockIsZeroClock = false;

    for (int i = 0; i < absEventNumIncrement; ++i)
    {
        if (!next_DVD_event_to_PAST_exists())
        {
            targetClockIsZeroClock = true;
            break;
        }

        reversely_handle_next_DVD_event_to_PAST();
    }


    double targetVoronoiClock = DBL_MAX;
    if (targetClockIsZeroClock)
    {
        targetVoronoiClock = 0.0;
    }
    else
    {
        targetVoronoiClock = get_most_recent_event_of_DVD()->get_occurring_clock();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVD2DSimulator::handle_next_DVD_event_to_FUTURE()
{
    EventOfDynamicVD2D* nextFutureEvent = get_next_DVD_event_to_FUTURE();

	switch (nextFutureEvent->get_event_type())
	{
	    case EventOfDynamicVD2D::EDGE_FLIP: //flip this VEdge
	    {
		    EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(nextFutureEvent);
		    handle_EDGE_FLIP_event_to_FUTURE(edgeFlippingEvent);
		    break;
	    }

	    case EventOfDynamicVD2D::DISK_COLLISION: //move disks and change motion vectors 
	    {
		    DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(nextFutureEvent);
		    handle_DISK_COLLISION_event_to_FUTURE(diskCollisionEvent);
		    break;
	    }

	    case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
	    {
		    DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(nextFutureEvent);
		    handle_DISK_VELOCITY_CHANGE_event_to_FUTURE(diskVelocityChangeEvent);
		    break;
	    }

	    case EventOfDynamicVD2D::DISK_HOPPING:
	    {
		    DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(nextFutureEvent);
		    handle_DISK_HOPPING_event_to_FUTURE(diskHoppingEvent);
		    break;
	    }
	}

    set_current_Voronoi_clock(nextFutureEvent->get_occurring_clock());

    increment_index_of_most_recent_event();
}


void DynamicVD2DSimulator::reversely_handle_next_DVD_event_to_PAST()
{
    EventOfDynamicVD2D* nextPastEvent = get_next_DVD_event_to_PAST();

    switch (nextPastEvent->get_event_type())
	{
	    case EventOfDynamicVD2D::EDGE_FLIP: //flip this VEdge
	    {
		    EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(nextPastEvent);
		    handle_EDGE_FLIP_event_to_PAST(edgeFlippingEvent);
		    break;
	    }

	    case EventOfDynamicVD2D::DISK_COLLISION: //move disks and change motion vectors 
	    {
		    DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(nextPastEvent);
		    handle_DISK_COLLISION_event_to_PAST(diskCollisionEvent);
		    break;
	    }

	    case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
	    {
		    DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(nextPastEvent);
		    handle_DISK_VELOCITY_CHANGE_event_to_PAST(diskVelocityChangeEvent);
		    break;
	    }

	    case EventOfDynamicVD2D::DISK_HOPPING:
	    {
		    DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(nextPastEvent);
		    handle_DISK_HOPPING_event_to_PAST(diskHoppingEvent);
		    break;
	    }
	}

    set_current_Voronoi_clock(nextPastEvent->get_occurring_clock());

    decrement_index_of_most_recent_event();
}




void DynamicVD2DSimulator::get_dynamic_disks(list<DynamicDisk*>& dynamicDisks) const
{
    dynamicDisks = m_DynamicDisks;
}

void DynamicVD2DSimulator::get_containers(list<DynamicDisk*>& containers) const
{
    containers = m_AllContainers;
}

const DynamicDisk* DynamicVD2DSimulator::get_outer_most_container() const
{
    return m_AllContainers.front();
}

const DynamicVoronoiDiagramCIC& DynamicVD2DSimulator::get_dynamic_VD() const
{
    return m_DynamicVD;
}

const TimeStatistics<std::milli>& DynamicVD2DSimulator::get_DVDSIM_time_statistics() const
{
    return m_CompTimeStatistics;
}

string DynamicVD2DSimulator::get_disk_file_name_with_path() const
{
    return m_DiskFileNameWthPath;
}

double  DynamicVD2DSimulator::get_length_of_Voronoi_clock() const
{
    return m_DynamicVD.get_length_of_Voronoi_clock();
}

double DynamicVD2DSimulator::get_coefficient_of_restitution() const
{
    return m_DynamicVD.get_coefficient_of_restitution();
}

double DynamicVD2DSimulator::get_current_Voronoi_clock() const
{
    return m_DynamicVD.get_current_Voronoi_clock();
}

int DynamicVD2DSimulator::get_num_of_disks() const
{
    return m_DynamicDisks.size();
}

int DynamicVD2DSimulator::get_num_of_containers() const
{
    return m_AllContainers.size();
}

const EventOfDynamicVD2D* DynamicVD2DSimulator::get_most_recent_event_of_DVD() const
{
    return m_DynamicVD.get_most_recent_event();
}


int DynamicVD2DSimulator::get_total_num_of_events_of_DVD() const
{
    return m_DynamicVD.get_total_event_num();
}

const DiskCollisionWatcher& DynamicVD2DSimulator::get_disk_collision_watcher() const
{
    return m_DynamicVD.get_disk_collision_watcher();
}

void DynamicVD2DSimulator::set_current_Voronoi_clock(const double& currentClock)
{
    m_DynamicVD.set_current_Voronoi_clock(currentClock);
}


void DynamicVD2DSimulator::set_length_of_Voronoi_clock(const double& lengthOfVoronoiClock)
{
    m_DynamicVD.set_length_of_Voronoi_clock(lengthOfVoronoiClock);
}


DynamicDisk* DynamicVD2DSimulator::create_disk(const DynamicDisk& disk) const
{
    DynamicDisk* newDynamicDisk = new DynamicDisk(disk);
    return newDynamicDisk;
}

DynamicDisk* DynamicVD2DSimulator::create_dynamic_object(DynamicDisk* dynamicObject) const
{
    return create_disk(*dynamicObject);
}

void DynamicVD2DSimulator::delete_disk(DynamicDisk*& disk) const
{
    delete disk;
    disk = NULL;
}

bool DynamicVD2DSimulator::this_clock_is_future(const double& targetVoronoiClock) const
{
    return m_DynamicVD.this_clock_is_future(targetVoronoiClock);
}

bool DynamicVD2DSimulator::this_clock_is_past(const double& targetVoronoiClock) const
{
    return m_DynamicVD.this_clock_is_past(targetVoronoiClock);
}



void DynamicVD2DSimulator::write_hexa_characters(ofstream& fout, const int& size, unsigned char* chars) const
{
    for (int i = 0; i < size; ++i)
    {
        fout << hex << setw(2) << setfill('0') << int(chars[i]);
    }
}


void DynamicVD2DSimulator::get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(const double& value, unsigned char* chars) const
{
    memcpy(chars, &value, sizeof(double));
}


double DynamicVD2DSimulator::get_double_type_value_from_hexa_characters_in_order_of_little_endian(char* charsOfValue) const
{
    double eventTime = 0.0;

    const int sizeOfDouble = sizeof(double);

    unsigned char charsOfValueInMemory[sizeOfDouble];

    for (int i = 0; i < sizeOfDouble; ++i)
    {
        unsigned int char1 = char_to_val(charsOfValue[2 * i]) * 16;
        unsigned int char2 = char_to_val(charsOfValue[2 * i + 1]);

        unsigned char newChar = char1 + char2;
        charsOfValueInMemory[i] = newChar;
    }

    memcpy(&eventTime, charsOfValueInMemory, sizeof(double));

    return eventTime;
}


unsigned int DynamicVD2DSimulator::char_to_val(char aCharacter) const
{
    unsigned int value = aCharacter;

    if (value >= 48 && value < 58)
    {
        value = value - 48;
    }
    else if (value >= 65 && value < 91)
    {
        value = value - 55;
    }
    else if (value >= 97 && value < 123)
    {
        value = value - 87;
    }

    return value;
}

void DynamicVD2DSimulator::increment_index_of_most_recent_event()
{
    m_DynamicVD.increment_index_of_most_recent_event();
}

void DynamicVD2DSimulator::decrement_index_of_most_recent_event()
{
    m_DynamicVD.decrement_index_of_most_recent_event();
}

void DynamicVD2DSimulator::handle_EDGE_FLIP_event_to_FUTURE(EdgeFlipEvent* edgeFlipEvent)
{
    m_DynamicVD.update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_FUTURE(edgeFlipEvent);
}

void DynamicVD2DSimulator::handle_DISK_COLLISION_event_to_FUTURE(DiskCollisionEvent* diskCollisionEvent)
{
    m_DynamicVD.relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_FUTURE(diskCollisionEvent);
}

void DynamicVD2DSimulator::handle_DISK_VELOCITY_CHANGE_event_to_FUTURE(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    m_DynamicVD.relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_FUTURE(diskVelocityChangeEvent);
}

void DynamicVD2DSimulator::handle_DISK_HOPPING_event_to_FUTURE(DiskHoppingEvent* diskHoppingEvent)
{
    m_DynamicVD.relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_FUTURE(diskHoppingEvent);
}

void DynamicVD2DSimulator::handle_EDGE_FLIP_event_to_PAST(EdgeFlipEvent* edgeFlipEvent)
{
    m_DynamicVD.update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_PAST(edgeFlipEvent);
}

void DynamicVD2DSimulator::handle_DISK_COLLISION_event_to_PAST(DiskCollisionEvent* diskCollisionEvent)
{
    m_DynamicVD.relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_PAST(diskCollisionEvent);
}

void DynamicVD2DSimulator::handle_DISK_VELOCITY_CHANGE_event_to_PAST(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    m_DynamicVD.relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_PAST(diskVelocityChangeEvent);
}

void DynamicVD2DSimulator::handle_DISK_HOPPING_event_to_PAST(DiskHoppingEvent* diskHoppingEvent)
{
    m_DynamicVD.relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_PAST(diskHoppingEvent);
}

EventOfDynamicVD2D* DynamicVD2DSimulator::get_next_DVD_event_to_FUTURE() const
{
    return m_DynamicVD.get_next_DVD_event_to_FUTURE();
}

EventOfDynamicVD2D* DynamicVD2DSimulator::get_next_DVD_event_to_PAST() const
{
    return m_DynamicVD.get_next_DVD_event_to_PAST();
}

bool DynamicVD2DSimulator::next_DVD_event_to_FUTURE_exists() const
{
    return m_DynamicVD.next_DVD_event_to_FUTURE_exists();
}

bool DynamicVD2DSimulator::next_DVD_event_to_PAST_exists() const
{
    return m_DynamicVD.next_DVD_event_to_PAST_exists();
}

bool DynamicVD2DSimulator::next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(const double& targetVoronoiClock) const
{
    return m_DynamicVD.next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(targetVoronoiClock);
}

bool DynamicVD2DSimulator::next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(const double& targetVoronoiClock) const
{
    return m_DynamicVD.next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(targetVoronoiClock);
}

bool DynamicVD2DSimulator::next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock() const
{
    return m_DynamicVD.next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock();
}

void DynamicVD2DSimulator::___start_to_track_statistics_in_playing_simulation(const double& targetVoronoiClock)
{
    m_DynamicVD.___start_to_track_statistics_of_playing_DVD(targetVoronoiClock);
}

void DynamicVD2DSimulator::___stop_to_track_statistics_in_playing_simulation()
{
    m_DynamicVD.___stop_to_track_statistics_of_playing_DVD();
}



bool DynamicVD2DSimulator::compare_collision_clock(const pair<int, double>& collidedDiskNClock1, const pair<int, double>& collidedDiskNClock2)
{
    if (collidedDiskNClock1.second < collidedDiskNClock2.second)
    {
        return true;
    }
    else
    {
        return false;
    }
}

