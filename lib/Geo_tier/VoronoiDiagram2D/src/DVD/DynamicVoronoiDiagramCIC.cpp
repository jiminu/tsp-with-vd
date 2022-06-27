#include "DynamicVoronoiDiagramCIC.h"
#include "VoronoiDiagram2DC.h"
#include "BucketForDisk.h"
#include "EdgeFlipTimeCalculator.h"
#include "DiskCollisionTimeCalculator.h"
using namespace V::GeometryTier;

#include "rg_Point2D.h"

#include <set>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <map>
using namespace std;

DynamicVoronoiDiagramCIC::DynamicVoronoiDiagramCIC()
{
    initialize_parameters();
	prepare_statistics_in_generating_EVENTSEQ();
}


DynamicVoronoiDiagramCIC::DynamicVoronoiDiagramCIC(const DynamicVoronoiDiagramCIC& dynamicVD)
{
	copy_from(dynamicVD);
}


DynamicVoronoiDiagramCIC::DynamicVoronoiDiagramCIC(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks)
{
    initialize_parameters();
	prepare_statistics_in_generating_EVENTSEQ();

    construct_initial_VD(container, dynamicDisks);
}


DynamicVoronoiDiagramCIC::DynamicVoronoiDiagramCIC(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks)
{
    initialize_parameters();
	prepare_statistics_in_generating_EVENTSEQ();

    construct_initial_VD(container, dynamicDisks);
}


DynamicVoronoiDiagramCIC::~DynamicVoronoiDiagramCIC()
{
    clear_EVENTSEQ();
}



DynamicVoronoiDiagramCIC& DynamicVoronoiDiagramCIC::operator=(const DynamicVoronoiDiagramCIC& dynamicVD)
{
	if (this != &dynamicVD)
	{
        clear();
        copy_from(dynamicVD);
	}

	return *this;
}



void DynamicVoronoiDiagramCIC::clear()
{
    clear_input();

    clear_EVENTSEQ();

    clear_mapper_from_id_to_generator();

    clear_VD_reinstatement_objects();

    clear_event_generation_related_objects();

    clear_time_statistics();
    clear_count_statistics();

    initialize_parameters();
}


void DynamicVoronoiDiagramCIC::clear_EVENTSEQ_related_objects()
{
    clear_EVENTSEQ();
    clear_event_generation_related_objects();

    clear_time_statistics();
    clear_count_statistics();
}


void DynamicVoronoiDiagramCIC::clear_EVENTSEQ()
{
    if (!m_EventSequence.empty())
    {
        for (vector<EventOfDynamicVD2D*>::const_iterator it_Event = m_EventSequence.begin();
            it_Event != m_EventSequence.end();
            it_Event++)
        {
            delete (*it_Event);
        }

        m_EventSequence.clear();
    }

}

void DynamicVoronoiDiagramCIC::clear_time_statistics()
{
    for (int i = 0; i < DVD_TIME_STATISTICS_SIZE; ++i)
    {
        m_TimeStatistics[i].reset(0);
    }
}


void DynamicVoronoiDiagramCIC::clear_count_statistics()
{
    for (int i = 0; i < DVD_COUNT_STATISTICS_SIZE; ++i)
    {
        m_CountStatistics[i].reset(0);
    }
}


void DynamicVoronoiDiagramCIC::clear_event_generation_related_objects()
{
    m_PriorityQForEdgeFlipping.clear();
    m_PriorityQForDiskCollision.clear();
    m_PriorityQForDiskVelocityChange.clear();

    m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.clear();
    m_VEdgesDefinedByIdentical4DiskWithFlippingVEdge.clear();
    m_NewVEdgesByDiskHopping.clear();
    m_DeletedVEdgesByDiskHopping.clear();
    m_RedefinedVEdgesByDiskHopping.clear();

    m_InThisContainerDiskWillBeInserted.clear();
}


void DynamicVoronoiDiagramCIC::clear_VD_reinstatement_objects()
{
    m_InitialVDForReinstatement.clear();
    m_MapFromCopiedInitialStateGeneratorsToOriginalGenerators.clear();
    m_MapFromOriginalDisksToCopiedInitialStateDisks.clear();
}


void DynamicVoronoiDiagramCIC::clear_input()
{    
    m_DiskHoppingProbabilityMatrix.clear();
    m_LengthOfVoronoiClock  = 0.0;
    m_CollisionHandlingType = DiskCollisionHandler();
    m_CollisionWatcher      = DiskCollisionWatcher();
    m_NumOfContainers       = 0;
}


void DynamicVoronoiDiagramCIC::copy_from(const DynamicVoronoiDiagramCIC& fromThisDVD)
{
    /************************************************************************/
    /*  Objects using event generation like prioirity queue are not copied. */
    /*  Dynamic disk in user data of Generator2D are shallowly copied.      */        
    /************************************************************************/

    VoronoiDiagramCIC::copyFrom(fromThisDVD);

    copy_input(fromThisDVD);
    copy_DVD_current_state_variables(fromThisDVD);
    copy_time_statistics(fromThisDVD);
    copy_count_statistics(fromThisDVD);

    create_event_sequence(fromThisDVD);

    // ** Generator ����, VEdge ������ This VD�� Other VD�� �������� ���� 
    map<Generator2D*, Generator2D*> generatorLinkFromOtherVDToThisVD;
    map<DynamicDisk*, DynamicDisk*> diskLinkFromOtherVDToThisVD;

    //0) create linke btw generators and btw moving disks
    make_link_btw_generators_from_other_VD_to_this_VD(fromThisDVD, generatorLinkFromOtherVDToThisVD);
    make_link_btw_disks_and_containers_from_other_VD_to_this_VD(fromThisDVD, diskLinkFromOtherVDToThisVD);

    modify_generators_in_event_sequence(generatorLinkFromOtherVDToThisVD);

    link_dynamic_disk_id_to_its_generator_pointer();

    copy_initial_VD_and_containers_and_disks_for_reinstatement();
}


void DynamicVoronoiDiagramCIC::copy_input(const DynamicVoronoiDiagramCIC& fromThisDVD)
{
    m_LengthOfVoronoiClock          = fromThisDVD.m_LengthOfVoronoiClock;
    m_CollisionHandlingType         = fromThisDVD.m_CollisionHandlingType;
    m_CollisionWatcher              = fromThisDVD.m_CollisionWatcher;
    m_DiskHoppingProbabilityMatrix  = fromThisDVD.m_DiskHoppingProbabilityMatrix;
    m_NumOfContainers               = fromThisDVD.m_NumOfContainers;
}


void DynamicVoronoiDiagramCIC::copy_DVD_current_state_variables(const DynamicVoronoiDiagramCIC& fromThisDVD)
{
    m_CurrentVoronoiClock                             = fromThisDVD.m_CurrentVoronoiClock;
    m_IndexOfMostRecentEvent                          = fromThisDVD.m_IndexOfMostRecentEvent;

    m_DuringEVENTSEQGeneration                        = fromThisDVD.m_DuringEVENTSEQGeneration;
    m_Vers1p0EVENTSEQFileUsed                         = fromThisDVD.m_Vers1p0EVENTSEQFileUsed; 
    m_InitialPQIsConstructed                          = fromThisDVD.m_InitialPQIsConstructed;
    m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ = fromThisDVD.m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ;
}


void DynamicVoronoiDiagramCIC::create_event_sequence(const DynamicVoronoiDiagramCIC& fromThisDVD)
{
    const vector<EventOfDynamicVD2D*>& eventSequenceInFromDVD = fromThisDVD.get_EVENTSEQ();

    const int eventSequenceSize = eventSequenceInFromDVD.size();
  
    m_EventSequence.resize(eventSequenceSize);

    for (int i = 0; i < eventSequenceInFromDVD.size(); ++i)
    {
        EventOfDynamicVD2D* currEventInFromDVD = eventSequenceInFromDVD[i];
        EventOfDynamicVD2D* newlyGeneratedEvent  = NULL;

        switch (currEventInFromDVD->get_event_type())
        {
            case EventOfDynamicVD2D::EDGE_FLIP:
            {
                EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(currEventInFromDVD);
                newlyGeneratedEvent = new EdgeFlipEvent(*edgeFlippingEvent);
                break;
            }

            case EventOfDynamicVD2D::DISK_COLLISION:
            {
                DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(currEventInFromDVD);
                newlyGeneratedEvent = new DiskCollisionEvent(*diskCollisionEvent);
                break;
            }

            case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
            {
                DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(currEventInFromDVD);
                newlyGeneratedEvent = new DiskVelocityChangeEvent(*diskVelocityChangeEvent);
                break;
            }

            case EventOfDynamicVD2D::DISK_HOPPING:
            {
                DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(currEventInFromDVD);
                newlyGeneratedEvent = new DiskHoppingEvent(*diskHoppingEvent);
                break;
            }
        }

        m_EventSequence[i] = newlyGeneratedEvent;
    }
}



void DynamicVoronoiDiagramCIC::copy_time_statistics(const DynamicVoronoiDiagramCIC& fromThisDVD)
{
    for (int i = 0; i < DVD_TIME_STATISTICS_SIZE; ++i)
    {
        m_TimeStatistics[i] = fromThisDVD.m_TimeStatistics[i];
    }
}


void DynamicVoronoiDiagramCIC::copy_count_statistics(const DynamicVoronoiDiagramCIC& fromThisDVD)
{
    for (int i = 0; i < DVD_COUNT_STATISTICS_SIZE; ++i)
    {
        m_CountStatistics[i] = fromThisDVD.m_CountStatistics[i];
    }
}

Generator2D* DynamicVoronoiDiagramCIC::get_generator_by_its_input_disk_id(const int& inputDiskId) const
{
    // DVD doesn't have this generator id
    if (m_LinkFromInputIDToGenerator.find(inputDiskId) == m_LinkFromInputIDToGenerator.end())
    {
        return NULL;
    }
    else
    {
        return m_LinkFromInputIDToGenerator.at(inputDiskId);
    }
}


void DynamicVoronoiDiagramCIC::link_dynamic_disk_id_to_its_generator_pointer()
{
    list<Generator2D*> allGeneratorsIncludingEnclosingCircleGenerator;
    getAllGenerators_hierarchy(allGeneratorsIncludingEnclosingCircleGenerator);

    for (list<Generator2D*>::const_iterator it_Generator = allGeneratorsIncludingEnclosingCircleGenerator.begin();
        it_Generator != allGeneratorsIncludingEnclosingCircleGenerator.end();
        it_Generator++)
    {
        Generator2D* currGenerator = *it_Generator;
        DynamicDisk* dynamicDisk = static_cast<DynamicDisk*>(currGenerator->getUserData());
        m_LinkFromInputIDToGenerator.insert(make_pair(dynamicDisk->getID(), currGenerator));
    }
}


VEdge2D* DynamicVoronoiDiagramCIC::get_VEdge_defined_by_these_four_input_disks(
    int* fourDiskIds, const bool& VEdgeDefinedByQuadrupletIsAlreadyFlipped /* = false*/) const
{
    EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = get_generator_quadruplet_by_its_id(fourDiskIds);
    VEdge2D* outputVEdge = get_VEdge_with_this_generator_quadruplet(generatorPtrQuadruplet, VEdgeDefinedByQuadrupletIsAlreadyFlipped);

    return outputVEdge;
}


void DynamicVoronoiDiagramCIC::initialize_parameters()
{
    m_DuringEVENTSEQGeneration                         = false;
    m_Vers1p0EVENTSEQFileUsed                          = false;
    m_InitialPQIsConstructed                           = false;
    m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ  = false;
    m_CurrentVoronoiClock                              = 0.0;
    m_LengthOfVoronoiClock                             = 0.0;
    m_IndexOfMostRecentEvent                           = -1;
    m_LastReservedNumber                               = 0;

    m_CollisionHandlingType.set_collision_type(DiskCollisionHandler::PHYSICAL_COLLISION);
    m_CollisionHandlingType.set_coefficient_of_restitution(1.0);
}


void DynamicVoronoiDiagramCIC::generate_EVENTSEQ(const double& lengthOfVoronoiClock, 
                                                 const bool& reinstateInitialVD /*= true*/)      
{
    m_DuringEVENTSEQGeneration = true;

    set_length_of_Voronoi_clock(lengthOfVoronoiClock);

    construct_initial_priority_queues();

    propagate_DVD_to_the_last_moment_of_Voronoi_clock();

    if (reinstateInitialVD)
    {
        clear_priority_queues_for_collision_and_flip();

        reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
    }

    complete_statistics_in_generating_EVENTSEQ(); 

    m_DuringEVENTSEQGeneration = false;
}



void DynamicVoronoiDiagramCIC::generate_EVENTSEQ(const int& numOfCollisions, 
                                                 const bool& reinstateInitialVD /*= true*/)
{
    m_DuringEVENTSEQGeneration = true;

    set_length_of_Voronoi_clock (0.0);

    construct_initial_priority_queues();

    propagate_DVD_until_target_collision_numbers_btw_disks(numOfCollisions);

	if (reinstateInitialVD)
    {
        clear_priority_queues_for_collision_and_flip();

		reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
	}

    
    complete_statistics_in_generating_EVENTSEQ();

    m_DuringEVENTSEQGeneration = false;
}


void DynamicVoronoiDiagramCIC::generate_EVENTSEQ(
    const DiskCollisionWatcher& collisionWatcher, 
    const bool& reinstateInitialVD /*= true*/)
{
    m_DuringEVENTSEQGeneration = true;

    DiskCollisionHandler collisionHandler;
    collisionHandler.set_collision_type(DiskCollisionHandler::PHYSICAL_COLLISION);
    collisionHandler.set_coefficient_of_restitution(1.0);

    set_parameters(0.0, collisionWatcher, collisionHandler);

    construct_initial_priority_queues();

    propagate_DVD_until_target_collision_numbers_of_watch_objects();

    if (reinstateInitialVD)
    {
        clear_priority_queues_for_collision_and_flip();

        reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
    }

    complete_statistics_in_generating_EVENTSEQ();

    m_DuringEVENTSEQGeneration = false;
}



void DynamicVoronoiDiagramCIC::generate_EVENTSEQ(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType, const bool& reinstateInitialVD /*= true*/)
{
    m_DuringEVENTSEQGeneration = true;

    set_parameters(lengthOfVoronoiClock, collisionHandlingType);

    construct_initial_priority_queues();

    propagate_DVD_to_the_last_moment_of_Voronoi_clock();

    if (reinstateInitialVD)
    {
        clear_priority_queues_for_collision_and_flip();

        reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
    }

    complete_statistics_in_generating_EVENTSEQ();

    m_DuringEVENTSEQGeneration = false;
}


DiskCollisionEvent* DynamicVoronoiDiagramCIC::generate_EVENTSEQ_by_next_collision(const double& endVoronoiClock, const bool& reinstateInitialVD /*= true*/)
{
    m_DuringEVENTSEQGeneration = true;

    if (!m_InitialPQIsConstructed)
    {
        construct_initial_priority_queues();
    }

    DiskCollisionEvent* output 
        = propagate_DVD_to_the_next_collision(endVoronoiClock);

    if (reinstateInitialVD)
    {
        clear_priority_queues_for_collision_and_flip();

        reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
    }

    complete_statistics_in_generating_EVENTSEQ();

    m_DuringEVENTSEQGeneration = false;

    return output;
}


void DynamicVoronoiDiagramCIC::generate_EVENTSEQ_collision_avoidance(const double& lengthOfVoronoiClock, const double& clockRelativeToCollisionClock, const bool& reinstateInitialVD /*= true*/)
{
    m_DuringEVENTSEQGeneration = true;

    initialize_statistics_for_collision_avoidance();

    ___start_clock_in_collision_avoidance(DVD_COL_AVD_TOTAL);

	set_length_of_Voronoi_clock(lengthOfVoronoiClock);

    while (1)
    {
		construct_initial_priority_queues();

        tuple<bool, double, DynamicDisk*, DynamicDisk*> lastOccuredDiskCollsionEvent =  propagate_DVD_until_next_collision();
        bool diskCollisionEventExist = get<0>(lastOccuredDiskCollsionEvent);
       
        if (diskCollisionEventExist)
        {
            double diskCollisionClock  = get<1>(lastOccuredDiskCollsionEvent);
            double targetVoronoiClock  = get_next_target_time_to_go_back(clockRelativeToCollisionClock, diskCollisionClock);

			go_back_to_target_Voronoi_clock_and_remove_events_after_the_target_Voronoi_clock(targetVoronoiClock);

            DynamicDisk* disk1 = get<2>(lastOccuredDiskCollsionEvent);
			DynamicDisk* disk2 = get<3>(lastOccuredDiskCollsionEvent);

            tuple<DynamicDisk*, rg_Point2D, DynamicDisk*, rg_Point2D> outputVelocityVectors = compute_velocity_vectors_by_brute_force_search_for_collision_avoidance(disk1, disk2, targetVoronoiClock);

			DiskVelocityChangeEvent* disk1VelocityChangeEvent  = create_velocity_vector_change_event(targetVoronoiClock, get<0>(outputVelocityVectors), get<1>(outputVelocityVectors));
			DiskVelocityChangeEvent* disk2VelocityChangeEvent  = create_velocity_vector_change_event(targetVoronoiClock, get<2>(outputVelocityVectors), get<3>(outputVelocityVectors));
		   
            change_to_new_velocity_vector_and_store_old_one(disk1VelocityChangeEvent);
			change_to_new_velocity_vector_and_store_old_one(disk2VelocityChangeEvent);

            store_this_VELOCITY_VECTOR_CHANGE_event_in_the_event_history_vector(disk1VelocityChangeEvent);
			store_this_VELOCITY_VECTOR_CHANGE_event_in_the_event_history_vector(disk2VelocityChangeEvent);

			clear_priority_queues_for_collision_and_flip();
        }
        else
        {
            clear_priority_queues_for_collision_and_flip();
            break;
        }
    }

	if (reinstateInitialVD)
	{
		reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
	}

	clear_priority_queues_for_collision_and_flip();


	___end_clock_in_collision_avoidance(DVD_COL_AVD_TOTAL);

    finalize_statistics_for_collision_avoidance();

    m_DuringEVENTSEQGeneration = false;
}



list<Generator2D*> DynamicVoronoiDiagramCIC::find_neighbor_generators(
    Generator2D* anchorGenerator, const double& startVoronoiClock, const double& endVoronoiClock, const bool& thisIsContainer /*=false*/)
{
    if (get_current_Voronoi_clock() != startVoronoiClock)
    {
        construct_VD_at_target_Voronoi_clock(startVoronoiClock);
    }

    list<Generator2D*> neighbors;
    anchorGenerator->getNeighborGenerators(neighbors, thisIsContainer);

    double currentVoronoiClock = get_current_Voronoi_clock();

    while (next_DVD_event_to_FUTURE_exists())
    {
        go_to_next_event_Voronoi_clock();

        currentVoronoiClock = get_current_Voronoi_clock();

        if (currentVoronoiClock > endVoronoiClock)
        {
            break;
        }

        const EventOfDynamicVD2D* event = get_most_recent_event();

        if (event->get_event_type() == EventOfDynamicVD2D::EDGE_FLIP)
        {
            EdgeFlipEvent* edgeFlipEvent= dynamic_cast<EdgeFlipEvent*>(const_cast<EventOfDynamicVD2D*>(event));

            EdgeFlipEvent::GeneratorPtrQuadruplet generatorQuadruplet = edgeFlipEvent->get_generator_quadruplet_before_flip();
            Generator2D* headGenerator = generatorQuadruplet.at(2);
            Generator2D* tailGenrator = generatorQuadruplet.at(3);
            
            if (headGenerator == anchorGenerator)
            {
                neighbors.push_back(tailGenrator);
            }
            else if (tailGenrator == anchorGenerator)
            {
                 neighbors.push_back(headGenerator);
            }
        }
    }

    if (get_current_Voronoi_clock() != startVoronoiClock)
    {
        construct_VD_at_target_Voronoi_clock(startVoronoiClock);
    }


    neighbors.sort();
    neighbors.unique();

    return neighbors;
}



void DynamicVoronoiDiagramCIC::go_to_target_Voronoi_clock(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed /*= false*/)
{
    double targetClock = targetVoronoiClock;

    if (targetClock == get_current_Voronoi_clock())
    {
        return;
    }

    if (targetClock > get_length_of_Voronoi_clock())
    {
        targetClock = get_length_of_Voronoi_clock();
    }
    else if (targetClock < 0.0)
    {
        targetClock = 0.0;
    }


    ___start_to_track_statistics_of_playing_DVD(targetClock);

    construct_VD_at_target_Voronoi_clock(targetClock, VEdgeGeometryIsComputed);

    ___stop_to_track_statistics_of_playing_DVD();
}



void DynamicVoronoiDiagramCIC::increase_Voronoi_clock(const double& VoronoiClockIncrement /*= 1.0*/, const bool& VEdgeGeometryIsComputed /*= false*/)
{
    double targetVoronoiClock = get_current_Voronoi_clock() + VoronoiClockIncrement;

    go_to_target_Voronoi_clock(targetVoronoiClock, VEdgeGeometryIsComputed);
}



void DynamicVoronoiDiagramCIC::go_to_next_event_Voronoi_clock(const int& eventNumIncrement /*= 1*/, const bool& VEdgeGeometryIsComputed /*= false*/)
{
    if (eventNumIncrement == 0)
    {
        return;
    }
    else if (m_IndexOfMostRecentEvent == -1 && eventNumIncrement < 0)
    {
        return;
    }
    else if (m_IndexOfMostRecentEvent == m_EventSequence.size() && eventNumIncrement > 0)
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
        compute_all_VEdges_equation_of_VD();
    }
}



void DynamicVoronoiDiagramCIC::write_statistics_for_generating_EVENTSEQ(const string& fileNameWithPath, const bool& appendingMode/*=false*/) const
{
    ofstream fout;
    
    if (appendingMode)
    {
        fout.open(fileNameWithPath, ofstream::out | ofstream::app);
    }
    else
    {
        fout.open(fileNameWithPath);
    }

    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(3);

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    DynamicDisk* containerDisk = get_most_outer_container();

    fout << "[ INPUT ]" << endl << endl;

    fout << "*. Num of disks                  :" << "\t" << dynamicDisks.size()           << endl;
    fout << "*. Simulation duration           :" << "\t" << get_length_of_Voronoi_clock()  << endl;
    fout << "*. Coefficient of restitution    :" << "\t" << m_CollisionHandlingType.get_coefficient_of_restitution()  << endl;
                                              
    fout << "*. Container center              :" << "\t" << "( " << containerDisk->getCircle().getCenterPt().getX() << ", "
                                                                 << containerDisk->getCircle().getCenterPt().getY() << " ), " 
         << "*. Radius : "                                       << containerDisk->getCircle().getRadius()          << endl;         
    fout << endl << endl;


    fout << "[ PHASE 1 - OUTPUT ]" << endl << endl;
    fout << "*. TIME UNIT:"            << "\t" << "Second" << endl;
    fout << "*. Computation time   :"  << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PHASE_ONE_TOTAL)  << endl;
    fout << "*. Number of events   :"   << "\t" << m_EventSequence.size()                                                                   << endl;
    fout << "   - Flipping   events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_EDGE_FLIPPING_IN_EVENT_HISTORY)          << endl;
    fout << "   - Collision  events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_COLLISION_IN_EVENT_HISTORY)         << endl;
    fout << "        - Btw disks (not with container) :" << "\t" <<  m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_COLLISION_BTW_INPUT_DISKS_NOT_WITH_CONTAINER_IN_EVENT_HISTORY) << endl;
    fout << "   - Velocity C events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_VELOCITY_CHANGE_IN_EVENT_HISTORY)   << endl;
    fout << "   - Del/Insert events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_HOPPING_IN_EVENT_HISTORY)   << endl;
    fout << "   - (R) Del/Insert events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_REAL_OCCURING_DELETE_N_INSERT_DISK_IN_EVENT_HISTORY) << endl;

    fout << endl << endl;


    fout << "[ WHOLE PROCESS TIME ]" << endl << endl;

    fout << "*. TOTAL                                                                        :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PHASE_ONE_TOTAL)                                                 << endl;
    fout << "1. Construct initial CIC VD                                                     :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_CONSTRUCT_INITIAL_CIC_VD)                                        << endl;
    fout << "2. Check whether disks are intersected by VD                                    :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_CHECK_WHETHER_DISKS_ARE_INTERSECTED)                             << endl;
    fout << "3. Construct initial priority Qs                                                :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS) << endl;
    fout << "   3.1. Construct initial flipping Q                                            :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME)             << endl;
    fout << "   3.2. Construct initial collision Q                                           :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME)             << endl;
    fout << "4. Generate events during simulation time window                                :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK)               << endl;
    fout << "5. Make VD and disks be at time zero                                            :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_REINSTATE_THE_INITIAL_VORONOI_DIAGRAM)                                 << endl;
    fout << endl << endl;

    fout << "[ TIME OF INTEREST ]" << endl << endl;

    fout << "1. Find possible flipping time of VEdges :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_FIND_POSSIBLE_FLIPPING_TIME_OF_EDGE)      << endl;
    fout << "   1.1. Get 4 defining disks             :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_GET_4DISKS_DEFINING_VEDGE)                 << endl;
    fout << "   1.2. Sort 4 disks in non-increasing   :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_SORTING_GENERATORS_IN_NON_INCRESING_ORDER) << endl;
    fout << "   1.3. Make polynomials                 :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_MAKING_EDGE_FLIPPING_POLYNOMIAL_EQ)       << endl;
    fout << "   1.4. Solve polynomials                :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_SOLVING_EDGE_FLIPPING_POLYNOMIAL_EQ)       << endl;
    fout << "2. Find possible collision time of disks :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_FIND_POSSIBLE_COLLISION_TIME_OF_DISKS)     << endl;
    fout << "3. Do Priority Q bookkeeping             :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PRIORITYS_BOOKKEEPING)                    << endl;
    fout << "   2.1. Do flipping Q bookkeeping        :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING)    << endl;
    fout << "   2.2. Do collision Q boookeeping       :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION)   << endl;
	fout << "   2.3. Do velocity   Q bookkeeping      :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_VELOCITY_CHANGE)<< endl;
    fout << "4. Move disks                            :" << "\t"<< m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_DISK_MOVING)                              << endl;
    fout << "5. Delete/insert disks                   :" << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_DISK_HOPPING) << endl;
    fout << endl << endl;

    fout << "[ COUNTS OF EVENT COMPUTATION ]" << endl << endl;

    fout << "1. Edge flipping  :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_EDGE_FLIPPING_IN_COMPUTATION)  << endl;
    fout << "2. Disk collision :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_COLLISION_IN_COMPUTATION) << endl;
    fout << "3. Velocity C     :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_VELOCITY_CHANGE_IN_COMPUTATION) << endl;
    fout << "4. Del/Insert     :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DELETE_N_INSERT_DISK_IN_COMPUTATION) << endl;

    fout << endl << endl;



    fout << "[ SPECIAL CASE OF IDENTICAL FLIPPING TIME ]" << endl << endl;

    fout << "1. 2-SAME :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_2_IDENTICAL_FLIPPING_TIME) << endl;
    fout << "2. 3-SAME :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_3_IDENTICAL_FLIPPING_TIME) << endl;
    fout << "3. 4-SAME :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_4_IDENTICAL_FLIPPING_TIME) << endl;
    fout << "4. 5-SAME :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_5_IDENTICAL_FLIPPING_TIME) << endl;

    fout << endl << endl;



    fout << "[ FLIPPING TYPE ]" << endl << endl;

    fout << "1. ORDINARY :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_ORDINARY_FLIPPING) << endl;
    fout << "2. 3_TO_2   :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_3_TO_2_FLIPPING)   << endl;
    fout << "3. 24_TO_33 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_24_TO_33_FLIPPING) << endl;
    fout << "4. 33_TO_24 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_33_TO_24_FLIPPING) << endl;
    fout << "5. 33_TO_22 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_33_TO_22_FLIPPING) << endl;
    fout << "6. SHADOW   :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_SHADOW_FLIPPING)   << endl;

    fout << endl << endl;



    fout << "[ EDGE FLIPPING POLYNOMIAL ]" << endl << endl;
    
    fout << "[ 1. DEGREE COUNT ]" << endl << endl;

    fout << "1. Degree 0 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE0_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "2. Degree 1 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE1_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "3. Degree 2 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE2_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "4. Degree 3 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE3_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "5. Degree 4 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE4_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "6. Degree 5 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE5_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "7. Degree 6 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE6_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "8. Degree 7 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE7_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;
    fout << "9. Degree 8 :" << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DEGREE8_POLYNOMIAL_FOR_EDGE_FLIPPING) << endl;

    fout << endl << endl;

    fout.close();
    
    return;
}



void DynamicVoronoiDiagramCIC::write_statistics_in_playing_DVD(const string& fileNameWithPath) const
{
    ofstream fout(fileNameWithPath.c_str());

    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(3);

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    DynamicDisk* containerDisk = get_most_outer_container();

    fout << "[ INPUT ]" << endl << endl;

    fout << "*. Num of disks                  :" << "\t" << dynamicDisks.size()          << endl;
    fout << "*. Simulation duration           :" << "\t" << get_length_of_Voronoi_clock()<< endl;
    fout << "*. Coefficient of restitution    :" << "\t" << m_CollisionHandlingType.get_coefficient_of_restitution()  << endl;
                                              
    fout << "*. Container center              :" << "\t" << "( " << containerDisk->getCircle().getCenterPt().getX() << ", "
                                                                 << containerDisk->getCircle().getCenterPt().getY() << " ), " 
         << "*. Radinus : "                                      << containerDisk->getCircle().getRadius()          << endl;         
    fout << endl << endl;

    fout << "[ PHASE 1 - OUTPUT ]" << endl << endl;
    fout << "*. TIME UNIT:"            << "\t" << "Second" << endl;
    fout << "*. Computation time   :"  << "\t" << m_TimeStatistics[PHASE1_COMP_TIME].timeInSec(DVD_TIME_PHASE_ONE_TOTAL) << endl;
    fout << "*. Number of events   :"   << "\t" << m_EventSequence.size()                                                                            << endl;
    fout << "   - Flipping   events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_EDGE_FLIPPING_IN_EVENT_HISTORY)          << endl;
    fout << "   - Collision  events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_COLLISION_IN_EVENT_HISTORY)         << endl;
    fout << "        - Btw disks (not with container) :" << "\t" <<  m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_COLLISION_BTW_INPUT_DISKS_NOT_WITH_CONTAINER_IN_EVENT_HISTORY) << endl;
    fout << "   - Velocity C events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_VELOCITY_CHANGE_IN_EVENT_HISTORY)   << endl;
    fout << "   - Del/Insert events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_HOPPING_IN_EVENT_HISTORY)   << endl;
    fout << "   - (R) Del/Insert events :"  << "\t" << m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_REAL_OCCURING_DELETE_N_INSERT_DISK_IN_EVENT_HISTORY) << endl;

    fout << endl << endl;


    fout << "[ PHASE 2 - OUTPUT  ]" << endl << endl;

    double cumulativeTime         = 0.0;
    double currTime               = 0.0;

    fout << "Number"                         << "\t";
    fout << "Time increment"                 << "\t";
    fout << "Curr time after time increment" << "\t";
    fout << "# Total Event"                  << "\t";
    fout << "# Flipping"                     << "\t";
    fout << "# Collision"                    << "\t";
    fout << "# Collision btw disks"          << "\t";
    fout << "Computation Time"               << "\t";
    fout << "Cumulative Comp. Time"          << "\t";
    fout << "Disk Moving Time"               << "\t";
    fout << "VD Topology Update Time"        << "\t";
    fout << "VD Geometry Update Time"        << endl;

    if (m_TimeStatistics[PHASE2_MOVING_DISKS_TIME].size() == 0)
    {
        for (int i = 0; i < m_TimeStatistics[PHASE2_COMP_TIME].size(); i++)
        {
            cumulativeTime += m_TimeStatistics[PHASE2_COMP_TIME].time(i);
            currTime       += m_TimeStatistics[PHASE2_SIMULATION_TIME_INCREMNT].time(i);

            fout << i + 1 << "\t";
            fout << m_TimeStatistics[PHASE2_SIMULATION_TIME_INCREMNT].time(i) << "\t";
            fout << currTime << "\t";
            fout << m_CountStatistics[PHASE2_FLIPPING_COUNT].count(i) + m_CountStatistics[PHASE2_COLLISION_COUNT].count(i) << "\t";
            fout << m_CountStatistics[PHASE2_FLIPPING_COUNT].count(i) << "\t";
            fout << m_CountStatistics[PHASE2_COLLISION_COUNT].count(i) << "\t";
            fout << m_CountStatistics[PHASE2_COLLISION_BTW_DISKS_COUNT].count(i) << "\t";
            fout << m_TimeStatistics[PHASE2_COMP_TIME].time(i) << "\t";
            fout << cumulativeTime << "\t";
            fout << " "  << "\t";
            fout << " "  << "\t";
            fout << " "  << endl;
        }
    }
    else
    {
        for (int i = 0; i < m_TimeStatistics[PHASE2_COMP_TIME].size(); i++)
        {
            cumulativeTime += m_TimeStatistics[PHASE2_COMP_TIME].time(i);
            currTime       += m_TimeStatistics[PHASE2_SIMULATION_TIME_INCREMNT].time(i);


            fout << i + 1 << "\t";
            fout << m_TimeStatistics[PHASE2_SIMULATION_TIME_INCREMNT].time(i) << "\t";
            fout << currTime << "\t";
            fout << m_CountStatistics[PHASE2_FLIPPING_COUNT].count(i) + m_CountStatistics[PHASE2_COLLISION_COUNT].count(i) << "\t";
            fout << m_CountStatistics[PHASE2_FLIPPING_COUNT].count(i) << "\t";
            fout << m_CountStatistics[PHASE2_COLLISION_COUNT].count(i) << "\t";
			fout << m_CountStatistics[PHASE2_COLLISION_BTW_DISKS_COUNT].count(i) << "\t";
            fout << m_TimeStatistics[PHASE2_COMP_TIME].time(i) << "\t";
            fout << cumulativeTime << "\t";
            fout << m_TimeStatistics[PHASE2_MOVING_DISKS_TIME].time(i) << "\t";
            fout << m_TimeStatistics[PHASE2_UPDATING_TOPOLOGY_TIME].time(i) << "\t";
            fout << m_TimeStatistics[PHASE2_UPDATING_GEOMETRY_TIME].time(i) << endl;
        }
    }

    fout.close();

    return;
}


void DynamicVoronoiDiagramCIC::clear_mapper_from_id_to_generator()
{
    m_LinkFromInputIDToGenerator.clear();
}


bool DynamicVoronoiDiagramCIC::these_two_VEdge_are_defined_by_two_same_disks(VEdge2D * VEdge1, VEdge2D * VEdge2) const
{
    bool areTwoVEdgesDefinedByTwoSameDisks = false;

    if (VEdge1->getRightFace() == VEdge2->getRightFace())
    {
        if (VEdge1->getLeftFace() == VEdge2->getLeftFace())
        {
            areTwoVEdgesDefinedByTwoSameDisks = true;
        }
        else
        {
            areTwoVEdgesDefinedByTwoSameDisks = false;
        }
    }
    else if (VEdge1->getLeftFace() == VEdge2->getRightFace())
    {
        if (VEdge1->getRightFace() == VEdge2->getLeftFace())
        {
            areTwoVEdgesDefinedByTwoSameDisks = true;
        }
        else
        {
            areTwoVEdgesDefinedByTwoSameDisks = false;
        }
    }
    else
    {
        areTwoVEdgesDefinedByTwoSameDisks = false;
    }

    return areTwoVEdgesDefinedByTwoSameDisks;
}


bool DynamicVoronoiDiagramCIC::this_VEdge_and_its_hands_and_legs_are_in_this_set(VEdge2D * VEdge, set<VEdge2D*>& VEdgeSet) const
{
    bool areThisVEdgeNItsHandsNLegsInThisSet = true;

    list<VEdge2D*> adjacentVEdges;
    VEdge->getAdjacentVEdges(adjacentVEdges);

    for (list<VEdge2D*>::const_iterator it_VEdge = adjacentVEdges.begin();
        it_VEdge != adjacentVEdges.end();
        it_VEdge++)
    {
        VEdge2D* currVEdge = *it_VEdge;

        if (!this_VEdge_is_in_this_set(currVEdge, VEdgeSet))
        {
            areThisVEdgeNItsHandsNLegsInThisSet = false;
            break;
        }
    }

    return areThisVEdgeNItsHandsNLegsInThisSet;
}


bool DynamicVoronoiDiagramCIC::these_three_VEdge_have_common_disk_in_right_or_left_VFace(const set<VEdge2D*>& VEdges) const
{
    if (VEdges.size()!= 3)
    {
        return false;
    }

    bool areThreeVEdgesDefinedByIdentical4Disks = false;

    set<VEdge2D*>::const_iterator it_VEdge = VEdges.begin();

    VEdge2D* firstVEdge = *it_VEdge;
    VEdge2D* secondVEdge = *(++it_VEdge);
    VEdge2D* thirdVEdge = *(++it_VEdge);

    if ((firstVEdge->getRightFace() == secondVEdge->getLeftFace())
        || (firstVEdge->getRightFace() == secondVEdge->getRightFace()))
    {
        if ((firstVEdge->getRightFace() == thirdVEdge->getLeftFace())
            || (firstVEdge->getRightFace() == thirdVEdge->getRightFace()))
        {
            areThreeVEdgesDefinedByIdentical4Disks = true;
        }
    }
    else if ((firstVEdge->getLeftFace() == secondVEdge->getLeftFace())
        || (firstVEdge->getLeftFace() == secondVEdge->getRightFace()))
    {
        if ((firstVEdge->getLeftFace() == thirdVEdge->getLeftFace())
            || (firstVEdge->getLeftFace() == thirdVEdge->getRightFace()))
        {
            areThreeVEdgesDefinedByIdentical4Disks = true;
        }
    }

    return  areThreeVEdgesDefinedByIdentical4Disks;
}



void DynamicVoronoiDiagramCIC::compute_velocity_vectors_by_scanning_EVENTSEQ()
{
    for (int i = 0; i < m_EventSequence.size(); i++)
    {
        EventOfDynamicVD2D* nextEvent = m_EventSequence.at(i);

        propagate_DVD_by_handling_and_storing_next_event(nextEvent);
    }
}


void DynamicVoronoiDiagramCIC::construct_VD_at_target_Voronoi_clock(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed)
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
		compute_all_VEdges_equation_of_VD();
    }
}



void DynamicVoronoiDiagramCIC::prepare_statistics_in_generating_EVENTSEQ()
{
    m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ = true;

    m_CountStatistics[PHASE1_COUNT].reset(DVD_COUNT_SIZE_PHASE_ONE);
    m_TimeStatistics[PHASE1_COMP_TIME].reset(DVD_TIME_SIZE_PHASE_ONE);
}


void DynamicVoronoiDiagramCIC::___start_to_track_statistics_of_playing_DVD( const double& targetVoronoiClock )
{
    m_CountStatistics[PHASE2_FLIPPING_COUNT].setCount(0);
    m_CountStatistics[PHASE2_COLLISION_COUNT].setCount(0);
    m_CountStatistics[PHASE2_COLLISION_BTW_DISKS_COUNT].setCount(0);
	m_CountStatistics[PHASE2_VELOCITY_CHANGE_COUNT].setCount(0);
	m_CountStatistics[PHASE2_DISK_HOPPING_COUNT].setCount(0);

    m_TimeStatistics[PHASE2_SIMULATION_TIME_INCREMNT].setTime(targetVoronoiClock - get_current_Voronoi_clock());

	___start_clock_in_phase2(PHASE2_COMP_TIME);
}


void DynamicVoronoiDiagramCIC::initialize_statistics_for_collision_avoidance()
{
    m_TimeStatistics[PHASE1_COLLISION_AVOIDANCE_TIME].reset(DVD_COL_AVD_SIZE);
}

void DynamicVoronoiDiagramCIC::complete_statistics_in_generating_EVENTSEQ(const bool& phase1WillBeContinued /*= false*/)
{
    if (!phase1WillBeContinued)
    {
        m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ = false;
    }
    
    m_TimeStatistics[PHASE1_COMP_TIME].setTime(DVD_TIME_PHASE_ONE_TOTAL, m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_CONSTRUCT_INITIAL_CIC_VD) 
                                                                      +  m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS) 
                                                                      +  m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK)
                                                                      +  m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_REINSTATE_THE_INITIAL_VORONOI_DIAGRAM));

    m_TimeStatistics[PHASE1_COMP_TIME].setTime(DVD_TIME_PRIORITYS_BOOKKEEPING, m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING)
                                                                            +  m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION)
                                                                            +  m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_VELOCITY_CHANGE));
}



void DynamicVoronoiDiagramCIC::___stop_to_track_statistics_of_playing_DVD()
{
	___end_clock_in_phase2(PHASE2_COMP_TIME);
}


void DynamicVoronoiDiagramCIC::finalize_statistics_for_collision_avoidance()
{
    complete_statistics_in_generating_EVENTSEQ();

    int totalTime = m_TimeStatistics[PHASE1_COLLISION_AVOIDANCE_TIME].time(DVD_COL_AVD_TOTAL);
    m_TimeStatistics[PHASE1_COMP_TIME].setTime(DVD_TIME_PHASE_ONE_TOTAL, totalTime);
}


void DynamicVoronoiDiagramCIC::___reset_statistics_of_playing_DVD()
{
    m_CountStatistics[PHASE2_FLIPPING_COUNT].reset(0);
    m_CountStatistics[PHASE2_COLLISION_COUNT].reset(0);
    m_CountStatistics[PHASE2_COLLISION_BTW_DISKS_COUNT].reset(0);
    m_CountStatistics[PHASE2_VELOCITY_CHANGE_COUNT].reset(0);
    m_CountStatistics[PHASE2_DISK_HOPPING_COUNT].reset(0);
    m_TimeStatistics[PHASE2_SIMULATION_TIME_INCREMNT].reset(0);
}


void DynamicVoronoiDiagramCIC::set_parameters(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType)
{
    DiskCollisionWatcher diskCollisionWatcher;
    diskCollisionWatcher.insert_watch_object_and_target_collision_number(DiskCollisionWatcher::WATCH_OBJ_FOR_ALL_PAIR_WISE_DISKS, UINT_MAX);

    set_parameters(lengthOfVoronoiClock, diskCollisionWatcher, collisionHandlingType);
}


void DynamicVoronoiDiagramCIC::set_parameters(
    const double& lengthOfVoronoiClock, 
    const DiskCollisionWatcher& collisionWatcher, 
    const DiskCollisionHandler& collisionHandlingType)
{
    set_length_of_Voronoi_clock(lengthOfVoronoiClock);
    set_collision_watcher(collisionWatcher);
    set_collision_handling_type(collisionHandlingType);
}


void DynamicVoronoiDiagramCIC::set_length_of_Voronoi_clock(const double& lengthOfVoronoiClock)
{
    m_LengthOfVoronoiClock = lengthOfVoronoiClock;
}


void DynamicVoronoiDiagramCIC::set_collision_handling_type(const DiskCollisionHandler& collisionHandlingType)
{
    m_CollisionHandlingType = collisionHandlingType;
}


void DynamicVoronoiDiagramCIC::set_collision_watcher(const DiskCollisionWatcher& collisionWatcher)
{
    m_CollisionWatcher = collisionWatcher;
}


void DynamicVoronoiDiagramCIC::set_offset_of_disks(const double& offsetOfDisk)
{
    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    for (list<DynamicDisk*>::iterator it_Disk = dynamicDisks.begin(); it_Disk != dynamicDisks.end(); ++it_Disk)
    {
        DynamicDisk* currDisk = *it_Disk;
        currDisk->setOffset(offsetOfDisk);
    }
}


void DynamicVoronoiDiagramCIC::set_disk_hopping_probability_matrix(const vector<float>& DHPM)
{
	m_DiskHoppingProbabilityMatrix = DHPM;
}



void DynamicVoronoiDiagramCIC::convert_input_for_constructing_VD(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks, list< pair< rg_Circle2D, void* > >& diskNUserDataPairs)
{
    //1. Make and insert constainer pair at the front of a list
    pair< rg_Circle2D, void* > containerDisksNUserDataPair;

    containerDisksNUserDataPair.first = container->getCircle();
    containerDisksNUserDataPair.second = static_cast<void*>(container);

    diskNUserDataPairs.push_back(containerDisksNUserDataPair);

    //2. Make and insert moving disk pair after the container pair in the list
    for (list<DynamicDisk*>::const_iterator it_DynamicDisk = dynamicDisks.begin(); it_DynamicDisk != dynamicDisks.end(); it_DynamicDisk++)
    {
        DynamicDisk* currMovingDisk = (*it_DynamicDisk);

        pair< rg_Circle2D, void* > currDisksNUserDataPair;

        currDisksNUserDataPair.first = currMovingDisk->getCircle();
        currDisksNUserDataPair.second = static_cast<void*>(currMovingDisk);

        diskNUserDataPairs.push_back(currDisksNUserDataPair);
    }
}


void DynamicVoronoiDiagramCIC::convert_input_for_constructing_VD(const list<DynamicDisk*>& containers, const list<DynamicDisk*>& dynamicDisks, list<pair<rg_Circle2D, void*>>& diskNUserDataPairs)
{
    m_NumOfContainers = containers.size();

    //1. Make and insert constainer pair at the front of a list
    for (list<DynamicDisk*>::const_iterator it_Container = containers.begin(); it_Container != containers.end(); it_Container++)
    {
        DynamicDisk* currContainer = (*it_Container);
        currContainer->set_this_disk_is_container(true);

        pair< rg_Circle2D, void* > currContainerNUserDataPair;

        currContainerNUserDataPair.first  = currContainer->getCircle();
        currContainerNUserDataPair.second = static_cast<void*>(currContainer);

        diskNUserDataPairs.push_back(currContainerNUserDataPair);
    }

    //2. Make and insert moving disk pair after the container pair in the list
    for (list<DynamicDisk*>::const_iterator it_DynamicDisk = dynamicDisks.begin(); it_DynamicDisk != dynamicDisks.end(); it_DynamicDisk++)
    {
        DynamicDisk* currMovingDisk = (*it_DynamicDisk);
        currMovingDisk->set_this_disk_is_container(false);

        pair< rg_Circle2D, void* > currDisksNUserDataPair;

        currDisksNUserDataPair.first = currMovingDisk->getCircle();
        currDisksNUserDataPair.second = static_cast<void*>(currMovingDisk);

        diskNUserDataPairs.push_back(currDisksNUserDataPair);
    }
}


void DynamicVoronoiDiagramCIC::copy_initial_VD_and_containers_and_disks_for_reinstatement(const VoronoiDiagramCIC& initialVD, const list<DynamicDisk*>& inputContainers, const list<DynamicDisk*>& inputDisks)
{
    copy_initial_VD_and_make_map_from_copied_generators_to_original_generators(initialVD);
    make_map_from_original_disks_to_copied_disks(inputContainers, inputDisks);
}


void DynamicVoronoiDiagramCIC::copy_initial_VD_and_containers_and_disks_for_reinstatement(const VoronoiDiagramCIC& initialVD, DynamicDisk* const inputContainer, const list<DynamicDisk*>& inputDisks)
{
    list<DynamicDisk*> inputContainers;
    inputContainers.push_back(inputContainer);
    copy_initial_VD_and_containers_and_disks_for_reinstatement(initialVD, inputContainers, inputDisks);
}


void DynamicVoronoiDiagramCIC::copy_initial_VD_and_containers_and_disks_for_reinstatement()
{
    list<DynamicDisk*> containers;
    get_containers(containers);

    list<DynamicDisk*> dynamicDisks;
    get_dynamic_disks(dynamicDisks);

    copy_initial_VD_and_containers_and_disks_for_reinstatement(*this, containers, dynamicDisks);
}


void DynamicVoronoiDiagramCIC::copy_initial_VD_and_make_map_from_copied_generators_to_original_generators(const VoronoiDiagramCIC& initialVD)
{
    m_InitialVDForReinstatement = initialVD;
    
    m_MapFromCopiedInitialStateGeneratorsToOriginalGenerators.clear();

    list<Generator2D*> allGeneratorsOfCopiedVD;
    m_InitialVDForReinstatement.getAllGenerators_hierarchy(allGeneratorsOfCopiedVD);

    list<Generator2D*> allGeneratorsOfOriginalVD;
    initialVD.getAllGenerators_hierarchy(allGeneratorsOfOriginalVD);

    list<Generator2D*>::const_iterator it_CopiedGenerator   = allGeneratorsOfCopiedVD.begin();
    list<Generator2D*>::const_iterator it_OriginalGenerator = allGeneratorsOfOriginalVD.begin();

    for (; it_CopiedGenerator != allGeneratorsOfCopiedVD.end(); ++it_CopiedGenerator, ++it_OriginalGenerator)
    {
        m_MapFromCopiedInitialStateGeneratorsToOriginalGenerators.insert(make_pair(*it_CopiedGenerator, *it_OriginalGenerator));
    }
}


void DynamicVoronoiDiagramCIC::make_map_from_original_disks_to_copied_disks(const list<DynamicDisk*>& inputContainers, const list<DynamicDisk*>& inputDisks)
{
    m_MapFromOriginalDisksToCopiedInitialStateDisks.clear();

    for (list<DynamicDisk*>::const_iterator it_Container = inputContainers.begin();
        it_Container != inputContainers.end();
        ++it_Container)
    {
        DynamicDisk* currContainer = *it_Container;
        m_MapFromOriginalDisksToCopiedInitialStateDisks.insert(make_pair(currContainer, *currContainer));
    }

	for (list<DynamicDisk*>::const_iterator it_Disk = inputDisks.begin();
		it_Disk != inputDisks.end();
		++it_Disk)
	{
		DynamicDisk* currDisk = *it_Disk;
		m_MapFromOriginalDisksToCopiedInitialStateDisks.insert(make_pair(currDisk, *currDisk));
	}
}


void DynamicVoronoiDiagramCIC::link_Voronoi_diagrams(
    const list<VoronoiDiagramCIC*>& allVDs_Copied, 
    const list<VoronoiDiagramCIC*>& allVDs_Current, 
    unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromCopiedToCurrent)
{
    list<VoronoiDiagramCIC*>::const_iterator  i_Copied_VD  = allVDs_Copied.begin();
    list<VoronoiDiagramCIC*>::const_iterator  i_Current_VD = allVDs_Current.begin();

    for (; i_Copied_VD != allVDs_Copied.end(); ++i_Copied_VD, ++i_Current_VD)
    {
        VoronoiDiagramCIC* copiedVD  = *i_Copied_VD;
        VoronoiDiagramCIC* currentVD = *i_Current_VD;
       
        VDMapFromCopiedToCurrent.insert(pair<VoronoiDiagramCIC*, VoronoiDiagramCIC*>(copiedVD, currentVD));
    }
    VDMapFromCopiedToCurrent.insert(pair<VoronoiDiagramCIC*, VoronoiDiagramCIC*>(NULL, NULL));
}


void DynamicVoronoiDiagramCIC::delete_all_VEntities_in_lists_and_clear_the_lists_in_hierarchy(const list<VoronoiDiagramCIC*>& allVDsInHierarchy)
{
    for (list<VoronoiDiagramCIC*>::const_iterator it_VD = allVDsInHierarchy.begin();
        it_VD != allVDsInHierarchy.end();
        ++it_VD)
    {
        (*it_VD)->clearVEntities();
    }
}


void DynamicVoronoiDiagramCIC::clear_all_generator_lists_in_hierarchy_while_the_dynamic_memory_of_the_generators_remains(const list<VoronoiDiagramCIC*>& allVDsInHierarchy)
{
    for (list<VoronoiDiagramCIC*>::const_iterator it_VD = allVDsInHierarchy.begin();
        it_VD != allVDsInHierarchy.end();
        ++it_VD)
    {
        (*it_VD)->clearGeneratorsWithoutDeletingDynamicMemoryOfTheGenerators();
    }
}


void DynamicVoronoiDiagramCIC::get_dynamic_disks_pointers(const list<DynamicDisk>& dynamicDisks, list<DynamicDisk*>& dynamicDiskPointers)
{
    for (list<DynamicDisk>::const_iterator it_DynamicDisk = dynamicDisks.begin();
         it_DynamicDisk != dynamicDisks.end();
         it_DynamicDisk++)
    {
        const DynamicDisk* currMovingDisk = &(*it_DynamicDisk);

        dynamicDiskPointers.push_back(const_cast<DynamicDisk*>(currMovingDisk));
    }
}


void DynamicVoronoiDiagramCIC::construct_initial_VD(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks, const bool& VEdgeGeometryIsComputed /*= false*/)
{
    DynamicDisk* containerPointer = const_cast<DynamicDisk*>(&container);
    list<DynamicDisk*> dynamicDiskPointers;
    get_dynamic_disks_pointers(dynamicDisks, dynamicDiskPointers);

    construct_initial_VD(containerPointer, dynamicDiskPointers,VEdgeGeometryIsComputed);
}


void DynamicVoronoiDiagramCIC::construct_initial_VD(const list<DynamicDisk>& containers, const list<DynamicDisk>& dynamicDisks, const bool& VEdgeGeometryIsComputed /*= false*/)
{
    list<DynamicDisk*> containerPointers;
    get_dynamic_disks_pointers(containers, containerPointers);

    list<DynamicDisk*> dynamicDiskPointers;
    get_dynamic_disks_pointers(dynamicDisks, dynamicDiskPointers);

    construct_initial_VD(containerPointers, dynamicDiskPointers, VEdgeGeometryIsComputed);
}


void DynamicVoronoiDiagramCIC::construct_initial_VD(const list<DynamicDisk*>& containers, const list<DynamicDisk*>& dynamicDisks, const bool& VEdgeGeometryIsComputed /*= false*/)
{
	___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_CIC_VD);

    m_NumOfContainers    = containers.size();
    
    list< pair< rg_Circle2D, void* > > diskNUserDataPairs;
    convert_input_for_constructing_VD(containers, dynamicDisks, diskNUserDataPairs);

    if (m_NumOfContainers > 1)
    {
        constructVoronoiDiagramCIC_hierarchy(diskNUserDataPairs);
    }
    else //(m_NumOfContainers == 1 )
    {
        constructVoronoiDiagramCIC(diskNUserDataPairs);
    }
    
    if (VEdgeGeometryIsComputed)
    {
		compute_all_VEdges_equation_of_VD();
    }
    
    link_dynamic_disk_id_to_its_generator_pointer();

	copy_initial_VD_and_containers_and_disks_for_reinstatement(*this, containers, dynamicDisks);

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_CIC_VD);
}


void DynamicVoronoiDiagramCIC::construct_initial_VD(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks, const bool& VEdgeGeometryIsComputed /*= false*/)
{
    list<DynamicDisk*> containerIncludingOne;
    containerIncludingOne.push_back(container);
    construct_initial_VD(containerIncludingOne, dynamicDisks);
}


void DynamicVoronoiDiagramCIC::construct_initial_priority_Qs_for_both_flipping_and_collision_events()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS);


    construct_initial_priorityQ_of_disk_collision_time();

    construct_initial_priorityQ_of_vedge_flipping_time();
    

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS);
}




void DynamicVoronoiDiagramCIC::construct_initial_priorityQ_of_vedge_flipping_time()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME);


    list<VEdge2D*> allVEdges;
    getVoronoiEdges(allVEdges);

    for (list<VEdge2D*>::const_iterator it_VEdge = allVEdges.begin();
        it_VEdge != allVEdges.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleEdgeFlipingTime(false, DBL_MAX);

        if (VEdge->isInfinite() || VEdge->isUnBounded()) 
        {
            possibleEdgeFlipingTime.first = false;
        }
        else if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
        {
            possibleEdgeFlipingTime.first = false;
        }
        else
        {
            possibleEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(m_CurrentVoronoiClock, VEdge);
        }

        do_bookKeeping_on_priorityQ_for_edge_flipping(possibleEdgeFlipingTime, VEdge, DURING_INITIAL_PQ_CONSTRUCTION);
    }


    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME);
}




void DynamicVoronoiDiagramCIC::construct_initial_priorityQ_of_disk_collision_time()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME);

    list<VEdge2D*> allVEdges;
    getVoronoiEdges(allVEdges);

    for (list<VEdge2D*>::const_iterator it_VEdge = allVEdges.begin();
        it_VEdge != allVEdges.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleDiskCollisionTime(false, DBL_MAX);

        if (   VEdge->isInfinite() || VEdge->isUnBounded()
            || VEdge->getLeftFace()->isUnBounded() || VEdge->getRightFace()->isUnBounded()) // dynamicDisk ���忡���� phantom�� �𸣴� ��
        {
            possibleDiskCollisionTime.first = false;
        }
        else
        {
            DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(VEdge->getLeftFace()->getGenerator()->getUserData());
            DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(VEdge->getRightFace()->getGenerator()->getUserData());

            possibleDiskCollisionTime = find_possible_collision_time_of_disks_after(m_CurrentVoronoiClock, movingDisk1, movingDisk2);
        }

        do_bookkeeping_on_priorityQ_for_disk_collision(possibleDiskCollisionTime, VEdge, DURING_INITIAL_PQ_CONSTRUCTION);
    }


    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME);
}


void DynamicVoronoiDiagramCIC::clear_priority_queues_for_collision_and_flip()
{
    m_PriorityQForDiskCollision.clear();
    m_PriorityQForEdgeFlipping.clear();

    m_InitialPQIsConstructed = false;
}


ReservationNumber DynamicVoronoiDiagramCIC::insert_disk_velocity_change_event_into_priorityQ(const double& eventTime, Generator2D* generator, const rg_Point2D& velocityVector)
{
    ReservationNumber reservedNumber = generate_reservation_number_for_velocity_vector_change();

    m_ReservedVelocityVectorChange[reservedNumber] = make_pair(generator, velocityVector);

    m_PriorityQForDiskVelocityChange.push(reservedNumber, eventTime);

    return reservedNumber;
}


bool DynamicVoronoiDiagramCIC::disk_velocity_change_priorityQ_is_empty() const
{
    return m_ReservedVelocityVectorChange.empty();
}


void DynamicVoronoiDiagramCIC::remove_disk_velocity_change_event_from_priorityQ(const ReservationNumber& reservationNum)
{
    m_PriorityQForDiskVelocityChange.killNodeRelatedWith(reservationNum);
    m_ReservedVelocityVectorChange.erase(reservationNum);
}


bool DynamicVoronoiDiagramCIC::this_number_is_reserved(const ReservationNumber& reservationNum) const
{
    if (m_ReservedVelocityVectorChange.find(reservationNum) != m_ReservedVelocityVectorChange.end())
    {
        return true;
    }
    else
    {
        return false;
    }
}


ReservationNumber DynamicVoronoiDiagramCIC::generate_reservation_number_for_velocity_vector_change()
{
    while (this_number_is_reserved(m_LastReservedNumber))
    {
        if (m_LastReservedNumber == UINT_MAX)
        {
            m_LastReservedNumber = 0;
        }
        else
        {
            ++m_LastReservedNumber;
        }
    }

    return m_LastReservedNumber;
}


void DynamicVoronoiDiagramCIC::propagate_DVD_to_the_last_moment_of_Voronoi_clock()
{
#ifdef WRITING_EVENT_TIME
    ___initialize_writing_event_time();
#endif //WRITING_EVENT_TIME

    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);

    while ( get_next_event_clock_in_priorityQ() <= get_length_of_Voronoi_clock() )
	{
        EventOfDynamicVD2D* nextEvent = get_next_event_from_priorityQs();

        propagate_DVD_by_handling_and_storing_next_event(nextEvent);
    }

    synchronize_disk_locations_and_current_clock_to_VD( get_length_of_Voronoi_clock() );

    compute_all_VVertex_coordinates_of_VD();

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);
}


void DynamicVoronoiDiagramCIC::propagate_DVD_until_target_collision_numbers_btw_disks(const int& numOfCollisions)
{
#ifdef WRITING_EVENT_TIME
    ___initialize_writing_event_time();
#endif //WRITING_EVENT_TIME

    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);

    while ( ___count_collisions_btw_disks() < numOfCollisions)
    {
        EventOfDynamicVD2D* nextEvent = get_next_event_from_priorityQs();

        propagate_DVD_by_handling_and_storing_next_event(nextEvent);
    }

    synchronize_disk_locations_and_current_clock_to_VD(get_current_Voronoi_clock());
    
    compute_all_VVertex_coordinates_of_VD();

    set_length_of_Voronoi_clock(get_current_Voronoi_clock());

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);
}


void DynamicVoronoiDiagramCIC::propagate_DVD_until_target_collision_numbers_of_watch_objects()
{
#ifdef WRITING_EVENT_TIME
    ___initialize_writing_event_time();
#endif //WRITING_EVENT_TIME

    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);

    while ( !m_CollisionWatcher.all_of_target_collision_numbers_are_satisfied() )
    {
        EventOfDynamicVD2D* nextEvent = get_next_event_from_priorityQs();

        propagate_DVD_by_handling_and_storing_next_event(nextEvent);
    }

    synchronize_disk_locations_and_current_clock_to_VD(get_current_Voronoi_clock());
    
    compute_all_VVertex_coordinates_of_VD();

    set_length_of_Voronoi_clock(get_current_Voronoi_clock());

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);
}



DiskCollisionEvent* DynamicVoronoiDiagramCIC::propagate_DVD_to_the_next_collision(const double& endVoronoiClock)
{
#ifdef WRITING_EVENT_TIME
    ___initialize_writing_event_time();
#endif //WRITING_EVENT_TIME

    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);

    DiskCollisionEvent* output = NULL;

    bool collisionIsOccured = false;

    while (!collisionIsOccured && get_next_event_clock_in_priorityQ() <= endVoronoiClock)
    {
        EventOfDynamicVD2D* nextEvent = get_next_event_from_priorityQs();

        propagate_DVD_by_handling_and_storing_next_event(nextEvent);

        if (nextEvent->get_event_type() == EventOfDynamicVD2D::DISK_COLLISION)
        {
            DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(const_cast<EventOfDynamicVD2D*>(nextEvent));
            DynamicDisk* disk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
            DynamicDisk* disk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

            if (!disk1->this_disk_is_container() && !disk2->this_disk_is_container())
            {
                collisionIsOccured = true;
                output = diskCollisionEvent;
            }
        }
    }

    double currentVoronoiClock = get_current_Voronoi_clock();

    if (!collisionIsOccured)
    {   //constrcut DVD by endVoronoiClock
        while (get_next_event_clock_in_priorityQ() <= endVoronoiClock)
        {
            EventOfDynamicVD2D* nextEvent = get_next_event_from_priorityQs();

            propagate_DVD_by_handling_and_storing_next_event(nextEvent);
        }

        currentVoronoiClock = endVoronoiClock;
    }

    synchronize_disk_locations_and_current_clock_to_VD(currentVoronoiClock);

    compute_all_VVertex_coordinates_of_VD();

    set_length_of_Voronoi_clock(currentVoronoiClock);

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);

    return output;
}



tuple<bool, double, DynamicDisk*, DynamicDisk*> DynamicVoronoiDiagramCIC::propagate_DVD_until_next_collision()
{
#ifdef WRITING_EVENT_TIME
    ___initialize_writing_event_time();
#endif //WRITING_EVENT_TIME

    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);

	tuple<bool, double, DynamicDisk*, DynamicDisk*> outputCollidingDisks = make_tuple(false, -DBL_MAX, (DynamicDisk*)NULL, (DynamicDisk*)NULL);
    bool collsionIsOccured = false;

    while (get_next_event_clock_in_priorityQ() <= get_length_of_Voronoi_clock() )
    {
        EventOfDynamicVD2D* nextEvent = get_next_event_from_priorityQs();

        propagate_DVD_by_handling_and_storing_next_event(nextEvent);

        if (nextEvent->get_event_type() == EventOfDynamicVD2D::DISK_COLLISION)
        {
			DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(const_cast<EventOfDynamicVD2D*>(nextEvent));
			DynamicDisk* disk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
			DynamicDisk* disk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

			if (!disk1->this_disk_is_container() && !disk2->this_disk_is_container())
			{
				collsionIsOccured    = true;
				outputCollidingDisks = make_tuple(true, diskCollisionEvent->get_occurring_clock(), disk1, disk2);
			}
        }

        if (collsionIsOccured)
        {
            break;
        }
    }


    return outputCollidingDisks;

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PROPAGATE_DVD_OVER_LENGTH_OF_VORONOI_CLOCK);
}


void V::GeometryTier::DynamicVoronoiDiagramCIC::propagate_DVD_by_handling_and_storing_next_event(EventOfDynamicVD2D* nextEvent)
{  
#ifdef WRITING_EVENT_TIME
    double startTime = ((double)clock() / CLOCKS_PER_SEC);
#endif //WRITING_EVENT_TIME

    switch (nextEvent->get_event_type())
    {
        case EventOfDynamicVD2D::EDGE_FLIP:
        {
            EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(nextEvent);
            handle_and_store_EDGE_FLIP_event(edgeFlippingEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_COLLISION:
        {
            DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(nextEvent);
            handle_and_store_DISK_COLLISION_event(diskCollisionEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
        {
            DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(nextEvent);
            handle_and_store_VELOCITY_VECTOR_CHANGE_event(diskVelocityChangeEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_HOPPING:
        {
            DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(nextEvent);
            handle_and_store_DISK_HOPPING_event(diskHoppingEvent);
            break;
        }
    }

    increment_index_of_most_recent_event();
    set_current_Voronoi_clock(nextEvent->get_occurring_clock());

#ifdef WRITING_EVENT_TIME
    double endTime = ((double)clock() / CLOCKS_PER_SEC);
    ___write_event_times(event->get_occurring_clock(), endTime - startTime);
#endif //WRITING_EVENT_TIME
}


void DynamicVoronoiDiagramCIC::reinstate_the_initial_Voronoi_diagram()
{
    reinstate_the_initial_VD_topology();

    reinstate_the_initial_disk_locations_and_velocity_vectors();

    set_index_of_most_recent_event(-1);

    set_current_Voronoi_clock(0.0);

    clear_current_collision_numbers_of_interest();

    ___reset_statistics_of_playing_DVD();
}


void DynamicVoronoiDiagramCIC::reinstate_the_initial_VD_topology()
{
    list<VoronoiDiagramCIC*> allVDs_Current;
    getAllVoronoiDiagrams_hierarchy(allVDs_Current);

    delete_all_VEntities_in_lists_and_clear_the_lists_in_hierarchy(allVDs_Current);

    clear_all_generator_lists_in_hierarchy_while_the_dynamic_memory_of_the_generators_remains(allVDs_Current);

    list<VoronoiDiagramCIC*> allVDs_Copied;
    list<Generator2D*>       allGeneratorsInHierarcicalVD_Copied;
    list<VFace2D*>           allFacesInHierarcicalVD_Copied;
    list<VEdge2D*>           allEdgesInHierarcicalVD_Copied;
    list<VVertex2D*>         allVerticesInHierarcicalVD_Copied;

    getAllVEntities(m_InitialVDForReinstatement,
        allVDs_Copied,
        allGeneratorsInHierarcicalVD_Copied,
        allFacesInHierarcicalVD_Copied,
        allEdgesInHierarcicalVD_Copied,
        allVerticesInHierarcicalVD_Copied);

    unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*> VDMapFromCopiedToCurrent;
    unordered_map<VFace2D*, VFace2D*>                     faceMapFromCopiedToNew;
    unordered_map<VEdge2D*, VEdge2D*>                     edgeMapFromCopiedToNew;
    unordered_map<VVertex2D*, VVertex2D*>                 vertexMapFromCopiedToNew;


    link_Voronoi_diagrams(allVDs_Copied, allVDs_Current, VDMapFromCopiedToCurrent);
    createNLinkFaces(allFacesInHierarcicalVD_Copied, faceMapFromCopiedToNew);
    createNLinkEdges(allEdgesInHierarcicalVD_Copied, edgeMapFromCopiedToNew);
    createNLinkVertices(allVerticesInHierarcicalVD_Copied, vertexMapFromCopiedToNew);


    insertCreatedEntitiesToVDs(allVDs_Copied,
        VDMapFromCopiedToCurrent,
        m_MapFromCopiedInitialStateGeneratorsToOriginalGenerators,
        faceMapFromCopiedToNew,
        edgeMapFromCopiedToNew,
        vertexMapFromCopiedToNew);


    connectCreatedEntities(allGeneratorsInHierarcicalVD_Copied,
        allFacesInHierarcicalVD_Copied,
        allEdgesInHierarcicalVD_Copied,
        allVerticesInHierarcicalVD_Copied,
        VDMapFromCopiedToCurrent,
        m_MapFromCopiedInitialStateGeneratorsToOriginalGenerators,
        faceMapFromCopiedToNew,
        edgeMapFromCopiedToNew,
        vertexMapFromCopiedToNew);
}


void DynamicVoronoiDiagramCIC::reinstate_the_initial_disk_locations_and_velocity_vectors()
{
    list<Generator2D*> allGenerators;
    getAllGenerators_hierarchy(allGenerators);

    for (list<Generator2D*>::const_iterator it_Generator = allGenerators.begin();
         it_Generator != allGenerators.end();
         ++it_Generator)
    {
        Generator2D* currGenerator   = *it_Generator;
        DynamicDisk* currDynamicDisk = static_cast<DynamicDisk*>(currGenerator->getUserData());

        DynamicDisk  initialStateDisk = m_MapFromOriginalDisksToCopiedInitialStateDisks.at(currDynamicDisk);
      
        *currDynamicDisk = initialStateDisk;
        currGenerator->setDisk(initialStateDisk);
    }
}


void DynamicVoronoiDiagramCIC::clear_current_collision_numbers_of_interest()
{
    m_CollisionWatcher.clear_current_collision_number();
}

void DynamicVoronoiDiagramCIC::find_influenced_disks_by_this_FLIP_event(set<Generator2D*>& influencedDiskGenerators, EdgeFlipEvent* anEvent)
{
    //Include all generators defining 5 influenced VEdges. 
    //These all VEdges also cover an VEdge by 1 influenced collision event, all kind of same flipping VEdges.

    VEdge2D* targetVEdge = get_VEdge_with_this_generator_quadruplet(anEvent->get_generator_quadruplet_before_flip(), false);

    set<VEdge2D*> influecedVEdges;
    collect_influenced_5VEdges_by_flipping_target_vedge_including_itself(targetVEdge, influecedVEdges);

    for (set<VEdge2D*>::const_iterator it_VEdge = influecedVEdges.begin();
        it_VEdge != influecedVEdges.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;
        find_4_definining_disks_of_VEdge(influencedDiskGenerators, VEdge);
    }
}



void DynamicVoronoiDiagramCIC::find_influenced_disks_by_this_DISK_COLLISION_event(set<Generator2D*>& influencedDiskGenerators, DiskCollisionEvent* anEvent)
{
    //Include all generators defining quill VEdges of VCells of colling two dynamicDisks.
    //These all VEdges also cover VEdges by several influenced collision event.

    Generator2D* diskGenerator1 = anEvent->get_disk_generator1();
    Generator2D* diskGenerator2 = anEvent->get_disk_generator2();

    set<VEdge2D*> quillVEdgesOfVCellsOfCollidingTwoDisks;
    collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(diskGenerator1, diskGenerator2,
        quillVEdgesOfVCellsOfCollidingTwoDisks);

    for (set<VEdge2D*>::const_iterator it_VEdge = quillVEdgesOfVCellsOfCollidingTwoDisks.begin();
        it_VEdge != quillVEdgesOfVCellsOfCollidingTwoDisks.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;
        find_4_definining_disks_of_VEdge(influencedDiskGenerators, VEdge);
    }
}


void DynamicVoronoiDiagramCIC::find_influenced_disks_by_this_VELOCITY_VECTOR_CHANGE_event(set<Generator2D*>& influencedDiskGenerators, DiskVelocityChangeEvent* anEvent)
{
    Generator2D* diskGenerator = anEvent->get_disk_generator();

    set<VEdge2D*> quillVEdgesOfVCell;
    collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_VEdge_flipping_times(diskGenerator, quillVEdgesOfVCell);

    for (set<VEdge2D*>::const_iterator it_VEdge = quillVEdgesOfVCell.begin();
        it_VEdge != quillVEdgesOfVCell.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;
        find_4_definining_disks_of_VEdge(influencedDiskGenerators, VEdge);
    }
}


void DynamicVoronoiDiagramCIC::find_influenced_disks_by_this_DISK_HOPPING_event(set<Generator2D*>& influencedDiskGenerators, DiskHoppingEvent* anEvent)
{
    Generator2D* fromThisContainer = anEvent->get_in_what_container_this_generator_is_deleted();
    Generator2D* toThisContainer   = anEvent->get_in_what_container_this_generator_is_inserted();

    list<Generator2D*> generatorsInFromContainer;
    fromThisContainer->getInnerGens(generatorsInFromContainer);

    list<Generator2D*> generatorsInToContainer;
    toThisContainer->getInnerGens(generatorsInToContainer);

    influencedDiskGenerators.insert(generatorsInFromContainer.begin(), generatorsInFromContainer.end());
    influencedDiskGenerators.insert(generatorsInToContainer.begin(), generatorsInToContainer.end());
}



void DynamicVoronoiDiagramCIC::update_the_event_times_in_the_priorityQs_influenced_by_this_FLIP_event(EdgeFlipEvent* edgeFlippingEvent)
{
    double   flipTime      = edgeFlippingEvent->get_occurring_clock();
    VEdge2D* flippingVEdge = NULL;

    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        flippingVEdge = get_VEdge_with_this_generator_quadruplet(edgeFlippingEvent->get_generator_quadruplet_before_flip(), true);
    }
    else
    {
        flippingVEdge = get_VEdge_with_this_generator_quadruplet(edgeFlippingEvent->get_generator_quadruplet_before_flip(), false);
    }

    set<VEdge2D*> influencedVEdges;

#ifndef PQ_STORING_ALL_VEDGES
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        //1. find colllision time : 1 VEdges
        influencedVEdges.insert(flippingVEdge);
        find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
            flipTime, influencedVEdges, REDEFINED_EDGE);

        //2. find flipping time : 5 VEdges + same flipping time VEdges
        collect_influenced_5VEdges_by_flipping_target_vedge_including_itself(flippingVEdge, influencedVEdges);
        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(flippingVEdge, flipTime, influencedVEdges, REDEFINED_EDGE);
    }
    else
    {
        //1. find colllision time : NONE

        //2. find flipping time : 1 VEdge
        influencedVEdges.insert(flippingVEdge);
        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(flippingVEdge, flipTime, influencedVEdges, REDEFINED_EDGE);
    }
#else
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        //1. find colllision time : 1 VEdges
        influencedVEdges.insert(flippingVEdge);
        find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
            flipTime, influencedVEdges);

        //2. find flipping time : 5 VEdges + same flipping time VEdges
        collect_influenced_5VEdges_by_flipping_target_vedge_including_itself(flippingVEdge, influencedVEdges);
        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(flippingVEdge, flipTime, influencedVEdges);
    }
    else
    {
        //1. find colllision time : NONE

        //2. find flipping time : 1 VEdge
        influencedVEdges.insert(flippingVEdge);
        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(flippingVEdge, flipTime, influencedVEdges);
    }
#endif
}


void DynamicVoronoiDiagramCIC::update_the_event_times_in_the_priorityQs_influenced_by_this_DISK_COLLISION_event(DiskCollisionEvent* diskCollisionEvent)
{
    Generator2D* diskGenerator1 = diskCollisionEvent->get_disk_generator1();
    Generator2D* diskGenerator2 = diskCollisionEvent->get_disk_generator2();
    double        collisionTime = diskCollisionEvent->get_occurring_clock();

#ifndef PQ_STORING_ALL_VEDGES
    //1. find colllision time : scale VEdges
    set<VEdge2D*> scaleVEdges;
    collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_PHYSICAL_COLLISION(diskGenerator1, diskGenerator2, scaleVEdges);
    find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
        collisionTime, scaleVEdges, NONE_BUT_DISK_VEl_VEC_CHANGED);

    //2. find flipping time : quill VEdges
    set<VEdge2D*> quillVEdges;
    collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(diskGenerator1, diskGenerator2, quillVEdges);
    find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(collisionTime, quillVEdges, NONE_BUT_DISK_VEl_VEC_CHANGED);
#else
    //1. find colllision time : scale VEdges
    switch (m_CollisionType)
    {
        case AVOIDANCE: /*COVIDSIM CODE IS NOT REFLECTED IN HERE*/
        {
            set<VEdge2D*> sharedVEdges;
            set<VEdge2D*> theOtherVEdges;

            collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_AVOIDANCE(diskGenerator1, diskGenerator2, sharedVEdges, theOtherVEdges);
            find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(collisionTime, theOtherVEdges);
            find_possible_collision_time_of_shared_VEdges_after_input_time_N_store_in_priorityQ_AVOIDANCE(collisionTime, sharedVEdges);
            break;
        }

        default:
        {
            set<VEdge2D*> scaleVEdges;
            collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_PHYSICAL_COLLISION(diskGenerator1, diskGenerator2, scaleVEdges);
            find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(collisionTime, scaleVEdges);
            break;
        }
    }

    //2. find flipping time : quill VEdges
    set<VEdge2D*> quillVEdges;
    collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(diskGenerator1, diskGenerator2, quillVEdges);
    find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(collisionTime, quillVEdges);
#endif

}


void DynamicVoronoiDiagramCIC::update_the_event_times_in_the_priorityQs_influenced_by_this_VELOCITY_VECTOR_CHANGE_event(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    Generator2D* diskGenerator = diskVelocityChangeEvent->get_disk_generator();
    double       collisionTime = diskVelocityChangeEvent->get_occurring_clock();

#ifndef PQ_STORING_ALL_VEDGES
    // 1. find colllision time : scale VEdges
    set<VEdge2D*> scaleVEdges;
    collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_collision_times(diskGenerator, scaleVEdges);
    find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
        collisionTime, scaleVEdges, NONE_BUT_DISK_VEl_VEC_CHANGED);

    //2. find flipping time : quill VEdges
    set<VEdge2D*> quillVEdges;
    collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_VEdge_flipping_times(diskGenerator, quillVEdges);
    find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(collisionTime, quillVEdges, NONE_BUT_DISK_VEl_VEC_CHANGED);

#else
    //1. find colllision time : scale VEdges
    set<VEdge2D*> scaleVEdges;
    collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_collision_times(diskGenerator, scaleVEdges);
    find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(collisionTime, scaleVEdges);

    //2. find flipping time : quill VEdges
    set<VEdge2D*> quillVEdges;
    collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_VEdge_flipping_times(diskGenerator, quillVEdges);
    find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(collisionTime, quillVEdges);
#endif
}


void DynamicVoronoiDiagramCIC::update_the_event_times_in_the_priorityQs_influenced_by_this_DISK_HOPPING_event(DiskHoppingEvent * diskHoppingEvent)
{
    /*Redefined edge�� ���ϰ� �������� ����. collision�� ���ֺ��� disk�� �ٲ��������? Ȯ��!! */
    double  diskHoppingTime = diskHoppingEvent->get_occurring_clock();

    if (diskHoppingEvent->get_does_this_event_happen())
    {
        find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
            diskHoppingTime,
            m_NewVEdgesByDiskHopping,
            NEW_EDGE);

        find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
            diskHoppingTime,
            m_DeletedVEdgesByDiskHopping,
            DELETED_EDGE);

        find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
            diskHoppingTime,
            m_RedefinedVEdgesByDiskHopping,
            REDEFINED_EDGE);

        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
            diskHoppingTime,
            m_NewVEdgesByDiskHopping,
            NEW_EDGE);

        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
            diskHoppingTime,
            m_DeletedVEdgesByDiskHopping,
            DELETED_EDGE);

        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
            diskHoppingTime,
            m_RedefinedVEdgesByDiskHopping,
            REDEFINED_EDGE);
    }
    else
    {
        Generator2D* targetDisk = diskHoppingEvent->get_target_generator();
        Generator2D* container = diskHoppingEvent->get_in_what_container_this_generator_is_deleted();

        //1. find colllision time : scale VEdges
        set<VEdge2D*> scaleVEdges;
        collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_PHYSICAL_COLLISION(targetDisk, container, scaleVEdges);
        find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
            diskHoppingTime, scaleVEdges, NONE_BUT_DISK_VEl_VEC_CHANGED);

        //2. find flipping time : quill VEdges
        set<VEdge2D*> quillNScaleVEdges;
        collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(targetDisk, container, quillNScaleVEdges);
        find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(diskHoppingTime, quillNScaleVEdges, NONE_BUT_DISK_VEl_VEC_CHANGED);
    }
}


void DynamicVoronoiDiagramCIC::compute_velocity_vectors_of_disks_after_collided(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2, rg_Point2D& outputVector1, rg_Point2D& outputVector2) const
{
    m_CollisionHandlingType.compute_velocity_vectors_of_disks_after_collided(movingDisk1, movingDisk2, outputVector1, outputVector2);
}


void DynamicVoronoiDiagramCIC::watch_and_record_collision_to_FUTURE(DiskCollisionEvent* diskCollisionEvent)
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    if (m_CollisionWatcher.all_collision_btw_disks_are_being_watched())
    {
        if(!movingDisk1->this_disk_is_container() && !movingDisk2->this_disk_is_container())
        {
            m_CollisionWatcher.increase_current_collision_number(DiskCollisionWatcher::WATCH_OBJ_FOR_ALL_PAIR_WISE_DISKS);
        }
    }
    else
    {
        watch_and_record_collision_for_group_to_FUTURE(movingDisk1, movingDisk2);
        watch_and_record_collision_for_disk_to_FUTURE(movingDisk1, movingDisk2);
    }
}


void DynamicVoronoiDiagramCIC::watch_and_record_collision_for_group_to_FUTURE(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2)
{
    DiskCollisionWatcher::WatchObject watchObjectCandidate 
        = DiskCollisionWatcher::WatchObject(DiskCollisionWatcher::WatchType::GROUP_ID_PAIR, 
                                            make_pair(dynamicDisk1->getGroupId(), dynamicDisk2->getGroupId()));

    if (m_CollisionWatcher.this_is_watch_object(watchObjectCandidate))
    {
        m_CollisionWatcher.increase_current_collision_number(watchObjectCandidate);
    }

/*
    if(dynamicDisk1->getGroupId() == 2 
    && dynamicDisk2->getGroupId() == 2)
    {
        cout << "here: both group 2" << endl;
    }
*/

}


void DynamicVoronoiDiagramCIC::watch_and_record_collision_for_disk_to_FUTURE(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2)
{
    DiskCollisionWatcher::WatchObject watchObjectCandidate 
        = DiskCollisionWatcher::WatchObject(DiskCollisionWatcher::WatchType::DISK_ID_PAIR, 
                                            make_pair(dynamicDisk1->getID(), dynamicDisk2->getID()));

    if (m_CollisionWatcher.this_is_watch_object(watchObjectCandidate))
    {
        m_CollisionWatcher.increase_current_collision_number(watchObjectCandidate);
    }
}


void DynamicVoronoiDiagramCIC::watch_and_record_collision_to_PAST(DiskCollisionEvent* diskCollisionEvent)
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    if (m_CollisionWatcher.all_collision_btw_disks_are_being_watched())
    {
        m_CollisionWatcher.decrease_current_collision_number(DiskCollisionWatcher::WATCH_OBJ_FOR_ALL_PAIR_WISE_DISKS);
    }
    else
    {
        watch_and_record_collision_for_group_to_PAST(movingDisk1, movingDisk2);
        watch_and_record_collision_for_disk_to_PAST(movingDisk1, movingDisk2);
    }
}

void DynamicVoronoiDiagramCIC::watch_and_record_collision_for_group_to_PAST(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2)
{
    DiskCollisionWatcher::WatchObject watchObjectCandidate
        = DiskCollisionWatcher::WatchObject(DiskCollisionWatcher::WatchType::GROUP_ID_PAIR,
            make_pair(dynamicDisk1->getGroupId(), dynamicDisk2->getGroupId()));

    if (m_CollisionWatcher.this_is_watch_object(watchObjectCandidate))
    {
        m_CollisionWatcher.decrease_current_collision_number(watchObjectCandidate);
    }
}


void DynamicVoronoiDiagramCIC::watch_and_record_collision_for_disk_to_PAST(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2)
{
    DiskCollisionWatcher::WatchObject watchObjectCandidate
        = DiskCollisionWatcher::WatchObject(DiskCollisionWatcher::WatchType::DISK_ID_PAIR,
            make_pair(dynamicDisk1->getID(), dynamicDisk2->getID()));

    if (m_CollisionWatcher.this_is_watch_object(watchObjectCandidate))
    {
        m_CollisionWatcher.decrease_current_collision_number(watchObjectCandidate);
    }
}


std::pair<bool, double> DynamicVoronoiDiagramCIC::find_possible_flipping_time_of_an_vedge_after(
    const double& mostRecentlyEventOccurringTime, 
    const VEdge2D* targetVEdge, 
    const bool& isThisVEdgeFlipping /*= false*/)
{	
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_FIND_POSSIBLE_FLIPPING_TIME_OF_EDGE);

    vector<DynamicDisk*> movingDisksForVEdge;
    get_4Disks_defining_this_VEdge_and_its_end_points(targetVEdge, movingDisksForVEdge);
    sort_generators_in_non_increasing_order(movingDisksForVEdge);

    rg_Polynomial polynomialEq = make_polynomial_equation_for_vedge_flipping(movingDisksForVEdge);
	
    pair<bool, double> possibleEdgeFlipTime = solve_polynomial_equation_for_vedge_flipping(polynomialEq,
                                                                                           mostRecentlyEventOccurringTime,
                                                                                           movingDisksForVEdge,
                                                                                           isThisVEdgeFlipping);

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_FIND_POSSIBLE_FLIPPING_TIME_OF_EDGE);

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_EDGE_FLIPPING_IN_COMPUTATION);

    return possibleEdgeFlipTime;
}



rg_Polynomial DynamicVoronoiDiagramCIC::make_polynomial_equation_for_vedge_flipping(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_MAKING_EDGE_FLIPPING_POLYNOMIAL_EQ);

    rg_Polynomial result = EdgeFlipTimeCalculator::make_polynomial_equation_for_vedge_flipping(fourDiskGeneratorsInNonIncreasingOrder);

    int newDegree = compute_polynomial_degree_by_disk_info(fourDiskGeneratorsInNonIncreasingOrder);

    if (newDegree != result.getDegree())
    {
        adjust_polynomial_degree(newDegree, result);
    }

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_MAKING_EDGE_FLIPPING_POLYNOMIAL_EQ);

    ___check_polynomial_degree_for_statistics(result);

	return result;
}


void DynamicVoronoiDiagramCIC::adjust_polynomial_degree(const int& newDegree, rg_Polynomial& polynomial)
{
    double newCoeff[9];

    for (int i = 0; i <= newDegree; i++)
    {
        newCoeff[i] = polynomial.getCoefficient(i);
    }

    polynomial.setPolynomial(newDegree, newCoeff);

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE_REDUCING);
}


std::pair<bool, double> DynamicVoronoiDiagramCIC::solve_polynomial_equation_for_vedge_flipping(
    const rg_Polynomial& polynomialEq, 
    const double& mostRecentlyEventOccurringTime, 
    const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
    const bool& isThisVEdgeFlipping)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_SOLVING_EDGE_FLIPPING_POLYNOMIAL_EQ);

   pair<bool, double> nextEdgeFlipTime = EdgeFlipTimeCalculator::solve_polynomial_equation_for_vedge_flipping(
                                                                    POLYNOMIAL_SOLVER, 
                                                                    mostRecentlyEventOccurringTime,
                                                                    polynomialEq,
                                                                    fourDiskGeneratorsInNonIncreasingOrder,
                                                                    isThisVEdgeFlipping);
   
    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_SOLVING_EDGE_FLIPPING_POLYNOMIAL_EQ);

    return nextEdgeFlipTime;
}



bool DynamicVoronoiDiagramCIC::this_VEdge_is_NOT_shadow_flipping(VEdge2D* VEdge, const double& eventTime, const rg_REAL resolution /*= rg_MATH_RES*/)
{
    set<Generator2D*> fourDefiningGenerators;
    find_4_definining_disks_of_VEdge(fourDefiningGenerators, VEdge);
    move_disks_to_target_Voronoi_clock_without_changing_VD(fourDefiningGenerators, eventTime);

    VVertex2D* startVtx = VEdge->getStartVertex();
    VVertex2D* endVtx = VEdge->getEndVertex();

    update_geometry_of_VVertex(startVtx);
    update_geometry_of_VVertex(endVtx);
  
    rg_Point2D startPt(startVtx->getLocation());
    rg_Point2D endPt(endVtx->getLocation());

    double distance = startPt.distance(endPt);

    if (distance <= resolution)
    {
        return true;
    }
    else
    {
        return false;
    }
}


VEdge2D* DynamicVoronoiDiagramCIC::get_next_flipping_edge_from_flipping_priorityQ()
{
	VEdge2D* toBeFlippedVEdge = m_PriorityQForEdgeFlipping.topNode()->getEntity();
	double   flippingTime     = m_PriorityQForEdgeFlipping.topNode()->getKey();

	//if there are a lot of vedges defined by same disks, find a real flipping vedge.
	list<VEdge2D*> candidates;
	m_PriorityQForEdgeFlipping.findEntitiesHavingSamekey(flippingTime, candidates, TOLERANCE_OF_TOP_ELEMENT_WITH_IDENTICAL_KEY_IN_PQ);
	
    set<VEdge2D*> candidateSet;
    candidateSet.insert(candidates.begin(), candidates.end());
    find_all_VEdges_defined_by_identical_disks_TOPOLOGICAL_COMPUTATION(m_PriorityQForEdgeFlipping.topNode()->getEntity(), candidateSet, m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ);
   
    //find_all_VEdges_defined_by_identical_disks(m_PriorityQForEdgeFlipping.topNode()->getEntity(), candidates, m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ);

    switch (m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.size())
    {
    case 1:
    {
        toBeFlippedVEdge = (*m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin());
        break;
    }

    default:
    {
        // added - 180306
        // 0. move related disks
        VEdge2D* oneOfSameFlippingTimeVEdge = *m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin();

        set<Generator2D*> generatorsForFindingFlippingVEdge;
        find_4_definining_disks_of_VEdge(generatorsForFindingFlippingVEdge, oneOfSameFlippingTimeVEdge);
        move_disks_to_target_Voronoi_clock_without_changing_VD(generatorsForFindingFlippingVEdge, flippingTime);


        // 1. find VVertex
        set<VVertex2D*> VVerticesExceptInfiniteVVertices;

        for (set<VEdge2D*>::const_iterator it_VEdge = m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin();
            it_VEdge != m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.end();
            it_VEdge++)
        {
            VVertex2D* currVertex = (*it_VEdge)->getStartVertex();

            if (!currVertex->isInfinite())
            {
                VVerticesExceptInfiniteVVertices.insert(const_cast<VVertex2D*>(currVertex));
            }

            currVertex = (*it_VEdge)->getEndVertex();

            if (!currVertex->isInfinite())
            {
                VVerticesExceptInfiniteVVertices.insert(const_cast<VVertex2D*>(currVertex));
            }
        }

        // 2. compute VVertex Coordinates
        for (set<VVertex2D*>::const_iterator it_VVertex = VVerticesExceptInfiniteVVertices.begin();
            it_VVertex != VVerticesExceptInfiniteVVertices.end();
            it_VVertex++)
        {
            update_geometry_of_VVertex(*it_VVertex);
        }


        // 2. compute VEdge distance
        // 3. choose the smallest distance VEdge

        toBeFlippedVEdge = find_next_flipping_vedge_when_vedges_have_same_flipping_time();


        // for stat
        ___check_number_of_VEdges_defined_by_identical_4disks_for_statistics();

        break;
    }
    }

	return  toBeFlippedVEdge;
}



void DynamicVoronoiDiagramCIC::find_all_VEdges_defined_by_identical_disks(const VEdge2D* VEdge, const list<VEdge2D*>& candidates, set<VEdge2D*>& outputs)
{
	vector<DynamicDisk*> generatorsForCriterion;
    get_4Disks_defining_this_VEdge_and_its_end_points(VEdge, generatorsForCriterion);
    sort(generatorsForCriterion.begin(), generatorsForCriterion.end(), DynamicVoronoiDiagramCIC::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order);

	for (list<VEdge2D*>::const_iterator it_VEdge = candidates.begin();
		it_VEdge != candidates.end();
		it_VEdge++)
	{
        VEdge2D* currVEdge = *it_VEdge;

        vector<DynamicDisk*>  generators;
        get_4Disks_defining_this_VEdge_and_its_end_points(currVEdge, generators);
        sort(generators.begin(), generators.end(), DynamicVoronoiDiagramCIC::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order);

		int count = 0;

        for (int i = 0; i < generatorsForCriterion.size(); i++)
		{
			if (generatorsForCriterion.at(i) == generators.at(i))
			{
				count++;
			}
			else
			{
				break;
			}
		}

		if (count == 4)
		{
			outputs.insert(*it_VEdge);
		}
	}
}

void DynamicVoronoiDiagramCIC::find_all_VEdges_defined_by_identical_disks_TOPOLOGICAL_COMPUTATION(const VEdge2D * VEdge, const set<VEdge2D*>& candidates, set<VEdge2D*>& outputs) const
{
    vector<Generator2D*> generatorsForCriterion;
    get_4DiskGenerators_defining_this_VEdge_and_its_end_points(VEdge, generatorsForCriterion);
    sort(generatorsForCriterion.begin(), generatorsForCriterion.end());

    for (set<VEdge2D*>::const_iterator it_VEdge = candidates.begin();
         it_VEdge != candidates.end();
         it_VEdge++)
    {
        VEdge2D* currVEdge = *it_VEdge;

        vector<Generator2D*>  generators;
        get_4DiskGenerators_defining_this_VEdge_and_its_end_points(currVEdge, generators);
        sort(generators.begin(), generators.end());

        int count = 0;

        for (int i = 0; i < generatorsForCriterion.size(); i++)
        {
            if (generatorsForCriterion.at(i) == generators.at(i))
            {
                count++;
            }
            else
            {
                break;
            }
        }

        if (count == 4)
        {
            outputs.insert(*it_VEdge);
        }
    }
}



void DynamicVoronoiDiagramCIC::collect_influenced_5VEdges_by_flipping_target_vedge_including_itself(const VEdge2D* targetVEdge, set<VEdge2D*>& influencedVEdgeList) const
{
	influencedVEdgeList.insert(const_cast<VEdge2D*>(targetVEdge->getRightHand()));
    influencedVEdgeList.insert(const_cast<VEdge2D*>(targetVEdge->getLeftHand()));
    influencedVEdgeList.insert(const_cast<VEdge2D*>(targetVEdge->getRightLeg()));
    influencedVEdgeList.insert(const_cast<VEdge2D*>(targetVEdge->getLeftLeg()));

    influencedVEdgeList.insert(const_cast<VEdge2D*>(targetVEdge));
}



void DynamicVoronoiDiagramCIC::find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
    const double& mostRecentlyEventOccurringTime, 
    const set<VEdge2D*>& influencedEdges)
{
	for (set<VEdge2D*>::const_iterator it_VEdge = influencedEdges.begin();
		it_VEdge != influencedEdges.end();
		it_VEdge++)
	{
		const VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleEdgeFlipingTime(false, DBL_MAX);

		if (VEdge->isInfinite() || VEdge->isUnBounded()) //ie, if this vedge is a virtual vedge : ������
		{
			possibleEdgeFlipingTime.first = false;
		}
        else if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
		{
            possibleEdgeFlipingTime.first = false;
		}
		else
		{
			possibleEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(mostRecentlyEventOccurringTime, VEdge);

		}

#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::EDGE_FLIP, const_cast<VEdge2D*>(VEdge), possibleEdgeFlipingTime.second);
#endif //DEBUG_DVD

		do_bookKeeping_on_priorityQ_for_edge_flipping(possibleEdgeFlipingTime, VEdge, AFTER_PQ_CONSTRUCTION);
	}
}

void DynamicVoronoiDiagramCIC::find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(const VEdge2D * targetVEdge, 
                                                                                                                      const double & mostRecentlyEventOccurringTime, 
                                                                                                                      const set<VEdge2D*>& influencedEdges)
{
    pair<bool, double> possibleTargetEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(mostRecentlyEventOccurringTime, targetVEdge, true);


    for (set<VEdge2D*>::const_iterator it_VEdge = influencedEdges.begin();
		it_VEdge != influencedEdges.end();
		it_VEdge++)
	{
		const VEdge2D* VEdge = *it_VEdge;

		pair<bool, double> possibleEdgeFlipingTime(false, DBL_MAX);

		if (VEdge->isInfinite() || VEdge->isUnBounded()) //ie, if this vedge is a virtual vedge : ������
		{
			possibleEdgeFlipingTime.first = false;
		}
        else if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
		{
            possibleEdgeFlipingTime.first = false;
		}
        //else if (this_VEdge_is_in_this_set(VEdge, m_VEdgesDefinedByIdentical4DiskWithFlippingVEdge))
        else if (this_VEdge_is_in_this_set(VEdge, m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ))
        {
            possibleEdgeFlipingTime = possibleTargetEdgeFlipingTime;
        }
		else
		{
            possibleEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(mostRecentlyEventOccurringTime, VEdge);
		}


#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::EDGE_FLIP, const_cast<VEdge2D*>(VEdge), possibleEdgeFlipingTime.second);
#endif //DEBUG_DVD


		do_bookKeeping_on_priorityQ_for_edge_flipping(possibleEdgeFlipingTime, VEdge, AFTER_PQ_CONSTRUCTION);
	}
}



void DynamicVoronoiDiagramCIC::do_bookKeeping_on_priorityQ_for_edge_flipping(const pair<bool, double>& newEdgeFlipTime, const VEdge2D* vEdge, const PriorityQBookkeepingType& PQBookkeepingType)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING);

    bool   doesHaveNextEdgeFlip = newEdgeFlipTime.first;
    double nextEdgeFlipTime     = newEdgeFlipTime.second;

    switch (PQBookkeepingType)
    {
    case DURING_INITIAL_PQ_CONSTRUCTION:
    {
        if (doesHaveNextEdgeFlip)
        {
            m_PriorityQForEdgeFlipping.push(const_cast<VEdge2D*>(vEdge), nextEdgeFlipTime);
        }
        else
        {
            m_PriorityQForEdgeFlipping.push(const_cast<VEdge2D*>(vEdge), DBL_MAX);
        }

        break;
    }

    case  AFTER_PQ_CONSTRUCTION:
    {
        if (doesHaveNextEdgeFlip)
        {
            m_PriorityQForEdgeFlipping.replaceKey(const_cast<VEdge2D*>(vEdge), nextEdgeFlipTime);
        }
        else
        {
            m_PriorityQForEdgeFlipping.replaceKey(const_cast<VEdge2D*>(vEdge), DBL_MAX);
        }

        break;
    }
    default:
        break;
    }

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING);
}

bool DynamicVoronoiDiagramCIC::this_VEdge_is_in_this_set(const VEdge2D * VEdge, const set<VEdge2D*>& VEdgeSet) const
{
    set<VEdge2D*>::const_iterator it_TargetVEdge = VEdgeSet.find(const_cast<VEdge2D*>(VEdge));

    if (it_TargetVEdge == VEdgeSet.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}


double DynamicVoronoiDiagramCIC::get_next_target_time_to_go_back(const double& timeRelativeToCollisionTime, const double& currentVoronoiClock) const
{
    double lastVelocityVectorChangeTime = 0.0;
    for (int i = m_IndexOfMostRecentEvent - 1; i >= 0; --i)
    {
        EventOfDynamicVD2D* currEvent = m_EventSequence[i];
        if (currEvent->get_event_type() == EventOfDynamicVD2D::DISK_VELOCITY_CHANGE)
        {
            lastVelocityVectorChangeTime = currEvent->get_occurring_clock();
            break;
        }
    }

    double timeDiffereceFromCurrTimeToLastVelocityChange = lastVelocityVectorChangeTime - currentVoronoiClock;
    if (abs(timeDiffereceFromCurrTimeToLastVelocityChange) < abs(timeRelativeToCollisionTime))
    {
        return lastVelocityVectorChangeTime;
    }
    else
    {
        return currentVoronoiClock + timeRelativeToCollisionTime;
    }
}


void DynamicVoronoiDiagramCIC::go_back_to_target_Voronoi_clock_and_remove_events_after_the_target_Voronoi_clock(const double& targetVoronoiClock)
{
    progress_DVD_to_PAST(targetVoronoiClock);

    remove_events_in_EVENTSEQ_after_current_Voronoi_clock();
}


void DynamicVoronoiDiagramCIC::remove_disk(DynamicDisk* disk)
{
    //about VD
    Generator2D* targetGenerator = get_generator_by_its_input_disk_id(disk->getID());

    list<VEdge2D*> newVEdges;
    list<VEdge2D*> deletedVEdges;
    list<VEdge2D*> redefinedVEdges;

    updateVoronoiDiagram_removal_deleteFromGeneratorList(
        targetGenerator, newVEdges, deletedVEdges, redefinedVEdges);

    Generator2D* containerGenerator = const_cast<Generator2D*>(getContainerGenerator());
    containerGenerator->removeOneInnerGen(targetGenerator);

    clear_delete_RED_VEdges_and_VVertices();

    //about DVD
    m_LinkFromInputIDToGenerator.erase(disk->getID());
}



void DynamicVoronoiDiagramCIC::remove_events_in_EVENTSEQ_after_current_Voronoi_clock()
{
    const int eventSeqSize = m_EventSequence.size();

    int count = 0;

	for (int i = m_IndexOfMostRecentEvent + 1; i < eventSeqSize; ++i)
	{
		delete m_EventSequence[i];
		++count;
	}

    /*
    for (int i = 0; i < count; ++i)
    {
        m_EventSequence.pop_back();
    }
    */

    int targetSize = eventSeqSize - count;
    m_EventSequence.resize(targetSize);
}


rg_Point2D DynamicVoronoiDiagramCIC::rotate_velocity_vector(const rg_Point2D& inputVector, const double& radian) const
{
	double v_x = cos(radian) * inputVector.getX() - sin(radian) * inputVector.getY();
	double v_y = sin(radian) * inputVector.getX() + cos(radian) * inputVector.getY();

    return rg_Point2D(v_x, v_y);
}


rg_Point2D DynamicVoronoiDiagramCIC::compute_first_disk_velocity_vector_for_collision_avoidance(DynamicDisk* disk1, DynamicDisk* disk2) const
{
    rg_Point2D disk1VelocityVector = disk1->getVelocityVector();

    double half_pi = 0.5 * rg_PI;
   
    return  rotate_velocity_vector(disk1VelocityVector, half_pi);
}



tuple<DynamicDisk*, rg_Point2D, DynamicDisk*, rg_Point2D> DynamicVoronoiDiagramCIC::compute_velocity_vectors_by_brute_force_search_for_collision_avoidance(DynamicDisk* disk1, DynamicDisk* disk2, const double& currentVoronoiClock)
{
    DynamicDisk* fasterDisk = NULL;
    DynamicDisk* slowerDisk = NULL;

    if (disk1->getVelocityVector().magnitude() > disk2->getVelocityVector().magnitude())
    {
        fasterDisk = disk1;
        slowerDisk = disk2;
    }
    else
    {
		fasterDisk = disk2;
		slowerDisk = disk1;
    }

    rg_Point2D originalVelocityVectorOfFasterDisk = fasterDisk->getVelocityVector();
	rg_Point2D originalVelocityVectorOfSlowerDisk = slowerDisk->getVelocityVector();


    const int    NUM_OF_DIRECTION = 30;
	const int    NUM_OF_SPEED     = 20;

    const double ANGLE_INCREMENT  = 2.0 * rg_PI / NUM_OF_DIRECTION;

    
	Generator2D* generatorOfFasterDisk = get_generator_by_its_input_disk_id(fasterDisk->getID());

	list<Generator2D*> firstNeighborGenerators;
	generatorOfFasterDisk->getNeighborGenerators(firstNeighborGenerators, false);

    pair<rg_Point2D, double> goodVelocityVectorNCollisionTimeWithNeighbors( rg_Point2D(0.0, 0.0), -DBL_MAX );

    for (int i = 0; i < NUM_OF_DIRECTION - 1; ++i)
    {
        rg_Point2D candidateVelocityVector = rotate_velocity_vector(fasterDisk->getVelocityVector(), ANGLE_INCREMENT);
   
        for (int j = 0; j < NUM_OF_SPEED; ++j)
		{
			//rg_Point2D speedAdjustedCandidateVelocityVector = (0.5 + 0.1 * j) * candidateVelocityVector;
			rg_Point2D speedAdjustedCandidateVelocityVector = 0.1 * (j + 1) * candidateVelocityVector;

			fasterDisk->setVelocityVector(speedAdjustedCandidateVelocityVector);

			double possibleFirstOccurentCollisionWithNeighbors = DBL_MAX;

			for (list<Generator2D*>::const_iterator it_FirstNeighbor = firstNeighborGenerators.begin();
				it_FirstNeighbor != firstNeighborGenerators.end();
				++it_FirstNeighbor)
			{
				Generator2D* currGenerator = *it_FirstNeighbor;
				DynamicDisk* currDisk = static_cast<DynamicDisk*>(currGenerator->getUserData());

				pair<bool, double> possibleCollsionTime = find_possible_collision_time_of_disks_after(currentVoronoiClock, fasterDisk, currDisk);
				if (possibleCollsionTime.first)
				{
					double possibleCollsionTimeWithOne = possibleCollsionTime.second;

					if (possibleCollsionTimeWithOne < possibleFirstOccurentCollisionWithNeighbors)
					{
						possibleFirstOccurentCollisionWithNeighbors = possibleCollsionTimeWithOne;
					}
				}
			}

			if (possibleFirstOccurentCollisionWithNeighbors > goodVelocityVectorNCollisionTimeWithNeighbors.second)
			{
				goodVelocityVectorNCollisionTimeWithNeighbors = make_pair(speedAdjustedCandidateVelocityVector, possibleFirstOccurentCollisionWithNeighbors);
			}
        }  
    }


    fasterDisk->setVelocityVector(originalVelocityVectorOfFasterDisk);
    slowerDisk->setVelocityVector(originalVelocityVectorOfSlowerDisk);


    rg_Point2D outputVelocityVectorForFasterDisk = goodVelocityVectorNCollisionTimeWithNeighbors.first;
    return make_tuple(fasterDisk, outputVelocityVectorForFasterDisk, slowerDisk, slowerDisk->getVelocityVector());
}


DiskVelocityChangeEvent* DynamicVoronoiDiagramCIC::create_velocity_vector_change_event(const double& eventTime, DynamicDisk* disk, const rg_Point2D& newVelocityVector) const
{
    Generator2D* generatorDisk = get_generator_by_its_input_disk_id(disk->getID());

	DiskVelocityChangeEvent* diskVelocityChangeEvent = new DiskVelocityChangeEvent(eventTime, generatorDisk, newVelocityVector);

    return diskVelocityChangeEvent;
}



VEdge2D* DynamicVoronoiDiagramCIC::find_next_flipping_vedge_when_vedges_have_same_flipping_time()
{
	VEdge2D* outputVEdge = nullptr;
	double minDistance = DBL_MAX;

	for (set<VEdge2D*>::const_iterator it_VEdge = m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin();
		it_VEdge != m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.end();
		it_VEdge++)
	{
		VEdge2D* VEdge = *it_VEdge;

		rg_Point2D startPt = VEdge->getStartVertex()->getLocation();
		rg_Point2D endPt = VEdge->getEndVertex()->getLocation();

		double currDistance = startPt.distance(endPt);
        
		if (currDistance < minDistance)
		{
			minDistance = currDistance;
			outputVEdge = VEdge;
		}
	}

	return outputVEdge;
}




void DynamicVoronoiDiagramCIC::clear_info_about_vedges_defined_by_same_disks()
{
	m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.clear();
    m_VEdgesDefinedByIdentical4DiskWithFlippingVEdge.clear();
}


int DynamicVoronoiDiagramCIC::compute_polynomial_degree_by_disk_info(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder)
{
    int sumOfSameVelocityVectorCount = 0;

    for (int i = 0; i < 4; i++)
    {
        const rg_Point2D& motionVector1 = fourDiskGeneratorsInNonIncreasingOrder.at(i)->getVelocityVector();
        
        int numOfSameVelocityVectors = 0;

        for (int j = 1; j <= 3; j++)
        {
            const rg_Point2D& motionVector2 = fourDiskGeneratorsInNonIncreasingOrder.at((i + j) % 4)->getVelocityVector();

            //if(motionVector1 == motionVector2)
            if (  rg_EQ(motionVector1.getX(), motionVector2.getX(), TOLERANCE_OF_EQUAL_VELOCITY_VEC)
               && rg_EQ(motionVector1.getY(), motionVector2.getY(), TOLERANCE_OF_EQUAL_VELOCITY_VEC))
            {
                numOfSameVelocityVectors++;
            }
        }

        sumOfSameVelocityVectorCount += numOfSameVelocityVectors;
    }


    int polynomialDegree = 0;

    switch (sumOfSameVelocityVectorCount)
    {
        case 0: //All Different Velocity Vectors
        {
            polynomialDegree = 8;
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_ALL_DIFFERENT_VELOCITY_VECTOR);
            break;
        }

        case 2: //2 Same Velocity Vectors
        {
            polynomialDegree = 6;
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_2_SAME_VELOCITY_VECTOR);
            break;
        }

        case 4:
        {
            polynomialDegree = 4;
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_2_2_SAME_VELOCITY_VECTOR);
            break;
        }

        case 6:
        {
            polynomialDegree = 4;
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_3_SAME_VELOCITY_VECTOR);
            break;
        }

        case 12:
        {
            polynomialDegree = 0;
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_4_SAME_VELOCITY_VECTOR);
            break;
        }
    }

    if (EdgeFlipTimeCalculator::all_four_disks_have_same_radius(fourDiskGeneratorsInNonIncreasingOrder))
    {
        polynomialDegree = polynomialDegree / 2;
    }

    return polynomialDegree;
}





void DynamicVoronoiDiagramCIC::collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(const Generator2D* diskGenerator1, 
                                                                                                               const Generator2D* diskGenerator2,
                                                                                                               set<VEdge2D*>& influencedVEdges)
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskGenerator1->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskGenerator2->getUserData());

    list<VEdge2D*> boundaryOfVCellsOfTwoCollidingDisks;

    if (movingDisk1->getCircle().isIncludedIn(movingDisk2->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        diskGenerator1->getOuterFace()->getBoundaryVEdges(boundaryOfVCellsOfTwoCollidingDisks);
    }
    else if (movingDisk2->getCircle().isIncludedIn(movingDisk1->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        diskGenerator2->getOuterFace()->getBoundaryVEdges(boundaryOfVCellsOfTwoCollidingDisks);
    }
    else // There is no intersection, so this case is two disks
    {
        diskGenerator1->getOuterFace()->getBoundaryVEdges(boundaryOfVCellsOfTwoCollidingDisks);
        diskGenerator2->getOuterFace()->getBoundaryVEdges(boundaryOfVCellsOfTwoCollidingDisks);
    }

    //for erase same vEdges
    boundaryOfVCellsOfTwoCollidingDisks.sort();
    boundaryOfVCellsOfTwoCollidingDisks.unique();

    for (list<VEdge2D*>::iterator it_VEdge = boundaryOfVCellsOfTwoCollidingDisks.begin();
        it_VEdge != boundaryOfVCellsOfTwoCollidingDisks.end();
        it_VEdge++)
    {
        VEdge2D* anVEdge = (*it_VEdge);

        influencedVEdges.insert(anVEdge->getRightHand());
        influencedVEdges.insert(anVEdge->getRightLeg());
        influencedVEdges.insert(anVEdge->getLeftHand());
        influencedVEdges.insert(anVEdge->getLeftLeg());
        influencedVEdges.insert(anVEdge);
    }
}


pair<bool, double> DynamicVoronoiDiagramCIC::find_possible_collision_time_of_disks_after( const double& mostRecentlyEventOccurringTime, const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2 )
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_FIND_POSSIBLE_COLLISION_TIME_OF_DISKS);

    pair<bool, double> nextCollisionTime = DiskCollisionTimeCalculator::find_possible_collision_time_of_disks_after(mostRecentlyEventOccurringTime, movingDisk1, movingDisk2);

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_FIND_POSSIBLE_COLLISION_TIME_OF_DISKS);

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DISK_COLLISION_IN_COMPUTATION);

    return nextCollisionTime;
}


void DynamicVoronoiDiagramCIC::collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_PHYSICAL_COLLISION(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& influencedVEdges)
{
    list<VEdge2D*> boundaryVEdges;

    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskGenerator1->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskGenerator2->getUserData());

    if (movingDisk1->getCircle().isIncludedIn(movingDisk2->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        diskGenerator1->getOuterFace()->getBoundaryVEdges(boundaryVEdges);
    }
    else if (movingDisk2->getCircle().isIncludedIn(movingDisk1->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        diskGenerator2->getOuterFace()->getBoundaryVEdges(boundaryVEdges);
    }
    else // There is no intersection, so this case is two disks
    {
        diskGenerator1->getOuterFace()->getBoundaryVEdges(boundaryVEdges);
        diskGenerator2->getOuterFace()->getBoundaryVEdges(boundaryVEdges);
    }

    influencedVEdges.insert(boundaryVEdges.begin(), boundaryVEdges.end());
}



void DynamicVoronoiDiagramCIC::collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_AVOIDANCE(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& sharedVEdges, set<VEdge2D*>& theOtherVEdges)
{
    list<VEdge2D*> boundaryVEdge1;
    list<VEdge2D*> boundaryVEdge2;

    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskGenerator1->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskGenerator2->getUserData());

    if (get_generator_by_its_input_disk_id(movingDisk1->getID()) != getContainerGenerator())
    {
        diskGenerator1->getOuterFace()->getBoundaryVEdges(boundaryVEdge1);
    }

    if (get_generator_by_its_input_disk_id(movingDisk2->getID()) != getContainerGenerator())
    {
        diskGenerator2->getOuterFace()->getBoundaryVEdges(boundaryVEdge2);
    }

    set<VEdge2D*> setOfBoundaryVEdges1;
    setOfBoundaryVEdges1.insert(boundaryVEdge1.begin(), boundaryVEdge1.end());

    for (list<VEdge2D*>::const_iterator it_VEdge = boundaryVEdge2.begin(); it_VEdge != boundaryVEdge2.end(); ++it_VEdge)
    {
        VEdge2D* currVEdge = (*it_VEdge);
        if (setOfBoundaryVEdges1.find(currVEdge) == setOfBoundaryVEdges1.end())
        {
            theOtherVEdges.insert(currVEdge);
        }
        else
        {
            sharedVEdges.insert(currVEdge);
        }
    }

    for (list<VEdge2D*>::const_iterator it_VEdge = boundaryVEdge1.begin(); it_VEdge != boundaryVEdge1.end(); ++it_VEdge)
    {
        VEdge2D* currVEdge = (*it_VEdge);
        if (sharedVEdges.find(currVEdge) == sharedVEdges.end())
        {
            theOtherVEdges.insert(currVEdge);
        }
    }
}


void DynamicVoronoiDiagramCIC::find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(const double& mostRecentlyEventOccurringTime, set<VEdge2D*>& influencedVEdges)
{
    for (set<VEdge2D*>::iterator it_VEdge = influencedVEdges.begin();
        it_VEdge != influencedVEdges.end();
        it_VEdge++)
    {
        const VEdge2D* VEdge = *it_VEdge;
        
        pair<bool, double> possibleCollisionTime(false, DBL_MAX);

        if (VEdge->isInfinite() || VEdge->isUnBounded()
         || VEdge->getLeftFace()->isUnBounded() || VEdge->getRightFace()->isUnBounded()) // dynamicDisk ���忡���� phantom�� �𸣴� ��

        {
            possibleCollisionTime.first = false;
        }
        else
        {
            DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(VEdge->getLeftFace()->getGenerator()->getUserData());
            DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(VEdge->getRightFace()->getGenerator()->getUserData());

            possibleCollisionTime = find_possible_collision_time_of_disks_after(mostRecentlyEventOccurringTime, movingDisk1, movingDisk2);
        }


#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::DISK_COLLISION, const_cast<VEdge2D*>(VEdge), possibleCollisionTime.second);
#endif //DEBUG_DVD


        do_bookkeeping_on_priorityQ_for_disk_collision(possibleCollisionTime, VEdge, AFTER_PQ_CONSTRUCTION);
    }
}




void DynamicVoronoiDiagramCIC::find_possible_collision_time_of_shared_VEdges_after_input_time_N_store_in_priorityQ_AVOIDANCE(const double& mostRecentlyEventOccurringTime, set<VEdge2D*>& sharedVEdges)
{
    for (set<VEdge2D*>::iterator it_VEdge = sharedVEdges.begin();
         it_VEdge != sharedVEdges.end();
         it_VEdge++)
    {
        const VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleCollisionTime(false, DBL_MAX);

#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::DISK_COLLISION, const_cast<VEdge2D*>(VEdge), possibleCollisionTime.second);
#endif //DEBUG_DVD
        
        do_bookkeeping_on_priorityQ_for_disk_collision(possibleCollisionTime, VEdge, AFTER_PQ_CONSTRUCTION);
    }
}


void DynamicVoronoiDiagramCIC::do_bookkeeping_on_priorityQ_for_disk_collision(const pair<bool, double>& newCollisionTime, const VEdge2D* vEdge, const PriorityQBookkeepingType& PQBookkeepingType)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION);


    bool   doesHaveNextDiskCollision = newCollisionTime.first;
    double nextDiskCollisionTime     = newCollisionTime.second;

    switch (PQBookkeepingType)
    {
    case DURING_INITIAL_PQ_CONSTRUCTION:
    {
        if (doesHaveNextDiskCollision)
        {
            m_PriorityQForDiskCollision.push(const_cast<VEdge2D*>(vEdge), nextDiskCollisionTime);
        }
        else{
            m_PriorityQForDiskCollision.push(const_cast<VEdge2D*>(vEdge), DBL_MAX);
        }


        break;
    }

    case  AFTER_PQ_CONSTRUCTION:
    {
        if (doesHaveNextDiskCollision)
        {
            m_PriorityQForDiskCollision.replaceKey(const_cast<VEdge2D*>(vEdge), nextDiskCollisionTime);
        }
        else
        {
            m_PriorityQForDiskCollision.replaceKey(const_cast<VEdge2D*>(vEdge), DBL_MAX);
        }

        break;
    }
    default:
        break;
    }


    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION);
}


void DynamicVoronoiDiagramCIC::collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_VEdge_flipping_times(const Generator2D* diskGenerator, set<VEdge2D*>& influencedVEdges)
{
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskGenerator->getUserData());

    list<VEdge2D*> boundaryOfVCellOfDisk;

    if (get_generator_by_its_input_disk_id(movingDisk->getID()) != getContainerGenerator())
    {
        diskGenerator->getOuterFace()->getBoundaryVEdges(boundaryOfVCellOfDisk);
    }

    for (list<VEdge2D*>::iterator it_VEdge = boundaryOfVCellOfDisk.begin();
        it_VEdge != boundaryOfVCellOfDisk.end();
        it_VEdge++)
    {
        VEdge2D* anVEdge = (*it_VEdge);

        influencedVEdges.insert(anVEdge->getRightHand());
        influencedVEdges.insert(anVEdge->getRightLeg());
        influencedVEdges.insert(anVEdge->getLeftHand());
        influencedVEdges.insert(anVEdge->getLeftLeg());
        influencedVEdges.insert(anVEdge);
    }
}


void DynamicVoronoiDiagramCIC::collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_collision_times(const Generator2D* diskGenerator, set<VEdge2D*>& influencedVEdges)
{
    list<VEdge2D*> boundaryVEdges;

    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskGenerator->getUserData());

    if (get_generator_by_its_input_disk_id(movingDisk->getID()) != getContainerGenerator())
    {
        diskGenerator->getOuterFace()->getBoundaryVEdges(boundaryVEdges);
    }

    influencedVEdges.insert(boundaryVEdges.begin(), boundaryVEdges.end());
}


void V::GeometryTier::DynamicVoronoiDiagramCIC::handle_and_store_EDGE_FLIP_event(EdgeFlipEvent* edgeFlippingEvent)
{
    if (m_DuringEVENTSEQGeneration || m_Vers1p0EVENTSEQFileUsed)
    {
        set<Generator2D*> influencedDiskGenerators;
        find_influenced_disks_by_this_FLIP_event(influencedDiskGenerators, edgeFlippingEvent);
        move_disks_to_target_Voronoi_clock_without_changing_VD(influencedDiskGenerators, edgeFlippingEvent->get_occurring_clock());
    }

    VEdge2D* targetVEdge = get_VEdge_with_this_generator_quadruplet(edgeFlippingEvent->get_generator_quadruplet_before_flip(), false);
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        flip_Vedge_to_FUTURE(targetVEdge);
    }

    if (m_DuringEVENTSEQGeneration)
    {
        update_the_event_times_in_the_priorityQs_influenced_by_this_FLIP_event(edgeFlippingEvent);

        store_this_FLIP_event_in_the_event_history_vector(edgeFlippingEvent);
    }
}


void V::GeometryTier::DynamicVoronoiDiagramCIC::handle_and_store_DISK_COLLISION_event(DiskCollisionEvent* diskCollisionEvent)
{
    if (m_DuringEVENTSEQGeneration || m_Vers1p0EVENTSEQFileUsed)
    {
        set<Generator2D*> influencedDiskGenerators;
        find_influenced_disks_by_this_DISK_COLLISION_event(influencedDiskGenerators, diskCollisionEvent);
        move_disks_to_target_Voronoi_clock_without_changing_VD(influencedDiskGenerators, diskCollisionEvent->get_occurring_clock());
        
        store_current_positions_of_two_colliing_disks_to_event(diskCollisionEvent);
    }
    else
    {
        change_positions_of_two_colliding_disks_from_event(diskCollisionEvent);
    }

    compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(diskCollisionEvent);

    watch_and_record_collision_to_FUTURE(diskCollisionEvent);

    if (m_DuringEVENTSEQGeneration)
    {
        update_the_event_times_in_the_priorityQs_influenced_by_this_DISK_COLLISION_event(diskCollisionEvent);

        store_this_DISK_COLLISION_event_in_the_event_history_vector(diskCollisionEvent);
    }
}


void V::GeometryTier::DynamicVoronoiDiagramCIC::handle_and_store_VELOCITY_VECTOR_CHANGE_event(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    if (m_DuringEVENTSEQGeneration || m_Vers1p0EVENTSEQFileUsed)
    {
        set<Generator2D*> influencedDiskGenerators;
        find_influenced_disks_by_this_VELOCITY_VECTOR_CHANGE_event(influencedDiskGenerators, diskVelocityChangeEvent);
        move_disks_to_target_Voronoi_clock_without_changing_VD(influencedDiskGenerators, diskVelocityChangeEvent->get_occurring_clock());
        
        store_current_position_of_disk_to_event(diskVelocityChangeEvent);
    }
    else
    {
        change_position_of_disk_from_event(diskVelocityChangeEvent);
    }
    
    change_to_new_velocity_vector_and_store_old_one(diskVelocityChangeEvent);

    if (m_DuringEVENTSEQGeneration)
    {
        update_the_event_times_in_the_priorityQs_influenced_by_this_VELOCITY_VECTOR_CHANGE_event(diskVelocityChangeEvent);

        store_this_VELOCITY_VECTOR_CHANGE_event_in_the_event_history_vector(diskVelocityChangeEvent);
    }
}


void V::GeometryTier::DynamicVoronoiDiagramCIC::handle_and_store_DISK_HOPPING_event(DiskHoppingEvent* diskHoppingEvent)
{
    set<Generator2D*> influencedDiskGenerators;
    find_influenced_disks_by_this_DISK_HOPPING_event(influencedDiskGenerators, diskHoppingEvent);
    //move_disks_in_two_containers_to_the_locations_of_next_event_time_without_changing_VD
    move_disks_to_target_Voronoi_clock_without_changing_VD(influencedDiskGenerators, diskHoppingEvent->get_occurring_clock());

//     if (m_DuringEVENTSEQGeneration || m_Vers1p0EVENTSEQFileUsed)
//     {
//         set<Generator2D*> influencedDiskGenerators;
//         find_influenced_disks_by_this_DISK_HOPPING_event(influencedDiskGenerators, diskHoppingEvent);
//         //move_disks_in_two_containers_to_the_locations_of_next_event_time_without_changing_VD
//         move_disks_to_target_Voronoi_clock_without_changing_VD(influencedDiskGenerators, diskHoppingEvent->get_occurring_clock());
//     }
//     else
//     {
//         change_position_of_disk_from_event(diskHoppingEvent);
//     }


    if (m_DuringEVENTSEQGeneration)
    {
        pair<bool, rg_Point2D> canBeHoppedNTargetPosition
            = check_two_conditions_for_hopping_and_get_inserted_position_if_satisfies_the_conditions(diskHoppingEvent);

        bool thisDiskCanBeHopped = canBeHoppedNTargetPosition.first;
        switch (thisDiskCanBeHopped)
        {
            case true:
            {
				rg_Point2D insertedPosition = canBeHoppedNTargetPosition.second;

				diskHoppingEvent->set_does_this_event_happen(true);

				store_current_and_target_positions_for_disk_hopping(diskHoppingEvent, insertedPosition);
				store_velocity_vectors_before_and_after_disk_hopping_by_current_velocity_vector(diskHoppingEvent);
				relocate_a_disk_btw_containers_by_DISK_HOPPING_to_FUTURE(diskHoppingEvent, true);

				___count_in_GeneratingEVENTSEQ(DVD_COUNT_REAL_OCCURING_DELETE_N_INSERT_DISK_IN_EVENT_HISTORY);
                break;
            }

            case false:
            {
				diskHoppingEvent->set_does_this_event_happen(false);

				store_current_position_of_disk_to_event(diskHoppingEvent);
				compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(diskHoppingEvent);
                break;
            }
        }

        update_the_event_times_in_the_priorityQs_influenced_by_this_DISK_HOPPING_event(diskHoppingEvent);

        store_this_DISK_HOPPING_event_in_the_event_history_vector(diskHoppingEvent);
    }
    else
    {
        Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
        Generator2D* toContainer     = diskHoppingEvent->get_in_what_container_this_generator_is_inserted();

        switch (diskHoppingEvent->get_does_this_event_happen())
        {
            case true:
            {
				store_velocity_vectors_before_and_after_disk_hopping_by_current_velocity_vector(diskHoppingEvent);
				update_circumcircles_of_VVertices_in_this_container(toContainer);
				relocate_a_disk_btw_containers_by_DISK_HOPPING_to_FUTURE(diskHoppingEvent,  false);

				___count_in_GeneratingEVENTSEQ(DVD_COUNT_REAL_OCCURING_DELETE_N_INSERT_DISK_IN_EVENT_HISTORY);
                break;
            }

            case false:
            {
                compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(diskHoppingEvent);
                break;
            }
        }
    }
}


EventOfDynamicVD2D* DynamicVoronoiDiagramCIC::get_next_event_from_priorityQs()
{
    EventOfDynamicVD2D* anEvent = NULL;

    switch (get_next_event_type_in_priorityQ()) 
    {
        case EventOfDynamicVD2D::EDGE_FLIP:
        {
            anEvent = generate_next_edge_flip_event_from_priorityQ();

            ___check_flipping_type_of_VEdge_for_statistics(dynamic_cast<EdgeFlipEvent*>(anEvent));

            break;
        }

        case EventOfDynamicVD2D::DISK_COLLISION:
        {
            if (this_is_collision_event_NOT_disk_hopping_event(m_PriorityQForDiskCollision.topNode()->getEntity()))
            {
                anEvent = generate_next_disk_collision_event_from_priorityQ();
            }
            else
            {
                anEvent = generate_next_disk_hopping_event_from_priorityQ();
            }

            break;
        }

    case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
    {
        anEvent = pop_and_generate_next_disk_velocity_change_event_from_priorityQ();
        
        break;
    }

    default:
        break;
    }

    return anEvent;
}


EdgeFlipEvent* DynamicVoronoiDiagramCIC::generate_next_edge_flip_event_from_priorityQ()
{
    EdgeFlipEvent* edgeFlipEvent = NULL;

    // CASE : SAME FLIPPING TIME IN PQ : 2, 3, 5 // 
    VEdge2D* VEdge = get_next_flipping_edge_from_flipping_priorityQ();
    double   eventTime = m_PriorityQForEdgeFlipping.topNode()->getKey();

    EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = get_generator_quadruplet_by_its_VEdge(VEdge);

    if (this_VEdge_is_NOT_shadow_flipping(VEdge, eventTime, TOLERANCE_OF_SHADOW_FLIP))  //length
    {
        edgeFlipEvent = new EdgeFlipEvent(eventTime, generatorPtrQuadruplet, true);
    }
    else
    {
        //  CASE : SHADOW FLIPPING  //
        edgeFlipEvent = new EdgeFlipEvent(eventTime, generatorPtrQuadruplet, false);
    }

    return edgeFlipEvent;
}


DiskCollisionEvent* DynamicVoronoiDiagramCIC::generate_next_disk_collision_event_from_priorityQ()
{
    DiskCollisionEvent* diskCollisionEvent = NULL;

    VEdge2D* aVEdge = m_PriorityQForDiskCollision.topNode()->getEntity();
    double   eventTime = m_PriorityQForDiskCollision.topNode()->getKey();

    Generator2D* diskGenerator1 = aVEdge->getLeftFace()->getGenerator();
    Generator2D* diskGenerator2 = aVEdge->getRightFace()->getGenerator();

    diskCollisionEvent = new DiskCollisionEvent(eventTime, diskGenerator1, diskGenerator2);

    return diskCollisionEvent;
}


DiskVelocityChangeEvent* DynamicVoronoiDiagramCIC::pop_and_generate_next_disk_velocity_change_event_from_priorityQ()
{
    DiskVelocityChangeEvent* diskVelocityChangeEvent = NULL;

    double eventTime = m_PriorityQForDiskVelocityChange.topNode()->getKey();

    pair<Generator2D*, rg_Point2D> pairOfDiskGeneratorNNewVelocityVector;

    ReservationNumber reservedNum         = m_PriorityQForDiskVelocityChange.pop();
    pairOfDiskGeneratorNNewVelocityVector = m_ReservedVelocityVectorChange.at(reservedNum);
    m_ReservedVelocityVectorChange.erase(reservedNum);

    Generator2D* generatorDisk     = pairOfDiskGeneratorNNewVelocityVector.first;
    rg_Point2D   newVelocityVector = pairOfDiskGeneratorNNewVelocityVector.second;

    diskVelocityChangeEvent = new DiskVelocityChangeEvent(eventTime, generatorDisk, newVelocityVector);

    return diskVelocityChangeEvent;
}



void DynamicVoronoiDiagramCIC::synchronize_disk_locations_and_current_clock_to_VD(const double& targetVoronoiClock)
{
    list<Generator2D*> allGenerator;
    getAllGenerators_hierarchy(allGenerator);

    set<Generator2D*> allGeneratorSet;
    allGeneratorSet.insert(allGenerator.begin(), allGenerator.end());

    move_disks_to_target_Voronoi_clock_without_changing_VD(allGeneratorSet, targetVoronoiClock);

    set_current_Voronoi_clock(targetVoronoiClock);
}


void DynamicVoronoiDiagramCIC::move_disks_to_target_Voronoi_clock_without_changing_VD(set<Generator2D*>& diskGenerators, const double& targetVoronoiClock)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_DISK_MOVING);


    for (set<Generator2D*>::iterator it_Generator = diskGenerators.begin(); it_Generator != diskGenerators.end(); it_Generator++)
    {
        Generator2D* diskGenerator         = *it_Generator;
        DynamicDisk*  dynamicDisk          = static_cast<DynamicDisk*>(diskGenerator->getUserData());
        double       mostRecentUpdateTime  = dynamicDisk->getMostRecentLocationUpdateTime();
        double       timeIncrement         = targetVoronoiClock - mostRecentUpdateTime;

        if (timeIncrement == 0.0)
        {
            continue;
        }

//         float newXCoord = dynamicDisk->getCircle().getCenterPt().getX() + dynamicDisk->getVelocityVectorX()*timeIncrement;
//         float newYCoord = dynamicDisk->getCircle().getCenterPt().getY() + dynamicDisk->getVelocityVectorY()*timeIncrement;

        double newXCoord = dynamicDisk->getCircle().getCenterPt().getX() + dynamicDisk->getVelocityVectorX()*timeIncrement;
        double newYCoord = dynamicDisk->getCircle().getCenterPt().getY() + dynamicDisk->getVelocityVectorY()*timeIncrement;

        rg_Circle2D newCircle(newXCoord, newYCoord, dynamicDisk->getCircle().getRadius());

        diskGenerator->setDisk(newCircle);
        
        dynamicDisk->setCircle(newCircle);
        dynamicDisk->setMostRecentLocationUpdateTime(targetVoronoiClock);
    }


    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_DISK_MOVING);
}



void DynamicVoronoiDiagramCIC::make_link_btw_generators_from_other_VD_to_this_VD(const DynamicVoronoiDiagramCIC& fromThisVD, map<Generator2D*, Generator2D*>& generatorLinkFromOtherVDToThisVD)
{
    list<Generator2D*> fromGenerators;
    list<Generator2D*> toGenerators;

    fromThisVD.getAllGenerators_hierarchy(fromGenerators);
    getAllGenerators_hierarchy(toGenerators);

    list<Generator2D*>::iterator it_from_generator = fromGenerators.begin();
    list<Generator2D*>::iterator it_to_generator   = toGenerators.begin();

    for (; it_from_generator != fromGenerators.end(); it_from_generator++, it_to_generator++)
    {
        pair<Generator2D*, Generator2D*> generatorPair(*it_from_generator, *it_to_generator);

        generatorLinkFromOtherVDToThisVD.insert(generatorPair);
    }
}


void DynamicVoronoiDiagramCIC::make_link_btw_VEdges(list<VEdge2D*>& fromVEdges, list<VEdge2D*>& toVEdges, map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD)
{
    list<VEdge2D*>::iterator it_from_VEdges = fromVEdges.begin();
    list<VEdge2D*>::iterator it_to_VEdges   = toVEdges.begin();

    for (; it_from_VEdges != fromVEdges.end(); it_from_VEdges++, it_to_VEdges++)
    {
        pair<VEdge2D*, VEdge2D*> VEdgePair(*it_from_VEdges, *it_to_VEdges);

        VEdgeLinkFromOtherVDToThisVD.insert(VEdgePair);
    }
}


void DynamicVoronoiDiagramCIC::make_link_btw_disks_and_containers_from_other_VD_to_this_VD(const DynamicVoronoiDiagramCIC& fromThisVD, map<DynamicDisk*, DynamicDisk*>& diskLinkFromOtherVDToThisVD)
{
    list<DynamicDisk*> fromDisks;
    list<DynamicDisk*> toDisks;

    fromThisVD.get_dynamic_disks_and_containers(fromDisks);
    get_dynamic_disks_and_containers(toDisks);

    list<DynamicDisk*>::iterator it_fromDisks = fromDisks.begin();
    list<DynamicDisk*>::iterator it_toDisks   = toDisks.begin();

    for (; it_fromDisks != fromDisks.end(); it_fromDisks++, it_toDisks++)
    {
        pair<DynamicDisk*, DynamicDisk*> diskPair(*it_fromDisks, *it_toDisks);

        diskLinkFromOtherVDToThisVD.insert(diskPair);
    }
}



void DynamicVoronoiDiagramCIC::modify_user_data_pointing_to_disk_in_generators(map<DynamicDisk*, DynamicDisk*>& diskLinkFromOtherVDToThisVD)
{
    list<Generator2D*> allGeneratorsFromThisVD;
    getAllGenerators_hierarchy(allGeneratorsFromThisVD);

    for (list<Generator2D*>::iterator it_Generator = allGeneratorsFromThisVD.begin();
        it_Generator != allGeneratorsFromThisVD.end();
        it_Generator++)
    {
        Generator2D*  currGenerator        = *it_Generator;
        DynamicDisk*  movingDiskFromThisVD = diskLinkFromOtherVDToThisVD.at(static_cast<DynamicDisk*>(currGenerator->getUserData()));

        currGenerator->setUserData(movingDiskFromThisVD);
    }
}



void DynamicVoronoiDiagramCIC::modify_generators_in_event_sequence(const map<Generator2D*, Generator2D*>& generatorLinkFromOtherVDToThisVD)
{
    const int eventSequenceSize = m_EventSequence.size();

    for (int i = 0; i < eventSequenceSize; i++)
    {
        EventOfDynamicVD2D* event = m_EventSequence.at(i);

        switch (event->get_event_type())
        {
            case EventOfDynamicVD2D::EDGE_FLIP:
            {
                EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(event);

                EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = edgeFlippingEvent->get_generator_quadruplet_before_flip();

                Generator2D* leftDiskGeneratorInThisVD  = generatorLinkFromOtherVDToThisVD.at(generatorPtrQuadruplet[0]);
                Generator2D* rightDiskGeneratorInThisVD = generatorLinkFromOtherVDToThisVD.at(generatorPtrQuadruplet[1]);
                Generator2D* headDiskGeneratorInThisVD  = generatorLinkFromOtherVDToThisVD.at(generatorPtrQuadruplet[2]);
                Generator2D* tailDiskGeneratorInThisVD  = generatorLinkFromOtherVDToThisVD.at(generatorPtrQuadruplet[3]);

                EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadrupletOfThisVD =
                         { leftDiskGeneratorInThisVD, rightDiskGeneratorInThisVD, headDiskGeneratorInThisVD, tailDiskGeneratorInThisVD };

                edgeFlippingEvent->set_generator_quadruplet_before_flip(generatorPtrQuadrupletOfThisVD);

                break;
            }

            case EventOfDynamicVD2D::DISK_COLLISION:
            {
                DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(event);
            
                Generator2D* diskGenerator1InThisVD = generatorLinkFromOtherVDToThisVD.at(diskCollisionEvent->get_disk_generator1());
                Generator2D* diskGenerator2InThisVD = generatorLinkFromOtherVDToThisVD.at(diskCollisionEvent->get_disk_generator2());

                diskCollisionEvent->set_two_collided_generators(diskGenerator1InThisVD, diskGenerator2InThisVD);

                break;
            }

            case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
            {
                DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(event);

                Generator2D* diskGeneratorInThisVD = generatorLinkFromOtherVDToThisVD.at(diskVelocityChangeEvent->get_disk_generator());

                diskVelocityChangeEvent->set_velocity_changed_generator(diskGeneratorInThisVD);

                break;
            }

            case EventOfDynamicVD2D::DISK_HOPPING:
            {
                DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(event);

                Generator2D* targetGeneratorInThisVD = generatorLinkFromOtherVDToThisVD.at(diskHoppingEvent->get_target_generator());
                Generator2D* fromContainerInThisVD   = generatorLinkFromOtherVDToThisVD.at(diskHoppingEvent->get_in_what_container_this_generator_is_deleted());
                Generator2D* toContainerInThisVD     = generatorLinkFromOtherVDToThisVD.at(diskHoppingEvent->get_in_what_container_this_generator_is_inserted());

                diskHoppingEvent->set_target_generator(targetGeneratorInThisVD);
                diskHoppingEvent->set_in_what_container_this_generator_is_deleted_n_inserted(fromContainerInThisVD, toContainerInThisVD);
 
                break;
            }
        }
    }
}



void DynamicVoronoiDiagramCIC::modify_VEdges_in_priorityQs_from_other_VDs_one_to_this_VDs_one(const map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD)
{ 
    if (m_PriorityQForEdgeFlipping.size() != 0)
    {
        for (map<VEdge2D*, VEdge2D*>::const_iterator it_pair = VEdgeLinkFromOtherVDToThisVD.begin();
            it_pair != VEdgeLinkFromOtherVDToThisVD.end();
            it_pair++)
        {
            VEdge2D* fromVEdge = (*it_pair).first;
            VEdge2D* toVEdge = (*it_pair).second;

            m_PriorityQForEdgeFlipping.replaceEntity(fromVEdge, toVEdge);
        }
    }

    if (m_PriorityQForDiskCollision.size() != 0)
    {
        for (map<VEdge2D*, VEdge2D*>::const_iterator it_pair = VEdgeLinkFromOtherVDToThisVD.begin();
            it_pair != VEdgeLinkFromOtherVDToThisVD.end();
            it_pair++)
        {
            VEdge2D* fromVEdge = (*it_pair).first;
            VEdge2D* toVEdge = (*it_pair).second;

            m_PriorityQForDiskCollision.replaceEntity(fromVEdge, toVEdge);
        }
    }
}


void DynamicVoronoiDiagramCIC::rewind_EVENTSEQ(const double& targetVoronoiClock)
{
    construct_VD_at_target_Voronoi_clock(targetVoronoiClock);

    remove_events_in_EVENTSEQ_after_current_Voronoi_clock();

    set_length_of_Voronoi_clock(get_current_Voronoi_clock());

    if (m_InitialPQIsConstructed)
    {
        //This part can be improved.
        clear_priority_queues_for_collision_and_flip();
    }
}


void DynamicVoronoiDiagramCIC::set_EVENTSEQ(const vector<EventOfDynamicVD2D*>& eventSequence, const bool& reinstateInitialVD/*= true*/)
{
    if (!m_EventSequence.empty())
    {
        for (vector<EventOfDynamicVD2D*>::const_iterator it_Event = m_EventSequence.begin();
            it_Event != m_EventSequence.end();
            it_Event++)
        {
            delete (*it_Event);
        }

        m_EventSequence.clear();
    }

    m_EventSequence = eventSequence;

    compute_velocity_vectors_by_scanning_EVENTSEQ();

    if (reinstateInitialVD)
    {
        reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
    }
}


void DynamicVoronoiDiagramCIC::get_dynamic_disks_and_containers(list<DynamicDisk*>& dynamicDisksNContainers) const
{
    list<Generator2D*> allGenerators;
    getAllGenerators_hierarchy(allGenerators);

    for (list<Generator2D*>::iterator it_Generator = allGenerators.begin();
        it_Generator != allGenerators.end();
        it_Generator++)
    {
        Generator2D* currGeneratorDisk = (*it_Generator);
        DynamicDisk* currDynamicDisk = static_cast<DynamicDisk*>(currGeneratorDisk->getUserData());

        dynamicDisksNContainers.push_back(currDynamicDisk);
    }
}


void DynamicVoronoiDiagramCIC::get_dynamic_disks(list<DynamicDisk*>& dynamicDisks) const
{
    list<Generator2D*> allGenerators;
    getAllGenerators_hierarchy(allGenerators);
   
    for (list<Generator2D*>::iterator it_Generator = allGenerators.begin();
        it_Generator != allGenerators.end();
        it_Generator++)
    {
        Generator2D* currGeneratorDisk = (*it_Generator);
        DynamicDisk* currDynamicDisk   = static_cast<DynamicDisk*>(currGeneratorDisk->getUserData());

        if (!currDynamicDisk->this_disk_is_container())
        {
            dynamicDisks.push_back(currDynamicDisk);
        }
    }
}


void DynamicVoronoiDiagramCIC::get_containers(list<DynamicDisk*>& containers) const
{
    list<Generator2D*> allGenerators;
    getAllGenerators_hierarchy(allGenerators);

    for (list<Generator2D*>::iterator it_Generator = allGenerators.begin();
        it_Generator != allGenerators.end();
        it_Generator++)
    {
        Generator2D* currGeneratorDisk = (*it_Generator);
        DynamicDisk* currDynamicDisk = static_cast<DynamicDisk*>(currGeneratorDisk->getUserData());

        if (currDynamicDisk->this_disk_is_container())
        {
            containers.push_back(currDynamicDisk);
        }
    }
}


DynamicDisk* DynamicVoronoiDiagramCIC::get_most_outer_container() const
{
    const Generator2D* containerGenerator = getContainerGenerator();

    DynamicDisk* containerDisk = static_cast<DynamicDisk*>(containerGenerator->getUserData());

    return containerDisk;
}


void DynamicVoronoiDiagramCIC::progress_DVD_using_EVENTNUM_to_FUTURE(const int& eventNumIncrement)
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
        targetVoronoiClock = get_most_recent_event()->get_occurring_clock();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}

 
void DynamicVoronoiDiagramCIC::progress_DVD_using_EVENTNUM_to_PAST(const int& eventNumIncrement)
{
    int absEventNumIncrement = abs(eventNumIncrement);

    const double mostRecentEventTime = get_most_recent_event()->get_occurring_clock();
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
        targetVoronoiClock = get_most_recent_event()->get_occurring_clock();
    }

	synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVoronoiDiagramCIC::progress_DVD_to_FUTURE(const double& targetVoronoiClock)
{
    while (next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(targetVoronoiClock))
    {
        handle_next_DVD_event_to_FUTURE();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}


void DynamicVoronoiDiagramCIC::progress_DVD_to_PAST(const double& targetVoronoiClock)
{
    while (next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(targetVoronoiClock) )
        //|| next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock())
    {
        reversely_handle_next_DVD_event_to_PAST();
    }

    synchronize_disk_locations_and_current_clock_to_VD(targetVoronoiClock);
}



void DynamicVoronoiDiagramCIC::handle_next_DVD_event_to_FUTURE()
{
    EventOfDynamicVD2D* nextFutureEvent = get_next_DVD_event_to_FUTURE();

    switch (nextFutureEvent->get_event_type())
    {
        case EventOfDynamicVD2D::EDGE_FLIP: //flip this VEdge
        {
            EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(nextFutureEvent);
            update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_FUTURE(edgeFlippingEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_COLLISION: //move disks and change motion vectors 
        {
            DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(nextFutureEvent);
            relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_FUTURE(diskCollisionEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
        {
            DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(nextFutureEvent);
            relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_FUTURE(diskVelocityChangeEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_HOPPING:
        {
            DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(nextFutureEvent);
            relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_FUTURE(diskHoppingEvent);
            break;
        }
    }

    set_current_Voronoi_clock(nextFutureEvent->get_occurring_clock());

    increment_index_of_most_recent_event();
}


void DynamicVoronoiDiagramCIC::reversely_handle_next_DVD_event_to_PAST()
{
    EventOfDynamicVD2D* nextPastEvent = get_next_DVD_event_to_PAST();

    switch (nextPastEvent->get_event_type())
    {
        case EventOfDynamicVD2D::EDGE_FLIP: //flip this VEdge
        {
            EdgeFlipEvent* edgeFlippingEvent = dynamic_cast<EdgeFlipEvent*>(nextPastEvent);
            update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_PAST(edgeFlippingEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_COLLISION: //move disks and change motion vectors 
        {
            DiskCollisionEvent* diskCollisionEvent = dynamic_cast<DiskCollisionEvent*>(nextPastEvent);
            relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_PAST(diskCollisionEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
        {
            DiskVelocityChangeEvent* diskVelocityChangeEvent = dynamic_cast<DiskVelocityChangeEvent*>(nextPastEvent);
            relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_PAST(diskVelocityChangeEvent);
            break;
        }

        case EventOfDynamicVD2D::DISK_HOPPING:
        {
            DiskHoppingEvent* diskHoppingEvent = dynamic_cast<DiskHoppingEvent*>(nextPastEvent);
            relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_PAST(diskHoppingEvent);
            break;
        }
    }

    set_current_Voronoi_clock(nextPastEvent->get_occurring_clock());

    decrement_index_of_most_recent_event();
}


void DynamicVoronoiDiagramCIC::update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_FUTURE(EdgeFlipEvent* edgeFlippingEvent)
{  
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        VEdge2D* targetVEdge = get_VEdge_with_this_generator_quadruplet(edgeFlippingEvent->get_generator_quadruplet_before_flip(), false);
        flip_Vedge_to_FUTURE(targetVEdge);
        ___flipping_count_in_phase2();
    }
}



void DynamicVoronoiDiagramCIC::update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_PAST(EdgeFlipEvent* edgeFlippingEvent)
{
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        VEdge2D* targetVEdge = get_VEdge_with_this_generator_quadruplet(edgeFlippingEvent->get_generator_quadruplet_before_flip(), true);
        flip_Vedge_to_PAST(targetVEdge); //[2020.03.17]
    }
}



void DynamicVoronoiDiagramCIC::relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_FUTURE(DiskCollisionEvent* diskCollisionEvent)
{
    rg_Point2D newVelocityVector1 = diskCollisionEvent->get_future_velocity_vector_of_disk1();
    rg_Point2D newVelocityVector2 = diskCollisionEvent->get_future_velocity_vector_of_disk2();

    relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION(diskCollisionEvent, newVelocityVector1, newVelocityVector2);

    DynamicDisk* dynamicDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
	DynamicDisk* dynamicDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    watch_and_record_collision_to_FUTURE(diskCollisionEvent);

    if (!dynamicDisk1->this_disk_is_container() && !dynamicDisk2->this_disk_is_container())
    {
        ___collision_btw_disks_count_in_phase2();
    }

    ___collision_count_in_phase2();
}


void DynamicVoronoiDiagramCIC::relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION(DiskCollisionEvent* diskCollisionEvent, const rg_Point2D& newVelocityVector1, const rg_Point2D& newVelocityVector2)
{
    Generator2D* diskGenerator1 = diskCollisionEvent->get_disk_generator1();
    Generator2D* diskGenerator2 = diskCollisionEvent->get_disk_generator2();

    const double currEventTime = diskCollisionEvent->get_occurring_clock();

    const rg_Point2D& disk1CenterAtCollidingTime = diskCollisionEvent->get_center_of_disk1();
    const rg_Point2D& disk2CenterAtCollidingTime = diskCollisionEvent->get_center_of_disk2();

    relocate_two_disks_at_colliding_time(currEventTime, diskGenerator1, diskGenerator2, disk1CenterAtCollidingTime, disk2CenterAtCollidingTime);

    update_two_disks_velocity_vectors_by_collision(diskGenerator1, diskGenerator2, newVelocityVector1, newVelocityVector2);
}


void DynamicVoronoiDiagramCIC::relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_PAST(DiskCollisionEvent* diskCollisionEvent)
{
    rg_Point2D newVelocityVector1 = diskCollisionEvent->get_past_velocity_vector_of_disk1();
    rg_Point2D newVelocityVector2 = diskCollisionEvent->get_past_velocity_vector_of_disk2();

    relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION(diskCollisionEvent, newVelocityVector1, newVelocityVector2);

    watch_and_record_collision_to_PAST(diskCollisionEvent);
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_FUTURE(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    rg_Point2D newVelocityVector = diskVelocityChangeEvent->get_future_velocity_vector_of_disk();

    relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE(diskVelocityChangeEvent, newVelocityVector);

    ___velocity_change_count_in_phase2();
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE(DiskVelocityChangeEvent* diskVelocityChangeEvent, const rg_Point2D& newVelocityVector)
{
    Generator2D* diskGenerator                  = diskVelocityChangeEvent->get_disk_generator();
    const rg_Point2D& diskCenterAtCollidingTime = diskVelocityChangeEvent->get_center_of_disk();
    const double currEventTime                  = diskVelocityChangeEvent->get_occurring_clock();

    relocate_a_disk_at_velocity_vector_change_time(currEventTime, diskGenerator, diskCenterAtCollidingTime);

    update_a_disk_velocity_vector(diskGenerator, newVelocityVector);
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_PAST(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    rg_Point2D newVelocityVector = diskVelocityChangeEvent->get_past_velocity_vector_of_disk();
  
    relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE(diskVelocityChangeEvent, newVelocityVector);
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_FUTURE(DiskHoppingEvent* diskHoppingEvent)
{
    if (diskHoppingEvent->get_does_this_event_happen())
    {
        move_disks_and_update_VVertices_related_to_disk_hopping_event_to_FUTURE(diskHoppingEvent);
        relocate_a_disk_btw_containers_by_DISK_HOPPING_to_FUTURE(diskHoppingEvent, false);
        ___disk_hopping_in_phase2();
    }
    else //handle like velocity change
	{
        rg_Point2D newVelocityVector = diskHoppingEvent->get_future_velocity_vector_of_disk();

        relocate_a_disk_and_update_velocity_vector_by_FAKE_DISK_HOPPING(diskHoppingEvent, newVelocityVector);

        ___collision_count_in_phase2();
    }
}

void DynamicVoronoiDiagramCIC::relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_PAST(DiskHoppingEvent* diskHoppingEvent)
{
    if (diskHoppingEvent->get_does_this_event_happen())
    {
        move_disks_and_update_VVertices_related_to_disk_hopping_event_to_PAST(diskHoppingEvent);
        relocate_a_disk_btw_containers_by_DISK_HOPPING_to_PAST(diskHoppingEvent, false);
    }
    else //handle like velocity change
	{
        rg_Point2D  newVelocityVector = diskHoppingEvent->get_past_velocity_vector_of_disk();

        relocate_a_disk_and_update_velocity_vector_by_FAKE_DISK_HOPPING(diskHoppingEvent, newVelocityVector);
    }
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_btw_containers_by_DISK_HOPPING(
    DiskHoppingEvent* diskHoppingEvent, Generator2D* toContainer, Generator2D* fromContainer, const rg_Point2D& newPosition, const bool& duringComputation /*= true*/)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_DISK_HOPPING);

    Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
    DynamicDisk* targetDisk      = static_cast<DynamicDisk*>(targetGenerator->getUserData());

    list<VEdge2D*> newVEdges;
    list<VEdge2D*> deletedVEdges;
    list<VEdge2D*> redefinedVEdges;

    fromContainer->getInnerVD()->updateVoronoiDiagram_removal_deleteFromGeneratorList(
        targetGenerator, newVEdges, deletedVEdges, redefinedVEdges);
    fromContainer->removeOneInnerGen(targetGenerator);

    rg_Circle2D newDisk(newPosition, targetDisk->getCircle().getRadius());
    targetDisk->setCircle(newDisk);
    targetGenerator->setDisk(newDisk);

    toContainer->getInnerVD()->updateVoronoiDiagram_insertion_addToGeneratorList(
        targetGenerator, newVEdges, deletedVEdges, redefinedVEdges);
    toContainer->addOneInnerGen(targetGenerator);

    if (duringComputation)
    {
        m_NewVEdgesByDiskHopping.insert(newVEdges.begin(), newVEdges.end());
        m_DeletedVEdgesByDiskHopping.insert(deletedVEdges.begin(), deletedVEdges.end());
        m_RedefinedVEdgesByDiskHopping.insert(redefinedVEdges.begin(), redefinedVEdges.end());
    }

    clear_delete_RED_VEdges_and_VVertices();

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_DISK_HOPPING);
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_btw_containers_by_DISK_HOPPING_to_FUTURE(DiskHoppingEvent* diskHoppingEvent, const bool& duringComputation /*= true*/)
{
    relocate_a_disk_btw_containers_by_DISK_HOPPING(diskHoppingEvent,
                                                   diskHoppingEvent->get_in_what_container_this_generator_is_inserted(), 
                                                   diskHoppingEvent->get_in_what_container_this_generator_is_deleted(),
                                                   diskHoppingEvent->get_inserted_position(),
                                                   duringComputation);
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_btw_containers_by_DISK_HOPPING_to_PAST(DiskHoppingEvent* diskHoppingEvent, const bool& duringComputation /*= true*/)
{
    relocate_a_disk_btw_containers_by_DISK_HOPPING(diskHoppingEvent,
                                                   diskHoppingEvent->get_in_what_container_this_generator_is_deleted(), 
                                                   diskHoppingEvent->get_in_what_container_this_generator_is_inserted(),
                                                   diskHoppingEvent->get_deleted_position(),
                                                   duringComputation);
}


void DynamicVoronoiDiagramCIC::relocate_a_disk_and_update_velocity_vector_by_FAKE_DISK_HOPPING(DiskHoppingEvent* diskHoppingEvent, const rg_Point2D& newVelocityVector)
{
	Generator2D* targetGenerator        = diskHoppingEvent->get_target_generator();
	const double currEventTime          = diskHoppingEvent->get_occurring_clock();
    const rg_Point2D& deletedPosition   = diskHoppingEvent->get_deleted_position();

    relocate_a_disk_at_velocity_vector_change_time(currEventTime, targetGenerator, deletedPosition);
    update_a_disk_velocity_vector(targetGenerator, newVelocityVector);
}


void DynamicVoronoiDiagramCIC::do_bookkeeping_for_disk_hopping(
    const double& mostRecentlyEventOccurringTime,
    const VEdge2D* vEdgeForDiskNContainer,
    const DynamicDisk* movingDisk,
    const DynamicDisk* container,
    const DVD_VEdgeStatus& VDEdgeStatus,
    const bool& thereExistNextCollisionTime)
{
    switch (VDEdgeStatus)
    {
        case NEW_EDGE:
        {
            if (thereExistNextCollisionTime)
            {
                if (movingDisk->getCircle().isIncludedIn(container->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
                {
                    Generator2D* containerToWhichDiskIsInserted = choose_next_container_by_transition_probability(container->getID());
                    DynamicDisk* dynamicContainerToWhichDiskIsInserted = static_cast<DynamicDisk*>((containerToWhichDiskIsInserted)->getUserData());

                    if (dynamicContainerToWhichDiskIsInserted != container)
                    {
                        m_InThisContainerDiskWillBeInserted.insert(make_pair(const_cast<VEdge2D*>(vEdgeForDiskNContainer), containerToWhichDiskIsInserted));

                        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DELETE_N_INSERT_DISK_IN_COMPUTATION);
                    }
                }
            }

            break;
        }

        case DELETED_EDGE:
        {
            if (!this_is_collision_event_NOT_disk_hopping_event(vEdgeForDiskNContainer))
            {
                m_InThisContainerDiskWillBeInserted.erase(const_cast<VEdge2D*>(vEdgeForDiskNContainer));
            }

            break;
        }

        case REDEFINED_EDGE:
        {
            //flip�Ǽ� container�� ���ֺ��� ���?: 
            
            //vedge�� �þ�ų� �پ��? ���?
            if (movingDisk->getCircle().isIncludedIn(container->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
            {
                if (this_is_collision_event_NOT_disk_hopping_event(vEdgeForDiskNContainer) && thereExistNextCollisionTime)
                {
                    Generator2D* containerToWhichDiskIsInserted = choose_next_container_by_transition_probability(container->getID());
                    DynamicDisk* dynamicContainerToWhichDiskIsInserted = static_cast<DynamicDisk*>((containerToWhichDiskIsInserted)->getUserData());

                    if (dynamicContainerToWhichDiskIsInserted != container)
                    {
                        m_InThisContainerDiskWillBeInserted.insert(make_pair(const_cast<VEdge2D*>(vEdgeForDiskNContainer), containerToWhichDiskIsInserted));
                    
                        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DELETE_N_INSERT_DISK_IN_COMPUTATION);
                    }
                }
                else if (!this_is_collision_event_NOT_disk_hopping_event(vEdgeForDiskNContainer) && !thereExistNextCollisionTime)
                {
                    m_InThisContainerDiskWillBeInserted.erase(const_cast<VEdge2D*>(vEdgeForDiskNContainer));
                }
                else if (!this_is_collision_event_NOT_disk_hopping_event(vEdgeForDiskNContainer) && thereExistNextCollisionTime)
                {
                    //do nothing
                }
                else //if (this_event_is_collision_not_deleteNInsert(vEdgeForDiskNContainer) && !thereExistNextCollisionTime)
                {
                    //do nothing
                }
            }
            else
            {
                if (!this_is_collision_event_NOT_disk_hopping_event(vEdgeForDiskNContainer))
                {
                    m_InThisContainerDiskWillBeInserted.erase(const_cast<VEdge2D*>(vEdgeForDiskNContainer));
                }
            }

            break;
        }

        case NONE_BUT_DISK_VEl_VEC_CHANGED: 
        {
            //deleteNInsert�� �ִٸ� update ���ٸ� insert
            if (movingDisk->getCircle().isIncludedIn(container->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
            {
                if (!thereExistNextCollisionTime)
                {
                    if (!this_is_collision_event_NOT_disk_hopping_event(vEdgeForDiskNContainer))
                    {
                        m_InThisContainerDiskWillBeInserted.erase(const_cast<VEdge2D*>(vEdgeForDiskNContainer));
                    }
                }
            }

            break;
        }
    }
}


void DynamicVoronoiDiagramCIC::construct_initial_priority_queues()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS);

#ifndef PQ_STORING_ALL_VEDGES
    construct_initial_priorityQ_of_disk_collision_time_and_mark_hopping_disk_among_them();

    construct_initial_priorityQ_of_vedge_flipping_time_vers2();
#else
    construct_initial_priorityQ_of_disk_collision_time();

    construct_initial_priorityQ_of_vedge_flipping_time();
#endif

    m_InitialPQIsConstructed = true;

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS);
}


void DynamicVoronoiDiagramCIC::construct_initial_priorityQ_of_vedge_flipping_time_vers2()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME);

    list<VEdge2D*> allVEdges;
    if (m_NumOfContainers > 1)
    {
        getAllEdges_hierarchy(allVEdges);
    }
    else
    {
        getVoronoiEdges(allVEdges);
    }

    for (list<VEdge2D*>::const_iterator it_VEdge = allVEdges.begin();
        it_VEdge != allVEdges.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleEdgeFlipingTime(false, DBL_MAX);

        if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
        {
            possibleEdgeFlipingTime.first = false;
        }
        else
        {
            possibleEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(m_CurrentVoronoiClock, VEdge);
        }

        do_bookkeeping_on_priorityQ_for_edge_flipping(possibleEdgeFlipingTime, VEdge, NEW_EDGE);
    }

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME);
}


void DynamicVoronoiDiagramCIC::construct_initial_priorityQ_of_disk_collision_time_and_mark_hopping_disk_among_them()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME);

    list<VEdge2D*> allVEdges;
    getAllEdges_hierarchy(allVEdges);

    for (list<VEdge2D*>::const_iterator it_VEdge = allVEdges.begin();
        it_VEdge != allVEdges.end();
        it_VEdge++)
    {
        VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleDiskCollisionTime(false, DBL_MAX);
        DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(VEdge->getLeftFace()->getGenerator()->getUserData());
        DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(VEdge->getRightFace()->getGenerator()->getUserData());

        possibleDiskCollisionTime = find_possible_collision_time_of_disks_after(m_CurrentVoronoiClock, movingDisk1, movingDisk2);
        do_bookkeeping_on_priorityQ_for_disk_collision(possibleDiskCollisionTime, VEdge, NEW_EDGE);

        if (m_NumOfContainers > 1)
        {
            DVD_StatusOfTwoDisks stateOfTwoDisks = find_state_of_two_disks(movingDisk1, movingDisk2);

            switch (stateOfTwoDisks)
            {
                case DISK_N_CONTAINER:
                {
                    do_bookkeeping_for_disk_hopping(
                        m_CurrentVoronoiClock, VEdge, movingDisk1, movingDisk2, NEW_EDGE, possibleDiskCollisionTime.first);
                    break;
                }

                case CONTAINER_N_DISK:
                {
                    do_bookkeeping_for_disk_hopping(
                        m_CurrentVoronoiClock, VEdge, movingDisk2, movingDisk1, NEW_EDGE, possibleDiskCollisionTime.first);
                    break;
                }
            }
        }
    }

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME);
}


void DynamicVoronoiDiagramCIC::reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ()
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_REINSTATE_THE_INITIAL_VORONOI_DIAGRAM);

    reinstate_the_initial_Voronoi_diagram();

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_REINSTATE_THE_INITIAL_VORONOI_DIAGRAM);
}


void DynamicVoronoiDiagramCIC::do_bookkeeping_on_priorityQ_for_disk_collision(
    const pair<bool, double>& newCollisionTime,
    const VEdge2D* vEdge, 
    const DVD_VEdgeStatus& VEdgeStatus)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION);

    switch (VEdgeStatus)
    {
        case NEW_EDGE:
        {
            if (newCollisionTime.first)
            {
                m_PriorityQForDiskCollision.push(const_cast<VEdge2D*>(vEdge), newCollisionTime.second);
            }
            break;
        }

        case DELETED_EDGE:
        {
            m_PriorityQForDiskCollision.killNodeRelatedWith(const_cast<VEdge2D*>(vEdge));
            break;
        }

        case REDEFINED_EDGE:
        case NONE_BUT_DISK_VEl_VEC_CHANGED:
        {
            if (newCollisionTime.first)
            {
                if (m_PriorityQForDiskCollision.doesHave(const_cast<VEdge2D*>(vEdge)))
                {
                    m_PriorityQForDiskCollision.replaceKey(const_cast<VEdge2D*>(vEdge), newCollisionTime.second);
                }
                else
                {
                    m_PriorityQForDiskCollision.push(const_cast<VEdge2D*>(vEdge), newCollisionTime.second);
                }
            }
            else
            {
                if (m_PriorityQForDiskCollision.doesHave(const_cast<VEdge2D*>(vEdge)))
                {
                    m_PriorityQForDiskCollision.killNodeRelatedWith(const_cast<VEdge2D*>(vEdge));
                }
                else
                {
                    //do nothing
                }
            }

            break;
        }
    }


    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION);
}

void DynamicVoronoiDiagramCIC::do_bookkeeping_on_priorityQ_for_edge_flipping(const pair<bool, double>& newCollisionTime, const VEdge2D * vEdge, const DVD_VEdgeStatus & DVD_VEdgeStatus)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING);

    switch (DVD_VEdgeStatus)
    {
        case NEW_EDGE:
        {
            if (newCollisionTime.first)
            {
                m_PriorityQForEdgeFlipping.push(const_cast<VEdge2D*>(vEdge), newCollisionTime.second);
            }
            break;
        }

        case DELETED_EDGE:
        {
            m_PriorityQForEdgeFlipping.killNodeRelatedWith(const_cast<VEdge2D*>(vEdge));
            break;
        }

        case REDEFINED_EDGE:
        case NONE_BUT_DISK_VEl_VEC_CHANGED:
        {
            if (newCollisionTime.first)
            {
                if (m_PriorityQForEdgeFlipping.doesHave(const_cast<VEdge2D*>(vEdge)))
                {
                    m_PriorityQForEdgeFlipping.replaceKey(const_cast<VEdge2D*>(vEdge), newCollisionTime.second);
                }
                else
                {
                    m_PriorityQForEdgeFlipping.push(const_cast<VEdge2D*>(vEdge), newCollisionTime.second);
                }
            }
            else
            {
                if (m_PriorityQForEdgeFlipping.doesHave(const_cast<VEdge2D*>(vEdge)))
                {
                    m_PriorityQForEdgeFlipping.killNodeRelatedWith(const_cast<VEdge2D*>(vEdge));
                }
                else
                {
                    //do nothing
                }  
            }

            break;
        }
    }

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING);
}

void DynamicVoronoiDiagramCIC::compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(DiskHoppingEvent * diskHoppingEvent) const
{
    DynamicDisk* targetDisk = static_cast<DynamicDisk*>(diskHoppingEvent->get_target_generator()->getUserData());
    DynamicDisk* container  = static_cast<DynamicDisk*>(diskHoppingEvent->get_in_what_container_this_generator_is_deleted()->getUserData());

    rg_Point2D newVelocityVector1, newVelocityVector2;

    compute_velocity_vectors_of_disks_after_collided(targetDisk, container, newVelocityVector1, newVelocityVector2);

    store_velocity_vector_before_and_after_collision(diskHoppingEvent, newVelocityVector1);

    update_velocity_vector(targetDisk, newVelocityVector1);
}



DiskHoppingEvent* DynamicVoronoiDiagramCIC::generate_next_disk_hopping_event_from_priorityQ() const
{
    DiskHoppingEvent* diskHoppingEvent = NULL;

    VEdge2D*  VEdgeOfDiskNContainer  = m_PriorityQForDiskCollision.topNode()->getEntity();
    double           eventTime       = m_PriorityQForDiskCollision.topNode()->getKey();

    Generator2D* disk1 = VEdgeOfDiskNContainer->getRightFace()->getGenerator();
    Generator2D* disk2 = VEdgeOfDiskNContainer->getLeftFace()->getGenerator();

    Generator2D* disk      = NULL;
    Generator2D* container = NULL;

    if (disk1->getDisk().getRadius() < disk2->getDisk().getRadius())
    {
        disk      = disk1;
        container = disk2;
    }
    else
    {
        disk      = disk2;
        container = disk1;
    }

    Generator2D* containerInWhichDiskWillBeInserted = m_InThisContainerDiskWillBeInserted.at(VEdgeOfDiskNContainer);

    diskHoppingEvent = new DiskHoppingEvent(eventTime, 
                                                    disk, 
                                                    container,
                                                    containerInWhichDiskWillBeInserted);

    return diskHoppingEvent;
}



pair<bool, rg_Point2D> DynamicVoronoiDiagramCIC::find_insert_position_of_disk_in_this_container_using_VVertices(Generator2D* movingDisk, Generator2D * container)
{
    bool       thereIsPositionToInsert = false;
    rg_Point2D insertPosition(DBL_MAX, DBL_MAX);

    list<VVertex2D*> VVerticesInThisContainer;
    container->getInnerVD()->getVoronoiVertices(VVerticesInThisContainer);

    for (list<VVertex2D*>::const_iterator it_Vertex = VVerticesInThisContainer.begin();
         it_Vertex != VVerticesInThisContainer.end();
         ++it_Vertex)
    {
        VVertex2D* currVVertex = *it_Vertex;

        if (currVVertex->getStatus() == RED_V)
        {
            continue;
        }

        if (movingDisk->getDisk().getRadius() < currVVertex->computeRadiusOfTangentCircle())
        {
            thereIsPositionToInsert = true;
            insertPosition = currVVertex->getCircumcircle().getCenterPt();
            break;
        }
    }

    return make_pair(thereIsPositionToInsert, insertPosition);
}


void DynamicVoronoiDiagramCIC::update_circumcircles_of_VVertices_in_this_container(Generator2D * container)
{
    list<VVertex2D*> VVerticesInThisContainer;
    container->getInnerVD()->getVoronoiVertices(VVerticesInThisContainer);
    container->getInnerVD()->updateCircumcirclesOfNewVVertices(VVerticesInThisContainer);
}


void DynamicVoronoiDiagramCIC::find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
    const double & mostRecentlyEventOccurringTime, 
    const set<VEdge2D*>& influencedVEdges, 
    const DVD_VEdgeStatus& DVD_VEdgeStatus)
{
    for (set<VEdge2D*>::iterator it_VEdge = influencedVEdges.begin();
         it_VEdge != influencedVEdges.end();
         it_VEdge++)
    {
        const VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleCollisionTime(false, DBL_MAX);

        if (DVD_VEdgeStatus != DELETED_EDGE)
        {
            DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(VEdge->getLeftFace()->getGenerator()->getUserData());
            DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(VEdge->getRightFace()->getGenerator()->getUserData());

            possibleCollisionTime = find_possible_collision_time_of_disks_after(mostRecentlyEventOccurringTime, movingDisk1, movingDisk2);
            do_bookkeeping_on_priorityQ_for_disk_collision(possibleCollisionTime, VEdge, DVD_VEdgeStatus);

            if (m_NumOfContainers > 1)
            {
                DVD_StatusOfTwoDisks statusOfTwoDisks = find_state_of_two_disks(movingDisk1, movingDisk2);
                switch (statusOfTwoDisks)
                {
                    case DISK_N_DISK:
                    case CONTAINER_N_CONTAINER:
                    {
                        if (!this_is_collision_event_NOT_disk_hopping_event(VEdge))
                        {
                            m_InThisContainerDiskWillBeInserted.erase(const_cast<VEdge2D*>(VEdge));
                        }

                        break;
                    }

                    case DISK_N_CONTAINER:
                    {
                        do_bookkeeping_for_disk_hopping(
                            mostRecentlyEventOccurringTime, VEdge, movingDisk1, movingDisk2, DVD_VEdgeStatus, possibleCollisionTime.first);
                        break;
                    }

                    case CONTAINER_N_DISK:
                    {
                        do_bookkeeping_for_disk_hopping(
                            mostRecentlyEventOccurringTime, VEdge, movingDisk2, movingDisk1, DVD_VEdgeStatus, possibleCollisionTime.first);
                        break;
                    }
                }
            }
        }
        else
        {
            possibleCollisionTime = make_pair(false, DBL_MAX);
            do_bookkeeping_on_priorityQ_for_disk_collision(possibleCollisionTime, VEdge, DELETED_EDGE);

            if (m_NumOfContainers > 1)
            {
                do_bookkeeping_for_disk_hopping(
                    mostRecentlyEventOccurringTime, VEdge, NULL, NULL, DELETED_EDGE, possibleCollisionTime.first);
            }
        }

#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::DISK_COLLISION, const_cast<VEdge2D*>(VEdge), possibleCollisionTime.second);
#endif //DEBUG_DVD
    }


}


void DynamicVoronoiDiagramCIC::find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
    const double & mostRecentlyEventOccurringTime, const set<VEdge2D*>& influencedEdges, const DVD_VEdgeStatus & DVD_VEdgeStatus)
{
    for (set<VEdge2D*>::const_iterator it_VEdge = influencedEdges.begin();
        it_VEdge != influencedEdges.end();
        it_VEdge++)
    {
        const VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleEdgeFlipingTime(false, DBL_MAX);


        if (DVD_VEdgeStatus != DELETED_EDGE)
        {
            if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
            {
                possibleEdgeFlipingTime.first = false;
            }
            else
            {
                possibleEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(mostRecentlyEventOccurringTime, VEdge);

            }
        }
        else
        {
            possibleEdgeFlipingTime.first = false;
        }

#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::EDGE_FLIP, const_cast<VEdge2D*>(VEdge), possibleEdgeFlipingTime.second);
#endif //DEBUG_DVD

        do_bookkeeping_on_priorityQ_for_edge_flipping(possibleEdgeFlipingTime, VEdge, DVD_VEdgeStatus);
    }
}

void DynamicVoronoiDiagramCIC::find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
    const VEdge2D * targetVEdge, const double & mostRecentlyEventOccurringTime, 
    const set<VEdge2D*>& influencedEdges, const DVD_VEdgeStatus & DVD_VEdgeStatus)
{
    pair<bool, double> possibleTargetEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(mostRecentlyEventOccurringTime, targetVEdge, true);


    for (set<VEdge2D*>::const_iterator it_VEdge = influencedEdges.begin();
        it_VEdge != influencedEdges.end();
        it_VEdge++)
    {
        const VEdge2D* VEdge = *it_VEdge;

        pair<bool, double> possibleEdgeFlipingTime(false, DBL_MAX);

        if (VEdge->isInfinite() || VEdge->isUnBounded()) //ie, if this vedge is a virtual vedge : ������
        {
            possibleEdgeFlipingTime.first = false;
        }
        else if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
        {
            possibleEdgeFlipingTime.first = false;
        }
        //else if (this_VEdge_is_in_this_set(VEdge, m_VEdgesDefinedByIdentical4DiskWithFlippingVEdge))
        else if (this_VEdge_is_in_this_set(VEdge, m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ))
        {
            possibleEdgeFlipingTime = possibleTargetEdgeFlipingTime;
        }
        else
        {
            possibleEdgeFlipingTime = find_possible_flipping_time_of_an_vedge_after(mostRecentlyEventOccurringTime, VEdge);
        }


#ifdef DEBUG_DVD
        ___write_update_event_info(EventOfDynamicVD2D::EDGE_FLIP, const_cast<VEdge2D*>(VEdge), possibleEdgeFlipingTime.second);
#endif //DEBUG_DVD


        do_bookkeeping_on_priorityQ_for_edge_flipping(possibleEdgeFlipingTime, VEdge, DVD_VEdgeStatus);
    }
}



VEdge2D* DynamicVoronoiDiagramCIC::get_VEdge_with_this_generator_quadruplet(
    const EdgeFlipEvent::GeneratorPtrQuadruplet& generatorQuadruplet, const bool& VEdgeDefinedByQuadrupletIsAlreadyFlipped /* = false  */) const
{
    VEdge2D* outputVEdge = NULL;

    //index 0 : left, 1 : right, 2 : head, 3: tail
    list<VEdge2D*> boundingVEdgesOfTargetGenerator;

    Generator2D* leftGeneratorOfQaudruplet  = generatorQuadruplet[0];
    Generator2D* rightGeneratorOfQaudruplet = generatorQuadruplet[1];
    Generator2D* headGeneratorOfQaudruplet  = generatorQuadruplet[2];
    Generator2D* tailGeneratorOfQaudruplet  = generatorQuadruplet[3];

    if (!VEdgeDefinedByQuadrupletIsAlreadyFlipped)
    {
        DynamicDisk* leftDisk  = static_cast<DynamicDisk*>(leftGeneratorOfQaudruplet->getUserData());
        DynamicDisk* rightDisk = static_cast<DynamicDisk*>(rightGeneratorOfQaudruplet->getUserData());

        if (leftDisk->this_disk_is_container())
        {
            if (rightDisk->getCircle().isIncludedIn(leftDisk->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
            {
                leftGeneratorOfQaudruplet->getInnerFace()->getBoundaryVEdges(boundingVEdgesOfTargetGenerator);
            }
            else
            {
                leftGeneratorOfQaudruplet->getOuterFace()->getBoundaryVEdges(boundingVEdgesOfTargetGenerator);
            }
        }
        else
        {
            leftGeneratorOfQaudruplet->getOuterFace()->getBoundaryVEdges(boundingVEdgesOfTargetGenerator);
        }
    }
    else
    {
        DynamicDisk* headDisk = static_cast<DynamicDisk*>(headGeneratorOfQaudruplet->getUserData());

        DynamicDisk* rightDisk = static_cast<DynamicDisk*>(rightGeneratorOfQaudruplet->getUserData());

        if (headDisk->this_disk_is_container())
        {
            if (rightDisk->getCircle().isIncludedIn(headDisk->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
            {
                headGeneratorOfQaudruplet->getInnerFace()->getBoundaryVEdges(boundingVEdgesOfTargetGenerator);
            }
            else
            {
                headGeneratorOfQaudruplet->getOuterFace()->getBoundaryVEdges(boundingVEdgesOfTargetGenerator);
            }
        }
        else
        {
            headGeneratorOfQaudruplet->getOuterFace()->getBoundaryVEdges(boundingVEdgesOfTargetGenerator);
        }
    }


    for (list<VEdge2D*>::const_iterator it_VEdge = boundingVEdgesOfTargetGenerator.begin();
        it_VEdge != boundingVEdgesOfTargetGenerator.end();
        ++it_VEdge)
    {
        VEdge2D* currVEdge = *it_VEdge;

        Generator2D* leftGeneratorOfCurrVEdge  = currVEdge->getLeftFace()->getGenerator();
        Generator2D* rightGeneratorOfCurrVEdge = currVEdge->getRightFace()->getGenerator();
        Generator2D* headGeneratorOfCurrVEdge  = currVEdge->getMateFace(currVEdge->getStartVertex())->getGenerator();
        Generator2D* tailtGeneratorOfCurrVEdge = currVEdge->getMateFace(currVEdge->getEndVertex())->getGenerator();

        if (!VEdgeDefinedByQuadrupletIsAlreadyFlipped)
        {
            if (leftGeneratorOfCurrVEdge == leftGeneratorOfQaudruplet
                && rightGeneratorOfCurrVEdge == rightGeneratorOfQaudruplet
                && headGeneratorOfCurrVEdge == headGeneratorOfQaudruplet
                && tailtGeneratorOfCurrVEdge == tailGeneratorOfQaudruplet)
            {
                outputVEdge = currVEdge;
                break;
            }
            else if (leftGeneratorOfCurrVEdge == rightGeneratorOfQaudruplet
                && rightGeneratorOfCurrVEdge == leftGeneratorOfQaudruplet
                && headGeneratorOfCurrVEdge == tailGeneratorOfQaudruplet
                && tailtGeneratorOfCurrVEdge == headGeneratorOfQaudruplet)
            {
                outputVEdge = currVEdge;
                break;
            }

        }
        else
        {
             if (leftGeneratorOfCurrVEdge    == headGeneratorOfQaudruplet
                && rightGeneratorOfCurrVEdge == tailGeneratorOfQaudruplet
                && headGeneratorOfCurrVEdge  == rightGeneratorOfQaudruplet
                && tailtGeneratorOfCurrVEdge == leftGeneratorOfQaudruplet)
            {
                outputVEdge = currVEdge;
                break;
            }
            else if (leftGeneratorOfCurrVEdge == tailGeneratorOfQaudruplet
                && rightGeneratorOfCurrVEdge  == headGeneratorOfQaudruplet
                && headGeneratorOfCurrVEdge   == leftGeneratorOfQaudruplet
                && tailtGeneratorOfCurrVEdge  == rightGeneratorOfQaudruplet)
            {
                outputVEdge = currVEdge;
                break;
            }

        } 
    }

    return outputVEdge;
}


EdgeFlipEvent::GeneratorPtrQuadruplet DynamicVoronoiDiagramCIC::get_generator_quadruplet_by_its_id(int * generatorIDs) const
{
    EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = { get_generator_by_its_input_disk_id(generatorIDs[0]),
                                                                     get_generator_by_its_input_disk_id(generatorIDs[1]),
                                                                     get_generator_by_its_input_disk_id(generatorIDs[2]),
                                                                     get_generator_by_its_input_disk_id(generatorIDs[3]) };

    return generatorPtrQuadruplet;
}


EdgeFlipEvent::GeneratorPtrQuadruplet DynamicVoronoiDiagramCIC::get_generator_quadruplet_by_its_VEdge(VEdge2D * VEdge) const
{
    EdgeFlipEvent::GeneratorPtrQuadruplet generatorPtrQuadruplet = { VEdge->getLeftFace()->getGenerator(),
                                                                     VEdge->getRightFace()->getGenerator(),
                                                                     VEdge->getMateFace(VEdge->getStartVertex())->getGenerator(),
                                                                     VEdge->getMateFace(VEdge->getEndVertex())->getGenerator() };

    return generatorPtrQuadruplet;
}




void DynamicVoronoiDiagramCIC::move_disks_and_update_VVertices_related_to_disk_hopping_event(DiskHoppingEvent * diskHoppingEvent, Generator2D* containerIncludingVVerticesToBeUpdated)
{
    move_disks_related_to_disk_hopping_event(diskHoppingEvent);

    update_circumcircles_of_VVertices_in_this_container(containerIncludingVVerticesToBeUpdated);
}


void DynamicVoronoiDiagramCIC::move_disks_related_to_disk_hopping_event(DiskHoppingEvent* diskHoppingEvent)
{
    Generator2D* fromThisContainer = diskHoppingEvent->get_in_what_container_this_generator_is_deleted();
    Generator2D* toThisContainer   = diskHoppingEvent->get_in_what_container_this_generator_is_inserted();

    list<Generator2D*> generatorsInFromContainer;
    fromThisContainer->getInnerGens(generatorsInFromContainer);

    list<Generator2D*> generatorsInToContainer;
    toThisContainer->getInnerGens(generatorsInToContainer);

    set<Generator2D*> influencedDiskGenerators;
    influencedDiskGenerators.insert(generatorsInFromContainer.begin(), generatorsInFromContainer.end());
    influencedDiskGenerators.insert(generatorsInToContainer.begin(), generatorsInToContainer.end());
    move_disks_to_target_Voronoi_clock_without_changing_VD(influencedDiskGenerators, diskHoppingEvent->get_occurring_clock());
}


void DynamicVoronoiDiagramCIC::move_disks_and_update_VVertices_related_to_disk_hopping_event_to_FUTURE(DiskHoppingEvent* diskHoppingEvent)
{
    Generator2D* toThisContainer = diskHoppingEvent->get_in_what_container_this_generator_is_inserted();
    move_disks_and_update_VVertices_related_to_disk_hopping_event(diskHoppingEvent, toThisContainer);
}


void DynamicVoronoiDiagramCIC::move_disks_and_update_VVertices_related_to_disk_hopping_event_to_PAST(DiskHoppingEvent* diskHoppingEvent)
{
    Generator2D* fromThisContainer = diskHoppingEvent->get_in_what_container_this_generator_is_deleted();
    move_disks_and_update_VVertices_related_to_disk_hopping_event(diskHoppingEvent, fromThisContainer);
}


void DynamicVoronoiDiagramCIC::update_geometry_of_VVertex(VVertex2D* VVertex)
{
    if (m_NumOfContainers == 1)
	{
        updateCircumcircle(VVertex);
        return;
    }

    list<Generator2D*> definingForVVertex;
    VVertex->getDefining3Generators(definingForVVertex);

    vector<DynamicDisk*> movingDisksForVVertex;
    for (list<Generator2D*>::const_iterator it_Generator = definingForVVertex.begin();
        it_Generator != definingForVVertex.end();
        ++it_Generator)
    {
        DynamicDisk* currDisk = static_cast<DynamicDisk*>((*it_Generator)->getUserData());
        movingDisksForVVertex.push_back(currDisk);
    }

    sort_generators_in_non_increasing_order(movingDisksForVVertex);

    DynamicDisk* largestMovingDisk       = movingDisksForVVertex[0];
    DynamicDisk* secondLargestMovingDisk = movingDisksForVVertex[1];

    if (largestMovingDisk->this_disk_is_container() && secondLargestMovingDisk->getCircle().isIncludedIn(largestMovingDisk->getCircle(), -TOLERANCE_OF_DISK_INTERSECTION))
    {
        Generator2D* container = get_generator_by_its_input_disk_id(largestMovingDisk->getID());

        list<VVertex2D*> oneVertexInList;
        oneVertexInList.push_back(VVertex);
        container->getInnerVD()->updateCircumcirclesOfNewVVertices(oneVertexInList);
    }
    else
    {
        updateCircumcircle(VVertex);
    }
}



void DynamicVoronoiDiagramCIC::clear_delete_RED_VEdges_and_VVertices()
{
    list<Generator2D*> containerGenerators;
    collectContainerGenerators(containerGenerators);

    for (list<Generator2D*>::const_iterator it_Generator = containerGenerators.begin();
        it_Generator != containerGenerators.end();
        ++it_Generator)
    {
        Generator2D* currContainer = *it_Generator;
        currContainer->getInnerVD()->removeREDVVerticesAndVEdges();
    }
}





bool DynamicVoronoiDiagramCIC::validate_DVD_by_empty_circle_test_at_vertices()
{
    compute_all_VVertex_coordinates_of_VD();

    list<VVertex2D*> VVertices;

    if (m_NumOfContainers > 1)
    {
        getAllVertices_hierarchy(VVertices);
    }
    else
    {
        getVoronoiVertices(VVertices);
    }

    for (list<VVertex2D*>::const_iterator iter_VVertex = VVertices.begin(); iter_VVertex != VVertices.end(); iter_VVertex++)
    {
        VVertex2D*  currVtx = *iter_VVertex;
        
        list<Generator2D*> defining3Generators;
        currVtx->getDefining3Generators(defining3Generators);

        DynamicDisk* defining3Disks[3];
        int count = 0;
        for (list<Generator2D*>::iterator it_Generator = defining3Generators.begin();
            it_Generator != defining3Generators.end();
            ++it_Generator)
        {
            DynamicDisk* dynamicDisk = static_cast<DynamicDisk*>((*it_Generator)->getUserData());
            defining3Disks[count++] = dynamicDisk;
        }

        rg_Circle2D currCircumCircle = currVtx->getCircumcircle();

        list<DynamicDisk*> allDynamicDisks;
        get_dynamic_disks(allDynamicDisks);

        for (list<DynamicDisk*>::const_iterator it_DynamicDisk = allDynamicDisks.begin(); it_DynamicDisk != allDynamicDisks.end(); it_DynamicDisk++)
        {
            const DynamicDisk* currMovingDisk = (*it_DynamicDisk);

            if (currMovingDisk == defining3Disks[0]
             || currMovingDisk == defining3Disks[1]
             || currMovingDisk == defining3Disks[2])
            {
                continue;
            }

            const rg_Circle2D& currDisk = currMovingDisk->getCircle();
            /*
            if (currDisk.isIntersectWith(currCircumCircle) || currDisk.isIncludedIn(currCircumCircle))
            {
            return false;
            }
            */

            if (there_is_intersection_btw(currCircumCircle, currDisk, TOLERANCE_OF_DISK_INTERSECTION) || currDisk.isIncludedIn(currCircumCircle))
            {
                return false;
            }
        }
    }

    return true;
}


void DynamicVoronoiDiagramCIC::get_VVertices_except_infinite_ones(list<VVertex2D*>& VVerticesExceptInfiniteOnes)
{
    list<VVertex2D*> allVVertices;
    getVoronoiVertices(allVVertices);

    for (list<VVertex2D*>::iterator it_VVertex = allVVertices.begin(); it_VVertex != allVVertices.end(); it_VVertex++)
    {
        VVertex2D* currVertex = *it_VVertex;

        if (!currVertex->isInfinite())
        {
            VVerticesExceptInfiniteOnes.push_back(currVertex);
        }
    }
}



bool DynamicVoronoiDiagramCIC::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2)
{
    bool isLargerThan = false;

    rg_Circle2D disk1(movingDisk1->getCircle());
    rg_Circle2D disk2(movingDisk2->getCircle());

    if (disk1.getRadius() > disk2.getRadius())
    {
        isLargerThan = true;
    }
    else if (disk1.getRadius() == disk2.getRadius())
    {

        if (disk1.getCenterPt().getX() > disk2.getCenterPt().getX())
        {
            isLargerThan = true;
        }
        else if (disk1.getCenterPt().getX() == disk2.getCenterPt().getX())
        {


            if (disk1.getCenterPt().getY() > disk2.getCenterPt().getY())
            {
                isLargerThan = true;
            }
        }
    }

    return isLargerThan;
}

bool DynamicVoronoiDiagramCIC::compare_disk_with_ID_in_non_decreasing_order(const DynamicDisk * movingDisk1, const DynamicDisk * movingDisk2)
{
    if (movingDisk1->getID() < movingDisk2->getID())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void DynamicVoronoiDiagramCIC::compute_all_VEdges_equation_of_VD()
{
    if (m_NumOfContainers > 1)
    {
        list<Generator2D*> containerGenerators;
        collectContainerGenerators(containerGenerators);

        for (list<Generator2D*>::const_iterator it_Generator = containerGenerators.begin();
            it_Generator != containerGenerators.end();
            ++it_Generator)
        {
            Generator2D* currContainer = *it_Generator;
            currContainer->getInnerVD()->updateGeometry();
        }
    }
    else
    {
        updateGeometry();
    }
}


#ifdef DEBUG_DVD
void DynamicVoronoiDiagramCIC::___initialize_debug() const
{
    string fileName(debugFileNameWithPath);
    ofstream fout(fileName.c_str(), ios_base::out);

    fout << "[Debug in event generation]" << endl;
    fout << " * File Name : " << "\t" << fileName << endl;

    fout.close();
}


void DynamicVoronoiDiagramCIC::___write_current_event_info(EventOfDynamicVD2D * anEvent) const
{
    string fileName(debugFileNameWithPath);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);
    fout.precision(6);

    fout << "****************************************" << endl;
    fout << "***          CURR EVENT              ***" << endl;
    fout << "****************************************" << endl;
    fout << endl;
   
    char type;


    switch (anEvent->get_event_type())
    {
        case EventOfDynamicVD2D::EDGE_FLIP:
        {
            EdgeFlipEvent* edgeFlippingEvent = static_cast<EdgeFlipEvent*>(anEvent);

            DynamicDisk* leftDisk  = static_cast<DynamicDisk*>(edgeFlippingEvent->get_generator_quadruplet_before_flip().at(0)->getUserData());

            DynamicDisk* rightDisk = static_cast<DynamicDisk*>(edgeFlippingEvent->get_generator_quadruplet_before_flip().at(1)->getUserData());

            DynamicDisk* headDisk  = static_cast<DynamicDisk*>(edgeFlippingEvent->get_generator_quadruplet_before_flip().at(2)->getUserData());

            DynamicDisk* tailDisk  = static_cast<DynamicDisk*>(edgeFlippingEvent->get_generator_quadruplet_before_flip().at(3)->getUserData());


            int IDOfLeftGenerator  = leftDisk->getID();
            int IDOfRightGenerator = rightDisk->getID();
            int IDOfHeadGenerator  = headDisk->getID();
            int IDOfTailGenerator  = tailDisk->getID();


            if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
            {
                type = 'E';
            }
            else
            {
                type = 'S';
            }

            fout << "*. Type : " << type                                    << "\t"; 
            fout << "*. Time : " << edgeFlippingEvent->get_occurring_clock() << "\t";
            //fout << "*. Edge : " << edgeFlippingEvent->get_VEdge()->getID() << "\t";
            fout << "*. 4Gen(LRHT) : " << IDOfLeftGenerator << "\t" << IDOfRightGenerator << "\t" << IDOfHeadGenerator << "\t" << IDOfTailGenerator << endl;

            break;
        }

        case EventOfDynamicVD2D::DISK_COLLISION:
        {
            DiskCollisionEvent* diskCollisionEvent = static_cast<DiskCollisionEvent*>(anEvent);

            int IDOfGenerator1 = diskCollisionEvent->get_disk_generator1()->getID();
            int IDOfGenerator2 = diskCollisionEvent->get_disk_generator2()->getID();

            rg_Point2D newVelocityOfGen1 = diskCollisionEvent->get_future_velocity_vector_of_disk1();
            rg_Point2D newVelocityOfGen2 = diskCollisionEvent->get_future_velocity_vector_of_disk2();

            type = 'C';

            fout << "*. Type : " << type                                     << "\t"; 
            fout << "*. Time : " << diskCollisionEvent->get_occurring_clock() << "\t";
            fout << "*. Gen(1) : " << IDOfGenerator1 << "\t" << "(" << newVelocityOfGen1.getX() << "," << newVelocityOfGen1.getY() << ")" << "\t ,";
            fout << "*. Gen(2) : " << IDOfGenerator2 << "\t" << "(" << newVelocityOfGen2.getX() << "," << newVelocityOfGen2.getY() << ")" << endl;

            break;
        }

        case EventOfDynamicVD2D::DISK_VELOCITY_CHANGE:
        {
            DiskVelocityChangeEvent* diskVelocityChangeEvent = static_cast<DiskVelocityChangeEvent*>(anEvent);
           
            int IDOfGenerator      = diskVelocityChangeEvent->get_disk_generator()->getID();
            rg_Point2D newVelocity = diskVelocityChangeEvent->get_future_velocity_vector_of_disk();
            
            type = 'V';
            fout << "*. Type : " << type << "\t";
            fout << "*. Time : " << diskVelocityChangeEvent->get_occurring_clock() << "\t";
            fout << "*. Gen  : " << IDOfGenerator      << "\t" << "("
                                 << newVelocity.getX() << "," 
                                 << newVelocity.getY() << ")" << endl;

            break;
        }
    }
}


void DynamicVoronoiDiagramCIC::___write_update_event_info(const EventOfDynamicVD2D::EventType & eventType, VEdge2D * VEdge, const double newEventTime) const
{
    string fileName(debugFileNameWithPath);
    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);
    fout.precision(6);

    switch (eventType)
    {
        case EventOfDynamicVD2D::EDGE_FLIP:
        {
            int IDOfLeftGenerator = VEdge->getLeftFace()->getGenerator()->getID();
            int IDOfRightGenerator = VEdge->getRightFace()->getGenerator()->getID();
            int IDOfHeadGenerator = VEdge->getMateFace(VEdge->getStartVertex())->getGenerator()->getID();
            int IDOfTailGenerator = VEdge->getMateFace(VEdge->getEndVertex())->getGenerator()->getID();

            double pastEventTime = m_PriorityQForEdgeFlipping.findKey(VEdge);
            string edgeType;

            if (VEdge->isInfinite() || VEdge->isUnBounded()) //ie, if this vedge is a virtual vedge : ������
            {
                edgeType = "Infinite";
            }
            else if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
            {
                edgeType = "Anomaly ";
            }
            else
            {
                int IDOfLeftGenerator = VEdge->getLeftFace()->getGenerator()->getID();
                int IDOfRightGenerator = VEdge->getRightFace()->getGenerator()->getID();
                int IDOfHeadGenerator = VEdge->getMateFace(VEdge->getStartVertex())->getGenerator()->getID();
                int IDOfTailGenerator = VEdge->getMateFace(VEdge->getEndVertex())->getGenerator()->getID();

                double pastEventTime = m_PriorityQForEdgeFlipping.findKey(VEdge);
                string edgeType;

                if (VEdge->isInfinite() || VEdge->isUnBounded()) //ie, if this vedge is a virtual vedge : ������
                {
                    edgeType = "Infinite";
                }
                else if (this_VEdge_is_boundary_of_anomalizing_cell(VEdge)) //(ie, if this edge is a boundary edge of an anomaly cell)
                {
                    edgeType = "Anomaly ";
                }
                else
                {
                    edgeType = "Normal  ";
                }

                fout << "\t[" << VEdge->getID() << "]\t" << "E" << "\t" << edgeType << "\t" << pastEventTime << "\t" << newEventTime << "\t";
                fout << "G : { " << IDOfLeftGenerator << "\t" << IDOfRightGenerator
                    << "\t" << IDOfHeadGenerator << "\t" << IDOfTailGenerator << " }" << endl;
            }

            break;
        }

        case EventOfDynamicVD2D::DISK_COLLISION:
        {
            int IDOfGenerator1 = VEdge->getLeftFace()->getGenerator()->getID();
            int IDOfGenerator2 = VEdge->getRightFace()->getGenerator()->getID();

            double pastEventTime = m_PriorityQForDiskCollision.findKey(VEdge);
            string edgeType;

            if (VEdge->isInfinite() || VEdge->isUnBounded())
            {
                edgeType = "Inf. VEdge";
            }
            else if(VEdge->getLeftFace()->isUnBounded() || VEdge->getRightFace()->isUnBounded()) // dynamicDisk ���忡���� phantom�� �𸣴� ��
            {
                edgeType = "Inf. VFace";
            }
            else
            {
                edgeType = "Normal    ";
            }

            fout << "\t[" << VEdge->getID() << "]\t" << "C" << "\t" << edgeType << "\t" << pastEventTime << "\t" << newEventTime << "\t";
            fout << "G : { " << IDOfGenerator1 << "\t" << IDOfGenerator2 << " }" << endl;
            break;
        }
    }
}
#endif //DEBUG_DVD


#ifdef WRITING_EVENT_TIME

void DynamicVoronoiDiagramCIC::___write_event_times(const double& eventTime, const double& computationTime)
{
    string fileName(m_DiskFileNameWthPath + "-" + eventTimeFileNameWithPath);

    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);
    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(10);

    accumulatedTimeInEventGeneration += computationTime;

    fout << eventTime  << "\t" << accumulatedTimeInEventGeneration << endl;
    fout.close();
}


void DynamicVoronoiDiagramCIC::___initialize_writing_event_time()
{
    string fileName(m_DiskFileNameWthPath + "-" + eventTimeFileNameWithPath);

    ofstream fout(fileName.c_str(), ios_base::out | ios_base::app);
    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(10);

    accumulatedTimeInEventGeneration =  m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_CONSTRUCT_INITIAL_CIC_VD)
                                   + m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_DUPLICATE_INITIAL_VD_N_DISKS_AND_DEFINE_CORRESPONDANCE_MAP_WITH_ORIGINAL_ONES)
                                   + m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS)
                                   + m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME)
                                   + m_TimeStatistics[PHASE1_COMP_TIME].time(DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME);

    fout << "Event Time" << "\t" << "Comp. Time" << endl;
    fout.close();
}

#endif //WRITING_EVENT_TIME



void DynamicVoronoiDiagramCIC::compute_all_VVertex_coordinates_of_VD()
{
    if (m_NumOfContainers > 1)
    {
        list<Generator2D*> containers;
        collectContainerGenerators(containers);

        for (list<Generator2D*>::iterator i_gen = containers.begin(); i_gen != containers.end(); ++i_gen) {
            Generator2D* currContainer = *i_gen;
            VoronoiDiagramCIC* innerVD = const_cast<VoronoiDiagramCIC*>(currContainer->getInnerVD());

            list<VVertex2D*> vertices;
            innerVD->getVoronoiVertices(vertices);
            innerVD->updateCircumcirclesOfNewVVertices(vertices);
        }
    }
    else
    {
        list<VVertex2D*> allVVerticesExceptInfiniteVVertices;

        get_VVertices_except_infinite_ones(allVVerticesExceptInfiniteVVertices);

        compute_circum_circles_Of_VVertices(allVVerticesExceptInfiniteVVertices);
    }
}



void DynamicVoronoiDiagramCIC::compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(DiskCollisionEvent * diskCollisionEvent) const
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    rg_Point2D newVelocityVector1, newVelocityVector2;

    compute_velocity_vectors_of_disks_after_collided(movingDisk1, movingDisk2, newVelocityVector1, newVelocityVector2);

    store_two_disks_velocity_vectors_before_and_after_collision(diskCollisionEvent, newVelocityVector1, newVelocityVector2);

    update_velocity_vector(movingDisk1, newVelocityVector1);

    update_velocity_vector(movingDisk2, newVelocityVector2);
}



void DynamicVoronoiDiagramCIC::store_current_positions_of_two_colliing_disks_to_event(DiskCollisionEvent* diskCollisionEvent) const
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    diskCollisionEvent->set_disk_centers_at_collision_time(movingDisk1->getCircle().getCenterPt(), movingDisk2->getCircle().getCenterPt());
}



void DynamicVoronoiDiagramCIC::change_positions_of_two_colliding_disks_from_event(DiskCollisionEvent* diskCollisionEvent) const
{
    Generator2D* generator1 = diskCollisionEvent->get_disk_generator1();
    Generator2D* generator2 = diskCollisionEvent->get_disk_generator2();

    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(generator1->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(generator2->getUserData());

    rg_Circle2D newCircle1(diskCollisionEvent->get_center_of_disk1(), movingDisk1->getCircle().getRadius());
    rg_Circle2D newCircle2(diskCollisionEvent->get_center_of_disk2(), movingDisk2->getCircle().getRadius());

    generator1->setDisk(newCircle1);
    generator2->setDisk(newCircle2);

    movingDisk1->setCircle(newCircle1);
    movingDisk2->setCircle(newCircle2);

    const double targetVoronoiClock = diskCollisionEvent->get_occurring_clock();
    movingDisk1->setMostRecentLocationUpdateTime(targetVoronoiClock);
    movingDisk2->setMostRecentLocationUpdateTime(targetVoronoiClock);
}


void DynamicVoronoiDiagramCIC::store_current_position_of_disk_to_event(DiskVelocityChangeEvent* diskVelocityChangeEvent) const
{
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskVelocityChangeEvent->get_disk_generator()->getUserData());

    diskVelocityChangeEvent->set_disk_center_at_velocity_change_time(movingDisk->getCircle().getCenterPt());
}


void DynamicVoronoiDiagramCIC::change_position_of_disk_from_event(DiskVelocityChangeEvent* diskVelocityChangeEvent) const
{
    Generator2D* generator  = diskVelocityChangeEvent->get_disk_generator();
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(generator->getUserData());

    rg_Circle2D newCircle(diskVelocityChangeEvent->get_center_of_disk(), movingDisk->getCircle().getRadius());

    generator->setDisk(newCircle);
    movingDisk->setCircle(newCircle);

    const double targetVoronoiClock = diskVelocityChangeEvent->get_occurring_clock();
    movingDisk->setMostRecentLocationUpdateTime(targetVoronoiClock);
}



void DynamicVoronoiDiagramCIC::change_position_of_disk_from_event(DiskHoppingEvent* diskVelocityChangeEvent) const
{
    Generator2D* generator = diskVelocityChangeEvent->get_target_generator();
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(generator->getUserData());

    rg_Circle2D newCircle(diskVelocityChangeEvent->get_deleted_position(), movingDisk->getCircle().getRadius());

    generator->setDisk(newCircle);
    movingDisk->setCircle(newCircle);

    const double targetVoronoiClock = diskVelocityChangeEvent->get_occurring_clock();
    movingDisk->setMostRecentLocationUpdateTime(targetVoronoiClock);
}


void DynamicVoronoiDiagramCIC::store_this_FLIP_event_in_the_event_history_vector(EdgeFlipEvent* edgeFlippingEvent)
{
    m_EventSequence.push_back(edgeFlippingEvent);

    clear_info_about_vedges_defined_by_same_disks();

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_EDGE_FLIPPING_IN_EVENT_HISTORY);
}


void DynamicVoronoiDiagramCIC::store_this_DISK_COLLISION_event_in_the_event_history_vector(DiskCollisionEvent* diskCollisionEvent)
{
    m_EventSequence.push_back(diskCollisionEvent);

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DISK_COLLISION_IN_EVENT_HISTORY);

    DynamicDisk* disk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* disk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    if (!disk1->this_disk_is_container() && !disk2->this_disk_is_container())
    {
         ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DISK_COLLISION_BTW_INPUT_DISKS_NOT_WITH_CONTAINER_IN_EVENT_HISTORY);
    }
}


void DynamicVoronoiDiagramCIC::store_this_VELOCITY_VECTOR_CHANGE_event_in_the_event_history_vector(DiskVelocityChangeEvent* diskVelocityChangeEvent)
{
    m_EventSequence.push_back(diskVelocityChangeEvent);

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DISK_VELOCITY_CHANGE_IN_EVENT_HISTORY);
}


void DynamicVoronoiDiagramCIC::store_this_DISK_HOPPING_event_in_the_event_history_vector(DiskHoppingEvent* diskHoppingEvent)
{
    m_EventSequence.push_back(diskHoppingEvent);

    clear_info_about_vedges_related_to_disk_hopping();

    ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DISK_HOPPING_IN_EVENT_HISTORY);
}


void DynamicVoronoiDiagramCIC::sort_generators_in_non_increasing_order(vector<DynamicDisk*>& disks)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_SORTING_GENERATORS_IN_NON_INCRESING_ORDER);

    sort(disks.begin(), disks.end(), DynamicVoronoiDiagramCIC::compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order);

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_SORTING_GENERATORS_IN_NON_INCRESING_ORDER);
}




Generator2D* DynamicVoronoiDiagramCIC::choose_next_container_by_transition_probability(const int& currContainerID) const
{
    int currContainerIndexOfCHTM = convert_container_id_to_transition_probability_index(currContainerID);

    Generator2D* nextContainer = NULL;

    float valueBtwZeroToOne = (double)(rand()) / (double)(RAND_MAX);

    float minInterval = 0.0;

#ifdef NEW_CONTAINER_ID_MAPPING
    for (int i = 0; i < m_NumOfContainers -1; ++i)
#else
    for (int i = 0; i < m_NumOfContainers; ++i)
#endif
    {
        float maxInterval = minInterval + get_transition_probability(currContainerIndexOfCHTM, i);

        if (valueBtwZeroToOne >= minInterval && valueBtwZeroToOne <= maxInterval)
        {
            int containerId = convert_transition_probability_index_to_container_id(i);
            nextContainer   = get_generator_by_its_input_disk_id(containerId);
            break;
        }

        minInterval = maxInterval;
    }

    return nextContainer;
}



std::pair<bool, rg_Point2D> DynamicVoronoiDiagramCIC::check_two_conditions_for_hopping_and_get_inserted_position_if_satisfies_the_conditions(DiskHoppingEvent* diskHoppingEvent)
{
    Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
    Generator2D* fromContainer   = diskHoppingEvent->get_in_what_container_this_generator_is_deleted();
    Generator2D* toContainer     = diskHoppingEvent->get_in_what_container_this_generator_is_inserted();

    // condition 1.
    int numOfDisksInFromContainer = find_number_of_generator_disks_in_this_container(fromContainer);

    // condition 2.
    pair<bool, rg_Point2D> newPosition(false, rg_Point2D(DBL_MAX, DBL_MAX));
    if (toContainer->getInnerVD() != NULL) // if toContainer works as virtual container.
    {
        update_circumcircles_of_VVertices_in_this_container(toContainer);
        newPosition = find_insert_position_of_disk_in_this_container_using_VVertices(targetGenerator, toContainer);
    }


    if (numOfDisksInFromContainer >= 3 && newPosition.first)
    {
        return make_pair(true, newPosition.second);
    }
    else
    {
        return make_pair(false, rg_Point2D(DBL_MAX, DBL_MAX));
    }
}



void DynamicVoronoiDiagramCIC::store_velocity_vectors_before_and_after_disk_hopping_by_current_velocity_vector(
    DiskHoppingEvent* diskHoppingEvent) const
{
    Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
    DynamicDisk* targetDisk      = static_cast<DynamicDisk*>(targetGenerator->getUserData());

    diskHoppingEvent->set_velocity_vector_after_collision(targetDisk->getVelocityVector());
    diskHoppingEvent->set_velocity_vector_before_collision(targetDisk->getVelocityVector());
    diskHoppingEvent->set_velocity_vector_setting(true);
}


void DynamicVoronoiDiagramCIC::store_current_and_target_positions_for_disk_hopping(DiskHoppingEvent* diskHoppingEvent, const rg_Point2D& targetPosition) const
{
    Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
    rg_Point2D currPosition      = targetGenerator->getDisk().getCenterPt();

    diskHoppingEvent->set_deleted_n_inserted_position(currPosition, targetPosition);
}


void DynamicVoronoiDiagramCIC::store_current_position_of_disk_to_event(DiskHoppingEvent* diskHoppingEvent) const
{
    Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
    DynamicDisk* targetDisk      = static_cast<DynamicDisk*>(targetGenerator->getUserData());
    rg_Point2D currPosition      = targetGenerator->getDisk().getCenterPt();
    diskHoppingEvent->set_velocity_vector_change_position(currPosition);
}


void DynamicVoronoiDiagramCIC::store_velocity_vector_before_and_after_collision(
    DiskHoppingEvent* diskHoppingEvent, const rg_Point2D& velocityVectorAfterCollding) const
{
    Generator2D* targetGenerator = diskHoppingEvent->get_target_generator();
    DynamicDisk* targetDisk      = static_cast<DynamicDisk*>(targetGenerator->getUserData());

    diskHoppingEvent->set_velocity_vector_before_collision(targetDisk->getVelocityVector());
    diskHoppingEvent->set_velocity_vector_after_collision(velocityVectorAfterCollding);
    diskHoppingEvent->set_velocity_vector_setting(true);
}










/********************************************************/
/***************     TO BE REMOVED     ******************/
/********************************************************/

void DynamicVoronoiDiagramCIC::find_VEdges_defined_by_identical_4disks_from_PQ_N_flipping_VEdge(const VEdge2D * flippingVEdge, const set<VEdge2D*> upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ, set<VEdge2D*>& VEdgesDefinedByIdentical4Disk) const
{
    set<VEdge2D*> totalSet;

    list<VEdge2D*> adjacentVEdges;
    flippingVEdge->getAdjacentVEdges(adjacentVEdges);

    totalSet.insert(adjacentVEdges.begin(), adjacentVEdges.end());
    totalSet.insert(upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin(), upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.end());

    find_all_VEdges_defined_by_identical_disks_TOPOLOGICAL_COMPUTATION(flippingVEdge, totalSet, VEdgesDefinedByIdentical4Disk);
}
