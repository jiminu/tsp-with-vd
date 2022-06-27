#ifndef _DYNAMIC_VORONOIDIAGRAM_2D_
#define _DYNAMIC_VORONOIDIAGRAM_2D_

//#define USING_BUCKET_FOR_CHECKING_ALL_PAIR_WISE_DISK_INTERSECTION
#define WRITE_HEXA_VALUES_FOR_DOUBLE_TYPE
//#define WRITING_POLYNOMIAL_RELATED_INFORMATION
//#define WRITING_EVENT_TIME
//#define DEBUG_DVD
//#define USING_BIG_VD
//#define PQ_STORING_ALL_VEDGES
#define NEW_CONTAINER_ID_MAPPING

#include "VoronoiDiagramCIC.h"

#include "SimpleTypesForDVD2.h"
#include "DefaultValuesForDVD2.h"

#include "EventOfDynamicVD2D.h"
#include "EdgeFlipEvent.h"
#include "DiskCollisionEvent.h"
#include "DiskVelocityChangeEvent.h"
#include "DiskHoppingEvent.h"
#include "EntityAccessiblePriorityQ.h"
#include "DynamicDisk.h"
#include "DiskCollisionHandler.h"
#include "DiskCollisionWatcher.h"
#include "TimeStatistics.h"
#include "CountStatistics.h"
#include "rg_Polynomial.h"
#include "rg_RelativeOp.h"
#include "rg_Triplet.h"
#include <random>
#include <chrono>

#include <list>
#include <array>
#include <vector>
#include <map>
using namespace std;


namespace V {
namespace GeometryTier {


class DynamicVoronoiDiagramCIC : public VoronoiDiagramCIC
{
public: 
    friend class DynamicVD2DSimulator;

private:
    //For input parameter               
    double                              m_LengthOfVoronoiClock;
    DiskCollisionHandler                m_CollisionHandlingType;
    DiskCollisionWatcher                m_CollisionWatcher;

    vector<float>                       m_DiskHoppingProbabilityMatrix;
    int                                 m_NumOfContainers;

    //For output
    vector<EventOfDynamicVD2D*>         m_EventSequence;

    //For current state
    int                                 m_IndexOfMostRecentEvent;
    double                              m_CurrentVoronoiClock;

    bool                                m_DuringEVENTSEQGeneration;
    bool                                m_Vers1p0EVENTSEQFileUsed;
    bool                                m_InitialPQIsConstructed;
    bool                                m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ;

    //For simulation
    unordered_map<int, Generator2D*>    m_LinkFromInputIDToGenerator; //including container


    //For reinstatement
	VoronoiDiagramCIC                         m_InitialVDForReinstatement;
    unordered_map<Generator2D*, Generator2D*> m_MapFromCopiedInitialStateGeneratorsToOriginalGenerators;
    unordered_map<DynamicDisk*, DynamicDisk>  m_MapFromOriginalDisksToCopiedInitialStateDisks;

    //******************************** For event generation ********************************//
        //After phase1, the elements of these data structure will be cleared.

        //For calculating events
        EntityAccessiblePriorityQ<VEdge2D*>                        m_PriorityQForEdgeFlipping;
        EntityAccessiblePriorityQ<VEdge2D*>                        m_PriorityQForDiskCollision;
        EntityAccessiblePriorityQ<ReservationNumber>               m_PriorityQForDiskVelocityChange;
        unordered_map<ReservationNumber, VelocityVectorChange>     m_ReservedVelocityVectorChange;
        ReservationNumber                                          m_LastReservedNumber;

        //EntityAccessiblePriorityQ_<pair<Generator2D*, rg_Point2D>> m_PriorityQForDiskVelocityChange;
        

        unordered_map<VEdge2D*, Generator2D*>                      m_InThisContainerDiskWillBeInserted;

        //For special case
        set<VEdge2D*>                       m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ;
        set<VEdge2D*>                       m_VEdgesDefinedByIdentical4DiskWithFlippingVEdge;

        set<VEdge2D*>                       m_NewVEdgesByDiskHopping;
        set<VEdge2D*>                       m_DeletedVEdgesByDiskHopping;
        set<VEdge2D*>                       m_RedefinedVEdgesByDiskHopping;

    //**************************************************************************************//

    //For statistics                    
    TimeStatistics<milli>               m_TimeStatistics[DVD_TIME_STATISTICS_SIZE];
    CountStatistics                     m_CountStatistics[DVD_COUNT_STATISTICS_SIZE];
    
public:
	//Constructor
    DynamicVoronoiDiagramCIC();
    DynamicVoronoiDiagramCIC(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks);
    DynamicVoronoiDiagramCIC(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks);
	DynamicVoronoiDiagramCIC(const DynamicVoronoiDiagramCIC& dynamicVD);

	//Deconstructor
	~DynamicVoronoiDiagramCIC();

	//Getter
    inline double get_current_Voronoi_clock()             const { return m_CurrentVoronoiClock; };
    inline double get_length_of_Voronoi_clock()           const { return m_LengthOfVoronoiClock; };

	inline void                                get_EVENTSEQ(vector<EventOfDynamicVD2D*>& eventSequence) const { eventSequence = m_EventSequence; };
	inline const vector<EventOfDynamicVD2D*>&  get_EVENTSEQ()                                            const { return m_EventSequence; };
	
	inline EventOfDynamicVD2D*                 get_event(const int& idx)                const { return (idx >= 0 && idx < m_EventSequence.size()) ? m_EventSequence.at(idx) : NULL; };
	inline double                              get_event_occuring_clock(const int& idx) const { return (idx >= 0) ? m_EventSequence.at(idx)->get_occurring_clock() : DBL_MAX; };
    inline const EventOfDynamicVD2D*           get_most_recent_event()                  const { return (m_IndexOfMostRecentEvent == -1) ? NULL : m_EventSequence.at(m_IndexOfMostRecentEvent); };
	inline int                                 get_most_recent_event_idx()              const { return m_IndexOfMostRecentEvent; };
	inline int                                 get_total_event_num()                    const { return m_EventSequence.size(); };

	inline DiskCollisionHandler get_collision_handling_type()  const { return m_CollisionHandlingType; };
	inline double get_coefficient_of_restitution()             const { return m_CollisionHandlingType.get_coefficient_of_restitution(); };
	inline bool   get_collision_type()                         const { return m_CollisionHandlingType.get_collision_type(); };

    inline const DiskCollisionWatcher& get_disk_collision_watcher() const { return m_CollisionWatcher; };

    inline const TimeStatistics<milli>&            get_time_statistics(const int& index)        const { return m_TimeStatistics[index]; };
	inline const TimeStatistics<milli>&            get_time_statistics_in_GeneratingEVENTSEQ()  const { return m_TimeStatistics[PHASE1_COMP_TIME]; };
    inline const TimeStatistics<milli>&            get_time_statistics_in_PlayingDVD()          const { return m_TimeStatistics[PHASE2_COMP_TIME]; };
    inline const CountStatistics&                  get_count_statistics_in_GeneratingEVENTSEQ() const { return m_CountStatistics[PHASE1_COUNT]; };
    inline const CountStatistics&                  get_count_statistics(const int& index)       const { return m_CountStatistics[index]; };


	Generator2D* get_generator_by_its_input_disk_id(const int& inputDiskId) const;
	VEdge2D*     get_VEdge_defined_by_these_four_input_disks(
                     int* fourDiskIds, const bool& VEdgeDefinedByQuadrupletIsAlreadyFlipped = false) const;
	VEdge2D*     get_VEdge_with_this_generator_quadruplet(
                     const EdgeFlipEvent::GeneratorPtrQuadruplet& generatorQuadruplet, 
                     const bool& VEdgeDefinedByQuadrupletIsAlreadyFlipped = false) const;

	EdgeFlipEvent::GeneratorPtrQuadruplet get_generator_quadruplet_by_its_id(int* generatorIDs) const;
	EdgeFlipEvent::GeneratorPtrQuadruplet get_generator_quadruplet_by_its_VEdge(VEdge2D* VEdge) const;

    //Setter
    inline void set_current_Voronoi_clock(const double& currentVoronoiClock);
    void        set_parameters(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType);
    void        set_length_of_Voronoi_clock(const double& lengthOfVoronoiClock);
    void        set_collision_handling_type(const DiskCollisionHandler& collisionHandlingType);
    void        set_collision_watcher(const DiskCollisionWatcher& collisionWatcher);
    void        set_offset_of_disks(const double& offsetOfDisk);
    void        set_disk_hopping_probability_matrix(const vector<float>& DHPM);

    void        set_parameters(const double& lengthOfVoronoiClock, 
                               const DiskCollisionWatcher& collisionWatcher,
                               const DiskCollisionHandler& collisionHandlingType);


    //Boolean
    inline bool this_clock_is_future(const double& targetVoronoiClock) const;
    inline bool this_clock_is_past(const double& targetVoronoiClock)   const;


    //Operator
    DynamicVoronoiDiagramCIC& operator=(const DynamicVoronoiDiagramCIC& dynamicVD);


    //Construct VD by dynamic disks.
    void construct_initial_VD(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks, const bool& VEdgeGeometryIsComputed = false);
    void construct_initial_VD(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks, const bool& VEdgeGeometryIsComputed = false);
    void construct_initial_VD(const list<DynamicDisk>& containers, const list<DynamicDisk>& dynamicDisks, const bool& VEdgeGeometryIsComputed = false);
    void construct_initial_VD(const list<DynamicDisk*>& containers, const list<DynamicDisk*>& dynamicDisks, const bool& VEdgeGeometryIsComputed = false);


    //Main function 1: Generate an event history and Go to VD at any target time in the input time horizon. 
    
    void generate_EVENTSEQ(const double& lengthOfVoronoiClock, 
                           const bool& reinstateInitialVD = true);

    void generate_EVENTSEQ(const int& numOfCollisions,
                           const bool& reinstateInitialVD = true);

    void generate_EVENTSEQ(const DiskCollisionWatcher& collisionWatcher,
                           const bool& reinstateInitialVD = true);

    DiskCollisionEvent* generate_EVENTSEQ_by_next_collision(
                           const double& endVoronoiClock, 
                           const bool& reinstateInitialVD = true);

    void rewind_EVENTSEQ(const double& targetVoronoiClock);


    // for temporary.... (not used)
    void generate_EVENTSEQ(const double& lengthOfVoronoiClock, 
                           const DiskCollisionHandler& collisionHandlingType,
                           const bool& reinstateInitialVD = true);

    void generate_EVENTSEQ_collision_avoidance(const double& lengthOfVoronoiClock,
                                               const double& clockRelativeToCollisionClock,
                                               const bool& reinstateInitialVD = true);



    void go_to_target_Voronoi_clock(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed = false);
	void increase_Voronoi_clock(const double& VoronoiClockIncrement = 1.0, const bool& VEdgeGeometryIsComputed = false);
    void go_to_next_event_Voronoi_clock(const int& eventNumIncrement = 1, const bool& VEdgeGeometryIsComputed = false);
    void reinstate_the_initial_Voronoi_diagram();

    //Velocity vector change
    ReservationNumber insert_disk_velocity_change_event_into_priorityQ(const double& eventTime, Generator2D* generator, const rg_Point2D& velocityVector);
    void remove_disk_velocity_change_event_from_priorityQ(const ReservationNumber& reservationNum);
    bool disk_velocity_change_priorityQ_is_empty() const;


    //VD 
    list<Generator2D*> find_neighbor_generators(
        Generator2D* anchorGenerator, const double& startVoronoiClock, const double& endVoronoiClock, const bool& thisIsContainer =false);
    void compute_all_VEdges_equation_of_VD();
    void compute_all_VVertex_coordinates_of_VD();


    //Etc.
    bool validate_DVD_by_empty_circle_test_at_vertices();
    void clear();
    void clear_EVENTSEQ_related_objects();

    //For statistics
    void write_statistics_for_generating_EVENTSEQ(const string& fileNameWithPath, const bool& appendingMode=false) const;
    void write_statistics_in_playing_DVD(const string& fileNameWithPath) const;


    /********************************************************************************************/
    /*           From here, the functions below can be used by one understanding DVD well.      */
    /*            By these funtions, you can apply many functions consist of DVD to your prob.  */
    /********************************************************************************************/
    void set_EVENTSEQ(const vector<EventOfDynamicVD2D*>& eventSequence, const bool& reinstateInitialVD= true);

    void update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_FUTURE(EdgeFlipEvent* edgeFlippingEvent);
    void update_VD_by_EDGE_FLIP_if_NOT_shadow_flip_to_PAST(EdgeFlipEvent* edgeFlippingEvent);
  
    void relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_FUTURE(DiskCollisionEvent* diskCollisionEvent);
    void relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION(DiskCollisionEvent* diskCollisionEvent, const rg_Point2D& newVelocityVector1, const rg_Point2D& newVelocityVector2);
    void relocate_two_disks_and_update_velocity_vectors_by_DISK_COLLISION_to_PAST(DiskCollisionEvent* diskCollisionEvent);

    void relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_FUTURE(DiskVelocityChangeEvent* diskVelocityChangeEvent);
    void relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE(DiskVelocityChangeEvent* diskVelocityChangeEvent,  const rg_Point2D& newVelocityVector);
    void relocate_a_disk_and_update_velocity_vector_by_DISK_VELOCITY_CHANGE_to_PAST(DiskVelocityChangeEvent* diskVelocityChangeEvent);
    
    void relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_FUTURE(DiskHoppingEvent* diskHoppingEvent);
    void relocate_a_disk_btw_containers_by_DISK_HOPPING_or_handle_like_velocity_change_event_by_FAKE_DISK_HOPPING_to_PAST(DiskHoppingEvent* diskHoppingEvent);
    void relocate_a_disk_and_update_velocity_vector_by_FAKE_DISK_HOPPING(DiskHoppingEvent* diskHoppingEvent, const rg_Point2D& newVelocityVector);

    inline EventOfDynamicVD2D* get_next_DVD_event_to_FUTURE()    const;
    inline EventOfDynamicVD2D* get_next_DVD_event_to_PAST()      const;
    inline bool                next_DVD_event_to_FUTURE_exists() const;
    inline bool                next_DVD_event_to_PAST_exists()   const;
    inline bool                next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(const double& targetVoronoiClock) const;
    inline bool                next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(const double& targetVoronoiClock)   const;
    inline bool                next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock() const;

    inline void increment_index_of_most_recent_event() { ++m_IndexOfMostRecentEvent; };
    inline void decrement_index_of_most_recent_event() { --m_IndexOfMostRecentEvent; };

    void   go_back_to_target_Voronoi_clock_and_remove_events_after_the_target_Voronoi_clock(const double& targetVoronoiClock);

    void  remove_disk(DynamicDisk* disk);

    void  synchronize_disk_locations_and_current_clock_to_VD(const double& targetVoronoiClock);

    void ___reset_statistics_of_playing_DVD();
    void ___start_to_track_statistics_of_playing_DVD(const double& targetVoronoiClock);
    void ___stop_to_track_statistics_of_playing_DVD();


private:
    //For V-edge flip event
    pair<bool, double> find_possible_flipping_time_of_an_vedge_after(const double& mostRecentlyEventOccurringTime, const VEdge2D* targetVEdge, const bool& isThisVEdgeFlipping = false);
    void collect_influenced_5VEdges_by_flipping_target_vedge_including_itself(const VEdge2D* targetVEdge, set<VEdge2D*>& influencedVEdgeList) const;
    void find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(const double& mostRecentlyEventOccurringTime, const set<VEdge2D*>& influencedEdges);
    void find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(const VEdge2D* targetVEdge, const double& mostRecentlyEventOccurringTime, const set<VEdge2D*>& influencedEdges);
    void do_bookKeeping_on_priorityQ_for_edge_flipping(const pair<bool, double>& newEdgeFlipTime, const VEdge2D* vEdge, const PriorityQBookkeepingType& PQBookkeepingType);

    void find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
        const double& mostRecentlyEventOccurringTime,
        const set<VEdge2D*>& influencedEdges,
        const DVD_VEdgeStatus& DVD_VEdgeStatus);
    void find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(
        const VEdge2D* targetVEdge,
        const double& mostRecentlyEventOccurringTime,
        const set<VEdge2D*>& influencedEdges,
        const DVD_VEdgeStatus& DVD_VEdgeStatus);
    void do_bookkeeping_on_priorityQ_for_edge_flipping(const pair<bool, double>& newCollisionTime, const VEdge2D* vEdge, const DVD_VEdgeStatus& DVD_VEdgeStatus);

    inline void flip_Vedge_to_FUTURE(VEdge2D* targetVEdge);
    inline void flip_Vedge_to_PAST(VEdge2D* targetVEdge);

    //For disk collision event
    pair<bool, double> find_possible_collision_time_of_disks_after(const double& mostRecentlyEventOccurringTime, const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    void collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& influencedVEdges);
    void collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_PHYSICAL_COLLISION(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& influencedVEdges);
    void collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_AVOIDANCE(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& sharedVEdges, set<VEdge2D*>& theOtherVEdges);
    void find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(const double& mostRecentlyEventOccurringTime, set<VEdge2D*>& influencedVEdges);
    void find_possible_collision_time_of_shared_VEdges_after_input_time_N_store_in_priorityQ_AVOIDANCE(const double& mostRecentlyEventOccurringTime, set<VEdge2D*>& sharedVEdges);
    void do_bookkeeping_on_priorityQ_for_disk_collision(const pair<bool, double>& newCollisionTime, const VEdge2D* vEdge, const PriorityQBookkeepingType& PQBookkeepingType);

    void find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(
        const double& mostRecentlyEventOccurringTime,
        const set<VEdge2D*>& influencedVEdges,
        const DVD_VEdgeStatus& DVD_VEdgeStatus);
    void do_bookkeeping_on_priorityQ_for_disk_collision(const pair<bool, double>& newCollisionTime, const VEdge2D* vEdge, const DVD_VEdgeStatus& VEdgeStatus);


    //For velocity change: this event is from user
    void collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_VEdge_flipping_times(const Generator2D* diskGenerator, set<VEdge2D*>& influencedVEdges);
    void collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_collision_times(const Generator2D* diskGenerator, set<VEdge2D*>& influencedVEdges);
    bool this_number_is_reserved(const ReservationNumber& reservationNum) const;
    ReservationNumber generate_reservation_number_for_velocity_vector_change();


    //For disk hopping
    //void relocate_a_disk_btw_containers_by_DISK_HOPPING(DiskHoppingEvent* diskHoppingEvent, const TimeProgressDirectionType& timeProgressDirection = TO_FUTURE, const bool& duringComputation = true);
    void relocate_a_disk_btw_containers_by_DISK_HOPPING(DiskHoppingEvent* diskHoppingEvent, Generator2D* toContainer, Generator2D* fromContainer, const rg_Point2D& newPosition, const bool& duringComputation = true);
    void relocate_a_disk_btw_containers_by_DISK_HOPPING_to_FUTURE(DiskHoppingEvent* diskHoppingEvent, const bool& duringComputation = true);
    void relocate_a_disk_btw_containers_by_DISK_HOPPING_to_PAST(DiskHoppingEvent* diskHoppingEvent, const bool& duringComputation = true);

    pair<bool, rg_Point2D> find_insert_position_of_disk_in_this_container_using_VVertices(Generator2D* movingDisk, Generator2D* container);
    inline bool            this_is_collision_event_NOT_disk_hopping_event(const VEdge2D* vEdgeForDiskNContainer) const;
    void                   do_bookkeeping_for_disk_hopping(const double& mostRecentlyEventOccurringTime,
                                                           const VEdge2D* vEdgeForDiskNContainer,
                                                           const DynamicDisk* movingDisk,
                                                           const DynamicDisk* container,
                                                           const DVD_VEdgeStatus& VDEdgeStatus,
                                                           const bool& thereExistNextCollisionTime);

#ifdef NEW_CONTAINER_ID_MAPPING
    inline int   convert_container_id_to_transition_probability_index(const int& containerID) const { return -(containerID + 2); };
    inline int   convert_transition_probability_index_to_container_id(const int& index) const { return -(index + 2); };
    inline float get_transition_probability(const int& fromIndex, const int& toIndex) const { return m_DiskHoppingProbabilityMatrix[fromIndex * (m_NumOfContainers-1) + toIndex]; };

#else
    inline int   convert_container_id_to_transition_probability_index(const int& containerID) const { return -(containerID + 1); };
    inline int   convert_transition_probability_index_to_container_id(const int& index) const { return -(index + 1); };
    inline float get_transition_probability(const int& fromIndex, const int& toIndex) const { return m_DiskHoppingProbabilityMatrix[fromIndex * m_NumOfContainers + toIndex]; };

#endif
    Generator2D* choose_next_container_by_transition_probability(const int& currContainerID) const;

    pair<bool, rg_Point2D> check_two_conditions_for_hopping_and_get_inserted_position_if_satisfies_the_conditions(DiskHoppingEvent* diskHoppingEvent);
    void store_velocity_vectors_before_and_after_disk_hopping_by_current_velocity_vector(DiskHoppingEvent* diskHoppingEvent) const;
    void store_current_and_target_positions_for_disk_hopping(DiskHoppingEvent* diskHoppingEvent, const rg_Point2D& targetPosition) const;
    void store_current_position_of_disk_to_event(DiskHoppingEvent* diskHoppingEvent) const;
    void store_velocity_vector_before_and_after_collision(DiskHoppingEvent* diskHoppingEvent, const rg_Point2D& velocityVectorAfterCollding) const;


    //For computing velocity vectors 
    void compute_velocity_vectors_of_disks_after_collided(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2, rg_Point2D& outputVector1, rg_Point2D& outputVector2) const;

    //For disk collision watching
    void watch_and_record_collision_to_FUTURE(DiskCollisionEvent* diskCollisionEvent);
    void watch_and_record_collision_for_group_to_FUTURE(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2);
    void watch_and_record_collision_for_disk_to_FUTURE(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2);
    void watch_and_record_collision_to_PAST(DiskCollisionEvent* diskCollisionEvent);
    void watch_and_record_collision_for_group_to_PAST(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2);
    void watch_and_record_collision_for_disk_to_PAST(DynamicDisk* dynamicDisk1, DynamicDisk* dynamicDisk2);

    void compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(DiskHoppingEvent* diskHoppingEvent) const;
    void compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(DiskCollisionEvent* diskCollisionEvent) const;
    
    void store_current_positions_of_two_colliing_disks_to_event(DiskCollisionEvent* diskCollisionEvent) const;
    void change_positions_of_two_colliding_disks_from_event(DiskCollisionEvent* diskCollisionEvent) const;
   
    void store_current_position_of_disk_to_event(DiskVelocityChangeEvent* diskVelocityChangeEvent) const;
    void change_position_of_disk_from_event(DiskVelocityChangeEvent* diskVelocityChangeEvent) const;

    void change_position_of_disk_from_event(DiskHoppingEvent* diskVelocityChangeEvent) const;

    inline void store_two_disks_velocity_vectors_before_and_after_collision(
                    DiskCollisionEvent* diskCollisionEvent, const rg_Point2D& disk1VelocityVectorAfterCollding, const rg_Point2D& disk2VelocityVectorAfterCollding) const;
    inline void update_velocity_vector(DynamicDisk* dynamicDisk, const rg_Point2D& newVelocityVector) const;
    inline void change_to_new_velocity_vector_and_store_old_one(DiskVelocityChangeEvent* diskVelocityChangeEvent) const;


    //For event queue(priority queue): V-edge flip and disk collision are always tracked by DVD so that we do not pop any events from the Qs.
    EventOfDynamicVD2D*  get_next_event_from_priorityQs();
    EdgeFlipEvent*       generate_next_edge_flip_event_from_priorityQ();
    DiskCollisionEvent*  generate_next_disk_collision_event_from_priorityQ();
    DiskVelocityChangeEvent*  pop_and_generate_next_disk_velocity_change_event_from_priorityQ();
    DiskHoppingEvent*    generate_next_disk_hopping_event_from_priorityQ() const;

    inline EventOfDynamicVD2D::EventType get_next_event_type_in_priorityQ() const;
    inline double        get_next_event_clock_in_priorityQ() const;
    VEdge2D*             get_next_flipping_edge_from_flipping_priorityQ();


    //For moving disks
    void move_disks_to_target_Voronoi_clock_without_changing_VD(set<Generator2D*>& diskGenerators, const double& targetVoronoiClock);

    //For VD queries
    inline void find_4_definining_disks_of_VEdge(set<Generator2D*>& fourDefiningGenerators, VEdge2D* targetVEdge);
    void        get_VVertices_except_infinite_ones(list<VVertex2D*>& VVerticesExceptInfiniteOnes);
    inline bool this_VEdge_is_boundary_of_anomalizing_cell(const VEdge2D* VEdge) const;
    inline int  find_number_of_generator_disks_in_this_container(const Generator2D* container) const;

    inline void compute_circum_circles_Of_VVertices(list<VVertex2D*>& VVertices) { updateCircumcirclesOfNewVVertices(VVertices); };

    void update_geometry_of_VVertex(VVertex2D* VVertex);
    void update_circumcircles_of_VVertices_in_this_container(Generator2D* container);
    void clear_delete_RED_VEdges_and_VVertices();


    //For sorting predicate
    static bool compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    static bool compare_disk_with_ID_in_non_decreasing_order(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    
    //For generating an event
    void handle_and_store_EDGE_FLIP_event(EdgeFlipEvent* edgeFlippingEvent);
    void handle_and_store_DISK_COLLISION_event(DiskCollisionEvent* diskCollisionEvent);
    void handle_and_store_VELOCITY_VECTOR_CHANGE_event(DiskVelocityChangeEvent* diskVelocityChangeEvent);
    void handle_and_store_DISK_HOPPING_event(DiskHoppingEvent* diskHoppingEvent);

	void find_influenced_disks_by_this_FLIP_event(set<Generator2D*>& influencedDiskGenerators, EdgeFlipEvent* anEvent);
	void find_influenced_disks_by_this_DISK_COLLISION_event(set<Generator2D*>& influencedDiskGenerators, DiskCollisionEvent* anEvent);
	void find_influenced_disks_by_this_VELOCITY_VECTOR_CHANGE_event(set<Generator2D*>& influencedDiskGenerators, DiskVelocityChangeEvent* anEvent);
	void find_influenced_disks_by_this_DISK_HOPPING_event(set<Generator2D*>& influencedDiskGenerators, DiskHoppingEvent* anEvent);

	void update_the_event_times_in_the_priorityQs_influenced_by_this_FLIP_event(EdgeFlipEvent* edgeFlippingEvent);
	void update_the_event_times_in_the_priorityQs_influenced_by_this_DISK_COLLISION_event(DiskCollisionEvent* diskCollisionEvent);
	void update_the_event_times_in_the_priorityQs_influenced_by_this_VELOCITY_VECTOR_CHANGE_event(DiskVelocityChangeEvent* diskVelocityChangeEvent);
	void update_the_event_times_in_the_priorityQs_influenced_by_this_DISK_HOPPING_event(DiskHoppingEvent* diskHoppingEvent);

	void store_this_FLIP_event_in_the_event_history_vector(EdgeFlipEvent* edgeFlippingEvent);
	void store_this_DISK_COLLISION_event_in_the_event_history_vector(DiskCollisionEvent* diskCollisionEvent);
	void store_this_VELOCITY_VECTOR_CHANGE_event_in_the_event_history_vector(DiskVelocityChangeEvent* diskVelocityChangeEvent);
	void store_this_DISK_HOPPING_event_in_the_event_history_vector(DiskHoppingEvent* diskHoppingEvent);

    //For handling special case where upto 5 V-edges can have an identical flip time.
    void clear_info_about_vedges_defined_by_same_disks();

    //For clearing V-edges influenced by current disk hopping event.
    void clear_info_about_vedges_related_to_disk_hopping();


    //For playing
    void handle_next_DVD_event_to_FUTURE();
    void reversely_handle_next_DVD_event_to_PAST();

    inline void set_index_of_most_recent_event(const int& index) { m_IndexOfMostRecentEvent = index; };

    //For etc.f
	void link_dynamic_disk_id_to_its_generator_pointer();
    void convert_input_for_constructing_VD(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks, list< pair< rg_Circle2D, void* > >& diskNUserDataPairs);
    void convert_input_for_constructing_VD(const list<DynamicDisk*>& containers, const list<DynamicDisk*>& dynamicDisks, list< pair< rg_Circle2D, void* > >& diskNUserDataPairs);

    void copy_initial_VD_and_containers_and_disks_for_reinstatement();
	void copy_initial_VD_and_containers_and_disks_for_reinstatement(const VoronoiDiagramCIC& initialVD, DynamicDisk* const inputContainer, const list<DynamicDisk*>& inputDisks);
    void copy_initial_VD_and_containers_and_disks_for_reinstatement(const VoronoiDiagramCIC& initialVD, const list<DynamicDisk*>& inputContainers, const list<DynamicDisk*>& inputDisks);
    void copy_initial_VD_and_make_map_from_copied_generators_to_original_generators(const VoronoiDiagramCIC& initialVD);
    void make_map_from_original_disks_to_copied_disks(const list<DynamicDisk*>& inputContainers, const list<DynamicDisk*>& inputDisks);
    
    void link_Voronoi_diagrams(const list<VoronoiDiagramCIC*>& allVDs_Copied, const list<VoronoiDiagramCIC*>& allVDs_Current, unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromCopiedToCurrent);
    void delete_all_VEntities_in_lists_and_clear_the_lists_in_hierarchy(const list<VoronoiDiagramCIC*>& allVDsInHierarchy);
    void clear_all_generator_lists_in_hierarchy_while_the_dynamic_memory_of_the_generators_remains(const list<VoronoiDiagramCIC*>& allVDsInHierarchy);

    void initialize_parameters();

    void get_dynamic_disks_pointers(const list<DynamicDisk>& dynamicDisks, list<DynamicDisk*>& dynamicDiskPointers);
    void propagate_DVD_to_the_last_moment_of_Voronoi_clock();
    void propagate_DVD_until_target_collision_numbers_btw_disks(const int& numOfCollisions);
    void propagate_DVD_until_target_collision_numbers_of_watch_objects();
    DiskCollisionEvent* propagate_DVD_to_the_next_collision(const double& endVoronoiClock);

    void propagate_DVD_by_handling_and_storing_next_event(EventOfDynamicVD2D* nextEvent);

    void construct_initial_priority_Qs_for_both_flipping_and_collision_events();
    void construct_initial_priorityQ_of_vedge_flipping_time();
    void construct_initial_priorityQ_of_disk_collision_time();
    
    void construct_initial_priority_queues();
    void construct_initial_priorityQ_of_vedge_flipping_time_vers2();
    void construct_initial_priorityQ_of_disk_collision_time_and_mark_hopping_disk_among_them();

    void reinstate_the_initial_Voronoi_diagram_after_generating_EVENTSEQ();
    void reinstate_the_initial_VD_topology();
    void reinstate_the_initial_disk_locations_and_velocity_vectors();
    void clear_current_collision_numbers_of_interest();
    void clear_priority_queues_for_collision_and_flip();

    bool this_VEdge_is_in_this_set(const VEdge2D* VEdge, const set<VEdge2D*>& VEdgeSet) const;
    inline void get_4Disks_defining_this_VEdge_and_its_end_points(const VEdge2D* VEdge, vector<DynamicDisk*>& disksForVEdge);
    inline void get_4DiskGenerators_defining_this_VEdge_and_its_end_points(const VEdge2D* VEdge, vector<Generator2D*>& generatorsForVEdge) const;

    //For collision avoiding ( temp...)
    tuple<bool, double, DynamicDisk*, DynamicDisk*> propagate_DVD_until_next_collision();
    double     get_next_target_time_to_go_back(const double& timeRelativeToCollisionTime, const double& currentVoronoiClock) const;
    
    rg_Point2D rotate_velocity_vector(const rg_Point2D& inputVector, const double& radian) const;

    rg_Point2D compute_first_disk_velocity_vector_for_collision_avoidance(DynamicDisk* disk1,  DynamicDisk* disk2) const;
    tuple<DynamicDisk*, rg_Point2D, DynamicDisk*, rg_Point2D> compute_velocity_vectors_by_brute_force_search_for_collision_avoidance(DynamicDisk* disk1, DynamicDisk* disk2, const double& currentVoronoiClock);

    DiskVelocityChangeEvent* create_velocity_vector_change_event(const double& eventTime, DynamicDisk* disk, const rg_Point2D& newVelocityVector) const;


	//For generating polynomial of V-edge flip and solving the polynomial
    rg_Polynomial      make_polynomial_equation_for_vedge_flipping(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);
    void               adjust_polynomial_degree(const int& newDegree, rg_Polynomial& polynomial);
    pair<bool, double> solve_polynomial_equation_for_vedge_flipping(
                          const rg_Polynomial& polynomialEq,
                          const double& mostRecentlyEventOccurringTime,
                          const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder,
                          const bool& isThisVEdgeFlipping);

    void sort_generators_in_non_increasing_order(vector<DynamicDisk*>& disks);
   
    bool this_VEdge_is_NOT_shadow_flipping(VEdge2D* VEdge, const double& eventTime, const rg_REAL resolution = rg_MATH_RES);
    int  compute_polynomial_degree_by_disk_info(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);
    


    //For selecting a V-edge to be flipped.
    VEdge2D* find_next_flipping_vedge_when_vedges_have_same_flipping_time();
    void find_all_VEdges_defined_by_identical_disks(const VEdge2D* VEdge, const list<VEdge2D*>& candidates, set<VEdge2D*>& outputs);
    void find_all_VEdges_defined_by_identical_disks_TOPOLOGICAL_COMPUTATION(const VEdge2D* VEdge, const set<VEdge2D*>& candidates, set<VEdge2D*>& outputs) const;


    //For playing DVD given event history 
    void construct_VD_at_target_Voronoi_clock(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed = false);
    void progress_DVD_using_EVENTNUM_to_FUTURE(const int& eventNumIncrement);
    void progress_DVD_using_EVENTNUM_to_PAST(const int& eventNumIncrement);
    void progress_DVD_to_FUTURE(const double& targetVoronoiClock);
    void progress_DVD_to_PAST(const double& targetVoronoiClock);
   
    inline void relocate_two_disks_at_colliding_time(const double& collidingTime, Generator2D* diskGenerator1, Generator2D* diskGenerator2, const rg_Point2D& disk1CenterAtColldingTime, const rg_Point2D& disk2CenterAtColldingTime);
    inline void relocate_a_disk_at_velocity_vector_change_time(const double& velocityChangeTime, Generator2D* diskGenerator, const rg_Point2D& diskCenterAtVelocityChangeTime);

    inline void update_two_disks_velocity_vectors_by_collision(Generator2D* diskGenerator1, Generator2D* diskGenerator2, const rg_Point2D& newVelocityVector1, const rg_Point2D& newVelocityVector2);
    inline void update_a_disk_velocity_vector(Generator2D* diskGenerator, const rg_Point2D& newVelocityVector);
    void        move_disks_related_to_disk_hopping_event(DiskHoppingEvent* diskHoppingEvent);
    void        move_disks_and_update_VVertices_related_to_disk_hopping_event_to_FUTURE(DiskHoppingEvent* diskHoppingEvent);
    void        move_disks_and_update_VVertices_related_to_disk_hopping_event_to_PAST(DiskHoppingEvent* diskHoppingEvent);
    void        move_disks_and_update_VVertices_related_to_disk_hopping_event(DiskHoppingEvent* diskHoppingEvent, Generator2D* containerIncludingVVerticesToBeUpdated);

    //For dynamicDisk
    void          get_dynamic_disks_and_containers(list<DynamicDisk*>& dynamicDisksNContainers) const;
    void          get_dynamic_disks(list<DynamicDisk*>& dynamicDisks)                           const;
    void          get_containers(list<DynamicDisk*>& containers)                                const;
    DynamicDisk*  get_most_outer_container() const;


    //For copying DVD
    void copy_from(const DynamicVoronoiDiagramCIC& fromThisDVD);
    void copy_input(const DynamicVoronoiDiagramCIC& fromThisDVD);
    void copy_DVD_current_state_variables(const DynamicVoronoiDiagramCIC& fromThisDVD);
 
    void copy_time_statistics(const DynamicVoronoiDiagramCIC& fromThisDVD);
    void copy_count_statistics(const DynamicVoronoiDiagramCIC& fromThisDVD);
  
    void create_event_sequence(const DynamicVoronoiDiagramCIC& fromThisDVD);

    void make_link_btw_VEdges(list<VEdge2D*>& fromVEdges, list<VEdge2D*>& toVEdges, map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD);
    void make_link_btw_generators_from_other_VD_to_this_VD(const DynamicVoronoiDiagramCIC& fromThisVD, map<Generator2D*, Generator2D*>& generatorLinkFromOtherVDToThisVD);
    void make_link_btw_disks_and_containers_from_other_VD_to_this_VD(const DynamicVoronoiDiagramCIC& fromThisVD, map<DynamicDisk*, DynamicDisk*>& diskLinkFromOtherVDToThisVD);
    void modify_user_data_pointing_to_disk_in_generators(map<DynamicDisk*, DynamicDisk*>& diskLinkFromOtherVDToThisVD);
    void modify_generators_in_event_sequence(const map<Generator2D*, Generator2D*>& generatorLinkFromOtherVDToThisVD);
    void modify_VEdges_in_priorityQs_from_other_VDs_one_to_this_VDs_one(const map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD);

    void clear_EVENTSEQ();
    void clear_time_statistics();
    void clear_count_statistics();
    void clear_event_generation_related_objects();
    void clear_VD_reinstatement_objects();
    void clear_input();
    void clear_mapper_from_id_to_generator();

    //For VD queries - special case detection
    bool these_two_VEdge_are_defined_by_two_same_disks(VEdge2D* VEdge1, VEdge2D* VEdge2) const;
    bool this_VEdge_and_its_hands_and_legs_are_in_this_set(VEdge2D* VEdge, set<VEdge2D*>& VEdgeSet) const;
    bool these_three_VEdge_have_common_disk_in_right_or_left_VFace(const set<VEdge2D*>& VEdges) const;
    void find_VEdges_defined_by_identical_4disks_from_PQ_N_flipping_VEdge(const VEdge2D* flippingVEdge,
                                                                          const set<VEdge2D*> upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ,
                                                                          set<VEdge2D*>& VEdgesDefinedByIdentical4Disk) const;

    //For EVENTSEQ 
    void compute_velocity_vectors_by_scanning_EVENTSEQ();
    void remove_events_in_EVENTSEQ_after_current_Voronoi_clock();


//public:
    //For statistics
    void prepare_statistics_in_generating_EVENTSEQ();
    void complete_statistics_in_generating_EVENTSEQ(const bool& phase1WillBeContinued = false);

    void initialize_statistics_for_collision_avoidance();
    void finalize_statistics_for_collision_avoidance();

    inline void ___start_clock_in_GeneratingEVENTSEQ(const unsigned int& it){ if (m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ) { m_TimeStatistics[PHASE1_COMP_TIME].start_clock(it); } };
    inline void ___end_clock_in_GeneratingEVENTSEQ(const unsigned int& it)  { if (m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ) { m_TimeStatistics[PHASE1_COMP_TIME].end_clock(it); } };
    inline void ___count_in_GeneratingEVENTSEQ(const unsigned int& it)      { if (m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ) { m_CountStatistics[PHASE1_COUNT].addCount(it); } };

    inline void ___start_clock_in_phase2(const unsigned int& typeOfStatistics) {  m_TimeStatistics[typeOfStatistics].start_clock(); };
    inline void ___end_clock_in_phase2(const unsigned int& typeOfStatistics)   { m_TimeStatistics[typeOfStatistics].end_clock();  };
    inline void ___flipping_count_in_phase2()        {  m_CountStatistics[PHASE2_FLIPPING_COUNT].addCount(); };
    inline void ___collision_count_in_phase2()       {  m_CountStatistics[PHASE2_COLLISION_COUNT].addCount(); };
	inline void ___collision_btw_disks_count_in_phase2() { m_CountStatistics[PHASE2_COLLISION_BTW_DISKS_COUNT].addCount(); };
    inline void ___velocity_change_count_in_phase2() {  m_CountStatistics[PHASE2_VELOCITY_CHANGE_COUNT].addCount(); };
	inline void ___disk_hopping_in_phase2()          {  m_CountStatistics[PHASE2_DISK_HOPPING_COUNT].addCount(); };

	inline void ___start_clock_in_collision_avoidance(const unsigned int& it) {  m_TimeStatistics[PHASE1_COLLISION_AVOIDANCE_TIME].start_clock(it);  };
	inline void ___end_clock_in_collision_avoidance(const unsigned int& it)    {  m_TimeStatistics[PHASE1_COLLISION_AVOIDANCE_TIME].end_clock(it);  };

    inline void        ___check_polynomial_degree_for_statistics(const rg_Polynomial& polynomial);
    inline void        ___check_number_of_VEdges_defined_by_identical_4disks_for_statistics();
    inline void        ___check_flipping_type_of_VEdge_for_statistics(EdgeFlipEvent* edgeFlippingEvent);

    inline int ___count_collisions_btw_disks() const;

#ifdef DEBUG_DVD
    const string debugFileNameWithPath = "DVD-debug.txt";

    void ___initialize_debug() const;
    void ___write_current_event_info(EventOfDynamicVD2D* anEvent) const;
    void ___write_update_event_info(const EventOfDynamicVD2D::EventType& eventType, VEdge2D* VEdge, const double newEventTime) const;
#endif //DEBUG_DVD

#ifdef WRITING_EVENT_TIME
    const string eventTimeFileNameWithPath         = "event_times.txt";
    static double accumulatedTimeInEventGeneration = 0.0;

    void ___write_event_times(const double& eventTime, const double& computationTime);
    void ___initialize_writing_event_time();

#endif //WRITING_EVENT_TIME

};


inline void DynamicVoronoiDiagramCIC::set_current_Voronoi_clock(const double& currentVoronoiClock)
{ 
    m_CurrentVoronoiClock = currentVoronoiClock; 
};


inline bool DynamicVoronoiDiagramCIC::this_is_collision_event_NOT_disk_hopping_event(const VEdge2D * vEdgeForDiskNContainer) const
{
    if (  m_InThisContainerDiskWillBeInserted.find(const_cast<VEdge2D*> (vEdgeForDiskNContainer)) 
       == m_InThisContainerDiskWillBeInserted.end())
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline int DynamicVoronoiDiagramCIC::find_number_of_generator_disks_in_this_container(const Generator2D * container) const
{
    list<Generator2D*> generatorDisksInContainer;
    container->getInnerGens(generatorDisksInContainer);

    return generatorDisksInContainer.size();
}


inline void DynamicVoronoiDiagramCIC::clear_info_about_vedges_related_to_disk_hopping()
{
    m_NewVEdgesByDiskHopping.clear();
    m_DeletedVEdgesByDiskHopping.clear();
    m_RedefinedVEdgesByDiskHopping.clear();
}


inline void DynamicVoronoiDiagramCIC::flip_Vedge_to_FUTURE(VEdge2D* targetVEdge)
{
    targetVEdge->flip();
}


inline void DynamicVoronoiDiagramCIC::flip_Vedge_to_PAST(VEdge2D* targetVEdge)
{
    targetVEdge->flip_reverse();
}


inline void DynamicVoronoiDiagramCIC::change_to_new_velocity_vector_and_store_old_one(DiskVelocityChangeEvent* diskVelocityChangeEvent) const
{
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskVelocityChangeEvent->get_disk_generator()->getUserData());
    
    diskVelocityChangeEvent->set_velocity_vector_before_change(movingDisk->getVelocityVector());
    diskVelocityChangeEvent->set_velocity_vector_setting(true);

    update_velocity_vector(movingDisk, diskVelocityChangeEvent->get_future_velocity_vector_of_disk());
}



inline void DynamicVoronoiDiagramCIC::store_two_disks_velocity_vectors_before_and_after_collision(DiskCollisionEvent* diskCollisionEvent,
                                                                                                               const rg_Point2D& disk1VelocityVectorAfterCollding,
                                                                                                               const rg_Point2D& disk2VelocityVectorAfterCollding) const
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    diskCollisionEvent->set_velocity_vectors_before_collision(movingDisk1->getVelocityVector(), movingDisk2->getVelocityVector());
    diskCollisionEvent->set_velocity_vectors_after_collision(disk1VelocityVectorAfterCollding, disk2VelocityVectorAfterCollding);
    diskCollisionEvent->set_velocity_vector_setting(true);
}


inline void DynamicVoronoiDiagramCIC::update_velocity_vector(DynamicDisk* dynamicDisk, const rg_Point2D& newVelocityVector) const
{
    dynamicDisk->setVelocityVector(newVelocityVector);
}

inline bool DynamicVoronoiDiagramCIC::this_clock_is_future(const double& targetVoronoiClock) const
{
    if (targetVoronoiClock > get_current_Voronoi_clock())
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool DynamicVoronoiDiagramCIC::this_clock_is_past(const double& targetVoronoiClock) const
{
    if (targetVoronoiClock < get_current_Voronoi_clock())
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline EventOfDynamicVD2D* DynamicVoronoiDiagramCIC::get_next_DVD_event_to_FUTURE() const
{
    return m_EventSequence[m_IndexOfMostRecentEvent + 1];
}


inline EventOfDynamicVD2D* DynamicVoronoiDiagramCIC::get_next_DVD_event_to_PAST() const
{
    return m_EventSequence[m_IndexOfMostRecentEvent];
}


inline bool DynamicVoronoiDiagramCIC::next_DVD_event_to_FUTURE_exists() const
{
    const int indexOfNextEventToFuture = m_IndexOfMostRecentEvent + 1;

    if (indexOfNextEventToFuture < m_EventSequence.size())
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool DynamicVoronoiDiagramCIC::next_DVD_event_to_PAST_exists() const
{
    const int indexOfNextEventToPast = m_IndexOfMostRecentEvent;

    if (indexOfNextEventToPast >= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool DynamicVoronoiDiagramCIC::next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(const double& targetVoronoiClock) const
{
    if (!next_DVD_event_to_FUTURE_exists())
    {
        return false;
    }

    EventOfDynamicVD2D* nextFutureEvent = get_next_DVD_event_to_FUTURE();
    if (nextFutureEvent->get_occurring_clock() <= targetVoronoiClock)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool DynamicVoronoiDiagramCIC::next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(const double& targetVoronoiClock) const
{
    if (!next_DVD_event_to_PAST_exists())
    {
        return false;
    }

    EventOfDynamicVD2D* nextPastEvent = get_next_DVD_event_to_PAST();
    if (nextPastEvent->get_occurring_clock() > targetVoronoiClock)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool DynamicVoronoiDiagramCIC::next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock() const
{
    if (!next_DVD_event_to_PAST_exists())
    {
        return false;
    }

    EventOfDynamicVD2D* nextPastEvent = get_next_DVD_event_to_PAST();
    if (nextPastEvent->get_occurring_clock() == 0.0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline void DynamicVoronoiDiagramCIC::relocate_a_disk_at_velocity_vector_change_time(const double& velocityChangeTime, Generator2D* diskGenerator, const rg_Point2D& diskCenterAtVelocityChangeTime)
{
    DynamicDisk* dynamicDisk = static_cast<DynamicDisk*>(diskGenerator->getUserData());

    dynamicDisk->setMostRecentLocationUpdateTime(velocityChangeTime);

    rg_Circle2D diskAtVelocityChangeTime(diskCenterAtVelocityChangeTime, dynamicDisk->getCircle().getRadius());

    dynamicDisk->setCircle(diskAtVelocityChangeTime);
    diskGenerator->setDisk(diskAtVelocityChangeTime);
}



inline void DynamicVoronoiDiagramCIC::relocate_two_disks_at_colliding_time(const double& collidingTime,
                                                                          Generator2D* diskGenerator1, Generator2D* diskGenerator2,
                                                                          const rg_Point2D& disk1CenterAtColldingTime, const rg_Point2D& disk2CenterAtColldingTime)
{
    DynamicDisk* disk1 = static_cast<DynamicDisk*>(diskGenerator1->getUserData());
    DynamicDisk* disk2 = static_cast<DynamicDisk*>(diskGenerator2->getUserData());

    disk1->setMostRecentLocationUpdateTime(collidingTime);
    disk2->setMostRecentLocationUpdateTime(collidingTime);

    rg_Circle2D disk1AtCollidingTime(disk1CenterAtColldingTime, disk1->getCircle().getRadius());
    rg_Circle2D disk2AtCollidingTime(disk2CenterAtColldingTime, disk2->getCircle().getRadius());

    disk1->setCircle(disk1AtCollidingTime);
    disk2->setCircle(disk2AtCollidingTime);

    diskGenerator1->setDisk(disk1AtCollidingTime);
    diskGenerator2->setDisk(disk2AtCollidingTime);
}


inline void DynamicVoronoiDiagramCIC::update_two_disks_velocity_vectors_by_collision(Generator2D* diskGenerator1, Generator2D* diskGenerator2, const rg_Point2D& newVelocityVector1, const rg_Point2D& newVelocityVector2)
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskGenerator1->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskGenerator2->getUserData());

    movingDisk1->setVelocityVector(newVelocityVector1);
    movingDisk2->setVelocityVector(newVelocityVector2);
}


inline void DynamicVoronoiDiagramCIC::update_a_disk_velocity_vector(Generator2D* diskGenerator, const rg_Point2D& newVelocityVector)
{
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskGenerator->getUserData());
    movingDisk->setVelocityVector(newVelocityVector);
}


inline void DynamicVoronoiDiagramCIC::get_4Disks_defining_this_VEdge_and_its_end_points(const VEdge2D* VEdge, vector<DynamicDisk*>& disksForVEdge)
{
    ___start_clock_in_GeneratingEVENTSEQ(DVD_TIME_GET_4DISKS_DEFINING_VEDGE);

    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getRightFace()->getGenerator()->getUserData()));
    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getLeftFace()->getGenerator()->getUserData()));
    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getMateFace(VEdge->getStartVertex())->getGenerator()->getUserData()));
    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getMateFace(VEdge->getEndVertex())->getGenerator()->getUserData()));

    ___end_clock_in_GeneratingEVENTSEQ(DVD_TIME_GET_4DISKS_DEFINING_VEDGE);

}


inline EventOfDynamicVD2D::EventType DynamicVoronoiDiagramCIC::get_next_event_type_in_priorityQ() const
{
    EventOfDynamicVD2D::EventType minKeyType = EventOfDynamicVD2D::EDGE_FLIP;
    double    minKey     = DBL_MAX;

    if (!m_PriorityQForEdgeFlipping.empty())
    {
        if (m_PriorityQForEdgeFlipping.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2D::EDGE_FLIP;
            minKey = m_PriorityQForEdgeFlipping.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskCollision.empty())
    {
        if (m_PriorityQForDiskCollision.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2D::DISK_COLLISION;
            minKey = m_PriorityQForDiskCollision.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskVelocityChange.empty())
    {
        if (m_PriorityQForDiskVelocityChange.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2D::DISK_VELOCITY_CHANGE;
            minKey = m_PriorityQForDiskVelocityChange.topNode()->getKey();
        }
    }


    return minKeyType;
}


inline double DynamicVoronoiDiagramCIC::get_next_event_clock_in_priorityQ() const
{
    EventOfDynamicVD2D::EventType minKeyType = EventOfDynamicVD2D::EDGE_FLIP;
    double    minKey     = DBL_MAX;

    if (!m_PriorityQForEdgeFlipping.empty())
    {
        if (m_PriorityQForEdgeFlipping.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2D::EDGE_FLIP;
            minKey = m_PriorityQForEdgeFlipping.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskCollision.empty())
    {
        if (m_PriorityQForDiskCollision.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2D::DISK_COLLISION;
            minKey = m_PriorityQForDiskCollision.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskVelocityChange.empty())
    {
        if (m_PriorityQForDiskVelocityChange.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2D::DISK_VELOCITY_CHANGE;
            minKey = m_PriorityQForDiskVelocityChange.topNode()->getKey();
        }
    }

    return minKey;
}


inline void DynamicVoronoiDiagramCIC::get_4DiskGenerators_defining_this_VEdge_and_its_end_points(const VEdge2D * VEdge, vector<Generator2D*>& generatorsForVEdge ) const
{
    generatorsForVEdge.push_back(VEdge->getRightFace()->getGenerator());
    generatorsForVEdge.push_back(VEdge->getLeftFace()->getGenerator());
    generatorsForVEdge.push_back(VEdge->getMateFace(VEdge->getStartVertex())->getGenerator());
    generatorsForVEdge.push_back(VEdge->getMateFace(VEdge->getEndVertex())->getGenerator());
}



inline void DynamicVoronoiDiagramCIC::find_4_definining_disks_of_VEdge(set<Generator2D*>& fourDefiningGenerators, VEdge2D* targetVEdge)
{
    Generator2D* rightGenerator = targetVEdge->getRightFace()->getGenerator();
    Generator2D* leftGenerator = targetVEdge->getLeftFace()->getGenerator();
    Generator2D* headGenerator = targetVEdge->getMateFace(targetVEdge->getEndVertex())->getGenerator();
    Generator2D* legGenerator = targetVEdge->getMateFace(targetVEdge->getStartVertex())->getGenerator();

    fourDefiningGenerators.insert(rightGenerator);
    fourDefiningGenerators.insert(leftGenerator);
    fourDefiningGenerators.insert(headGenerator);
    fourDefiningGenerators.insert(legGenerator);
}



inline bool DynamicVoronoiDiagramCIC::this_VEdge_is_boundary_of_anomalizing_cell(const VEdge2D* VEdge) const
{
    Generator2D* headGenerator = VEdge->getMateFace(VEdge->getStartVertex())->getGenerator();
    Generator2D* tailGenerator = VEdge->getMateFace(VEdge->getEndVertex())->getGenerator();

    if (headGenerator == tailGenerator)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline void DynamicVoronoiDiagramCIC::___check_number_of_VEdges_defined_by_identical_4disks_for_statistics()
{
    if (!m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ)
    {
        return;
    }

    switch (m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.size())
    {
    case 2:
        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_2_IDENTICAL_FLIPPING_TIME);
        break;

    case 3:
        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_3_IDENTICAL_FLIPPING_TIME);
        break;

    case 4:
        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_4_IDENTICAL_FLIPPING_TIME);
        break;

    case 5:
        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_5_IDENTICAL_FLIPPING_TIME);
        break;
    }
}



inline void DynamicVoronoiDiagramCIC::___check_flipping_type_of_VEdge_for_statistics(EdgeFlipEvent* edgeFlippingEvent)
{
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        switch (m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.size())
        {
        case 1:
        {            
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_ORDINARY_FLIPPING);
            break;
        }

        case 2:
        {      
            VEdge2D* VEdge1 = *m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin();
            VEdge2D* VEdge2 = *(++m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin());

            if (these_two_VEdge_are_defined_by_two_same_disks(VEdge1, VEdge2))
            {
                ___count_in_GeneratingEVENTSEQ(DVD_COUNT_24_TO_33_FLIPPING);
            }

            break;
        }

        case 3:
        {
            if (these_three_VEdge_have_common_disk_in_right_or_left_VFace(m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ))
            {
                ___count_in_GeneratingEVENTSEQ(DVD_COUNT_3_TO_2_FLIPPING);
            }

            break;
        }

        case 5:
        {
            VEdge2D* flippingVEdge = get_VEdge_with_this_generator_quadruplet(edgeFlippingEvent->get_generator_quadruplet_before_flip());

            if (this_VEdge_and_its_hands_and_legs_are_in_this_set(flippingVEdge, m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ))
            {
                ___count_in_GeneratingEVENTSEQ(DVD_COUNT_33_TO_22_FLIPPING);
            }
            else
            {
                ___count_in_GeneratingEVENTSEQ(DVD_COUNT_33_TO_24_FLIPPING);
            }

            break;
        }
        }
    }
    else
    {
        ___count_in_GeneratingEVENTSEQ(DVD_COUNT_SHADOW_FLIPPING);
    }
}



inline void DynamicVoronoiDiagramCIC::___check_polynomial_degree_for_statistics(const rg_Polynomial& polynomial)
{
    if (!m_StatisticsTrackingIsOn_Phase_GeneratingEVENTSEQ)
    {
        return;
    }

    switch (polynomial.getDegree())
    {
        case 0:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE0_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 1:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE1_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 2:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE2_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 3:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE3_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 4:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE4_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 5:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE5_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 6:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE6_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 7:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE7_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        case 8:
            ___count_in_GeneratingEVENTSEQ(DVD_COUNT_DEGREE8_POLYNOMIAL_FOR_EDGE_FLIPPING);
            break;

        default:
            break;
    }
}


inline int DynamicVoronoiDiagramCIC::___count_collisions_btw_disks() const
{
    return m_CountStatistics[PHASE1_COUNT].count(DVD_COUNT_DISK_COLLISION_BTW_INPUT_DISKS_NOT_WITH_CONTAINER_IN_EVENT_HISTORY);
}




}
}

#endif











