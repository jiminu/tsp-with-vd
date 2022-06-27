#ifndef _DYNAMIC_VORONOIDIAGRAM_2D_
#define _DYNAMIC_VORONOIDIAGRAM_2D_

//#define USING_BUCKET_FOR_CHECKING_ALL_PAIR_WISE_DISK_INTERSECTION
#define WRITING_EVENT_HISTORY_FILE_BY_HEXA_FLIPPING_TIME
//#define TEST_FOR_TOPOLOGY_RESTORATION
//#define WRITING_POLYNOMIAL_RELATED_INFORMATION
//#define WRITING_EVENT_TIME
//#define DEBUG_DVD

#include "VoronoiDiagramCIC.h"
#include "EventOfDynamicVD2DC.h"
#include "EdgeFlipEvent.h"
#include "DiskCollisionEvent.h"
#include "DiskVelocityChangeEvent.h"
#include "EntityAccessiblePriorityQ.h"
#include "DynamicDisk.h"
#include "TimeStatistics.h"
#include "CountStatistics.h"
#include "rg_Polynomial.h"
#include "rg_RelativeOp.h"
#include "rg_Triplet.h"

#include <list>
#include <array>
#include <vector>
#include <map>
using namespace std;


namespace V {
namespace GeometryTier {


class DynamicVoronoiDiagram2DC : public VoronoiDiagramCIC
{
public: 
    template <class T>
    static inline void hash_combine(std::size_t& seed, const T& v)
    {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    struct Pointer_N_Point2D_pair_hash {
        template <class T>
        std::size_t operator () (const std::pair<T, rg_Point2D>&p) const {
            size_t seed = 0;
            DynamicVoronoiDiagram2DC::hash_combine(seed, p.first);
            DynamicVoronoiDiagram2DC::hash_combine(seed, p.second.getX());
            DynamicVoronoiDiagram2DC::hash_combine(seed, p.second.getY());
            return seed;
        }
    };

    struct Pointer_N_Point2D_pair_Equality_Pred {
        template <class T>
        bool operator () (const std::pair<T, rg_Point2D>&p, const std::pair<T, rg_Point2D>&q) const 
        {
            if (p.first == q.first && p.second == q.second)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    };

    template<typename T> using EntityAccessiblePriorityQ_ = EntityAccessiblePriorityQ<T, less<PQ_Node<T>>, Pointer_N_Point2D_pair_hash, Pointer_N_Point2D_pair_Equality_Pred>;

    enum EdgeFlippingFormulaGeneratingType { FOR_FINDING_TIME, FOR_VERIFING_TIME };
    enum PriorityQBookkeepingType { DURING_INITIAL_PQ_CONSTRUCTION, AFTER_PQ_CONSTRUCTION };
    enum SimulationDirectionType { TO_FUTURE, TO_PAST, NONE_DIRECTION };
    enum TimeStatisticsType { PHASE1_COMP_TIME, PHASE2_COMP_TIME, PHASE2_MOVING_DISKS_TIME, PHASE2_UPDATING_TOPOLOGY_TIME, PHASE2_UPDATING_GEOMETRY_TIME, PHASE2_SIMULATION_TIME_INCREMNT, TIME_STATISTICS_SIZE };
    enum CountStatisticsType { PHASE1_COUNT, PHASE2_FLIPPING_COUNT, PHASE2_COLLISION_COUNT, COUNT_STATISTICS_SIZE };
    enum PolynomialSolverType { LABSRC, JENKINS_TRAUB, BBPR, STURMS_SEQ };
    enum CollisionType { PHYSICAL_COLLISION, INITIALLY_MOVED_OR_FIXED_COLLISION, AVOIDANCE };
    enum DynamicVDAlgorithmPhase1
    {
        DVD_TIME_PHASE_ONE_TOTAL,
        DVD_TIME_CONSTRUCT_INITIAL_CIC_VD,
        DVD_TIME_CHECK_WHETHER_DISKS_ARE_INTERSECTED,
        DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_QS_FOR_FLIPPING_AND_COLLISION_EVENTS,
        DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_VEDGE_FLIPPING_TIME,
        DVD_TIME_CONSTRUCT_INITIAL_PRIORITY_Q_OF_DISK_COLLISION_TIME,
        DVD_TIME_GENERATE_ALL_EVENTS_DURING_SIMULATION_TIME_WINDOW,
        DVD_TIME_MAKE_VD_N_DISKS_BE_AT_TIME_ZERO,
        /*Below this line, these are not directly related to the algorithm steps*/
        DVD_TIME_FIND_POSSIBLE_FLIPPING_TIME_OF_EDGE,
        DVD_TIME_GET_4DISKS_DEFINING_VEDGE,
        DVD_TIME_SORTING_GENERATORS_IN_NON_INCRESING_ORDER,
        DVD_TIME_MAKING_EDGE_FLIPPING_POLYNOMIAL_EQ,
        DVD_TIME_SOLVING_EDGE_FLIPPING_POLYNOMIAL_EQ,
        DVD_TIME_FIND_POSSIBLE_COLLISION_TIME_OF_DISKS,
        DVD_TIME_PRIORITYS_BOOKKEEPING,
        DVD_TIME_PRIORITY_BOOKKEEPING_FOR_EDGE_FLIPPING,
        DVD_TIME_PRIORITY_BOOKKEEPING_FOR_DISK_COLLISION,
        DVD_TIME_DISK_MOVING,
        DVD_TIME_WRITE_POLYNOMIALS_TO_FILE_IN_MAIN_FUNC,
        DVD_TIME_WRITE_POLYNOMIALS_TO_FILE_OUT_OF_MAIN_FUNC,
        DVD_TIME_SIZE_PHASE_ONE
    };
    enum DynamicVDInterestingSubject
    {
        DVD_COUNT_EDGE_FLIPPING_IN_EVENT_HISTORY,
        DVD_COUNT_DISK_COLLISION_IN_EVENT_HISTORY,
        DVD_COUNT_DISK_COLLISION_BTW_INPUT_DISKS_NOT_WITH_CONTAINER_IN_EVENT_HISTORY,
        DVD_COUNT_DISK_VELOCITY_CHANGE_IN_EVENT_HISTORY,
        DVD_COUNT_EDGE_FLIPPING_IN_COMPUTATION,
        DVD_COUNT_DISK_COLLISION_IN_COMPUTATION,
        DVD_COUNT_DISK_VELOCITY_CHANGE_IN_COMPUTATION,
        DVD_COUNT_2_IDENTICAL_FLIPPING_TIME,
        DVD_COUNT_3_IDENTICAL_FLIPPING_TIME,
        DVD_COUNT_4_IDENTICAL_FLIPPING_TIME,
        DVD_COUNT_5_IDENTICAL_FLIPPING_TIME,
        DVD_COUNT_SHADOW_FLIPPING,
        DVD_COUNT_ORDINARY_FLIPPING,
        DVD_COUNT_3_TO_2_FLIPPING,
        DVD_COUNT_24_TO_33_FLIPPING,
        DVD_COUNT_33_TO_24_FLIPPING,
        DVD_COUNT_33_TO_22_FLIPPING,
        DVD_COUNT_DEGREE0_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE1_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE2_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE3_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE4_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE5_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE6_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE7_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_DEGREE8_POLYNOMIAL_FOR_EDGE_FLIPPING,
        DVD_COUNT_NUM_OF_REAL_ROOT_0,
        DVD_COUNT_NUM_OF_REAL_ROOT_1,
        DVD_COUNT_NUM_OF_REAL_ROOT_2,
        DVD_COUNT_NUM_OF_REAL_ROOT_3,
        DVD_COUNT_NUM_OF_REAL_ROOT_4,
        DVD_COUNT_NUM_OF_REAL_ROOT_5,
        DVD_COUNT_NUM_OF_REAL_ROOT_6,
        DVD_COUNT_NUM_OF_REAL_ROOT_7,
        DVD_COUNT_NUM_OF_REAL_ROOT_8,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_0,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_1,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_2,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_3,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_4,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_5,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_6,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_7,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_8,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_0,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_1,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_2,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_3,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_4,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_5,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_6,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_7,
        DVD_COUNT_NUM_OF_POSITIVE_REAL_ROOT_SATIFYING_INEQUALITY_8,
        DVD_COUNT_ALL_DIFFERENT_VELOCITY_VECTOR,
        DVD_COUNT_2_SAME_VELOCITY_VECTOR,
        DVD_COUNT_3_SAME_VELOCITY_VECTOR,
        DVD_COUNT_4_SAME_VELOCITY_VECTOR,
        DVD_COUNT_2_2_SAME_VELOCITY_VECTOR,
        DVD_COUNT_DEGREE_REDUCING,
        DVD_COUNT_SIZE_PHASE_ONE
    };

    //const PolynomialSolverType POLYNOMIAL_SOLVER = STURMS_SEQ;
    //const PolynomialSolverType POLYNOMIAL_SOLVER = JENKINS_TRAUB;
    const PolynomialSolverType POLYNOMIAL_SOLVER   = PolynomialSolverType::LABSRC;

    /*******************           USED TOLERANCE IN DVD                ******************/
        const double TOLERANCE_OF_EDGE_FLIP_ROOT                       = resNeg4;
        const double TOLERANCE_OF_GAVRILOVA_INEQUALITY                 = resNeg4;
        const double TOLERANCE_OF_SHADOW_FLIP                          = resNeg4;
        const double TOLERANCE_OF_DISK_INTERSECTION                    = resNeg6;
        const double TOLERANCE_OF_EQUAL_RADIUS                         = resNeg4;
        const double TOLERANCE_OF_ZERO_VELOCITY                        = resNeg6;
        const double TOLERANCE_OF_COLLISION_PREDICTION                 = resNeg6;
        const double TOLERANCE_OF_EQUAL_VELOCITY_VEC                   = resNeg6;
        const double TOLERANCE_OF_TOP_ELEMENT_WITH_IDENTICAL_KEY_IN_PQ = resNeg4;
        const int    MAX_ITERATION_OF_FALSE_POSITION                   = 100;
        const double TOLERANCE_OF_FALSE_POSITION                       = resNeg12;
    /*************************************************************************************/

private:
    //For input parameter               
    double                              m_UpperBoundOfTimeHorizon;
    double                              m_CoefficientOfRestitution;
    CollisionType                       m_CollisionType;

    //For output
    vector<EventOfDynamicVD2DC*>        m_EventHistory;
    int                                 m_IndexOfMostRecentEvent;
                                        
    //For simulation
    unordered_map<int, Generator2D*>    m_LinkFromInputIDToGenerator; //including container
    double                              m_CurrentTime;
   
    //******************************** For event generation ********************************//
        //After phase1, the elements of these data structure will be cleared.

        //For calculating events
        EntityAccessiblePriorityQ<VEdge2D*>                        m_PriorityQForEdgeFlipping;
        EntityAccessiblePriorityQ<VEdge2D*>                        m_PriorityQForDiskCollision;
        
        EntityAccessiblePriorityQ_<pair<Generator2D*, rg_Point2D>> m_PriorityQForDiskVelocityChange;
 
        //For special case
        set<VEdge2D*>                       m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ;
        set<VEdge2D*>                       m_VEdgesDefinedByIdentical4DiskWithFlippingVEdge;
    //**************************************************************************************//

    //For statistics                    
    TimeStatistics<milli>               m_TimeStatistics[TIME_STATISTICS_SIZE];
    CountStatistics                     m_CountStatistics[COUNT_STATISTICS_SIZE];
    
    bool                                m_SetStatisticsForPhase1;
    bool                                m_SetStatisticsForPhase2;

#ifdef TEST_FOR_TOPOLOGY_RESTORATION
    //************************** TEST For RESTORING wrong topology *****************************//
        list<pair<VEdge2D*, double>>   m_FourOldEdgeNFlippingTimeIncidentToCurrentFlipEdge;
        double                         m_OldCollisionTimeRelatedToFlipEdge;

        bool restore_VD_topology_if_flipping_time_is_wrong(EventOfDynamicVD2DC* event);
        void store_disk_location_at_current_time(const set<Generator2D*>& diskGenerators, const double& targetTime, list<rg_Triplet<Generator2D*, double, rg_Point2D>>& diskCentersItsUpdatedTimes) const;
        void make_generators_as_in_original_status(const list<rg_Triplet<Generator2D*, double, rg_Point2D>>& diskCentersItsUpdatedTimes);
    //**************************************************************************************//
#endif //TEST_FOR_TOPOLOGY_RESTORATION

public:
	//Constructor
    DynamicVoronoiDiagram2DC();
    DynamicVoronoiDiagram2DC(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks);
    DynamicVoronoiDiagram2DC(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks);
	DynamicVoronoiDiagram2DC(const DynamicVoronoiDiagram2DC& dynamicVD);

	//Deconstructor
	~DynamicVoronoiDiagram2DC();

	//Getter
    inline double get_current_time()                      const { return m_CurrentTime; };
    inline double get_coefficient_of_restitution()        const { return m_CoefficientOfRestitution; };
    inline double get_simulation_terminating_time()       const { return m_UpperBoundOfTimeHorizon; };
    inline bool   get_collision_type() const { return m_CollisionType; };

    inline const unordered_map<int, Generator2D*>& get_linker_from_input_id_to_generator() const { return m_LinkFromInputIDToGenerator; };
    inline void                        get_event_history(vector<EventOfDynamicVD2DC*>& eventHistory) { eventHistory = m_EventHistory; };
    inline const EventOfDynamicVD2DC*  get_most_recent_event()        const { return (m_IndexOfMostRecentEvent == -1) ? NULL : m_EventHistory.at(m_IndexOfMostRecentEvent); };
    inline int                         get_most_recent_event_idx()    const { return m_IndexOfMostRecentEvent; };
    inline int                         get_total_event_num()          const { return m_EventHistory.size(); };
    inline double                      get_event_time(const int& idx) const { return (idx >= 0) ? m_EventHistory.at(idx)->get_occurring_time() : DBL_MAX; };
    void                               get_dynamic_disks(list<DynamicDisk*>& dynamicDisks) const;
    DynamicDisk*                       get_container() const;

    inline const TimeStatistics<milli>&       get_time_statistics_in_phase1()  const { return m_TimeStatistics[PHASE1_COMP_TIME]; };
    inline const CountStatistics&             get_count_statistics_in_phase1() const { return m_CountStatistics[PHASE1_COUNT]; };
    inline const TimeStatistics<milli>&       get_time_statistics_int_phase2() const { return  m_TimeStatistics[PHASE2_COMP_TIME]; };

    //Setter
    inline void set_current_time(const double& currentTime) { m_CurrentTime = currentTime; };
    void        set_simulation_parameter(const double& upperBoundOfTimeHorizon, const double& coefficientOfRestitution, const CollisionType& collisionType);
    void        set_simulation_terminating_time(const double& upperBoundOfTimeHorizon);
    void        set_coefficient_of_restitution(const double& coefficientOfRestitution);
    void        set_collision_type(const CollisionType& collisionType);
    void        set_offset_of_disks(const double& offsetOfDisk);

    //Operator
    DynamicVoronoiDiagram2DC& operator=(const DynamicVoronoiDiagram2DC& dynamicVD);

    //Filltering
    inline bool initial_DVD_is_constructed() const { return ((m_generators.size() > 0) ? true : false); };
    bool        there_is_a_pair_of_disks_intersected_by_VD(const double& offset = 0.0, const double& tolerance = rg_MATH_RES);
    bool        all_disks_are_stopped() const;

    //For link input disk to generators
    VEdge2D*     get_VEdge_defined_by_these_four_input_disks(int* fourDiskIds) const;
    Generator2D* get_generator_by_its_input_disk_id(const int& inputDiskId) const;
    void         link_dynamic_disk_id_to_its_generator_pointer();

   
    //Function: Construct VD by dynamic disks.
    void construct_initial_CIC_VD(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks);

    void construct_initial_CIC_VD(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks);


    //Main function 1: Generate an event history and Go to VD at any target time in the input time horizon. 
      //PHASE ONE
    void generate_event_history_PHASE_ONE(const double& inputUpperBoundOfTimeHorizon, 
                                          const double& coefficientOfRestitution,
                                          const CollisionType& collisionType);

    void generate_event_history_ALLOWING_initial_intersection_btw_disks_and_MOVING_only_initially_dynamic_disks_PHASE_ONE(
                                          const double& inputUpperBoundOfTimeHorizon,
                                          const double& coefficientOfRestitution);

    void generate_event_history_ALLOWING_disk_offset_PHASE_ONE(
                                          const double& inputUpperBoundOfTimeHorizon,
                                          const double& coefficientOfRestitution,
                                          const double& offsetOfDisk);

    //Use event history which is calculated outside of this class.
    void insert_already_computed_events_into_event_history_N_set_DVD_at_time_zero_PHASE_ONE(const vector<EventOfDynamicVD2DC*>& eventHistory);

      //PHASE TWO
    void construct_VD_at_target_time_by_scanning_event_history_PHASE_TWO(const double& inputTargetTime);

    void construct_VD_at_target_time_by_scanning_event_history_phase_two_with_moving_dynamic_disks_N_updating_topology_SEPARATELY_PHASE_TWO(const double& inputTargetTime);

    void construct_VD_with_time_increment_by_scanning_event_history_PHASE_TWO(const double& timeIncrement = 1.0);
   
    void construct_VD_with_event_increment_by_scanning_event_history_PHASE_TWO(const int& eventNumIncrement = 1);


    //Main function 2: Go to VD at a target time with calculating events at the same time.
    void preprocess_dynamic_simulation_of_DVD(const DynamicDisk& container, const list<DynamicDisk>& dynamicDisks, const double& coefficientOfRestitution, const CollisionType& collisionType);

    void preprocess_dynamic_simulation_of_DVD(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks, const double& coefficientOfRestitution, const CollisionType& collisionType);

    void construct_VD_at_target_time_by_computing_events_occuring_until_the_target_time(const double& inputTargetTime);


    //Etc.
    void clear();


    //For statistics
    void write_statistics_in_phase_one(const string& fileNameWithPath) const;
    void write_statistics_in_phase_two(const string& fileNameWithPath) const;


    /********************************************************************************************/
    /*           From here, the functions below can be used by one understanding DVD well.      */
    /*            By these funtions, you can apply many functions consist of DVD to your prob.  */
    /********************************************************************************************/
public:
    //For V-edge flip event
    pair<bool, double> find_possible_flipping_time_of_an_vedge_after(const double& mostRecentlyEventOccurringTime, const VEdge2D* targetVEdge, const bool& isThisVEdgeFlipping = false);
    void collect_influenced_5VEdges_by_flipping_target_vedge_including_itself(const VEdge2D* targetVEdge, set<VEdge2D*>& influencedVEdgeList) const;
    void find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(const double& mostRecentlyEventOccurringTime, const set<VEdge2D*>& influencedEdges);
    void find_possible_flipping_time_of_influenced_VEdges_after_input_time_N_store_in_priorityQ(const VEdge2D* targetVEdge, const double& mostRecentlyEventOccurringTime, const set<VEdge2D*>& influencedEdges);
    void do_bookKeeping_on_priorityQ_for_edge_flipping(const pair<bool, double>& newEdgeFlipTime, const VEdge2D* vEdge, const PriorityQBookkeepingType& PQBookkeepingType);

    inline void flip_target_vedge(VEdge2D* targetVEdge);
    inline void flip_reversely_target_vedge(VEdge2D* targetVEdge);

    //For disk collision event
    pair<bool, double> find_possible_collision_time_of_disks_after(const double& mostRecentlyEventOccurringTime, const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    void collect_influenced_VEdges_by_a_disk_collision_for_updating_VEdge_flipping_times(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& influencedVEdges);
    void collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_PHYSICAL_COLLISION(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& influencedVEdges);
    void collect_influenced_VEdges_by_a_disk_collision_for_updating_collision_times_AVOIDANCE(const Generator2D* diskGenerator1, const Generator2D* diskGenerator2, set<VEdge2D*>& sharedVEdges, set<VEdge2D*>& theOtherVEdges);
    void find_possible_collision_time_of_influeced_VEdges_after_input_time_N_store_in_priorityQ_PHYSICAL_COLLISION(const double& mostRecentlyEventOccurringTime, set<VEdge2D*>& influencedVEdges);
    void find_possible_collision_time_of_shared_VEdges_after_input_time_N_store_in_priorityQ_AVOIDANCE(const double& mostRecentlyEventOccurringTime, set<VEdge2D*>& sharedVEdges);
    void do_bookkeeping_on_priorityQ_for_disk_collision(const pair<bool, double>& newCollisionTime, const VEdge2D* vEdge, const PriorityQBookkeepingType& PQBookkeepingType);


    //For velocity change: this event is from user
    void collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_VEdge_flipping_times(const Generator2D* diskGenerator, set<VEdge2D*>& influencedVEdges);
    void collect_influenced_VEdges_by_a_disk_velocity_change_for_updating_collision_times(const Generator2D* diskGenerator, set<VEdge2D*>& influencedVEdges);
    void insert_disk_velocity_change_event_into_priorityQ(const double& eventTime, DynamicDisk* dynamicDisk, const rg_Point2D& velocityVector);
    bool this_disk_velocity_change_event_is_in_priorityQ(DynamicDisk* dynamicDisk, const rg_Point2D& velocityVector) const;
    double find_the_occuring_time_of_disk_velocity_change_event_in_priorityQ(DynamicDisk* dynamicDisk, const rg_Point2D& velocityVector) const;
    bool disk_velocity_change_priorityQ_is_empty() const;
    void delete_disk_velocity_change_event_from_priorityQ(DynamicDisk* dynamicDisk, const rg_Point2D& velocityVector);


    //For computing velocity vectors 
    void compute_velocity_vectors_of_disks_after_collided_PHYSICAL_COLLISION(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2, rg_Point2D& outputVector1, rg_Point2D& outputVector2) const;
    void compute_velocity_vectors_of_disks_after_collided_for_INTIALLY_MOVED_OR_FIXED_COLLISION(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2, rg_Point2D& outputVector1, rg_Point2D& outputVector2) const;
    void compute_velocity_vectors_of_disks_after_collided_for_AVOIDANCE(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2, rg_Point2D& outputVector1, rg_Point2D& outputVector2) const;
    
    void compute_N_update_with_new_velocity_vectors_and_store_old_N_new_ones(DiskCollisionEvent* diskCollisionEvent);
    inline void store_two_colliding_disk_location_N_velocity_vectors_right_after_colliding(DiskCollisionEvent* diskCollisionEvent, const rg_Point2D& disk1VelocityVectorAfterCollding, const rg_Point2D& disk2VelocityVectorAfterCollding) const;
    inline void update_velocity_vector(DynamicDisk* dynamicDisk, const rg_Point2D& newVelocityVector);
    inline void change_to_new_velocity_vector_and_store_old_one(DiskVelocityChangeEvent* diskVelocityChangeEvent) const;


    //For event queue(priority queue): V-edge flip and disk collision are always tracked by DVD so that we do not pop any events from the Qs.
    EventOfDynamicVD2DC* get_next_event_from_priorityQs();
    EdgeFlipEvent*   generate_next_edge_flip_event_from_priorityQ();
    DiskCollisionEvent*  generate_next_disk_collision_event_from_priorityQ();
    DiskVelocityChangeEvent*  pop_and_generate_next_disk_velocity_change_event_from_priorityQ();
    inline EventOfDynamicVD2DC::EventType get_next_event_type() const;
    inline double        get_next_event_time() const;
    VEdge2D*             get_next_flipping_edge_from_flipping_priorityQ();


    //For moving disks
    void move_disks_to_this_target_time_without_changing_VD(const double& targetTime);
    void move_disks_to_this_target_time_without_changing_VD(set<Generator2D*>& diskGenerators, const double& targetTime);


    //For VD queries
    void update_VD_geometry();
    bool this_VD_is_valid();
    inline void find_4_definining_disks_of_VEdge(set<Generator2D*>& fourDefiningGenerators, VEdge2D* targetVEdge);

    //For VD queries -CHECK!!!!
    
    void compute_all_VVertice_coordinates_of_current_VD();
    void get_VVertices_except_infinite_ones(list<VVertex2D*>& VVerticesExceptInfiniteOnes);
    inline bool this_VEdge_is_boundary_of_anomalizing_cell(const VEdge2D* VEdge) const;
    static bool compare_disk_with_radius_N_XCoor_N_YCoord_in_non_increasing_order(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    static bool compare_disk_with_ID_in_non_decreasing_order(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2);
    inline void compute_circum_circles_Of_VVertices(list<VVertex2D*>& VVertices) { updateCircumcirclesOfNewVVertices(VVertices); };
    inline void compute_all_VEdges_equation_of_current_VD() { updateGeometry(); };
   
    //For generating an event
    void set_NULL_for_recent_update_time_of_dynamic_disks() const;
    void find_influenced_disks_by_next_event(set<Generator2D*>& influencedDiskGenerators, EventOfDynamicVD2DC* anEvent);
    void handle_this_event(EventOfDynamicVD2DC* event);
    void update_the_event_times_in_the_priorityQs_influenced_by_this_event(EventOfDynamicVD2DC* anEvent);
    void store_this_event_in_the_event_history_vector(EventOfDynamicVD2DC* event);


    //For handling special case where upto 5 V-edges can have an identical flip time.
    void clear_info_about_vedges_defined_by_same_disks();


    //For etc.
    void initialize_simulation_variable();
    void prepare_data_set_for_CIC(DynamicDisk* const container, const list<DynamicDisk*>& dynamicDisks, list< pair< rg_Circle2D, void* > >& diskNUserDataPairs);
    void get_dynamic_disks_pointers(const list<DynamicDisk>& dynamicDisks, list<DynamicDisk*>& dynamicDiskPointers);
    void generate_event_history_for_dynamic_disks();
    void generate_all_events_during_simulation_time_window();
    void construct_initial_priority_Qs_for_both_flipping_and_collision_events();
    void construct_initial_priorityQ_of_vedge_flipping_time();
    void construct_initial_priorityQ_of_disk_collision_time();
    void make_VD_N_disks_be_at_time_zero();
    void clear_priorityQs();

    bool this_VEdge_is_in_this_set(const VEdge2D* VEdge, const set<VEdge2D*>& VEdgeSet) const;
    inline bool there_is_intersection_btw(const rg_Circle2D& disk1, const rg_Circle2D& disk2);
    inline void get_4Disks_defining_this_VEdge_and_its_end_points(const VEdge2D* VEdge, vector<DynamicDisk*>& disksForVEdge);
    inline void get_4DiskGenerators_defining_this_VEdge_and_its_end_points(const VEdge2D* VEdge, vector<Generator2D*>& generatorsForVEdge) const;


protected:
	//For generating polynomial of V-edge flip and solving the polynomial
    rg_Polynomial make_polynomial_equation_for_vedge_flipping(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);
	pair<bool, double> solve_polynomial_equation_for_vedge_flipping(const rg_Polynomial& polynomialEq, const double& flippingTimeOfMostRecentlyFlippedEdge, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);
    pair<bool, double> solve_polynomial_equation_for_vedge_flipping_by_LABSRC(const rg_Polynomial& polynomialEq, const double& flippingTimeOfMostRecentlyFlippedEdge, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, const bool& isThisVEdgeFlipping);
    pair<bool, double> solve_polynomial_equation_for_vedge_flipping_by_JENKINS_TRAUB(const rg_Polynomial& polynomialEq, const double& flippingTimeOfMostRecentlyFlippedEdge, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, const bool& isThisVEdgeFlipping);
    pair<bool, double> solve_polynomial_equation_for_vedge_flipping_by_STURMS_SEQUENCE(const rg_Polynomial& polynomialEq, const double& flippingTimeOfMostRecentlyFlippedEdge, const VEdge2D* targetVEdge, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, const bool& isThisVEdgeFlipping);

    void fill_elements_of_all_polynomial_XY_matrix(EdgeFlippingFormulaGeneratingType mode, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, rg_Polynomial(*XYMat)[3]) const;
    void fill_elements_of_all_polynomial_matrices(EdgeFlippingFormulaGeneratingType mode, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, rg_Polynomial(*XMat)[3], rg_Polynomial(*YMat)[3], rg_Polynomial(*XYMat)[3]) const;
    void fill_elements_for_polynomial_matrix(const vector<DynamicDisk*>& movingDiskSet,
                                             rg_Polynomial& x1_x4, rg_Polynomial& x2_x4, rg_Polynomial& x3_x4, 
                                             rg_Polynomial& y1_y4, rg_Polynomial& y2_y4, rg_Polynomial& y3_y4,
                                             rg_Polynomial& r1_r4, rg_Polynomial& r2_r4, rg_Polynomial& r3_r4,
                                             rg_Polynomial& p1,    rg_Polynomial& p2,    rg_Polynomial& p3) const;
    
    void sort_generators_in_non_increasing_order(vector<DynamicDisk*>& disks);
   
	rg_Polynomial make_determinant_polynomial(rg_Polynomial(*matrix)[3]) const;
    void   find_positive_real_roots(rg_ComplexNumber* roots, const int& numOfRoots, list<double>& realRoots) const;
    double find_min_absolute_real_root(rg_ComplexNumber* roots, const int& numOfRoots) const;
    void   erase_numerically_zero_root_from_positive_real_roots(list<double>& realRoots, const double& minAbsoluteRealNumber) const;
    bool   four_disks_are_circumscribing(const double& flippingTime, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder) const;
	bool   all_four_disks_have_same_radius(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder) const;
	bool   this_VEdge_is_NOT_shadow_flipping(VEdge2D* VEdge, const double& eventTime, const rg_REAL resolution = rg_MATH_RES);
    int    compute_polynomial_degree_by_disk_info(const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);
    
    //For STURMS_SEQUENCE solver which is not good 
    double do_false_position_method(const rg_Polynomial& polynomial, 
                                    const double& lowerBound,
                                    const double& upperBound,
                                    const int& maxIteration,
                                    const double& tolerance) const;
    int find_number_of_roots_in_interval_using_sturms_chain(const list<rg_Polynomial>& sturmsChain, const double& intervalBegin, const double& intervalEnd) const;
    list<rg_Polynomial> make_sturms_chain(const rg_Polynomial& polynomial) const;
    pair<rg_Polynomial, rg_Polynomial> divide_polynomial(const rg_Polynomial& dividend, const rg_Polynomial& divider) const;
    int count_sign_transition_in_sturms_chain(const list<rg_Polynomial>& sturmsChain, const double& valueForEvaluation) const;
    double find_lower_bound_for_flipping_VEdge(const rg_Polynomial& polynomial, const list<rg_Polynomial>& sturmsChain, const double& lowerInterval, const double& upperInterval) const;
   

    //For checking whether two disks will be collided. 
    bool two_disks_are_moving_to_opposite_side_relatively(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2) const;
    bool two_disks_will_be_collided(const DynamicDisk* movingDisk1, const DynamicDisk* movingDisk2) const;


    //For selecting a V-edge to be flipped.
    VEdge2D* find_next_flipping_vedge_when_vedges_have_same_flipping_time();
    void find_all_VEdges_defined_by_identical_disks(const VEdge2D* VEdge, const list<VEdge2D*>& candidates, set<VEdge2D*>& outputs);
    void find_all_VEdges_defined_by_identical_disks_TOPOLOGICAL_COMPUTATION(const VEdge2D* VEdge, const set<VEdge2D*>& candidates, set<VEdge2D*>& outputs) const;


    //For playing DVD given event history 
    void compute_VD_at_target_time(const double& targetTime);
    void update_VD_topology_N_move_disks_by_scanning_event_history(const double& targetTime);
    inline SimulationDirectionType find_simulation_direction(const double& targetTime) const;
    inline SimulationDirectionType find_simulation_direction(const int& eventNumIncrement) const;
    void update_VD_N_move_disks_to_future(const double& targetTime);
    void update_VD_N_move_disks_to_past(const double& targetTime);
    void update_VD_N_move_disks_to_future_by_event_increment(const int& eventNumIncrement);
    void update_VD_N_move_disks_to_past_by_event_increment(const int& eventNumIncrement);

    void compute_VD_at_target_time_by_dynamic_disks_N_updating_topology_SEPARATELY(const double& targetTime);
    void update_VD_topology(const double& targetTime, int& updatedIndexOfMostRecentEvent);
    void update_VD_topology_to_future(const double& targetTime, int& updatedIndexOfMostRecentEvent);
    void update_VD_topology_to_past(const double& targetTime, int& updatedIndexOfMostRecentEvent);
    void move_disks(const double& targetTime, int& updatedIndexOfMostRecentEvent);
    void move_disks_to_future(const double& targetTime, int& updatedIndexOfMostRecentEvent);
    void move_disks_to_past(const double& targetTime, int& updatedIndexOfMostRecentEvent);

    inline void move_two_disks_until_colliding_time(const double& collidingTime, Generator2D* diskGenerator1, Generator2D* diskGenerator2);
    inline void move_two_disks_until_colliding_time(const double& collidingTime, Generator2D* diskGenerator1, Generator2D* diskGenerator2, const rg_Point2D& disk1CenterAtColldingTime, const rg_Point2D& disk2CenterAtColldingTime);
    inline void move_a_disk_until_velocity_change_time(const double& velocityChangeTime, Generator2D* diskGenerator);
    inline void move_a_disk_until_velocity_change_time(const double& velocityChangeTime, Generator2D* diskGenerator, const rg_Point2D& diskCenterAtVelocityChangeTime);

    inline void update_two_disks_velocity_vectors_by_collision(Generator2D* diskGenerator1, Generator2D* diskGenerator2, const rg_Point2D& newVelocityVector1, const rg_Point2D& newVelocityVector2);
    inline void update_a_disk_velocity_vector_by_velocity_change(Generator2D* diskGenerator, const rg_Point2D& newVelocityVector);


    //For copying DVD
    void copy_from(const DynamicVoronoiDiagram2DC& dynamicVD);
    void make_link_btw_VEdges(list<VEdge2D*>& fromVEdges, list<VEdge2D*>& toVEdges, map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD);
    void make_link_btw_generators(list<Generator2D*>& fromGenerators, list<Generator2D*>& toGenerators, map<Generator2D*, Generator2D*>& generatorLinkFromOtherVDToThisVD);
    void make_link_btw_disks(list<DynamicDisk*>& fromDisks, list<DynamicDisk*>& toDisks, map<DynamicDisk*, DynamicDisk*>& diskLinkFromOtherVDToThisVD);
    void modify_user_data_moving_disk_of_generators_in_this_VD_from_other_VDs_one_to_this_VDs_one(list<Generator2D*>& allGeneratorsFromThisVD, map<DynamicDisk*, DynamicDisk*>& diskLinkFromOtherVDToThisVD);
    void modify_generators_and_VEdges_in_event_history_from_other_VDs_one_to_this_VDs_one(const map<Generator2D*, Generator2D*>& generatorLinkFromOtherVDToThisVD, const map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD);
    void modify_VEdges_in_priorityQs_from_other_VDs_one_to_this_VDs_one(const map<VEdge2D*, VEdge2D*>& VEdgeLinkFromOtherVDToThisVD);


    //For VD queries - special case detection
    bool these_two_VEdge_are_defined_by_two_same_disks(VEdge2D* VEdge1, VEdge2D* VEdge2) const;
    bool this_VEdge_and_its_hands_and_legs_are_in_this_set(VEdge2D* VEdge, set<VEdge2D*>& VEdgeSet) const;
    bool these_three_VEdge_have_common_disk_in_right_or_left_VFace(const set<VEdge2D*>& VEdges) const;
    void find_VEdges_defined_by_identical_4disks_from_PQ_N_flipping_VEdge(const VEdge2D* flippingVEdge,
                                                                          const set<VEdge2D*> upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ,
                                                                          set<VEdge2D*>& VEdgesDefinedByIdentical4Disk) const;

    //For finding velocity vector in event history 
    void  compute_motion_vectors_by_scanning_event_history();


    //For statistics
    void initialize_statistics_in_phase_one();
    void initialize_statistics_in_phase_two(const double& targetTime);
    void finalize_statistics_in_phase_one(const bool& phase1WillBeContinued = false);
    void finalize_statistics_in_phase_two();

    inline void ___start_clock_in_phase1(const unsigned int& it){ if (m_SetStatisticsForPhase1) { m_TimeStatistics[PHASE1_COMP_TIME].start_clock(it); } };
    inline void ___end_clock_in_phase1(const unsigned int& it)  { if (m_SetStatisticsForPhase1) { m_TimeStatistics[PHASE1_COMP_TIME].end_clock(it); } };
    inline void ___count_in_phase1(const unsigned int& it)      { if (m_SetStatisticsForPhase1) { m_CountStatistics[PHASE1_COUNT].addCount(it); } };

    inline void ___start_clock_in_phase2(const unsigned int& typeOfStatistics) { if (m_SetStatisticsForPhase2) { m_TimeStatistics[typeOfStatistics].start_clock(); } };
    inline void ___end_clock_in_phase2(const unsigned int& typeOfStatistics)   { if (m_SetStatisticsForPhase2) { m_TimeStatistics[typeOfStatistics].end_clock(); } };
    inline void ___flipping_count_in_phase2() { if (m_SetStatisticsForPhase2) { m_CountStatistics[PHASE2_FLIPPING_COUNT].addCount(); } };
    inline void ___collision_count_in_phase2(){ if (m_SetStatisticsForPhase2){ m_CountStatistics[PHASE2_COLLISION_COUNT].addCount(); } };
 
    inline void        ___check_polynomial_degree_for_statistics(const rg_Polynomial& polynomial);
    inline void        ___check_number_of_VEdges_defined_by_identical_4disks_for_statistics();
    inline void        ___check_flipping_type_of_VEdge_for_statistics(EdgeFlipEvent* edgeFlippingEvent);
    void        ___check_roots_of_polynomial_for_statistics(const int& numOfRoots, rg_ComplexNumber* roots, const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder);


#ifdef DEBUG_DVD
    const string debugFileNameWithPath = "DVD-debug.txt";

    void ___initialize_debug() const;
    void ___write_current_event_info(EventOfDynamicVD2DC* anEvent) const;
    void ___write_update_event_info(const EventOfDynamicVD2DC::EventType& eventType, VEdge2D* VEdge, const double newEventTime) const;
#endif //DEBUG_DVD

#ifdef WRITING_EVENT_TIME
    const string eventTimeFileNameWithPath         = "event_times.txt";
    static double accumulatedTimeInEventGeneration = 0.0;

    void ___write_event_times(const double& eventTime, const double& computationTime);
    void ___initialize_writing_event_time();

#endif //WRITING_EVENT_TIME

#ifdef WRITING_POLYNOMIAL_RELATED_INFORMATION
    const int LIMIT_NUM_OF_POLYNOMIALS_ON_MEMORY       = 2000;
    static int   numOfWriting                          = 0;
    const string polynomialCoefficientFileNameWithPath = "polynomial_coefficients.txt";
    const string movingDiskFileNameWithPath            = "dynamic_disks.txt";
    const string realRootFileNameWithPath              = "real_roots.txt";
    const string realRootEvaluationFileNameWithPath    = "real_roots_evalution.txt";
    const string VEdgeFlippingStatus                   = "flipping_status.txt";

    //For getting disk and polynomial information during simulation by writing file 
    //Need about 200 byte per 1 polynomial, so write a file when 5,000 polynomial (1 MB) accumulated. 
    vector<double>                      m_SimulationTimesAtMakingPolynomial;
    vector<vector<DynamicDisk>>         m_QuadrupletsOfMovingDisks;
    vector<rg_Polynomial>               m_Polynomials;
    vector<char>                        m_FlippingStatusForVEdge;
    vector<vector<double>>              m_RealRoots;
    
    //for collecting polynomial info during phase1
    void ___store_this_polynomial_related_info(const double& currTime, 
                                               const vector<DynamicDisk*>& fourDiskGeneratorsInNonIncreasingOrder, 
                                               const rg_Polynomial& polynomial,
                                               const rg_ComplexNumber* roots,
                                               const char& isThisFlippingVEdge);

    void ___clear_all_polynomial_related_info();
    void ___write_dynamic_disks_info();
    void ___write_polynomial_coefficient_info();
    void ___write_real_roots_info();
    void ___write_roots_evaluation_info();
    void ___write_VEdge_flipping_status();
    char ___change_bool_type_of_flipping_status_to_char(const bool& type);
#endif //WRITING_POLYNOMIAL_RELATED_INFORMATION


    /******************************************* To be removed *************************************************/
    void collect_upto_5VEdges_of_same_flipping_time_by_identical_4disks_from_PQ(set<VEdge2D*>& influencedVEdges) const;

    //For debugging
    bool four_generator_ids_are_same_as(list<Generator2D*>& generators, const int& id1, const int& id2, const int& id3, const int& id4);
    bool this_VEdge_is_defined_by_these_two_disks(const VEdge2D* VEdge, const rg_Circle2D& circle1, const rg_Circle2D& circle2);
    /***********************************************************************************************************/

};


inline void DynamicVoronoiDiagram2DC::flip_target_vedge(VEdge2D* targetVEdge)
{
    targetVEdge->flip();
}


inline void DynamicVoronoiDiagram2DC::flip_reversely_target_vedge(VEdge2D* targetVEdge)
{
    targetVEdge->flip_reverse();
}



inline bool DynamicVoronoiDiagram2DC::there_is_intersection_btw(const rg_Circle2D& disk1, const rg_Circle2D& disk2)
{
    double distanceOfCircleCenters = disk1.getCenterPt().distance(disk2.getCenterPt());
    double sumOfRadius = disk1.getRadius() + disk2.getRadius();

    if (distanceOfCircleCenters < sumOfRadius - TOLERANCE_OF_DISK_INTERSECTION)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline void DynamicVoronoiDiagram2DC::change_to_new_velocity_vector_and_store_old_one(DiskVelocityChangeEvent* diskVelocityChangeEvent) const
{
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskVelocityChangeEvent->get_disk_generator()->getUserData());
    
    diskVelocityChangeEvent->set_disk_center_at_velocity_change_time(movingDisk->getCircle().getCenterPt());
    diskVelocityChangeEvent->set_velocity_vector_before_change(movingDisk->getVelocityVector());
    diskVelocityChangeEvent->set_velocity_vector_setting(true);

    movingDisk->setVelocityVector(diskVelocityChangeEvent->get_future_velocity_vector_of_disk());
}



inline void DynamicVoronoiDiagram2DC::store_two_colliding_disk_location_N_velocity_vectors_right_after_colliding(DiskCollisionEvent* diskCollisionEvent,
                                                                                                               const rg_Point2D& disk1VelocityVectorAfterCollding,
                                                                                                               const rg_Point2D& disk2VelocityVectorAfterCollding) const
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator1()->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskCollisionEvent->get_disk_generator2()->getUserData());

    diskCollisionEvent->set_motion_vectors_before_collision(movingDisk1->getVelocityVector(), movingDisk2->getVelocityVector());
    diskCollisionEvent->set_motion_vectors_after_collision(disk1VelocityVectorAfterCollding, disk2VelocityVectorAfterCollding);
    diskCollisionEvent->set_disk_centers_at_collision_time(movingDisk1->getCircle().getCenterPt(), movingDisk2->getCircle().getCenterPt());
    diskCollisionEvent->set_motion_vector_setting(true);
}


inline void DynamicVoronoiDiagram2DC::update_velocity_vector(DynamicDisk* dynamicDisk, const rg_Point2D& newVelocityVector)
{
    dynamicDisk->setVelocityVector(newVelocityVector);
}


inline DynamicVoronoiDiagram2DC::SimulationDirectionType DynamicVoronoiDiagram2DC::find_simulation_direction(const double& targetTime) const
{
    SimulationDirectionType simulationFlowDirection;

    if (targetTime > get_current_time())
    {
        simulationFlowDirection = TO_FUTURE;
    }
    else if (targetTime < get_current_time())
    {
        simulationFlowDirection = TO_PAST;
    }
    else
    { 
        simulationFlowDirection = NONE_DIRECTION;
    }

    return simulationFlowDirection;
}


inline DynamicVoronoiDiagram2DC::SimulationDirectionType DynamicVoronoiDiagram2DC::find_simulation_direction(const int& eventNumIncrement) const
{
    SimulationDirectionType simulationFlowDirection;

    if (eventNumIncrement > 0)
    {
        simulationFlowDirection = TO_FUTURE;
    }
    else if (eventNumIncrement < 0)
    {
        simulationFlowDirection = TO_PAST;
    }
    else
    {
        simulationFlowDirection = NONE_DIRECTION;
    }

    return simulationFlowDirection;
}



inline void DynamicVoronoiDiagram2DC::move_two_disks_until_colliding_time(const double& collidingTime, Generator2D* diskGenerator1, Generator2D* diskGenerator2)
{
    set<Generator2D*> collidedDiskGenerators;

    if (diskGenerator1 != NULL)
    {
        collidedDiskGenerators.insert(diskGenerator1);
    }

    if (diskGenerator2 != NULL)
    {
        collidedDiskGenerators.insert(diskGenerator2);
    }

    move_disks_to_this_target_time_without_changing_VD(collidedDiskGenerators, collidingTime);
}


inline void DynamicVoronoiDiagram2DC::move_a_disk_until_velocity_change_time(const double& velocityChangeTime, Generator2D* diskGenerator)
{
    set<Generator2D*> collidedDiskGenerators;

    if (diskGenerator != NULL)
    {
        collidedDiskGenerators.insert(diskGenerator);
    }


    move_disks_to_this_target_time_without_changing_VD(collidedDiskGenerators, velocityChangeTime);
}


inline void DynamicVoronoiDiagram2DC::move_a_disk_until_velocity_change_time(const double& velocityChangeTime, Generator2D* diskGenerator, const rg_Point2D& diskCenterAtVelocityChangeTime)
{
    DynamicDisk* dynamicDisk = static_cast<DynamicDisk*>(diskGenerator->getUserData());

    dynamicDisk->setMostRecentLocationUpdateTime(velocityChangeTime);

    rg_Circle2D diskAtVelocityChangeTime(diskCenterAtVelocityChangeTime, dynamicDisk->getCircle().getRadius());

    dynamicDisk->setCircle(diskAtVelocityChangeTime);
    diskGenerator->setDisk(diskAtVelocityChangeTime);
}



inline void DynamicVoronoiDiagram2DC::move_two_disks_until_colliding_time(const double& collidingTime,
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


inline void DynamicVoronoiDiagram2DC::update_two_disks_velocity_vectors_by_collision(Generator2D* diskGenerator1, Generator2D* diskGenerator2, const rg_Point2D& newVelocityVector1, const rg_Point2D& newVelocityVector2)
{
    DynamicDisk* movingDisk1 = static_cast<DynamicDisk*>(diskGenerator1->getUserData());
    DynamicDisk* movingDisk2 = static_cast<DynamicDisk*>(diskGenerator2->getUserData());

    movingDisk1->setVelocityVector(newVelocityVector1);
    movingDisk2->setVelocityVector(newVelocityVector2);
}


inline void DynamicVoronoiDiagram2DC::update_a_disk_velocity_vector_by_velocity_change(Generator2D* diskGenerator, const rg_Point2D& newVelocityVector)
{
    DynamicDisk* movingDisk = static_cast<DynamicDisk*>(diskGenerator->getUserData());
    movingDisk->setVelocityVector(newVelocityVector);
}


inline void DynamicVoronoiDiagram2DC::get_4Disks_defining_this_VEdge_and_its_end_points(const VEdge2D* VEdge, vector<DynamicDisk*>& disksForVEdge)
{
    ___start_clock_in_phase1(DVD_TIME_GET_4DISKS_DEFINING_VEDGE);

    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getRightFace()->getGenerator()->getUserData()));
    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getLeftFace()->getGenerator()->getUserData()));
    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getMateFace(VEdge->getStartVertex())->getGenerator()->getUserData()));
    disksForVEdge.push_back(static_cast<DynamicDisk*>(VEdge->getMateFace(VEdge->getEndVertex())->getGenerator()->getUserData()));

    ___end_clock_in_phase1(DVD_TIME_GET_4DISKS_DEFINING_VEDGE);

}


inline EventOfDynamicVD2DC::EventType DynamicVoronoiDiagram2DC::get_next_event_type() const
{
    EventOfDynamicVD2DC::EventType minKeyType = EventOfDynamicVD2DC::EDGE_FLIP;
    double    minKey     = DBL_MAX;

    if (!m_PriorityQForEdgeFlipping.empty())
    {
        if (m_PriorityQForEdgeFlipping.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2DC::EDGE_FLIP;
            minKey = m_PriorityQForEdgeFlipping.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskCollision.empty())
    {
        if (m_PriorityQForDiskCollision.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2DC::DISK_COLLISION;
            minKey = m_PriorityQForDiskCollision.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskVelocityChange.empty())
    {
        if (m_PriorityQForDiskVelocityChange.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2DC::DISK_VELOCITY_CHANGE;
            minKey = m_PriorityQForDiskVelocityChange.topNode()->getKey();
        }
    }


    return minKeyType;
}


inline double DynamicVoronoiDiagram2DC::get_next_event_time() const
{
    EventOfDynamicVD2DC::EventType minKeyType = EventOfDynamicVD2DC::EDGE_FLIP;
    double    minKey     = DBL_MAX;

    if (!m_PriorityQForEdgeFlipping.empty())
    {
        if (m_PriorityQForEdgeFlipping.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2DC::EDGE_FLIP;
            minKey = m_PriorityQForEdgeFlipping.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskCollision.empty())
    {
        if (m_PriorityQForDiskCollision.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2DC::DISK_COLLISION;
            minKey = m_PriorityQForDiskCollision.topNode()->getKey();
        }
    }

    if (!m_PriorityQForDiskVelocityChange.empty())
    {
        if (m_PriorityQForDiskVelocityChange.topNode()->getKey() <= minKey)
        {
            minKeyType = EventOfDynamicVD2DC::DISK_VELOCITY_CHANGE;
            minKey = m_PriorityQForDiskVelocityChange.topNode()->getKey();
        }
    }

    return minKey;
}


inline void DynamicVoronoiDiagram2DC::get_4DiskGenerators_defining_this_VEdge_and_its_end_points(const VEdge2D * VEdge, vector<Generator2D*>& generatorsForVEdge ) const
{
    generatorsForVEdge.push_back(VEdge->getRightFace()->getGenerator());
    generatorsForVEdge.push_back(VEdge->getLeftFace()->getGenerator());
    generatorsForVEdge.push_back(VEdge->getMateFace(VEdge->getStartVertex())->getGenerator());
    generatorsForVEdge.push_back(VEdge->getMateFace(VEdge->getEndVertex())->getGenerator());
}



inline void DynamicVoronoiDiagram2DC::find_4_definining_disks_of_VEdge(set<Generator2D*>& fourDefiningGenerators, VEdge2D* targetVEdge)
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



inline bool DynamicVoronoiDiagram2DC::this_VEdge_is_boundary_of_anomalizing_cell(const VEdge2D* VEdge) const
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


inline void DynamicVoronoiDiagram2DC::___check_number_of_VEdges_defined_by_identical_4disks_for_statistics()
{
    if (!m_SetStatisticsForPhase1)
    {
        return;
    }

    switch (m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.size())
    {
    case 2:
        ___count_in_phase1(DVD_COUNT_2_IDENTICAL_FLIPPING_TIME);
        break;

    case 3:
        ___count_in_phase1(DVD_COUNT_3_IDENTICAL_FLIPPING_TIME);
        break;

    case 4:
        ___count_in_phase1(DVD_COUNT_4_IDENTICAL_FLIPPING_TIME);
        break;

    case 5:
        ___count_in_phase1(DVD_COUNT_5_IDENTICAL_FLIPPING_TIME);
        break;
    }
}



inline void DynamicVoronoiDiagram2DC::___check_flipping_type_of_VEdge_for_statistics(EdgeFlipEvent* edgeFlippingEvent)
{
    if (edgeFlippingEvent->get_whether_this_vedge_to_be_flipped())
    {
        switch (m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.size())
        {
        case 1:
        {            
            ___count_in_phase1(DVD_COUNT_ORDINARY_FLIPPING);
            break;
        }

        case 2:
        {      
            VEdge2D* VEdge1 = *m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin();
            VEdge2D* VEdge2 = *(++m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ.begin());

            if (these_two_VEdge_are_defined_by_two_same_disks(VEdge1, VEdge2))
            {
                ___count_in_phase1(DVD_COUNT_24_TO_33_FLIPPING);
            }

            break;
        }

        case 3:
        {
            if (these_three_VEdge_have_common_disk_in_right_or_left_VFace(m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ))
            {
                ___count_in_phase1(DVD_COUNT_3_TO_2_FLIPPING);
            }

            break;
        }

        case 5:
        {
            VEdge2D* flippingVEdge = edgeFlippingEvent->get_VEdge();

            if (this_VEdge_and_its_hands_and_legs_are_in_this_set(flippingVEdge, m_Upto5VEdgesOfSameFlippingTimeByIdentical4DisksInPQ))
            {
                ___count_in_phase1(DVD_COUNT_33_TO_22_FLIPPING);
            }
            else
            {
                ___count_in_phase1(DVD_COUNT_33_TO_24_FLIPPING);
            }

            break;
        }
        }
    }
    else
    {
        ___count_in_phase1(DVD_COUNT_SHADOW_FLIPPING);
    }
}



inline void DynamicVoronoiDiagram2DC::___check_polynomial_degree_for_statistics(const rg_Polynomial& polynomial)
{
    if (!m_SetStatisticsForPhase1)
    {
        return;
    }

    switch (polynomial.getDegree())
    {
    case 0:
        ___count_in_phase1(DVD_COUNT_DEGREE0_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 1:
        ___count_in_phase1(DVD_COUNT_DEGREE1_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 2:
        ___count_in_phase1(DVD_COUNT_DEGREE2_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 3:
        ___count_in_phase1(DVD_COUNT_DEGREE3_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 4:
        ___count_in_phase1(DVD_COUNT_DEGREE4_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 5:
        ___count_in_phase1(DVD_COUNT_DEGREE5_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 6:
        ___count_in_phase1(DVD_COUNT_DEGREE6_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 7:
        ___count_in_phase1(DVD_COUNT_DEGREE7_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    case 8:
        ___count_in_phase1(DVD_COUNT_DEGREE8_POLYNOMIAL_FOR_EDGE_FLIPPING);
        break;

    default:
        break;
    }
}

} // GeometryTier
} // V

#endif











