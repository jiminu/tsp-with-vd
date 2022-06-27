#ifndef _DYNAMICVD2D_SIMULATOR_
#define _DYNAMICVD2D_SIMULATOR_

#define DEBUG_DVD_SIM

#include <cstring>
#include <iomanip>
#include "DefaultValuesForDVD2.h"
#include "SimpleTypesForDVD2.h"
#include "DynamicVoronoiDiagramCIC.h"

namespace V {
namespace GeometryTier {


class DynamicVD2DSimulator
{
protected:
    list<DynamicDisk*>               m_AllContainers;
    list<DynamicDisk*>               m_DynamicDisks;
    DynamicVoronoiDiagramCIC	     m_DynamicVD;

    string                           m_DiskFileNameWthPath;
    TimeStatistics<milli>            m_CompTimeStatistics;

    bool                             m_DynamicDisksAreGeneratedInsideOfThisObject;

public:
    DynamicVD2DSimulator();
    DynamicVD2DSimulator(const DynamicVD2DSimulator& DVDSimulator);

    virtual ~DynamicVD2DSimulator();

    DynamicVD2DSimulator& operator=(const DynamicVD2DSimulator& dvdSimulator);

     void                        get_dynamic_disks(list<DynamicDisk*>& dynamicDisks) const;
     void                        get_containers   (list<DynamicDisk*>& containers)   const;
     const DynamicDisk*          get_outer_most_container()                          const;
     const DynamicVoronoiDiagramCIC& get_dynamic_VD()                                const;
     string                      get_disk_file_name_with_path()                      const;
     double                      get_length_of_Voronoi_clock()                       const;
     double                      get_coefficient_of_restitution()                    const;
     double                      get_current_Voronoi_clock()                         const;
     int                         get_num_of_containers()                             const;
     int                         get_num_of_disks()                                  const;


	 const EventOfDynamicVD2D*   get_most_recent_event_of_DVD()     const;
	 int                         get_total_num_of_events_of_DVD()   const;

     const DiskCollisionWatcher& get_disk_collision_watcher()       const;

     const TimeStatistics<milli>& get_DVDSIM_time_statistics()      const;


     void set_current_Voronoi_clock(const double& currentClock);
     void set_length_of_Voronoi_clock(const double& lengthOfVoronoiClock);
     void set_disk_file_name(const string& diskFileName) { m_DiskFileNameWthPath = diskFileName; };


    //////////////////////////////////////////////////////////////////////////
    // PART 1. File Load/Save
    //////////////////////////////////////////////////////////////////////////
    
    void load_disks(const string& fileNameWithFullPath);
    void load_disks(const list<DynamicDisk*>& containers, const list<DynamicDisk*>& dynamicDisks);
    void load_disks(const list<DynamicDisk>& containers, const list<DynamicDisk>& dynamicDisks);

    virtual void load_EVENTSEQ(const string& fileNameWithFullPath, const bool& reinstateInitialVD = true);
    virtual void load_EVENTSEQ(const string& fileNameWithFullPath, const DiskCollisionWatcher& diskCollisionWatcher, const bool& reinstateInitialVD = true);
    void load_disk_hopping_probability_matrix(const string& fileNameWithFullPath);
    
    virtual void generate_EVENTSEQ(const double& lengthOfVoronoiClock, const bool& reinstateInitialVD = true);
    virtual void generate_EVENTSEQ(const int& numOfCollisions, const bool& reinstateInitialVD = true);
    virtual void generate_EVENTSEQ(const DiskCollisionWatcher& collisionWatcher,
                                   const bool& reinstateInitialVD = true);

    DiskCollisionEvent* generate_EVENTSEQ_by_next_collision(const double& endVoronoiClock = DBL_MAX, 
                                                            const bool& reinstateInitialVD = true);

    void rewind_EVENTSEQ(const double& targetVoronoiClock);

    ReservationNumber insert_disk_velocity_change_event_into_priorityQ(const double& eventTime, const int& inputDiskId, const rg_Point2D& velocityVector);
    void remove_disk_velocity_change_event_from_priorityQ(const ReservationNumber& reservationNum);



    //test code...
    virtual void generate_EVENTSEQ(const double& lengthOfVoronoiClock, const DiskCollisionHandler& collisionHandlingType, const bool& reinstateInitialVD = true);
    virtual void generate_EVENTSEQ_collision_avoidance(const double& lengthOfVoronoiClock, const double& timeRelativeToCollisionTime, const bool& reinstateInitialVD = true);


    void write_EVENTSEQ();
    void write_EVENTSEQ(const string& fileNameWithFullPath);

    void write_average_distances_for_one_disks_to_first_neighbors(const string& fileNameWithPath) const;
    void write_speed_of_all_disks(const string& fileNameWithPath) const;

    void write_disk_contact_numbers_by_time_increment(const string& fileNameWithPath) const;

    void write_statistics_for_generating_EVENTSEQ()                                const;
	void write_statistics_for_generating_EVENTSEQ(const string& fileNameWithPath)    const;
    void write_data_loading_time_statistics(const string& fileNameWithPath) const;
  


    //////////////////////////////////////////////////////////////////////////
    // PART 2. Simulation Control
    //////////////////////////////////////////////////////////////////////////
    void go_to_target_clock_either_to_FUTURE_or_to_PAST(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed = false);
    void increase_Voronoi_clock_either_to_FUTURE_or_to_PAST(const double& VoronoiClockIncrement = 1.0, const bool& VEdgeGeometryIsComputed = false);
    void go_to_next_event_Voronoi_clock_either_to_FUTURE_or_to_PAST(const int& eventNumIncrement = 1, const bool& VEdgeGeometryIsComputed = false);
    void reinstate_the_initial_Voronoi_diagram();
    void compute_VEdges_geometry_of_VD();


    //////////////////////////////////////////////////////////////////////////
    // PART 3. VD queries
    //////////////////////////////////////////////////////////////////////////
    list<int> find_neighbor_disk_ids(const int& anchorDiskId, const double& startVoronoiClock, const double& endVoronoiClock, const bool& thisIsContainer /*= false*/);


    //////////////////////////////////////////////////////////////////////////
    // PART 4. Collision related APIs
    //////////////////////////////////////////////////////////////////////////
    list<pair<int, double>> find_future_collisions_if_target_disk_moves_differently(
        const int& targetDiskId, const rg_Point2D& velocityVector, const list<int>& otherDiskIds);


    // DO NOT WORK !!! NEED NEW DVD
    list<pair<int, double>> find_collisions_with_neighbors_if_target_disk_moves_differently_at_current_clock(
        const int& targetDiskId, const rg_Point2D& velocityVector, const double& endVoronoiClock);



    //////////////////////////////////////////////////////////////////////////
    // PART 5. Statistics
    //////////////////////////////////////////////////////////////////////////
    double      compute_density()                            const;
    double      compute_current_average_velocity_magnitude() const;
    double      compute_current_energy()                     const;
    rg_Point2D  compute_current_momentum()                   const;

	 void ___start_to_track_statistics_in_playing_simulation(const double& targetVoronoiClock);
	 void ___stop_to_track_statistics_in_playing_simulation();

     void ___start_clock_in_DVD_simulator(const unsigned int& it) {  m_CompTimeStatistics.start_clock(it); };
     void ___end_clock_in_DVD_simulator(const unsigned int& it)   {  m_CompTimeStatistics.end_clock(it);  };


    //////////////////////////////////////////////////////////////////////////
    // PART 6. Data set validation test
    //////////////////////////////////////////////////////////////////////////
    bool       data_set_is_valid();
    bool       there_is_a_pair_of_disks_intersected_by_VD(const double& offset = 0.0, const double& tolerance = rg_MATH_RES);
    bool       all_disks_are_stopped() const;



protected:
    //For virtual functions
	virtual void progress_DVD_to_FUTURE(const double& targetVoronoiClock);
	virtual void progress_DVD_to_PAST(const double& targetVoronoiClock);
	virtual void progress_DVD_using_EVENTNUM_to_FUTURE(const int& eventNumIncrement);
	virtual void progress_DVD_using_EVENTNUM_to_PAST(const int& eventNumIncrement);

    virtual  void handle_EDGE_FLIP_event_to_FUTURE(EdgeFlipEvent* edgeFlipEvent);
	virtual  void handle_DISK_COLLISION_event_to_FUTURE(DiskCollisionEvent* diskCollisionEvent);
	virtual  void handle_DISK_VELOCITY_CHANGE_event_to_FUTURE(DiskVelocityChangeEvent* diskVelocityChangeEvent);
	virtual  void handle_DISK_HOPPING_event_to_FUTURE(DiskHoppingEvent* diskHoppingEvent);
   
    virtual  void handle_EDGE_FLIP_event_to_PAST(EdgeFlipEvent* edgeFlipEvent);
    virtual  void handle_DISK_COLLISION_event_to_PAST(DiskCollisionEvent* diskCollisionEvent);
    virtual  void handle_DISK_VELOCITY_CHANGE_event_to_PAST(DiskVelocityChangeEvent* diskVelocityChangeEvent);
    virtual  void handle_DISK_HOPPING_event_to_PAST(DiskHoppingEvent* diskHoppingEvent);

    void handle_next_DVD_event_to_FUTURE();
    void reversely_handle_next_DVD_event_to_PAST();
    void synchronize_disk_locations_and_current_clock_to_VD(const double& targetVoronoiClock);

     EventOfDynamicVD2D* get_next_DVD_event_to_FUTURE()    const;
     EventOfDynamicVD2D* get_next_DVD_event_to_PAST()      const;
     bool         next_DVD_event_to_FUTURE_exists() const;
     bool         next_DVD_event_to_PAST_exists()   const;
     bool         next_DVD_event_to_FUTURE_occurs_before_or_at_this_Voronoi_clock(const double& targetVoronoiClock) const;
     bool         next_DVD_event_to_PAST_occurs_before_this_Voronoi_clock(const double& targetVoronoiClock)   const;
     bool         next_DVD_event_to_PAST_occurs_at_ZERO_Voronoi_clock() const;

    //For creating and deleting dynamic disks
     DynamicDisk* create_disk(const DynamicDisk& disk) const;
     void delete_disk(DynamicDisk*& disk) const;

    void                create_disks(const list<DynamicDisk>& input, list<DynamicDisk*>& output);
    void                create_dynamic_disks(const list<DynamicDisk>& dynamicDisks);
    void                create_containers(const list<DynamicDisk>& containers);
    void                delete_disks(list<DynamicDisk*>& disks);
    void                delete_dynamic_disks();
    void                delete_containers();

    void                initialize_dynamic_disks_for_simulation();         
    void                set_containers(const list<DynamicDisk*>& containers);
    void                set_dynamic_disks(const list<DynamicDisk*>& disks);

    void                remove_disk(DynamicDisk* disk);

    //For handling priority queue
    void construct_initial_priority_queues();
    void clear_priority_queues_for_collision_and_flip();
    tuple<bool, double, DynamicDisk*, DynamicDisk*> propagate_DVD_until_next_collision();


    //For time progress
    void   go_back_to_target_Voronoi_clock_and_remove_events_after_the_target_Voronoi_clock(const double& targetVoronoiClock);

    
     bool this_clock_is_future(const double& targetVoronoiClock) const;
     bool this_clock_is_past(const double& targetVoronoiClock)   const;

    virtual void parse_input_file_and_create_dynamic_objects(const string& fileNameWithFullPath);

private:
    void construct_initial_CIC_VD();

	void construct_VD_at_target_Voronoi_clock(const double& targetVoronoiClock, const bool& VEdgeGeometryIsComputed = false);
    void compute_all_VVertex_coordinates_of_VD();
     void increment_index_of_most_recent_event();
	 void decrement_index_of_most_recent_event();
 
	//For reading and writing history file
	void write_edge_flipping_event(EventOfDynamicVD2D* anEvent, ofstream& fileStream)  const;
	void write_disk_collision_event(EventOfDynamicVD2D* anEvent, ofstream& fileStream) const;
	void write_disk_velocity_change_event(EventOfDynamicVD2D* anEvent, ofstream& fileStream) const;
    void write_disk_hopping_event(EventOfDynamicVD2D* anEvent, ofstream& fileStream)   const;

     void write_hexa_characters(ofstream& fout, const int& size, unsigned char* chars) const;

	 void         get_hexa_characters_of_DOUBLE_type_by_order_of_little_endian(const double& value, unsigned char* chars) const;
	 double       get_double_type_value_from_hexa_characters_in_order_of_little_endian(char* charsOfValue) const;
	 unsigned int char_to_val(char aCharacter) const;

	void extract_DVD_data_from_EVENTSEQ_file(
        const string& fileNameWithFullPath, 
        vector<EventOfDynamicVD2D*>& eventSequence, 
        double& lengthOfVoronoiClock, 
        double& eventSeqFileVers, 
        DiskCollisionHandler& collisionHandlingType) const;

	void set_input_parameters_to_DVD(
        const double& lengthOfVoronoiClock, 
        const DiskCollisionHandler& collisionHandlingType);

    void set_input_parameters_to_DVD(
        const double& lengthOfVoronoiClock,
        const DiskCollisionHandler& collisionHandlingType,
        const DiskCollisionWatcher& collisionWatcher);

    void set_EVENTSEQ_to_DVD(
        const vector<EventOfDynamicVD2D*>& EVENTSEQ, const double eventSeqFileVers, const bool& reinstateInitialVD= true);


    void parse_comment_parts_of_EVENTSEQ_file(
        ifstream& fin, double& lengthOfVoronoiClock, double& eventSeqFileVers, DiskCollisionHandler& collisionHandlingType) const;
  
        void                parse_event_definitions_of_of_EVENTSEQ_file(ifstream& fin, vector<EventOfDynamicVD2D*>& eventSequence) const;
        EventOfDynamicVD2D* parse_event_type_and_create_event(const char* token) const;
        double               parse_to_double(const char* token) const;
        void                 parse_and_set_edge_flip_event(EdgeFlipEvent* edgeFlipEvent, const vector<string>& tokens) const;
        void                 parse_and_set_disk_collision_event(DiskCollisionEvent* diskCollisionEvent, const vector<string>& tokens) const;
        void                 parse_and_set_disk_velocity_change_event(DiskVelocityChangeEvent* diskVelocityChangeEvent, const vector<string>& tokens) const;
		void                 parse_and_set_disk_hopping_event(DiskHoppingEvent* diskHoppingEvent, const vector<string>& tokens) const;

    string create_EVENTSEQ_file_name() const;
    string create_EVNETSEQ_generation_statistics_file_name() const;

    list<DynamicDisk*>          create_dynamic_objects_from_current_dynamic_objects() const;
    list<DynamicDisk*>          create_containers_from_current_containers()           const;
    list<DynamicDisk*>          create_from(const list<DynamicDisk*>& input)          const;
    virtual DynamicDisk*        create_dynamic_object(DynamicDisk* dynamicObject)     const;

    void change_user_data_of_generators_in_DVD_to_point_current_dynamic_objects_or_containers();
    void change_user_data_of_generators_in_DVD_to_point_current_one(const list<DynamicDisk*>& input);
    void make_map_from_new_disks_to_copied_disks_at_time_ZERO_in_DVD();

    static bool compare_collision_clock(const pair<int,double>& collidedDiskNClock1, const pair<int,double>& collidedDiskNClock2);
    	 
protected:
    void clear();
    void clear_EVENTSEQ_related_objects();
    void copy_from(const DynamicVD2DSimulator& dvdSimulator);
};


}
}

#endif //_DYNAMICVD2DC_SIMULATOR_
