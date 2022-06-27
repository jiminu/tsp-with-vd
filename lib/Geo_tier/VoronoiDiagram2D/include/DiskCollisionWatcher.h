#ifndef _DISK_COLLISION_WATCHER_
#define _DISK_COLLISION_WATCHER_

#include "SimpleTypesForDVD2.h"
#include <unordered_map>
#include <tuple>
#include <utility>
using namespace std;


namespace V {
namespace GeometryTier {


class DiskCollisionWatcher
{
public:
    enum WatchType   {DISK_ID_PAIR, GROUP_ID_PAIR, ALL_COLLISION};

    typedef pair<int, int>                     WatchIDPair;
    typedef pair<WatchType, WatchIDPair>       WatchObject;
    typedef unsigned int                       COLLISION_NUM;
    typedef bool                               TARGET_COLLISION_NUM_IS_SATISFIED;

    static const WatchObject WATCH_OBJ_FOR_ALL_PAIR_WISE_DISKS;

    struct HashForWatchObject 
    {
        std::size_t operator () (const WatchObject& watchObject) const 
        {
            size_t seed = 0;
            hash_combiner(seed, watchObject.first);
            hash_combiner(seed, watchObject.second.first);
            hash_combiner(seed, watchObject.second.second);
            return seed;
        }
    };

    struct EqualityPredicateForWatchObject 
    {
        bool operator () (const WatchObject& p, const WatchObject& q) const
        {
            if (p.first != q.first)
            {
                return false;
            }
            else
            {
                const WatchIDPair& idPairFor_p = p.second;
                const WatchIDPair& idPairFor_q = q.second;

                if (idPairFor_p.first == idPairFor_q.first && idPairFor_p.second == idPairFor_q.second)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
    };

    typedef unordered_map<WatchObject, COLLISION_NUM, HashForWatchObject, EqualityPredicateForWatchObject> UMapFromWatchObjToColNum;
    typedef unordered_map<WatchObject, TARGET_COLLISION_NUM_IS_SATISFIED, HashForWatchObject, EqualityPredicateForWatchObject> UMapFromWatchObjToColNumIsSatisfied;

private:
    UMapFromWatchObjToColNum             m_MapFromWatchObjectTo_TARGET_CollisionNumber;
    UMapFromWatchObjToColNum             m_MapFromWatchObjectTo_CURRENT_CollisionNumber;
    UMapFromWatchObjToColNumIsSatisfied  m_MapFromWatchObjectTo_CollisionNumberIsSatisfied;
    int                                  m_NumOfWatchObjectsSatisfyingTargetCollisionNumber;
    bool                                 m_AllCollisionsBtwDisksAreBeingWatched;

public:
    DiskCollisionWatcher();
    DiskCollisionWatcher(const DiskCollisionWatcher& diskCollisionWatcher);
    DiskCollisionWatcher(DiskCollisionWatcher&& diskCollisionWatcher);
    virtual ~DiskCollisionWatcher();

    DiskCollisionWatcher& operator =(const DiskCollisionWatcher& diskCollisionWatcher);
    DiskCollisionWatcher& operator =(DiskCollisionWatcher&& diskCollisionWatcher);

    void insert_watch_object_and_target_collision_number(const WatchObject& watchObj, const COLLISION_NUM& targetNumber);
    void insert_watch_object_and_target_collision_number(const WatchType& watchType, const int& id1, const int& id2, const COLLISION_NUM& targetNumber);

    void increase_current_collision_number(const WatchObject& watchObj);
    void increase_current_collision_number(const WatchType& watchType, const int& id1, const int& id2);

    void decrease_current_collision_number(const WatchObject& watchObj);
    void decrease_current_collision_number(const WatchType& watchType, const int& id1, const int& id2);

    bool this_is_watch_object(const WatchObject& watchObj)                                const; 
    bool this_is_watch_object(const WatchType& watchType, const int& id1, const int& id2) const;

    bool target_collision_number_is_satisfied(const WatchObject& watchObj)                                const;
    bool target_collision_number_is_satisfied(const WatchType& watchType, const int& id1, const int& id2) const;

    bool all_of_target_collision_numbers_are_satisfied()                   const;
    bool all_collision_btw_disks_are_being_watched()                       const;

    int  find_target_collision_number(const WatchObject& watchObj)                                 const; 
    int  find_target_collision_number(const WatchType& watchType, const int& id1, const int& id2)  const;

    int  find_current_collision_number(const WatchObject& watchObj)                                const;  
    int  find_current_collision_number(const WatchType& watchType, const int& id1, const int& id2) const;

    void clear_current_collision_number();

private:
    void prepare_storage_for_tracking_this_watch_object(const WatchObject& orderedWatchObject);

    void prepare_storage_for_tracking_current_collision_number(const WatchObject& orderedWatchObject);
    void prepare_storage_for_tracking_whether_target_collision_number_is_satisfied(const WatchObject& orderedWatchObject);

    bool current_collision_number_is_equal_to_target_collision_number(const WatchObject& watchObj) const;

    WatchObject make_ordered_watch_object(const WatchObject& watchObject) const;
    WatchObject make_watch_object(const WatchType& watchType, const int& id1, const int& id2) const;

    void copy_from(const DiskCollisionWatcher& diskCollisionWatcher);
    void move_from(DiskCollisionWatcher&& diskCollisionWatcher);
};

}
}


#endif // !_DISK_COLLISION_WATCHER_


