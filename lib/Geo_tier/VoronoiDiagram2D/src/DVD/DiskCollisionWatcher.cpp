#include "DiskCollisionWatcher.h"
using namespace V::GeometryTier;


#include <stdexcept>
#include <iostream>
using namespace std;

const DiskCollisionWatcher::WatchObject DiskCollisionWatcher::WATCH_OBJ_FOR_ALL_PAIR_WISE_DISKS 
                    = DiskCollisionWatcher::WatchObject(DiskCollisionWatcher::ALL_COLLISION, DiskCollisionWatcher::WatchIDPair(-1,-1));

DiskCollisionWatcher::DiskCollisionWatcher()
{
    m_NumOfWatchObjectsSatisfyingTargetCollisionNumber = 0;
    m_AllCollisionsBtwDisksAreBeingWatched             = false;
}


DiskCollisionWatcher::DiskCollisionWatcher(const DiskCollisionWatcher& diskCollisionWatcher)
{
    copy_from(diskCollisionWatcher);
}


DiskCollisionWatcher::DiskCollisionWatcher(DiskCollisionWatcher&& diskCollisionWatcher)
{
    move_from(move(diskCollisionWatcher));
}


DiskCollisionWatcher::~DiskCollisionWatcher()
{
}


DiskCollisionWatcher& DiskCollisionWatcher::operator=(const DiskCollisionWatcher& diskCollisionWatcher)
{
    if (this != &diskCollisionWatcher)
    {
        copy_from(diskCollisionWatcher);
    }

    return *this;
}


DiskCollisionWatcher& DiskCollisionWatcher::operator=(DiskCollisionWatcher&& diskCollisionWatcher)
{
    if (this != &diskCollisionWatcher)
    {
        move_from(move(diskCollisionWatcher));
    }

    return *this;
}


void DiskCollisionWatcher::copy_from(const DiskCollisionWatcher& diskCollisionWatcher)
{
    m_MapFromWatchObjectTo_TARGET_CollisionNumber     
        = diskCollisionWatcher.m_MapFromWatchObjectTo_TARGET_CollisionNumber;

    m_MapFromWatchObjectTo_CURRENT_CollisionNumber   
        = diskCollisionWatcher.m_MapFromWatchObjectTo_CURRENT_CollisionNumber;

    m_MapFromWatchObjectTo_CollisionNumberIsSatisfied 
        = diskCollisionWatcher.m_MapFromWatchObjectTo_CollisionNumberIsSatisfied;

    m_NumOfWatchObjectsSatisfyingTargetCollisionNumber
        = diskCollisionWatcher.m_NumOfWatchObjectsSatisfyingTargetCollisionNumber;

    m_AllCollisionsBtwDisksAreBeingWatched
        = diskCollisionWatcher.m_AllCollisionsBtwDisksAreBeingWatched;
}


void DiskCollisionWatcher::move_from(DiskCollisionWatcher&& diskCollisionWatcher)
{
    m_MapFromWatchObjectTo_TARGET_CollisionNumber
        = move(diskCollisionWatcher.m_MapFromWatchObjectTo_TARGET_CollisionNumber);

    m_MapFromWatchObjectTo_CURRENT_CollisionNumber
        = move(diskCollisionWatcher.m_MapFromWatchObjectTo_CURRENT_CollisionNumber);

    m_MapFromWatchObjectTo_CollisionNumberIsSatisfied
        = move(diskCollisionWatcher.m_MapFromWatchObjectTo_CollisionNumberIsSatisfied);

    m_NumOfWatchObjectsSatisfyingTargetCollisionNumber
        = diskCollisionWatcher.m_NumOfWatchObjectsSatisfyingTargetCollisionNumber;

    m_AllCollisionsBtwDisksAreBeingWatched
        = diskCollisionWatcher.m_AllCollisionsBtwDisksAreBeingWatched;
}


void DiskCollisionWatcher::insert_watch_object_and_target_collision_number(const WatchObject& watchObj, const COLLISION_NUM& targetNumber)
{
    WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);

    if (orderedWatchObject.first == ALL_COLLISION)
    {
        m_AllCollisionsBtwDisksAreBeingWatched = true;
    }

    if (   m_MapFromWatchObjectTo_TARGET_CollisionNumber.find(orderedWatchObject) 
        != m_MapFromWatchObjectTo_TARGET_CollisionNumber.end() )
    {
        m_MapFromWatchObjectTo_TARGET_CollisionNumber.at(orderedWatchObject) = targetNumber;
    }
    else
    {
        m_MapFromWatchObjectTo_TARGET_CollisionNumber.insert(make_pair(orderedWatchObject, targetNumber));
        prepare_storage_for_tracking_this_watch_object(orderedWatchObject);
    }
}


void DiskCollisionWatcher::insert_watch_object_and_target_collision_number(
    const WatchType& watchType, const int& id1, const int& id2, const COLLISION_NUM& targetNumber)
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    insert_watch_object_and_target_collision_number(watchObject, targetNumber);
}


void DiskCollisionWatcher::increase_current_collision_number(const WatchObject& watchObj)
{
    WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);

    try
    {
        m_MapFromWatchObjectTo_CURRENT_CollisionNumber.at(orderedWatchObject) += 1;

        if (current_collision_number_is_equal_to_target_collision_number(orderedWatchObject))
        {
            ++m_NumOfWatchObjectsSatisfyingTargetCollisionNumber;

            m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.at(orderedWatchObject) = true;
        }
    }
    catch (const out_of_range& oor) {
        cerr << "In increase_current_collision_number(watchObj): " << oor.what() << endl;
    }
}


void DiskCollisionWatcher::increase_current_collision_number(
    const WatchType& watchType, const int& id1, const int& id2)
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    increase_current_collision_number(watchObject);
}


void DiskCollisionWatcher::decrease_current_collision_number(const WatchObject& watchObj)
{
    WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);

    try
    {
        if (current_collision_number_is_equal_to_target_collision_number(orderedWatchObject))
        {
            --m_NumOfWatchObjectsSatisfyingTargetCollisionNumber;

            m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.at(orderedWatchObject) = false;
        }

        m_MapFromWatchObjectTo_CURRENT_CollisionNumber.at(orderedWatchObject) -= 1;
    }
    catch (const out_of_range& oor) {
        cerr << "In increase_current_collision_number(watchObj): " << oor.what() << endl;
    }
}

void DiskCollisionWatcher::decrease_current_collision_number(const WatchType& watchType, const int& id1, const int& id2)
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    decrease_current_collision_number(watchObject);
}


bool DiskCollisionWatcher::this_is_watch_object(const WatchObject& watchObj) const
{
    WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);

    if (   m_MapFromWatchObjectTo_TARGET_CollisionNumber.find(orderedWatchObject) 
        != m_MapFromWatchObjectTo_TARGET_CollisionNumber.end() )
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool DiskCollisionWatcher::this_is_watch_object(
    const WatchType& watchType, const int& id1, const int& id2) const
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    return this_is_watch_object(watchObject);
}


bool DiskCollisionWatcher::target_collision_number_is_satisfied(const WatchObject& watchObj) const
{
    bool collisionNumberIsSatisfied = false;

    try
    {
        WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);
        collisionNumberIsSatisfied     = m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.at(orderedWatchObject);
    }
    catch (const out_of_range& oor) {
        cerr << "In target_collision_number_is_satisfied(watchObj): " << oor.what() << endl;
    }

    return collisionNumberIsSatisfied;
}


bool DiskCollisionWatcher::target_collision_number_is_satisfied(
    const WatchType& watchType, const int& id1, const int& id2) const
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    return target_collision_number_is_satisfied(watchObject);
}


bool DiskCollisionWatcher::all_of_target_collision_numbers_are_satisfied() const
{
    if (m_NumOfWatchObjectsSatisfyingTargetCollisionNumber == m_MapFromWatchObjectTo_TARGET_CollisionNumber.size())
    {
        return true;
    }
    else
    {
        return false;
    }

    /*
        bool allOfTargetCollisionNumbersAreSatisfied = true;

        for (UMapFromWatchObjToColNumIsSatisfied::const_iterator it_Satisfied = m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.begin();
             it_Satisfied != m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.end();
             ++it_Satisfied)
        {
            if ((*it_Satisfied).second == false)
            {
                allOfTargetCollisionNumbersAreSatisfied = false;
                break;
            }
        }

        return allOfTargetCollisionNumbersAreSatisfied;
    */
}


bool DiskCollisionWatcher::all_collision_btw_disks_are_being_watched() const
{
    return m_AllCollisionsBtwDisksAreBeingWatched;
}


int DiskCollisionWatcher::find_target_collision_number(const WatchObject& watchObj) const
{
    int numOfTargetCollisionNumber = -1;
    try
    {
        WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);
        numOfTargetCollisionNumber     = m_MapFromWatchObjectTo_TARGET_CollisionNumber.at(orderedWatchObject);
    }
    catch (const out_of_range& oor) {
        cerr << "In find_target_collision_number(watchObj): " << oor.what() << endl;
    }

    return numOfTargetCollisionNumber;
}


int DiskCollisionWatcher::find_target_collision_number(
    const WatchType& watchType, const int& id1, const int& id2) const
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    return find_target_collision_number(watchObject);
}


int DiskCollisionWatcher::find_current_collision_number(const WatchObject& watchObj) const
{
    int numOfCurrentCollisionNumber = -1;
    try
    {
        WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);
        numOfCurrentCollisionNumber    = m_MapFromWatchObjectTo_CURRENT_CollisionNumber.at(orderedWatchObject);
    }
    catch (const out_of_range& oor) {
        cerr << "In find_current_collision_number(watchObj): " << oor.what() << endl;
    }

    return numOfCurrentCollisionNumber;
}


int DiskCollisionWatcher::find_current_collision_number(
    const WatchType& watchType, const int& id1, const int& id2) const
{
    WatchObject watchObject = make_watch_object(watchType, id1, id2);
    return find_current_collision_number(watchObject);
}


void DiskCollisionWatcher::clear_current_collision_number()
{
    m_NumOfWatchObjectsSatisfyingTargetCollisionNumber = 0;

    for (UMapFromWatchObjToColNumIsSatisfied::iterator it_Satisfied = m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.begin();
        it_Satisfied != m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.end();
        ++it_Satisfied)
    {
        (*it_Satisfied).second = false;
    }

    for (UMapFromWatchObjToColNum::iterator it_CurrentNum = m_MapFromWatchObjectTo_CURRENT_CollisionNumber.begin();
        it_CurrentNum != m_MapFromWatchObjectTo_CURRENT_CollisionNumber.end();
        ++it_CurrentNum)
    {
        (*it_CurrentNum).second = 0;
    }
}

bool DiskCollisionWatcher::current_collision_number_is_equal_to_target_collision_number(const WatchObject& watchObj) const
{
    unsigned int currentCollisionNumber = 0;
    unsigned int targetCollisionNumber  = 0;

    try
    {
        WatchObject orderedWatchObject = make_ordered_watch_object(watchObj);
        currentCollisionNumber = m_MapFromWatchObjectTo_CURRENT_CollisionNumber.at(orderedWatchObject);
        targetCollisionNumber  = m_MapFromWatchObjectTo_TARGET_CollisionNumber.at(orderedWatchObject);
    }
    catch (const out_of_range& oor) {
        cerr << "In current_collision_number_is_less_than_target_collision_number(watchObj): " << oor.what() << endl;
    }

    if (currentCollisionNumber == targetCollisionNumber)
    {
        return true;
    }
    else
    {
        return false;
    }
}


DiskCollisionWatcher::WatchObject DiskCollisionWatcher::make_ordered_watch_object(const WatchObject& watchObject) const
{
    WatchObject orderedWatchObject = watchObject;

    if (watchObject.second.first > watchObject.second.second)
    {
        orderedWatchObject.second.first  = watchObject.second.second;
        orderedWatchObject.second.second = watchObject.second.first;
    }

    return orderedWatchObject;
}


DiskCollisionWatcher::WatchObject DiskCollisionWatcher::make_watch_object(
    const WatchType& watchType, const int& id1, const int& id2) const
{
    DiskCollisionWatcher::WatchObject watchObject;
    watchObject.first         = watchType;
    watchObject.second.first  = id1;
    watchObject.second.second = id2;
    return watchObject;
}


void DiskCollisionWatcher::prepare_storage_for_tracking_this_watch_object(const WatchObject& orderedWatchObject)
{
    prepare_storage_for_tracking_current_collision_number(orderedWatchObject);
    prepare_storage_for_tracking_whether_target_collision_number_is_satisfied(orderedWatchObject);
}


void DiskCollisionWatcher::prepare_storage_for_tracking_current_collision_number(const WatchObject& orderedWatchObject)
{
    m_MapFromWatchObjectTo_CURRENT_CollisionNumber.insert(make_pair(orderedWatchObject, 0));
}


void DiskCollisionWatcher::prepare_storage_for_tracking_whether_target_collision_number_is_satisfied(const WatchObject& orderedWatchObject)
{
    m_MapFromWatchObjectTo_CollisionNumberIsSatisfied.insert(make_pair(orderedWatchObject, false));
}

