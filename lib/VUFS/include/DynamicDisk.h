#ifndef _DYNAMIC_DISK_H
#define _DYNAMIC_DISK_H


#include "Disk.h"

class DynamicDisk : public Disk
{

protected:
    rg_Point2D  m_VelocityVector;

    double      m_Offset;
    double      m_MostRecentLocationUpdateTime;
    bool        m_ThisIsContainer;
    int         m_GroupId;

public:
    DynamicDisk();
    DynamicDisk(const int& ID, 
                const rg_Circle2D& circle, 
                const rg_Point2D& velocityVector,
                const double& mostRecentLocationUpdateTime);

    DynamicDisk(const int& ID, 
                const rg_Circle2D& circle, 
                const rg_Point2D& velocityVector, 
                const double& mostRecentLocationUpdateTime,
                const double& offset);

    DynamicDisk(const DynamicDisk& dynamicDisk);
    
    virtual ~DynamicDisk();
    virtual void copy_from(const DynamicDisk& dynamicDisk);

    inline double             getOffset()          const { return m_Offset; };
    inline double             getVelocityVectorX() const { return m_VelocityVector.getX(); };
    inline double             getVelocityVectorY() const { return m_VelocityVector.getY(); };
    inline const rg_Point2D&  getVelocityVector()  const { return m_VelocityVector; };
    inline const rg_Circle2D& getCircle()          const { return *this; };
    inline double             getMostRecentLocationUpdateTime() const { return m_MostRecentLocationUpdateTime; };
    inline int                getGroupId()                       const { return m_GroupId; };

    inline void               setOffset(const double& offset) { m_Offset = offset; };
    inline void               setVelocityVector(const double& vecX, const double& vecY) { m_VelocityVector.setPoint(vecX, vecY); };
    inline void               setVelocityVector(const rg_Point2D& velocityVector)       { m_VelocityVector = velocityVector; };
    inline void               setMostRecentLocationUpdateTime(const double& locationUpdateTime){ m_MostRecentLocationUpdateTime = locationUpdateTime; };
    inline void               setGroupId(const int& diskGroupId) { m_GroupId = diskGroupId; };
    inline void               set_this_disk_is_container(const bool& thisDiskIsContainer) { m_ThisIsContainer = thisDiskIsContainer; };

    inline bool               this_disk_is_container() const { return m_ThisIsContainer; };


    DynamicDisk& operator =(const DynamicDisk& dynamicDisk);
    bool        operator<(const DynamicDisk& dynamicDisk) const;
    bool        operator>(const DynamicDisk& dynamicDisk) const;
};


#endif

