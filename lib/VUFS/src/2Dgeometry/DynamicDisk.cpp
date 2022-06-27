#include "DynamicDisk.h"

DynamicDisk::DynamicDisk()
    :Disk(-1)
{
    m_Offset                       = 0.0;
    m_MostRecentLocationUpdateTime = 0.0;  
    m_VelocityVector               = rg_Point2D(0.0, 0.0);
    m_ThisIsContainer              = false;
    m_GroupId                      = 0;
}


DynamicDisk::DynamicDisk(const int& ID, const rg_Circle2D& circle, const rg_Point2D& velocityVector, const double& mostRecentLocationUpdateTime)
    :Disk(ID, circle), m_VelocityVector(velocityVector)
{
    m_Offset                       = 0.0;
    m_MostRecentLocationUpdateTime = mostRecentLocationUpdateTime;
    m_ThisIsContainer              = false;
    m_GroupId                      = 0;
}


DynamicDisk::DynamicDisk(const int& ID, const rg_Circle2D& circle, const rg_Point2D& velocityVector, const double& mostRecentLocationUpdateTime, const double& offset)
    :Disk(ID, circle), m_VelocityVector(velocityVector)
{
    m_Offset                       = offset;
    m_MostRecentLocationUpdateTime = mostRecentLocationUpdateTime;
    m_ThisIsContainer              = false;
    m_GroupId                  = 0;
}


DynamicDisk::DynamicDisk(const DynamicDisk& dynamicDisk)
{
    copy_from(dynamicDisk);
}



DynamicDisk::~DynamicDisk()
{
}


void DynamicDisk::copy_from(const DynamicDisk& dynamicDisk)
{
     setID(dynamicDisk.getID());
     setCircle(dynamicDisk.getCircle());

     m_Offset                       = dynamicDisk.m_Offset;
     m_VelocityVector               = dynamicDisk.m_VelocityVector;
     m_MostRecentLocationUpdateTime = dynamicDisk.m_MostRecentLocationUpdateTime;
     m_ThisIsContainer              = dynamicDisk.m_ThisIsContainer;
     m_GroupId                  = dynamicDisk.m_GroupId;
}


DynamicDisk& DynamicDisk::operator=(const DynamicDisk& dynamicDisk)
{
    if (this != &dynamicDisk) 
    {
        copy_from(dynamicDisk);
    }

    return *this;
}


bool DynamicDisk::operator<(const DynamicDisk& dynamicDisk) const
{
    if (getRadius() < dynamicDisk.getRadius())
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool DynamicDisk::operator>(const DynamicDisk& dynamicDisk) const
{
    if (getRadius() > dynamicDisk.getRadius())
    {
        return true;
    }
    else
    {
        return false;
    }
}

#ifdef PYVORONOI
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
namespace py = pybind11;
void init_DynamicDisk(py::module& m) {
    py::class_<DynamicDisk, std::shared_ptr<DynamicDisk>>(m, "DynamicDisk")
        .def(py::init([]() { return new DynamicDisk(); }))
        .def(py::init([](const DynamicDisk& generator) { return new DynamicDisk(generator); }))
        //please check
        .def(py::init([](const int& ID, const rg_Circle2D& circle, const rg_Point2D& velocityVector, const double& mostRecentLocationUpdateTime) { return new DynamicDisk(ID, circle, velocityVector, mostRecentLocationUpdateTime); }))
        .def(py::init([](const int& ID, const rg_Circle2D& circle, const rg_Point2D& velocityVector, const double& mostRecentLocationUpdateTime, const double& offset) { return new DynamicDisk(ID, circle, velocityVector, mostRecentLocationUpdateTime, offset); }))
        
        .def("copy_from", &DynamicDisk::copy_from)

        .def("get_offset", &DynamicDisk::getOffset)
        .def("get_velocity_vectorX", &DynamicDisk::getVelocityVectorX)
        .def("get_velocity_vectorY", &DynamicDisk::getVelocityVectorY)
        .def("get_velocity_vector", &DynamicDisk::getVelocityVector)
        .def("get_circle", &DynamicDisk::getCircle)
        .def("get_most_recent_location_update_time", &DynamicDisk::getMostRecentLocationUpdateTime)
        .def("get_group_id", &DynamicDisk::getGroupId)
        
        .def("set_offset", &DynamicDisk::setOffset)
        //please check
        .def("set_velocity_vector", static_cast< void (DynamicDisk::*)(const double& vecX, const double& vecY) > (&DynamicDisk::setVelocityVector))
        .def("set_velocity_vector", static_cast< void (DynamicDisk::*)(const rg_Point2D& velocityVector) > (&DynamicDisk::setVelocityVector))

        .def("set_most_recent_location_update_time", &DynamicDisk::setMostRecentLocationUpdateTime)
        .def("set_group_id", &DynamicDisk::setGroupId)
        .def("set_this_disk_is_container", &DynamicDisk::set_this_disk_is_container)

        .def("this_disk_is_container", &DynamicDisk::this_disk_is_container)

        .def(py::self == py::self)
        .def(py::self > py::self)
        .def(py::self < py::self);
}
#endif