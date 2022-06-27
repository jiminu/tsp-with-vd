#ifndef _OPERATIONSFORPOLYGONMESHMODEL_H
#define _OPERATIONSFORPOLYGONMESHMODEL_H

#include "PolygonMeshModel.h"

class OperationsForPolygonMeshModel
{
public:
    struct VertexWithProperty {
        PolygonMeshModel::Vertex*   m_vertex;
        double                      m_property;
    };

    struct ShellWithVolume {
        PolygonMeshModel::Shell*    m_shell;
        double                      m_volume;
    };


public:
    OperationsForPolygonMeshModel(void);
    ~OperationsForPolygonMeshModel(void);

    static bool LessThanVertexWithProperty( const VertexWithProperty& v1, const VertexWithProperty& v2); 
    static bool GreaterThanShellWithVolume( const ShellWithVolume& s1, const ShellWithVolume& s2); 
};



inline bool OperationsForPolygonMeshModel::LessThanVertexWithProperty( const VertexWithProperty& v1, const VertexWithProperty& v2) {
    return (v1.m_property < v2.m_property) ? true : false;
}

inline bool OperationsForPolygonMeshModel::GreaterThanShellWithVolume( const ShellWithVolume& s1, const ShellWithVolume& s2) {
    return (s1.m_volume > s2.m_volume ) ? true : false;
}

#endif


