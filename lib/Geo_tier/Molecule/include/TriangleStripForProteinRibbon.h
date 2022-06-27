#ifndef _TRIANGLESTRIPFORPROTEINRIBBON_H
#define _TRIANGLESTRIPFORPROTEINRIBBON_H

#include "rg_Const.h"
#include "rg_Point3D.h"


namespace V {
namespace GeometryTier {



class Chain;

class TriangleStripForProteinRibbon
{
private:
    rg_INT      m_numVertices;    
    rg_Point3D* m_vertex[2];
    rg_Point3D* m_normal[2];


public:
    TriangleStripForProteinRibbon();
    TriangleStripForProteinRibbon(const TriangleStripForProteinRibbon& tStrip);
    ~TriangleStripForProteinRibbon();

    rg_INT       getNumVertices() const;
    rg_Point3D** getVertices();
    rg_Point3D** getNormals();

    void clean();
    void makeProteinRibbonByCarsonBugg86(Chain* currChain);


    TriangleStripForProteinRibbon& operator =(const TriangleStripForProteinRibbon& tStrip);
};



} // namespace GeometryTier
} // namespace V


#endif


