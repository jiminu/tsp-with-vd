#ifndef DISKGENERATOR2D_H
#define DISKGENERATOR2D_H

#include "rg_Point2D.h"
#include "rg_Circle2D.h"
#include "Generator2D.h"

#include <list>
using namespace std;

namespace V {
namespace GeometryTier {

class VFace2D;

class DiskGenerator2D : public Generator2D
{
public:
	DiskGenerator2D();
	DiskGenerator2D(const int& ID);
	DiskGenerator2D(const int& ID, const rg_Circle2D& disk);
	DiskGenerator2D(const int& ID, const rg_Circle2D& disk, VFace2D* const face);
	DiskGenerator2D(const int& ID, const rg_Circle2D& disk, void* const userData);
	DiskGenerator2D(const DiskGenerator2D& generator);
	~DiskGenerator2D();

	DiskGenerator2D&    operator=(const DiskGenerator2D& generator);
	bool                operator==(const DiskGenerator2D& generator) const;

	void                getNeighborGeneratorsInThisVoronoiDiagram(list<DiskGenerator2D*>& neighborGeneratorsList, VoronoiDiagram2DC* const VD);
	void                getNeighborGenerators(list<DiskGenerator2D*>& neighborGeneratorsList, const bool& isThisContainer) const;

};


} // GeometryTier
} // V

#endif


