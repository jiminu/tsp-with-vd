#ifndef _FUNCTIONS_FOR_QUASITRIANGULATION_H_
#define _FUNCTIONS_FOR_QUASITRIANGULATION_H_

#include "ComponentBSFace.h"

namespace V {
namespace GeometryTier {

class FunctionsForQuasiTriangulation2D
{
public:
	//static rg_INT compute_Euler_characteristic(rg_dList<ComponentBSFace>& components);
	static rg_INT number_of_incident_regular_edges_of_vertex_on_Beta_shape_boundary(VertexBU2D* vertex, const rg_REAL& betaValue);
	static rg_INT incident_regular_edges_of_vertex_on_Beta_shape_boundary(VertexBU2D* vertex, const rg_REAL& betaValue, list<EdgeBU2D*>& incidentRegularEdges);
};

}
}

#endif


