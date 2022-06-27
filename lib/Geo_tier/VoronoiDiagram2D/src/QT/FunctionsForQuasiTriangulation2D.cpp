#include "FunctionsForQuasiTriangulation2D.h"
using namespace V::GeometryTier;
//rg_INT FunctionsForQuasiTriangulation2D::compute_Euler_characteristic(rg_dList<ComponentBSFace>& components)
//{
//    rg_INT V = 0;
//    rg_INT E = 0;
//    rg_INT F = 0;
//    rg_INT EulerCharacteristic = V + E - F;
//    return EulerCharacteristic;
//}

rg_INT FunctionsForQuasiTriangulation2D::number_of_incident_regular_edges_of_vertex_on_Beta_shape_boundary(VertexBU2D* vertex, const rg_REAL& betaValue)
{
    rg_dList<EdgeBU2D*> incidentEdges;
    rg_INT numEdges = vertex->getIncidentEdges(incidentEdges);

    rg_INT degreeOfVertex = 0;
    incidentEdges.reset4Loop();
    while (incidentEdges.setNext4Loop())
    {
        EdgeBU2D* currEdge = incidentEdges.getEntity();
        rg_INT boundingState = currEdge->getBoundingState(betaValue);
        if (boundingState == REGULAR_SIMPLEX)
            degreeOfVertex++;
    }

    return degreeOfVertex;
}


rg_INT FunctionsForQuasiTriangulation2D::incident_regular_edges_of_vertex_on_Beta_shape_boundary(VertexBU2D* vertex, const rg_REAL& betaValue, list<EdgeBU2D*>& incidentRegularEdges)
{
    rg_dList<EdgeBU2D*> incidentEdges;
    rg_INT numEdges = vertex->getIncidentEdges(incidentEdges);

    rg_INT degreeOfVertex = 0;
    incidentEdges.reset4Loop();
    while (incidentEdges.setNext4Loop())
    {
        EdgeBU2D* currEdge = incidentEdges.getEntity();
        rg_INT boundingState = currEdge->getBoundingState(betaValue);
        if (boundingState == REGULAR_SIMPLEX)
        {
            incidentRegularEdges.push_back(currEdge);
            degreeOfVertex++;
        }
    }

    return degreeOfVertex;
}