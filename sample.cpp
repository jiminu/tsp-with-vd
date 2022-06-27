#include <iostream>
#include <rg_Circle2D.h>
#include <VoronoiDiagramCIC.h>
#include <list>

using namespace V::GeometryTier;

int main(int, char**) {
    rg_Circle2D circleOne(1,1,0);
    rg_Circle2D circleTwe(1,0,0);
    rg_Circle2D circleThree(2,1,0);
    std::list<rg_Circle2D> circles;
    circles.push_back(circleOne);
    circles.push_back(circleTwe);
    circles.push_back(circleThree);
    
    VoronoiDiagramCIC VD;
    
    VD.constructVoronoiDiagramCIC_noContainerInInput(circles);
    
    list<VEdge2D*> edges;
    VD.getVoronoiEdges(edges);
    std::cout<< edges.size() <<endl;
    std::cout << "Hello, world!\n";
}
