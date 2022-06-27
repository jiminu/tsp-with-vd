#ifndef _EDGE_FLIPPING_EVENT_H_
#define _EDGE_FLIPPING_EVENT_H_

#include <array>
using namespace std;

#include "EventOfDynamicVD2D.h"
#include "VEdge2D.h"
#include "rg_Point2D.h"
#include "Generator2D.h"

namespace V {
namespace GeometryTier {


class EdgeFlipEvent : public EventOfDynamicVD2D
{
public:
    typedef array<Generator2D*, 4> GeneratorPtrQuadruplet;    //index 0 : left, 1 : right, 2 : head, 3: tail

private:
    GeneratorPtrQuadruplet   m_GeneratorQuadrupletBeforeFlip; 
	bool					 m_DoesThisVEdgeToBeFlipped;

public:
	//constructor
	EdgeFlipEvent();
    EdgeFlipEvent(const double& occuringTime);
    EdgeFlipEvent(const double& occuringTime, const GeneratorPtrQuadruplet& generatorQuadruplet);
    EdgeFlipEvent(const double& occuringTime, const GeneratorPtrQuadruplet& generatorQuadruplet, const bool& toBeFlipped);
	EdgeFlipEvent(const EdgeFlipEvent& edgeFlippingEvent);

	//deconstructor
	~EdgeFlipEvent();

    //operator
    EdgeFlipEvent& operator=(const EdgeFlipEvent& edgeFlippingEvent);

	//setter
    inline void set_generator_quadruplet_before_flip(const GeneratorPtrQuadruplet& generatorQuadruplet) { m_GeneratorQuadrupletBeforeFlip = generatorQuadruplet; };
    inline void set_whether_this_vedge_to_be_flipped(const bool& thisEdgeIsToBeFlipped) { m_DoesThisVEdgeToBeFlipped = thisEdgeIsToBeFlipped; };

	//getter
    inline const GeneratorPtrQuadruplet& get_generator_quadruplet_before_flip() const { return m_GeneratorQuadrupletBeforeFlip; };
    inline bool  get_whether_this_vedge_to_be_flipped() const { return m_DoesThisVEdgeToBeFlipped; };

    inline void initialize_generator_quadruplet_as_NULL() { m_GeneratorQuadrupletBeforeFlip = { NULL, NULL, NULL, NULL }; };

private:
    void copy_from(const EdgeFlipEvent& edgeFlipEvent);
};


}
}


/*

class EdgeFlipEvent : public EventOfDynamicVD2D
{

private:
    VEdge2D*                 m_VEdge;
    bool					 m_doesThisVEdgeToBeFlipped;

public:
    //constructor
    EdgeFlipEvent();
    EdgeFlipEvent(const double& occuringTime);
    EdgeFlipEvent(const double& occuringTime, VEdge2D* vEdge);
    EdgeFlipEvent(const double& occuringTime, VEdge2D* vEdge, const bool& toBeFlipped);
    EdgeFlipEvent(const EdgeFlipEvent& edgeFlippingEvent);

    //deconstructor
    ~EdgeFlipEvent();

    //operator
    EdgeFlipEvent& operator=(const EdgeFlipEvent& edgeFlippingEvent);

    //setter
    inline void set_vedge(VEdge2D* VEdge) { m_VEdge = VEdge; };

    //getter
    inline VEdge2D* get_VEdge() const { return m_VEdge; };
    inline bool     get_whether_this_vedge_to_be_flipped() const { return m_doesThisVEdgeToBeFlipped; };
};


*/

#endif

