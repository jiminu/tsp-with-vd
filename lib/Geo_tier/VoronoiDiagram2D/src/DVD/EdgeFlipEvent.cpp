#include "EdgeFlipEvent.h"
using namespace V::GeometryTier;

EdgeFlipEvent::EdgeFlipEvent()
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = -1.0;
    m_DoesThisVEdgeToBeFlipped = true;

    initialize_generator_quadruplet_as_NULL();
}


EdgeFlipEvent::EdgeFlipEvent(const double& occuringTime)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = occuringTime;
    m_DoesThisVEdgeToBeFlipped = true;

    initialize_generator_quadruplet_as_NULL();
}


EdgeFlipEvent::EdgeFlipEvent(const double & occuringTime, 
                             const GeneratorPtrQuadruplet & generatorQuadruplet)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = occuringTime;
    m_DoesThisVEdgeToBeFlipped = true;
    m_GeneratorQuadrupletBeforeFlip      = generatorQuadruplet;

}


EdgeFlipEvent::EdgeFlipEvent(const double & occuringTime, 
                             const GeneratorPtrQuadruplet & generatorQuadruplet, const bool & toBeFlipped)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = occuringTime;
    m_DoesThisVEdgeToBeFlipped = toBeFlipped;
    m_GeneratorQuadrupletBeforeFlip      = generatorQuadruplet;
}


EdgeFlipEvent::EdgeFlipEvent(const EdgeFlipEvent& edgeFlippingEvent)
{
    copy_from(edgeFlippingEvent);
}



EdgeFlipEvent::~EdgeFlipEvent()
{

}



void EdgeFlipEvent::copy_from(const EdgeFlipEvent& edgeFlipEvent)
{
    EventOfDynamicVD2D::copy_from(edgeFlipEvent);

    m_DoesThisVEdgeToBeFlipped       = edgeFlipEvent.m_DoesThisVEdgeToBeFlipped;
    m_GeneratorQuadrupletBeforeFlip  = edgeFlipEvent.m_GeneratorQuadrupletBeforeFlip;
}


EdgeFlipEvent& EdgeFlipEvent::operator=(const EdgeFlipEvent& edgeFlippingEvent)
{
    if (this != &edgeFlippingEvent)
    {
        copy_from(edgeFlippingEvent);
    }

    return *this;
}







/*
EdgeFlipEvent::EdgeFlipEvent()
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = -1.0;
    m_VEdge                    = NULL;
    m_doesThisVEdgeToBeFlipped = true;
}



EdgeFlipEvent::EdgeFlipEvent(const double& occuringTime)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = occuringTime;
    m_VEdge                    = NULL;
    m_doesThisVEdgeToBeFlipped = true;
}



EdgeFlipEvent::EdgeFlipEvent(const double& occuringTime, VEdge2D* vEdge)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = occuringTime;
    m_VEdge                    = vEdge;
    m_doesThisVEdgeToBeFlipped = true;
}


EdgeFlipEvent::EdgeFlipEvent(const double& occuringTime, VEdge2D* vEdge, const bool& toBeFlipped)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = occuringTime;
    m_VEdge                    = vEdge;
    m_doesThisVEdgeToBeFlipped = toBeFlipped;
}



EdgeFlipEvent::EdgeFlipEvent(const EdgeFlipEvent& edgeFlippingEvent)
{
    m_eventType                = EDGE_FLIP;
    m_occurringTime            = edgeFlippingEvent.m_occurringTime;
    m_VEdge                    = edgeFlippingEvent.m_VEdge;
    m_doesThisVEdgeToBeFlipped = edgeFlippingEvent.m_doesThisVEdgeToBeFlipped;
}



EdgeFlipEvent::~EdgeFlipEvent()
{
}



EdgeFlipEvent& EdgeFlipEvent::operator=(const EdgeFlipEvent& edgeFlippingEvent)
{
    if (this == &edgeFlippingEvent)
    {
        return *this;
    }

     m_eventType                = EDGE_FLIP;
     m_occurringTime            = edgeFlippingEvent.m_occurringTime;
     m_VEdge                    = edgeFlippingEvent.m_VEdge;
     m_doesThisVEdgeToBeFlipped = edgeFlippingEvent.m_doesThisVEdgeToBeFlipped;

     return *this;
}



*/
