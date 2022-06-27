#include "rg_QTGateList.h"
using namespace V::GeometryTier;

QTGateList::QTGateList()
{
}



QTGateList::QTGateList(QTGateList& gateList)
{
    m_gates.duplicateList(gateList.m_gates);
}



QTGateList::~QTGateList()
{
}




rg_INT QTGateList::getSize() const
{
    return m_gates.getSize();
}



rg_dList<QTGate>* QTGateList::getGates()
{
    return &m_gates;
}




QTGate* QTGateList::addGate(const QTGate& gate)
{
    return m_gates.add(gate);
}



QTGate* QTGateList::find(const QTGate& gate)
{
    QTGate* currGate = rg_NULL;

    m_gates.reset4Loop();
    while ( m_gates.setNext4Loop() )  {
        currGate = m_gates.getpEntity();

        if ( *currGate == gate )
            return currGate;
    }
    return rg_NULL;
}




QTGateList& QTGateList::operator =(QTGateList& gateList)
{
    if ( this == &gateList )
        return *this;

    m_gates.duplicateList(gateList.m_gates);

    return *this;
}



