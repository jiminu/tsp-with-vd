#ifndef _QTGATELIST_H
#define _QTGATELIST_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "rg_QTGate.h"

namespace V {

namespace GeometryTier {


class QTGateList
{
private:
    rg_dList<QTGate> m_gates;

public:
    QTGateList();
    QTGateList(QTGateList& gateList);
    ~QTGateList();
    
    rg_INT            getSize() const;
    rg_dList<QTGate>* getGates();

    QTGate* addGate(const QTGate& gate);
    QTGate* find(const QTGate& gate);

    QTGateList& operator =(QTGateList& gateList);
};

} // namespace GeometryTier

} // namespace V


#endif
