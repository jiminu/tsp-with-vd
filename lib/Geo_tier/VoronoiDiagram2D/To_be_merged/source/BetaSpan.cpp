#include "BetaSpan.h"
using namespace BULL2D::GeometryTier;

BetaSpan::BetaSpan()
: m_numBetaInterval(0), m_boundingState(rg_NULL), m_betaInterval(rg_NULL)
{

}



BetaSpan::BetaSpan(const rg_INT& numInterval)
: m_numBetaInterval(numInterval)
{
    if ( m_numBetaInterval > 0 ) {
        m_boundingState = new rg_INT[m_numBetaInterval];
        rg_INT i = 0;
        for ( i=0; i<m_numBetaInterval; i++ ) {
            m_boundingState[i] = EXTRANEOUS_SIMPLEX;
        }
        m_betaInterval  = new rg_REAL[m_numBetaInterval+1];
        for ( i=0; i<=m_numBetaInterval; i++ ) {
            m_betaInterval[i] = REAL_MINUS_INFINITE;
        }
    }
    else {
        m_boundingState = rg_NULL;
        m_betaInterval  = rg_NULL;
    }
}



BetaSpan::BetaSpan(const BetaSpan& span)
: m_numBetaInterval(span.m_numBetaInterval)
{
    if ( m_numBetaInterval > 0 ) {
        m_boundingState = new rg_INT[m_numBetaInterval];
        rg_INT i = 0;
        for ( i=0; i<m_numBetaInterval; i++ ) {
            m_boundingState[i] = span.m_boundingState[i];
        }
        m_betaInterval  = new rg_REAL[m_numBetaInterval+1];
        for ( i=0; i<=m_numBetaInterval; i++ ) {
            m_betaInterval[i] = span.m_betaInterval[i];
        }
    }
    else {
        m_boundingState = rg_NULL;
        m_betaInterval  = rg_NULL;
    }
}



BetaSpan::~BetaSpan()
{
    if ( m_boundingState != rg_NULL ) {
        delete [] m_boundingState;
    }
    if ( m_betaInterval != rg_NULL ) {
        delete [] m_betaInterval;
    }
}




rg_INT   BetaSpan::getNumBetaInterval() const
{
    return m_numBetaInterval;
}



rg_INT   BetaSpan::getBoundingState(const rg_INT& i) const
{
    if ( i>=0 && i<m_numBetaInterval ) {
        return m_boundingState[i];
    }
    else {
        return EXTRANEOUS_SIMPLEX;
    }

}



rg_INT*  BetaSpan::getBoundingStates() const
{
    return m_boundingState;
}



rg_REAL* BetaSpan::getBetaIntervals() const
{
    return m_betaInterval;
}



Interval_REAL BetaSpan::getBetaInterval(const rg_INT& i) const
{
    if ( i>=0 && i<m_numBetaInterval ) {
        return Interval_REAL(m_betaInterval[i], m_betaInterval[i+1]);
    }
    else{
        return Interval_REAL();
    }
}



rg_REAL BetaSpan::getLowerValueOfInterval(const rg_INT& i) const
{
    if ( i>=0 && i<m_numBetaInterval ) {
        return m_betaInterval[i];
    }
    else {
        return REAL_MINUS_INFINITE;
    }
}



rg_REAL BetaSpan::getUpperValueOfInterval(const rg_INT& i) const
{
    if ( i>=0 && i<m_numBetaInterval ) {
        return m_betaInterval[i+1];
    }
    else {
        return REAL_MINUS_INFINITE;
    }
}



Interval_REAL BetaSpan::getBetaIntervalOfSingularRegularState() const
{
    rg_INT startInterval = -1;
    rg_INT endInterval   = -1;
    for ( rg_INT i=0; i<m_numBetaInterval; i++ ) {
        if ( m_boundingState[i] == SINGULAR_SIMPLEX || m_boundingState[i] == REGULAR_SIMPLEX ) {
            startInterval = i;
            break;
        }
    }

    for ( rg_INT j=m_numBetaInterval-1; j>=0; j-- ) {
        if ( m_boundingState[j] == SINGULAR_SIMPLEX || m_boundingState[j] == REGULAR_SIMPLEX ) {
            endInterval = j;
            break;
        }
    }

    Interval_REAL betaIntervalOfSingularRegular;
    if ( startInterval != -1 && endInterval != -1 ) {
        betaIntervalOfSingularRegular.setLowerValue( m_betaInterval[startInterval] );
        betaIntervalOfSingularRegular.setUpperValue( m_betaInterval[endInterval+1] );
    }
    else {
        betaIntervalOfSingularRegular.setLowerValue( REAL_MINUS_INFINITE );
        betaIntervalOfSingularRegular.setUpperValue( REAL_MINUS_INFINITE );
    }

    return betaIntervalOfSingularRegular;
}



void BetaSpan::setNumBetaInterVal(const rg_INT& numInterval)
{
    m_numBetaInterval = numInterval;
    if ( m_boundingState != rg_NULL ) {
        delete [] m_boundingState;
        m_boundingState = rg_NULL;
    }
    if ( m_betaInterval != rg_NULL ) {
        delete [] m_betaInterval;
        m_betaInterval = rg_NULL;
    }

    if ( m_numBetaInterval > 0 ) {
        m_boundingState = new rg_INT[m_numBetaInterval];
        rg_INT i = 0;
        for ( i=0; i<m_numBetaInterval; i++ ) {
            m_boundingState[i] = EXTRANEOUS_SIMPLEX;
        }
        m_betaInterval  = new rg_REAL[m_numBetaInterval+1];
        for ( i=0; i<=m_numBetaInterval; i++ ) {
            m_betaInterval[i] = REAL_MINUS_INFINITE;
        }
    }
}



void BetaSpan::setBetaInterval(const rg_INT& i, const rg_REAL& lower, const rg_REAL& upper)
{
    if ( i>=0 && i<m_numBetaInterval ) {
        m_betaInterval[i]   = lower;
        m_betaInterval[i+1] = upper;
    }
}



void BetaSpan::setBetaInterval(const rg_INT& i, const Interval_REAL& interval)
{
    if ( i>=0 && i<m_numBetaInterval ) {
        m_betaInterval[i]   = interval.getLowerValue();
        m_betaInterval[i+1] = interval.getUpperValue();
    }
}



void BetaSpan::setBetaInterval(const rg_INT& i, const rg_INT& boundingState, const rg_REAL& lower, const rg_REAL& upper)
{
    if ( i>=0 && i<m_numBetaInterval ) {
        m_boundingState[i]  = boundingState;

        m_betaInterval[i]   = lower;
        m_betaInterval[i+1] = upper;
    }
}

void BetaSpan::setBetaInterval(const rg_INT& i, const rg_INT& boundingState, const Interval_REAL& interval)
{
    if ( i>=0 && i<m_numBetaInterval ) {
        m_boundingState[i]  = boundingState;

        m_betaInterval[i]   = interval.getLowerValue();
        m_betaInterval[i+1] = interval.getUpperValue();
    }
}


void BetaSpan::shiftBetaSpan(const rg_REAL& delta)
{
    for (rg_INT i=0; i<=m_numBetaInterval; i++ ) 
	{
        m_betaInterval[ i ] = m_betaInterval[ i ] + delta;
	}
}


BetaSpan& BetaSpan::operator=(const BetaSpan& span)
{
    if ( this == &span ) {
        return *this;
    }

    m_numBetaInterval = span.m_numBetaInterval;
    if ( m_boundingState != rg_NULL ) {
        delete [] m_boundingState;
    }
    if ( m_betaInterval != rg_NULL ) {
        delete [] m_betaInterval;
    }
    m_boundingState = rg_NULL;
    m_betaInterval = rg_NULL;


    if ( m_numBetaInterval > 0 ) {
        m_boundingState = new rg_INT[m_numBetaInterval];
        rg_INT i = 0;
        for ( i=0; i<m_numBetaInterval; i++ ) {
            m_boundingState[i] = span.m_boundingState[i];
        }
        m_betaInterval  = new rg_REAL[m_numBetaInterval+1];
        for ( i=0; i<=m_numBetaInterval; i++ ) {
            m_betaInterval[i] = span.m_betaInterval[i];
        }
    }
    
    return *this;
}




rg_INT   BetaSpan::findBoundingState(const rg_REAL& beta) const
{
    rg_INT boundingState = EXTRANEOUS_SIMPLEX;
    if ( m_numBetaInterval <= 0 || beta < m_betaInterval[0] ) {
        return EXTRANEOUS_SIMPLEX;
    }
    if ( beta >= m_betaInterval[m_numBetaInterval] ) {
        return INTERIOR_SIMPLEX;
    }

    for (rg_INT i=0; i<m_numBetaInterval; i++ ) {
        if ( m_betaInterval[i] <= beta && beta < m_betaInterval[i+1] ) {
            boundingState = m_boundingState[i];
            break;
        }
    }

    return boundingState;
}



BetaSpan BetaSpan::merge(const BetaSpan& thatSpan) const
{
    rg_INT   numDividedBetaIntervals = m_numBetaInterval + thatSpan.m_numBetaInterval;
    rg_INT*  dividedBoundingState    = new rg_INT[numDividedBetaIntervals];
    rg_REAL* dividedBetaInterval     = new rg_REAL[numDividedBetaIntervals+1];


    //  divide beta-intervals by lower and upper values 
    //  of each beta-interval in this and that span.
    rg_INT   i = 0;
    rg_INT   i_thisSpan = 0;
    rg_INT   i_thatSpan = 0;
    for ( i=0; i<=numDividedBetaIntervals; i++ ) {
        if ( i_thisSpan<m_numBetaInterval && i_thatSpan<thatSpan.m_numBetaInterval ) {
            if ( m_betaInterval[i_thisSpan] == thatSpan.m_betaInterval[i_thatSpan] ) {
                dividedBetaInterval[i] = m_betaInterval[i_thisSpan];
                i_thisSpan++;
                i_thatSpan++;
            }
            else if ( m_betaInterval[i_thisSpan] > thatSpan.m_betaInterval[i_thatSpan] ) {
                dividedBetaInterval[i] = thatSpan.m_betaInterval[i_thatSpan];
                i_thatSpan++;
            }
            else { //if ( m_betaInterval[i_thisSpan] < thatSpan.m_betaInterval[i_thatSpan] ) {
                dividedBetaInterval[i] = m_betaInterval[i_thisSpan];
                i_thisSpan++;
            }
        }
        else if ( i_thisSpan<m_numBetaInterval && i_thatSpan>=thatSpan.m_numBetaInterval ) {
            dividedBetaInterval[i] = m_betaInterval[i_thisSpan];
            i_thisSpan++;
        }
        else if ( i_thisSpan>=m_numBetaInterval && i_thatSpan<thatSpan.m_numBetaInterval ) {
            dividedBetaInterval[i] = thatSpan.m_betaInterval[i_thatSpan];
            i_thatSpan++;
        }
        else {
            dividedBetaInterval[i] = REAL_PLUS_INFINITE;
        }
    }

    //  determine the bounding states of divided beta-intervals.
    for ( i=0; i<numDividedBetaIntervals; i++ ) {
        rg_INT thisBoundingState = INTERIOR_SIMPLEX;
        rg_INT thatBoundingState = INTERIOR_SIMPLEX;
        rg_INT j=0;
        for ( j=0; j<m_numBetaInterval; j++ ) {
            if ( m_betaInterval[j] <= dividedBetaInterval[i] && dividedBetaInterval[i+1] <= m_betaInterval[j+1] ) {
                thisBoundingState = m_boundingState[j];
                break;
            }
        }
        for ( j=0; j<thatSpan.m_numBetaInterval; j++ ) {
            if ( thatSpan.m_betaInterval[j] <= dividedBetaInterval[i] && dividedBetaInterval[i+1] <= thatSpan.m_betaInterval[j+1] ) {
                thatBoundingState = thatSpan.m_boundingState[j];
                break;
            }
        }



        if ( thisBoundingState == thatBoundingState ) {
            dividedBoundingState[i] = thisBoundingState;
        }
        else {
            if ( thisBoundingState == REGULAR_SIMPLEX || thatBoundingState == REGULAR_SIMPLEX ) {
                dividedBoundingState[i] = REGULAR_SIMPLEX;
            }
            else {
                if ( thisBoundingState == EXTRANEOUS_SIMPLEX ) { 
                    if ( thatBoundingState == SINGULAR_SIMPLEX ) {
                        dividedBoundingState[i] = EXTRANEOUS_SIMPLEX;
                    }
                    else { // if ( thatBoundingState == INTERIOR_SIMPLEX ) {
                        dividedBoundingState[i] = INTERIOR_SIMPLEX;
                    }
                }
                else if ( thisBoundingState == SINGULAR_SIMPLEX ) {
                    if ( thatBoundingState == EXTRANEOUS_SIMPLEX ) {
                        dividedBoundingState[i] = EXTRANEOUS_SIMPLEX;
                    }
                    else { // if ( thatBoundingState[i] == INTERIOR_SIMPLEX ) {
                        dividedBoundingState[i] = SINGULAR_SIMPLEX;
                    }
                }
                else { // if ( thisBoundingState == INTERIOR_SIMPLEX ) {
                    if ( thatBoundingState == EXTRANEOUS_SIMPLEX ) {
                        dividedBoundingState[i] = INTERIOR_SIMPLEX;
                    }
                    else { // if ( thatBoundingState[i] == SINGULAR_SIMPLEX ) {
                        dividedBoundingState[i] = SINGULAR_SIMPLEX;
                    }
                }
            }
        }
    }



    //  make merged beta-span.
    //  1. count the number of beta-intervals of merged beta-span.
    rg_INT numNewBetaIntervals = 1;
    for ( i=0; i<numDividedBetaIntervals-1; i++ ) {
        if ( dividedBoundingState[i] != dividedBoundingState[i+1] ) {
            numNewBetaIntervals++;
        }
    }

    //  2. determine the beta-intervals and their bounding states of merged beta-span.
    rg_INT*  newBoundingState = new rg_INT[numNewBetaIntervals];
    rg_REAL* newBetaInterval  = new rg_REAL[numNewBetaIntervals+1];

    rg_INT currInterval = 0;
    newBoundingState[0] = dividedBoundingState[0];
    newBetaInterval[0]  = dividedBetaInterval[0];
    for ( i=0; i<numDividedBetaIntervals; i++ ) {
        newBetaInterval[currInterval+1] = dividedBetaInterval[i+1];
        if ( newBoundingState[currInterval] != dividedBoundingState[i+1] ) {
            newBoundingState[currInterval+1] = dividedBoundingState[i+1];
            currInterval++;
        }
    }

    delete [] dividedBoundingState;
    delete [] dividedBetaInterval;

    BetaSpan newBetaSpan;
    newBetaSpan.m_numBetaInterval = numNewBetaIntervals;
    newBetaSpan.m_boundingState   = newBoundingState;
    newBetaSpan.m_betaInterval    = newBetaInterval;

    return newBetaSpan;
}



void     BetaSpan::unify(const BetaSpan& span)
{
    BetaSpan newBetaSpan = merge(span);

    (*this) = newBetaSpan;
}



void     BetaSpan::report(ofstream& fout)
{
    if ( m_numBetaInterval <= 0 ) {
        return;
    }

    for ( rg_INT i=0; i<m_numBetaInterval; i++ ) {
        if ( m_boundingState[i] == EXTRANEOUS_SIMPLEX ) {
            fout << "EXT" << "\t";
        }
        else if ( m_boundingState[i] == SINGULAR_SIMPLEX ) {
            fout << "SIN" << "\t";
        }
        else if ( m_boundingState[i] == REGULAR_SIMPLEX ) {
            fout << "REG" << "\t";
        }
        else { //if ( m_boundingState[i] == INTERIOR_SIMPLEX ) {
            fout << "INT" << "\t";
        }
    }

    for ( rg_INT j=0; j<=m_numBetaInterval; j++ ) {
        fout << m_betaInterval[j] << "\t";
    }
    fout << endl;
}


