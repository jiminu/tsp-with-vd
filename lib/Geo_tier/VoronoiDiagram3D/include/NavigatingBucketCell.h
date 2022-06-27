#ifndef _NAVIGATINGBUCKETCELL_H
#define _NAVIGATINGBUCKETCELL_H

#include "rg_Const.h"
#include "BucketCellIndex.h"
#include "rg_Point3D.h"


namespace V {

namespace GeometryTier {


const rg_INT I_MINUS = 1;
const rg_INT J_MINUS = 2;
const rg_INT K_MINUS = 3;
const rg_INT I_PLUS  = 4;
const rg_INT J_PLUS  = 5;
const rg_INT K_PLUS  = 6;

const rg_INT PREV_BUCKETCELL[6]    = { I_PLUS, J_PLUS, K_PLUS, I_MINUS, J_MINUS, K_MINUS };
const rg_INT NEXT_BUCKETCELL[6]    = { I_MINUS, J_MINUS, K_MINUS, I_PLUS, J_PLUS, K_PLUS };
const rg_INT BUCKETCELL_FACE[6][4] = { {0, 3, 4, 7},    //  i- of curr cell
                                       {0, 1, 4, 5},    //  j- of curr cell
                                       {0, 1, 2, 3},    //  k- of curr cell
                                       {1, 2, 5, 6},    //  i+ of curr cell
                                       {2, 3, 6, 7},    //  j+ of curr cell
                                       {4, 5, 6, 7} };  //  k+ of curr cell

// cellMask&BUCKETCELL_FACE_FILTER[i] (bit operation)
//   -> mask of face[i]
const rg_INT BUCKETCELL_FACE_FILTER[6] = { 50115, 3855, 255, 15420, 61680, 65280 };

//  sign of distances of all vertices on face[i] >0,    <0,    =0
const rg_INT BUCKETCELL_FACE_STATUS[6][3] = {   { 0, 16705, 33410},
                                                { 0,  1285,  2570},
                                                { 0,    85,   170},
                                                { 0,  5140, 10280},
                                                { 0, 20560, 41120},
                                                { 0, 21760, 43520}  };
//  mask values for dist(Pg, vi)           >0,     <0,     =0
const rg_INT BUCKETVERTEX_MASK[8][3] = { { 0,     1,     2}, 
                                         { 0,     4,     8}, 
                                         { 0,    16,    32}, 
                                         { 0,    64,   128},
                                         { 0,   256,   512}, 
                                         { 0,  1024,  2048}, 
                                         { 0,  4096,  8192}, 
                                         { 0, 16384, 32768} };





class NavigatingBucketCell
{
private:
public:
    BucketCellIndex m_index;
    rg_FLAG         m_status;
    rg_REAL         m_distBetPandVi;
    rg_FLAG         m_prevCell;
    rg_INT          m_knownFaceMask[4];

    NavigatingBucketCell();
    NavigatingBucketCell(const BucketCellIndex& index, 
                         const rg_FLAG&         status); 
    NavigatingBucketCell(const BucketCellIndex& index, 
                         const rg_FLAG&         status, 
                         const rg_REAL&         dist  );
    NavigatingBucketCell(const BucketCellIndex& index, 
                         const rg_FLAG&         status, 
                         const rg_REAL&         dist,
                         const rg_FLAG&         prevCell,
                               rg_INT*          knownFaceMask);
    NavigatingBucketCell(const NavigatingBucketCell& naviCell);
    ~NavigatingBucketCell();

    BucketCellIndex getIndex() const;

    rg_FLAG getStatus() const;
    rg_REAL getDistanceVi() const;
    rg_FLAG getPrevCell() const;
    rg_INT  getKnownFaceMask(const rg_INT& i) const;

    void    setIndex(const BucketCellIndex& index);
    void    setStatus(const rg_FLAG& status);
    void    setDistanceVi(const rg_REAL& dist);
    void    setPrevCell(const rg_FLAG& prevCell);
    void    setKnownFaceMask(rg_INT* knownFaceMask);

    void    setCell(const BucketCellIndex& index, 
                    const rg_FLAG&         status, 
                    const rg_REAL&         dist,
                    const rg_FLAG&         prevCell,
                          rg_INT*          knownFaceMask);


    NavigatingBucketCell& operator =(const NavigatingBucketCell& naviCell);


    void computeVisFromPrevCell(rg_REAL* distVis, const rg_Point3D& unitLength) const;
};

} // namespace GeometryTier

} // namespace V


#endif

