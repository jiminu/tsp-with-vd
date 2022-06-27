#ifndef _BUCKET_CELL_FOR_DISKS_H_
#define _BUCKET_CELL_FOR_DISKS_H_

#include "rg_Circle2D.h"
#include <list>
using namespace std;

class BucketCellForDisks
{
private:
	int m_Xindex;
	int m_Yindex;

	list<rg_Circle2D*> m_ContainedDisks;

public:
	BucketCellForDisks();
	BucketCellForDisks(const BucketCellForDisks& bucketCellForDisk);
	 
	~BucketCellForDisks();

	BucketCellForDisks& operator=(const BucketCellForDisks& bucketCellForDisk);

    inline int					        get_Xindex() const { return m_Xindex; };
    inline int					        get_Yindex() const { return m_Yindex; };
    inline const list<rg_Circle2D*>&	get_contained_disks()const { return m_ContainedDisks; } ;

    inline void set_Xindex(const int& xindex) { m_Xindex = xindex; };
    inline void set_Yindex(const int& yindex) { m_Yindex = yindex; };
    inline void set_contained_disks(const list<rg_Circle2D*>& movingDisks) { m_ContainedDisks = movingDisks; };

	inline void insert_disk(rg_Circle2D* disk) { m_ContainedDisks.push_back(disk); };
	inline void remove_disk(rg_Circle2D* disk) { m_ContainedDisks.remove(disk); };
};

#endif //_BUCKET_CELL_FOR_RSS_H_





