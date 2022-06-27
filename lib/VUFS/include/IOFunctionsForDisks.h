#ifndef _IO_FUNCTIONS_FOR_DISKS_H
#define _IO_FUNCTIONS_FOR_DISKS_H

#include "rg_Const.h"
#include "DynamicDisk.h"
#include <string>
#include <list>
#include <fstream>
using namespace std;

class IOFunctionsForDisks
{
public:  
    static void read(const string& fileNameWithPath, list<DynamicDisk>& movingDisks);
    static void read(const string& fileNameWithPath, DynamicDisk& container, list<DynamicDisk>& movingDisks);
    static void read(const string& fileNameWithPath, list<DynamicDisk>& container, list<DynamicDisk>& movingDisks);

    static void read(const string& fileNameWithPath, list<rg_Circle2D>& disks);
    static void read(const string& fileNameWithPath, rg_Circle2D& container, list<rg_Circle2D>& disks);

    static void write(const string& fileNameWithPath, list<Disk>& disks);
    static void write(const string& fileNameWithPath, Disk& container, list<Disk>& disks);
        
    static void write(const string& fileNameWithPath, list<DynamicDisk>& dynamicDisks);
};


#endif

