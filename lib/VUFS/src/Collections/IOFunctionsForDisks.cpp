#include "IOFunctionsForDisks.h"
#include "string.h"

void IOFunctionsForDisks::read(const string& fileNameWithPath, list<DynamicDisk>& movingDisks)
{
    list<DynamicDisk> containers;
    IOFunctionsForDisks::read(fileNameWithPath, containers, movingDisks);

    if (!containers.empty())
    {
        movingDisks.insert(movingDisks.begin(), containers.begin(), containers.end());
    }
}


void IOFunctionsForDisks::read(const string& fileNameWithPath, DynamicDisk& container, list<DynamicDisk>& movingDisks)
{
    list<DynamicDisk> containers;
    IOFunctionsForDisks::read(fileNameWithPath, containers, movingDisks);

    if (containers.size() == 0)
    {
        container = DynamicDisk();
    }
    else if (containers.size() == 1)
    {
        container = containers.front();
    }
    else
    {
        DynamicDisk containerHavingLargestNegativeId = containers.front();
        for (list<DynamicDisk>::const_iterator it_Container = next(containers.begin());
            it_Container != containers.end();
            ++it_Container)
        {
            const DynamicDisk& currContainer = *it_Container;
            if (currContainer.getID() > containerHavingLargestNegativeId.getID())
            {
                containerHavingLargestNegativeId = currContainer;
            }
        }

        container = containerHavingLargestNegativeId;
    }
}

void IOFunctionsForDisks::read(const string& fileNameWithPath, list<DynamicDisk>& containers, list<DynamicDisk>& movingDisks)
{
    ifstream fin;
    fin.open(fileNameWithPath.c_str());

    if (!fin.is_open())
    {
        return;
    }

    char  seps[4]    = { ' ', '\t', '\n', '\0' };
    char* token      = NULL;
    char* next_token = NULL;
    char buffer[300];

    while (fin.getline(buffer, 300)) 
    {
        if (buffer[0] != '#')
        {
            break;
        }
    };

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    rg_INT numMovingDisks = atoi(strtok_s(buffer, seps, &next_token));

    for (rg_INT i = 0; i < numMovingDisks; i++) {
        fin.getline(buffer, 300);
        rg_INT   ID = atoi(strtok_s(buffer, seps, &next_token));
        rg_REAL  x = atof(strtok_s(NULL, seps, &next_token));
        rg_REAL  y = atof(strtok_s(NULL, seps, &next_token));
        rg_REAL  radius = atof(strtok_s(NULL, seps, &next_token));


        char* redString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT red = 0.5f;
        rg_FLOAT green = 0.5f;
        rg_FLOAT blue = 0.5f;
        if (redString != rg_NULL) {
            red = atof(redString);
            green = atof(strtok_s(NULL, seps, &next_token));
            blue = atof(strtok_s(NULL, seps, &next_token));
        }


        char* vecXString = strtok_s(NULL, seps, &next_token);
        rg_REAL  vecX;
        rg_REAL  vecY;

        if (/*BCFVersion == 2 && +*/vecXString != rg_NULL) {
            vecX = atof(vecXString);
            vecY = atof(strtok_s(NULL, seps, &next_token));
        }

        char* redBoundayString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT redBoundary = 0.0f;
        rg_FLOAT greenBoundary = 0.0f;
        rg_FLOAT blueBoundary = 0.0f;;
        if (/*BCFVersion == 2 && +*/redBoundayString != rg_NULL) {
            redBoundary = atof(redBoundayString);
            greenBoundary = atof(strtok_s(NULL, seps, &next_token));
            blueBoundary = atof(strtok_s(NULL, seps, &next_token));
        }

        char* redCentPTString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT redCenterPT = 0.0f;
        rg_FLOAT greenCenterPT = 0.0f;
        rg_FLOAT blueCenterPT = 0.0f;;
        if (/*BCFVersion == 2 && +*/redCentPTString != rg_NULL) {
            redCenterPT = atof(redCentPTString);
            greenCenterPT = atof(strtok_s(NULL, seps, &next_token));
            blueCenterPT = atof(strtok_s(NULL, seps, &next_token));
        }

        rg_Circle2D disk(x, y, radius);
        rg_Point2D  motionVector(vecX, vecY);

        DynamicDisk dynamicDisk(ID, disk, motionVector, 0.0);

        if (ID < 0)
        {
            containers.push_back(dynamicDisk);
        }
        else
        {
            movingDisks.push_back(dynamicDisk);
        }
    }

    fin.close();
#else
    rg_INT numMovingDisks = atoi(strtok_r(buffer, seps, &next_token));

    for (rg_INT i = 0; i < numMovingDisks; i++) {
        fin.getline(buffer, 300);
        rg_INT   ID = atoi(strtok_r(buffer, seps, &next_token));
        rg_REAL  x = atof(strtok_r(NULL, seps, &next_token));
        rg_REAL  y = atof(strtok_r(NULL, seps, &next_token));
        rg_REAL  radius = atof(strtok_r(NULL, seps, &next_token));


        char* redString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT red = 0.5f;
        rg_FLOAT green = 0.5f;
        rg_FLOAT blue = 0.5f;
        if (redString != rg_NULL) {
            red = atof(redString);
            green = atof(strtok_r(NULL, seps, &next_token));
            blue = atof(strtok_r(NULL, seps, &next_token));
        }


        char* vecXString = strtok_r(NULL, seps, &next_token);
        rg_REAL  vecX;
        rg_REAL  vecY;

        if (/*BCFVersion == 2 && +*/vecXString != rg_NULL) {
            vecX = atof(vecXString);
            vecY = atof(strtok_r(NULL, seps, &next_token));
        }

        char* redBoundayString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT redBoundary = 0.0f;
        rg_FLOAT greenBoundary = 0.0f;
        rg_FLOAT blueBoundary = 0.0f;;
        if (/*BCFVersion == 2 && +*/redBoundayString != rg_NULL) {
            redBoundary = atof(redBoundayString);
            greenBoundary = atof(strtok_r(NULL, seps, &next_token));
            blueBoundary = atof(strtok_r(NULL, seps, &next_token));
        }

        char* redCentPTString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT redCenterPT = 0.0f;
        rg_FLOAT greenCenterPT = 0.0f;
        rg_FLOAT blueCenterPT = 0.0f;;
        if (/*BCFVersion == 2 && +*/redCentPTString != rg_NULL) {
            redCenterPT = atof(redCentPTString);
            greenCenterPT = atof(strtok_r(NULL, seps, &next_token));
            blueCenterPT = atof(strtok_r(NULL, seps, &next_token));
        }

        rg_Circle2D disk(x, y, radius);
        rg_Point2D  motionVector(vecX, vecY);

        DynamicDisk dynamicDisk(ID, disk, motionVector, 0.0);

        if (ID < 0)
        {
            containers.push_back(dynamicDisk);
        }
        else
        {
            movingDisks.push_back(dynamicDisk);
        }
    }

    fin.close();

#endif
}





void IOFunctionsForDisks::read(const string& fileNameWithPath, list<rg_Circle2D>& disks)
{
 ifstream fin;
    fin.open(fileNameWithPath.c_str());

    char  seps[4]    = { ' ', '\t', '\n', '\0' };
    char* token      = NULL;
    char* next_token = NULL;
    char buffer[300];

    while (fin.getline(buffer, 300)) 
    {
        if (buffer[0] != '#') 
        {
            break;
        }
    };

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    rg_INT numMovingDisks = atoi(strtok_s(buffer, seps, &next_token));

    for (rg_INT i = 0; i < numMovingDisks; i++) {
        fin.getline(buffer, 300);
        rg_INT   ID     = atoi(strtok_s(buffer, seps, &next_token));
        rg_REAL  x      = atof(strtok_s(NULL, seps, &next_token));
        rg_REAL  y      = atof(strtok_s(NULL, seps, &next_token));
        rg_REAL  radius = atof(strtok_s(NULL, seps, &next_token));

        char* redString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT red = 0.5f;
        rg_FLOAT green = 0.5f;
        rg_FLOAT blue = 0.5f;
        if (redString != rg_NULL) {
            red   = atof(redString);
            green = atof(strtok_s(NULL, seps, &next_token));
            blue  = atof(strtok_s(NULL, seps, &next_token));
        }


        char* vecXString = strtok_s(NULL, seps, &next_token);
        rg_REAL  vecX;
        rg_REAL  vecY;

        if (/*BCFVersion == 2 && +*/vecXString != rg_NULL) {
            vecX = atof(vecXString);
            vecY = atof(strtok_s(NULL, seps, &next_token));
        }

        char* redBoundayString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT redBoundary   = 0.0f;
        rg_FLOAT greenBoundary = 0.0f;
        rg_FLOAT blueBoundary  = 0.0f;;
        if (/*BCFVersion == 2 && +*/redBoundayString != rg_NULL) {
            redBoundary   = atof(redBoundayString);
            greenBoundary = atof(strtok_s(NULL, seps, &next_token));
            blueBoundary  = atof(strtok_s(NULL, seps, &next_token));
        }

        char* redCentPTString  = strtok_s(NULL, seps, &next_token);
        rg_FLOAT redCenterPT   = 0.0f;
        rg_FLOAT greenCenterPT = 0.0f;
        rg_FLOAT blueCenterPT  = 0.0f;;
        if (/*BCFVersion == 2 && +*/redCentPTString != rg_NULL) {
            redCenterPT   = atof(redCentPTString);
            greenCenterPT = atof(strtok_s(NULL, seps, &next_token));
            blueCenterPT  = atof(strtok_s(NULL, seps, &next_token));
        }

        rg_Circle2D disk(x, y, radius);
 
        disks.push_back(disk);
    }

    fin.close();
#else
     rg_INT numMovingDisks = atoi(strtok_r(buffer, seps, &next_token));

    for (rg_INT i = 0; i < numMovingDisks; i++) {
        fin.getline(buffer, 300);
        rg_INT   ID     = atoi(strtok_r(buffer, seps, &next_token));
        rg_REAL  x      = atof(strtok_r(NULL, seps, &next_token));
        rg_REAL  y      = atof(strtok_r(NULL, seps, &next_token));
        rg_REAL  radius = atof(strtok_r(NULL, seps, &next_token));

        char* redString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT red = 0.5f;
        rg_FLOAT green = 0.5f;
        rg_FLOAT blue = 0.5f;
        if (redString != rg_NULL) {
            red   = atof(redString);
            green = atof(strtok_r(NULL, seps, &next_token));
            blue  = atof(strtok_r(NULL, seps, &next_token));
        }


        char* vecXString = strtok_r(NULL, seps, &next_token);
        rg_REAL  vecX;
        rg_REAL  vecY;

        if (/*BCFVersion == 2 && +*/vecXString != rg_NULL) {
            vecX = atof(vecXString);
            vecY = atof(strtok_r(NULL, seps, &next_token));
        }

        char* redBoundayString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT redBoundary   = 0.0f;
        rg_FLOAT greenBoundary = 0.0f;
        rg_FLOAT blueBoundary  = 0.0f;;
        if (/*BCFVersion == 2 && +*/redBoundayString != rg_NULL) {
            redBoundary   = atof(redBoundayString);
            greenBoundary = atof(strtok_r(NULL, seps, &next_token));
            blueBoundary  = atof(strtok_r(NULL, seps, &next_token));
        }

        char* redCentPTString  = strtok_r(NULL, seps, &next_token);
        rg_FLOAT redCenterPT   = 0.0f;
        rg_FLOAT greenCenterPT = 0.0f;
        rg_FLOAT blueCenterPT  = 0.0f;;
        if (/*BCFVersion == 2 && +*/redCentPTString != rg_NULL) {
            redCenterPT   = atof(redCentPTString);
            greenCenterPT = atof(strtok_r(NULL, seps, &next_token));
            blueCenterPT  = atof(strtok_r(NULL, seps, &next_token));
        }

        rg_Circle2D disk(x, y, radius);
 
        disks.push_back(disk);
    }

    fin.close();
#endif
}


void IOFunctionsForDisks::read(const string& fileNameWithPath, rg_Circle2D& container, list<rg_Circle2D>& disks)
{
    ifstream fin;
    fin.open(fileNameWithPath.c_str());

    char  seps[3]    = { ' ', '\t', '\n' };
    char* token      = NULL;
    char* next_token = NULL;
    char buffer[300];

    while (fin.getline(buffer, 300)) 
    {
        if (buffer[0] != '#') 
        {
            break;
        }
    };

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    rg_INT numMovingDisks = atoi(strtok_s(buffer, seps, &next_token));

    for (rg_INT i = 0; i < numMovingDisks; i++) {
        fin.getline(buffer, 300);
        rg_INT   ID = atoi(strtok_s(buffer, seps, &next_token));
        rg_REAL  x = atof(strtok_s(NULL, seps, &next_token));
        rg_REAL  y = atof(strtok_s(NULL, seps, &next_token));
        rg_REAL  radius = atof(strtok_s(NULL, seps, &next_token));

        char* redString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT red = 0.5f;
        rg_FLOAT green = 0.5f;
        rg_FLOAT blue = 0.5f;
        if (redString != rg_NULL) {
            red = atof(redString);
            green = atof(strtok_s(NULL, seps, &next_token));
            blue = atof(strtok_s(NULL, seps, &next_token));
        }


        char* vecXString = strtok_s(NULL, seps, &next_token);
        rg_REAL  vecX;
        rg_REAL  vecY;

        if (/*BCFVersion == 2 && +*/vecXString != rg_NULL) {
            vecX = atof(vecXString);
            vecY = atof(strtok_s(NULL, seps, &next_token));
        }

        char* redBoundayString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT redBoundary = 0.0f;
        rg_FLOAT greenBoundary = 0.0f;
        rg_FLOAT blueBoundary = 0.0f;;
        if (/*BCFVersion == 2 && +*/redBoundayString != rg_NULL) {
            redBoundary = atof(redBoundayString);
            greenBoundary = atof(strtok_s(NULL, seps, &next_token));
            blueBoundary = atof(strtok_s(NULL, seps, &next_token));
        }

        char* redCentPTString = strtok_s(NULL, seps, &next_token);
        rg_FLOAT redCenterPT = 0.0f;
        rg_FLOAT greenCenterPT = 0.0f;
        rg_FLOAT blueCenterPT = 0.0f;;
        if (/*BCFVersion == 2 && +*/redCentPTString != rg_NULL) {
            redCenterPT = atof(redCentPTString);
            greenCenterPT = atof(strtok_s(NULL, seps, &next_token));
            blueCenterPT = atof(strtok_s(NULL, seps, &next_token));
        }

        rg_Circle2D disk(x, y, radius);
        
        if (ID == -1)
        {
            container = disk;
        }
        else
        {
            disks.push_back(disk);
        }
    }

    fin.close();

#else
    rg_INT numMovingDisks = atoi(strtok_r(buffer, seps, &next_token));

    for (rg_INT i = 0; i < numMovingDisks; i++) {
        fin.getline(buffer, 300);

        rg_INT   ID     = atoi(strtok_r(buffer, seps, &next_token));
        rg_REAL  x      = atof(strtok_r(NULL, seps, &next_token));
        rg_REAL  y      = atof(strtok_r(NULL, seps, &next_token));
        rg_REAL  radius = atof(strtok_r(NULL, seps, &next_token));

        char* redString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT red = 0.5f;
        rg_FLOAT green = 0.5f;
        rg_FLOAT blue = 0.5f;
        if (redString != rg_NULL) {
            red = atof(redString);
            green = atof(strtok_r(NULL, seps, &next_token));
            blue = atof(strtok_r(NULL, seps, &next_token));
        }

        char* vecXString = strtok_r(NULL, seps, &next_token);
        rg_REAL  vecX;
        rg_REAL  vecY;

        if (/*BCFVersion == 2 && +*/vecXString != rg_NULL) {
            vecX = atof(vecXString);
            vecY = atof(strtok_r(NULL, seps, &next_token));
        }

        char* redBoundayString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT redBoundary = 0.0f;
        rg_FLOAT greenBoundary = 0.0f;
        rg_FLOAT blueBoundary = 0.0f;;
        if (/*BCFVersion == 2 && +*/redBoundayString != rg_NULL) {
            redBoundary = atof(redBoundayString);
            greenBoundary = atof(strtok_r(NULL, seps, &next_token));
            blueBoundary = atof(strtok_r(NULL, seps, &next_token));
        }

        char* redCentPTString = strtok_r(NULL, seps, &next_token);
        rg_FLOAT redCenterPT = 0.0f;
        rg_FLOAT greenCenterPT = 0.0f;
        rg_FLOAT blueCenterPT = 0.0f;;
        if (/*BCFVersion == 2 && +*/redCentPTString != rg_NULL) {
            redCenterPT   = atof(redCentPTString);
            greenCenterPT = atof(strtok_r(NULL, seps, &next_token));
            blueCenterPT  = atof(strtok_r(NULL, seps, &next_token));
        }


        rg_Circle2D disk(x, y, radius);

        if (ID == -1)
        {
            container = disk;
        }
        else
        {
            disks.push_back(disk);
        }
    }

    fin.close();
#endif
}


void IOFunctionsForDisks::write(const string& fileNameWithPath, list<Disk>& disks)
{
    ofstream fout(fileNameWithPath);

    fout << disks.size() << endl;

    list<Disk>::iterator i_disk = disks.begin();
    for (; i_disk != disks.end(); ++i_disk)
    {
        rg_INT id = i_disk->getID();
        rg_Point2D center = i_disk->getCenterPt();
        rg_REAL radius = i_disk->getRadius();
        fout << id++ << " " << center.getX() << " " << center.getY() << " " << radius << endl;
    }

    fout.close();
}


void IOFunctionsForDisks::write(const string& fileNameWithPath, Disk& container, list<Disk>& disks)
{
    ofstream fout(fileNameWithPath);

    fout << disks.size() + 1 << endl;

    rg_INT id = container.getID();
    rg_Point2D center = container.getCenterPt();
    rg_REAL radius = container.getRadius();
    fout << id++ << " " << center.getX() << " " << center.getY() << " " << radius << endl;
        
    list<Disk>::iterator i_disk = disks.begin();
    for (; i_disk != disks.end(); ++i_disk)
    {
        rg_INT id = i_disk->getID();
        rg_Point2D center = i_disk->getCenterPt();
        rg_REAL radius = i_disk->getRadius();
        fout << id << " " << center.getX() << " " << center.getY() << " " << radius << endl;
    }

    fout.close();
}


void IOFunctionsForDisks::write(const string& fileNameWithPath, list<DynamicDisk>& dynamicDisks)
{
    ofstream fout(fileNameWithPath);
    fout << "# BCF Ver 2.1" << endl;
    fout << "# The first line is the number of circles." << endl;
    fout << "# The lines from the second to the last are the data for each circle where each has 15 columns as follows:" << endl;
    fout << "# 1. Circle ID, 2. X-coordinate, 3. Y-coordinate, 4. radius," << endl;
    fout << "# 5. R-coord, 6. G-coord, 7. B-coord of color for face," << endl;
    fout << "# 8. X-component of motion vector, 9. Y-component of motion vector," << endl;
    fout << "# 10. R-coord, 11. G-coord, 12. B-coord of color for boundary," << endl;
    fout << "# 13. R-coord, 14. G-coord, and 15. B-coord of color for center point" << endl;
    fout << "# The first four fields are required and the rest eleven fields are optional." << endl;
    fout << "# If the R, G, and B-coords of color for face are not defined, a default color (grey: R=0.686275, G=0.686275, B=0.686275) is assumed." << endl;
    fout << "# If the X and Y components of motion vector are not defined, the circle is static without any motion." << endl;
    fout << "# If a motion is assigned to a circle, the R, G, B coordinates of color for face should also be defined." << endl;
    fout << "# If the R, G, and B-coords of color for boundary/center point are not defined, a default color (grey: R=1.0, G=1.0, B=1.0) is assumed." << endl;
    fout << "# If the R, G, and B-coords of color center point, the R, G, B coordinates of color for boundary should also be defined." << endl;
    fout << "# Either \"tab\" or \"blank\" can be used as a delimiter. " << endl;
    fout << "# " << endl;
    fout << "# If Id is -1, the circle is container including the others" << endl;


    fout << dynamicDisks.size() << endl;

    for (list<DynamicDisk>::iterator it_DynamicDisk = dynamicDisks.begin(); 
         it_DynamicDisk != dynamicDisks.end(); 
         ++it_DynamicDisk)
    {
        const DynamicDisk& dynamicDisk = *it_DynamicDisk;

        int     id     = dynamicDisk.getID();
        double  xCoord = dynamicDisk.getCenterPt().getX();
        double  yCoord = dynamicDisk.getCenterPt().getY();
        double  radius = dynamicDisk.getRadius();
        double  faceColor_R = 0.686275;
        double  faceColor_G = 0.686275;
        double  faceColor_B = 0.686275;
        double  xSpeed = dynamicDisk.getVelocityVectorX();
        double  ySpeed = dynamicDisk.getVelocityVectorY();
        double  vertexColor_R = 0.0;
        double  vertexColor_G = 0.0;
        double  vertexColor_B = 0.0;
        double  edgeColor_R   = 0.0;
        double  edgeColor_G   = 0.0;
        double  edgeColor_B   = 0.0;

        fout << id     << "\t"
             << xCoord << "\t"
             << yCoord << "\t"
             << radius << "\t"
             << faceColor_R << "\t"
             << faceColor_G << "\t"
             << faceColor_B << "\t"
             << xSpeed      << "\t"
             << ySpeed      << "\t"
             << vertexColor_R << "\t"
             << vertexColor_G << "\t"
             << vertexColor_B << "\t"
             << edgeColor_R   << "\t"
             << edgeColor_G   << "\t"
             << edgeColor_B   << endl;
    }

    fout.close();
}

