#pragma once

class City {
public :
    // member variable
    float x;
    float y;

    // method
    City();
    City(const float& x, const float& y);
    ~City() {};
    
    float distance_to(const City& c1);
    
    bool operator==(const City& rhs) const;
};