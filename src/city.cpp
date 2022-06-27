#include"city.h"
#include<cmath>

City::City() {
    this->x = 0;
    this->y = 0;
}

City::City(const float& x, const float& y) {
    this->x = x;
    this->y = y;
}

float City::distance_to(const City& p1) {
    float distance;

    distance = sqrt(pow(this->x - p1.x, 2) + pow(this->y - p1.y, 2));

    return distance;
}

bool City::operator==(const City& rhs) const {
    return (x == rhs.x && y == rhs.y);
}