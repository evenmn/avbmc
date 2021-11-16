#pragma once
#include <string>
#include <valarray>

struct Particle {
    std::valarray<double> r;
    std::valarray<double> v;
    std::valarray<double> a;

    int type;
    double mass;
    std::string label;
};
