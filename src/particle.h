#pragma once
#include <string>
#include <valarray>


class Particle
{
public:
    Particle(std::string, std::valarray<double>);

    std::valarray<double> r;
    std::valarray<double> v;
    std::valarray<double> f;

    int type;
    double mass;
    std::string label;
};
