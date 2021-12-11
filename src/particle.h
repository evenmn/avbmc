#pragma once
#include <string>
#include <valarray>

// Because the memory of structs has to be
// allocated manually, it is easier to 
// use a class in this particular case
/*
struct Particle 
{
    std::valarray<double> r;
    std::valarray<double> v;
    std::valarray<double> a;

    int type;
    double mass;
    std::string label;
};
*/

class Particle
{
public:
    Particle(std::string, std::valarray<double>);

    std::valarray<double> r;
    std::valarray<double> v;
    std::valarray<double> a;

    int type;
    double mass;
    std::string label;
};
