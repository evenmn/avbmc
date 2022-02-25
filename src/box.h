#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "particle.h"

class Box
{
public:
    Box(class System *); 

    // methods
    void set_boundary(class Boundary *);
    void add_particle(class Particle);
    void add_particle(std::string, std::valarray<double>);
    void add_particles(std::vector<class Particle>);

    std::string file_marking();
    void snapshot(std::string, bool = true);
    void set_dump(int, std::string, std::vector<std::string>, bool = true);
    void set_thermo(int, std::string, std::vector<std::string>, bool = true);
    std::vector<int> build_neigh_list(int, double);
    std::vector<int> build_neigh_list(int, double**);
    void write_nsystemsize(std::string);
    ~Box();

    // variables
    class Dump* dump = nullptr;
    class Thermo* thermo = nullptr;
    class Boundary* boundary = nullptr;
    //class Velocity* velocity = nullptr;
    class System* system = nullptr;
    class DistanceManager* distance_manager = nullptr;

    unsigned int npar, step, ntype, nmove, box_id;
    double poteng, time;

    std::vector<int> nsystemsize;
    std::vector<class Particle> particles;
};
