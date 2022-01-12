#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

#include "moves.h"
#include "../particle.h"
#include "../system.h"
#include "../rng/rng.h"


/* --------------------------------------------------------
   Moves base class constructor
----------------------------------------------------------- */

Moves::Moves(System* system_in)
{
    system = system_in;
    rng = system->rng;
}


/* --------------------------------------------------------
   Rotate molecule by a random angle in all dimensions. Might
   send in positions as a reference in the future, even
   though we still may need to copy the positions
----------------------------------------------------------- */

std::vector<std::valarray<double> > Moves::rotate_molecule(std::vector<std::valarray<double> > positions_in)
{
    std::vector<std::valarray<double> > positions_out;
    if(system->ndim == 1){
        return positions_in;
    }
    else if(system->ndim == 2){
        double angle = 2 * pi * rng->next_double();
        for(std::valarray<double> position_in : positions_in){
            std::valarray<double> position_out = {
                position_in[0] * std::cos(angle) - position_in[1] * std::sin(angle),
                position_in[0] * std::sin(angle) + position_in[1] * std::cos(angle)
            };
            positions_out.push_back(position_out);
        }
        return positions_out;
    }
    else{
        double anglex = 2 * pi * rng->next_double();
        double angley = 2 * pi * rng->next_double();
        double anglez = 2 * pi * rng->next_double();
        for(std::valarray<double> position_in : positions_in){
            std::valarray<double> position_out = {
                position_in[0] * (1 + std::cos(angley) + std::cos(anglez))
                    - position_in[2] * std::sin(anglez)
                    + position_in[1] * std::sin(angley),
                position_in[1] * (1 + std::cos(anglex) + std::cos(anglez))
                    + position_in[0] * std::sin(angley)
                    - position_in[2] * std::sin(anglex),
                position_in[2] * (1 + std::cos(anglex) + std::cos(angley))
                    - position_in[0] * std::sin(angley)
                    + position_in[1] * std::sin(anglex)
            };
            positions_out.push_back(position_out);
        }
        return positions_out;
    }
}


/* -------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
---------------------------------------------------------------- */

double Moves::norm(std::valarray<double> array)
{
    double normsq = 0.;
    for (unsigned int i=0; i < array.size(); i++)
    {
        normsq += array[i] * array[i];
    }
    return normsq;
}


/* -------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor
   distance squared 'rsq'
---------------------------------------------------------------- */

std::vector<int> Moves::build_neigh_list(std::vector<Particle> particles, const int i, const double rsq)
{
    unsigned int npar = particles.size();
    double rijsq;
    std::valarray<double> ri = particles[i].r;
    std::vector<int> neigh_list;
    for(unsigned int j=0; j<i; j++){
        rijsq = std::pow(particles[j].r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for(unsigned int j=i+1; j<npar; j++){
        rijsq = std::pow(particles[j].r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}



/* ---------------------------------------------------------------
   Check if particles match molecule type recursively.
------------------------------------------------------------------ */

void Moves::check_neighbors(const int k, std::vector<Particle> molecule, unsigned int elm_count,
                            std::vector<int> &elm_idx, std::vector<Particle> particles, double rc) {
    if (elm_count <= molecule.size()) {  // ensure that recursion stops when molecule has correct size
        if (particles[k].type == molecule[elm_count].type) {
            elm_idx.push_back(k);
            elm_count ++;
            std::vector<int> neigh_list = build_neigh_list(particles, k, rc);
            for (int neigh : neigh_list) {
                check_neighbors(neigh, molecule, elm_count, elm_idx, particles, rc);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule'. randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
------------------------------------------------------------------------------- */

std::vector<int> Moves::detect_molecule(std::vector<Particle> particles,
                                        std::vector<Particle> molecule,
                                        bool &detected, double rc)
{
    std::vector<int> elm_idx;
    unsigned int count = 0;
    while (count < particles.size() || !detected)
    {
        elm_idx.clear();
        int k = system->rng->next_int(particles.size());     // pick initial particle
        check_neighbors(k, molecule, 0, elm_idx, particles, rc);
        if (elm_idx.size() == particles.size()) {
            detected = true;
        }
        count ++;
    }
    if (!detected)
    {
        elm_idx.clear();
    }
    return elm_idx;
}


