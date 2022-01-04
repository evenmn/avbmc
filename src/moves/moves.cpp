#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

#include "moves.h"
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
    for (int i=0; i < array.size(); i++)
    {
        normsq += array[i] * array[i];
    }
    return normsq;
}


/*
Moves::~Moves()
{
    delete system;
    delete rng;
}
*/
