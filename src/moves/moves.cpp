#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

#include "moves.h"
#include "../box.h"


/* --------------------------------------------------------
   Moves base class constructor
----------------------------------------------------------- */

Moves::Moves(Box* box_in)
{
    box = box_in;
    rng = box->rng;
}


/* --------------------------------------------------------
   Rotate molecule by a random angle in all dimensions. Might
   send in positions as a reference in the future, even
   though we still may need to copy the positions
----------------------------------------------------------- */

std::vector<std::valarray<double> > Moves::rotate_molecule(std::vector<std::valarray<double> > positions_in)
{
    std::vector<std::valarray<double> > positions_out;
    if(box->ndim == 1){
        return positions_in;
    }
    else if(box->ndim == 2){
        double angle = 2 * pi * box->rng->next_double();
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
        double anglex = 2 * pi * box->rng->next_double();
        double angley = 2 * pi * box->rng->next_double();
        double anglez = 2 * pi * box->rng->next_double();
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

