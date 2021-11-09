#include "init_velocity.h"
#include "../box.h"


Velocity::Velocity(class Box* box_in)
{
    box = box_in;
}

Gauss::Gauss(class Box* box_in, const double mean_in, const double var_in) 
    : Velocity(box_in)
{
    mean = mean_in;
    var = var_in;
}

mat Gauss::get_velocity(const int npar, const int ndim)
{
    mat velocity(npar, ndim);
    for(int i=0; i<npar; i++){
        for(int j=0; j<ndim; j++){
            velocity(i, j) = box->rng->next_gaussian(mean, var);
        }
    }
    return velocity;
}

Temp::Temp(class Box* box_in, const double temp_in)
    : Velocity(box_in)
{
    temp = temp_in;
}

mat Temp::get_velocity(const int npar, const int ndim)
{
    double mean = 0.;
    double var = sqrt(temp);
    Gauss gauss(box, mean, var);
    return gauss.get_velocity(npar, ndim);
}
