#include "init_velocity.h"
#include "rng/rng.h"


Velocity::Velocity() {}

Zero::Zero() : Velocity() {}

mat Zero::get_velocity(const int npar, const int ndim)
{
    return zeros(npar, ndim);
}

Gauss::Gauss(class RandomNumberGenerator* rng_in, const double mean_in, const double var_in) 
    : Velocity()
{
    rng = rng_in;
    mean = mean_in;
    var = var_in;
}

mat Gauss::get_velocity(const int npar, const int ndim)
{
    mat velocity(npar, ndim);
    for(int i=0; i<npar; i++){
        for(int j=0; j<ndim; j++){
            velocity(i, j) = rng->next_gaussian(mean, var);
        }
    }
    return velocity;
}

Temp::Temp(class RandomNumberGenerator* rng_in, const double temp_in)
    : Gauss(rng_in, 0., sqrt(temp_in)) {}

/*
mat Temp::get_velocity(const int npar, const int ndim)
{
    double mean = 0.;
    double var = sqrt(temp);
    Gauss gauss(box, mean, var);
    return gauss.get_velocity(npar, ndim);
}
*/
