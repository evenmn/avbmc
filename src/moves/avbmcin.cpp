#include "avbmcin.h"
#include "../box.h"
#include "../particle.h"


AVBMCIn::AVBMCIn(class Box* box_in, const double r_below_in, const double r_above_in)
    : Moves(box_in)
{
    /* Aggregate-volume-biased out move. 
     * This is a fictitious inter-box move, and work is the
     * grand canonical essemble only
     */

    r_below = r_below_in;
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    r_below2 = r_below * r_below;
    v_in = 4 * pi * std::pow(r_above, 3)/3;

    type = 0;
    label = "Ar";
}

void AVBMCIn::perform_move()
{
    /* Insert particle into the bonded region
     * of particle i
     */

    // create local neighbor list of particle i
    int i = box->rng->next_int(box->npar);
    std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
    n_in = neigh_listi.size();

    // compute norm
    auto norm = [] (std::valarray<double> x) -> double { 
        double sqrd_sum = 0.;
        for(double x_ : x){
            sqrd_sum += x_ * x_;
        }
        return sqrd_sum;
    };

    // construct new particle
    std::valarray<double> dr(box->ndim);
    double norm_ = norm(dr);
    while(norm_ > r_above2 || norm_ < r_below2){
        for(double &d : dr){
            d = 2 * box->rng->next_double() - 1;
        }
        norm_ = norm(dr);
    }
    Particle *particle_in = new Particle(label, box->particles[i]->r + dr);
    particle_in->type = type;
    box->particles.push_back(particle_in);
    box->npar ++;

    // compute du
    box->sampler->du = box->forcefield->comp_energy_par(box->particles, box->npar - 1);
    box->poteng += box->sampler->du;
}

double AVBMCIn::accept(double temp, double chempot)
{
    /* Get pre-exponential factor of Boltzmann
     * ratio
     */
    return (v_in * box->npar) / ((n_in + 1) * (box->npar + 1)) * std::exp(chempot/temp);
}

void AVBMCIn::reset()
{
    box->npar --;
    box->poteng -= box->sampler->du;
    box->particles.erase(box->particles.begin() + box->npar);
}
