#include "avbmcin.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"


AVBMCIn::AVBMCIn(Box* box_in, const double r_below_in, const double r_above_in)
    : Moves(box_in)
{
    /* Aggregate-volume-biased in move. 
     * This is a fictitious inter-box move, and works in the
     * grand canonical essemble only
     */

    r_below = r_below_in;
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    r_below2 = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    type = 0;      // type and label of inserted particle
    label = "Ar";  // this has to be generalized
}

std::vector<std::valarray<double> > AVBMCIn::rotate_molecule(std::vector<std::valarray<double> > positions_in)
{
    /* Rotate rotate by a random angle in all dimensions. Might send in positions
     * as a reference in the future, even though we still may need to copy the positions
     */
    std::vector<std::valarray<double> > positions_out;
    if(box->ndim == 1){
        return positions_in;
    }
    else if(box->ndim == 2){
        double angle = 2 * pi * box->rng->next_double();
        for(std::valarray<double> position_in : positions_in){
            std::valarray<double> position_out = {position_in[0] * std::cos(angle) - position_in[1] * std::sin(angle), position_in[0] * std::sin(angle) + position_in[1] * std::cos(angle)};
            positions_out.push_back(position_out);
        }
        return positions_out;
    }
    else{
        double anglex = 2 * pi * box->rng->next_double();
        double angley = 2 * pi * box->rng->next_double();
        double anglez = 2 * pi * box->rng->next_double();
        for(std::valarray<double> position_in : positions_in){
            std::valarray<double> position_out = {position_in[0] * (1 + std::cos(angley) + std::cos(anglez)) - position_in[2] * std::sin(anglez) + position_in[1] * std::sin(angley), position_in[0] * std::sin(angley) + position_in[1] * (1 + std::cos(anglex) + std::cos(anglez)) - position_in[2] * std::sin(anglex), - position_in[0] * std::sin(angley) + position_in[1] * std::sin(anglex) + position_in[2] * (1 + std::cos(anglex) + std::cos(angley))};
            positions_out.push_back(position_out);
        }
        return positions_out;
    }
}





void AVBMCIn::perform_move()
{
    /* Insert particle into the bonded region
     * of particle i
     */

    // pick molecule type to be inserted
    int mol_idx = box->rng->choice(box->molecule_types->molecule_probs);
    std::vector<std::valarray<double> > positions = box->molecule_types->default_mols[mol_idx];

    // positions = rotate_molecule(positions)
    
    // pick particle among COM particles
    int type, i;
    int count = 0;
    //while(type != box->particles[box->molecule_types->molecule_types[box->molecule_types->coms[mol_idx]]]->type && count < box->npar){
    //    i = box->rng->next_int(box->npar);
    //    type = box->particles[i]->type;
    //}
    type = 0;
    i = 0;
    std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
    n_in = neigh_listi.size();
    

    // pick particle i and create local neighbor list of particle
    //int i = box->rng->next_int(box->npar);
    //std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
    //n_in = neigh_listi.size();

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
    for(std::valarray<double> &position : positions){
        position += box->particles[i]->r + dr;
    }
    //Molecule* molecule_in = new Molecule(positions, box->molecule_types->coms[mol_idx]);

    //du = - box->forcefield->comp_energy_mol(box->particles, molecule_in);
    du = 1e10;

    //Particle *particle_in = new Particle(label, box->particles[i]->r + dr);
    //particle_in->type = type;
    //box->particles.push_back(particle_in);
    //box->npar += molecule_in->natom;

    // compute du
    //du = box->forcefield->comp_energy_par(box->particles, box->npar - 1);
    //box->poteng += du;
}

double AVBMCIn::accept(double temp, double chempot)
{
    /* Get acceptance probability of move
     */
    double dw = box->sampler->w(box->npar) - box->sampler->w(box->npar-1);
    return (v_in * box->npar) / ((n_in + 1) * (box->npar + 1)) * std::exp(-(du-chempot+dw)/temp);
}

void AVBMCIn::reset()
{
    /* Set back to old state is move
     * is rejected
     */
    //box->npar --;
    //box->poteng -= du;
    //box->particles.erase(box->particles.begin() + box->npar);
}
