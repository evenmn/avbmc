#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcmoloutres.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../distance_manager.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../constraint/constraint.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Aggregate-volume-biased molecule out move. This is a fictitous inter-box
   move, and works in the grand canonical ensemble only. 'molecule_in'
   specifies what molecule to remove. 'r_above_in' is the maximum distance 
   between target molecule (main atom) and molecule to remove (main atom). 
   'r_max_inner_in' is the maximum distance between main atom and other atoms
   to be considered a molecule. 'energy_bias_in' is a boolean controlling
   whether or not energy biasing should be used. 'target_mol_in' is a boolean
   defining the target molecule as the entire molecule or main atom.
------------------------------------------------------------------------------- */

AVBMCMolOutRes::AVBMCMolOutRes(System* system_in, Box* box_in,
         std::vector<Particle> molecule_in, const double r_above_in,
         const double r_inner_in, const bool energy_bias_in,
         const bool target_mol_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_inner = r_inner_in;
    energy_bias = energy_bias_in;
    target_mol = target_mol_in;
    molecule = molecule_in;
    natom = molecule.size();
    natom_inv = 1. / natom;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    cum_time = 0.;
    label = "AVBMCMolOut";

    neigh_id_above = box->distance_manager->add_cutoff(r_above);
    neigh_id_inner = box->distance_manager->add_cutoff(r_inner);

    for (Particle &particle : molecule) {
        particle.type = box->forcefield->label2type.at(particle.label);
    }
}


/* ----------------------------------------------------------------------------
   Remove a random molecule from the bonded region of another similar molecule.
------------------------------------------------------------------------------- */

void AVBMCMolOutRes::perform_move()
{
    bool detected_out, detected_target;
    unsigned int count, i, j, n_in;
    std::vector<int> target_molecule, neigh_listi;
    std::vector<std::vector<int> > neigh_list_inner;

    reject_move = true;
    if (box->npar > 2 * natom - 1) {
        count = 0;
        while (count < box->npar && reject_move) {  // do maximum npar attempts to detect target molecule
            count ++;
            detected_target = false;
            if (target_mol) {
                neigh_list_inner = box->distance_manager->neigh_lists[neigh_id_inner];
                //target_molecule = detect_molecule(neigh_list_inner, molecule, detected_target);
                target_molecule = detect_molecule(box->particles, molecule, detected_target, r_inner);
                //target_molecule = detect_molecule(neigh_list_inner, box->particles, molecule, detected_target);
                i = target_molecule[0];
            }
            else {
                i = rng->next_int(box->npar);
                detected_target = (box->particles[i].type == molecule[0].type);
            }
            //std::cout << "FJDSJGSRJGJSGJ " << detected_target << std::endl;
            if (detected_target) {  // target molecule detected
                //neigh_listi = box->build_neigh_list(i, r_abovesq);
                neigh_listi = box->distance_manager->neigh_lists[neigh_id_above][i];
                n_in = neigh_listi.size();
                if (n_in >= natom) {  // ensure that there is a least one molecule left
                    std::vector<Particle> particles_tmp;
                    for (j=0; j < n_in; j++){
                        particles_tmp.push_back(box->particles[neigh_listi[j]]);
                    }
                    detected_out = false;

                    //std::vector<int> molecule_out = detect_molecule(
                    //neigh_list_inner = box->distance_manager->neigh_lists[neigh_id_inner];
                    //std::vector<int> molecule_out = detect_molecule(neigh_list_inner, particles_tmp, molecule, detected_out);
                    std::vector<int> molecule_out = detect_molecule(particles_tmp, molecule, detected_out, r_inner);
                    //std::cout << "detected_out: " << detected_out << " " << molecule_out.size() << std::endl;
                    if (detected_out) {
                        std::cout << "DETECTED OUT" << std::endl;
                        std::vector<int> molecule_out2;
                        for (int idx : molecule_out) {
                            molecule_out2.push_back(neigh_listi[idx]);
                        }
                        reject_move = false;
                        // compute change of energy when removing molecule
                        du = 0.;
                        if (box->store_energy) {
                            for (int j : molecule_out2) {
                                du -= box->forcefield->poteng_vec[j];
                            }
                        }
                        else {
                            for (int j : molecule_out2) {
                                du -= box->forcefield->comp_energy_par_force0(j);
                            }
                        }

                        // remove molecule
                        npartype_old = box->npartype;
                        particles_old = box->particles;
                        std::sort(molecule_out2.begin(), molecule_out2.end(), std::greater<int>()); // sort in descending order
                        box->distance_manager->set();
                        std::cout << "NPAR: " << box->npar << std::endl;
                        for (int j : molecule_out2){
                            box->npar --;
                            box->npartype[box->particles[j].type] --;
                            //box->particles[j] = box->particles.back();
                            //box->particles.pop_back();
                            //box->distance_manager->update_remove(j);

                            //box->particles.erase(box->particles.begin() + j);
                        }
                        std::cout << std::endl;
                        box->poteng += du;
                        //box->npar -= natom;
                        nmolavg = n_in * natom_inv;
                    }  // end 
                }  // end if n_in
            }  // end if constructed_target
        }  // end while
    }  // end if npar
}  // end perform_move


/* ----------------------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
------------------------------------------------------------------------------- */

double AVBMCMolOutRes::accept(double temp, double chempot)
{
    bool constr_satis = true;
    for (Constraint* constraint : box->constraints) {
        constr_satis *= constraint->verify();
    }
    if (reject_move || !constr_satis){
        return 0.;
    }
    else {
        //box->boundary->correct_position();
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar + natom);
        double prefactor = nmolavg * box->npar / (v_in * (box->npar - natom)); 
        return prefactor * std::exp(-(du+chempot+dw)/temp);
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state if move is rejected
------------------------------------------------------------------------------- */

void AVBMCMolOutRes::reset()
{
    if (!reject_move) {
        box->distance_manager->reset();
        box->forcefield->reset();
        box->npar += natom;
        box->npartype = npartype_old;
        box->poteng -= du;
        box->particles = particles_old;
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------------------------- */

void AVBMCMolOutRes::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string AVBMCMolOutRes::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Energy bias:            " + std::to_string(energy_bias) + "\n";
    move_info += "    Search target molecule: " + std::to_string(target_mol) + "\n";
    return move_info;
}
