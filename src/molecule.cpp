#include <vector>
#include <valarray>
#include <string>
#include <memory>

#include "box.h"
#include "system.h"
#include "molecule.h"
#include "particle.h"
#include "rng/rng.h"


/* ---------------------------------------------------------
   Construct molecule, given atom indices 'atoms_in'. The
   indices refer to particles in the official particle list,
   box->particles. Second argument is center (of mass) atom,
   which is 0 by default (and should probably always be 0).
------------------------------------------------------------ */

Molecule::Molecule(std::vector<int> atoms_in)
{
    atoms_idx = atoms_in;
    natom = atoms_idx.size();
}


/* ---------------------------------------------------------
   Construct object containing all molecule types. This is
   needed when performing AVBMC types of moves with
   molecules. 
------------------------------------------------------------ */

MoleculeTypes::MoleculeTypes(System* system_in)
{
    system = system_in;
    configured = false;
    ntype = 0;
}


/* ---------------------------------------------------------
   Add known molecule type, defined by elements, maximum
   distance between center (of mass) atom and all other
   atoms. Probability of picking this molecule relative 
   to the other molecule types is also needed, and default
   molecule configuration is needed for the input moves.

   Arguments:
     elements :
         vector of molecule atoms/elements
     rc:
         maximum distance between atoms to be considered as molecule
     com_atom:
         index of center atom, needed by AVBMC type of moves
     molecule_prob:
         probability of inserting or removing this molecule
     default_mol:
         molecule to be inserted by AVBMC insertation moves
------------------------------------------------------------ */

void MoleculeTypes::add_molecule_type(std::vector<std::string> elements, const double rc,
                                      const double molecule_prob, std::vector<std::valarray<double> > default_mol)
{
    molecule_elements.push_back(elements);
    rcs.push_back(rc);
    molecule_probs.push_back(molecule_prob);
    default_mols.push_back(default_mol);
    configured = true;
    ntype ++;
}


/* -------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor
   distance squared 'rsq'
---------------------------------------------------------------- */

std::vector<int> MoleculeTypes::build_neigh_list(std::vector<std::shared_ptr<Particle> > particles, const int i, const double rsq)
{
    int npar = particles.size();
    double rijsq;
    std::valarray<double> ri = particles[i]->r;
    std::vector<int> neigh_list;
    for(int j=0; j<i; j++){
        rijsq = std::pow(particles[j]->r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for(int j=i+1; j<npar; j++){
        rijsq = std::pow(particles[j]->r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


/* ---------------------------------------------------------------
   Check if particles match molecule type recursively.
------------------------------------------------------------------ */

void MoleculeTypes::check_neighbors(const int k, const int i, int elm_count,
                                    std::vector<int> &elm_idx, std::vector<std::shared_ptr<Particle> > particles){
    if(elm_count <= molecule_elements[i].size()){  // ensure that recursion stops when molecule has correct size
        if(particles[k]->type == molecule_types[i][elm_count]){  // check if element is matching
            elm_idx.push_back(k);  // add atom to molecule atom idxs
            elm_count ++;
            std::vector<int> neigh_list = build_neigh_list(particles, k, rcs[i]);
            for(int neigh : neigh_list){
                check_neighbors(neigh, i, elm_count, elm_idx, particles);
            }
        }
    }
}


/* ---------------------------------------------------------------
   Construct molecule of type 'i' randomly by picking a random 
   atom among the elements and checking the neighbor list
------------------------------------------------------------------ */

Molecule* MoleculeTypes::construct_molecule(std::vector<std::shared_ptr<Particle> > particles,
                                            const int i, bool &constructed)
{
    int count;
    std::vector<int> elm_idx;
    count = 0;
    while (count < particles.size())
    {
        elm_idx.clear();
        int k = system->rng->next_int(particles.size());     // pick initial particle
        check_neighbors(k, i, 0, elm_idx, particles);
        if(elm_idx.size() == molecule_types[i].size()){
            constructed = true;
            break;
        }
        count ++;
    }
    Molecule* molecule = nullptr; 
    if (!constructed)
    {
        std::vector<int> atoms;
        molecule = new Molecule(atoms);
    }
    else {
        molecule = new Molecule(elm_idx);
    }
    return molecule;
}


