#include <vector>
#include <valarray>
#include <string>

#include "box.h"
#include "molecule.h"
#include "particle.h"


/* ---------------------------------------------------------
   Construct molecule, given atom indices 'atoms_in'. The
   indices refer to particles in the official particle list,
   box->particles. Second argument is center (of mass) atom,
   which is 0 by default (and should probably always be 0).
------------------------------------------------------------ */

Molecule::Molecule(std::vector<int> atoms_in, const int com_atom_in)
{
    atoms_idx = atoms_in;
    com_atom = com_atom_in;
    natom = atoms_idx.size();
}


/* ---------------------------------------------------------
   Construct object containing all molecule types. This is
   needed when performing AVBMC types of moves with
   molecules. 
------------------------------------------------------------ */

MoleculeTypes::MoleculeTypes(Box* box_in)
{
    box = box_in;
    configured = false;
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

void MoleculeTypes::add_molecule_type(std::vector<std::string> elements, const double rc, const int com_atom,
                                      const double molecule_prob, std::vector<std::valarray<double> > default_mol)
{
    molecule_elements.push_back(elements);
    rcs.push_back(rc);
    coms.push_back(com_atom);
    molecule_probs.push_back(molecule_prob);
    default_mols.push_back(default_mol);
    configured = true;
}


/* ---------------------------------------------------------------
   Check if particles match molecule type recursively.
------------------------------------------------------------------ */

void MoleculeTypes::check_neighbors(const int k, const int i, int elm_count, std::vector<int> &elm_idx){
    if(elm_count <= molecule_elements[i].size()){  // ensure that recursion stops when molecule has correct size
        if(box->particles[k]->type == molecule_types[i][elm_count]){  // check if element is matching
            elm_idx.push_back(k);  // add atom to molecule atom idxs
            elm_count ++;
            std::vector<int> neigh_list = box->forcefield->build_neigh_list(k, rcs[i]);
            for(int neigh : neigh_list){
                check_neighbors(neigh, i, elm_count, elm_idx);
            }
        }
    }
}


/* ---------------------------------------------------------------
   Construct molecule of type 'i' randomly by picking a random 
   atom among the elements and checking the neighbor list
------------------------------------------------------------------ */

Molecule* MoleculeTypes::construct_molecule(const int i)
{
    std::vector<int> elm_idx;
    int count = 0;
    while(count < box->npar){
        elm_idx.clear();
        int k = box->rng->next_int(box->npar);     // pick particle
        check_neighbors(k, i, 0, elm_idx);
        if(elm_idx.size() == molecule_types[i].size()){
            break;
        }
        count ++;
    }
    Molecule* molecule = nullptr; 
    if(count == box->npar - 1){
        std::vector<int> atoms;
        molecule = new Molecule(atoms, 0);
    }
    else{
        molecule = new Molecule(elm_idx, coms[i]);
    }
    return molecule;
}


