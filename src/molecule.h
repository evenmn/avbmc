#pragma once
#include <vector>
#include <valarray>
#include <string>

class Molecule
{
public:
    Molecule(std::vector<int>);
    std::vector<int> atoms_idx;
    int natom;
};


class MoleculeTypes
{
public:
    MoleculeTypes(class System *);
    void add_molecule_type(std::vector<std::string>, double, double,
                           std::vector<std::valarray<double> >);
    std::vector<int> build_neigh_list(std::vector<class Particle>, int, double);
    void check_neighbors(int, int, unsigned int, std::vector<int>&,
                         std::vector<class Particle>);
    Molecule* construct_molecule(std::vector<class Particle>, int, bool&);

    std::vector<std::vector<std::string> > molecule_elements;
    std::vector<std::vector<std::valarray<double> > > default_mols;
    std::vector<std::vector<int> > molecule_types;
    std::vector<double> molecule_probs, rcs;
    bool configured;
    int ntype;

private:
    class System* system = nullptr;
};


