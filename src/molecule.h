#pragma once
#include <vector>
#include <valarray>
#include <string>

class Molecule
{
public:
    Molecule(std::vector<int>, int);

    std::vector<int> atoms_idx;
    int com_atom, natom;
};


class MoleculeTypes
{
public:
    MoleculeTypes(class Box*);
    void add_molecule_type(std::vector<std::string>, double, int, double, std::vector<std::valarray<double> >);
    void check_neighbors(int, int, int, std::vector<int> &);
    Molecule* construct_molecule(int, std::vector<class Particle *>, bool &);

    std::vector<std::vector<std::string> > molecule_elements;
    std::vector<std::vector<std::valarray<double> > > default_mols;
    std::vector<std::vector<int> > molecule_types;
    std::vector<double> molecule_probs, rcs;
    std::vector<int> coms;
    bool configured;
    int ntype;

private:
    class Box* box = nullptr;
};


