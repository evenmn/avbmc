#pragma once
#include <iostream>
#include <string>
#include <vector>


class Constraint
{
public:
    Constraint(class Box *);
    Constraint(const Constraint &);
    virtual bool verify() = 0;
    virtual ~Constraint() = default;

    std::string label;

protected:
    unsigned int cutoff_id, type1, type2;
    std::vector<std::vector<unsigned int> > neigh_list;
    class Box* box = nullptr;
};
