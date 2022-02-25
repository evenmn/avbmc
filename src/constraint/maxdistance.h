#include <iostream>
#include <string>
#include <vector>

#include "constraint.h"


class MaxDistance : public Constraint
{
public:
    MaxDistance(class Box *, std::string, std::string, double);
    bool verify();

private:
    unsigned int cutoff_id, type1, type2;
    std::vector<std::vector<int> > neigh_list;
};
