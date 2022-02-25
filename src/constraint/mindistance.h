#include <iostream>
#include <string>
#include <vector>

#include "constraint.h"


class MinDistance : public Constraint
{
public:
    MinDistance(class Box *, std::string, std::string, double);
    bool verify();
};
