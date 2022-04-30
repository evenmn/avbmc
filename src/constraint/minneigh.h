#pragma once
#include <iostream>
#include <string>
#include <vector>

#include "constraint.h"


class MinNeigh : public Constraint
{
public:
    MinNeigh(class Box *, std::string, std::string, double, int);
    bool verify() override;

private:
    unsigned int nc;
};
