#pragma once
#include <iostream>
#include <string>
#include <vector>

#include "constraint.h"


class MaxDistance : public Constraint
{
public:
    MaxDistance(class Box *, std::string, std::string, double);
    bool verify() override;
};
